/*
 Copyright (C) 2001 EBI, GRL

 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.

 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
package org.ensembl.driver.plugin.variation;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.ensembl.datamodel.variation.Allele;
import org.ensembl.datamodel.variation.ValidationState;
import org.ensembl.datamodel.variation.Variation;
import org.ensembl.datamodel.variation.impl.AlleleImpl;
import org.ensembl.datamodel.variation.impl.VariationImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.plugin.standard.InternalIDOrderComparator;
import org.ensembl.driver.plugin.standard.MySQLDriver;
import org.ensembl.driver.variation.VariationAdaptor;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.util.StringUtil;
import org.ensembl.util.Util;

/**
 * Provides access to variations in an ensembl variation database.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 * @see MySQLVariationDriver parent driver this adaptor usually belongs to.
 */
public class MySQLVariationAdaptor implements VariationAdaptor {

  private VariationDriver vdriver;

  public MySQLVariationAdaptor(VariationDriver vdriver) {
    this.vdriver = vdriver;

  }

  /**
   * @see org.ensembl.driver.variation.VariationAdaptor#fetch(long)
   */
  public Variation fetch(long internalID) throws AdaptorException {

    return fetchByConstraint("v.variation_id = " + internalID);
  }

  private Variation fetchByConstraint(String constraint)
      throws AdaptorException {
    Variation v = null;
    List tmp = fetchListByConstraint(constraint);
    if (tmp.size() > 0)
      v = (Variation) tmp.get(0);
    return v;
  }

  /**
   * Creates a variation object from the current and next N rows.
   * 
   * @param rs
   *          resultset with next() called at least once.
   * @return a variation or null if after last row in rs.
   */
  private Variation createObject(ResultSet rs) throws AdaptorException,
      SQLException {

    if (rs.isAfterLast())
      return null;

    final long internalID = rs.getLong("v.variation_id");
    Variation v = new VariationImpl(vdriver);
    v.setInternalID(internalID);
    v.setName(rs.getString("v.name"));
    String stateString = rs.getString("v.validation_status");
    String[] states = stateString.split(",");
    for (int i = 0; i < states.length; i++) {
      String code = states[i];
      if (code.length() == 0)
        continue;
      ValidationState s = ValidationState.createValidationState(code);
      v.addValidationState(s);
    }

    Set alleleSet = new HashSet();
    Set synonymSet = new HashSet();

    // TODO v.setSource(rs.getString());
    do {

      long alleleID = rs.getLong("a.allele_id");
      Long alleleIDKey = new Long(alleleID);
      if (!alleleSet.contains(alleleIDKey)) {
        alleleSet.add(alleleIDKey);
        Allele a = new AlleleImpl(vdriver);
        v.addAllele(a);
        a.setInternalID(rs.getLong("a.allele_id"));
        a.setAlleleString(rs.getString("a.allele"));
        a.setFrequency(rs.getDouble("a.frequency"));
        a.setPopulation(vdriver.getPopulationAdaptor().fetch(
            rs.getLong("a.population_id")));
      }

      String synonym = rs.getString("vs.name") + ":" + rs.getString("s2.name");
      if (!synonymSet.contains(synonym)) {
        synonymSet.add(synonym);
        v.addSynonym(synonym);
      }

    } while (rs.next() && rs.getLong("v.variation_id") == internalID);
    return v;
  }

  /**
   * @return VariationAdaptor.TYPE
   * @see org.ensembl.driver.Adaptor#getType()
   */
  public String getType() throws AdaptorException {
    return TYPE;
  } /*
     * (non-Javadoc)
     * 
     * @see org.ensembl.driver.Adaptor#closeAllConnections()
     */

  public void closeAllConnections() throws AdaptorException { // TODO
    // Auto-generated
    // method stub

  } /*
     * (non-Javadoc)
     * 
     * @see org.ensembl.driver.Adaptor#clearCache()
     */

  public void clearCache() throws AdaptorException { // TODO Auto-generated
    // method stub

  }

  /**
   * @see org.ensembl.driver.variation.VariationAdaptor#fetch(java.lang.String)
   */
  public Variation fetch(String name) throws AdaptorException {
    return fetchByConstraint("v.name = \"" + name + "\"");
  }

  /**
   * @see org.ensembl.driver.variation.VariationAdaptor#fetch(long[])
   */
  public List fetch(long[] internalIDs) throws AdaptorException {
    List r = new ArrayList();
    // TODO optimise with batch fetch
    //    for (int i = 0, n = internalIDs.length; i < n; i++) {
    //      Variation v = fetch(internalIDs[i]);
    //      if (v!=null)
    //        r.add(v);
    //    }

    long[][] idBatches = Util.batch(internalIDs, 100);
    for (int i = 0; i < idBatches.length; i++) {
      r.addAll(fetchBatch(idBatches[i]));
    }

    Collections.sort(r, new InternalIDOrderComparator(internalIDs));
    return r;
  }

  private List fetchBatch(long[] internalIDs) throws AdaptorException {
    if (internalIDs.length == 0)
      return Collections.EMPTY_LIST;

    StringBuffer sql = new StringBuffer();
    sql.append("v.variation_id IN (");
    sql.append(StringUtil.toString(internalIDs));
    sql.append(")");
    return fetchListByConstraint(sql.toString());
  }

  /**
   * @param string
   * @return
   * @throws AdaptorException
   */
  private List fetchListByConstraint(String constraint) throws AdaptorException {

    List r = new ArrayList();

    String sql = "SELECT v.variation_id, v.name, v.validation_status, s1.name,"
        + " a.allele_id, a.allele, a.frequency, a.population_id,"
        + " vs.name, s2.name"
        + " FROM   variation v, source s1, source s2, allele a, variation_synonym vs"
        + " WHERE  v.variation_id = a.variation_id"
        + " AND    v.variation_id = vs.variation_id"
        + " AND    v.source_id = s1.source_id"
        + " AND    vs.source_id = s2.source_id" + " AND " + constraint
        + " ORDER BY v.variation_id ";

    Connection conn = null;
    try {
      conn = vdriver.getConnection();
      ResultSet rs = conn.createStatement().executeQuery(sql);
      rs.next();
      Object o = null;
      while ((o = createObject(rs)) != null) {
        r.add(o);
      }
    } catch (SQLException e) {
      throw new AdaptorException("Failed to fetch Variation with constriant: "
          + constraint + ":"+sql, e);
    } finally {
      MySQLDriver.close(conn);
    }

    return r;
  }

}
