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
import java.util.List;

import org.ensembl.datamodel.variation.Population;
import org.ensembl.datamodel.variation.impl.PopulationImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.plugin.standard.MySQLDriver;
import org.ensembl.driver.variation.PopulationAdaptor;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.util.LruCache;

/**
 * Implementation of Population adaptor.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public class MySQLPopulationAdaptor implements PopulationAdaptor {

  private VariationDriver vdriver;

  /**
   * LRU cache is used to store populations and thereby
   * reduce database retrievals.
   */
  private LruCache cache = new LruCache(1000);

  /**
   * @param driver
   */
  public MySQLPopulationAdaptor(VariationDriver vdriver) {
    this.vdriver = vdriver;
  }

  /**
   * @see org.ensembl.driver.variation.PopulationAdaptor#fetch(long)
   */
  public Population fetch(long internalID) throws AdaptorException {

    Population p = (Population) cache.get(internalID);
    if (p == null)
      p = fetchByConstraint("population_id = " + internalID);
    return p;
  }

  /**
   * @return PopulationAdaptor.TYPE
   * @see org.ensembl.driver.Adaptor#getType()
   * @see org.ensembl.driver.variation.PopulationAdaptor#TYPE
   */
  public String getType() throws AdaptorException {
    return TYPE;
  }

  /**
   * @see org.ensembl.driver.Adaptor#closeAllConnections()
   */
  public void closeAllConnections() throws AdaptorException {
    // we just use the standard connections from the driver
    // so close it there.
    vdriver.closeAllConnections();
  }

  /**
   * @see org.ensembl.driver.Adaptor#clearCache()
   */
  public void clearCache() throws AdaptorException {
    cache.clear();
  }

  /**
   * @see org.ensembl.driver.variation.PopulationAdaptor#fetch(java.lang.String)
   */
  public Population fetch(String name) throws AdaptorException {

    Population p = (Population) cache.get(name);
    if (p == null)
      p = fetchByConstraint("name = \"" + name + "\"");
    return p;
  }

  /**
   * @see org.ensembl.driver.variation.PopulationAdaptor#fetchSuperPopulations(org.ensembl.datamodel.variation.Population)
   */
  public List fetchSuperPopulations(Population subPopulation)
    throws AdaptorException {

    String sql =
      "SELECT p.population_id, p.name, p.size,  p.description"
        + " FROM   population p, population_structure ps"
        + " WHERE  p.population_id = ps.super_population_id"
        + " AND ps.sub_population_id = "
        + subPopulation.getInternalID();

    return fetchListByQuery(sql);
  }

  /**
   * @see org.ensembl.driver.variation.PopulationAdaptor#fetchSubPopulations(org.ensembl.datamodel.variation.Population)
   */
  public List fetchSubPopulations(Population superPopulation)
    throws AdaptorException {

    String sql =
      "SELECT p.population_id, p.name, p.size,  p.description"
        + " FROM   population p, population_structure ps"
        + " WHERE  p.population_id = ps.super_population_id"
        + " AND ps.super_population_id = "
        + superPopulation.getInternalID();

    return fetchListByQuery(sql);
  }

  private List fetchListByConstraint(String constraint)
    throws AdaptorException {
    String sql =
      "SELECT population_id, name, size, description"
        + " FROM   population"
        + " WHERE "
        + constraint;
    return fetchListByQuery(sql);
  }

  private List fetchListByQuery(String sql) throws AdaptorException {
    List r = new ArrayList();

    Connection conn = null;
    try {
      conn = vdriver.getConnection();
      ResultSet rs = conn.createStatement().executeQuery(sql);
      if (rs.next()) {
        Population p = null;
        while ((p = createObject(rs)) != null)
          r.add(p);
      }
    } catch (SQLException e) {
      throw new AdaptorException(
        "Failed to fetch populations with query: " + sql,
        e);
    } finally {
      MySQLDriver.close(conn);
    }

    return r;
  }

  private Population fetchByConstraint(String constraint)
    throws AdaptorException {

    List r = fetchListByConstraint(constraint);
    return (Population) (r.size() > 0 ? r.get(0) : null);

  }

  /**
   * Creates a population object from the current and next N rows.
   * @param rs resultset with next() called at least once.
   * @return a population or null if after last row in rs.
   */

  private Population createObject(ResultSet rs)
    throws AdaptorException, SQLException {

    if (rs.isAfterLast())
      return null;

    final long internalID = rs.getLong("population_id");
    Population p = new PopulationImpl(vdriver);

    do {
      p.setInternalID(internalID);
      p.setName(rs.getString("name"));
      p.setDescription(rs.getString("description"));
      p.setSize(rs.getInt("size"));
    } while (rs.next() && rs.getLong("population_id") != internalID);

    cache.put(p, new Long(p.getInternalID()), p.getName());

    return p;
  }
}
