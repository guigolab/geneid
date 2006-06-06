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
import java.util.List;

import org.ensembl.datamodel.InvalidLocationException;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.variation.Variation;
import org.ensembl.datamodel.variation.VariationFeature;
import org.ensembl.datamodel.variation.impl.VariationFeatureImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor;
import org.ensembl.driver.plugin.standard.MySQLDriver;
import org.ensembl.driver.variation.VariationFeatureAdaptor;

/**
 * Adaptor for accessing variation features from database.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 * @see org.ensembl.datamodel.variation.VariationFeature
 *
 */
public class MySQLVariationFeatureAdaptor
  extends MySQLBaseFeatureAdaptor
  implements VariationFeatureAdaptor {

  // TODO use vdriver.getCoreDriver() to get on demand because it can't be passed
  // in via the generic DriverManager.loadVariationAdaptor() mechanism.

  private MySQLVariationDriver vdriver;

  public MySQLVariationFeatureAdaptor(MySQLVariationDriver vdriver) {
    // TODO consider changing baseadaptor to handle Driver instead of MySQLDriver
    super((MySQLDriver)vdriver.getCoreDriver(), TYPE);
    this.vdriver = vdriver;
  }

  /**
   * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#tables()
   */
  public String[][] tables() {
    String[][] tables = {
      {"variation_feature", "vf"}
    };

    return tables;

  }

  /**
   * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#columns()
   */
  public String[] columns() {
    String[] cols =
  {
    "vf.variation_feature_id",
    "vf.seq_region_id",
    "vf.seq_region_start",
    "vf.seq_region_end",
    "vf.seq_region_strand",
    "vf.variation_id",
    "vf.allele_string",
    "vf.variation_name",
    "vf.map_weight"
  };
    return cols;
  }
    
  /**
   * Creates a VariationFeature using the next row from the result set.
   * @return a VariationFeature if there is another row, otherwise null.
   * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
   */
  public Object createObject(ResultSet rs) throws AdaptorException {
    VariationFeature vf = null;

    try {
      if (rs.next()) {
        Location loc =
          vdriver.getCoreDriver().getLocationConverter().idToLocation(
            rs.getLong("seq_region_id"),
            rs.getInt("seq_region_start"),
            rs.getInt("seq_region_end"),
            rs.getInt("seq_region_strand"));

        vf = new VariationFeatureImpl(vdriver);
        vf.setLocation(loc);
        vf.setInternalID(rs.getLong("variation_feature_id"));
        vf.setAlleleString(rs.getString("allele_string"));
        vf.setVariationName(rs.getString("variation_name"));
        vf.setMapWeight(rs.getInt("map_weight"));
        vf.setVariationInternalID(rs.getLong("variation_id"));
        
      }

    } catch (InvalidLocationException ee) {
      throw new AdaptorException("Error when building Location", ee);
    } catch (SQLException se) {
      throw new AdaptorException("SQL error when building object", se);
    }

    return vf;

  }

  /**
   * @see org.ensembl.driver.variation.VariationFeatureAdaptor#fetch(long)
   */
  public VariationFeature fetch(long internalID) throws AdaptorException {
    return (VariationFeature) super.fetchByInternalID(internalID);
  }

  /**
   * @see org.ensembl.driver.variation.VariationFeatureAdaptor#fetch(org.ensembl.datamodel.variation.Variation)
   */
  public List fetch(Variation variation) throws AdaptorException {
    // TODO assign features to variation?
    return genericFetch("vf.variation_id = "+variation.getInternalID(),null);
  }

  



  /**
   * @see org.ensembl.driver.plugin.standard.BaseAdaptor#getConnection()
   */
  public Connection getConnection() throws AdaptorException {
    return vdriver.getConnection();
  }

}
