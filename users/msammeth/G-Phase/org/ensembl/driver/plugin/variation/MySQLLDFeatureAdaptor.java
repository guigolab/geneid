/*
	Copyright (C) 2003 EBI, GRL

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
import org.ensembl.datamodel.variation.LDFeature;
import org.ensembl.datamodel.variation.LDFeatureContainer;
import org.ensembl.datamodel.variation.VariationFeature;
import org.ensembl.datamodel.variation.impl.LDFeatureImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor;
import org.ensembl.driver.plugin.standard.MySQLDriver;
import org.ensembl.driver.variation.LDFeatureAdaptor;
import org.ensembl.driver.variation.VariationDriver;

/**
 * Fetches LDFeatures from an ensembl database as either a simple
 * List or an LDFeatureContainer.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class MySQLLDFeatureAdaptor extends MySQLBaseFeatureAdaptor implements
		LDFeatureAdaptor {

	private VariationDriver vdriver;

	/**
	 * @param vdriver parent driver
	 */
	public MySQLLDFeatureAdaptor(VariationDriver vdriver) {
		super((MySQLDriver)vdriver.getCoreDriver(),TYPE);
    this.vdriver = vdriver;
	}

  /**
   * @see org.ensembl.driver.plugin.standard.BaseAdaptor#getConnection()
   */
  public Connection getConnection() throws AdaptorException {
    return vdriver.getConnection();
  }

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#tables()
	 */
	public String[][] tables() {
    final String[][] tables = new String[][]{{"pairwise_ld", "pl"}};
      return tables;
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#columns()
	 */
	public String[] columns() {
		final String[] columns = {"pl.variation_feature_id_1", "pl.variation_feature_id_2", "pl.population_id",
                              "pl.seq_region_id", "pl.seq_region_start", "pl.seq_region_end", 
                              "pl.snp_distance_count", "pl.r2", "pl.d_prime", "pl.sample_count"};
		return columns;
	}

	/**
   * @return an LDFeature created from next row in rs if available, null if no more rows.
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
	 */
	public Object createObject(ResultSet rs) throws AdaptorException {
    LDFeature vf = null;

    try {
      if (rs.next()) {
        Location loc =
          vdriver.getCoreDriver().getLocationConverter().idToLocation(
            rs.getLong("seq_region_id"),
            rs.getInt("seq_region_start"),
            rs.getInt("seq_region_end"),
            0); // unstranded information

        vf = new LDFeatureImpl(vdriver,
						loc,
						rs.getLong("pl.variation_feature_id_1"),
						rs.getLong("pl.variation_feature_id_2"),
						rs.getLong("pl.population_id"),
						rs.getInt("pl.snp_distance_count"),
						rs.getDouble("pl.r2"),
						rs.getDouble("pl.d_prime"),
						rs.getInt("pl.sample_count")
						);
        
      }

    } catch (InvalidLocationException ee) {
      throw new AdaptorException("Error when building Location", ee);
    } catch (SQLException se) {
      throw new AdaptorException("SQL error when building object", se);
    }

    return vf;
	}

	/**
	 * @see org.ensembl.driver.variation.LDFeatureAdaptor#fetch(org.ensembl.datamodel.variation.VariationFeature)
	 */
	public List fetch(VariationFeature variationFeature)
			throws AdaptorException {
		return fetchByNonLocationConstraint("pl.variation_feature_id_1 = "+variationFeature.getInternalID());
	}

	/**
	 * @see org.ensembl.driver.variation.LDFeatureAdaptor#fetchLDFeatureContainer(org.ensembl.datamodel.variation.VariationFeature)
	 */
	public LDFeatureContainer fetchLDFeatureContainer(
			VariationFeature variationFeature) throws AdaptorException {
    List l = fetch(variationFeature);
		return new LDFeatureContainer(l);
	}

	/**
	 * @see org.ensembl.driver.variation.LDFeatureAdaptor#fetchLDFeatureContainer(org.ensembl.datamodel.Location)
	 */
	public LDFeatureContainer fetchLDFeatureContainer(Location location)
			throws AdaptorException {
    List l = fetch(location);
		return new LDFeatureContainer(l);
	}

}
