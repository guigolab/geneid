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

package org.ensembl.driver.plugin.standard;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.AffyProbe;
import org.ensembl.datamodel.impl.AffyProbeImpl;
import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AffyProbeAdaptor;

/**
 * The point of this class is....
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class MySQLAffyProbeAdaptor extends MySQLBaseFeatureAdaptor implements
		Adaptor, AffyProbeAdaptor {

	/**
	 * @param driver
	 */
	public MySQLAffyProbeAdaptor(MySQLDriver driver) {
		super(driver, TYPE);
	}

	/**
	 * @see org.ensembl.driver.AffyProbeAdaptor#fetch(long)
	 */
	public AffyProbe fetch(long internalID) throws AdaptorException {
		return (AffyProbe) fetchByInternalID(internalID);
	}

	/**
	 * @see org.ensembl.driver.AffyProbeAdaptor#fetch(org.ensembl.datamodel.AffyArray)
	 */
	public List fetch(AffyArray array) throws AdaptorException {
		return fetchByNonLocationConstraint(" ap.affy_array_id = "
				+ array.getInternalID());
	}

	/**
	 * @see org.ensembl.driver.AffyProbeAdaptor#fetch(java.lang.String)
	 */
	public List fetch(String probeSetName) throws AdaptorException {
		return fetchByNonLocationConstraint(" ap.probeset = \"" + probeSetName + "\"");
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#tables()
	 */
	public String[][] tables() {
		final String[][] tables = { { "affy_probe", "ap" } };
		return tables;
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#columns()
	 */
	public String[] columns() {
		final String[] columns = { "ap.affy_probe_id", "ap.affy_array_id",
				"ap.probeset", "ap.name" };
		return columns;
	}

	public String finalClause() {
		return " order by ap.affy_probe_id, ap.affy_array_id";
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
	 */
public Object createObject(ResultSet rs) throws AdaptorException {

		AffyProbe p = null;

		try {
			
			while(rs.next()) {
        long internalID = rs.getLong(1);
        
        if (p==null) {
          p = new AffyProbeImpl(driver, 
          			internalID,
								rs.getString(3));
        }
				
        if (internalID==p.getInternalID())
          p.addArrayWithProbeName(driver.getAffyArrayAdaptor().fetch(rs.getLong(2)), 
          		rs.getString(4));
			}
      
      // reset the row cursor for when this method is called again for the
      // next probe
			if (p!=null)
				rs.previous();

		} catch (SQLException e) {
			throw new AdaptorException("Problem loading AffyProbe from database.",e); 
		}
		return p;
	}
}
