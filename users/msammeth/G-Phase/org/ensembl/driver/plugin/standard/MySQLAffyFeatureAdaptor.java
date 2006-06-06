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
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.AffyFeature;
import org.ensembl.datamodel.AffyProbe;
import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.impl.AffyFeatureImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AffyFeatureAdaptor;

/**
 * The point of this class is....
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class MySQLAffyFeatureAdaptor extends MySQLBaseFeatureAdaptor implements
		AffyFeatureAdaptor {

	/**
	 * @param driver
	 * @param type
	 */
	public MySQLAffyFeatureAdaptor(MySQLDriver driver) {
		super(driver, TYPE);
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#tables()
	 */
	public String[][] tables() {
		final String[][] tables = { { "affy_feature", "af" },
				{ "affy_probe", "ap" }, { "affy_array", "aa" } };
		return tables;
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#columns()
	 */
	public String[] columns() {
		final String[] columns = { "af.affy_feature_id", "af.seq_region_id",
				"af.seq_region_start", "af.seq_region_end",
				"af.seq_region_strand", "af.mismatches", "af.affy_probe_id",
				"af.analysis_id", "aa.name", "ap.probeset", "ap.name" };
		return columns;
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#finalWhereClause()
	 */
	public String finalWhereClause() {
		return "af.affy_probe_id = ap.affy_probe_id and ap.affy_array_id=aa.affy_array_id";
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
	 */
	public Object createObject(ResultSet rs) throws AdaptorException {
		try {
			if (!rs.next())
				return null;

			Location loc = locationConverter
					.idToLocation(rs.getLong("seq_region_id"), rs
							.getInt("seq_region_start"), rs
							.getInt("seq_region_end"), rs
							.getInt("seq_region_strand"));

			return new AffyFeatureImpl(driver, rs.getLong(1), loc, rs
					.getString(10), rs.getLong(7), rs.getInt(6));

		} catch (SQLException e) {
			throw new AdaptorException("Failed to create AffyFeature.", e);
		}
	}

	/**
	 * @throws AdaptorException
	 * @see org.ensembl.driver.AffyFeatureAdaptor#fetch(org.ensembl.datamodel.AffyProbe)
	 */
	public List fetch(AffyProbe affyProbe) throws AdaptorException {
		return fetchByNonLocationConstraint("af.affy_probe_id = "
				+ affyProbe.getInternalID());
	}

	/**
	 * @throws AdaptorException
	 * @see org.ensembl.driver.AffyFeatureAdaptor#fetch(org.ensembl.datamodel.AffyArray)
	 */
	public List fetch(AffyArray array) throws AdaptorException {
		return fetchByNonLocationConstraint("aa.affy_array_id = "
				+ array.getInternalID());
	}

	/**
	 * @throws AdaptorException
	 * @see org.ensembl.driver.AffyFeatureAdaptor#fetch(org.ensembl.datamodel.Location,
	 *      org.ensembl.datamodel.AffyArray)
	 */
	public List fetch(Location loc, AffyArray array) throws AdaptorException {
		return fetchAllByConstraint(loc, "aa.affy_array_id = "
				+ array.getInternalID());
	}

	/**
	 * @see org.ensembl.driver.AffyFeatureAdaptor#fetch(long)
	 */
	public AffyFeature fetch(long internalID) throws AdaptorException {
		return (AffyFeature) fetchByInternalID(internalID);
	}

	/**
	 * @throws AdaptorException
	 * @see org.ensembl.driver.AffyFeatureAdaptor#fetchUniqueHits(org.ensembl.datamodel.Location)
	 */
	public List fetchUniqueProbeAndLocation(Location location)
			throws AdaptorException {

		// "Squash" duplicate affy_features that
		// differ only by their microarray 
		
		// Note: possible optimisation of squashing via SQL "group by"
		// proved too slow for big regions e.g. chr1 in homo_sapiens_core_27_35a.
		
		Location requestLoc = driver.getLocationConverter().fetchComplete(
				location);
		if (requestLoc == null)
			return Collections.EMPTY_LIST;

		List r = new ArrayList();
		
		CoordinateSystem[] relevantCSs = coordinateSystemAdaptor
				.fetchAllByFeatureTable(getPrimaryTableName());
		
		for (int i = 0; i < relevantCSs.length; i++) {

			// TODO - dereference
			// vvvvvvvvvv This should be dereferencedLoc
			Location loc = locationConverter
					.convert(requestLoc, relevantCSs[i]);

			// Load data for one seq region at a time
			// in an effort to reduce total memory footprint.
			// We can discard duplicate affy features on per region
			// basis
			for (Location node = loc; node != null; node = node
					.next()) {
				
				if (node.getSeqRegionName()==null) continue;
						
				// ensure only a single seq region is used as filter per query
				Location query = node.copy();
				query.setNext(null);

				// we rely on the fact that the affy features are sorted by location
				// and then probeInternalID for the filtering we do next.
				List tmp = fetch(query);
				AffyFeature previous = null;
				for (int j = 0, n = tmp.size(); j < n; j++) {
					AffyFeature af = (AffyFeature) tmp.get(j);
					
					if (previous!=null && af.getProbeInternalID()==previous.getProbeInternalID() 
							&& af.getLocation().compareTo(previous.getLocation())==0 )
						continue;
					
					r.add(af);
					previous = af;
				}
			
			}
		}

		return r;
	}
}
