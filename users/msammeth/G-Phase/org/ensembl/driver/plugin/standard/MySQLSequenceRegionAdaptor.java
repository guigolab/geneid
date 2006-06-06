/*
 * Copyright (C) 2002 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at your
 * option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package org.ensembl.driver.plugin.standard;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.logging.Logger;

import org.ensembl.datamodel.Attribute;
import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.datamodel.impl.AttributeImpl;
import org.ensembl.datamodel.impl.SequenceRegionImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoordinateSystemAdaptor;
import org.ensembl.driver.SequenceRegionAdaptor;

/**
 * Implementation of SequenceRegionAdaptor interface.
 */

public class MySQLSequenceRegionAdaptor extends BaseAdaptor implements SequenceRegionAdaptor {

	private static final Logger logger = Logger.getLogger(MySQLSequenceRegionAdaptor.class.getName());
	private CoordinateSystemAdaptor coordSystemAdaptor = null;
	public MySQLSequenceRegionAdaptor(MySQLDriver driver) throws AdaptorException {

		super(driver);
		if (coordSystemAdaptor == null) {
			try {
				coordSystemAdaptor = driver.getCoordinateSystemAdaptor();
			} catch (AdaptorException ae) {
				logger.severe("Error creating CoordinateSystemAdaptor " + ae.getMessage());
			}
		}

	}

  
  
  
	public SequenceRegion fetch(long internalID) throws AdaptorException {

		String sql = "SELECT * FROM seq_region WHERE seq_region_id=" + internalID;

		return fetchOne(sql);
	}

	public SequenceRegion fetch(String name, CoordinateSystem cs) throws AdaptorException {

		if (cs.getInternalID() == 0) {
			cs = coordSystemAdaptor.fetch(cs.getName(), cs.getVersion());
		}

		String sql = "SELECT * FROM seq_region WHERE name='" + name + "' AND coord_system_id=" + cs.getInternalID();

		return fetchOne(sql);

	}

	public SequenceRegion fetch(Location loc) throws AdaptorException {

		return fetch(loc.getSeqRegionName(), loc.getCoordinateSystem());

	}

	public SequenceRegion[] fetchAllByCoordinateSystem(CoordinateSystem cs) throws AdaptorException {

		ArrayList result = new ArrayList();

		if (cs.getInternalID() == 0) {
			cs = coordSystemAdaptor.fetch(cs.getName(), cs.getVersion());
		}

		String sql = "SELECT * FROM seq_region WHERE coord_system_id=" + cs.getInternalID();

		return fetchAll(sql);

	}

	public SequenceRegion[] fetchAllByAttributeCode(String code) throws AdaptorException {

		String sql =
			"SELECT sr.seq_region_id, sr.name, sr.length, sr.coord_system_id "
				+ "FROM seq_region sr, attrib_type at, seq_region_attrib sra "
				+ "WHERE sr.seq_region_id=sra.seq_region_id AND sra.attrib_type_id=at.attrib_type_id AND at.code='"
				+ code
				+ "' ";

		return fetchAll(sql);

	}

	public SequenceRegion[] fetchAllByAttributeValue(String code, String value) throws AdaptorException {

		String sql =
			"SELECT sr.seq_region_id, sr.name, sr.length, sr.coord_system_id "
				+ "FROM seq_region sr, attrib_type at, seq_region_attrib sra "
				+ "WHERE sr.seq_region_id=sra.seq_region_id AND sra.attrib_type_id=at.attrib_type_id AND at.code='"
				+ code
				+ "' AND sra.value='"
				+ value
				+ "' ";

		return fetchAll(sql);
	}


  public void fetchComplete(SequenceRegion seqRegion) throws AdaptorException {
    
    Connection conn = null;
    conn = getConnection();
    fetchAttributes(seqRegion, conn);
    close(conn);
  }

	public String getType() throws AdaptorException {

		return TYPE;
	}

	/**
	 * Create a sequence region from a ResultSet row; does fetching of co-ord system and attributes.
	 * @param rs result set from which we use the current row
	 * @return sequence region
	 * @throws SQLException
	 */
	private SequenceRegion createSequenceRegion(ResultSet rs) throws SQLException, AdaptorException {

		SequenceRegion result = new SequenceRegionImpl(getDriver());
		result.setInternalID(rs.getLong("seq_region_id"));
		result.setName(rs.getString("name"));
		result.setLength(rs.getLong("length"));

		// Coordinate system
		long csID = rs.getLong("coord_system_id");
		if (coordSystemAdaptor == null) {
			try {
				coordSystemAdaptor = driver.getCoordinateSystemAdaptor();
			} catch (AdaptorException ae) {
				logger.severe("Error creating CoordinateSystemAdaptor " + ae.getMessage());
			}
		}

		result.setCoordinateSystem(coordSystemAdaptor.fetch(csID));

		return result;
	}

	/**
	* Fetch the attributes for a particular SequenceRegion object from the database.
	 */
	private void fetchAttributes(SequenceRegion sr, Connection con) throws AdaptorException {

		ResultSet rs;

		String sql =
			"SELECT * FROM seq_region_attrib sra, attrib_type at "
				+ "WHERE at.attrib_type_id=sra.attrib_type_id AND sra.seq_region_id="
				+ sr.getInternalID();

		try {

			con = getConnection();
			rs = con.createStatement().executeQuery(sql);

			while (rs.next()) {

				Attribute sra =
					new AttributeImpl(
						rs.getString("code"),
						rs.getString("name"),
						rs.getString("description"),
						rs.getString("value"));

				sr.addAttribute(sra);
			}

		} catch (Exception e) {
			throw new AdaptorException("Rethrow + stacktrace" + sql, e);
		} finally {
			close(con);
		}

	}

	private SequenceRegion fetchOne(String sql) throws AdaptorException {

		SequenceRegion result = null;
		Connection con = null;

		try {

			con = getConnection();
			ResultSet rs = con.createStatement().executeQuery(sql);

			if (rs.next()) {
				result = createSequenceRegion(rs); // populates basic fields and coord-system
				fetchAttributes(result, con); // populates attributes
			}

		} catch (Exception e) {
			throw new AdaptorException("Rethrow + stacktrace" + sql, e);
		} finally {
			close(con);
		}

		return result;

	}

	private SequenceRegion[] fetchAll(String sql) throws AdaptorException {

		ArrayList result = new ArrayList();

		ResultSet rs;
		Connection con = null;

		try {

			con = getConnection();
			rs = con.createStatement().executeQuery(sql);

			while (rs.next()) {
				SequenceRegion sr = createSequenceRegion(rs); // populates basic fields and coord-system
//				fetchAttributes(sr, con); // populates attributes
				result.add(sr);
			}

		} catch (Exception e) {
			throw new AdaptorException("Rethrow + stacktrace" + sql, e);
		} finally {
			close(con);
		}

		return (SequenceRegion[])result.toArray(new SequenceRegion[result.size()]);

	}


}
