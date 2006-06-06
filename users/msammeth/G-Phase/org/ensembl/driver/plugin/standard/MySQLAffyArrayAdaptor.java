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
import org.ensembl.datamodel.impl.AffyArrayImpl;
import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AffyArrayAdaptor;

/**
 * Implemementation of the AffyArrayAdaptor that works with standard ensembl
 * myssql databases.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class MySQLAffyArrayAdaptor extends MySQLBaseFeatureAdaptor implements
		Adaptor, AffyArrayAdaptor {

	/**
	 * @param driver
	 */
	public MySQLAffyArrayAdaptor(MySQLDriver driver) {
		super(driver, TYPE);
	}

	/**
	 * @see org.ensembl.driver.AffyArrayAdaptor#fetch(long)
	 */
	public AffyArray fetch(long internalID) throws AdaptorException {
		return (AffyArray) fetchByInternalID(internalID);
	}

	/**
	 * @see org.ensembl.driver.AffyArrayAdaptor#fetch()
	 */
	public List fetch() throws AdaptorException {
		return fetchByNonLocationConstraint("");
	}

	/**
	 * @see org.ensembl.driver.AffyArrayAdaptor#fetch(java.lang.String)
	 */
	public AffyArray fetch(String name) throws AdaptorException {
		List tmp = fetchByNonLocationConstraint(" name = \"" + name + "\"");
		return (AffyArray) (tmp.size() == 0 ? null : tmp.get(0));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#tables()
	 */
	public String[][] tables() {
		final String[][] tables = { { "affy_array", "aa" } };
		return tables;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#columns()
	 */
	public String[] columns() {
		final String[] columns = { "aa.affy_array_id", "aa.parent_array_id",
				"aa.probe_setsize", "aa.name" };
		return columns;
	}

	/**
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
	 */
	public Object createObject(ResultSet rs) throws AdaptorException {
		try {
			if (!rs.next())
				return null;

			return new AffyArrayImpl(driver, rs.getLong(1), rs.getString(4), rs
					.getInt(3));
		} catch (SQLException e) {
			throw new AdaptorException(
					"Failed to create AffyArray from database", e);
		}
	}

}
