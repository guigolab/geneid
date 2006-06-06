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

import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.plugin.standard.MySQLDriver;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.util.LruCache;

/**
 * Base class for non-location adaptors.
 */
public abstract class BasePersistentAdaptor {

	protected LruCache cache;

	public BasePersistentAdaptor(VariationDriver vdriver) {
		this.vdriver = vdriver;
	}

	public BasePersistentAdaptor(VariationDriver vdriver, int cacheSize) {
		this.vdriver = vdriver;
		cache = new LruCache(cacheSize);
	}

	protected List fetchListByQuery(String sql) throws AdaptorException {
		List r = new ArrayList();

		Connection conn = null;
		try {
			conn = vdriver.getConnection();
			ResultSet rs = conn.createStatement().executeQuery(sql);
			if (rs.next()) {
				Object o = null;
				while ((o = createObject(rs)) != null)
					r.add(o);
			}
		} catch (SQLException e) {
			throw new AdaptorException("Failed to fetch items with query: "
					+ sql, e);
		} finally {
			MySQLDriver.close(conn);
		}

		return r;

	}

	protected Object fetchByQuery(String sql) throws AdaptorException {
		List r = fetchListByQuery(sql);
		return r.size() > 0 ? r.get(0) : null;
	}

	/**
	 * @see org.ensembl.driver.Adaptor#clearCache()
	 */
	public void clearCache() throws AdaptorException {
		if (cache != null)
			cache.clear();
	}

	/**
	 * @see org.ensembl.driver.Adaptor#closeAllConnections()
	 */
	public void closeAllConnections() throws AdaptorException {
		vdriver.closeAllConnections();
	}

	protected VariationDriver vdriver;

	/**
	 * Returns an object of the correct type for the implementing adaptor or
	 * null if there are no more in the result set.
	 * 
	 * @param rs
	 * @return
	 * @throws SQLException
	 * @throws AdaptorException
	 */
	protected abstract Object createObject(ResultSet rs) throws SQLException,
			AdaptorException;

}
