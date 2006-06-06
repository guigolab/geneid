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

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.logging.Logger;

import javax.sql.DataSource;

import org.ensembl.datamodel.CloneFragmentLocation;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Persistent;
import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.util.ConnectionPoolDataSource;
import org.ensembl.util.LruCache;
import org.ensembl.util.Warnings;

abstract public class BaseAdaptor implements Adaptor {
  private static final Logger logger =
    Logger.getLogger(BaseAdaptor.class.getName());

  final String NULL = "NULL";

  protected LruCache cache;
  protected MySQLDriver driver;
  protected DataSource dataSource = null;

  protected BaseAdaptor(MySQLDriver driver) {
    this.driver = driver;
  }

  /**
   * Creates driver member with specified driver attached and an
   * internalID->Persistent cache. See Cache for meaning of cache parameters.
   * @see LruCache org.ensembl.util.LruCache
   */
  protected BaseAdaptor(MySQLDriver driver, int cachemaxCapacity) {
    this(driver);
    if (cachemaxCapacity > 0)
      cache = new LruCache(cachemaxCapacity);
  }

  /**
   * Closes any open datasource connections if the underlying datasource is a 
   * ConnectionPoolDataSource (which it is by default). Prints a warning if the 
   * underlying connection is not a ConnectionPoolDataSource.
   * 
   * Each adaptor can have a separate datasource so it is necessary to ensure all connections
   * are closed.
   * 
   * @throws AdaptorException if a problem occurs closing the connections.
   */
  public void closeAllConnections() throws AdaptorException {
    ConnectionPoolDataSource.closeAllConnections(dataSource);
  }

  /**
   * Clears cache if cache is used. Does nothing otherwise.
   *
   */
  public void clearCache() {
    if (cache != null)
      cache.clear();
  }

  public final org.ensembl.driver.Driver getDriver() {
    return driver;
  }

  public final boolean supportsMap(String map) {
    return false;
  }

  /**
   * Configure the adaptor.
   */
  void configure() throws ConfigurationException, AdaptorException {
    // do nothing
  }

  /**
   * Convenience method for closing a connection. Also sets
   * conn.setAutoCommit(true) before closing, this is useful if conn is returned
   * to a connection pool. Prints logger.warninging if exception occurs.
   */
  public static void close(Connection conn) {

    MySQLDriver.close(conn);
  }

  /**
   * Convenience method for closing a connection on a ResultSet. Also sets
   * conn.setAutoCommit(true) before closing, this is useful if conn is returned
   * to a connection pool. Prints logger.warninging if exception occurs.
   */
  public static void close(ResultSet rs) {

    MySQLDriver.close(rs);
  }

  public static void rollback(Connection conn) {
    try {
      if (conn != null) {
        conn.rollback();
      }
    } catch (SQLException e) {
      logger.warning("Failed to rollback transaction. " + e.getMessage());
    }
  }

  /**
   * Convenience method which wraps conn.createStatement().executeUpdate(sql).
   * Sql is included in exception.message if exception thrown.
   * @return number of rows affected.
   */
  public static int executeUpdate(Connection conn, String sql)
    throws AdaptorException {
    try {
      return conn.createStatement().executeUpdate(sql);
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql, e);
    }
  }

  /**
   * Convenience method which wraps ps.executeUpdate() in try/catch and
   * includes SQL in thrown exception.
   * @return number of rows affected.
   */
  public static int executeUpdate(PreparedStatement ps, String sql)
    throws AdaptorException {
    try {
      return ps.executeUpdate();
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql, e);
    }
  }

  /**
   * Convenience method which wraps conn.createStatement().executeQuery(sql).
   * Sql is included in exception.message if exception thrown.
   * @return ResultSet generated by executing the query.
   */
  public static ResultSet executeQuery(Connection conn, String sql)
    throws AdaptorException {
    try {
      return conn.createStatement().executeQuery(sql);
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql, e);
    }
  }

  /**
   * Convenience method which wraps ps.executeQuery().
   * Sql is included in exception.message if exception thrown.
   * @return ResultSet generated by executing the query.
   */
  public static ResultSet executeQuery(PreparedStatement ps, String sql)
    throws AdaptorException {
    try {
      return ps.executeQuery();
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql, e);
    }
  }

  /**
   * Executes the sql which should include an autoincrement element.
   * @return internalID internalID auto generated by database.
   */
  public static long executeAutoInsert(Connection conn, String sql)
    throws AdaptorException {

    long internalID = 0;
    String sql2 = sql;

    try {

      int nRows = conn.createStatement().executeUpdate(sql2);
      if (nRows != 1)
        throw new AdaptorException("Failed to insert to database: " + sql2);

      sql2 = "select last_insert_id()";
      ResultSet rs = conn.createStatement().executeQuery(sql2);
      rs.next();
      internalID = rs.getLong(1);
      if (internalID <= 0)
        throw new AdaptorException(
          "Auto increment generated an unacceptable internalID: "
            + internalID
            + " : "
            + sql2);
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql2, e);
    }

    return internalID;
  }

  /**
   * Executes the PreparedStatement which should include an autoincrement element.
   * @return internalID auto generated by database (last autoincremented value).
   */
  public static long executeAutoInsert(PreparedStatement ps, String sql)
    throws AdaptorException {

    long internalID = 0;
    String sql2 = sql;

    try {

      int nRows = ps.executeUpdate();
      if (nRows != 1)
        throw new AdaptorException("Failed to insert to database: " + sql2);

      sql2 = "select last_insert_id()";
      ResultSet rs = ps.getConnection().createStatement().executeQuery(sql2);
      rs.next();
      internalID = rs.getLong(1);
      if (internalID <= 0)
        throw new AdaptorException(
          "Auto increment generated an unacceptable internalID: "
            + internalID
            + " : "
            + sql2);
    } catch (SQLException e) {
      throw new AdaptorException("Failed to execute sql:" + sql2, e);
    }

    return internalID;
  }

  /**
   * @return location as a clonefragmentLocation, if location is an
   * AssemblyLocation it is converted to a CloneFragmentLocation without gaps.
   */
  protected CloneFragmentLocation getAsCloneFragmentLocation(Location loc)
    throws AdaptorException {

    Warnings.deprecated(
      "CloneFrangmentLocations no longer supported - returning null");

    return null;
  }

  /**
   * Convenience clears buffer.
   */
  final void clear(StringBuffer buf) {
    buf.delete(0, Integer.MAX_VALUE);
  }

  /**
   * Add item to cache.
   * @param persistent object to be stored in cache
   */
  protected void addToCache(Persistent persistent) {
    if (cache != null && persistent != null)
      cache.put(persistent, new Long(persistent.getInternalID()));
  }

  /**
   * @return object from cache with specified internalID if present,
   * otherwise null
   */
  Persistent fetchFromCache(long internalID) {

    if (cache == null || internalID < 1)
      return null;
    else {
      Object o = cache.get(new Long(internalID));
      if (o == null)
        return null;
      else
        return (Persistent) o;
    }

  }

  Persistent deleteFromCache(long internalID) {
    if (cache == null || internalID < 1)
      return null;
    else {
      return (Persistent) cache.removeValueByKey(new Long(internalID));
    }
  }

  /**
   * Returns a connection to either the generic (driver) 
   * database or a database specified for this adaptor.
   * The specific database is only used if one is specified 
   * in the driver configuration.
   * @return connection that is relevant to this adaptor.
   */
  public Connection getConnection() throws AdaptorException {

    Connection conn = null;
    try {
      conn = getDataSource().getConnection();
    } catch (SQLException e) {
      throw new AdaptorException("", e);
    }

    return conn;
  }

  /**
   * Returns the datasource for this adaptor.
   * 
   * @return datasource that is relevant to this adaptor.
   * @throws AdaptorException
   */
  public DataSource getDataSource() throws AdaptorException {

    if (dataSource == null)
      dataSource = driver.getDataSource(getType());

    return dataSource;
  }

}
