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
package org.ensembl.driver.plugin.standard;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import javax.sql.DataSource;

import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.EnsemblDriver;
import org.ensembl.util.ConnectionPoolDataSource;
import org.ensembl.util.StringUtil;

/**
 * Base class for all Ensembl Drivers providing database connection 
 * support and adaptor management.
 * 
 * Derived classes should implement loadAdaptors().
 * 
 * 
 * @see #loadAdaptors()
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public abstract class BaseDriver implements EnsemblDriver {

  private final static Logger logger =
    Logger.getLogger(BaseDriver.class.getName());

  private HashMap adaptors = new HashMap();

  protected Configuration configuration;

  private String[] databaseNames = null;

  /**
   * Creates a drive initialised with the specified configuration.
   * @param configuration this should be a Properties instance.
   * @throws AdaptorException
   * @see #initialise(Object) for configuration parameters.
   */
  public BaseDriver(Properties configuration) throws AdaptorException {
    initialise(configuration);
  }

  /**
   * Constructs a driver pointing at the specified database. Assumes no
   * password and port = 3306.
   * 
   * @param host
   *            computer hosting mysqld database
   * @param database
   *            database name
   * @param user
   *            user name
   * @param databaseIsPrefix true is database is to be used as a prefix
   * or false if it is to be used unmodified as a database name.
   */
  public BaseDriver(
    String host,
    String database,
    String user,
    boolean databaseIsPrefix)
    throws AdaptorException {
    this(host, database, user, null, null, databaseIsPrefix);
  }

  /**
   * Constructs a driver pointing at the specified database. Assumes port =
   * 3306.
   * 
   * @param host
   *            computer hosting mysqld database
   * @param database
   *            database name
   * @param user
   *            user name
   * @param password
   *            password
   * @param databaseIsPrefix true is database is to be used as a prefix
   * or false if it is to be used unmodified as a database name.
   */
  public BaseDriver(
    String host,
    String database,
    String user,
    String password,
    boolean databaseIsPrefix)
    throws AdaptorException {
    this(host, database, user, null, null, databaseIsPrefix);
  }

  /**
   * Creates an unitialised driver with no adaptor.
   * 
   * Call initialise(Properties) to initialise this driver.
   * 
   * @see #initialise(Properties)
   */
  public BaseDriver() {
  }

  /**
   * Constructs a driver pointing at the specified database.
   * 
   * @param host
   *            computer hosting mysqld database
   * @param database
   *            database name
   * @param user
   *            user name
   * @param password
   *            password
   * @param port
   *            port on host computer that mysqld is running on
   */
  public BaseDriver(
    String host,
    String database,
    String user,
    String password,
    String port,
    boolean databaseIsPrefix)
    throws AdaptorException {
    Properties p = new Properties();

    if (host == null)
      throw new AdaptorException("host can not be null");
    else
      p.setProperty("host", host);

    if (user == null)
      throw new AdaptorException("user can not be null");
    else
      p.setProperty("user", user);

    if (password != null && !"".equals(password))
      p.setProperty("password", password);

    if (port != null)
      p.setProperty("port", port);

    if (database != null) {
      if (databaseIsPrefix)
        p.setProperty("database_prefix", database);
      else
        p.setProperty("database", database);
    }

    try {
      initialise(p);
    } catch (ConfigurationException e) {
      throw new AdaptorException("Failed to configure driver with.", e);
    }

  }

  /**
     * Attempts to connects to the database if not already connected.
     *
     * Useful for checking if a driver has been correctly configured.
     * 
     * @return true if driver can connect to it's persistent store, otherwise
     *         false.
     */
  public synchronized boolean isConnected() {

    boolean connected = false;
    Connection conn = null;
    try {
      conn = getConnection();
    } catch (AdaptorException e) {
      logger.warning(e.getMessage());
    } finally {
      if (conn != null) {
        connected = true;
        close(conn);
      }
    }
    return connected;
  }

  public Adaptor addAdaptor(Adaptor adaptor) throws AdaptorException {
    BaseAdaptor previous =
      (BaseAdaptor) adaptors.put(adaptor.getType(), adaptor);
    if (previous != null) {
      previous.driver = null;
      logger.fine(
        "Adaptor "
          + previous.getClass().getName()
          + " replaced by "
          + previous.getClass().getName());
    } else {
      logger.fine("Added " + adaptor.getClass().getName() + " to Driver");
    }

    return adaptor;
  }

  public synchronized void removeAdaptor(Adaptor adaptor)
    throws AdaptorException {
    removeAdaptor(adaptor.getType());
  }

  public synchronized void removeAdaptor(String type) {
    BaseAdaptor adaptor = (BaseAdaptor) adaptors.remove(type);
    if (adaptor != null) {
      adaptor.driver = null;
      logger.fine("Removed " + adaptor.getClass().getName() + " from Driver");
    }
  }

  public synchronized void removeAllAdaptors() throws AdaptorException {
    for (Iterator iter = adaptors.values().iterator(); iter.hasNext();)
      removeAdaptor((Adaptor) iter.next());

  }

  public synchronized Connection getConnection() throws AdaptorException {

    DataSource ds = getDataSource();
    if (ds == null)
      return null;

    Connection conn = null;
    try {

      logger.fine("Getting connection ... ");
      conn = ds.getConnection();
      logger.fine("Got connection.");
    } catch (Exception e) {
      throw new AdaptorException(
        "Failed to initialise database connection pool : ",
        e);
    }

    return conn;
  }

  /**
     * Closes any open datasource connections if the underlying datasource is a 
     * ConnectionPoolDataSource (which it is by default). Prints a warning if the 
     * underlying connection is not a ConnectionPoolDataSource.
     * 
     * @throws AdaptorException if a problem occurs closing the connections.
     */
  public synchronized void closeAllConnections() throws AdaptorException {
    for (Iterator iter = dataSource.values().iterator(); iter.hasNext();) {
      DataSource ds = (DataSource) iter.next();
      ConnectionPoolDataSource.closeAllConnections(ds);
      dataSource.remove(ds);
    }

    // each adaptor might have it's own datasource so we need
    // to make sure those are closed as well.
    Adaptor[] adaptors = getAdaptors();
    for (int i = 0; i < adaptors.length; i++)
      adaptors[i].closeAllConnections();

  }

  /**
     * Clears all caches.
     * 
     * @throws AdaptorException if a problem occurs closing the connections.
     */
  public synchronized void clearAllCaches() throws AdaptorException {

    // each adaptor might have it's own datasource so we need
    // to make sure those are closed as well.
    Adaptor[] adaptors = getAdaptors();
    for (int i = 0; i < adaptors.length; i++)
      adaptors[i].clearCache();

    databaseNames = null;
  }

  /**
     * Convenience method for closing a connection. Also sets
     * conn.setAutoCommit(true) before closing, this is useful if conn is
     * returned to a connection pool. Prints logger.warninging if exception
     * occurs.
     */
  public static void close(Connection conn) {

    try {
      if (conn != null) {
        if (conn.isClosed()) {
          logger.warning("connection already closed, can't close again!");
        } else {
          conn.setAutoCommit(true);
          conn.close();
        }
      }
    } catch (SQLException e) {
      logger.warning(
        "Failed to set connection.autocommit=true OR close connection. "
          + e.getMessage());
    }
  }

  /**
       * Convenience method for closing a connection on a result set. 
       * It can not setAutoCommit(true) before closing. Prints logger.warninging if exception
       * occurs.
       */
  public static void close(ResultSet rs) {

    try {
      if (rs != null) {
        rs.close();
      }
    } catch (SQLException e) {
      logger.warning(
        "Failed to set connection.autocommit=true OR close connection. "
          + e.getMessage());
    }
  }

  /**
     * Creates a connection pool datasource pointing at the database specified in properties.
     * If there are adaptor specific settings then these override the default (driver)
     * settings.
     * @param type adaptor.getType() is used to identify adaptor specific datasources.
     * @return new connection pool.
     */
  protected DataSource createDataSource(String type) throws AdaptorException {

    StringBuffer err = new StringBuffer();

    int poolSize = -1;
    String poolSizeStr =
      configuration.getProperty(type, "connection_pool_size");
    if (poolSizeStr == null)
      err.append("connection_pool_size, ");
    else
      poolSize = Integer.parseInt(poolSizeStr);

    String jdbc_driver = configuration.getProperty(type, "jdbc_driver");
    if (jdbc_driver == null)
      err.append("jdbc_driver, ");

    String connection_string =
      configuration.getProperty(type, "connection_string");
    if (connection_string == null)
      err.append("connection_string, ");

    String databaseName = configuration.getProperty(type, "database");
    String databaseNamePrefix =
      configuration.getProperty(type, "database_prefix");
    if (databaseName == null)
      err.append("database, ");

    String user = configuration.getProperty(type, "user");
    if (user == null)
      err.append("user, ");

    String password = configuration.getProperty(type, "password");

    if (err.length() > 0)
      throw new AdaptorException(
        "Failed to initialise driver for "
          + type
          + " because missing parameters: "
          + err.toString());

    try {

      String fullConnectionString =
        connection_string + "/" + databaseName + "?autoReconnect=true";

      ConnectionPoolDataSource ds =
        new ConnectionPoolDataSource(
          jdbc_driver,
          fullConnectionString,
          user,
          password,
          poolSize);

      // attempt to resolve database name if necessary
      if (databaseName.equals("") && databaseNamePrefix != null) {
        databaseName = resolveDatabaseName(ds, type, databaseNamePrefix);
        if (type != null && type.length() > 0)
          configuration.putProperty(type, "database", databaseName);
        else
          configuration.put("database", databaseName);

        fullConnectionString =
          connection_string + "/" + databaseName + "?autoReconnect=true";
        ds =
          new ConnectionPoolDataSource(
            jdbc_driver,
            fullConnectionString,
            user,
            password,
            poolSize);
      }

      return ds;

    } catch (ClassNotFoundException e) {
      throw new AdaptorException("", e);
    } catch (SQLException e) {
      throw new AdaptorException("", e);
    }

  }

  public synchronized Adaptor getAdaptor(String type) throws AdaptorException {
    return (Adaptor) adaptors.get(type);
  }

  public synchronized Adaptor[] getAdaptors() throws AdaptorException {
    int len = adaptors.size();
    Adaptor[] adaptorArray = new Adaptor[len];
    adaptors.values().toArray(adaptorArray);
    return adaptorArray;
  }

  /**
   * Initialises the driver using the parameters in the
   * configuration.
   * 
   * If "database_prefix" is specified then this is used to 
   * specify the latest version of the database beginning
   * with the prefix. In this case "database" should not be specified.
   * For example "database_prefix"="homo_sapiens_core" will
   * be resolved to "database"="homo_sapiens_core_HIGHEST_VERSION".
   * 
   * All caches are cleared, connections closed and
   * adaptors removed before calling 
   * <code>processConfiguration(Properties)</code> and then
   * <code>loadAdaptors()</code>. 
   * 
   * Derived classes can implement there own 
   * processConfiguration(Configuration) if they want to modify
   * the configuration before 
   * <code>loadAdaptors()</code> is called.
   * 
   * @param config configuration parameters.
   * @throws ConfigurationException
   * @throws AdaptorException
   */
  public synchronized void initialise(Properties config)
    throws ConfigurationException, AdaptorException {

    clearAllCaches();
    closeAllConnections();
    removeAllAdaptors();

    // copy to avoid modifying the parameter passed in
    configuration = new Configuration(config);

    // amend the configuration, should be overridden by 
    // derived classes for special behaviour
    processConfiguration(configuration);

    // load the adaptors, should be overridden by implementing
    // classes/
    loadAdaptors();

  }

  /**
   * Modifies properties if needed and validates it.
   * 
   * Prints warnings or
   * throws an exception if the configuration is invalid.
   * 
   * Derived classes should override this method if they want 
   * different behaviour.
   * 
   * <ul>Default properties inserted if the key is missing:
   * <li>"jdbc_driver" : "org.gjt.mm.mysql.Driver"
   * <li>"connection_pool_size" : "10"
   * <li>"connection_string" : "jdbc:mysql://" + host + <port>:
   * <li>"password" : ""
   * <li>"database" : database
   * </ul>
   *  
   * @param properties object to be modified if necessary.
   * @throws ConfigurationException
   */
  protected void processConfiguration(Properties properties)
    throws ConfigurationException {

    String databaseNamePrefix = properties.getProperty("database_prefix");
    String databaseName = properties.getProperty("database");
    String user = properties.getProperty("user");
    String host = properties.getProperty("host");
    String port = properties.getProperty("port");
    String connStr = properties.getProperty("connection_string");

    if (user == null)
      throw new ConfigurationException("user is not set.");
    if (host == null && connStr == null)
      throw new ConfigurationException("host is not set.");

    if (databaseName == null) {
      // default to empty string, this allows the user to connect to a db via
      // a driver and retrieve database names without having to know any
      // legal database names before hand.
      databaseName = "";
      properties.put("database", databaseName);
    }

    if (!databaseName.equals("") && databaseNamePrefix != null)
      logger.warning(
        "Ignoring \"database\" because \"database_prefix\" is set. "
          + "Remove one of these parameters from the configuration.");

    if (!properties.containsKey("jdbc_driver"))
      properties.put("jdbc_driver", "org.gjt.mm.mysql.Driver");
    if (!properties.containsKey("connection_pool_size"))
      properties.put("connection_pool_size", "10");
    if (connStr == null)
      properties.put(
        "connection_string",
        "jdbc:mysql://" + host + ((port != null) ? (":" + port) : ""));
    if (!properties.containsKey("password")) {
      properties.put("password", "");
    }
    // name can be used by a DriverManager instance to identify
    // this driver instance.
    if (!properties.containsKey("name"))
      if (databaseNamePrefix != null)
        properties.put("name", databaseNamePrefix);
      else
        properties.put("name", databaseName);

    if (logger.isLoggable(Level.FINE)) {
      logger.fine("Configuring driver : " + properties);
      java.util.Enumeration names = this.configuration.elements();
      while (names.hasMoreElements()) {
        String key = (String) names.nextElement();
        String value = properties.getProperty(key);
        logger.fine("name: " + key + " value: " + value);
      }
    }

  }

  /**
   * This method is to be overridden by implementing classes.
   * 
   * It is a call back method called by initialise(Object).
   * 
   * @see #initialise(Properties)
   */
  protected abstract void loadAdaptors()
    throws AdaptorException, ConfigurationException;

  public synchronized Properties getConfiguration() {
    return configuration;
  }

  /**
     * Lists databases available on this server (if is a database server).
     * 
     * @return names of zero or more databases on the same server as this
     *         driver instance.
     */
  public synchronized String[] fetchDatabaseNames() throws AdaptorException {
    return fetchDatabaseNames(getDataSource());
  }

  private synchronized String[] fetchDatabaseNames(DataSource ds)
    throws AdaptorException {

    if (databaseNames == null) {
      List dbNames = null;
      Connection conn = null;
      String sql = "show databases;";
      try {
        conn = ds.getConnection();
        ResultSet rs = conn.createStatement().executeQuery(sql);
        if (rs.next()) {
          dbNames = new ArrayList();
          do {
            dbNames.add(rs.getString(1));
          } while (rs.next());
        }
      } catch (SQLException e) {
        throw new AdaptorException("Failed to read database names" + sql, e);
      } finally {
        BaseAdaptor.close(conn);
      }

      databaseNames = new String[] {
      };
      if (dbNames != null)
        databaseNames = (String[]) dbNames.toArray(databaseNames);

    }

    return databaseNames;
  }

  /**
     * Returns the default datasource.
     * @return default datasource.
     */
  protected DataSource getDataSource() throws AdaptorException {
    return getDataSource("");
  }

  /**
     * Returns the datasource for the adaptor type. This could be the default datasource
     * or one specific to the adaptor.
     * @param adaptorType adaptor type.
     * @return the datasource for the adaptor.
     * @throws AdaptorException
     */
  protected DataSource getDataSource(String adaptorType)
    throws AdaptorException {
    if (adaptorType == null)
      throw new NullPointerException("adaptorType can not be null");
    DataSource ds = (DataSource) dataSource.get(adaptorType);
    if (ds == null) {
      ds = createDataSource(adaptorType);
      dataSource.put(adaptorType, ds);
    }
    return ds;
  }

  /**
   * Returns the name of the "latest" database that begins with prefix.
   * 
   * The latest database name has the pattern prefix + "_"  + DIGITS + chars.
   * If more than one database name match this pattern then the one with
   * the highest DIGITS is chosen. e.g. converts "homo_sapiens_core" into 
   * "homo_sapiens_core_24_34e" if 24_34e is the latest 
   * version of the database.
   * 
   * @param ds datasource to retrieve data from.
   * @return real database name begining with prefix.
   * @throws AdaptorException if problem occurs resolving database
   * name or if no database corresponding to the prefix exists
   * in the datasource.
   */
  private String resolveDatabaseName(
    DataSource ds,
    String adaptorType,
    String databaseNamePrefix)
    throws AdaptorException {

    String[] all = fetchDatabaseNames(ds);

    // filter irrelevant database names:
    // we are only interested in databases that match
    // prefix + "_" + digits + other chars.
  	List dbs = new ArrayList();    
    Pattern p = Pattern.compile(databaseNamePrefix+"_\\d+.*");
    for (int i = 0; i < all.length; i++) 
    	if (p.matcher(all[i]).find())
    		dbs.add(all[i]);
	
    Collections.sort(dbs);
    String realDatabaseName = (String) (dbs.size()==0 ? null : dbs.get(dbs.size()-1));
    
    if (realDatabaseName == null)
      throw new AdaptorException(
        "Failed to match database name to any  on databaser server. database="
          + databaseNamePrefix
          + "\tavailable on server="
          + StringUtil.toString(dbs));

    return realDatabaseName;

  }

  public String toString() {
    String s = null;
    try {
      Connection conn = getConnection();
      s = conn.getMetaData().getURL();
      conn.close();
    } catch (Exception e) {
      s = "ERROR: " + e.getMessage();
    }
    return s;
  }

  protected Map dataSource = new HashMap();

}
