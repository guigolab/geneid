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
package org.ensembl.driver;

import java.sql.Connection;
import java.util.Properties;


/**
 * An ensembl driver provides access to an ensembl database
 * throw it's adaptors.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public interface EnsemblDriver {
  
  
  /**
   * Initialises the driver and it's adaptor with the 
   * configuration. 
   * @param configuration parameters specifying the database
   * and any additional datasource specific settings. 
   * @throws ConfigurationException
   * @throws AdaptorException
   */
  void initialise(Properties configuration) throws ConfigurationException, AdaptorException;

  /**
   * 
   * @return configuration provided to initialise(Properties) plus
   * any locally made modifications.
   * @see #initialise(Properties)
   */
  Properties getConfiguration();


  /**
   * @return Adaptor of the specified type, or
   * null if no such driver available.
   */
  Adaptor getAdaptor(String type) throws AdaptorException;


  /**
   * @return array of zero or more Adaptors.
   */
  Adaptor[] getAdaptors() throws AdaptorException;


  /** @return whether the driver is connected to it's persistent store */
  boolean isConnected();

  
  /**
   * Lists databases available on this server.
   * 
   * This method does not require the database or 
   * database_prefix to specified in the configuration. 
   * This means
   * that you can connect to a database server and retrieve the 
   * database names without having to first know one of them.
   * 
   * @return zero or more database names available on the same server as this driver instance. 
   */
  String[] fetchDatabaseNames() throws AdaptorException;

  
  /**
   * Attempts to connects to the database if not already connected.
   *
   * Useful for checking if a driver has been correctly configured.
   * 
   * @return true if driver can connect to it's persistent store, otherwise
   *         false.
   */
  Adaptor addAdaptor(Adaptor adaptor) throws AdaptorException;
  
  void removeAdaptor(Adaptor adaptor) throws AdaptorException;
  
  void removeAdaptor(String type);
  
  void removeAllAdaptors() throws AdaptorException;

  Connection getConnection() throws AdaptorException;
  
  /**
   * Clears all caches.
   * 
   * @throws AdaptorException if a problem occurs closing the connections.
   */
  void clearAllCaches() throws AdaptorException;

  /**
   * Closes all connections opened by the driver.
   * @throws AdaptorException
   */
  void closeAllConnections() throws AdaptorException;

}