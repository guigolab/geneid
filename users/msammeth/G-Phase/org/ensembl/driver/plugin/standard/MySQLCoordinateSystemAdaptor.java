/*
 * Copyright (C) 2003 EBI, GRL
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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Logger;

import javax.sql.DataSource;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.CoordinateSystemMapping;
import org.ensembl.datamodel.InvalidLocationException;
import org.ensembl.datamodel.Location;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.CoordinateSystemAdaptor;

/**
 * Implementation of CoordinateSystemAdaptor interface.
 *  
 */
// TODO - Implement store/delete
public class MySQLCoordinateSystemAdaptor
  extends BaseAdaptor
  implements CoordinateSystemAdaptor {

  private class MaxLengthMap extends HashMap {

	private static final long serialVersionUID = 1L;

	private void put(CoordinateSystem cs, String tableName, int maxLength) {
      put(maxLengthKey(cs, tableName), new Integer(maxLength));
    }

    private int get(CoordinateSystem cs, String tableName) {
      Object o = get(maxLengthKey(cs, tableName));
      if (o == null)
        return -1;
      else
        return ((Integer) o).intValue();
    }

    private String maxLengthKey(CoordinateSystem cs, String tableName) {
      return cs.getName() + "." + cs.getVersion() + "." + tableName;
    }

  }

  private static final Logger logger =
    Logger.getLogger(MySQLCoordinateSystemAdaptor.class.getName());

  private List additionalDataSources = new ArrayList();
  private ArrayList mappings = null;
  private HashMap coordSystemCache = null;
  private HashMap featureTableCache = null;
  private MaxLengthMap maxLengthCache = null;

  public MySQLCoordinateSystemAdaptor(MySQLDriver driver)
    throws AdaptorException {

    super(driver);
  }

  public CoordinateSystem fetch(long internalID) throws AdaptorException {

    if (coordSystemCache == null) {
      buildCaches();
    }
    return (CoordinateSystem) coordSystemCache.get(new Long(internalID));

  }

  public CoordinateSystem fetch(String name, String version)
    throws AdaptorException {

    if (coordSystemCache == null) {
      buildCaches();
    }

    Iterator it = coordSystemCache.keySet().iterator();
    while (it.hasNext()) {
      Long id = (Long) it.next();
      CoordinateSystem cs = (CoordinateSystem) coordSystemCache.get(id);
      if (cs.getName().equals(name)
        && (cs.getVersion() == null || cs.getVersion().length() == 0)) {
        return cs;
      } else if (
        cs.getName().equals(name)
          && (cs.getVersion() != null || cs.getVersion().equals(version))) {
        return cs;
      }
    }

    return null;
  }

  public CoordinateSystem[] fetchAll() throws AdaptorException {

    if (coordSystemCache == null) {
      buildCaches();
    }
    ArrayList result = new ArrayList(coordSystemCache.values());
    // you have to sort the result by rank
    Collections.sort(result, new Comparator() {
      public int compare(Object a, Object b) {
        int rankA = ((CoordinateSystem) a).getRank();
        int rankB = ((CoordinateSystem) b).getRank();
        if (rankA == rankB) {
          return 0;
        }
        if (rankA < rankB) {
          return -1;
        } else {
          return 1;
        }
      }
    });

    return (CoordinateSystem[]) result.toArray(
      new CoordinateSystem[result.size()]);
  }

  public CoordinateSystem fetchSequenceLevel() throws AdaptorException {

    if (coordSystemCache == null) {
      buildCaches();
    }
    Iterator it = coordSystemCache.keySet().iterator();
    while (it.hasNext()) {
      Long id = (Long) it.next();
      CoordinateSystem cs = (CoordinateSystem) coordSystemCache.get(id);
      if (cs.isSequenceLevel()) {
        return cs;
      }
    }
    return null;
  }

  public long store(CoordinateSystem cs) throws AdaptorException {

    return -1;
  }

  /**
  * Check if a string containing one or more comma-separated attributes
  * contains a particular attribute.
  * 
  * @param attribStr
  *            The string to check.
  * @param attrib
  *            The attribute to look for,.
  * @return True if attrib appears in attribStr
  */
  private boolean hasAttrib(String attribStr, String attrib) {

    String[] attribs = attribStr.split(",");
    for (int i = 0; i < attribs.length; i++) {
      if (attribs[i].equalsIgnoreCase(attrib)) {
        return true;
      }
    }

    return false;
  }

  private CoordinateSystem createCoordSystemFromResultSetRow(ResultSet rs)
    throws SQLException {

    CoordinateSystem result = new CoordinateSystem();
    result.setInternalID(rs.getLong("coord_system_id"));
    result.setName(rs.getString("name"));
    result.setVersion(rs.getString("version"));
    String attribStr = rs.getString("attrib");
    result.setDefault(hasAttrib(attribStr, "default_version"));
    result.setSequenceLevel(hasAttrib(attribStr, "sequence_level"));
    result.setRank(Integer.parseInt(rs.getString("rank")));

    return result;

  }

  public String getType() throws AdaptorException {

    return TYPE;
  }

  public CoordinateSystem[] getMappingPath(
    CoordinateSystem cs1,
    CoordinateSystem cs2)
    throws AdaptorException {

    CoordinateSystem[] result = new CoordinateSystem[2];
    // load and cache mapping info from meta table if required
    // mappings holds a list of CoordinateSystemMapping objects
    if (mappings == null) {

      mappings = new ArrayList();
      Connection con = null;
      String sql =
        "SELECT meta_value FROM meta WHERE meta_key='assembly.mapping'";
      try {

        con = getConnection();
        ResultSet rs = con.createStatement().executeQuery(sql);
        while (rs.next()) {

          String row = rs.getString("meta_value");
          if (row != null) { // build CoordinateSystemMapping
            // object from database
            String[] mapping = row.split("\\|");
            CoordinateSystem[] mappingPath =
              new CoordinateSystem[mapping.length];
            int pathIndex = 0;
            for (int i = 0; i < mapping.length; i++) {
              String[] parts = mapping[i].split(":");
              String name = parts[0];
              // always present
              String version = parts.length > 1 ? parts[1] : "";
              // optional
              mappingPath[pathIndex++] = fetch(name, version);
            }
            CoordinateSystemMapping csm =
              new CoordinateSystemMapping(mappingPath);
            mappings.add(csm);
          }
        }

      } catch (Exception e) {
        logger.warning(sql);
        throw new AdaptorException("Rethrow + stacktrace", e);
      } finally {
        close(con);
      }

    } // --------------------------------
    // get the required mapping path from the list of mappings
    for (Iterator it = mappings.iterator(); it.hasNext();) {

      CoordinateSystemMapping mapping = (CoordinateSystemMapping) it.next();
      if (mapping.getFirst().equals(cs1)
        && mapping.getLast().equals(cs2)
        || mapping.getFirst().equals(cs2)
        && mapping.getLast().equals(cs1)) {

        return mapping.getPath();
      }
    }

    return null;
  }

  /**
   * Store all the entries in the database in a HashMap of CoordinateSystem
   * objects. The HashMap is keyed on internal ID, so retrieval by internal
   * ID will be very fast. Other retrieval methods will be slower but this
   * should not be a problem as the coord_system table is always likely to be
   * small. Also stores the feature table name / coordinate system mappings
   * in featureTableCache.
   */
  private void buildCaches() throws AdaptorException {

    Connection con = null;
    coordSystemCache = new HashMap();
    featureTableCache = new HashMap();
    maxLengthCache = new MaxLengthMap();
    String sql = "";

    try {
      con = getConnection();
      sql = "SELECT * FROM coord_system";
      ResultSet rs = con.createStatement().executeQuery(sql);
      while (rs.next()) {
        CoordinateSystem cs = createCoordSystemFromResultSetRow(rs);
        coordSystemCache.put(new Long(cs.getInternalID()), cs);
      }
    } catch (Exception e) {
      throw new AdaptorException("Rethrow + stacktrace" + sql, e);
    } finally {
      close(con);
    }


    // load meta_coord data from ALL available datasources.
    List dss = getAllDataSources();
    for (int i = 0, n = dss.size(); i < n; i++) {
      DataSource ds = (DataSource) dss.get(i);

      try {

        // feature tables
        // note key for featureTableCache = tablename, 
        // value = ArrayList of CoordinateSystem objects
        // deliberately do "SELECT *" because the code needs to
        // handle older schema (2 columns) and newer schema 
        // (3 columns)
        sql = "SELECT * FROM meta_coord";
        con = ds.getConnection();
        ResultSet rs = con.createStatement().executeQuery(sql);
        int nCols = 0;
        while (rs.next()) {
          String tableName = rs.getString("table_name").toLowerCase();
          // note case insensitive
          CoordinateSystem cs = fetch(rs.getLong("coord_system_id"));
          if (cs == null) {
            throw new AdaptorException(
              "meta_coord table refers to non-existant co-ordinate system with ID "
                + rs.getLong("coord_system_id"));
          }
          ArrayList csList = (ArrayList) featureTableCache.get(tableName);
          if (csList == null) {
            csList = new ArrayList();
          }
          csList.add(cs);
          featureTableCache.put(tableName, csList);

          if (nCols == 0)
            nCols = rs.getMetaData().getColumnCount();
          if (nCols > 2)
            maxLengthCache.put(cs, tableName, rs.getInt("max_length"));

        }

      } catch (Exception e) {
        throw new AdaptorException("Rethrow + stacktrace" + sql, e);
      } finally {
        close(con);
      }
    }
    //				Iterator it1 = featureTableCache.keySet().iterator();
    //				while (it1.hasNext()) {
    //					String tableName = (String)it1.next();
    //					System.out.print(tableName + ": ");
    //					ArrayList csList = (ArrayList)featureTableCache.get(tableName);
    //					Iterator it2 = csList.iterator();
    //					while (it2.hasNext()) {
    //					CoordinateSystem cs = (CoordinateSystem)it2.next();
    //						System.out.println("\t" + cs.getInternalID() + " " + cs.getName() + " " + cs.getVersion());
    //					}
    //				}

  }

  public CoordinateSystem[] fetchAllByFeatureTable(String featureTableName)
    throws AdaptorException {

    if (featureTableCache == null) {
      buildCaches();
    }
    String tableName = featureTableName.toLowerCase();
    ArrayList csList = (ArrayList) featureTableCache.get(tableName);
    if (csList == null) {
      throw new AdaptorException(
        "Cannot get coordinate system for " + tableName);
    }

    return (CoordinateSystem[]) csList.toArray(
      new CoordinateSystem[csList.size()]);

  }

  /**
   * This function supports old style Locations with Maps. Maps are essentially CoordinateSystems
   * so we retrieve one by mapname
   * @param mapName a Map identifier 
   * @return the equivalent CoordinateSystem object
   * @throws AdaptorException
   */

  public CoordinateSystem fetchByMap(String mapName) throws AdaptorException {
    return fetch(mapName, "");
  }

  public CoordinateSystem fetchComplete(CoordinateSystem skeletonCS)
    throws AdaptorException {

    return fetch(skeletonCS.getName(), skeletonCS.getVersion());

  }

  public List fetchTopLevelLocations() throws AdaptorException {

    List locs = new ArrayList();

    String sql =
      "SELECT coord_system_id, sr.name, length "
        + "FROM seq_region sr,seq_region_attrib sra, attrib_type at "
        + "WHERE sr.seq_region_id=sra.seq_region_id and sra.attrib_type_id=at.attrib_type_id and code=\"toplevel\"";

    Connection conn = null;
    CoordinateSystemAdaptor csAdaptor = driver.getCoordinateSystemAdaptor();
    conn = getConnection();
    ResultSet rs = executeQuery(conn, sql);
    try {
      while (rs.next()) {
        locs.add(
          new Location(
            csAdaptor.fetch(rs.getLong(1)),
            rs.getString(2),
            1,
            rs.getInt(3)));
      }
    } catch (InvalidLocationException e) {
      throw new AdaptorException("Problem constructing top level location", e);
    } catch (SQLException e) {
      throw new AdaptorException("Problem constructing top level location", e);
    }

    return locs;
  }

  /**
   * Returns the maximum length of the feature 
   * in the coordinate system or -1 if no max length is 
   * found.
   * @param cs coordinate system feature is in.
   * @param tableName table name feature is stored in.
   * @return max length, or -1 if no max length is available.
   */
  public int fetchMaxLength(CoordinateSystem cs, String tableName) {
    return maxLengthCache.get(cs, tableName);
  }

  /**
   * Clears this adaptors caches.
   */
  public void clearCache() {
    coordSystemCache = null;
    featureTableCache = null;
    mappings = null;
    maxLengthCache = null;
  }

  public void addDataSource(DataSource dataSource) {
    additionalDataSources.add(dataSource);
  }

  public boolean removeDataSource(DataSource dataSource) {
    return additionalDataSources.remove(dataSource);
  }

  public List getAllDataSources() throws AdaptorException {
    List dss = new ArrayList(additionalDataSources);
    dss.add(0, getDataSource());
    return dss;
  }
}
