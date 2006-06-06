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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.ensembl.datamodel.AssemblyException;
import org.ensembl.datamodel.AssemblyMapper;
import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AssemblyExceptionAdaptor;
import org.ensembl.driver.AssemblyMapperAdaptor;
import org.ensembl.driver.CoordinateSystemAdaptor;
import org.ensembl.driver.LocationConverter;
import org.ensembl.util.LruCache;
import org.ensembl.util.mapper.Coordinate;
import org.ensembl.util.mapper.Mapper;

/**
 * Location converter which satisfies convert() calls by querying the
 * database. Fast to initialise because no 'data loading' stage. Individual
 * calls to convert however require separate database calls.
 *
 * <p>Possible optimisations based on cacing assembly:
 * <ul>
 *   <li>Include 'gaps' in cache rather than dynamically generate them.
 *   <li>Hold 2 versions of assembly in memory for fast access:
 *       assemblyMap->cloneFragmentMap and cloneFragmentMap->assemblyMap.
 *   <li>To speed loading consider writing/readinf serialised version of cache.
 * </ul>
 * @version $Revision: 1.1 $
 * */
public class MySQLLocationConverter
  extends BaseAdaptor
  implements LocationConverter {
  private static final Logger logger =
    Logger.getLogger(MySQLLocationConverter.class.getName());
  public static final int SEQ_REGION_CACHE_SIZE = 200000;

  private LruCache seqRegionCache;

  private class RegionCacheElement {
    long id;
    CoordinateSystem cs;
    String seqRegionName;
    int regionLength;
  };

  public MySQLLocationConverter(MySQLDriver driver) {
    super(driver);
    seqRegionCache = new LruCache(SEQ_REGION_CACHE_SIZE);
  }

  public void cacheSeqRegion(
    String regionName,
    CoordinateSystem cs,
    long regionId,
    int regionLength) {

    RegionCacheElement rce = new RegionCacheElement();
    rce.cs = cs;
    rce.id = regionId;
    rce.seqRegionName = regionName;
    rce.regionLength = regionLength;

    seqRegionCache.put(
      rce,
      (regionName + ":" + cs.getInternalID()),
      new Long(regionId));
  }

  public long[] namesToIds(String[] names, CoordinateSystem cs)
    throws AdaptorException {
    long[] result = new long[names.length];
    try {
      for (int i = 0; i < names.length; i++) {
        result[i] = nameToId(names[i], cs);
      }
    } catch (Exception e) {
      // rethrow
      throw new AdaptorException("rethrow", e);
    }
    return result;
  }

  /**
   * Extract the internalIds from the given Location
   * Ignore gap bits. 
   */
  public long[] locationToIds(Location loc) throws AdaptorException {
    Set result = new HashSet();
    while (loc != null) {
      if (loc.isGap() && logger.isLoggable(Level.FINE)) {
        logger.fine(
          "MySQLLocatinConvertor.locationToIds ignoring gap for "
            + loc.toString());
      }
      if (!loc.isGap()) {
        result.add(
          new Long(
            nameToId(loc.getSeqRegionName(), loc.getCoordinateSystem())));
      }
      loc = loc.next();
    }

    long[] ids = new long[result.size()];
    int i = 0;
    Iterator iter = result.iterator();
    while (iter.hasNext())
      ids[i++] = ((Long) iter.next()).longValue();

    return ids;
  }

  public int getLengthByLocation(Location loc) throws AdaptorException {
    RegionCacheElement rce =
      (RegionCacheElement) seqRegionCache.get(
        loc.getSeqRegionName()
          + ":"
          + loc.getCoordinateSystem().getInternalID());
    if (rce != null) {
      return rce.regionLength;
    }

    String sql =
      " SELECT "
        + "    seq_region_id, length "
        + " FROM  "
        + "    seq_region "
        + " WHERE  "
        + "    name = ? AND "
        + "    coord_system_id = ? ";
    logger.fine(sql);
    Connection conn = null;
    try {
      conn = getConnection();
      PreparedStatement ps = conn.prepareStatement(sql);
      ps.setString(1, loc.getSeqRegionName());
      ps.setLong(2, loc.getCoordinateSystem().getInternalID());
      ResultSet result = ps.executeQuery();
      if (!result.next()) {
        return -1;
      }
      long seqRegionId = result.getLong(1);
      int regionLength = result.getInt(2);
      cacheSeqRegion(
        loc.getSeqRegionName(),
        loc.getCoordinateSystem(),
        seqRegionId,
        regionLength);
      return regionLength;
    } catch (Exception e) {
      throw new AdaptorException("rethrow", e);
    } finally {
      close(conn);
    }
  }

  public long nameToId(String seqRegionName, CoordinateSystem cs)
    throws AdaptorException {

    // SQL requires a fully-specified co-ordinate system INCLUDING internal ID
    // often locations will be passed in that just have the coordinate system name
    // set; if this is the case we need to get a "full" co-ordinate system
    if (cs.getInternalID() == 0) {
      cs = driver.getCoordinateSystemAdaptor().fetchComplete(cs);
    }

    RegionCacheElement rce =
      (RegionCacheElement) seqRegionCache.get(
        seqRegionName + ":" + cs.getInternalID());
    if (rce != null) {
      return rce.id;
    }

    String sql =
      " SELECT "
        + "    seq_region_id, length "
        + " FROM  "
        + "    seq_region "
        + " WHERE  "
        + "    name = ? AND "
        + "    coord_system_id = ? ";
    logger.fine(sql);
    Connection conn = null;
    try {
      conn = getConnection();
      PreparedStatement ps = conn.prepareStatement(sql);
      ps.setString(1, seqRegionName);
      ps.setLong(2, cs.getInternalID());
      ResultSet result = ps.executeQuery();
      if (!result.next()) {
        return -1;
      }
      long seqRegionId = result.getLong(1);
      int regionLength = result.getInt(2);
      cacheSeqRegion(seqRegionName, cs, seqRegionId, regionLength);
      return seqRegionId;
    } catch (Exception e) {
      throw new AdaptorException("rethrow", e);
    } finally {
      close(conn);
    }
  }

  public Location idToLocation(long id) throws AdaptorException {
    return idToLocation(id, 0, 0, 0);
  }

  public Location idToLocation(long id, int start, int end, int strand)
    throws AdaptorException {
    RegionCacheElement rce =
      (RegionCacheElement) seqRegionCache.get(new Long(id));
    if (rce == null) {
      String sql =
        " SELECT "
          + "    coord_system_id, name, length "
          + " FROM  "
          + "    seq_region "
          + " WHERE  "
          + "    seq_region_id = ?";

      Connection conn = null;
      try {
        conn = getConnection();
        PreparedStatement ps = conn.prepareStatement(sql);
        ps.setLong(1, id);
        ResultSet result = ps.executeQuery();
        if (!result.next()) {
          return null;
        }
        long coordSystemId = result.getLong(1);
        String regionName = result.getString(2);
        int regionLength = result.getInt(3);
        CoordinateSystem cs =
          driver.getCoordinateSystemAdaptor().fetch(coordSystemId);

        cacheSeqRegion(regionName, cs, id, regionLength);

        // default values if necessary
        if (start == 0)
          start = 1;
        if (end == 0)
          end = regionLength;

        Location loc = new Location(cs, regionName, start, end, strand, false);

        return loc;
      } catch (Exception e) {
        throw new AdaptorException("rethrow", e);
      } finally {
        close(conn);
      }
    } else {
      //		default values if necessary
      if (start == 0)
        start = 1;
      if (end == 0)
        end = rce.regionLength;

      Location loc =
        new Location(rce.cs, rce.seqRegionName, start, end, strand, false);
      return loc;
    }
  }

  /**
   * Note: converting an incomplete sourceLocation (missing either or both start and end) will
   * result in a new location with a start and end based on where the sourceLocation overlaps
   * the target coordinate system. Consequently, converting an incomplete location into 
   * another coordinate system and back again will not produce exactly the same location as the original. 
   */
  public Location convert(
    Location sourceLocation,
    CoordinateSystem targetCS,
    boolean includeGaps,
    boolean allList,
    boolean setSequenceRegion)
    throws AdaptorException {

    if (targetCS == null) {
      throw new IllegalArgumentException("Need target coordinate system to be set");
    }

    CoordinateSystem sourceCS = sourceLocation.getCoordinateSystem();
    if (sourceCS == null) {
      throw new IllegalArgumentException("Location needs valid coordinate system");
    }

    if (!sourceLocation.isSeqRegionNameSet())
      return new Location(targetCS);

    // we need to ensure the CoordinateSystems are complete before we can use
    // them.
    CoordinateSystem cs =
      driver.getCoordinateSystemAdaptor().fetchComplete(
        sourceLocation.getCoordinateSystem());
    if (cs == null)
      throw new AdaptorException(
        "No CoordinateSystem in database corresponding to: "
          + sourceLocation.getCoordinateSystem());
    sourceLocation = sourceLocation.copy();
    sourceLocation.setCoordinateSystem(cs);
    targetCS = driver.getCoordinateSystemAdaptor().fetchComplete(targetCS);

    // we need to know the start and end before we can do the conversion
    if (!sourceLocation.isStartSet() || !sourceLocation.isEndSet())
      sourceLocation = fetchComplete(sourceLocation);

    // can't convert the location if it doesn not correspond to entries 
    // in the database
    if (sourceLocation == null)
      return null;

    if (sourceLocation.getCoordinateSystem().equals(targetCS)) {

      // not 100 percent correct yet
      // should remove gaps when includeGpas is false
      // maybe should make a copy 
      return sourceLocation;

    }

    AssemblyMapperAdaptor ama = driver.getAssemblyMapperAdaptor();
    AssemblyMapper am =
      ama.fetchByCoordSystems(sourceLocation.getCoordinateSystem(), targetCS);
    Coordinate coords[];
    Location cLoc = sourceLocation;
    Location resultLocation = null;
    Location loc;

    do {
      coords = am.map(cLoc);
      for (int i = 0; i < coords.length; i++) {
        loc = null;
        if (!coords[i].isGap()) {
          loc =
            new Location(
              targetCS,
              coords[i].id,
              coords[i].start,
              coords[i].end,
              coords[i].strand);

        } else if (includeGaps) {
          loc =
            new Location(
              targetCS,
              coords[i].id,
              coords[i].start,
              coords[i].end,
              coords[i].strand,
              true);

        }

        if (loc != null) {

          if (setSequenceRegion) {
            loc.setSequenceRegion(
              driver.getSequenceRegionAdaptor().fetch(
                loc.getSeqRegionName(),
                loc.getCoordinateSystem()));
          }

          if (resultLocation == null) {
            resultLocation = loc;
          } else {
            resultLocation.append(loc);
          }
        }
      }
      cLoc = cLoc.next();
      loc = null;
    } while (allList && cLoc != null);

    // load the sequence region attribute of resultLocation if required
    if (setSequenceRegion) {

      SequenceRegion sr =
        driver.getSequenceRegionAdaptor().fetch(
          resultLocation.getSeqRegionName(),
          resultLocation.getCoordinateSystem());
      resultLocation.setSequenceRegion(sr);

    }

    return resultLocation;
  }

  /**
   * Delegates conversion to appropriate sister method.
   *
   * @return location translated into map coordinates, or null if no
   * equivalent exists.
   *
   * @throws AdaptorException if the location type is unsupported.  
   */
  public Location convert(
    Location sourceLocation,
    String targetMap,
    boolean includeGaps,
    boolean allList,
    boolean setSequenceRegion)
    throws AdaptorException {
    Location newLoc = null;
    CoordinateSystemAdaptor csa = driver.getCoordinateSystemAdaptor();
    CoordinateSystem tcs = csa.fetchByMap(targetMap);
    return convert(
      sourceLocation,
      tcs,
      includeGaps,
      allList,
      setSequenceRegion);
  }

  /**
   * Converts location and returned location does not have gaps.
   * @return location translated into map coordinates, or null if no
   * equivalent exists.
   */
  public Location convert(Location sourceLocation, String targetMap)
    throws AdaptorException {
    return convert(sourceLocation, targetMap, false, true, false);
  }

  /**
   * Standard conversion with gaps and all Location is converted.
   * SequenceRegion attrib in returned location is NOT filled.
   * @param loc
   * @param cs
   * @return location translated into map coordinates, or null if no
   * equivalent exists.
   * @throws AdaptorException
   */

  public Location convert(Location loc, CoordinateSystem cs)
    throws AdaptorException {
    return convert(loc, cs, true, true, false);
  }

  /**
   * Perform an in place edit where necessary. Does nothing if the location is complete or 
   * can not be made complete.
   */
  public Location fetchComplete(Location location) throws AdaptorException {

    // Note: we use recursive implementation rather than a loop so that we can rollback
    // changes (by not making them) if a problem occurs trying to edit either one of
    // the values in this location node or one of the later ones in a location list.

    // collect the values for cs, srt, start and end
    CoordinateSystem cs =
      driver.getCoordinateSystemAdaptor().fetchComplete(
        location.getCoordinateSystem());

    // if coordinate system is not in database we can't return a complete
    // location.
    if (cs == null)
      return null;

    String srName = location.getSeqRegionName();
    
    // If only the coordinate system is specified then we need to 
    // create a location list that includes all of the sequence regions 
    // in the coordinate system.
    if (srName==null) {
      
      SequenceRegion[] rs = driver.getSequenceRegionAdaptor().fetchAllByCoordinateSystem(cs);
      Location l = location;
      for(int i=0; i<rs.length; ++i) {
        l.setSeqRegionName(rs[i].getName());
        l.setStart(1);
        l.setEnd((int)rs[i].getLength());
        l.setStrand(0);
        l.setNext(new Location(cs));
        l = l.next();
      }
      
      return location;
    }
    
    int start = -1;
    int end = -1;
    if (srName != null) {

      start = location.getStart();
      if (start < 1)
        start = 1;

      end = location.getEnd();
      if (end < 1) {
        SequenceRegion sr = driver.getSequenceRegionAdaptor().fetch(location);
        // if sequence region is not in db we can't construct the complete
        // location
        if (sr == null)
          return null;
        if (end < 1)
          end = (int) sr.getLength();
      }
    }

    // update next node if necessary, return null if problem occured updating it
    Location next = location.next();
    if (next != null)
      if (fetchComplete(next) == null)
        return null;

    // store cs, sr, start, end if we get this far
    location.setCoordinateSystem(cs);
    if (srName!=null) {
      location.setSeqRegionName(srName);
      location.setStart(start);
      location.setEnd(end);
    }

    return location;
  }

  public Location dereference(Location loc) throws AdaptorException {
    Location deref;
    AssemblyExceptionAdaptor aea = driver.getAssemblyExceptionAdaptor();
    AssemblyException[] aexs = aea.fetchLinks(loc);

    if (aexs.length == 0) {
      return loc;
    }

    Mapper mapper = new Mapper("link", "org");

    for (int i = 0; i < aexs.length; i++) {
      if (aexs[i].type.equals("PAR")) {
        mapper.addMapCoordinates(
          aexs[i].linkedLocation.getSeqRegionName(),
          aexs[i].linkedLocation.getStart(),
          aexs[i].linkedLocation.getEnd(),
          aexs[i].linkedLocation.getStrand(),
          aexs[i].originalLocation.getSeqRegionName(),
          aexs[i].originalLocation.getStart(),
          aexs[i].originalLocation.getEnd());
      } else if (aexs[i].type.equals("HAP")) {
        mapper.addMapCoordinates(
          aexs[i].linkedLocation.getSeqRegionName(),
          1,
          aexs[i].linkedLocation.getStart() - 1,
          aexs[i].linkedLocation.getStrand(),
          aexs[i].originalLocation.getSeqRegionName(),
          1,
          aexs[i].originalLocation.getStart() - 1);
        // now the bit after the HAP, need to find the end coordinates
        int orgLength = getLengthByLocation(aexs[i].originalLocation);
        int linkedLength = getLengthByLocation(aexs[i].linkedLocation);
        mapper.addMapCoordinates(
          aexs[i].linkedLocation.getSeqRegionName(),
          aexs[i].linkedLocation.getEnd() + 1,
          linkedLength,
          aexs[i].linkedLocation.getStrand(),
          aexs[i].originalLocation.getSeqRegionName(),
          aexs[i].originalLocation.getEnd() + 1,
          orgLength);
      } else {
        throw new AdaptorException("Unknown AssemblyException type. Know HAP and PAR currently");
      }
    }
    // now mapper can do the translation
    Coordinate[] coords =
      mapper.mapCoordinate(
        loc.getSeqRegionName(),
        loc.getStart(),
        loc.getEnd(),
        loc.getStrand(),
        "link");

    // now map over to original location, coordinate system stays the same
    Location targetLoc = null;
    for (int i = 0; i < coords.length; i++) {
      if (i == 0) {
        targetLoc =
          new Location(
            loc.getCoordinateSystem(),
            coords[i].id,
            coords[i].start,
            coords[i].end,
            coords[i].strand,
            coords[i].isGap());
      } else {
        targetLoc.append(
          new Location(
            loc.getCoordinateSystem(),
            coords[i].id,
            coords[i].start,
            coords[i].end,
            coords[i].strand,
            coords[i].isGap()));
      }
    }
    return targetLoc;
  }

  public String getType() {
    return TYPE;
  }

  /**
   * Assemblies in dataset.
   * @return list of zero or more assemblies.
   * @deprecated should use CoordinateSystemAdaptor
   */
  public String[] fetchAssemblyNames() throws AdaptorException {
    CoordinateSystemAdaptor csa = driver.getCoordinateSystemAdaptor();
    CoordinateSystem cs[] = csa.fetchAll();
    ArrayList resultList = new ArrayList();
    for (int i = 0; i < cs.length; i++) {
      if (cs[i].getName().equals("Chromosome")) {
        resultList.add(cs[i].getVersion());
      }
    }
    String[] assembliesArr = null;
    assembliesArr = (String[]) resultList.toArray(assembliesArr);

    return assembliesArr;
  }

  /**
   * 
   * @deprecated please use locationToIds( Location loc) after you converted to your target 
   * coordinate system
   * @throws AdaptorException
   */
  public long[] listSeqRegionIds(
    Location loc,
    CoordinateSystem targetCoordinateSystem)
    throws AdaptorException {
    return new long[0];
  }

  /**
   * @see org.ensembl.driver.LocationConverter#convertToTopLevel(org.ensembl.datamodel.Location)
   */
  public Location convertToTopLevel(Location location)
    throws AdaptorException {

    Location loc = null;

    CoordinateSystem[] css = driver.getCoordinateSystemAdaptor().fetchAll();
    for (int i = 0; loc == null && i < css.length; i++) {
      CoordinateSystem cs = css[i];
      loc = convert(location, cs);
    }

    return loc;
  }

} // MySQLLocationConverter
