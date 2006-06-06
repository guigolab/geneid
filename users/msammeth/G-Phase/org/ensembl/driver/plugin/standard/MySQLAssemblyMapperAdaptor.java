/*
 * Created on 09-Oct-2003
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package org.ensembl.driver.plugin.standard;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.ensembl.datamodel.AssemblyMapper;
import org.ensembl.datamodel.ChainedAssemblyMapper;
import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.SimpleAssemblyMapper;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AssemblyMapperAdaptor;
import org.ensembl.driver.CoordinateSystemAdaptor;
import org.ensembl.util.mapper.Coordinate;
import org.ensembl.util.mapper.Mapper;
import org.ensembl.util.mapper.Range;
import org.ensembl.util.mapper.RangeRegistry;

/**
 * @author arne
 *
 */
public class MySQLAssemblyMapperAdaptor
	extends BaseAdaptor
	implements AssemblyMapperAdaptor {

	private HashMap seqRegionCache;
	private HashMap mapperCache;

	private static final int CHUNKFACTOR = 20;
	private static final int MAXPAIRCOUNT = 3000;

	public static void main(String[] args) {
	}

	public MySQLAssemblyMapperAdaptor(MySQLDriver driver) {
		super(driver);
		mapperCache = new HashMap();
		seqRegionCache = new HashMap();
	}

	public String getType() {
		return TYPE;
	}

	public AssemblyMapper fetchByCoordSystems(
		CoordinateSystem cs1,
		CoordinateSystem cs2)
		throws AdaptorException {

		AssemblyMapper resultMapper;
		resultMapper =
			(AssemblyMapper) mapperCache.get(
				cs1.toString() + " " + cs2.toString());
		if (resultMapper != null) {
			return resultMapper;
		}
		CoordinateSystemAdaptor csa = driver.getCoordinateSystemAdaptor();
		CoordinateSystem[] path = csa.getMappingPath(cs1, cs2);
		if (path == null) {
			throw new AdaptorException("Can't map between coordinate systems " + cs1 + " and " + cs2);
		}
		if (path.length == 2) {
			// simple mapping	
			resultMapper = new SimpleAssemblyMapper(driver, path[0], path[1]);
		} else if (path.length == 3) {
			// chained mapping
			resultMapper =
				new ChainedAssemblyMapper(driver, path[0], path[1], path[2]);
		} else {
			// no mapping
			throw new AdaptorException("Cant map between" + cs1 + " and " + cs2);
		}

		mapperCache.put((cs1.toString() + " " + cs2.toString()), resultMapper);
		mapperCache.put((cs2.toString() + " " + cs1.toString()), resultMapper);
		return resultMapper;
	}

	public void registerAssembled(
		SimpleAssemblyMapper asmMapper,
		String seqRegionName,
		int start,
		int end)
		throws AdaptorException {

		int startChunk = start >> CHUNKFACTOR;
		int endChunk = end >> CHUNKFACTOR;
		LinkedList regions = new LinkedList();

		Integer beginChunkRegion, endChunkRegion;
		beginChunkRegion = null;

		//	find regions of continuous unregistered chunks
		for (int i = startChunk; i <= endChunk; i++) {
			if (asmMapper.haveRegisteredAssembled(seqRegionName, i)) {
				if (beginChunkRegion != null) {
					// this is the end of an unregistered region.
					endChunkRegion = new Integer((i << CHUNKFACTOR) - 1);
					regions.add(beginChunkRegion);
					regions.add(endChunkRegion);
					beginChunkRegion = null;
				}
			} else {
				if (beginChunkRegion == null) {
					beginChunkRegion = new Integer((i << CHUNKFACTOR) + 1);
				}
				asmMapper.registerAssembled(seqRegionName, i);
			}
		}
		// the last part may have been an unregistered region too
		if (beginChunkRegion != null) {
			endChunkRegion = new Integer(((endChunk + 1) << CHUNKFACTOR) - 1);
			regions.add(beginChunkRegion);
			regions.add(endChunkRegion);
		}

		// nothing new needs registering
		if (regions.size() == 0) {
			return;
		}

		//	keep the Mapper to a reasonable size
		if (asmMapper.getSize() > MAXPAIRCOUNT) {
			asmMapper.flush();
			regions.clear();
			regions.add(new Integer((startChunk << CHUNKFACTOR) + 1));
			regions.add(new Integer(endChunk << CHUNKFACTOR));
			for (int i = startChunk; i <= endChunk; i++) {
				asmMapper.registerAssembled(seqRegionName, i);
			}
		}

		long seqRegionId;
		try {
			seqRegionId =
				driver.getLocationConverter().nameToId(
					seqRegionName,
					asmMapper.getAssembledCoordinateSystem());
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		String sql =
			"SELECT "
				+ "		asm.cmp_start,"
				+ "		asm.cmp_end,"
				+ "		asm.cmp_seq_region_id,"
				+ "	  sr.name,"
        + "   sr.length,"
				+ "		asm.ori,"
				+ "   asm.asm_start,"
				+ "		asm.asm_end"
				+ "	FROM "
				+ "		 assembly asm, seq_region sr "
				+ "	WHERE "
				+ "		 asm.asm_seq_region_id = ? AND "
				+ "		 ? <= asm.asm_end AND "
				+ "		 ? >= asm.asm_start AND "
				+ "		 asm.cmp_seq_region_id = sr.seq_region_id AND "
				+ "    sr.coord_system_id = ? ";

		//		Retrieve the description of how the assembled region is made from
		//		component regions for each of the continuous blocks of unregistered,
		//		chunked regions

		Iterator i = regions.iterator();
		Connection conn = null;
		try {
			conn = getConnection();
			PreparedStatement ps = conn.prepareStatement(sql);
			while (i.hasNext()) {
				int beginRegion = ((Integer) i.next()).intValue();
				int endRegion = ((Integer) i.next()).intValue();

				ps.setLong(1, seqRegionId);
				ps.setInt(2, beginRegion);
				ps.setInt(3, endRegion);
				ps.setLong(4, asmMapper.getComponentCoordinateSystem().getInternalID());
				ResultSet result = ps.executeQuery();
				while (result.next()) {
					String componentRegionName = result.getString(4);
					if (asmMapper.haveRegisteredComponent(componentRegionName)) {
						continue;
					}
					asmMapper.registerComponent(componentRegionName);
					asmMapper.mapper.addMapCoordinates(
						componentRegionName,
						result.getInt(1),
						result.getInt(2),
						result.getInt(6),
						seqRegionName,
						result.getInt(7),
						result.getInt(8));
					driver.getLocationConverter().cacheSeqRegion(
						componentRegionName,
						asmMapper.getComponentCoordinateSystem(),
						result.getLong(3),
            result.getInt(5));
				}
			}
		} catch (Exception e) {
			// something went wrong ...
			e.printStackTrace();
		} finally {
			close(conn);
		}
	}

	public void registerComponent(
		SimpleAssemblyMapper asmMapper,
		String componentSeqRegionName)
		throws AdaptorException {

		long componentCoordSystemId =
			asmMapper.getComponentCoordinateSystem().getInternalID();
		long asmCoordSystemId =
			asmMapper.getAssembledCoordinateSystem().getInternalID();

		if (asmMapper.haveRegisteredComponent(componentSeqRegionName)) {
			return;
		}

		long seqRegionId =
			driver.getLocationConverter().nameToId(
				componentSeqRegionName,
				asmMapper.getComponentCoordinateSystem());
		String sql =
			"	SELECT "
				+ "    asm.asm_start, "
				+ "    asm.asm_end, "
				+ "    asm.asm_seq_region_id, "
        + "    sr.name, "
        + "    sr.length "
				+ "	FROM "
				+ "    assembly asm, seq_region sr "
				+ " WHERE "
				+ "    asm.cmp_seq_region_id = ? AND "
				+ "    asm.asm_seq_region_id = sr.seq_region_id AND "
				+ "    sr.coord_system_id = ? ";

		Connection conn = null;

		try {
			conn = getConnection();
			PreparedStatement ps = conn.prepareStatement(sql);

			ps.setLong(1, seqRegionId);
			ps.setLong(2, asmCoordSystemId);

			ResultSet result = ps.executeQuery();
			if (!result.next()) {
				asmMapper.registerComponent(componentSeqRegionName);
				return;
			}
			String asmSeqRegion = result.getString(4);
			int asmStart = result.getInt(1);
			int asmEnd = result.getInt(2);
			long asmRegionId = result.getLong(3);
      int regionLength = result.getInt( 5 );
      
			if (result.next()) {
				throw (
					new Exception("Cant handle 2 assembly areas for same component"));
			}
			driver.getLocationConverter().cacheSeqRegion(
				asmSeqRegion,
				asmMapper.getAssembledCoordinateSystem(),
				asmRegionId,
        regionLength );
			registerAssembled(asmMapper, asmSeqRegion, asmStart, asmEnd);
		} catch (Exception e) {
			// some sql exception, maybe we should let it through
			e.printStackTrace();
		} finally {
			close(conn);
		}
	}

	public void registerChained(
		ChainedAssemblyMapper cam,
		String startTag,
		String seqRegionName,
		List ranges)
		throws AdaptorException {

		String asmsql =
			"SELECT "
				+ "  asm.cmp_start, "
				+ "  asm.cmp_end, "
				+ "  asm.cmp_seq_region_id, "
				+ "  sr.name, "
				+ "  asm.ori, "
				+ "  asm.asm_start, "
				+ "  asm.asm_end, "
        + "  sr.length "
				+ "FROM "
				+ "  assembly asm, seq_region sr "
				+ "WHERE "
				+ "  asm.asm_seq_region_id = ? AND"
				+ "  ? <= asm.asm_end AND "
				+ "  ? >= asm.asm_start AND "
				+ "  asm.cmp_seq_region_id = sr.seq_region_id AND "
				+ "  sr.coord_system_id = ?";

		String cmpsql =
			"SELECT "
				+ "  asm.asm_start, "
				+ "  asm.asm_end, "
				+ "  asm.asm_seq_region_id, "
				+ "  sr.name, "
				+ "  asm.ori, "
				+ "  asm.cmp_start, "
				+ "  asm.cmp_end, "
        + "  sr.length "
				+ "FROM "
				+ "  assembly asm, seq_region sr "
				+ "WHERE "
				+ "  asm.cmp_seq_region_id = ? AND"
				+ "  ? <= asm.cmp_end AND "
				+ "  ? >= asm.cmp_start AND "
				+ "  asm.asm_seq_region_id = sr.seq_region_id AND "
				+ "  sr.coord_system_id = ?";

		Mapper startMiddleMapper, endMiddleMapper, combinedMapper;
		CoordinateSystem startCS, endCS, midCS;
		RangeRegistry startReg, endReg;
		String endTag;

		if (startTag.equals("first")) {
			startMiddleMapper = cam.getFirstMiddleMapper();
			startCS = cam.getCsFirst();
			startReg = cam.getFirstReg();
			endMiddleMapper = cam.getLastMiddleMapper();
			endCS = cam.getCsLast();
			endReg = cam.getLastReg();
			endTag = "last";
		} else if (startTag.equals("last")) {
			startMiddleMapper = cam.getLastMiddleMapper();
			startCS = cam.getCsLast();
			startReg = cam.getLastReg();
			endMiddleMapper = cam.getFirstMiddleMapper();
			endCS = cam.getCsFirst();
			endReg = cam.getFirstReg();
			endTag = "first";
		} else {
			throw new AdaptorException("Wrong coord system tag"); //up
		}

		midCS = cam.getCsMiddle();
		combinedMapper = cam.getFirstLastMapper();
		CoordinateSystemAdaptor csa = driver.getCoordinateSystemAdaptor();
		CoordinateSystem[] path = csa.getMappingPath(startCS, midCS);
		// path[0] is assembled, path[1] component

		if (path.length != 2) {
			throw new AdaptorException(
				"should be able to go direct from " + startCS + " to " + midCS);
		}
		//			
		// obtain the first half of the mappings and load them into the start mapper
		//
		String sql;
		if (path[0].equals(startCS)) {
			sql = asmsql;
		} else {
			sql = cmpsql;
		}

		long seqRegionId =
			driver.getLocationConverter().nameToId(seqRegionName, startCS);
		long midCSId = midCS.getInternalID();

		LinkedList midRanges = new LinkedList();
		LinkedList startRanges = new LinkedList();

		Iterator i = ranges.iterator();

		Connection conn = getConnection();
		try {
			PreparedStatement ps = conn.prepareStatement(sql);

			while (i.hasNext()) {
				Range range = (Range) i.next();

				ps.setLong(1, seqRegionId);
				ps.setInt(2, range.start);
				ps.setInt(3, range.end);
				ps.setLong(4, midCSId);

				ResultSet result = executeQuery(ps, sql);
				while (result.next()) {
					int midStart = result.getInt(1);
					int midEnd = result.getInt(2);
					long midRegionId = result.getLong(3);
					String midRegionName = result.getString(4);
					int ori = result.getInt(5);
					int startStart = result.getInt(6);
					int startEnd = result.getInt(7);

					if (startMiddleMapper
						.addMapCoordinates(
							seqRegionName,
							startStart,
							startEnd,
							ori,
							midRegionName,
							midStart,
							midEnd)) {

						driver.getLocationConverter().cacheSeqRegion(
							midRegionName,
							midCS,
							midRegionId,
              result.getInt( 8 ));
              
						midRanges.add(new String(midRegionName));
						midRanges.add(new Long(midRegionId));
						midRanges.add(new Integer(midStart));
						midRanges.add(new Integer(midEnd));

						startRanges.add(new Range(startStart, startEnd));
						// there is a chance that we got more back from the query than what
						// we already registered, we register that, too
						if (startStart < range.start || startEnd > range.end) {
							startReg.checkAndRegister(seqRegionName, startStart, startEnd);
						}
					}
				}
			}

			//
			//			
			// now the second half of the mapping
			// perform another query and load the mid <-> end mapper using the mid cs
			// ranges

			path = csa.getMappingPath(midCS, endCS);
			// path[0] is assembled, path[1] component

			if (path.length != 2) {
				throw new Exception(
					"should be able to go direct from " + midCS + " to " + endCS);
			}

			if (path[0].equals(midCS)) {
				sql = asmsql;
			} else {
				sql = cmpsql;
			}

			long endCSId = endCS.getInternalID();

			i = midRanges.iterator();
			ps = conn.prepareStatement(sql);

			while (i.hasNext()) {
				String midRegionName = (String) i.next();
				long midRegionId = ((Long) i.next()).longValue();
				int midRegionStart = ((Integer) i.next()).intValue();
				int midRegionEnd = ((Integer) i.next()).intValue();

				ps.setLong(1, midRegionId);
				ps.setInt(2, midRegionStart);
				ps.setInt(3, midRegionEnd);
				ps.setLong(4, endCSId);

				ResultSet result = ps.executeQuery();
				while (result.next()) {
					int endStart = result.getInt(1);
					int endEnd = result.getInt(2);
					long endSeqRegionId = result.getLong(3);
					String endSeqRegionName = result.getString(4);
					int ori = result.getInt(5);
					int midStart = result.getInt(6);
					int midEnd = result.getInt(7);

					if (endMiddleMapper
						.addMapCoordinates(
							endSeqRegionName,
							endStart,
							endEnd,
							ori,
							midRegionName,
							midStart,
							midEnd)) {
						driver.getLocationConverter().cacheSeqRegion(
							endSeqRegionName,
							endCS,
							endSeqRegionId,
              result.getInt(8));
						endReg.checkAndRegister(endSeqRegionName, endStart, endEnd);
					}
				}
			}

			//		
			// Now that both halves are loaded
			// Do stepwise mapping using both of the loaded mappers to load
			// the final start <-> end mapper
			//		

			i = startRanges.iterator();
			while (i.hasNext()) {
				Range range = (Range) i.next();
				int sum = 0;
				Coordinate[] initialCoords =
					startMiddleMapper.mapCoordinate(
						seqRegionName,
						range.start,
						range.end,
						1,
						startTag);

				for (int coordIdx = 0; coordIdx < initialCoords.length; coordIdx++) {
					Coordinate coord = initialCoords[coordIdx];
					if (coord.isGap()) {
						sum += coord.length();
						continue;
					}

					Coordinate[] finalCoords =
						endMiddleMapper.mapCoordinate(
							coord.id,
							coord.start,
							coord.end,
							coord.strand,
							"middle");
					for (int finalCoordIdx = 0;
						finalCoordIdx < finalCoords.length;
						finalCoordIdx++) {
						Coordinate fcoord = finalCoords[finalCoordIdx];
						if (!fcoord.isGap()) {
							int totalStart = sum + range.start;
							int totalEnd = fcoord.length() + totalStart - 1;
							if (startTag.equals("first")) {
								combinedMapper.addMapCoordinates(
									seqRegionName,
									totalStart,
									totalEnd,
									fcoord.strand,
									fcoord.id,
									fcoord.start,
									fcoord.end);
							} else {
								combinedMapper.addMapCoordinates(
									fcoord.id,
									fcoord.start,
									fcoord.end,
									fcoord.strand,
									seqRegionName,
									totalStart,
									totalEnd);
							}
						}
						sum += fcoord.length();
					}
				}
			}
		} catch (SQLException e) {
			throw new AdaptorException("rethrow", e);
		} catch (Exception e) {
			throw new AdaptorException("rethrow", e);
		} finally {
			close(conn);
		}
	}
}
