/*
 * Copyright (C) 2003 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package org.ensembl.driver.plugin.standard;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.ensembl.datamodel.Analysis;
import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Feature;
import org.ensembl.datamodel.Locatable;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Persistent;
import org.ensembl.datamodel.impl.LocatableImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AnalysisAdaptor;
import org.ensembl.driver.CoordinateSystemAdaptor;
import org.ensembl.driver.FeatureAdaptor;
import org.ensembl.driver.FeatureIterator;
import org.ensembl.driver.LocationConverter;
import org.ensembl.util.LocatableList;
import org.ensembl.util.LongList;
import org.ensembl.util.StringUtil;

/**
 * This class is the base class for all adaptors that work with "features" in
 * the EnsEMBL core database. A feature may be a SimpleFeature, RepeatFeature,
 * Gene etc. This class handles, in as generic a way as is possible, the
 * generation of SQL and the retrieval of the features from the database.
 * 
 * Subclasses should only need to define the tables and columns in which the
 * feature information is stored in the database, how to create a feature object
 * from that data, and of course any feature-specific functionality.
 * 
 * Note that feature storage is not yet supported.
 * 
 * HOW TO IMPLEMENT AN ADAPTOR
 * 
 * 1. Create a class that extends MySQLBaseFeatureAdaptor 2. Implement
 * constructors as necessary 3. Implement the tables() and columns() methods -
 * see Javadoc. 4. Override leftJoin(), finalClause(), finalWhereClause() if
 * required. 5. Implement createObject() 6. (for backwards API compatibility)
 * Implement <type>fetch(internalID) </type>7. Implement any feature-specific
 * functionality
 * 
 * See MySQLSimpleFeatureAdaptor as an example.
 * 
 */
public abstract class MySQLBaseFeatureAdaptor extends BaseAdaptor implements
		FeatureAdaptor {

	private static final Logger logger = Logger
			.getLogger(MySQLBaseFeatureAdaptor.class.getName());

	private String[] logicNames;

	private String featureType;

	private boolean available;

	protected CoordinateSystemAdaptor coordinateSystemAdaptor = null;

	protected MySQLLocationConverter locationConverter = null;

	protected AnalysisAdaptor analysisAdaptor = null;

	public final static int DEFAULT_FEATURE_BUFFER_SIZE = 100;

	private final static int QUERY_THRESHOLD = 3;

	private boolean loadChildren = true;

	// -----------------------------------------------------------------

	public MySQLBaseFeatureAdaptor(MySQLDriver driver, String type) {
		super(driver);
		this.featureType = type;
		available = true;
		if (driver != null)
			try {
				coordinateSystemAdaptor = driver.getCoordinateSystemAdaptor();
				locationConverter = (MySQLLocationConverter) driver
						.getLocationConverter();
				analysisAdaptor = driver.getAnalysisAdaptor();
			} catch (AdaptorException e) {
				e.printStackTrace();
			}
	}

	// -----------------------------------------------------------------

	public MySQLBaseFeatureAdaptor(MySQLDriver driver, String logicName,
			String type) {
		this(driver, type);
		logicNames = new String[1];
		logicNames[0] = logicName;
		available = logicName != null;
	}

	// -----------------------------------------------------------------

	public MySQLBaseFeatureAdaptor(MySQLDriver driver, String[] logicNames,
			String type) {
		this(driver, type);
		this.logicNames = logicNames;
		available = false;
		for (int i = 0; !available && i < logicNames.length; ++i)
			if (logicNames[i] != null)
				available = true;
	}

	// -----------------------------------------------------------------

	protected AnalysisAdaptor getAnalysisAdaptor() throws AdaptorException {
		if (analysisAdaptor == null)
			return (AnalysisAdaptor) getDriver().getAdaptor("analysis");
		else
			return analysisAdaptor;
	}

	// -----------------------------------------------------------------

	public String getType() {
		return featureType;
	}

	// -----------------------------------------------------------------

	/**
	 * Generates a condition string for an sql statement.
	 * <p>
	 * <i>For example: </i> <br>
	 * <code>SELECT * FROM table WHERE &lt;condition&gt;</code><br>
	 * where &lt;condition&gt; is being generated by this method.
	 * 
	 * @return String containing sql condition, or null.
	 */
	protected String getAnalysisIDCondition() throws AdaptorException {

		int n = 0;
		if (logicNames != null)
			n = logicNames.length;

		switch (n) {

		case 0:
			return null;

		case 1:
			return getAnalysisIDCondition(logicNames[0]);

		default:
			return getAnalysisIDCondition(logicNames);

		}

	}

	// -----------------------------------------------------------------

	/**
	 * Adaptors available if they have been successfully initialised, this does
	 * not necessarilly mean that there is any data fro them to retrieve.
	 * 
	 * @return true if adaptor is available for use, false otherwise.
	 */
	public boolean isAvailable() {
		return available;
	}

	// -----------------------------------------------------------------

	/**
	 * Creates an sql condition (to be used in a WHERE clause).
	 * 
	 * @return string containing an SQL condition if logicalName is resolved to
	 *         an analysis, otherwise null.
	 */
	protected String getAnalysisIDCondition(String logicalName)
			throws AdaptorException {
		Analysis analysis = driver.getAnalysisAdaptor().fetchByLogicalName(
				logicalName);
		if (analysis == null) {
			return null;
		} else {
			return getPrimaryTableSynonym() + ".analysis_id = "
					+ Long.toString(analysis.getInternalID());
		}
	}

	// -----------------------------------------------------------------

	protected String getAnalysisIDCondition(Analysis analysis)
			throws AdaptorException {
		return getAnalysisIDCondition(analysis.getLogicalName());
	}

	// -----------------------------------------------------------------

	protected String getAnalysisIDCondition(String[] logicalNames)
			throws AdaptorException {
		StringBuffer buf = new StringBuffer();
		buf.append(" " + getPrimaryTableSynonym() + ".analysis_id IN ( ");
		AnalysisAdaptor aa = driver.getAnalysisAdaptor();
		for (int i = 0; i < logicalNames.length; ++i) {
			Analysis a = aa.fetchByLogicalName(logicalNames[i]);
			buf.append(Long.toString(a.getInternalID()));
			if (i + 1 < logicalNames.length) {
				buf.append(", ");
			}
		}
		buf.append(" ) ");
		return buf.toString();
	}

	// -----------------------------------------------------------------

	protected String getAnalysisIDCondition(Analysis[] analyses)
			throws AdaptorException {

		String[] logicalNames = new String[analyses.length];
		for (int i = 0; i < logicalNames.length; ++i) {
			logicalNames[i] = analyses[i].getLogicalName();
		}

		return getAnalysisIDCondition(logicalNames);
	}

	// -----------------------------------------------------------------

	// -----------------------------------------------------------------

	public String[] getLogicNames() {
		return logicNames;
	}

	// -----------------------------------------------------------------
	/**
	 * Fetch all features of this type, over all coordinate systems.
	 */
	public List fetchAll() throws AdaptorException {

		return fetchByNonLocationConstraint("");

	}

	/**
	 * Default implementation ignores the loadChildren hint and just fetches all
	 * the Features. Sub classes that wish to implement the prefetch behaviour
	 * should over ride this method.
	 * 
	 * @param loadChildren
	 *            ignored in this default implementation.
	 * @return zero or more Features with the specified internalIDs.
	 * @throws AdaptorException
	 */
	public List fetchAll(boolean loadChildren) throws AdaptorException {

		return fetchAll();

	}

	/**
	 * 
	 * @return List of all the internal IDs for this type. Use with start, end
	 *         if you want batch
	 * @throws AdaptorException
	 */
	public long[] fetchInternalIDs() throws AdaptorException {
		return fetchInternalIDs(-1, -1);
	}

	/**
	 * Fetch internalIDs in batches.
	 * 
	 * @return zero or more internal ids for this type you should overwrite this
	 *         function and return empty List if your feature doesnt support
	 *         internal ids.
	 * @throws AdaptorException
	 * @param start
	 *            batch start if you want to use batches. First id is 1
	 * @param end
	 *            batch end. Includes this row as well. Counting starts at 1
	 */
	public long[] fetchInternalIDs(int start, int end) throws AdaptorException {
		LongList tmp = new LongList();
		String query = "";
		Connection conn = null;

		try {

			conn = getConnection();
			String tablename = getPrimaryTableName();
			String limit;
			if (start >= 0 && end >= 0) {
				limit = " LIMIT " + (start - 1) + "," + (end - start + 1);
			} else {
				limit = "";
			}
			query = "SELECT " + tablename + "_id" + " FROM " + tablename
					+ limit;

			ResultSet rs = executeQuery(conn, query);
			while (rs.next())
				tmp.add(rs.getLong(1));

		} catch (Exception e) {
			throw new AdaptorException("QUERY " + query + " failed ", e);
		} finally {
			close(conn);
		}
		return tmp.toArray();
	}

	public long[] fetchInternalIDs(final Location location)
			throws AdaptorException {

		LocatableList tmp = new LocatableList();
		Location fullLoc = driver.getLocationConverter()
				.fetchComplete(location);

		if (fullLoc == null)
			return tmp.toSortedInternalIDArray();

		Connection conn = null;
		StringBuffer query = new StringBuffer();

		try {

			conn = getConnection();
			String tablename = getPrimaryTableName();

			// get all the coordinate systems that this type of feature can be
			// stored in
			CoordinateSystem[] featureCoordSystems = coordinateSystemAdaptor
					.fetchAllByFeatureTable(getPrimaryTableName());

			for (int i = 0; i < featureCoordSystems.length; i++) {

				CoordinateSystem relCoordSys = featureCoordSystems[i];

				// TODO - dereference
				// vvvvvvvvvv This should be dereferencedLoc
				Location queryLoc = locationConverter.convert(fullLoc,
						relCoordSys);
				if (queryLoc == null) {
					throw new AdaptorException(
							"LocationConvertor.convert returned null");
				}

				for (Location loc = queryLoc; loc != null; loc = loc.next()) {

					// can't search on gaps because they don't correspond to a
					// specific region.
					if (loc.isGap())
						continue;

					query.delete(0, query.length());
					query.append("SELECT ").append(tablename).append("_id")
							.append(", seq_region_id").append(
									", seq_region_start").append(
									", seq_region_end").append(
									", seq_region_strand").append(" FROM ")
							.append(tablename);
					
					
					boolean where = true;
					boolean and = false;

					if (loc.getSeqRegionName() != null) {
						if (where){
							query.append(" WHERE ");
							where = false;
						}
						long seqRegionID = locationConverter.nameToId(loc
								.getSeqRegionName(), featureCoordSystems[i]);
						if (and)
							query.append(" AND ");
						query.append(tablename).append(".seq_region_id = ")
								.append(seqRegionID);
						and = true;
					}
					if (loc.getEnd() > 0) {
						if (where){
							query.append(" WHERE ");
							where = false;
						}
						if (and)
							query.append(" AND ");
						query.append(tablename).append(".seq_region_start <= ")
								.append(loc.getEnd());
						and = true;
					}

					if (loc.getStart() > 0) {
						if (where){
							query.append(" WHERE ");
							where = false;
						}
						if (and)
							query.append(" AND ");
						query.append(tablename).append(".seq_region_end >= ")
								.append(loc.getStart());
						and = true;
					}

					if (loc.getStrand() != 0) {
						if (where){
							query.append(" WHERE ");
							where = false;
						}
						if (and)
							query.append(" AND ");
						query.append(tablename).append(".seq_region_strand = ")
								.append(loc.getStrand());
						and = true;
					}

					ResultSet rs = executeQuery(conn, query.toString());

					while (rs.next()) {
						long id = rs.getLong(1);
						Location relLoc = locationConverter.idToLocation(rs
								.getInt(2));
						relLoc.setStart(rs.getInt(3));
						relLoc.setEnd(rs.getInt(4));
						relLoc.setStrand(rs.getInt(5));
						tmp.add(new LocatableImpl(id, relLoc));

					}

				}

			}

		} catch (Exception e) {
			throw new AdaptorException(
					"QUERY " + query.toString() + " failed ", e);
		} finally {
			close(conn);
		}

		return tmp.toSortedInternalIDArray();
	}

	
	/**
	 * Convenience method that returns an array of internal IDs 
	 * created by executing the SQL.
	 * @param sql sql statement with internal Ids in the first column of the
	 * select statement.
	 * @return iterator with zero or more elements.
	 * @throws AdaptorException
	 */
	protected long[] fetchInternalIDsBySQL(String sql) throws AdaptorException {
		
		LongList buf = new LongList();
		
		Connection conn = null;
		try {
			conn = getConnection();
			ResultSet rs = executeQuery(conn, sql);
			while(rs.next())
				buf.add(rs.getLong(1));
			
		} catch (SQLException e) {
			throw new AdaptorException("Failed to load internal ids from sql: " + sql,e);
		}
			
		return buf.toArray();
	}
	
	
	/**
	 * Convenience method that returns an iterator over the internal IDs 
	 * returned by executing the SQL.
	 * @param sql sql statement with internal Ids in the first column of the
	 * select statement.
	 * @return iterator with zero or more elements.
	 * @throws AdaptorException
	 */
	protected Iterator fetchIteratorBySQL(String sql) throws AdaptorException {
		return fetchIterator(fetchInternalIDsBySQL(sql));
	}
	
	// ---------------------------------------------------------------------

	public long fetchCount() throws AdaptorException {
		long result = 0;
		String query = null;
		Connection conn = null;

		try {
			conn = getConnection();
			String tablename = getPrimaryTableName();
			query = "SELECT count(distinct(" + tablename + "_id))" + " FROM "
					+ tablename;
			ResultSet rs = executeQuery(conn, query);
			if (rs.next())
				result = rs.getLong(1);

			System.out.println(query);
			System.out.println(result);

		} catch (Exception e) {
			throw new AdaptorException("QUERY " + query + " failed ", e);
		} finally {
			close(conn);
		}
		return result;
	}

	// ---------------------------------------------------------------------

	/**
	 * @param location
	 *            Location in which to search for Protein Similarity Features
	 * @param analysis
	 *            Filter to search
	 * 
	 * @return List of features matching the location and filtered by analysis.
	 */
	public List fetch(Location location, Analysis analysis)
			throws AdaptorException {

		return fetchAllByConstraint(location, getAnalysisIDCondition(analysis
				.getLogicalName()));
	}

	// -----------------------------------------------------------------

	/**
	 * @param location
	 *            Location in which to search for features
	 * @param analyses
	 *            Filter to search
	 * @return List of features matching the location and filtered by analysis.
	 */
	public List fetch(Location location, Analysis[] analyses)
			throws AdaptorException {

		return fetchAllByConstraint(location, getAnalysisIDCondition(analyses));

	}

	// -----------------------------------------------------------------

	public List fetch(Location location, String logicalName)
			throws AdaptorException {

		return fetchAllByConstraint(location,
				getAnalysisIDCondition(logicalName));

	}

	// -----------------------------------------------------------------
	/**
	 * @return A list of features inside the Location, and filtered by the
	 *         Analyses. An empty List is returned if non found.
	 */
	public List fetch(Location location, String[] logicalNames)
			throws AdaptorException {

		return fetchAllByConstraint(location,
				getAnalysisIDCondition(logicalNames));

	}

	// -----------------------------------------------------------------
	/**
	 * Constructs and executes sql statement which removes the entry with
	 * _internalID_ from the specified table.
	 */
	void delete(Connection conn, String tableName, long internalID)
			throws AdaptorException {
		executeUpdate(conn, "delete from " + tableName + " where " + tableName
				+ "_id=" + internalID);
	}

	// ---------------------------------------------------------------------
	//
	// NEW BASEFEATUREADAPTOR FUNCTIONALITY FOR NEW SCHEMA DATABASE
	//
	// ---------------------------------------------------------------------

	// ---------------------------------------------------------------------
	/**
	 * This method must return the table names and synonyms to be used in the
	 * SQL statement that is generated to query features from the database.
	 * Hence all the tables that are joined on/to must appear here.
	 * 
	 * @return An array with one entry for each table. Each array entry for each
	 *         table has two elements - table name and table synonym. The
	 *         <em>first</em> name/synonym pair must be the "primary" table,
	 *         i.e. the name of the feature table and its synonym.
	 */
	public abstract String[][] tables();

	// ---------------------------------------------------------------------
	/**
	 * This method must return the list of column names to use in the SQL
	 * statement.
	 * 
	 * @return List of column names (in the form table.column)
	 */
	public abstract String[] columns();

	// ---------------------------------------------------------------------
	/**
	 * This method should create a feature object from one or more rows of a
	 * ResultSet. Note that it is up to the implementation of this method to
	 * call rs.next(). Calling methods should also check that the returned
	 * object is not null before using it.
	 * 
	 * @return A feature type object, or null if the object cannot be created
	 *         (e.g. if the end of the ResultSet has been reached, or the
	 *         ResultSet is empty).
	 */
	public abstract Object createObject(ResultSet rs) throws AdaptorException;

	// ---------------------------------------------------------------------
	/**
	 * Subclasses that need to add a final clause to the WHERE clause of the
	 * generated SQL should override this method.
	 * 
	 * @return The final WHERE clause; note that a prefix of <code>AND</code>
	 *         is <em>not</em> required; the SQL generator will add this if
	 *         required.
	 */
	public String finalWhereClause() {
		return "";
	}

	// ---------------------------------------------------------------------
	/**
	 * Subclasses that need to add a final clause to the generated SQL (e.g.
	 * ORDER BY) should override this method.
	 * 
	 * @return The final clause.
	 */
	public String finalClause() {
		return "";
	}

	// ---------------------------------------------------------------------
	/**
	 * Subclasses that need to add a left join clause should override this
	 * method.
	 * 
	 * @return The left join string.
	 */
	public String[][] leftJoin() {

		String[][] lj = { {} };
		return lj;

	}

	// ---------------------------------------------------------------------
	/**
	 * Get all the features corresponding to a particular location from the
	 * database. TODO - dereference/rereference symlinked regions e.g.
	 * 
	 * @param requestLoc
	 *            the area to fetch features from. Note that currently this must
	 *            be a fully specified Location object, e.g. with start, end,
	 *            strand set, and have only one node.
	 * @param constraint
	 *            An optional extra constraint that is used in the SQL.
	 */
	protected List fetchAllByConstraint(Location requestLoc, String constraint)
			throws AdaptorException {

		requestLoc = driver.getLocationConverter().fetchComplete(requestLoc);

		// return empty list if the request location could not be resolved to a
		// complete=
		// location (e.g. because the coordinate system or sequence region were
		// bogus.)
		if (requestLoc == null)
			return Collections.EMPTY_LIST;

		List features = new ArrayList();

		// get all the coordinate systems that this type of feature can be
		// stored in
		CoordinateSystem[] featureCoordSystems = coordinateSystemAdaptor
				.fetchAllByFeatureTable(getPrimaryTableName());

		for (int i = 0; i < featureCoordSystems.length; i++) {

			// TODO - dereference
			// vvvvvvvvvv This should be dereferencedLoc
			Location queryLoc = locationConverter.convert(requestLoc,
					featureCoordSystems[i]);

			if (logger.isLoggable(Level.FINE)) {
				logger.fine("request cs = " + featureCoordSystems[i]);
				logger.fine("requestLoc = " + requestLoc);
				logger.fine("queryLoc = " + queryLoc);
			}

			if (queryLoc == null) {
				throw new AdaptorException(
						"LocationConvertor.convert returned null");
			}

			StringBuffer query = new StringBuffer(constraint);

			// query should be built differently depending on size: use IN list
			// for queries with many
			// locations and separate queries if there are a few.

			if (queryLoc.size() >= QUERY_THRESHOLD) {

				// build IN clause for IDs
				String idString = buildIDString(queryLoc,
						featureCoordSystems[i]);
				if (idString.length() == 0) {
					logger.warning("empty idString for loc = " + requestLoc);
				}
				if (query.length() > 0)
					query.append(" AND ");
				query.append(getPrimaryTableSynonym() + ".seq_region_id IN (" + idString + ")");

				List tmp = genericFetch(query.toString(), queryLoc);
				convertLocationsToRequestCoordinateSystem(requestLoc, tmp);
				tmp = filterNonOverlapping(requestLoc, tmp);
				features.addAll(tmp); // XXX
				// should
				// be
				// dereferencedLoc
			} else {

				// execute a separate query for each node in location

				for (Location loc = queryLoc; loc != null; loc = loc.next()) {

					// can't search on gaps because they don't correspond to a
					// specific region.
					if (loc.isGap())
						continue;

					StringBuffer tmpQuery = new StringBuffer(query.toString());

					boolean and = constraint != null && constraint.length() > 0;

					and = location2PartialSQLWhereClause(loc,
							featureCoordSystems[i], and, tmpQuery);

					List tmp = genericFetch(tmpQuery.toString(), loc);
					convertLocationsToRequestCoordinateSystem(requestLoc, tmp);
					features.addAll(tmp); // XXX
					// should
					// be
					// dereferencedLoc

				}
			}
		}

		// TODO - rebuild references

		sortByLocation(features);

		return features;

	} // fetchAllByConstraint

	/**
	 * Appends a "location" filter consiting of several conditions connected by 
	 * "and" operators to the buffer. These are designed for inclusion in an SQL "where"
	 * clauses. Does nothing if the location is a gap.
	 * 
	 * @param loc
	 *            location to create create conditions for
	 * @param cs
	 *            coordinate system to create conditions for
	 * @param and
	 *            whether an "and" should be appended to the buffer before
	 * @param buffer
	 * @return whether an "and" value should be included before the next
	 *         condition in the where clause
	 * @throws AdaptorException
	 */
	protected boolean location2PartialSQLWhereClause(Location loc,
			CoordinateSystem cs, 
			boolean initialAnd, 
			StringBuffer buffer)
			throws AdaptorException {

		if (loc.isGap())
			return initialAnd;

		boolean and = initialAnd;
		String tab_syn = getPrimaryTableSynonym();

		if (loc.getSeqRegionName() != null) {
			long seqRegionID = locationConverter.nameToId(loc
					.getSeqRegionName(), cs);
			if (and)
				buffer.append(" AND ");
			buffer.append(tab_syn + ".seq_region_id = " + seqRegionID);
			and = true;
		}
		if (loc.getEnd() > 0) {
			if (and)
				buffer.append(" AND ");
			buffer.append(tab_syn).append(".seq_region_start <= ").append(
					loc.getEnd());
			and = true;
		}

		if (loc.getStart() > 0) {
			if (and)
				buffer.append(" AND ");
			buffer.append(tab_syn).append(".seq_region_end >= ").append(
					loc.getStart());
			and = true;

			// Optimisation for quickly finding features at
			// the "end" of a sequence region by quickly
			// filtering features before the location
			// (nearer the the beginning of the sequence region).
			//
			// Problem: The "start" column is indexed in the database
			// but the "end" column is not. Because of the
			// the way the above query constraints work we
			// can quickly filter features where
			// feature.start>location.end (after the location)
			// but we have to sequentially filter those where
			// feature.end<location.start (before the location).
			//
			// Solution: If the maximum length for the feature
			// is available we can use that to create a
			// further constraint on the "start"
			// column in the database. This means we can quickly
			// filter most features that appear before location.
			int max = getMaxFeatureLength(loc.getCoordinateSystem(),
					getPrimaryTableName());
			if (max > 0) {
				int minStart = loc.getStart() - max;
				buffer.append(" AND ").append(tab_syn).append(
						".seq_region_end >= ").append(minStart);
			}

		}
		if (loc.getStrand() != 0) {
			if (and)
				buffer.append(" AND ");
			buffer.append(tab_syn).append(".seq_region_strand = ").append(
					loc.getStrand());
			and = true;
		}

		return and;
	}

	/**
	 * Sorts objects by location if the objects are Features, otherwise does
	 * nothing.
	 * 
	 * EnsJ returns Features sorted by Location.
	 * 
	 * @param objects
	 *            zero or more objects of the same type
	 */
	private void sortByLocation(List features) {
		if (features.size() > 0 && features.get(0) instanceof Feature)
			Collections.sort(features);
	}

	/**
	 * Returns the max length of this feature type from the "meta_coord" table
	 * in the database if that information is avaiable, otherwise returns -1.
	 * 
	 * Derived classes should override this method if they want different
	 * behaviour.
	 * 
	 * @param cs
	 *            coordinate system feature appears in.
	 * @param tableName
	 *            table name used to store the feature.
	 * @return max feature length if available, otherwise -1.
	 */
	protected int getMaxFeatureLength(CoordinateSystem cs, String tableName)
			throws AdaptorException {

		return ((MySQLCoordinateSystemAdaptor) driver
				.getCoordinateSystemAdaptor()).fetchMaxLength(cs, tableName);

	}

	// ---------------------------------------------------------------------

	private List filterNonOverlapping(Location requestLoc, List original) {
		if (requestLoc == null || original.size() == 0)
			return original;

		CoordinateSystem sample = ((Locatable) original.get(0)).getLocation()
				.getCoordinateSystem();
		if (!sample.getName()
				.equals(requestLoc.getCoordinateSystem().getName()))
			throw new RuntimeException(
					"RequestLocation and original are in different coordinate systems");

		List r = new ArrayList();
		for (int i = 0, n = original.size(); i < n; i++) {
			Locatable l = (Locatable) original.get(i);
			if (requestLoc.overlaps(l.getLocation()))
				r.add(l);
		}
		;
		return r;
	}

	private void convertLocationsToRequestCoordinateSystem(Location requestLoc,
			List tmp) throws AdaptorException {
		if (requestLoc == null)
			return;

		LocationConverter lc = driver.getLocationConverter();
		CoordinateSystem cs = requestLoc.getCoordinateSystem();
		for (int i = 0, n = tmp.size(); i < n; i++) {
			Locatable l = (Locatable) tmp.get(i);
			// TODO - need to recurse location conversion down in
			// gene/transcript/translation/exon classes
			// OR create g/t/t/e.setCoordinateSystem() that does the conversion.
			l.setLocation(lc.convert(l.getLocation(), cs));
		}
	}

	/**
	 * Fetches the data matching the queries as separate SQL queries.
	 * 
	 * It deals with the queries in order querySQL[0], then querySQL[1], then
	 * querySQL[2] ...This enables us to reduce the size of SQL queries fired at
	 * the DB server. The results from each batch are aggregated.
	 * 
	 * <p>
	 * Can work in conjunction with createConstraintBatches.
	 * </p>
	 * 
	 * @param querySQL
	 *            array of filters (WHERE condition), each one defines a single
	 *            batch.
	 * @param loc
	 *            location filter. This applied in each batch.
	 * @return list of zero or more matching features
	 * @throws AdaptorException
	 *             if a problem occurs retrieving the data
	 * @see MySQLBaseFeatureAdaptor#createConstraintBatches(String, String)
	 */
	protected List genericFetch(String[] querySQL, Location loc)
			throws AdaptorException {
		List r = new ArrayList();
		for (int i = 0; i < querySQL.length; i++) {
			genericFetch(querySQL[i], loc, r);
		}
		sortByLocation(r);
		return r;
	}

	protected List genericFetch(String[] querySQL) throws AdaptorException {
		return genericFetch(querySQL, null);
	}

	protected List genericFetch(String querySQL, Location dereferencedLoc)
			throws AdaptorException {
		return genericFetch(querySQL, dereferencedLoc, new ArrayList());
	}


	protected List genericFetch(String querySQL, Location dereferencedLoc,
			List features) throws AdaptorException {
		return genericFetch(querySQL, dereferencedLoc, features, finalClause());
	}
	
	/**
	 * Executes the SQL and puts the resulting objects into features.
	 * 
	 * Note that derefencedLoc is NOT used as a filter.
	 * 
	 * @param querySQL
	 * @param dereferencedLoc
	 * @param features
	 * @param finalClause
	 * @return features.
	 * @throws AdaptorException
	 */
	protected List genericFetch(String querySQL, Location dereferencedLoc,
			List features, String finalClause) throws AdaptorException {

		String[][] tableArray = null;

		// --------------------------------
		// Build LEFT JOIN clause if one is defined
		// tableArray is built from the table names that are NOT in the left
		// join list

		String[][] leftJoinArray = leftJoin();
		StringBuffer leftJoin = new StringBuffer();

		if (leftJoinArray.length == 0) {

			tableArray = tables();

		} else {
			List tableList = new ArrayList();
			String[][] tmpTabs = tables();
			for (int i = 0; i < tmpTabs.length; i++) {
				int lja;
				if ((lja = leftJoinArrayContains(leftJoinArray, tmpTabs[i][0])) > -1) {
					String syn = tmpTabs[i][1];
					String ljTable = leftJoinArray[lja][0];
					String ljCondition = leftJoinArray[lja][1];
					leftJoin.append(" LEFT JOIN " + ljTable + " " + syn
							+ " ON " + ljCondition);
				} else {
					tableList.add(tmpTabs[i]); // table not used in left join
				}
			}

			// convert tableList to tableArray
			tableArray = new String[tableList.size()][2];
			Iterator it = tableList.iterator();
			int ta = 0;
			while (it.hasNext()) {
				String[] tablePair = (String[]) it.next();
				tableArray[ta++] = tablePair;
			}
		}

		// --------------------------------

		String cols = columnsToString(columns());
		String tabs = tablesToString(tableArray);

		StringBuffer query = new StringBuffer();
		query.append("SELECT " + cols + " FROM " + tabs + " "
				+ leftJoin.toString());

		// --------------------------------
		// append a WHERE clause if one was defined
		if (querySQL.length() > 0) {
			query.append(" WHERE " + querySQL);
			if (finalWhereClause().length() > 0) {
				query.append(" AND ").append(finalWhereClause());
			}
		} else if (finalWhereClause().length() > 0) {
			query.append("WHERE " + finalWhereClause());
		}

		// append additional clauses which may have been defined
		query.append(" ").append(finalClause);

		// System.out.println("About to execute: " + query.toString());
		Connection conn = null;

		try {

			conn = getConnection();
			ResultSet rs = executeQuery(conn, query.toString());

			Object feature;
			while ((feature = createObject(rs)) != null) {
				// TODO - convert back to dereferencedLoc coordinate system
				// TODO - filtering
				features.add(feature);
			} // while feature
		} finally {
			close(conn);
		}
		return features;

	}

	// ---------------------------------------------------------------------

	public String columnsToString(String[] columns) {

		StringBuffer result = new StringBuffer();
		for (int i = 0; i < columns.length; i++) {
			result.append(columns[i]);
			if (i < columns.length - 1) {
				result.append(", ");
			}
		}
		return result.toString();
	}

	// ---------------------------------------------------------------------

	public String tablesToString(String[][] tables) {

		StringBuffer result = new StringBuffer();

		for (int i = 0; i < tables.length; i++) {
			result.append(tables[i][0]); // table name
			result.append(" ");
			result.append(tables[i][1]); // synonum
			if (i < tables.length - 1) {
				result.append(", ");
			}
		}
		return result.toString();

	}

	// -----------------------------------------------------------------

	public String getPrimaryTableName() {

		String[][] tables = tables();
		return tables[0][0];

	}

	// -----------------------------------------------------------------

	public String getPrimaryTableSynonym() {

		String[][] tables = tables();
		return tables[0][1];

	}

	// ---------------------------------------------------------------------

	private String buildIDString(Location loc, CoordinateSystem cs)
			throws AdaptorException {

		StringBuffer result = new StringBuffer();
		long[] ids = locationConverter.locationToIds(loc);

		for (int j = 0; j < ids.length; j++) {
			result.append(ids[j]);
			if (j < ids.length - 1) {
				result.append(", ");
			}

		}

		return result.toString();

	}

	/**
	 * To be overridden by derived classes where necessary.
	 * 
	 * Default implementation ignores "loadChildren" and delegate to fetch(loc).
	 * 
	 * 
	 * @param loc
	 *            location filter
	 * @param loadChildren
	 *            ignored.
	 * @return zero or more features.
	 * @throws AdaptorException
	 * @see #fetch(Location)
	 */
	public List fetch(Location loc, boolean loadChildren)
			throws AdaptorException {
		return fetch(loc);
	}

	// ----------------------------------------------------------------

	/**
	 * Fetch features that overlap the location.
	 * 
	 * @param loc
	 *            location filter
	 * @return zero or more features.
	 * @throws AdaptorException
	 * @see #fetch(Location)
	 */
	public List fetch(Location loc) throws AdaptorException {

		// Optimisation: if the sequence region is not specified
		// it is generally faster to do a fetchAll()
		// and then remove irrelevant items that do not
		// map to loc.coordinateSyatem than doing many
		// fetch by sequence regions (which fetchAllByConstraint()
		// will resolve to.

		if (loc.getSeqRegionName() == null) {

			List r = fetchAll();
			for (int i = 0; i < r.size();) {
				Locatable f = (Locatable) r.get(i);
				Location tmp = locationConverter.convert(f.getLocation(), loc
						.getCoordinateSystem());
				if (tmp == null) {
					// remove feature if it can't be mapped to cs of interest
					r.remove(i);
				} else {
					f.setLocation(tmp);
					i++;
				}

			}
			return r;

		} else {
			return fetchAllByConstraint(loc, "");
		}

	}

	// -----------------------------------------------------------------

	/**
	 * Returns an iterator over all Features of this type in the database.
	 * 
	 * Returns the same Features as FetchAll() but uses much less memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @return iterator over all Features of this type in the database.
	 * @see #fetchAll() an alternative way of getting all the features.
	 * @see #DEFAULT_FEATURE_BUFFER_SIZE
	 * @see org.ensembl.driver.FeatureIterator
	 * @throws AdaptorException
	 */
	public Iterator fetchIterator() throws AdaptorException {
		return new FeatureIterator((FeatureAdaptor) this,
				DEFAULT_FEATURE_BUFFER_SIZE, false);
	}

	/**
	 * Returns the same Features as fetch(long[], boolean) but uses much less
	 * memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @see org.ensembl.driver.FeatureAdaptor#fetchIterator(long[], boolean)
	 */
	public Iterator fetchIterator(long[] internalIDs, boolean loadChildren)
			throws AdaptorException {
		return new FeatureIterator((FeatureAdaptor) this,
				DEFAULT_FEATURE_BUFFER_SIZE, loadChildren, internalIDs);
	}

	/**
	 * Returns the same Features as fetch(long[]) but uses much less memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @see org.ensembl.driver.FeatureAdaptor#fetchIterator(long[])
	 */
	public Iterator fetchIterator(long[] internalIDs) throws AdaptorException {
		return fetchIterator(internalIDs, false);
	}

	/**
	 * Returns an iterator over all Features of this type in the database.
	 * 
	 * Returns the same Features as FetchAll() but uses much less memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @param loadChildren
	 *            whether to preload child data if available.
	 * @return iterator over all Features of this type in the database.
	 * @see #fetchAll() an alternative way of getting all the features.
	 * @see #DEFAULT_FEATURE_BUFFER_SIZE
	 * @see org.ensembl.driver.FeatureIterator
	 * @throws AdaptorException
	 */
	public Iterator fetchIterator(boolean loadChildren) throws AdaptorException {
		return new FeatureIterator((FeatureAdaptor) this,
				DEFAULT_FEATURE_BUFFER_SIZE, loadChildren);
	}

	/**
	 * Returns an iterator over all Features in the specified location.
	 * 
	 * Returns the same Features as fetch(location) but uses potentially much
	 * less memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @param location
	 *            location filter.
	 * @return iterator over all Features that fall in the specified location.
	 * @see #fetch(Location) an alternative way of getting all the features.
	 * @see #DEFAULT_FEATURE_BUFFER_SIZE
	 * @see org.ensembl.driver.FeatureIterator *
	 * @throws AdaptorException
	 */
	public Iterator fetchIterator(Location location) throws AdaptorException {
		return new FeatureIterator((FeatureAdaptor) this,
				DEFAULT_FEATURE_BUFFER_SIZE, false, location);
	}

	/**
	 * Returns an iterator over all Features in the specified location.
	 * 
	 * Returns the same Features as fetch(location) but uses potentially much
	 * less memory.
	 * 
	 * The method uses a buffer of features with a size set to
	 * DEFAULT_FEATURE_BUFFER_SIZE. If you want to use a different buffer size
	 * you should create your own FeatureIterator.
	 * 
	 * @param location
	 *            location filter.
	 * @return iterator over all Features that fall in the specified location.
	 * @see #fetch(Location) an alternative way of getting all the features.
	 * @see #DEFAULT_FEATURE_BUFFER_SIZE
	 * @see org.ensembl.driver.FeatureIterator *
	 * @throws AdaptorException
	 */
	public Iterator fetchIterator(Location location, boolean loadChildren)
			throws AdaptorException {
		return new FeatureIterator((FeatureAdaptor) this,
				DEFAULT_FEATURE_BUFFER_SIZE, loadChildren, location);
	}

	// -----------------------------------------------------------------

	/**
	 * Fetch zero or more Features with the specified internalIDs
	 * 
	 * @param internalIDs
	 *            internal IDs
	 * @return zero or more Features with the specified internalIDs.
	 * @throws AdaptorException
	 */
	public List fetch(long[] internalIDs) throws AdaptorException {

		if (internalIDs.length == 0)
			return Collections.EMPTY_LIST;

		StringBuffer sqlConstraint = new StringBuffer();
		sqlConstraint.append(getPrimaryTableSynonym());
		sqlConstraint.append(".");
		sqlConstraint.append(getPrimaryTableName());
		sqlConstraint.append("_id IN (");
		sqlConstraint.append(StringUtil.toString(internalIDs));
		sqlConstraint.append(")");
		
		List r = fetchByNonLocationConstraint(sqlConstraint.toString());
		Collections.sort(r, new InternalIDOrderComparator(internalIDs));
		return r;

	}

	/**
	 * Default implementation ignores the loadChildren hint and just fetches the
	 * Features. Sub classes that wish to implement the prefetch behaviour
	 * should over ride this method.
	 * 
	 * @param internalIDs
	 *            feature internalIDs.
	 * @param loadChildren
	 *            ignored in this default implementation.
	 * @return zero or more Features with the specified internalIDs.
	 * @throws AdaptorException
	 */
	public List fetch(long[] internalIDs, boolean loadChildren)
			throws AdaptorException {
		return fetch(internalIDs);
	}

	// ---------------------------------------------------------------------

	protected Persistent fetchByInternalID(long internalID)
			throws AdaptorException {

		StringBuffer sqlConstraint = new StringBuffer();
		sqlConstraint.append(getPrimaryTableSynonym());
		sqlConstraint.append(".");
		sqlConstraint.append(getPrimaryTableName());
		sqlConstraint.append("_id=");
		sqlConstraint.append(internalID);

		List features = fetchByNonLocationConstraint(sqlConstraint.toString());

		if (features.size() == 0)
			return null;

		else if (features.size() == 1)
			return (Persistent) features.get(0);

		else
			throw new AdaptorException(
					"Unable to fetch by internal ID - expected 1 feature, got "
							+ features.size() + " internalID=" + internalID);

	}

	// ---------------------------------------------------------------------

	protected List fetchByNonLocationConstraint(String constraint)
			throws AdaptorException {

		List features = genericFetch(constraint, null);
		sortByLocation(features);

		return features;
	}

	// -----------------------------------------------------------------
	/**
	 * Check if a two-dimensional array contains a particular String in its
	 * first dimension.
	 * 
	 * @param arr
	 *            The array to search.
	 * @param str
	 *            The string to search for.
	 * @return index of location of str, or -1 if not found.
	 */
	private int leftJoinArrayContains(String[][] arr, String str) {
		if (arr.length > 0) {
			for (int i = 0; i < arr.length; i++) {
				if (arr[i].length > 0) {
					if (arr[i][0].equals(str)) {
						return i;
					}
				}
			}
		}

		return -1;

	}

	// -----------------------------------------------------------------

	/**
	 * @return true if the loadChildren flag is set.
	 */
	public boolean getIncludeChildren() {
		return loadChildren;
	}

	/**
	 * @param loadChildren
	 */
	public void setIncludeChildren(boolean loadChildren) {
		this.loadChildren = loadChildren;
	}

	/**
	 * Creates an array of constraints that can be fed into the WHERE clause of
	 * an SQL statement. Each constraint is of the form "column = X", or "column
	 * in (X,Y,Z,A,B....)" where X,Y,Z,A,B are values retrieved from the first
	 * column of the result set.
	 * 
	 * <p>
	 * It is useful to batch the ids into chunks to improve retrieval
	 * performance without making the (final) SQL queries too long for the db
	 * engine. This method can be used in conjunction with
	 * genericFetch(String[]) which can take the output of this method for
	 * input.
	 * </p>
	 * 
	 * @param column
	 *            column name, might need to be qualified with table name
	 * @param sql
	 *            SQL statement where the first column in the SELECT clause
	 *            contain the values to be inserted into the constraint after
	 *            the "=" or between the brackets of the IN list.
	 * @return an array containing zero or more constraints
	 * @see MySQLBaseFeatureAdaptor#genericFetch(String[])
	 */
	protected String[] createConstraintBatches(String column, String sql)
			throws AdaptorException {

		Connection conn = null;

		conn = getConnection();
		ResultSet rs = executeQuery(conn, sql);

		// read the ids from the db.
		List idBuf = new ArrayList();
		try {
			while (rs.next())
				idBuf.add(rs.getString(1));
		} catch (SQLException e) {
			throw new AdaptorException("Failed to get ids from SQL:" + sql, e);
		} finally {
			close(conn);
		}
		final int idBufSize = idBuf.size();

		// group the ids into chunks and then create a constraint for each chunk
		final int CHUNK_SIZE = 1000;
		final int extra = (idBufSize % CHUNK_SIZE == 0) ? 0 : 1;
		final int nConstraints = (int) (idBufSize / CHUNK_SIZE) + extra;
		String[] constraints = new String[nConstraints];
		int b = 0;
		int c = 0;
		while (b < idBufSize) {

			if (b + 1 == idBufSize) {
				// column = X
				constraints[c++] = column + " = " + (String) idBuf.get(b++);
				continue;

			} else {

				// column IN (X,Y,Z,A,B,C....)
				StringBuffer constraint = new StringBuffer();
				constraint.append(column).append(" IN (");
				int i = 0;
				while (b < idBufSize && i < CHUNK_SIZE) {
					if (i++ > 0)
						constraint.append(",");
					constraint.append(idBuf.get(b++));
				}
				constraint.append(")");
				constraints[c++] = constraint.toString();
			}

		}

		return constraints;
	}

}
