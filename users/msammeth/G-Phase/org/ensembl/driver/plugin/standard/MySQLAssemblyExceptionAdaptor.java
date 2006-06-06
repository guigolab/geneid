/*
 * Created on 09-Nov-2003
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package org.ensembl.driver.plugin.standard;

import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.ensembl.datamodel.AssemblyException;
import org.ensembl.datamodel.CloneFragmentLocation;
import org.ensembl.datamodel.Location;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.AssemblyExceptionAdaptor;
import org.ensembl.driver.LocationConverter;
/**
 * @author arne
 * This class delivers contents of the assembly_exception table to
 * the LocationConverter and to display code 
 */
public class MySQLAssemblyExceptionAdaptor
	extends MySQLBaseFeatureAdaptor
	implements AssemblyExceptionAdaptor {

	private List exceptionCache;

	public MySQLAssemblyExceptionAdaptor(MySQLDriver driver)
		throws AdaptorException {
		super(driver, "AssemblyExceptions");
	}

	private void loadExceptions() {

	}

	public String[] columns() {
		String columns[] =
			{
				"ae.seq_region_id",
				"ae.seq_region_start",
				"ae.seq_region_end",
				"ae.exc_seq_region_id",
				"ae.exc_seq_region_start",
				"ae.exc_seq_region_end",
				"ae.ori",
				"ae.exc_type" };
		return columns;
	}

	public String getType() {
		return TYPE;
	}

	public String[][] tables() {
		String[][] tables = { { "assembly_exception", "ae" }
		};
		return tables;
	}

	public AssemblyException[] fetchLinks(Location loc) throws AdaptorException {
		if (exceptionCache == null) {
			exceptionCache = genericFetch("", null);
		}

		// read result from cache 
		ArrayList result = new ArrayList();
		Iterator it = exceptionCache.iterator();
		while (it.hasNext()) {
			AssemblyException ae = (AssemblyException) it.next();
			if (ae.linkedLocation.getSeqRegionName().equals(loc.getSeqRegionName())
				&& ae.linkedLocation.getCoordinateSystem().equals(
					loc.getCoordinateSystem())) {
				result.add(ae);
			}
		}

		AssemblyException aes[] = new AssemblyException[result.size()];
		result.toArray(aes);
		return aes;
	}

	/* (non-Javadoc)
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#fetch(org.ensembl.datamodel.CloneFragmentLocation, java.lang.String)
	 */
	protected List fetch(
		CloneFragmentLocation location,
		String analysisIDCondition)
		throws AdaptorException {
		
		return null;
	}

	/* (non-Javadoc)
	 * @see org.ensembl.driver.plugin.standard.MySQLBaseFeatureAdaptor#createObject(java.sql.ResultSet)
	 */
	public Object createObject(ResultSet rs) throws AdaptorException {
		LocationConverter lc = driver.getLocationConverter();
		AssemblyException ae = null;
		try {

      if ( rs.next() ) {
      
			Location linkLocation =
				lc.idToLocation(rs.getLong(1), rs.getInt(2), rs.getInt(3), 1);
			Location orgLocation =
				lc.idToLocation(
					rs.getLong(4),
					rs.getInt(5),
					rs.getInt(6),
					rs.getInt(7));
			ae = new AssemblyException(orgLocation, linkLocation, rs.getString(8));
			ae.setDriver(getDriver());
      
      }
		} catch (SQLException e) {
			throw new AdaptorException("rethrow ", e);
		}
		return ae;
	}

}
