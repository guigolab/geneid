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

package org.ensembl.datamodel.impl;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Properties;
import java.util.Set;
import java.util.logging.Logger;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.ExternalDatabase;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.RuntimeAdaptorException;
import org.ensembl.util.PropertiesUtil;
import org.ensembl.util.StringUtil;

/**
 * Implementation of AffyArray.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class AffyArrayImpl extends PersistentImpl implements AffyArray {

	/**
	 * Used by the (de)serialization system to determine if the data in a
	 * serialized instance is compatible with this class.
	 * 
	 * It's presence allows for compatible serialized objects to be loaded when
	 * the class is compatible with the serialized instance, even if:
	 * 
	 * <ul>
	 * <li>the compiler used to compile the "serializing" version of the class
	 * differs from the one used to compile the "deserialising" version of the
	 * class.</li>
	 * 
	 * <li>the methods of the class changes but the attributes remain the same.
	 * </li>
	 * </ul>
	 * 
	 * Maintainers must change this value if and only if the new version of this
	 * class is not compatible with old versions. e.g. attributes change. See
	 * Sun docs for <a
	 * href="http://java.sun.com/j2se/1.4.2/docs/guide/serialization/"> details.
	 * </a>
	 * 
	 */
	private static final long serialVersionUID = 1L;

	private final static Logger logger = Logger.getLogger(AffyArrayImpl.class
			.getName());

	private String name;

	private int probeSetSize;

	private List affyProbes;

	private static Properties defaultProbeSetSizes;

	private ExternalDatabase externalDatabase;

	/**
	 * Set of array names that have default probe set sizes.
	 *  Used to prevent giving the same warning many times.
	 */ 
	private static Set warned = new HashSet();

	/**
	 * Creates an AffyArray with the specified fields.
	 * 
	 * Other fields will be lazy loaded from database.
	 * 
	 * @param internalID
	 *            internalID of this array in the database.
	 * @param driver
	 *            parent driver.
	 * @param name
	 *            name of the array.
	 * @param probeSetSize
	 *            number of probes in a typical probe set.
	 */
	public AffyArrayImpl(Driver driver, long internalID, String name,
			int probeSetSize) {
		super(driver);
		this.internalID = internalID;
		this.name = name;
		setProbeSetSize(probeSetSize, name);

	}

	/**
	 * Sets the probesetSize to size iof size>0, otherwise looks up this array
	 * in config file by name.
	 * 
	 * @param size
	 *            number of probes in each standard probeset in this array.
	 */
	private void setProbeSetSize(int size, String name) {
		if (size > 0) {
			probeSetSize = size;
		} else {
			String defaultSize = getDefaultProbeSetSizes().getProperty(
					formatName(name));
			if (defaultSize != null) {
				probeSetSize = Integer.parseInt(defaultSize);
				if (!warned .contains(name)) {
				logger
						.warning("Using default probe set size for microarray "
								+ name
								+ " ("
								+ probeSetSize
								+ "). The correct value should be set in 'affay_array.probeset_size' in the database.");
				warned.add(name);
				}
			}
		}
		if (probeSetSize < 1)
			throw new RuntimeAdaptorException(
					"Invalid probe set size for microarray: " + name);

	}

	/**
	 * We need to convert the microarray names to a standard 'format' so that we
	 * can do hash lookups.
	 * 
	 * @param name
	 * @return
	 */
	private String formatName(String name) {
		return name.replace('-', '_').toLowerCase();
	}

	private synchronized Properties getDefaultProbeSetSizes() {

		if (defaultProbeSetSizes == null) {

			defaultProbeSetSizes = PropertiesUtil
					.createProperties("resources/data/microarray_composite_size.properties");

			// add some aliases for each microarray name
			List names = new ArrayList(defaultProbeSetSizes.keySet());
			for (int i = 0, n = names.size(); i < n; i++) {

				String name = (String) names.get(i);
				String size = defaultProbeSetSizes.getProperty(name);

				String name2 = formatName(name.substring(5)); // remove AFFY_
				// prefix
				defaultProbeSetSizes.put(name2, size);

			}

		}

		return defaultProbeSetSizes;
	}

	/**
	 * @see org.ensembl.datamodel.AffyArray#getName()
	 */
	public String getName() {
		return name;
	}

	/**
	 * @see org.ensembl.datamodel.AffyArray#getProbeSetSize()
	 */
	public int getProbeSetSize() {
		return probeSetSize;
	}

	/**
	 * @see org.ensembl.datamodel.AffyArray#getAffyProbes()
	 * @throws RuntimeAdaptorException
	 *             if a problem occurs lazy loading the probes.
	 */
	public List getAffyProbes() {
		if (affyProbes == null)
			try {
				affyProbes = driver.getAffyProbeAdaptor().fetch(this);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException(
						"Failed to lazy load AffyProbes for AffyArray:" + this,
						e);
			}
		return affyProbes;
	}

	public String toString() {
		StringBuffer buf = new StringBuffer();

		buf.append("[").append(super.toString());
		buf.append(", name = ").append(name);
		buf.append(", probeSetSize = ").append(probeSetSize);
		buf.append(", affyProbes = ")
				.append(StringUtil.sizeOrUnset(affyProbes));
		buf.append(", externalDatabase = ").append(getExternalDatabase());
		buf.append("]");

		return buf.toString();
	}

	/**
	 * @see org.ensembl.datamodel.AffyArray#getExternalDatabase()
	 */
	public ExternalDatabase getExternalDatabase() {
		try {
			
			// We have to do some guessing and pattern matching to find the
			// external database corresponding to this AffyArray because the
			// names used in the relevant database tables do not quite match.
			
			String[] aliases = { name.toLowerCase(),
					name.toLowerCase().replace('-', '_'),
					("affy_" + name).toLowerCase().replace('-', '_'),
					("afyy_" + name).toLowerCase().replace('-', '_')
					};
			
			List xdbs = driver.getExternalDatabaseAdaptor().fetch();
			for (int i = 0, n = xdbs.size(); externalDatabase == null && i < n; i++) {
				ExternalDatabase xdb = (ExternalDatabase) xdbs.get(i);
				String xdbName = xdb.getName().toLowerCase();
				for (int j = 0; externalDatabase==null && j < aliases.length; j++) 
					if (aliases[j].equals(xdbName)) externalDatabase = xdb;
			}
			
			if (externalDatabase == null)
				logger
						.warning("Failed to lazy load External Database for AffyArray: "
								+ name);
		} catch (AdaptorException e) {
			throw new RuntimeAdaptorException(
					"Failed to lazy load External Database for AffyArray: "
							+ name, e);
		}
		return externalDatabase;
	}
}
