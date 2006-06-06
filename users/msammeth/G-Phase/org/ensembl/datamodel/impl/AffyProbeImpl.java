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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.AffyProbe;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.RuntimeAdaptorException;
import org.ensembl.util.StringUtil;

/**
 * Implementation of the AffyProbe type.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public class AffyProbeImpl extends PersistentImpl implements AffyProbe {

	private static final long serialVersionUID = 1L;
	
	private List arrays;
	private Map names = new HashMap();
	private String probeSetName;
	private List affyFeatures;

	
	
	/**
	 * Creates an AffyProbe instance that will lazy load the AffyArrays it 
	 * appears in and it's AffyFeatures.
	 * @param internalID internalID in database.
	 * @param driver source driver.
	 * @param probeSetName name of the probeset the probe belongs to.
	 */
	public AffyProbeImpl(Driver driver, long internalID, String probeSetName) {
		super(driver);
		this.internalID = internalID;
		this.probeSetName = probeSetName;
	}
	
	
	/**
	 * @see org.ensembl.datamodel.AffyProbe#getAffyArraysContainingThisProbe()
	 */
	public List getAffyArraysContainingThisProbe() {
		if (arrays==null) arrays = new ArrayList(names.keySet());
		return arrays;
	}

	/**
	 * @see org.ensembl.datamodel.AffyProbe#getQualifiedName(org.ensembl.datamodel.AffyArray)
	 */
	public String getQualifiedName(AffyArray array) {
		return array.getName()+":"+probeSetName+":"+names.get(array);
	}

	/**
	 * @see org.ensembl.datamodel.AffyProbe#getQualifiedNames()
	 */
	public String[] getQualifiedNames() {
		String[] qNames = new String[getAffyArraysContainingThisProbe().size()];
		for (int i = 0; i < qNames.length; i++) 
			qNames[i] = getQualifiedName((AffyArray) arrays.get(i));
		return qNames;
	}

	
	public boolean isProbeInArray(AffyArray array) {
		return names.containsKey(array);
	}
	
	/**
	 * @see org.ensembl.datamodel.AffyProbe#getNameInArray(AffyArray)
	 */
	public String getNameInArray(AffyArray array) {
		return (String) names.get(array);
	}

	/**
	 * @see org.ensembl.datamodel.AffyProbe#getProbeSetName()
	 */
	public String getProbeSetName() {
		return probeSetName;
	}

	/**
	 * @see org.ensembl.datamodel.AffyProbe#getAffyFeatures()
	 */
	public List getAffyFeatures() {
		if (affyFeatures==null)
			try {
				affyFeatures = driver.getAffyFeatureAdaptor().fetch(this);
			} catch (AdaptorException e) {
				throw new RuntimeAdaptorException("Failed to lazy load AffyFeatures for AffyProbe: "+this,e);
			}
		return affyFeatures;
	}

	
	public String toString() {
		StringBuffer buf = new StringBuffer();

		buf.append("[");
		buf.append(super.toString());
		buf.append(", names").append(StringUtil.toString(getQualifiedNames()));
		buf.append(", probeSetName").append(probeSetName);
		buf.append(", affyFeatures").append(StringUtil.sizeOrUnset(affyFeatures));
		buf.append(", parentAffyArrays").append(StringUtil.sizeOrUnset(arrays));
		buf.append("]");

		return buf.toString();
	}


	/**
	 * @see org.ensembl.datamodel.AffyProbe#addArrayWithProbeName(org.ensembl.datamodel.AffyArray, java.lang.String)
	 */
	public void addArrayWithProbeName(AffyArray array, String name) {
		names.put(array,name);		
		arrays = null;
	}
	
}
