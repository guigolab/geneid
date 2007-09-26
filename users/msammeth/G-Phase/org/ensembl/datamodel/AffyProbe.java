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

package org.ensembl.datamodel;

import java.util.List;

/**
 * An oligo spot from one or more Affymetric micro arrays. 
 * 
 * AffyProbes are grouped into probe sets and the same probe in the 
 * same set can appear in multiple arrays.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public interface AffyProbe extends Persistent {

	
	/**
	 * All affymetix microarrays that contain this probe.
	 * @return zero or more AffyArrays.
	 * @see AffyArray
	 */
	List getAffyArraysContainingThisProbe();
	
	/**
	 * Return the name of this probe when it's parent is _array_.
	 * 
	 * String format: arrayName:probesetName:probeName 
	 * 
	 * @param array parent array.
	 * @return fully qualified name of this probe or null if array is
	 * not one of it's parent arrays.
	 * @see #getQualifiedNames()
	 */
	String getQualifiedName(AffyArray array);
	
	/**
	 * All the possible names of this probe based on it's parent arrays.
	 * 
	 * For each parent array the string format is: arrayName:probesetName:probeName 
	 * 
	 * @return zero or more fully qualified names.
	 * @see #getQualifiedName(AffyArray)
	 */
	String[] getQualifiedNames();
	
	/**
	 * Name of this probe in specified microarray. 
	 * 
	 * This is last "component" of the probe's qualified names.
	 * 
	 * @param array microarray where this probe appears.
	 * @return name of this probe in the specified array, null if
	 * it does not appear in the array.
	 */
	String getNameInArray(AffyArray array);
	
	
	/**
	 * Name of the probeset this probe belongs to.
	 * @return Name of the probeset this probe belongs to.
	 */
	String getProbeSetName();
	
	/**
	 * The places where this probe appears on the genome.
	 * @return zero or more AffyFeatures.
	 * @see AffyFeature
	 */
	List getAffyFeatures();
	
	/**
	 * Adds the array-probename pair to the probe.
	 * 
	 * Probes can have different names in different microarrays.
	 * 
	 * @param array microarray.
	 * @param name probe name in the specified array.
	 */
	void addArrayWithProbeName(AffyArray array, String name);
	
	/**
	 * Whether this probe appears in the array.
	 * @param array microarray.
	 * @return true is this probe appears in the array, otherwise false.
	 */
	boolean isProbeInArray(AffyArray array);
}
