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

package org.ensembl.probemapping;

import java.util.ArrayList;
import java.util.List;

import org.ensembl.util.IDSet;

/**
 * A partial representation of a microarray probe set / composite
 * used for probeset to transcript mapping.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 * @see org.ensembl.datamodel.AffyProbe
 * @see org.ensembl.datamodel.AffyArray
 */
public class ProbeSet {

	public final String probeSetName;
	
	/**
	 * List of zero or more AffyArrays.
	 * @see org.ensembl.datamodel.AffyArray
	 */
	public final List arrays;
	
	/**
	 * List of zero or more MappableAffyFeatures that belong to this 
	 * probeset.
	 * @see MappableAffyFeature
	 */
	public final List affyFeatures = new ArrayList();

	private List mappedTranscripts;

	/**
	 * True if the probe set hits too many transcripts to be mapped to any 
	 * of them.
	 */
	public boolean tooManyTranscripts = false;

	/**
	 * @param probeSetName name of the probe sets
	 * @param arrays microarrays that contain this probeset.
	 */
	public ProbeSet(String probeSetName, List arrays) {
		this.probeSetName = probeSetName;
		this.arrays = arrays;
	}

	/**
	 * Adds a Mappable affy feature that belongs to this probe set.
	 * @param feature Mappable affy feature that belongs to this probe set.
	 */
	public void addMappableAffyFeature(MappableAffyFeature feature) {
		affyFeatures.add(feature);
		mappedTranscripts = null;
	}

	/**
	 * Returns a unique set of MappableTranscripts where 
	 * each one hits at least 
	 * one of the probes in this probe set.
	 * @return zero or more MappableTranscripts.
	 * @see MappableTranscript
	 */
	public List getOverlappingTranscripts() {
		if (mappedTranscripts == null) {
			IDSet uniqueIDs = new IDSet();
			mappedTranscripts = new ArrayList();
			for (int i = 0, n = affyFeatures.size(); i < n; i++) {
				MappableAffyFeature af = (MappableAffyFeature) affyFeatures.get(i);
				List transcripts = af.transcripts;
				for (int j = 0, z = transcripts.size(); j < z; j++) {
					MappableTranscript t = (MappableTranscript)transcripts.get(j);
					if (!uniqueIDs.contains(t)){
					mappedTranscripts.add(t);
					uniqueIDs.add(t);
					}
				}
			}
			
		}
		return mappedTranscripts;
	}
}
