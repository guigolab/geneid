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

package org.ensembl.driver;

import java.util.List;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.AffyFeature;
import org.ensembl.datamodel.AffyProbe;
import org.ensembl.datamodel.Location;

/**
 * Adaptor for retrieving AffyFeatures from a database.
 *
 * @see AffyFeature
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public interface AffyFeatureAdaptor extends FeatureAdaptor {

	final String TYPE = "affy_feature";

	/**
	 * Fetch AffyFeature with specified internalID.
	 * @param internalID internalID of the AffyFeature.
	 * @return AffyFeature with specified internalID or null if none found.
	 */
	AffyFeature fetch(long internalID) throws AdaptorException;

	/**
	 * Fetches AffyFeatures corresponding the the affyProbe.
	 * @param affyProbe AffyProbe of interest.
	 * @return zero or more AffyFeatures representing where the affyProbe hits 
	 * the genome.
	 * @throws AdaptorException
	 * @see org.ensembl.datamodel.AffyFeature
	 */
	List fetch(AffyProbe affyProbe) throws AdaptorException;

	/**
	 * Retrieve the genomic hits for the probes in the specified array.
	 * 
	 * This will return a large amount of data for some microarrays.
	 * 
	 * @param array micro array.
	 * @return zero or more AffyFeatures corresponding to the specified array.
	 * @throws AdaptorException
	 */
	List fetch(AffyArray array) throws AdaptorException;

	/**
	 * Fetches AffyFeatures that overlap with location and which
	 * appear in the specified array.
	 * 
	 * @param loc filter location condition.
	 * @param array filter microarray condition.
	 * @return zero or more AffyFeatures.
	 * @throws AdaptorException
	 * @see org.ensembl.datamodel.AffyFeature
	 */
	List fetch(Location loc, AffyArray array) throws AdaptorException;

	/**
	 * Like fetch(Location) but excludes duplicate AffyFeatures which share the
	 * same probe and location but differ by their microarray only.
	 * 
	 * @param locationFilter location of interest
	 * @return zero or more AffyFeatures.
	 * @throws AdaptorException
	 */
	List fetchUniqueProbeAndLocation(Location locationFilter) throws AdaptorException;

}
