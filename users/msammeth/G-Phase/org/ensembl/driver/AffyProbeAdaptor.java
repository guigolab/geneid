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
import org.ensembl.datamodel.AffyProbe;

/**
 * Adaptor for retrieving AffyProbes from the database.
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public interface AffyProbeAdaptor extends Adaptor {

	final String TYPE = "affy_probe";

	AffyProbe fetch(long internalID) throws AdaptorException;
	
	/**
	 * Returns AffyProbes that appear in the specified array.
	 * @param array micro array from which to fetch AffyProbes
	 * @return zero or more AffyProbes for the specified array.
	 */
	List fetch(AffyArray array) throws AdaptorException;

	/**
	 * Fetch Probes from the specified probe set.
	 * @param probeSetName name of the probe set.
	 * @return zero or more AffyProbes
	 * @throws AdaptorException
	 * @see AffyProbe
	 */
	List fetch(String probeSetName) throws AdaptorException;
	
}
