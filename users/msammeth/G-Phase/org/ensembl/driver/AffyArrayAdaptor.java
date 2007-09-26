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

/**
 * Adaptor for retrieving AffyArrays.
 *
 * @see org.ensembl.datamodel.AffyArray
 *
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 */
public interface AffyArrayAdaptor extends Adaptor {

	final String TYPE = "affy_array";
	
	/**
	 * Retrieve AffyArray with specified internal ID.
	 * @param internalID internal ID of AffyArray in database.
	 * @return AffyArray with specified internal ID, or null if
	 * none found.
	 * @throws AdaptorException
	 */
	AffyArray fetch(long internalID) throws AdaptorException;

	/**
	 * Retrieve AffyArray with specified name.
	 * @param name name of an AffyArray in database.
	 * @return AffyArray with specified name, or null if
	 * none found.
	 * @throws AdaptorException
	 */
	AffyArray fetch(String name) throws AdaptorException;
	
	/**
	 * Fetch all AffyArrays from the database.
	 * @return zero or more AffyArrays.
	 * @see AffyArray
	 */
	List fetch() throws AdaptorException;
	
	
}
