/*
    Copyright (C) 2001 EBI, GRL

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
package org.ensembl.util;


import java.util.Collection;
import java.util.Iterator;

import org.ensembl.datamodel.Persistent;

/**
 * Convenience class for managing sets of ids taken from 
 * Persistent item's internal ids.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public class IDSet extends LongSet {

  
	private static final long serialVersionUID = 1L;

	public IDSet() {
	}

	/**
	 * Creates a set from a collection of Persistent objects.
	 * @param c collection of Persistent objects.
	 */
	public IDSet(Collection c) {
		for (Iterator iter = c.iterator(); iter.hasNext();) 
		 add((Persistent) iter.next());
	}
	
  public boolean add(Persistent persistent) {
    return add(new Long(persistent.getInternalID()));
  }
  
  public boolean contains(Persistent persistent) {
    return contains(new Long(persistent.getInternalID()));
  }
}
