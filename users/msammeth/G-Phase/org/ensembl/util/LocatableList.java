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

import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;

import org.ensembl.datamodel.Locatable;

/**
 * Container for locatables that offers a convenient
 * means of retrieving an array of the locatable's internalIDs
 * where the locatables are sorted by genomic location.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public class LocatableList {


  private LinkedList list = new LinkedList();
  private int size = 0;
  private long[] result = null;


  public void add(Locatable l) {
    list.add(l);
    size++;
    result=null;
  }
  
  /**
   * Sorts locatable's by genomic location and returns
   * an array of the locatable's internalIDs.
   * @return zero or more internalIDs.
   */
  public long[] toSortedInternalIDArray() {
    if (result == null) {
      Collections.sort(list);
      result = new long[size];
      Iterator iter = list.listIterator();
      for (int i = 0; i < size; i++)
        result[i] = ((Locatable) iter.next()).getInternalID();
    }
    return result;
  }
  
  public int size() {
    return size;
  }
}
