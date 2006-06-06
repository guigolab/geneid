/*
 * Copyright (C) 2004 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License as published by the Free Software Foundation; either version
 * 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
 * even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 */

package org.ensembl.idmapping;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

/**
 * Envelope data structure to store a set of interlinked source and target items.
 */
public class Envelope {


    private List entries;

    //  -------------------------------------------------------------------------

    public Envelope() {
    		entries = new ArrayList();
    }

    public void addEntry(Entry e) {

            entries.add(e);

    }

    	public int size() {
    		return entries.size();
    	}
    	
    // -------------------------------------------------------------------------
    /**
     * Return the entries, sorted by score, <em>with highest score first</em>.
     */
    public List getEntries() {
    		Collections.sort( entries, new EntryScoreReverseComparator());
    		return entries;        
    }
    
    // -------------------------------------------------------------------------
    /**
     * Empty the old/new lists.
     */
    public void clear() {

        entries.clear();

    }

    // -------------------------------------------------------------------------
    /**
     * Return a String representation of the envelope; only the number of source/target items is
     * displayed.
     */
    public String toString() {

        StringBuffer buf = new StringBuffer();
        //        buf.append('[');
        //
        //        buf.append("sources=");
        //        buf.append(sourceObjects.size());
        //        buf.append(", ");
        //
        //        buf.append("targets=");
        //        buf.append(targetObjects.size());
        //
        //        buf.append(']');
        Iterator it = entries.iterator();
        while (it.hasNext()) {
            Entry e = (Entry) it.next();
            buf.append(e.toString() + "\n");
        }
        return buf.toString();
    }

    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------

}