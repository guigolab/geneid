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

import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Superclass for common functionality for transcript and gene score builders.
 */
public class ScoreBuilder {

    protected Cache cache;
    
    public ScoreBuilder(Cache cache) {
        
        this.cache = cache;
        
    }
    
    public ScoreBuilder() {

    }
    
    // -------------------------------------------------------------------------

    /**
     * Find the highest scoring entry in a list of entries, where the target exon is in a certain
     * list. Optionally remove highest-scoring exon from the Map.
     */
    protected Entry findHighestScoreTargetEntry(List entries, Map exons, boolean remove) {

        float max = -1.0f;

        Entry result = null;
	
	Long highestScoringExonID = new Long(-1);

        Iterator it = entries.iterator();
        while (it.hasNext()) {

            Entry entry = (Entry) it.next();
	    Long exonID = new Long(entry.getTarget());
            if (exons.containsKey(exonID) && entry.getScore() > max) {
                result = entry;
                max = entry.getScore();
		highestScoringExonID = exonID;
            }
        }

	if (remove && highestScoringExonID.longValue() > 0) {
	    exons.remove(highestScoringExonID);
	}

        return result;

    }

    // -------------------------------------------------------------------------

    /**
     * Find the highest scoring entry in a list of entries, where the source exon is in a certain
     * list.  Optionally remove highest-scoring exon from the Map.
     */
    protected Entry findHighestScoreSourceEntry(List entries, Map exons, boolean remove) {

        float max = -1.0f;

        Entry result = null;

	Long highestScoringExonID = new Long(-1);
        
	Iterator it = entries.iterator();
        while (it.hasNext()) {

            Entry entry = (Entry) it.next();
	    Long exonID = new Long(entry.getSource());
            if (exons.containsKey(exonID) && entry.getScore() > max) {
                result = entry;
                max = entry.getScore();
		highestScoringExonID = exonID;
            }
        }

	if (remove && highestScoringExonID.longValue() > 0) {
            exons.remove(highestScoringExonID);
        }

        return result;

    }
    
    // -------------------------------------------------------------------------
    
}
