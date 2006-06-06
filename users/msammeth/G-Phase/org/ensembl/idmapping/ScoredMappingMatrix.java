/*
 *  
 */
package org.ensembl.idmapping;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

import org.ensembl.util.Util;

/**
 * @author arne
 * 
 * Objects of this class should be able to take pairs of exons and their associated score. Access
 * functions should enable to get sorted lists and pairs out by various criteria. This will be used
 * to store Exon and Transcript scores.
 */
public class ScoredMappingMatrix implements Serializable {

    //private TreeSet sourceTree, targetTree;

  private static final long serialVersionUID = 1L;

	private Map sourceMap, targetMap;

    private HashMap combinedMap;

    public static void main(String[] args) {

    }

    public ScoredMappingMatrix() {

        sourceMap = new HashMap();
        targetMap = new HashMap();
        combinedMap = new HashMap();

    }

    /** 
     * Build a scored mapping matrix from a list of entries.
     *
     */
    public ScoredMappingMatrix(List entries) {

       this();
       Iterator it = entries.iterator();
       while (it.hasNext()) {
           
           Entry e = (Entry)it.next();

           combinedMap.put(new CombinedKey(e), e);

           Util.addToMapList(sourceMap, new Long(e.getSource()), e);
           Util.addToMapList(targetMap, new Long(e.getTarget()), e);

       }

    }
    
    /**
     * Use an Entry struct to contain all scored pair information.
     * 
     * @param source
     * @param target
     * @return an Entry object with source,target and score
     */
    public Entry getEntry(long source, long target) {

        CombinedKey key = new CombinedKey(source, target);
        if (combinedMap.containsKey(key)) {
            return (Entry) combinedMap.get(key);
        } else {
            return null;
        }
    }

    /**
     * Find if there was a score defined between given source and target object.
     * 
     * @param source
     * @param target
     * @return whether there was a score defined between given source and target object.
     */
    public boolean hasScore(long source, long target) {

        return (getEntry(source, target) != null);
    }

    /**
     * Return a score between the given objects or 0.0f if there is no score.
     * 
     * @param source
     * @param target
     * @return score between the given objects or 0.0f if there is no score.
     */
    public float getScore(long source, long target) {

        Entry entry = getEntry(source, target);
        if (entry != null) {
            return entry.score;
        } else {
            return 0.0f;
        }
    }

    /**
     * Put a score to the two given objects.
     * 
     * @param source
     * @param target
     * @param score
     */
    public void addScore(long source, long target, float score) {

        Entry newEntry = new Entry(source, target, score);
        CombinedKey key = new CombinedKey(source, target);

        if( combinedMap.containsKey( key )) {
        		Entry oldEntry = (Entry) combinedMap.get( key );
        		oldEntry.score = newEntry.score;
        } else {
    			combinedMap.put(key, newEntry);
    			Util.addToMapList(sourceMap, new Long(source), newEntry);
    			Util.addToMapList(targetMap, new Long(target), newEntry);
        }
    }

    /**
     * Gives back a List of Entry objects that have that source. If no Entry objects have that
     * source, an empty list (rather than null) is returned.
     */
    public List sourceEntries(long source) {

        List list = (List) sourceMap.get(new Long(source));

        return list != null ? list : new ArrayList();

    }

    /**
     * Gives back a List of Entry objects that have that target. If no Entry objects have that
     * target, an empty list (rather than null) is returned.
     */
    public List targetEntries(long target) {

        List list = (List) targetMap.get(new Long(target));

        return list != null ? list : new ArrayList();

    }

    /**
     * Get the targets that have a score with a given source. Note this only returns the target IDs,
     * not Entry objects, so no scores are returned.
     */
    public long[] getTargetsForSource(long source) {

        List sourceEntries = sourceEntries(source);
        long[] result = new long[sourceEntries.size()];
        Iterator it = sourceEntries.iterator();
        int i = 0;
        while (it.hasNext()) {
            result[i++] = ((Entry) it.next()).getTarget();
        }

        return result;

    }

    /**
     * Get the sources that have a score with a given target. Note this only returns the source IDs,
     * not Entry objects, so no scores are returned.
     */
    public long[] getSourcesForTarget(long target) {

        List targetEntries = targetEntries(target);
        long[] result = new long[targetEntries.size()];
        Iterator it = targetEntries.iterator();
        int i = 0;
        while (it.hasNext()) {
            result[i++] = ((Entry) it.next()).getSource();
        }

        return result;

    }

    /**
     * @return A list of all the sources which have entries in this matrix.
     */
    public long[] getAllSources() {

        return setToLongArray(sourceMap.keySet());

    }

    /**
     * @return A list of all the sources which have entries in this matrix.
     */
    public long[] getAllTargets() {

        return setToLongArray(targetMap.keySet());

    }

    private long[] setToLongArray(Set s) {

        long[] result = new long[s.size()];

        Iterator it = s.iterator();
        int i = 0;
        while (it.hasNext()) {
            result[i++] = ((Long) it.next()).longValue();
        }

        return result;

    }

    public void remove(long source, long target) {

        CombinedKey key = new CombinedKey(source, target);
        combinedMap.remove(key);
        sourceMap.remove(new Long(source));
        targetMap.remove(new Long(target));
        
    }

    public int getEntryCount() {

        return combinedMap.size();

    }

    /**
     * Get the minimum and maximum scores for this matrix.
     * 
     * @return A 2-element array, the first element being the lowest score, the second the highest.
     */
    public float[] getMinMaxScores() {

        float[] result = {Float.MAX_VALUE, Float.MIN_VALUE};
        Collection values = combinedMap.values();
        Iterator it = values.iterator();
        while (it.hasNext()) {
            Entry e = (Entry) it.next();
            result[0] = Math.min(result[0], e.getScore());
            result[1] = Math.max(result[1], e.getScore());
        }

        return result;

    }

    /**
     * Get the average scores for this matrix.
     * 
     * @return The average score.
     */
    public float getAverageScore() {

        float total = 0.0f;
        Collection values = combinedMap.values();
        Iterator it = values.iterator();
        while (it.hasNext()) {
            Entry e = (Entry) it.next();
            total += e.getScore();
        }

        return total / values.size();

    }
    /**
     * Get a list of all the entries in this matrix.
     *  
     */
    public List getAllEntries() {

        return new ArrayList(combinedMap.values());

    }

    public void dump() {

        Set entries = combinedMap.entrySet();
        Iterator it = entries.iterator();
        while (it.hasNext()) {

            Map.Entry e = (Map.Entry) it.next();
            CombinedKey key = (CombinedKey) e.getKey();
            Entry entry = (Entry) e.getValue();
            System.out.println("Key: " + key.source + "," + key.target + " Entry: " + entry.getSource() + "," + entry.getTarget()
                    + "," + entry.getScore());
        }
    }

    public void dumpToFile(String rootDir, String outputFileName) {

        try {

            Set entries = combinedMap.entrySet();
            Iterator it = entries.iterator();

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator + outputFileName));

            while (it.hasNext()) {

                Map.Entry e = (Map.Entry) it.next();
                CombinedKey key = (CombinedKey) e.getKey();
                Entry entry = (Entry) e.getValue();
                writer.write(entry.getSource() + "\t" + entry.getTarget() + "\t" + entry.getScore() + "\n");

            }

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    //  -------------------------------------------------------------------------

    /**
     * Find the "envelopes" from a mapping matrix. An envelope is a set of interconnected source and
     * target objects.
     * 
     * The algorithm used to construct each envelope is: start with the original object, find all
     * target objects it scores with, then find all the source objects <em>that</em> shares
     * objects with and repeat the process for any new source objects. The process terminates when
     * all objects have been added.
     * 
     * @return A List of Envelope objects.
     */
    public List buildEnvelopes() {

        List envelopes = new ArrayList();

        Set allSourcesDone = new HashSet(); // cumulative list of processed sources

        // iterate over _all_ sources, each time we find one that isn't done,
        // we build an envelope from the "seed" source transcript
        long[] allSources = getAllSources();

        for (int i = 0; i < allSources.length; i++) {

            if (allSourcesDone.contains(new Long(allSources[i]))) {
                continue;
            }

            // start building a new envelope
            Set sourcesDone = new HashSet(); // fully processed sources for this envelope
            Set targetsDone = new HashSet(); // fully processed targets for this envelope

            Stack sourcesToProcess = new Stack();
            sourcesToProcess.push(new Long(allSources[i])); // add "seed" source

            Envelope env = new Envelope();
            
            while (!sourcesToProcess.empty()) {

                Long source = (Long) sourcesToProcess.pop();
                
                // ensure that the same source is not processed more than once
                if (sourcesDone.contains(source)) {
                    continue;
                }
                
                sourcesDone.add(source);
                // find all targets for this source
                long[] targetsToProcess = getTargetsForSource(source.longValue());
                for (int t = 0; t < targetsToProcess.length; t++) {

                    env.addEntry(getEntry(source.longValue(), targetsToProcess[t]));

                    Long target = new Long(targetsToProcess[t]);
                    if (targetsDone.contains(target)) {
                        continue;
                    }
                    
                    targetsDone.add(target);
                    
                    // find all sources for this target and add them the list of sources to process
                    long[] tmpSourcesToProcess = getSourcesForTarget(targetsToProcess[t]);
                    for (int s = 0; s < tmpSourcesToProcess.length; s++) {

                        Long tmpSource = new Long(tmpSourcesToProcess[s]);
                        if (!sourcesDone.contains(tmpSource)) {
                            sourcesToProcess.push(tmpSource);
                        }

                    }
                }

            }

            // add the sources for this envelope to the cumulative list
            allSourcesDone.addAll(sourcesDone);

            envelopes.add(env);

            /*
             * if (sourcesDone.size() == 1 && targetsDone.size() == 1) { Long[] sA = (Long[])
             * sourcesDone.toArray(new Long[sourcesDone.size()]); Long[] tA = (Long[])
             * targetsDone.toArray(new Long[targetsDone.size()]); for (int j = 0; j < sA.length;
             * j++) { for (int k = 0; k < tA.length; k++) { System.out.println("source: " + sA[j] + "
             * target: " + tA[k]); } } }
             */

        } // for i in allSources

        System.out.println("Total number of envelopes: " + envelopes.size());
        System.out.print( "Total content of entries: ");
        Iterator i;
        int noEntries = 0;
        i= envelopes.iterator();
        while( i.hasNext()) {
        		noEntries += ((Envelope)i.next()).size();
        }
        System.out.println( noEntries);
        return envelopes;

    }

    // -------------------------------------------------------------------------
    /**
     * Produce a list of source/target mappings based on the scores in this matrix.
     * Does the mapping based <em>solely</em> on scores, no disambiguating is done.
     */
    public List doMapping() {
        
        List mappings = new ArrayList();
        
        Map sourcesDone = new HashMap();
        Map targetsDone = new HashMap();
        
        List envelopes = buildEnvelopes();
        
        Iterator envIt = envelopes.iterator();
        while (envIt.hasNext()) {
            
            Envelope env = (Envelope)envIt.next();
            List entries = env.getEntries();
            Iterator entryIt = entries.iterator();
            while (entryIt.hasNext()) {
                
                Entry e = (Entry)entryIt.next();
                Long sourceID = new Long(e.getSource());
                Long targetID = new Long(e.getTarget());
                if (!sourcesDone.containsKey(sourceID) && !targetsDone.containsKey(targetID)) {
                    
                    mappings.add(e);
                    sourcesDone.put(sourceID, sourceID);
                    targetsDone.put(targetID, targetID);
                    
                }
                
            }
        
        }
        
        return mappings;
        
    }
    
    // -------------------------------------------------------------------------
    /**
     * Combine this matrix with another to get the union. If a source/target pair has different
     * scores in each, the highest score is used.
     */
    public void combineWith(ScoredMappingMatrix smm) {

        List smmEntries = smm.getAllEntries();
        Iterator it = smmEntries.iterator();
        while (it.hasNext()) {
            Entry e = (Entry) it.next();
            Entry myEntry = getEntry(e.getSource(), e.getTarget());

            // if it's not in our matrix, add it
            if( myEntry == null ) {
            		addScore( e.getSource(), e.getTarget(), e.getScore());
            } else {
                // if it is in our matrix, set the score if required
                if (e.getScore() > myEntry.getScore()) {
                	  myEntry.score = e.score;
                }
            }
        }

    }

    // -------------------------------------------------------------------------
    
    public String toString() {
        
        float[] minMax = getMinMaxScores();
        return ("ScoredMappingMatrix: Size: " + getEntryCount() + " Min score: " + minMax[0] + " Max score: " + minMax[1] + " Average score: " + getAverageScore());
    
    }
    
    // -------------------------------------------------------------------------

}

// -------------------------------------------------------------------------

/**
 * Object containing a source and target object that can
 */

class CombinedKey implements Serializable {

  private static final long serialVersionUID = 1L;
	
	public long source, target;

    public CombinedKey(long source, long target) {

        this.source = source;
        this.target = target;
    }

    public CombinedKey(Entry e) {

        this.source = e.source;
        this.target = e.target;
        
    }
    
    // need to override equals and hashcode to allow sensible comparisons
    public boolean equals(Object o) {

        CombinedKey ck = (CombinedKey) o;

        return (ck.source == source && ck.target == target);

    }

    public int hashCode() {

        int hash = 17;

        hash = 37 * hash * (int) source;

        hash = 37 * hash * (int) target;

        return hash;
    }

}
// -------------------------------------------------------------------------

