/*
 * Created on 08-Mar-2004
 * 
 * To change the template for this generated file go to Window>Preferences>Java>Code Generation>Code and Comments
 */
package org.ensembl.idmapping;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.CoordinateSystem;
import org.ensembl.datamodel.Exon;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.SequenceRegion;
import org.ensembl.driver.Driver;
import org.ensembl.driver.LocationConverter;

/**
 * @author arne
 * 
 * To change the template for this generated type comment go to Window>Preferences>Java>Code Generation>Code and Comments
 */
public class ExonDirectMapper {

    private boolean debug = true;
    
    public static void main(String[] args) {

    }

    private Config conf;

    private Cache cache;

    private ScoredMappingMatrix exonScoringMatrix;

    private CoordinateSystem[] commonCoordSystems; 

    private int ees = 0;

    public ExonDirectMapper(Config conf, Cache cache) {

        this.conf = conf;
        this.cache = cache;
        exonScoringMatrix = new ScoredMappingMatrix();
        
        Map exons = cache.getSourceExonsByInternalID();
        //Exon exon = (Exon)exons.get(new Long(192948));
        //System.out.println("In edm constructor: " + exon.getLocation().toString());
    }

    /**
     * Given two coordinate systems, which have the same name and version, find out wether the contained SeqRegions are the same. At
     * least 50% should have the same name, otherwise they are not really usefull to compare things.
     *  
     */
    private boolean seqRegionsCompatible(CoordinateSystem source, CoordinateSystem target) throws Exception {

        // deal with "chunked" genomes where seq_region name and length being equal does not indicate they are the same
        if (source.getName().equals("chunk")) {
            return false;
        }

        HashMap sourceRegions = new HashMap();
        int equalCount = 0;

        SequenceRegion[] sRegs = conf.getSourceDriver().getSequenceRegionAdaptor().fetchAllByCoordinateSystem(source);
        for (int i = 0; i < sRegs.length; i++) {
            sourceRegions.put(sRegs[i].getName(), new Long(sRegs[i].getLength()));
        }

        SequenceRegion[] tRegs = conf.getTargetDriver().getSequenceRegionAdaptor().fetchAllByCoordinateSystem(target);
        for (int j = 0; j < tRegs.length; j++) {
            if (sourceRegions.containsKey(tRegs[j].getName())) {
                equalCount++;
                if (((Long) sourceRegions.get(tRegs[j].getName())).longValue() != tRegs[j].getLength()) {
                    // same name different length is ok criteria for uncompatible seqregions
                    return false;
                }
            }
        }
        if (((equalCount / (double) sRegs.length) > 0.5) && ((equalCount / (double) tRegs.length) > 0.5)) {
            return true;
        } else {
            System.out.println("Only " + equalCount + " equally named seqRegions - won't do exon overlap mapping");
            return false;
        }
    }

    /**
     * Finds out if source and target exons can map to each other. Assume that the highest ranking common coordinate system will do.
     */
    public boolean mappedComparable() throws Exception {

        Driver sourceDriver = conf.getSourceDriver();
        Driver targetDriver = conf.getTargetDriver();

        CoordinateSystem[] sCoords = sourceDriver.getCoordinateSystemAdaptor().fetchAll();
        CoordinateSystem[] tCoords = targetDriver.getCoordinateSystemAdaptor().fetchAll();

        List commonCoordSystemList = new ArrayList();

        for (int i = 0; i < sCoords.length; i++) {
            for (int j = 0; j < tCoords.length; j++) {
                if (sCoords[i].getName().equals(tCoords[j].getName())) {
                    if (sCoords[i].getVersion() != null) {
                        if (!sCoords[i].getVersion().equals(tCoords[j].getVersion())) {
                            continue;
                        }
                    }
                    if (seqRegionsCompatible(sCoords[i], tCoords[j])) {
                        commonCoordSystemList.add(sCoords[i]);
                    }
                }
            }
        }

        commonCoordSystems = (CoordinateSystem[])commonCoordSystemList.toArray(new CoordinateSystem[commonCoordSystemList.size()]);
        
        return commonCoordSystems.length > 0;
    }

    /**
     * Connects to source and target database and extracts exons. Converts them to common coordinate system. (Should be established
     * before calling this). Creates scoring matrix for overlaps.
     * 
     * @throws Exception
     */
    public void buildOverlapScoring() throws Exception {

        if (commonCoordSystems.length == 0) {
            return;
        }

        debug("Reading exons");
        Collection sourceExons = readExons(conf.getSourceDriver(), cache.getSourceExonsByInternalID().values());
        Collection targetExons = readExons(conf.getTargetDriver(), cache.getTargetExonsByInternalID().values());
        debug("Doing overlap scoring");

        HashSet sourceOverlap, targetOverlap;
        sourceOverlap = new HashSet();
        targetOverlap = new HashSet();
        ExonContainerComparator ekc = new ExonContainerComparator();

        Iterator sE = sourceExons.iterator();
        Iterator tE = targetExons.iterator();
        ExonSortContainer topSource, topTarget;
        topSource = (ExonSortContainer) sE.next();
        topTarget = (ExonSortContainer) tE.next();

        // iterate through source and target and find potetially
        // overlapping Exons
        while (topSource != null || topTarget != null) {
            boolean addSource = false;
            boolean addTarget = false;

            if (topSource != null && topTarget != null) {
                int cmp = ekc.compare(topSource, topTarget);
                if (cmp <= 0) {
                    addSource = true;
                }
                if (cmp >= 0) {
                    addTarget = true;
                }
            } else if (topSource != null) {
                addTarget = true;
            } else {
                addSource = true;
            }

            if (addSource) {
                if (sourceOverlap.contains(topSource.getExon())) {
                    sourceOverlap.remove(topSource.getExon());
                } else {
                    sourceOverlap.add(topSource.getExon());
                    scoreExon(topSource.exon, targetOverlap);
                }
                if (sE.hasNext()) {
                    topSource = (ExonSortContainer) sE.next();
                } else {
                    topSource = null;
                }
            }

            if (addTarget) {
                if (targetOverlap.contains(topTarget.getExon())) {
                    targetOverlap.remove(topTarget.getExon());
                } else {
                    targetOverlap.add(topTarget.getExon());
                    scoreExon(sourceOverlap, topTarget.getExon());
                }
                if (tE.hasNext()) {
                    topTarget = (ExonSortContainer) tE.next();
                } else {
                    topTarget = null;
                }
            }

        }
        debug("Finished overlap scoring");

        //System.out.println("Source exons: " + sourceExons.size() + " target exons: " + targetExons.size() + " exonExonScore
        // called: " + ees);
    }

    /**
     * Calculates overlap score between two exons. Its done by dividing overlap region by exons sizes. 1.0 is full overlap on both
     * exons.
     */
    private float exonExonScore(Exon a, Exon b) {

        float result = 0.0f;
        int overlap = 0;
        // find overlap from Exon a and Exon b
        Location node = a.getLocation();
        HashMap aNodes = new HashMap();

        String id;

        while (node != null) {

            if (node.getStrand() == -1) {
                id = "- " + node.getSeqRegionName();
            } else {
                id = node.getSeqRegionName();
            }
            aNodes.put(id, node);
            node = node.next();
        }

        node = b.getLocation();
        while (node != null) {

            if (node.getStrand() == -1) {
                id = "- " + node.getSeqRegionName();
            } else {
                id = node.getSeqRegionName();
            }
            if (aNodes.containsKey(id)) {
                Location aNode = (Location) aNodes.get(id);
                int start = node.getStart() > aNode.getStart() ? node.getStart() : aNode.getStart();
                int end = node.getEnd() < aNode.getEnd() ? node.getEnd() : aNode.getEnd();

                if (end >= start) {
                    overlap += (end - start + 1);
                }
            }
            node = node.next();
        }

        result = (overlap / (float) a.getLocation().getLength()) + (overlap / (float) b.getLocation().getLength());
        result /= 2;

        // if phases are different penalise by 10%
        if (a.getPhase() != b.getPhase()) {
            result *= 0.9f;
        }

        ees++;

        return result;
    }

    /**
     * score given source Exon against all target Exons in given Set Put the scores into the exonScoringMatrix.
     * 
     * @param sourceExons
     * @param targetExon
     */
    private void scoreExon(HashSet sourceExons, Exon targetExon) {

        //System.out.println("### sourceExons.size: " + sourceExons.size());
        Exon sE;

        Iterator i = sourceExons.iterator();
        while (i.hasNext()) {
            sE = (Exon) i.next();
            if (exonScoringMatrix.hasScore(sE.getInternalID(), targetExon.getInternalID())) {
                continue;
            }

            float score = exonExonScore(sE, targetExon);
            if (score < 0.5f) {
                continue;
            }
            exonScoringMatrix.addScore(sE.getInternalID(), targetExon.getInternalID(), score);
        }
    }

    /**
     * score given target Exon against all source Exons in given Set Put the scores into the exonScoringMatrix.
     * 
     * @param sourceExon
     * @param targetExons
     */
    private void scoreExon(Exon sourceExon, HashSet targetExons) {

        Exon tE;

        Iterator i = targetExons.iterator();
        while (i.hasNext()) {
            tE = (Exon) i.next();
            if (exonScoringMatrix.hasScore(sourceExon.getInternalID(), tE.getInternalID())) {
                continue;
            }

            float score = exonExonScore(sourceExon, tE);
            if (score < 0.5f) {
                continue;
            }
            exonScoringMatrix.addScore(sourceExon.getInternalID(), tE.getInternalID(), score);
        }
    }

    /**
     * Read exons from one Driver and put them into a sorted Tree. If an exon location has more than one seq_region the exon will be
     * in the tree more than once...
     * 
     * @param driver
     * @return the map of exon locations or null if we cant have it
     */
    private Collection readExons(Driver driver, Collection inExons) throws Exception {

        if (commonCoordSystems.length == 0) {
            return null;
        }
        
        // "fill in" all common coordinate systems
        for (int i = 0; i < commonCoordSystems.length; i++) {
            CoordinateSystem c = driver.getCoordinateSystemAdaptor().fetch(commonCoordSystems[i].getName(), commonCoordSystems[i].getVersion());
            commonCoordSystems[i] = c;
        }
       
        // sort exons by location to improve assembly mapper cache performance
        debug("Sorting exons by location");
        List exons = new ArrayList(inExons);
        Collections.sort(exons, new ExonLocationComparator());

        // convert exon locations to common coordinate system
        System.out.println("Converting exon locations to common co-ordinate system");
        LocationConverter lc = driver.getLocationConverter();
        Iterator it = exons.iterator();
        while (it.hasNext()) {
            Exon exon = (Exon) it.next();
            // check each common coord system in turn to see if this exon converts to it 
            // (i.e. LocationConvertor doesn't return a gap)
            for (int i = 0; i < commonCoordSystems.length; i++) {
                Location convertedLocation = lc.convert(exon.getLocation(), commonCoordSystems[i]);
                if (convertedLocation.getSeqRegionName() != null) {
                    exon.setLocation(convertedLocation);
                    break;
                }
            }
        }

        // enter each exon once for start-1 and once for end, this allows
        // for 1 length exon bits
        debug("Building sorted list of exons");
        List exonContainers = new ArrayList();
        Iterator i = exons.iterator();

        while (i.hasNext()) {
            Exon e = (Exon) i.next();
            Location node = e.getLocation();
            while (node != null) {
                // insert exon into list
                if (node.getSeqRegionName() == null) {
                    node = lc.fetchComplete(node);
                    node = lc.convert(node, new CoordinateSystem("chromosome"));
                    System.out.println(node.getSeqRegionName());
                }
//                if (node.getSeqRegionName() == null) {
//                    System.out.println("### Null seq region name for exon " + e.getInternalID() + " loc "
//                            + e.getLocation().toString() + " sr = " + e.getLocation().getSequenceRegion());
//                } else {
//                    System.out.println("### Non-null seq region name for exon " + e.getInternalID() + " loc "
//                            + e.getLocation().toString());
//                }
                exonContainers.add(new ExonSortContainer(node.getSeqRegionName(), node.getStart() - 1, e));
                exonContainers.add(new ExonSortContainer(node.getSeqRegionName(), node.getEnd(), e));

                node = node.next();
            }
        }

        Collections.sort(exonContainers, new ExonContainerComparator());

        return exonContainers;
    }

    public ScoredMappingMatrix getScoringMatrix() {

        return exonScoringMatrix;
    }

    // -------------------------------------------------------------------------

    private void debug(String s) {

        if (debug) {
            System.out.println(s);
        }

    }

    /**
     * Helper class to store Exons by their position and make their position part of the hashkey they are stored under.
     * 
     * @author arne
     * 
     * To change the template for this generated type comment go to Window>Preferences>Java>Code Generation>Code and Comments
     */
    private class ExonSortContainer {

        private String seqRegionName;

        private long position;

        private Exon exon;

        public ExonSortContainer(String regionName, long pos, Exon e) {

            this.seqRegionName = regionName;
            this.position = pos;
            this.exon = e;
        }

        public Exon getExon() {

            return exon;

        }

    }

    /**
     * Helper class so ExonKeys are comparable. Needed to make a tree from the keys.
     * 
     * @author arne
     * 
     * To change the template for this generated type comment go to Window>Preferences>Java>Code Generation>Code and Comments
     */
    private class ExonContainerComparator implements Comparator {

        public int compare(Object a, Object b) {

            int strComp = ((ExonSortContainer) a).seqRegionName.compareTo(((ExonSortContainer) b).seqRegionName);
            if (strComp == 0) {
                if (((ExonSortContainer) a).position < ((ExonSortContainer) b).position) {
                    return -1;
                }
                if (((ExonSortContainer) a).position > ((ExonSortContainer) b).position) {
                    return 1;
                }
                return 0;
            } else {
                return strComp;
            }
        }
    }

}

// -------------------------------------------------------------------------

class ExonLocationComparator implements Comparator {

    public int compare(Object o1, Object o2) {

        Exon e1 = (Exon) o1;
        Exon e2 = (Exon) o2;
        return e1.getLocation().compareTo(e2.getLocation());

    }

}

// -------------------------------------------------------------------------
