/*
 * Copyright (C) 2004 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License along with this library; if not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package org.ensembl.idmapping;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Statement;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.ensembl.datamodel.Accessioned;
import org.ensembl.datamodel.Exon;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Sequence;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.Translation;
import org.ensembl.driver.plugin.standard.BaseAdaptor;
import org.ensembl.util.AscendingInternalIDComparator;
import org.ensembl.util.StringUtil;
import org.ensembl.util.Util;

/**
 * Map stable IDs based on internal ID mappings, and generate stable ID events.
 */

public class StableIDMapper {

    private String rootDir;

    // Map of StableIDEventContainer objects created this session
    private Map createdStableIDEvents = new HashMap();

    // list of StableIDEventContainer objects to be propagated
    private Map propagatedStableIDEvents = new HashMap();

    // Map of Lists of old new stableID mappings for debugging
    // Keyed on type (exon, transcript, translation, gene)
    private Map debugMappings = new HashMap();

    // mapping session ID used for this session
    private long currentMappingSessionID = 0;

    private Config conf;

    public StableIDMapper(String rootDir, Config conf) {

        this.rootDir = rootDir;
        this.conf = conf;

    }

    // -------------------------------------------------------------------------
    /**
     * Generate stable IDs for targets. Where a source has been mapped to a target, the stable ID of the target is the same as the
     * stable ID of the source. For new targets, a new stable ID is assigned. Stable ID events for mappings, and for stable ID
     * creation and deletion are produced.
     * 
     * @param sourcesByInternalID All the source objects.
     * @param targetsByInternalID All the target objects.
     * @param mappings A List of targets with their accession IDs and versions filled in.
     */
    public List mapStableIDs(Map sourcesByInternalID, Map targetsByInternalID, List mappings) {

        List allSources = new ArrayList(sourcesByInternalID.values());

        String type = getTypeFromObject(allSources.get(0));

        int[] mappedCount = {0, 0}; // known, novel
        int newCount = 0;
        int[] lostCount = {0, 0}; // known, novel

        // cache mappings
        Map sourcesMapped = new HashMap();
        Map targetsMapped = new HashMap();
        Iterator mIt = mappings.iterator();
        while (mIt.hasNext()) {
            Entry entry = (Entry) mIt.next();
            sourcesMapped.put(new Long(entry.getSource()), new Long(entry.getTarget()));
            targetsMapped.put(new Long(entry.getTarget()), new Long(entry.getSource()));
        }

        // assign existing or new stable IDs to TARGET exons

        // transfer stable ID from source->target for mapped objects
        List allTargets = new ArrayList(targetsByInternalID.values());
        Iterator ateit = allTargets.iterator();
        while (ateit.hasNext()) {

            Accessioned target = (Accessioned) ateit.next();
            Long targetID = new Long(target.getInternalID());
            if (targetsMapped.containsKey(targetID)) {

                Long sourceID = (Long) targetsMapped.get(targetID);
                Accessioned source = (Accessioned) sourcesByInternalID.get(sourceID);

                // set target's stable ID to be the same as the source's stable
                // ID
                target.setAccessionID(source.getAccessionID());

                // version calculation is different for exons, transcripts,
                // genes
                int version = calculateNewVersion(source, target);
                target.setVersion(version);

                // generate a "mapped" stable ID event
                if (!type.equals("exon")) {
                    addCreatedStableIDEvent(new StableIDEventContainer(source.getAccessionID(), source.getVersion(), target
                            .getAccessionID(), target.getVersion(), getTypeFromObject(target), currentMappingSessionID));
                }

                // add this mapping to the list for debugging purposes
                addDebugMapping(type, source.getInternalID(), target.getInternalID(), target.getAccessionID());

                // update the statistics
                mappedCount = updateCount(source, mappedCount);

            }

        }

        // assign new stable IDs to targets that weren't mapped
        // find the highest mapped stable ID, and use it as the base for
        // assigning new ones
        String newID = findHighestStableID(allTargets);

        ateit = allTargets.iterator();
        while (ateit.hasNext()) {

            Accessioned target = (Accessioned) ateit.next();
            Long targetID = new Long(target.getInternalID());
            if (!targetsMapped.containsKey(targetID)) {
                // not mapped - generate new stable ID
                newID = incrementStableID(newID);
                target.setAccessionID(newID);

                // new stable ID, so version is 1
                target.setVersion(1);

                // generate a "new" stable ID event
                if (!type.equals("exon")) {
                    addCreatedStableIDEvent(new StableIDEventContainer(null, 0, target.getAccessionID(), target.getVersion(),
                            getTypeFromObject(target), currentMappingSessionID));
                }

                // update the statistics
                newCount++;
            }
        }

        // find sources that weren't mapped, mark as lost
        Iterator asit = allSources.iterator();
        while (asit.hasNext()) {
            Accessioned source = (Accessioned) asit.next();

            if (!sourcesMapped.containsKey(new Long(source.getInternalID()))) {

                // generate a "lost" stable ID event
                if (!type.equals("exon")) {
                    addCreatedStableIDEvent(new StableIDEventContainer(source.getAccessionID(), source.getVersion(), null, 0,
                            getTypeFromObject(source), currentMappingSessionID));
                }

                // update the statistics
                lostCount = updateCount(source, lostCount);

            }
        }

        Collections.sort(allTargets, new AscendingInternalIDComparator());

        String fileName = rootDir + File.separator + type + "_stable_id.txt";
        dumpStableIDsToFile(allTargets, fileName);
        System.out.println("Wrote " + allTargets.size() + " " + type + "s to " + fileName);

        String stats = generateMappingStatistics(type, mappedCount, lostCount, newCount);
        System.out.println("\n" + stats);
        writeStringToFile(stats, rootDir + File.separator + type + "_mapping_statistics.txt");

        return allTargets;

    }

    // -------------------------------------------------------------------------

    public void dumpLostGeneAndTranscripts(Cache cache) {

        Map deletedStableIDs = new HashMap();
        Map deletedGeneIDs = new HashMap();
        Map deletedTranscriptIDs = new HashMap();
	
	// iterate over all stable ID events
	// find any where target = null (deleted) - store in deletedStableIDs
	
	Iterator sideit = createdStableIDEvents.values().iterator();
	while (sideit.hasNext()) {
	    
	    StableIDEventContainer sidec = (StableIDEventContainer) sideit.next();
	    if (sidec.getNewStableID() == null) {
		if (sidec.getType().equals("gene")) {
		    deletedGeneIDs.put(sidec.getOldStableID(), "");
		} else if (sidec.getType().equals("transcript")) {
		    deletedTranscriptIDs.put(sidec.getOldStableID(), "");
		}
	    }
	}
	
	// iterate again over stable ID events
	// mark as merged those that have other non-"deleted" events
	
	sideit = createdStableIDEvents.values().iterator();
	while (sideit.hasNext()) {
	    
	    StableIDEventContainer sidec = (StableIDEventContainer) sideit.next();
	    if (sidec.getNewStableID() != null) {
		if (sidec.getType().equals("gene")) {
		    if (deletedGeneIDs.containsKey(sidec.getOldStableID())) {
			deletedGeneIDs.put(sidec.getOldStableID(), sidec.getNewStableID());
		    }
		} else if (sidec.getType().equals("transcript")) {
		    if (deletedTranscriptIDs.containsKey(sidec.getOldStableID())) {
			deletedTranscriptIDs.put(sidec.getOldStableID(), sidec.getNewStableID());
		    }
		}
	    }
	}

	// now dump them to appropriately-named files
	try {
	    
	    OutputStreamWriter deletedGenesWriter = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator + "genes_lost_deleted.txt"));
	    OutputStreamWriter mergedGenesWriter = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator + "genes_lost_merged.txt"));
	    OutputStreamWriter deletedTranscriptsWriter = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator + "transcripts_lost_deleted.txt"));
	    OutputStreamWriter mergedTranscriptsWriter = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator + "transcripts_lost_merged.txt"));
	    
	    Iterator deletedGenesIterator = deletedGeneIDs.keySet().iterator();
	    while (deletedGenesIterator.hasNext()) {
		String oldGeneID = (String)deletedGenesIterator.next();
		String newGeneID = (String)deletedGeneIDs.get(oldGeneID);
		Gene oldGene = (Gene)cache.getSourceGenesByStableID().get(oldGeneID);
		String knownNovel = oldGene.isKnown() ? "KNOWN" : "NOVEL";
		if (newGeneID.equals("")) {   // deleted
		    deletedGenesWriter.write(oldGeneID + "\t" + knownNovel + "\n");
		} else {                      // lost
		    mergedGenesWriter.write(oldGeneID + "\t" + newGeneID + "\t" + knownNovel + "\n");
		}
	    }
	    
	    Iterator deletedTranscriptsIterator = deletedTranscriptIDs.keySet().iterator();
            while (deletedTranscriptsIterator.hasNext()) {
                String oldTranscriptID = (String)deletedTranscriptsIterator.next();
                String newTranscriptID = (String)deletedTranscriptIDs.get(oldTranscriptID);
                Transcript oldTranscript = (Transcript)cache.getSourceTranscriptsByStableID().get(oldTranscriptID);
                String knownNovel = oldTranscript.isKnown() ? "KNOWN" : "NOVEL";
                if (newTranscriptID.equals("")) {   // deleted
                    deletedTranscriptsWriter.write(oldTranscriptID + "\t" + knownNovel + "\n");
                } else {                            // lost
                    mergedTranscriptsWriter.write(oldTranscriptID + "\t" + newTranscriptID + "\t" + knownNovel + "\n");
                }
            }

	    deletedGenesWriter.close();
	    mergedGenesWriter.close();
	    deletedTranscriptsWriter.close();
	    mergedTranscriptsWriter.close();
	    
        } catch (IOException e) {
	    
            e.printStackTrace();
        }

    }

    //---------------------------------------------------------------------

    /**
     * Find the highest currently-assigned stable ID.
     * 
     * @param existing The current list of {exons|transcripts etc} with assigned stable IDs; objects in list should be of type
     *            Accessioned.
     * @return The highest currently-assigned stable ID.
     */
    private String findHighestStableID(List existing) {

        String maxStableID = "";

        Iterator it = existing.iterator();
        while (it.hasNext()) {
            String id = ((Accessioned) it.next()).getAccessionID();
            if (id != null && id.compareTo(maxStableID) > 0) {
                maxStableID = id;
            }
        }

        //System.out.println("Max stable ID = " + maxStableID);

        return maxStableID;

    }

    //---------------------------------------------------------------------
    /**
     * Increment by 1 a stable ID, regardless of its prefix. Numeric part of accession ID will have exactly 11 digits.
     */
    private String incrementStableID(String stableID) {

        int prefixLength = org.ensembl.util.StringUtil.indexOfFirstDigit(stableID);
        String prefix = stableID.substring(0, prefixLength);
        long number = Long.parseLong(stableID.substring(prefixLength, stableID.length()));
        number++;

        // Pad with 0 if necessary
        DecimalFormat min11Digits = new DecimalFormat();
        min11Digits.setMinimumIntegerDigits(11);
        // Stop ',' being inserted between digits!
        min11Digits.setGroupingSize(min11Digits.getMinimumIntegerDigits() + 1);

        return prefix + min11Digits.format(number);

    }

    // -------------------------------------------------------------------------

    /**
     * Calculate new version for target object based on source object.
     */
    private int calculateNewVersion(Accessioned source, Accessioned target) {

        int version = source.getVersion();

        if (source instanceof Exon && target instanceof Exon) {

            // EXONS - increment version if sequence changed, otherwise keep the
            // same
            if (!((Exon) source).getSequence().getString().equals(((Exon) target).getSequence().getString())) {
                version++;
            }

        } else if ((source instanceof Transcript && target instanceof Transcript)) {

            // TRANSCRIPTS - if spliced sequence of exons changed, increment
            // version
            Sequence sourceSequence = ((Transcript) source).getSequence();
            Sequence targetSequence = ((Transcript) target).getSequence();
            if (sourceSequence != null && targetSequence != null) {
                if (!sourceSequence.getString().equals(targetSequence.getString())) {
                    version++;
                }
            }

        } else if (source instanceof Translation && target instanceof Translation) {

            // TRANSLATIONS - increment version if transcripts changed
            Transcript sourceTranscript = ((Translation) source).getTranscript();
            Transcript targetTranscript = ((Translation) target).getTranscript();
            Sequence sourceSequence = sourceTranscript.getSequence();
            Sequence targetSequence = targetTranscript.getSequence();
            if (sourceSequence != null && targetSequence != null) {
                if (!sourceSequence.getString().equals(targetSequence.getString())) {
                    version++;
                }
            }

        } else if (source instanceof Gene && target instanceof Gene) {

            // GENES - increment version if any transcript changes
            Set sourceTranscriptDescriptions = extractTranscriptAccessionAndVersions((Gene) source);
            Set targetTranscriptDescriptions = extractTranscriptAccessionAndVersions((Gene) target);
            if (!sourceTranscriptDescriptions.equals(targetTranscriptDescriptions)) {
                version++;
            }

        } else {
            System.err.println("Can't calculate version information for objects of type " + source.getClass().getName() + " "
                    + target.getClass().getName());
        }

        return version;

    }

    // -------------------------------------------------------------------------
    /**
     * Creates a set of strings where each string = transcript.accessionID + transcript.version.
     * 
     * @param gene The gene from which to extract the transcripts.
     * @return set of strings representing each transcript.
     */
    private Set extractTranscriptAccessionAndVersions(Gene gene) {

        Set descriptions = new HashSet();
        for (Iterator it = gene.getTranscripts().iterator(); it.hasNext();) {
            Transcript transcript = (Transcript) it.next();
            descriptions.add(transcript.getAccessionID() + transcript.getVersion());
        }

        return descriptions;

    }

    // -------------------------------------------------------------------------

    private String getTypeFromObject(Object o) {

        String type = "";

        if (o instanceof Exon) {

            type = "exon";

        } else if (o instanceof Transcript) {

            type = "transcript";

        } else if (o instanceof Translation) {

            type = "translation";

        } else if (o instanceof Gene) {

            type = "gene";

        } else {

            System.err.println("Cannot get type for " + o.toString());

        }

        return type;

    }

    // -------------------------------------------------------------------------

    public void dumpStableIDsToFile(List objects, String outputFileName) {

        try {

            Iterator it = objects.iterator();

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(outputFileName));

            while (it.hasNext()) {

                Accessioned obj = (Accessioned) it.next();
                writer.write(obj.getInternalID() + "\t" + obj.getAccessionID() + "\t" + obj.getVersion() + "\n");

            }

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }

    }

    //  -------------------------------------------------------------------------

    public void writeStringToFile(String str, String outputFileName) {

        try {

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(outputFileName));
            writer.write(str);

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    // -------------------------------------------------------------------------
    /**
     * Generate statistics for the mapping.
     * 
     * @param Type the type of object being considered (gene, transcript, translation, exon)
     * @param mappedCount Number of known (1st array entry) and novel (2nd array entry) objects.
     * @param lostCount Number of known (1st array entry) and novel (2nd array entry) objects.
     * @param newCount Number of new objects.
     * @return A formatted String of results, or an empty string for exons.
     */
    private String generateMappingStatistics(String type, int[] mappedCount, int[] lostCount, int newCount) {

        NumberFormat df = new DecimalFormat("##.##%");
        StringBuffer result = new StringBuffer();

        result.append(StringUtil.capitaliseFirstLetter(type) + " Mapping Statistics for "
                + System.getProperty("idmapping.source.database") + " -> " + System.getProperty("idmapping.target.database")
                + "\n\n");
        result.append("Type\tMapped\tLost\tPercentage\n");
        result.append("-----------------------------------\n");

        int totalMapped;
        int totalLost;

        // Exons aren't categorised into known/novel
        if (type.equalsIgnoreCase("exon")) {

            totalMapped = mappedCount[0];
            totalLost = lostCount[0];

        } else {
            // Genes, transcripts, translations can be known or novel
            for (int i = 0; i < 2; i++) {

                String s = (i == 0 ? "Known" : "Novel");
                float ratio = (float) mappedCount[i] / (float) (mappedCount[i] + lostCount[i]);
                result.append(s + "\t" + mappedCount[i] + "\t" + lostCount[i] + "\t" + df.format(ratio) + "\n");

            }

            totalMapped = mappedCount[0] + mappedCount[1];
            totalLost = lostCount[0] + lostCount[1];
        }

        float ratio = (float) totalMapped / (float) (totalMapped + totalLost);
        result.append("Total\t" + totalMapped + "\t" + totalLost + "\t" + df.format(ratio) + "\n");

        return result.toString();

    }

    // -------------------------------------------------------------------------

    private int[] updateCount(Accessioned object, int[] count) {

        int[] result = {count[0], count[1]};

        if (object instanceof Gene) {
            if (((Gene) object).isKnown()) {
                result[0]++;
            } else {
                result[1]++;
            }
        } else if (object instanceof Transcript) {
            if (((Transcript) object).isKnown()) {
                result[0]++;
            } else {
                result[1]++;
            }
        } else if (object instanceof Translation) {
            if (((Translation) object).isKnown()) {
                result[0]++;
            } else {
                result[1]++;
            }
        } else if (object instanceof Exon) {
            // Exons aren't known/novel, just use first array entry
            result[0]++;
        }

        return result;

    }

    //  -------------------------------------------------------------------------
    /**
     * Upload stable IDs and versions from a list of objects to a database.
     * 
     * @param con The database connection to use.
     * @param objects The list of Accessioned objects to get the stable IDs and versions from.
     */
    public void uploadStableIDs(Connection con, List objects) {

        String type = getTypeFromObject(objects.get(0));
        String table = type + "_stable_id";

        System.out.println("Uploading " + objects.size() + " " + type + " stable IDs to " + table);

        try {

            con.setAutoCommit(false); // required for batch updates
            PreparedStatement stmt = con.prepareStatement("INSERT INTO " + table + " VALUES(?, ?, ?)");
            Iterator it = objects.iterator();
            while (it.hasNext()) {
                Accessioned obj = (Accessioned) it.next();
                stmt.setLong(1, obj.getInternalID());
                stmt.setString(2, obj.getAccessionID());
                stmt.setInt(3, obj.getVersion());
                stmt.addBatch();
            }

            int[] numUpdates = stmt.executeBatch();
            con.commit();
            con.setAutoCommit(true);
            stmt.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    // -------------------------------------------------------------------------
    /**
     * Generate stable ID events for genes/transcripts that have been merged or split. Note that other stable ID events (for direct
     * mapping events, new and lost events) are created as part of the mapping process.
     */
    public void generateMergedSplitEvents(Cache cache, List exonMappings, List transcriptMappings, List geneMappings) {

        // for each source exon, check if it is part of a transcript that maps to the
        // transcript that the target exon maps to
        Set sourceExons = cache.getSourceExonsByInternalID().keySet();
        Iterator seit = sourceExons.iterator();
        while (seit.hasNext()) {

            Long sourceExonID = (Long) seit.next();
            List sourceTranscripts = (List) cache.getSourceTranscriptsByExonInternalID().get(sourceExonID);
            Long targetExonID = (Long) cache.getExonMappingsMap().get(sourceExonID);

            if (targetExonID != null) { // may not have mapped

                List targetTranscripts = (List) cache.getTargetTranscriptsByExonInternalID().get(targetExonID);

                Iterator ttit = targetTranscripts.iterator();

                while (ttit.hasNext()) {

                    Transcript targetTranscript = (Transcript) ttit.next();

                    Iterator stit = sourceTranscripts.iterator();

                    while (stit.hasNext()) {

                        Transcript sourceTranscript = (Transcript) stit.next();

                        // generate event
                        StableIDEventContainer sidec = new StableIDEventContainer(sourceTranscript.getAccessionID(),
                                sourceTranscript.getVersion(), targetTranscript.getAccessionID(), targetTranscript.getVersion(),
                                "transcript", currentMappingSessionID);
                        addCreatedStableIDEvent(sidec);
                    }

                }

                // now for gene
                Gene sourceGene = (Gene) cache.getSourceGeneByExonInternalID().get(sourceExonID);
                Gene targetGene = (Gene) cache.getTargetGeneByExonInternalID().get(targetExonID);

                // generate event
                StableIDEventContainer sidec = new StableIDEventContainer(sourceGene.getAccessionID(), sourceGene.getVersion(),
                        targetGene.getAccessionID(), targetGene.getVersion(), "gene", currentMappingSessionID);
                addCreatedStableIDEvent(sidec);
            }
        }

    }

    // -------------------------------------------------------------------------

    private void addCreatedStableIDEvent(StableIDEventContainer sidec) {

        createdStableIDEvents.put(sidec.getKey(), sidec);

    }

    // -------------------------------------------------------------------------

    public void writeNewStableIDEvents() {

        // make sure mapping session ID is set correctly
        Iterator it = createdStableIDEvents.values().iterator();
        while (it.hasNext()) {
            StableIDEventContainer sidec = (StableIDEventContainer) it.next();
            sidec.setMappingSessionID(currentMappingSessionID);
        }

        writeStableIDEvents(new ArrayList(createdStableIDEvents.values()), getCreatedStableIDEventFileName());

    }

    // -------------------------------------------------------------------------

    public void writePropagatedStableIDEvents() {

        writeStableIDEvents(new ArrayList(propagatedStableIDEvents.values()), getPropagatedStableIDEventFileName());

    }

    // -------------------------------------------------------------------------

    /**
     * Write stable ID events to a file.
     */
    private void writeStableIDEvents(List events, String fileName) {

        try {

            OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));
            Iterator it = events.iterator();
            while (it.hasNext()) {

                StableIDEventContainer sidec = (StableIDEventContainer) it.next();
                writeStableIDEvent(writer, sidec);

            }

            writer.close();

        } catch (IOException e) {

            e.printStackTrace();
        }

    }

    // -------------------------------------------------------------------------

    private void writeStableIDEvent(OutputStreamWriter writer, StableIDEventContainer sidec) throws IOException {

        String oldSID = (sidec.getOldStableID() == null) ? "\\N" : sidec.getOldStableID();
        String newSID = (sidec.getNewStableID() == null) ? "\\N" : sidec.getNewStableID();

        writer.write(oldSID + "\t" + sidec.getOldVersion() + "\t" + newSID + "\t" + sidec.getNewVersion() + "\t"
                + sidec.getMappingSessionID() + "\t" + sidec.getType() + "\n");

    }

    // -------------------------------------------------------------------------
    /**
     * Load the stable ID events into a database.
     */
    public void uploadNewStableIDEvents() {

        Config.uploadFromFile(getCreatedStableIDEventFileName(), "stable_id_event", conf.getTargetConnection(), false);

    }

    // -------------------------------------------------------------------------

    private String getCreatedStableIDEventFileName() {

        return rootDir + File.separator + "stable_id_event_new.txt";

    }

    private String getPropagatedStableIDEventFileName() {

        return rootDir + File.separator + "stable_id_event_propagated.txt";

    }

    // -------------------------------------------------------------------------

    /**
     * Generate the mapping session ID being used (stored internally); initialise it if necessary. The mapping_session table, with a
     * row for the new session, is written to mapping_session.txt
     */
    public void generateNewMappingSessionID() {

        try {

            // get the current mapping session ID - max existing + 1
            Connection sourceCon = conf.getSourceConnection();

            String sql = "SELECT MAX(mapping_session_id) FROM mapping_session";

            ResultSet rs = BaseAdaptor.executeQuery(sourceCon, sql);
            if (rs.next()) {

                currentMappingSessionID = rs.getInt(1) + 1;

            } else {

                System.out.println("Cannot get maximum mapping_session_id - mapping_session table empty?");
                currentMappingSessionID = 1;

            }
            System.out.println("Mapping session ID for this session is " + currentMappingSessionID);
            rs.close();

            // write the existing mapping session data and the new one to
            try {

                OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(rootDir + File.separator
                        + "mapping_session.txt"));

                // existing
                rs = sourceCon.createStatement().executeQuery("SELECT * FROM mapping_session");
                while (rs.next()) {
                    writer.write(rs.getInt(1) + "\t" + rs.getString(2) + "\t" + rs.getString(3) + "\t" + rs.getString(4) + "\n");
                }
                rs.close();

                // new
                SimpleDateFormat sf = new SimpleDateFormat("yyyyMMddHHmmss");
                writer.write(currentMappingSessionID + "\t" + System.getProperty("idmapping.source.database") + "\t"
                        + System.getProperty("idmapping.target.database") + "\t" + sf.format(new Date()) + "\n");

                writer.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

        } catch (Exception e) {

            System.err.println("Error getting new mapping session");
            e.printStackTrace();
        }

        System.out.println("Updated mapping_session table written to mapping_session.txt");
    }

    // -------------------------------------------------------------------------
    /**
     * @return A List of StableIDEventContainers for changed objects, i.e. where the new version != the old version. Note that this
     *         is only for existing objects, i.e. new ones are not counted.
     */
    public List getChanged(String type) {

        List result = new ArrayList();

        Iterator it = createdStableIDEvents.values().iterator();
        while (it.hasNext()) {
            StableIDEventContainer sidec = (StableIDEventContainer) it.next();
            if (sidec.getType().equalsIgnoreCase(type)) {
                if (sidec.getOldStableID() != null) { // ignore created
                    if (sidec.getOldVersion() != sidec.getNewVersion()) {
                        result.add(sidec);
                    }
                }
            }
        }

        return result;

    }

    // -------------------------------------------------------------------------
    /**
     * Propagate stable ID events. Write directly to file to avoid having to store lots of them in memory
     */

    public void propagateStableIDEvents() {

        // Find mapping_session_id of ALL/LATEST mapping session.
        long allLatest = getAllLatest();
        if (allLatest < 1) {
            System.err.println("Warning - could not get ALL/LATEST, some stable_id_events will have mapping_session=0");
        }
        System.out.println("Current mapping session ID is " + currentMappingSessionID + "; ALL/LATEST is " + allLatest);
        // Create an event with mapping session of allLatest for each recently
        // generated event
        Iterator it = createdStableIDEvents.values().iterator();
        while (it.hasNext()) {
            StableIDEventContainer sidec = (StableIDEventContainer) it.next();
            StableIDEventContainer pe = new StableIDEventContainer(sidec);
            pe.setMappingSessionID(allLatest);
            propagatedStableIDEvents.put(pe.getKey(), pe);
            sidec = null;
            pe = null;
        }

        System.out.println("Total propagated stable ID events now " + propagatedStableIDEvents.size());

        // Create records for chained mappings; for each previous
        // mapping_session_id *X*,
        // where *X* is not allLatest or currentMappingSessionID, do the
        // following query and insert all the results into stable_id_event

        long[] previousMappingSessions = getPreviousMappingSessionIDs(allLatest, currentMappingSessionID);

        Connection con = conf.getTargetConnection();

        for (int i = 0; i < previousMappingSessions.length; i++) {

            System.out.println("Propagating stable ID events for previous mapping session with ID " + previousMappingSessions[i]
                    + " (" + (i + 1) + " of " + previousMappingSessions.length + ")");
            String sql = "SELECT s1.old_stable_id, s1.old_version, s2.new_stable_id, s2.new_version, " + allLatest + ", s1.type "
                    + "FROM stable_id_event s1, stable_id_event s2 WHERE s1.mapping_session_id = " + previousMappingSessions[i]
                    + " AND s2.mapping_session_id = " + allLatest + " "
                    + "AND s1.new_stable_id is not null AND s2.old_stable_id is not null "
                    + "AND s1.new_stable_id = s2.old_stable_id " + "AND s1.new_version = s2.old_version "
                    + "AND s1.old_stable_id is not null";

            try {

                Statement stmt = con.createStatement(java.sql.ResultSet.TYPE_FORWARD_ONLY, java.sql.ResultSet.CONCUR_READ_ONLY);
                stmt.setFetchSize(100);
                //System.out.println("Executing " + sql);
                ResultSet rs = stmt.executeQuery(sql);
                //System.out.println("Results:");
                while (rs.next()) {
                    //System.out.println(rs.getString(1));
                    StableIDEventContainer sidec = new StableIDEventContainer(rs.getString(1), rs.getInt(2), rs.getString(3), rs
                            .getInt(4), rs.getString(6), rs.getLong(5));
                    propagatedStableIDEvents.put(sidec.getKey(), sidec);
                }
                rs.close();
                stmt.close();

            } catch (Exception e) {
                e.printStackTrace();
            }

            System.out.println("Total propagated stable ID events now " + propagatedStableIDEvents.size());

        } // foreach previous mapping session

        System.out.println("Propagated a total of " + propagatedStableIDEvents.size() + " stable ID events");

    }

    // -------------------------------------------------------------------------

    private long getAllLatest() {

        long allLatest = -1;

        String sql = "SELECT mapping_session_id FROM mapping_session WHERE old_db_name='ALL' AND new_db_name='LATEST'";

        Connection con = conf.getSourceConnection();

        try {

            Statement stmt = con.createStatement();
            ResultSet rs = stmt.executeQuery(sql);
            if (rs.next()) {
                allLatest = rs.getLong(1);
                System.out.println("ALL/LATEST mapping_session_id is " + allLatest);
            } else {

                // TODO create if required
                System.out.println("No ALL/LATEST mapping session found - create???");
            }
            rs.close();
            stmt.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        return allLatest;

    }

    // -------------------------------------------------------------------------

    public void uploadPropagatedStableIDEvents() {

        Config.uploadFromFile(getPropagatedStableIDEventFileName(), "stable_id_event", conf.getTargetConnection(), true);

    }

    //  -------------------------------------------------------------------------

    public void uploadExistingStableIDEvents() {

        Config.uploadFromFile(rootDir + File.separator + "stable_id_event_existing.txt", "stable_id_event", conf
                .getTargetConnection(), false);

    }

    //  -------------------------------------------------------------------------

    public void uploadMappingSession() {

        Config.uploadFromFile(rootDir + File.separator + "mapping_session.txt", "mapping_session", conf.getTargetConnection(),
                false);

    }

    // -------------------------------------------------------------------------

    public long getCurrentMappingSessionID() {

        if (currentMappingSessionID == 0) {

            generateNewMappingSessionID();

        }

        return currentMappingSessionID;

    }

    // -------------------------------------------------------------------------

    private long[] getPreviousMappingSessionIDs(long allLatest, long currentMappingSessionID) {

        String sql = "SELECT mapping_session_id FROM mapping_session WHERE mapping_session_id NOT IN (" + allLatest + ", "
                + currentMappingSessionID + ")";

        Connection con = conf.getTargetConnection();
        List ids = new ArrayList();

        try {

            Statement stmt = con.createStatement();
            ResultSet rs = stmt.executeQuery(sql);
            while (rs.next()) {
                ids.add(new Long(rs.getLong(1)));
            }
            rs.close();
            stmt.close();

        } catch (Exception e) {
            e.printStackTrace();
        }

        // ugh why doesn't toArray work with long[}??
        long[] result = new long[ids.size()];
        Iterator it = ids.iterator();
        int i = 0;
        while (it.hasNext()) {
            Long ll = (Long) it.next();
            result[i++] = ll.longValue();
        }

        return result;

    }

    // -------------------------------------------------------------------------

    private void addDebugMapping(String type, long sourceInternalID, long targetInternalID, String stableID) {

        Util.addToMapList(debugMappings, type, new MappingContainer(sourceInternalID, targetInternalID, stableID));

    }

    // -------------------------------------------------------------------------

    public void dumpDebugMappingsToFile() {

        try {

            Iterator kit = debugMappings.keySet().iterator();

            while (kit.hasNext()) {

                String type = (String) kit.next();
                String fileName = rootDir + File.separator + "debug" + File.separator + type + "_mappings.txt";
                OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));

                List typeMappings = (List) (debugMappings.get(type));
                Iterator it = typeMappings.iterator();

                while (it.hasNext()) {

                    MappingContainer mc = (MappingContainer) it.next();
                    writer.write(mc.toString() + "\n");

                }

                writer.close();

                System.out.println("Wrote " + type + " mappings to " + fileName);

            } // while type

        } catch (IOException e) {

            e.printStackTrace();
        }
    }

    // --------------------------------------------------------------------------------
    /**
     * Generate the SQL to update the created and modified date columns in the 
     * stable_id_event tables.
     * Note that this assumes stable_id_event has been updated, so should only
     * be run if this has been done.
     *
     * @param upload Whether or not to actually do the update on the database. 
     *               The SQL is generated regardless.
     */
    public void updateCreatedModified(boolean upload) {

	String[] types = { "gene", "transcript", "translation", "exon" };

	// make everything changed during this mapping session have the same timestamp.
	String ts = new SimpleDateFormat("yyyyMMddHHmmss").format(new Date());

	// update each stable_id table in turn
	for (int i = 0; i < types.length; i++) {

	    String table = types[i] + "_stable_id";

	    String sqlCreated = "UPDATE " + table + " si, stable_id_event sie SET si.created_date='" + ts + "', si.modified_date='" + ts + "' WHERE sie.new_stable_id=si.stable_id AND sie.type='" + types[i] +"' AND sie.old_stable_id IS NULL AND sie.mapping_session_id=" + currentMappingSessionID + ";";
	    
	    String sqlModified = "UPDATE " + table + " si, stable_id_event sie SET si.modified_date='" + ts + "' WHERE sie.new_stable_id=si.stable_id AND sie.type='" + types[i] + "' AND sie.old_stable_id IS NOT NULL AND sie.old_version != sie.new_version AND sie.mapping_session_id=" + currentMappingSessionID + ";";
	    
	    try {
		
		String fileName = rootDir + File.separator + "update_" + table + "_dates.sql";
		OutputStreamWriter writer = new OutputStreamWriter(new FileOutputStream(fileName));

		writer.write(sqlCreated + "\n\n" + sqlModified);
		System.out.println("Wrote timestamp update SQL for " + table + " to " + fileName);

		writer.close();

		// optionally execute
		if (upload) {
		    
		    Statement stmt = conf.getTargetConnection().createStatement();
		    
		    System.out.println("Setting created & modified dates for new " + types[i] + " stable IDs");
		    stmt.executeUpdate(sqlCreated);

		    System.out.println("Setting modified dates for new " + types[i] + " stable IDs");
		    stmt.executeUpdate(sqlModified);

		    stmt.close();

		}

	    } catch (Exception e) {
		e.printStackTrace();
	    }
	    
	}

    }


    // -------------------------------------------------------------------------

} // StableIDMapper

// -------------------------------------------------------------------------
/**
 * Container to represent a mapping - source & target internal IDs, stableID.
 */

class MappingContainer {

    private long sourceInternalID, targetInternalID;

    private String stableID;

    public MappingContainer(long sourceInternalID, long targetInternalID, String stableID) {

        this.sourceInternalID = sourceInternalID;
        this.targetInternalID = targetInternalID;
        this.stableID = stableID;

    }

    public String toString() {

        return sourceInternalID + "\t" + targetInternalID + "\t" + stableID;

    }
}

// -------------------------------------------------------------------------

/**
 * Container to represent a row in the stable_id_event_table.
 */

class StableIDEventContainer {

    private String oldStableID, newStableID;

    private int oldVersion, newVersion;

    private String type;

    private long mappingSessionID;

    /**
     * Create a new StableIDEventContainer. Note either oldStableID or newStableID may be null (but not both at the same time!)
     */
    public StableIDEventContainer(String oldStableID, int oldVersion, String newStableID, int newVersion, String type,
            long mappingSessionID) {

        this.oldStableID = oldStableID;
        this.oldVersion = oldVersion;
        this.newStableID = newStableID;
        this.newVersion = newVersion;
        this.type = type;
        this.mappingSessionID = mappingSessionID;

    }

    /**
     * Copy constructor.
     */
    public StableIDEventContainer(StableIDEventContainer sidec) {

        this.oldStableID = sidec.getOldStableID();
        this.oldVersion = sidec.getOldVersion();
        this.newStableID = sidec.getNewStableID();
        this.newVersion = sidec.getNewVersion();
        this.type = sidec.getType();
        this.mappingSessionID = sidec.getMappingSessionID();

    }

    /**
     * @return Returns the newStableID.
     */
    public String getNewStableID() {

        return newStableID;
    }

    /**
     * @param newStableID The newStableID to set.
     */
    public void setNewStableID(String newStableID) {

        this.newStableID = newStableID;
    }

    /**
     * @return Returns the newVersion.
     */
    public int getNewVersion() {

        return newVersion;
    }

    /**
     * @param newVersion The newVersion to set.
     */
    public void setNewVersion(int newVersion) {

        this.newVersion = newVersion;
    }

    /**
     * @return Returns the oldStableID.
     */
    public String getOldStableID() {

        return oldStableID;
    }

    /**
     * @param oldStableID The oldStableID to set.
     */
    public void setOldStableID(String oldStableID) {

        this.oldStableID = oldStableID;
    }

    /**
     * @return Returns the oldVersion.
     */
    public int getOldVersion() {

        return oldVersion;
    }

    /**
     * @param oldVersion The oldVersion to set.
     */
    public void setOldVersion(int oldVersion) {

        this.oldVersion = oldVersion;
    }

    /**
     * @return Returns the type.
     */
    public String getType() {

        return type;
    }

    /**
     * @param type The type to set.
     */
    public void setType(String type) {

        this.type = type;
    }

    /**
     * @return Returns the currentMappingSessionID.
     */
    public long getMappingSessionID() {

        return mappingSessionID;
    }

    /**
     * @param currentMappingSessionID The currentMappingSessionID to set.
     */
    public void setMappingSessionID(long mappingSessionID) {

        this.mappingSessionID = mappingSessionID;
    }

    public String toString() {

        return type + " " + oldStableID + ":" + oldVersion + " " + newStableID + ":" + newVersion + " " + mappingSessionID;

    }

    public String getKey() {

        return oldStableID + "." + oldVersion + ":" + newStableID + "." + newVersion + ":" + type + ":" + mappingSessionID;

    }
    
}

