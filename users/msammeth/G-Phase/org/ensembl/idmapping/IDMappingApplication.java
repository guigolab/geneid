/*
 * Copyright (C) 2003 EBI, GRL
 * 
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

package org.ensembl.idmapping;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.Exon;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.Translation;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.Driver;
import org.ensembl.util.NamedTimer;
import org.ensembl.util.ProgressPrinter;
import org.ensembl.util.SerialUtil;
import org.ensembl.util.StringUtil;
import org.ensembl.util.Util;

public class IDMappingApplication implements Runnable {

	private Config conf;

	private String rootDir = "";

	private boolean debug = true;

	// put this many deleted genes in the summary email
	private int SAMPLE_SIZE = 20;

	public IDMappingApplication(String configFile) {

		conf = new Config(configFile);

		// TODO - reinstate after debugging
		if (!conf.validateConfig()) {
			System.err.println("Configuration check failed");
			//System.exit(1);
		}

	}

	//---------------------------------------------------------------------

	/**
	 * Creates IDMappingLogs from specified source and target data sets.
	 */

	public static void main(String[] args) {

		if (args.length > 1 || (args.length == 1 && args[0].equals("-h"))) {
			System.out
					.println("Usage: IDMappingApplication {properties file}\n\nIf no properties file is specified, resources/data/idmapping.properties is used.");
			System.exit(1);
		}

		String configFile;

		if (args.length == 1) {
			configFile = args[0];
		} else {
			configFile = ".." + File.separator + "resources" + File.separator
					+ "data" + File.separator + "idmapping.properties";
		}

		// if the config validates do the id mapping
		System.out.println("\n----- Checking configuration -----");

		IDMappingApplication app = new IDMappingApplication(configFile);

		app.run();

	}

	// -------------------------------------------------------------------------
	/**
	 * Do the scoring, and internal and stable ID mapping.
	 */
	public void run() {

		NamedTimer timer = new NamedTimer();
		timer.start("all");

		rootDir = createWorkingDirectory();
		System.out.println("Using working directory " + rootDir);

		timer.start("caching");
		System.out
				.println("\n----- Reading and caching genes, transcripts, translations and exons -----");
		Cache cache = buildCaches();
		timer.stop("caching");
		InternalIDMapper internalIDMapper = new InternalIDMapper(rootDir, cache);

		// ----------------------------
		// SCORE BUILDING

		timer.start("scoring");

		// calculate scores for exons based on overlap and/or sequence
		// similarity
		timer.start("exonScoring");
		System.out.println("\n----- Generating exon scores -----");
		ScoredMappingMatrix exonScoringMatrix = internalIDMapper
				.scoreExons(conf);
		timer.stop("exonScoring");

		
		// calculate scores for transcripts based on exon scores
		timer.start("transcriptScoring");
		System.out.println("\n----- Generating transcript scores -----");
		ScoredMappingMatrix transcriptScoringMatrix = internalIDMapper
				.scoreTranscripts(exonScoringMatrix);
		timer.stop("transcriptScoring");
		// XXX comments
		// do transcript internal ID mapping based on transcript scores
		timer.start("transcriptInternalMapping");
		System.out.println("\n----- Mapping transcript internal IDs -----");
		List transcriptMappings = internalIDMapper
				.doTranscriptMapping(transcriptScoringMatrix);
		timer.stop("transcriptInternalMapping");
		
		// calculate scores for genes based on transcript scores
		timer.start("geneScoring");
		System.out.println("\n----- Generating gene scores -----");
		ScoredMappingMatrix geneScoringMatrix = internalIDMapper.scoreGenes(
				transcriptScoringMatrix, transcriptMappings);
		timer.stop("geneScoring");
		
		timer.stop("scoring");

		// ----------------------------
		// INTERNAL ID MAPPING
		timer.start("internalMapping");

		// do exon internal ID mapping based on exon scores and transcript
		// internal ID mappings
		timer.start("exonInternalMapping");
		System.out.println("\n----- Mapping exon internal IDs -----");
		List exonMappings = internalIDMapper.doExonMapping(exonScoringMatrix,
				transcriptMappings);
		timer.stop("exonInternalMapping");
		
		// do translation internal ID mapping based on transcript internal
		// IDmapping
		timer.start("translationInternalMapping");
		System.out.println("\n----- Mapping translation internal IDs -----");
		List translationMappings = internalIDMapper
				.doTranslationMapping(transcriptMappings);
		timer.stop("translationInternalMapping");

		// do gene internal ID mapping based on gene scores
		timer.start("geneInternalMapping");
		System.out.println("\n----- Mapping gene internal IDs -----");
		List geneMappings = internalIDMapper.doGeneMapping(geneScoringMatrix);
		timer.stop("geneInternalMapping");

		// build lookup tables of generated internal ID mappings
		debug("Caching internal ID mappings");
		cache.cacheMappings(exonMappings, transcriptMappings,
				translationMappings, geneMappings);

		// free up scoring matrices etc
		exonScoringMatrix = null;
		transcriptScoringMatrix = null;
		geneScoringMatrix = null;
		internalIDMapper = null;
		
		timer.stop("internalMapping");

		// ----------------------------
		// STABLE ID MAPPING

		timer.start("stableMapping");

		Driver sourceDriver = conf.getSourceDriver();
		Driver targetDriver = conf.getTargetDriver();

		StableIDMapper stableIDMapper = new StableIDMapper(rootDir, conf);

		timer.start("exonStableMapping");
		System.out.println("\n----- Mapping exon stable IDs -----");
		List newExons = stableIDMapper.mapStableIDs(cache
				.getSourceExonsByInternalID(), cache
				.getTargetExonsByInternalID(), exonMappings);
		timer.stop("exonStableMapping");

		timer.start("transcriptStableMapping");
		System.out.println("\n----- Mapping transcript Stable IDs -----");
		List newTranscripts = stableIDMapper.mapStableIDs(cache
				.getSourceTranscriptsByInternalID(), cache
				.getTargetTranscriptsByInternalID(), transcriptMappings);
		timer.stop("transcriptStableMapping");

		timer.start("translationStableMapping");
		System.out.println("\n----- Mapping translation stable IDs -----");
		List newTranslations = stableIDMapper.mapStableIDs(cache
				.getSourceTranslationsByInternalID(), cache
				.getTargetTranslationsByInternalID(), translationMappings);
		timer.stop("translationStableMapping");

		timer.start("geneStableMapping");
		System.out.println("\n----- Mapping gene stable IDs -----");
		List newGenes = stableIDMapper.mapStableIDs(cache
				.getSourceGenesByInternalID(), cache
				.getTargetGenesByInternalID(), geneMappings);
		timer.stop("geneStableMapping");

		stableIDMapper.dumpDebugMappingsToFile();

		if (doUpload("stableids")) {

			String[] tables = { "exon", "transcript", "translation", "gene" };
			for (int i = 0; i < tables.length; i++) {
				String table = tables[i] + "_stable_id";
				Config.uploadFromFile(
						rootDir + File.separator + table + ".txt", table, conf
								.getTargetConnection(), false);
			}

		}

		timer.stop("stableMapping");

		System.out.println("\n----- Stable ID event generation -----");

		// calculate merged/split stable ID events
		debug("Generating exon-jump stable ID events");
		stableIDMapper.generateMergedSplitEvents(cache, exonMappings,
				transcriptMappings, geneMappings);

		stableIDMapper.dumpLostGeneAndTranscripts(cache);

		// get (and write to file) new mapping session
		stableIDMapper.generateNewMappingSessionID();

		// write stable new ID events to file
		debug("Writing new stable ID events to file");
		stableIDMapper.writeNewStableIDEvents();

		// write existing stable ID events to file
		debug("Writing existing stable ID events to file");
		Util.dumpTableToFile(conf.getSourceConnection(), "stable_id_event",
				rootDir + File.separator + "stable_id_event_existing.txt");

		// upload stable ID events to database if required
		if (doUpload("events")) {

		    debug("Flushing cached exon sequences");
		    cache.flushSequences();

			stableIDMapper.uploadMappingSession();
			stableIDMapper.uploadExistingStableIDEvents();
			stableIDMapper.uploadNewStableIDEvents();
			// propagate stable ID events
			debug("Propagating stable ID events");
			stableIDMapper.propagateStableIDEvents();
			stableIDMapper.writePropagatedStableIDEvents();
			stableIDMapper.uploadPropagatedStableIDEvents();

		} else {

			System.out
					.println("Upload property not set, new mapping session and stable ID events not uploaded/propagated.");

		}
		
		// update created and/or modified date columns in stable_id tables
		stableIDMapper.updateCreatedModified(doUpload("events"));

		// ----------------------------
		// ARCHIVING

		timer.start("archiving");

		// update peptide_archive
		timer.start("peptideArchiving");
		System.out.println("\n----- Updating peptide archive -----");
		Archiver archiver = new Archiver(conf);
		Util.dumpTableToFile(conf.getSourceConnection(), "peptide_archive",
				rootDir + File.separator + "peptide_archive_existing.txt");
		archiver.writeNewPeptideArchiveToFile(cache, stableIDMapper
				.getChanged("translation"), rootDir + File.separator
				+ "peptide_archive_new.txt");
		timer.stop("peptideArchiving");

		// update gene_archive
		timer.start("geneArchiving");
		System.out.println("\n----- Updating gene archive -----");
		Util.dumpTableToFile(conf.getSourceConnection(), "gene_archive",
				rootDir + File.separator + "gene_archive_existing.txt");
		archiver.writeNewGeneArchiveToFile(cache, stableIDMapper
				.getChanged("gene"), stableIDMapper
				.getCurrentMappingSessionID(), rootDir + File.separator
				+ "gene_archive_new.txt");
		timer.stop("geneArchiving");

		if (doUpload("archive")) {

			Config.uploadFromFile(rootDir + File.separator
					+ "peptide_archive_existing.txt", "peptide_archive", conf
					.getTargetConnection(), true);
			Config.uploadFromFile(rootDir + File.separator
					+ "gene_archive_existing.txt", "gene_archive", conf
					.getTargetConnection(), true);
			Config.uploadFromFile(rootDir + File.separator
					+ "peptide_archive_new.txt", "peptide_archive", conf
					.getTargetConnection(), true);
			Config.uploadFromFile(rootDir + File.separator
					+ "gene_archive_new.txt", "gene_archive", conf
					.getTargetConnection(), true);

		}

		timer.stop("archiving");

		// ----------------------------

		timer.stop("all");
		printTimings(timer);

		createSummaryEmail(timer);

		System.out.println("\n----- ID Mapping finished -----");

	}

	// -------------------------------------------------------------------------

	private String createWorkingDirectory() {

		String dirName = "";
		if (System.getProperty("idmapping.base_directory") != null) {
			dirName = System.getProperty("idmapping.base_directory");
		}
		dirName += File.separator
				+ System.getProperty("idmapping.source.database") + "_"
				+ System.getProperty("idmapping.target.database");

		File f = new File(dirName);
		f.mkdir();

		File debug = new File(dirName + File.separator + "debug");
		debug.mkdir();

		return dirName;

	}

	//---------------------------------------------------------------------

	private void printTimings(NamedTimer nt) {

		System.out.println("\nCaching:       " + gf(nt, "caching"));

		System.out.println("\nScoring");
		System.out.println("    Exon:        " + gf(nt, "exonScoring"));
		System.out.println("    Transcript:  " + gf(nt, "transcriptScoring"));
		System.out.println("    Gene:        " + gf(nt, "geneScoring"));
		System.out.println("    Overall:     " + gf(nt, "scoring"));

		System.out.println("\nInternal ID Mapping");
		System.out.println("    Exon:        " + gf(nt, "exonInternalMapping"));
		System.out.println("    Transcript:  "
				+ gf(nt, "transcriptInternalMapping"));
		System.out.println("    Translation: "
				+ gf(nt, "translationInternalMapping"));
		System.out.println("    Gene:        " + gf(nt, "geneInternalMapping"));
		System.out.println("    Overall:     " + gf(nt, "internalMapping"));

		System.out.println("\nStable ID Mapping");
		System.out.println("    Exon:        " + gf(nt, "exonStableMapping"));
		System.out.println("    Transcript:  "
				+ gf(nt, "transcriptStableMapping"));
		System.out.println("    Translation: "
				+ gf(nt, "translationStableMapping"));
		System.out.println("    Gene:        " + gf(nt, "geneStableMapping"));
		System.out.println("    Overall:     " + gf(nt, "stableMapping"));

		System.out.println("\nArchiving:");
		System.out.println("    Peptide:      " + gf(nt, "peptideArchiving"));
		System.out.println("    Gene:         " + gf(nt, "geneArchiving"));
		System.out.println("    Overall:      " + gf(nt, "archiving"));

		System.out.println("\nTotal duration: " + gf(nt, "all"));

	}

	private String gf(NamedTimer nt, String s) {

		return nt.format(nt.getDuration(s));

	}

	// -------------------------------------------------------------------------

	private boolean doUpload(String prop) {

		String fullProp = System.getProperty("idmapping.upload." + prop);
		return (fullProp != null && fullProp.equals("yes"));

	}

	// -------------------------------------------------------------------------

	private void createSummaryEmail(NamedTimer nt) {

		StringBuffer message = new StringBuffer();

		message.append("ID mapping complete for "
				+ System.getProperty("idmapping.source.database") + " -> "
				+ System.getProperty("idmapping.target.database") + "\n");

		message.append("\nResults:\n\n");
		String[] types = { "exon", "transcript", "translation", "gene" };
		for (int i = 0; i < types.length; i++) {
			message.append(StringUtil.readTextFile(rootDir + File.separator
					+ types[i] + "_mapping_statistics.txt"));
			message.append("\n");
		}

		String[] uploads = { "stableids", "events", "archive" };
		String[] uploadDescriptions = { "Stable IDs",
				"Stable ID events and mapping session",
				"Gene and peptide archiving information" };

		for (int j = 0; j < uploads.length; j++) {

			String s = doUpload(uploads[j]) ? "" : "not ";
			message.append(uploadDescriptions[j] + " were " + s
					+ "uploaded to "
					+ System.getProperty("idmapping.target.database") + "\n");

		}

		message.append("\n");

		message.append("A sample of the first " + SAMPLE_SIZE
				+ " deleted known genes is at the end of this email.\n");

		String[] gt = { "genes", "transcripts" };
		for (int k = 0; k < gt.length; k++) {

			message.append("A full list of " + gt[k]
					+ " which were deleted is in " + rootDir + File.separator
					+ gt[k] + "_lost_deleted.txt" + "\n");
			message.append("A full list of " + gt[k]
					+ " which were lost due to merging, and the " + gt[k]
					+ " into which they were merged is in " + rootDir
					+ File.separator + gt[k] + "_lost_merged.txt" + "\n");

		}

		message.append("\nTotal run time: " + nt.format(nt.getDuration("all")));

		message.append("\n----------------------------------------------------------------------\n");
		
		String[] deletedGenes = StringUtil.readTextFile(
				rootDir + File.separator + "genes_lost_deleted.txt").split(
				"\\n");

		if (deletedGenes.length > 1) {

		    message.append("First " + SAMPLE_SIZE
				   + " known genes which were deleted:\n\n");

		    int knownCount = 0;
		    for (int i = 0; i < deletedGenes.length; i++) {
			String[] bits = deletedGenes[i].split("\\t");
			if (bits[1].equals("KNOWN")) {
			    message.append(bits[0]);
			    knownCount++;
			    String link = Util.makeEnsemblGeneLink(bits[0]);
			    if (link != null) {
				message.append("\t" + link);
			    }
			    message.append("\n");
			}
			if (knownCount > SAMPLE_SIZE) {
			    break;
			}
		    }
		    
		}

		String fileName = rootDir + File.separator + "summary_email.txt";

		try {

			OutputStreamWriter writer = new OutputStreamWriter(
					new FileOutputStream(fileName));
			writer.write(message.toString());
			writer.close();

		} catch (IOException e) {

			e.printStackTrace();
		}

		System.out
				.println("\nSummary information suitable for emailing written to "
						+ fileName);

	}

	// -------------------------------------------------------------------------

	private Cache buildCaches() {

		Cache cache = new Cache();

		// read from file if already done
		String fileName = rootDir + File.separator + "cache.ser";
		File f = new File(fileName);
		if (!f.exists()) {

			try {

				List allSourceGenes = conf.getSourceDriver().getGeneAdaptor()
						.fetchAll();
				System.out.println("Total of " + allSourceGenes.size()
						+ " source genes");
				ProgressPrinter spp = new ProgressPrinter(0, allSourceGenes
						.size(), "% of source genes read");
				int i = 0;

				Iterator sgit = allSourceGenes.iterator();
				while (sgit.hasNext()) {

					Gene gene = (Gene) sgit.next();
					gene.isKnown(); // force lazy-load of isKnown
					cache.getSourceGenesByInternalID().put(
							new Long(gene.getInternalID()), gene);
					cache.getSourceGenesByStableID().put(gene.getAccessionID(),
							gene);

					List transcripts = gene.getTranscripts();
					Iterator stit = transcripts.iterator();
					while (stit.hasNext()) {

						Transcript transcript = (Transcript) stit.next();
						transcript.isKnown(); // force lazy-load of isKnown
						Long transcriptID = new Long(transcript.getInternalID());
						cache.getSourceTranscriptsByInternalID().put(
								transcriptID, transcript);
						cache.getSourceGeneByTranscriptInternalID().put(
								transcriptID, gene);
						cache.getSourceTranscriptsByStableID().put(
								transcript.getAccessionID(), transcript);

						Translation translation = transcript.getTranslation();
						if (translation != null) { // ignore pseudogenes etc
							translation.isKnown(); // force lazy-load of isKnown
							cache.getSourceTranslationsByInternalID().put(
									new Long(translation.getInternalID()),
									translation);
							cache.getSourceTranslationsByTranscriptInternalID()
									.put(transcriptID, translation);
							cache.getSourceTranslationsByStableID().put(
									translation.getAccessionID(), translation);
						}

						List exons = transcript.getExons();
						Iterator seit = exons.iterator();
						while (seit.hasNext()) {

							Exon exon = (Exon) seit.next();
							exon.getSequence(); // forces sequence to be loaded
							Long exonID = new Long(exon.getInternalID());
							//if (exon.getInternalID() ==192948) {
							//System.out.println("In read: " +
							// exon.getLocation().toString());
							//}
							cache.getSourceExonsByInternalID()
									.put(exonID, exon);
							Util.addToMapList(cache
									.getSourceTranscriptsByExonInternalID(),
									exonID, transcript);
							cache.getSourceGeneByExonInternalID().put(exonID,
									gene);

						}
					}

					spp.printUpdate(i++);

				}
				spp.printUpdate(allSourceGenes.size());
				System.out.println("");

				// ------------------------
				List allTargetGenes = conf.getTargetDriver().getGeneAdaptor()
						.fetchAll();
				System.out.println("\nTotal of " + allTargetGenes.size()
						+ " target genes");
				ProgressPrinter tpp = new ProgressPrinter(0, allTargetGenes
						.size(), "% of target genes read");
				int j = 0;
				Iterator tgit = allTargetGenes.iterator();
				while (tgit.hasNext()) {

					Gene gene = (Gene) tgit.next();
					gene.isKnown();
					cache.getTargetGenesByInternalID().put(
							new Long(gene.getInternalID()), gene);

					List transcripts = gene.getTranscripts();
					Iterator ttit = transcripts.iterator();
					while (ttit.hasNext()) {

						Transcript transcript = (Transcript) ttit.next();
						transcript.isKnown();
						Long transcriptID = new Long(transcript.getInternalID());
						cache.getTargetTranscriptsByInternalID().put(
								transcriptID, transcript);
						cache.getTargetGeneByTranscriptInternalID().put(
								transcriptID, gene);

						Translation translation = transcript.getTranslation();
						if (translation != null) { // ignore pseudogenes etc
							translation.isKnown();
							cache.getTargetTranslationsByInternalID().put(
									new Long(translation.getInternalID()),
									translation);
							cache.getTargetTranslationsByTranscriptInternalID()
									.put(transcriptID, translation);
						}

						List exons = transcript.getExons();
						Iterator teit = exons.iterator();
						while (teit.hasNext()) {

							Exon exon = (Exon) teit.next();
							exon.getSequence(); // forces sequence to be loaded
							Long exonID = new Long(exon.getInternalID());
							cache.getTargetExonsByInternalID()
									.put(exonID, exon);
							Util.addToMapList(cache
									.getTargetTranscriptsByExonInternalID(),
									exonID, transcript);
							cache.getTargetGeneByExonInternalID().put(exonID,
									gene);

						}
					}

					tpp.printUpdate(j++);

				}
				tpp.printUpdate(allTargetGenes.size());
				System.out.println("");

			} catch (AdaptorException ae) {

				ae.printStackTrace();

			}

			debug("Writing cache to " + rootDir + File.separator + "cache.ser");
			serialiseCache(cache, rootDir + File.separator + "cache.ser");
			debug("Finished writing cache");

		} else {

			System.out
					.println("Using cached gene, transcript, translation and exon information in "
							+ fileName);
			cache = (Cache) SerialUtil.readObject(fileName);

		}

		//Exon exon = (Exon) cache.getSourceExonsByInternalID().get(new
		// Long(192948));
		//System.out.println("From cache: " + exon.getLocation().toString());

		return cache;

	}

	// -------------------------------------------------------------------------

	public void serialiseCache(Cache cache, String fileName) {

		SerialUtil.writeObject(cache, fileName);

	}

	// -------------------------------------------------------------------------

	private void debug(String s) {

		if (debug) {
			System.out.println(s);
		}

	}

	
	// -------------------------------------------------------------------------

} // IDMappingApplication

//---------------------------------------------------------------------

class Cache implements Serializable {
	
	private static final long serialVersionUID = 1L;

	private Map sourceGenesByInternalID = new HashMap();

	private Map targetGenesByInternalID = new HashMap();

	private Map sourceGenesByStableID = new HashMap();

	private Map sourceTranscriptsByStableID = new HashMap();

	private Map sourceExonsByInternalID = new HashMap();

	private Map targetExonsByInternalID = new HashMap();

	private Map sourceTranscriptsByInternalID = new HashMap();

	private Map targetTranscriptsByInternalID = new HashMap();

	private Map sourceTranslationsByInternalID = new HashMap();

	private Map targetTranslationsByInternalID = new HashMap();

	private Map sourceTranscriptsByExonInternalID = new HashMap(); // Map of
																   // lists

	private Map targetTranscriptsByExonInternalID = new HashMap(); // Map of
																   // lists

	private Map sourceTranslationsByTranscriptInternalID = new HashMap();

	private Map targetTranslationsByTranscriptInternalID = new HashMap();

	private Map sourceTranslationsByStableID = new HashMap();

	private Map sourceGeneByTranscriptInternalID = new HashMap();

	private Map targetGeneByTranscriptInternalID = new HashMap();

	private Map sourceGeneByExonInternalID = new HashMap();

	private Map targetGeneByExonInternalID = new HashMap();

	private Map exonMappingsMap = new HashMap();

	private Map transcriptMappingsMap = new HashMap();

	private Map translationMappingsMap = new HashMap();

	private Map geneMappingsMap = new HashMap();

	/**
	 * @return Returns the sourceExonsByInternalID.
	 */
	public Map getSourceExonsByInternalID() {

		return sourceExonsByInternalID;
	}

	/**
	 * @return Returns the sourceGenesByInternalID.
	 */
	public Map getSourceGenesByInternalID() {

		return sourceGenesByInternalID;
	}

	/**
	 * @return Returns the sourceGenesByTranscriptInternalID.
	 */
	public Map getSourceGeneByTranscriptInternalID() {

		return sourceGeneByTranscriptInternalID;
	}

	/**
	 * @return Returns the sourceTranscriptsByExonInternalID.
	 */
	public Map getSourceTranscriptsByExonInternalID() {

		return sourceTranscriptsByExonInternalID;
	}

	/**
	 * @return Returns the sourceTranscriptsByInternalID.
	 */
	public Map getSourceTranscriptsByInternalID() {

		return sourceTranscriptsByInternalID;
	}

	/**
	 * @return Returns the sourceTranslationsByInternalID.
	 */
	public Map getSourceTranslationsByInternalID() {

		return sourceTranslationsByInternalID;
	}

	/**
	 * @return Returns the targetTranslationsByInternalID.
	 */
	public Map getTargetTranslationsByInternalID() {

		return targetTranslationsByInternalID;
	}

	/**
	 * @return Returns the sourceTranslationsByTranscriptInternalID.
	 */
	public Map getSourceTranslationsByTranscriptInternalID() {

		return sourceTranslationsByTranscriptInternalID;
	}

	/**
	 * @return Returns the targetExonsByInternalID.
	 */
	public Map getTargetExonsByInternalID() {

		return targetExonsByInternalID;
	}

	/**
	 * @return Returns the targetGenesByInternalID.
	 */
	public Map getTargetGenesByInternalID() {

		return targetGenesByInternalID;
	}

	/**
	 * @return Returns the targetGenesByTranscriptInternalID.
	 */
	public Map getTargetGeneByTranscriptInternalID() {

		return targetGeneByTranscriptInternalID;
	}

	/**
	 * @return Returns the targetTranscriptsByExonInternalID.
	 */
	public Map getTargetTranscriptsByExonInternalID() {

		return targetTranscriptsByExonInternalID;
	}

	/**
	 * @return Returns the targetTranscriptsByInternalID.
	 */
	public Map getTargetTranscriptsByInternalID() {

		return targetTranscriptsByInternalID;
	}

	/**
	 * @return Returns the targetTranslationsByTranscriptInternalID.
	 */
	public Map getTargetTranslationsByTranscriptInternalID() {

		return targetTranslationsByTranscriptInternalID;
	}

	// -------------------------------------------------------------------------

	/**
	 * Build and store hashtables of source-target internal ID mappings from
	 * lists of Entry objects.
	 */
	public void cacheMappings(List exonMappings, List transcriptMappings,
			List translationMappings, List geneMappings) {

		exonMappingsMap = buildMapping(exonMappings);
		transcriptMappingsMap = buildMapping(transcriptMappings);
		translationMappingsMap = buildMapping(translationMappings);
		geneMappingsMap = buildMapping(geneMappings);

	}

	private Map buildMapping(List mappings) {

		Map result = new HashMap();
		Iterator it = mappings.iterator();
		while (it.hasNext()) {
			Entry entry = (Entry) it.next();
			result
					.put(new Long(entry.getSource()), new Long(entry
							.getTarget()));
		}

		return result;

	}

	// -------------------------------------------------------------------------

	/**
	 * @return Returns the exonMappingsMap.
	 */
	public Map getExonMappingsMap() {

		return exonMappingsMap;
	}

	/**
	 * @return Returns the geneMappingsMap.
	 */
	public Map getGeneMappingsMap() {

		return geneMappingsMap;
	}

	/**
	 * @return Returns the transcriptMappingsMap.
	 */
	public Map getTranscriptMappingsMap() {

		return transcriptMappingsMap;
	}

	/**
	 * @return Returns the translationMappingsMap.
	 */
	public Map getTranslationMappingsMap() {

		return translationMappingsMap;
	}

	/**
	 * @return Returns the sourceGenesByExonInternalID.
	 */
	public Map getSourceGeneByExonInternalID() {

		return sourceGeneByExonInternalID;
	}

	/**
	 * @return Returns the targetGenesByExonInternalID.
	 */
	public Map getTargetGeneByExonInternalID() {

		return targetGeneByExonInternalID;
	}

	/**
	 * @param sourceGenesByExonInternalID
	 *            The sourceGenesByExonInternalID to set.
	 */
	public void setSourceGeneByExonInternalID(Map sourceGenesByExonInternalID) {

		this.sourceGeneByExonInternalID = sourceGenesByExonInternalID;
	}

	/**
	 * @param targetGenesByExonInternalID
	 *            The targetGenesByExonInternalID to set.
	 */
	public void setTargetGeneByExonInternalID(Map targetGenesByExonInternalID) {

		this.targetGeneByExonInternalID = targetGenesByExonInternalID;
	}

	/**
	 * @return Returns the sourceTranslationsByStableID.
	 */
	public Map getSourceTranslationsByStableID() {

		return sourceTranslationsByStableID;
	}

	/**
	 * @return Returns the sourceGenesByStableID.
	 */
	public Map getSourceGenesByStableID() {

		return sourceGenesByStableID;
	}

	/**
	 * @return Returns the sourceTranscriptsByStableID.
	 */
	public Map getSourceTranscriptsByStableID() {

		return sourceTranscriptsByStableID;
	}

	/**
	 * Hopefully removes all sequences from cached Exons as we need to
	 * be careful with the space ..
	 *
	 */
	public void flushSequences() {
		Iterator i = sourceExonsByInternalID.values().iterator();
		while (i.hasNext()) {
			Exon e = (Exon) i.next();
			e.setSequence(null);
		}
		
		i = targetExonsByInternalID.values().iterator();
		while (i.hasNext()) {
			Exon e = (Exon) i.next();
			e.setSequence(null);
		}
	}
}

// -------------------------------------------------------------------------

