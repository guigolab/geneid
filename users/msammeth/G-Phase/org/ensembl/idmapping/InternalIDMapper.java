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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.Translation;
import org.ensembl.util.SerialUtil;

/**
 * Perform mapping of internal IDs based on scores and other criteria.
 */
public class InternalIDMapper {

	private String rootDir, debugDir;

	private boolean debug = true;

	private OutputStreamWriter ambiGeneOut, ambiTranscriptOut;

	// scores are considered the same if (2.0 * (s1-s2))/(s1 + s2) < this
	private static final float SIMILAR_SCORE_RATIO = 0.01f;

	/**
	 * The cache provides access to gene, transcripts and exons of hte source
	 * and target database
	 */
	private Cache cache;

	public InternalIDMapper(String rootDir, Cache cache) {

		this.rootDir = rootDir;
		this.debugDir = rootDir + File.separator + "debug";
		this.cache = cache;

	}

	// -------------------------------------------------------------------------
	/**
	 * Build exon scoring matrix based on overlap and/or sequence matching using
	 * exonerate.
	 */
	public ScoredMappingMatrix scoreExons(Config conf) {

		ScoredMappingMatrix exonScoringMatrix = new ScoredMappingMatrix();

		// read from file if already done
		String exonScoreFileName = rootDir + File.separator + "exon_scores.ser";
		File f = new File(exonScoreFileName);
		if (!f.exists()) {

			System.out
					.println("Did not find existing exon score matrix, will build a new one");
			ExonDirectMapper exonDirectMapper = new ExonDirectMapper(conf,
					cache);

			try {

				if (exonDirectMapper.mappedComparable()) {

					debug("Building exon overlap scores");
					exonDirectMapper.buildOverlapScoring();
					// Retrieve source exons, store projections to common
					// coordinate system. Retrieve target exons. Store
					// projection. Find overlaps.
					exonScoringMatrix = exonDirectMapper.getScoringMatrix();
					SerialUtil
							.writeObject(exonScoringMatrix, exonScoreFileName);
					debug("Wrote exon scoring matrix to " + exonScoreFileName);

				} else {

					// no direct exon mapping applicable
					System.out
							.println("No direct exon mapping - will only use exonerate (if enabled)");

				}

			} catch (Exception e) {

				System.err.println("Exception doing exon overlap scoring: ");
				e.printStackTrace();

			}

			// dump low-scoring exons to a FASTA file and map using exonerate
			if (useExonerate()) {

				System.out
						.println("\n----- Generating exon scores using exonerate -----");
				debug("Dumping exons to FASTA files for exonerate");
				new ExonDumper().dumpFilteredExons(cache, exonScoringMatrix,
						rootDir);
				debug("Running exonerate");
				ScoredMappingMatrix exonerateScoringMatrix = new ExonerateRunner()
						.run(rootDir);
				debug("Exonerate scoring matrix size: "
						+ exonerateScoringMatrix.getEntryCount());
				System.out
						.print("Combining direct and exonerate scores; size & average before "
								+ exonScoringMatrix.getEntryCount()
								+ " "
								+ exonScoringMatrix.getAverageScore());
				exonScoringMatrix.combineWith(exonerateScoringMatrix);
				exonerateScoringMatrix = null;
				System.out.println(" after "
						+ exonScoringMatrix.getEntryCount() + " "
						+ exonScoringMatrix.getAverageScore());
			}

		} else {

			System.out.println("Using existing exon score matrix in "
					+ exonScoreFileName);
			exonScoringMatrix = (ScoredMappingMatrix) SerialUtil
					.readObject(exonScoreFileName);

		}

		if (debug) {
			exonScoringMatrix.dumpToFile(debugDir, "exon_scores.txt");
		}

		debug(exonScoringMatrix.toString());

		return exonScoringMatrix;

	}

	// -------------------------------------------------------------------------

	public ScoredMappingMatrix scoreTranscripts(
			ScoredMappingMatrix exonScoringMatrix) {

		ScoredMappingMatrix transcriptScoringMatrix = null;

		// read from file if already done
		String transcriptScoreFileName = rootDir + File.separator
				+ "transcript_scores.ser";
		File f = new File(transcriptScoreFileName);
		if (!f.exists()) {

			System.out
					.println("Did not find existing transcript score matrix, will build a new one");
			// build a transcript score matrix from exon score matrix
			transcriptScoringMatrix = new TranscriptScoreBuilder(cache)
					.buildScoringMatrix(exonScoringMatrix);
			debug("Transcript scoring matrix has = "
					+ transcriptScoringMatrix.getEntryCount() + " entries");
			SerialUtil.writeObject(transcriptScoringMatrix,
					transcriptScoreFileName);
			System.out.println("Wrote transcript scoring matrix to "
					+ transcriptScoreFileName);

		} else {

			System.out.println("Using existing transcript score matrix in "
					+ transcriptScoreFileName);
			transcriptScoringMatrix = (ScoredMappingMatrix) SerialUtil
					.readObject(transcriptScoreFileName);

		}

		debug(transcriptScoringMatrix.toString());

		if (debug) {
			transcriptScoringMatrix.dumpToFile(debugDir,
					"transcript_scores.txt");
		}

		return transcriptScoringMatrix;

	}

	// -------------------------------------------------------------------------
	/**
	 * Score genes based on transcript scores.
	 */
	public ScoredMappingMatrix scoreGenes(ScoredMappingMatrix transcriptScores,
			List transcriptMappings) {

		ScoredMappingMatrix geneScoringMatrix = null;

		// read from file if already done
		String geneScoreFileName = rootDir + File.separator + "gene_scores.ser";
		File f = new File(geneScoreFileName);
		if (!f.exists()) {

			System.out
					.println("Did not find existing gene score matrix, will build a new one");

			// create a new scored mapping matrix based on transcript mappings
			// point transcriptScores at the new matrix
			transcriptScores = new ScoredMappingMatrix(transcriptMappings);

			// build a gene score matrix from transcript score matrix
			geneScoringMatrix = new GeneScoreBuilder(cache)
					.buildScoringMatrix(transcriptScores);
			debug("Gene scoring matrix has = "
					+ geneScoringMatrix.getEntryCount() + " entries");
			SerialUtil.writeObject(geneScoringMatrix, geneScoreFileName);
			System.out.println("Wrote gene scoring matrix to "
					+ geneScoreFileName);

		} else {

			System.out.println("Using existing gene score matrix in "
					+ geneScoreFileName);
			geneScoringMatrix = (ScoredMappingMatrix) SerialUtil
					.readObject(geneScoreFileName);

		}

		debug(geneScoringMatrix.toString());
		if (debug) {
			geneScoringMatrix.dumpToFile(debugDir, "gene_scores.txt");
		}

		return geneScoringMatrix;

	}

	// -------------------------------------------------------------------------
	/**
	 * Map transcript internal IDs based on transcript scores, rejecting
	 * mappings that are ambiguous.
	 * 
	 * @return A list of Entry objects containing the mappings.
	 */
	public List doTranscriptMapping(ScoredMappingMatrix transcriptScores) {

		List mappings = new ArrayList();
		;

		// read from file if already done
		String transcriptMappingFileName = rootDir + File.separator
				+ "transcript_mappings.ser";
		File f = new File(transcriptMappingFileName);
		if (!f.exists()) {

			int ambiguousMappings = 0;

			Map sourcesDone = new HashMap();
			Map targetsDone = new HashMap();

			List envelopes = transcriptScores.buildEnvelopes();

			Iterator envIt = envelopes.iterator();
			while (envIt.hasNext()) {

				Envelope env = (Envelope) envIt.next();
				List entries = env.getEntries();
				Iterator entryIt = entries.iterator();
				while (entryIt.hasNext()) {

					Entry e = (Entry) entryIt.next();
					Long sourceID = new Long(e.getSource());
					Long targetID = new Long(e.getTarget());
					if (!sourcesDone.containsKey(sourceID)
							&& !targetsDone.containsKey(targetID)) {

						sourcesDone.put(sourceID, sourceID);
						targetsDone.put(targetID, targetID);

						// only map if there are no other ambiguous mappings for
						// this source or target
						Entry ambiguousEntry;
						if ((ambiguousEntry = ambiguousMappingsExist(e,
								transcriptScores)) == null) {

							mappings.add(e);

						} else {
							ambiguousTranscript(e, ambiguousEntry);
							ambiguousMappings++;
						}

					}

				}

			}

			System.out.println("Got " + mappings.size()
					+ " transcript mappings, rejected " + ambiguousMappings
					+ " ambiguous mappings");

			SerialUtil.writeObject(mappings, transcriptMappingFileName);
			System.out.println("Wrote transcript internal ID mappings to "
					+ transcriptMappingFileName);

		} else {

			System.out
					.println("Using existing transcript internal ID mappings in "
							+ transcriptMappingFileName);
			mappings = (ArrayList) SerialUtil
					.readObject(transcriptMappingFileName);

		}
		if (ambiTranscriptOut != null) {
			try {
				ambiTranscriptOut.close();
			} catch (Exception x) {
			}
		}
		//if (debug) {
		//    dumpMappingsToFile(mappings, "transcript_mappings.txt");
		//}

		return mappings;

	}

	// -------------------------------------------------------------------------

	private Entry ambiguousMappingsExist(Entry e, ScoredMappingMatrix scores) {

		// iterate over related sources and targets
		// if any has a score similar (within MAPPING_THRESHOLD of e's score)
		// it is ambiguous and not considered.
		long[] targets = scores.getTargetsForSource(e.getSource());
		Entry otherEntry = null;
		for (int i = 0; i < targets.length; i++) {
			otherEntry = scores.getEntry(e.getSource(), targets[i]);
			if (targets[i] != e.getTarget()
					&& scoresSimilar(e.getScore(), otherEntry.getScore())) {
				return otherEntry;
			}
		}

		long[] sources = scores.getSourcesForTarget(e.getTarget());
		for (int i = 0; i < sources.length; i++) {
			otherEntry = scores.getEntry(sources[i], e.getTarget());
			if (sources[i] != e.getSource()
					&& scoresSimilar(e.getScore(), otherEntry.getScore())) {
				return otherEntry;
			}
		}

		return null;

	}

	// -------------------------------------------------------------------------

	private boolean scoresSimilar(float s1, float s2) {

		// always allow for exact match to be mapped in favor of very similar
		// match
		if (s1 == 1.0f && s2 < 1.0f) {
			return false;
		}

		float diff = Math.abs(s1 - s2);

		float pc = (2.0f * diff) / (s1 + s2);

		return pc < SIMILAR_SCORE_RATIO;

	}

	// -------------------------------------------------------------------------
	/**
	 * Map exon internal IDs based on exon scores and transcript mappings.
	 * 
	 * @return A list of Entry objects containing the mappings.
	 */
	public List doExonMapping(ScoredMappingMatrix exonScores,
			List transcriptMappings) {

		List mappings = new ArrayList();

		// read from file if already done
		String exonMappingFileName = rootDir + File.separator
				+ "exon_mappings.ser";
		File f = new File(exonMappingFileName);
		if (!f.exists()) {

			int ambiguousMappings = 0;

			// make a temporary map of transcript mappings.
			Map tmpTranscriptMappings = new HashMap();

			Iterator it = transcriptMappings.iterator();
			while (it.hasNext()) {
				Entry e = (Entry) it.next();
				tmpTranscriptMappings.put(new Long(e.getSource()), new Long(e
						.getTarget()));
			}

			Map sourcesDone = new HashMap();
			Map targetsDone = new HashMap();
			List envelopes = exonScores.buildEnvelopes();
			Iterator envIt = envelopes.iterator();
			while (envIt.hasNext()) {

				Envelope env = (Envelope) envIt.next();
				List entries = env.getEntries();
				Iterator entryIt = entries.iterator();
				while (entryIt.hasNext()) {

					Entry e = (Entry) entryIt.next();
					Long sourceID = new Long(e.getSource());
					Long targetID = new Long(e.getTarget());
					if (!sourcesDone.containsKey(sourceID)
							&& !targetsDone.containsKey(targetID)) {

						// increase scores of exons that are part of transcripts
						// that are already mapped
						if (exonsTranscriptsMap(e, tmpTranscriptMappings, cache
								.getSourceTranscriptsByExonInternalID(), cache
								.getTargetTranscriptsByExonInternalID())) {
							e.setScore(e.getScore()
									* (1.0f + (2.0f * SIMILAR_SCORE_RATIO)));
						}

						sourcesDone.put(sourceID, sourceID);
						targetsDone.put(targetID, targetID);

						// only map if there are no other ambiguous mappings for
						// this source or target
						if (ambiguousMappingsExist(e, exonScores) == null) {

							mappings.add(e);

						} else {

							ambiguousMappings++;

						}

					}

				}

			}

			System.out.println("Got " + mappings.size()
					+ " exon mappings, rejected " + ambiguousMappings
					+ " ambiguous mappings");

			SerialUtil.writeObject(mappings, exonMappingFileName);
			System.out.println("Wrote exon internal ID mappings to "
					+ exonMappingFileName);

		} else {

			System.out.println("Using existing exon internal ID mappings in "
					+ exonMappingFileName);
			mappings = (ArrayList) SerialUtil.readObject(exonMappingFileName);

		}

		//if (debug) {
		//    dumpMappingsToFile(mappings, "exon_mappings.txt");
		//}

		return mappings;

	}

	// -------------------------------------------------------------------------
	/**
	 * @return true if the source and target exons belong to transcripts that
	 *         have mappings.
	 */
	private boolean exonsTranscriptsMap(Entry e, Map transcriptMappings,
			Map sourceTranscriptsByExonID, Map targetTranscriptsByExonID) {

		List sourceTranscripts = (List) sourceTranscriptsByExonID.get(new Long(
				e.getSource()));
		List targetTranscripts = (List) targetTranscriptsByExonID.get(new Long(
				e.getTarget()));

		Iterator stit = sourceTranscripts.iterator();
		while (stit.hasNext()) {

			Transcript sourceTranscript = (Transcript) stit.next();
			Iterator ttit = targetTranscripts.iterator();
			while (ttit.hasNext()) {

				Transcript targetTranscript = (Transcript) ttit.next();
				Long target = (Long) transcriptMappings.get(new Long(
						sourceTranscript.getInternalID()));
				if (target != null) {
					if (target.longValue() == targetTranscript.getInternalID()) {
						return true;
					}
				}
			}
		}

		return false;
	}

	// -------------------------------------------------------------------------
	/**
	 * Map translation internal IDs based on transcript mappings.
	 * 
	 * @return A list of Entry objects containing the mappings.
	 */
	public List doTranslationMapping(List transcriptMappings) {

		List translationMappings = new ArrayList();

		// read from file if already done
		String translationMappingFileName = rootDir + File.separator
				+ "translation_mappings.ser";
		File f = new File(translationMappingFileName);
		if (!f.exists()) {

			int transcriptsWithoutTranslations = 0;

			Iterator it = transcriptMappings.iterator();
			while (it.hasNext()) {

				Entry transcriptEntry = (Entry) it.next();
				Long sourceTranscriptID = new Long(transcriptEntry.getSource());
				Translation sourceTranslation = (Translation) cache
						.getSourceTranslationsByTranscriptInternalID().get(
								sourceTranscriptID);
				Long targetTranscriptID = new Long(transcriptEntry.getTarget());
				Translation targetTranslation = (Translation) cache
						.getTargetTranslationsByTranscriptInternalID().get(
								targetTranscriptID);

				// avoid storing translation mappings for transcripts that have
				// no translation (e.g. pseudogenes)
				if (sourceTranslation != null && targetTranslation != null) {

					// note the score in the translation mapping Entry is the
					// corresponding TRANSCRIPT score
					Entry translationEntry = new Entry(sourceTranslation
							.getInternalID(),
							targetTranslation.getInternalID(), transcriptEntry
									.getScore());
					translationMappings.add(translationEntry);

				} else {

					transcriptsWithoutTranslations++;
				}

			}

			System.out.println("Skipped " + transcriptsWithoutTranslations
					+ " transcripts without translations");

			SerialUtil.writeObject(translationMappings,
					translationMappingFileName);
			System.out.println("Wrote translation internal ID mappings to "
					+ translationMappingFileName);

		} else {

			System.out
					.println("Using existing translation internal ID mappings in "
							+ translationMappingFileName);
			translationMappings = (ArrayList) SerialUtil
					.readObject(translationMappingFileName);

		}

		return translationMappings;

	}

	// -------------------------------------------------------------------------
	/**
	 * Map gene internal IDs based on gene scores.
	 * 
	 * @return A list of Entry objects containing the mappings.
	 */
	public List doGeneMapping(ScoredMappingMatrix geneScores) {

		List mappings = new ArrayList();
		;

		// read from file if already done
		String geneMappingFileName = rootDir + File.separator
				+ "gene_mappings.ser";
		File f = new File(geneMappingFileName);
		if (!f.exists()) {

			int ambiguousMappings = 0;

			Map sourcesDone = new HashMap();
			Map targetsDone = new HashMap();

			List envelopes = geneScores.buildEnvelopes();

			Iterator envIt = envelopes.iterator();
			while (envIt.hasNext()) {

				Envelope env = (Envelope) envIt.next();
				List entries = env.getEntries();
				Iterator entryIt = entries.iterator();
				while (entryIt.hasNext()) {

					Entry e = (Entry) entryIt.next();
					Long sourceID = new Long(e.getSource());
					Long targetID = new Long(e.getTarget());
					if (!sourcesDone.containsKey(sourceID)
							&& !targetsDone.containsKey(targetID)) {

						sourcesDone.put(sourceID, sourceID);
						targetsDone.put(targetID, targetID);

						// only map if there are no other ambiguous mappings for
						// this source or target
						Entry ambiguousEntry;
						if ((ambiguousEntry = ambiguousMappingsExist(e,
								geneScores)) == null) {

							mappings.add(e);

						} else {
							ambiguousGene(e, ambiguousEntry);
							ambiguousMappings++;

						}

					}

				}

			}

			System.out.println("Got " + mappings.size()
					+ " gene mappings, rejected " + ambiguousMappings
					+ " ambiguous mappings");

			SerialUtil.writeObject(mappings, geneMappingFileName);
			System.out.println("Wrote gene internal ID mappings to "
					+ geneMappingFileName);

		} else {

			System.out.println("Using existing gene internal ID mappings in "
					+ geneMappingFileName);
			mappings = (ArrayList) SerialUtil.readObject(geneMappingFileName);
		}

		if (ambiGeneOut != null) {
			try {
				ambiGeneOut.close();
			} catch (Exception x) {
			}
		}
		//if (debug) {
		//    dumpMappingsToFile(mappings, "gene_mappings.txt");
		//}

		return mappings;

	}

	// -------------------------------------------------------------------------

	public void ambiguousGene(Entry e, Entry other) {
		try {
			if (ambiGeneOut == null) {
				ambiGeneOut = new OutputStreamWriter(new FileOutputStream(
						debugDir + File.separator + "ambiguous_Genes.txt"));
			}
			ambiGeneOut.write("Gene ambiguous, lost" + "\n");
			ambiGeneOut.write("  " + e.toString() + "\n");
			ambiGeneOut.write("  " + other.toString() + "\n");
		} catch (IOException exc) {
			exc.printStackTrace();
		}
	}

	public void ambiguousTranscript(Entry e, Entry other) {
		try {
			if (ambiTranscriptOut == null) {
				ambiTranscriptOut = new OutputStreamWriter(
						new FileOutputStream(debugDir + File.separator
								+ "ambiguous_Transcripts.txt"));
			}
			ambiTranscriptOut.write("Transcript ambiguous, lost" + "\n");
			ambiTranscriptOut.write("  " + e.toString() + "\n");
			ambiTranscriptOut.write("  " + other.toString() + "\n");
		} catch (IOException exc) {
			exc.printStackTrace();
		}
	}

	public void dumpMappingsToFile(List mappings, String outputFileName) {

		try {

			Iterator it = mappings.iterator();

			OutputStreamWriter writer = new OutputStreamWriter(
					new FileOutputStream(rootDir + File.separator
							+ outputFileName));

			while (it.hasNext()) {

				Entry entry = (Entry) it.next();
				writer.write(entry.getSource() + "\t" + entry.getTarget()
						+ "\n");

			}

			writer.close();

		} catch (IOException e) {

			e.printStackTrace();
		}
	}

	//---------------------------------------------------------------------

	private boolean useExonerate() {

		String prop = System.getProperty("idmapping.use_exonerate");
		if (prop != null) {
			prop = prop.toLowerCase();
			if (prop.equals("yes") || prop.equals("true") || prop.equals("1")) {
				return true;
			}
		}

		return false;
	}

	//  -------------------------------------------------------------------------
	/**
	 * Prints argument string if the programms debug mode is enabled
	 *  
	 */
	private void debug(String s) {

		if (debug) {
			System.out.println(s);
		}

	}

	// -------------------------------------------------------------------------
}