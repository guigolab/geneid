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

package org.ensembl.probemapping;

import gnu.getopt.Getopt;
import gnu.getopt.LongOpt;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

import org.ensembl.datamodel.AffyArray;
import org.ensembl.datamodel.AffyFeature;
import org.ensembl.datamodel.ExternalDatabase;
import org.ensembl.datamodel.ExternalRef;
import org.ensembl.datamodel.Gene;
import org.ensembl.datamodel.Location;
import org.ensembl.datamodel.Transcript;
import org.ensembl.datamodel.impl.ExternalRefImpl;
import org.ensembl.driver.AdaptorException;
import org.ensembl.driver.ConfigurationException;
import org.ensembl.driver.Driver;
import org.ensembl.driver.DriverManager;
import org.ensembl.util.SimpleTimer;

/**
 * Program that maps MicroArray probesets to transcripts. The results are stored
 * as "XRefs" in an Ensembl database and a log file is produced that records why
 * all overlapping probeset-transcript pairs were or were not mapped.
 * 
 * <p>
 * The input data is loaded from source database(s) as AffyFeatures,AffyArrays,
 * AffyProbes and Transcripts and the mappings are saved as ExternalReferences
 * in the output databases.
 * </p>
 * 
 * @see #parseCommandLineInitParameters(String[]) for usage.
 * @see #map() for mapping algorithm.
 * @see org.ensembl.datamodel.AffyProbe
 * @see org.ensembl.datamodel.AffyFeature
 * @see org.ensembl.datamodel.AffArray
 * @see org.ensembl.datamodel.Transcript
 * @see org.ensembl.datamodel.ExternalRef
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 * 
 */
public class ProbeMapper {

	// DEFAULT PARAMETERS
	private static final double DEFAULT_THRESHOLD = 0.5;

	private static final int DEFAULT_TRANSCRIPTS_PER_PROBE_SET_THRESHOLD = 100;

	private static final int DEFAULT_DOWN_STREAM_FLANK = 2000;

	private static final int DEFAULT_MAX_TRANSCRIPTS_PER_COMPOSITE = 100;

	private static final String USER_DEFINED_FINISHED_LOCATION_FILE = "ignore.txt";

	private static final String DEFAULT_WORKING_DIRECTORY = ".";

	private static final String DEFAULT_LOG_FILENAME = "probeset2transcript.log";

	// PARAMETERS
	private File workingDirectory = new File(DEFAULT_WORKING_DIRECTORY);

	private int maxTranscriptsPerCompositeThreshold = DEFAULT_TRANSCRIPTS_PER_PROBE_SET_THRESHOLD;

	private Logger logger = Logger.getLogger(ProbeMapper.class.getName());

	private boolean useCache = false;

	private boolean verbose;

	private List locationFilterList = new ArrayList();

	private Location locationFilter;

	private Driver transcriptDriver;

	private Driver probeDriver;

	private Driver outputDriver;

	private Map transcript2ProbeSets = new HashMap();

	private int downStreamFlank = DEFAULT_DOWN_STREAM_FLANK;

	private double threshold = DEFAULT_THRESHOLD;

	private String logFilename = DEFAULT_LOG_FILENAME;

	// INTERNAL DATA STRUCTURES
	/**
	 * Affy features to map to transcripts;
	 */
	private MappableAffyFeature[] mappableAffyFeatures;

	/**
	 * Transcripts to be mapped to affy features and probesets.
	 */
	private MappableTranscript[] mappableTranscripts;

	/**
	 * Map from probset names to probesets.
	 * 
	 * @see ProbeSet
	 */
	private Map probeSets = new HashMap();

	/**
	 * probeset to transcript relationships. One will exist for every
	 * overlapping probeset and transcript.
	 */
	private List mappingStatuses = null;

	private Map xrefCache = new HashMap();

	/**
	 * Runs the application. Loads the data, performs the mapping and saves the
	 * results.
	 * 
	 * @param args
	 *            command line parameters, see commandLineInit(String[]) for
	 *            description.
	 * @see #parseCommandLineInitParameters(String[])
	 */
	public static void main(String[] args) throws ParseException, IOException {

		long time = System.currentTimeMillis();
		ProbeMapper app = new ProbeMapper();

		try {

			app.parseCommandLineInitParameters(args);
			app.run();

		} finally {

			if (app.verbose)
				System.out
						.println("ProbeMapper finished in = "
								+ (System.currentTimeMillis() - time) / 1000
								+ " secs.");
		}
	}

	/**
	 * For each of the locations specified load the probes and transcripts, map
	 * them and store the results.
	 */
	public void run() throws IOException {

		File file = new File(workingDirectory, logFilename);
		FileWriter logWriter = new FileWriter(file);
		
		
		for (int i = 0, n = locationFilterList.size(); i < n; i++) {

			// release memory and reduce search space by ensuring we only tyr to
			// map
			// probesets and transcripts
			// from the same location.
			mappingStatuses.clear();
			mappableTranscripts = null;
			mappableAffyFeatures = null;
			System.gc();

			locationFilter = (Location) locationFilterList.get(i);

			SimpleTimer timer = new SimpleTimer();

			timer.start();
			load();
			if (verbose)
				System.out.println("DURATION: load " + locationFilter + "= "
						+ timer.stop().getDurationInSecs() + "secs.");

			timer.start();
			mapAffyFeatures2Transcripts();
			markProbeSetsThatHitTooManyTranscripts(maxTranscriptsPerCompositeThreshold);
			findOverlappingProbeSetAndTranscripts();

			if (verbose)
				System.out.println("DURATION: mapping "
						+ timer.stop().getDurationInSecs() + "secs.");

			timer.start();
			// Possible optimisation: log and store could be parrallelised
			writeLog(logWriter);
			store();
			if (verbose)
				System.out.println("DURATION: store and log "
						+ timer.stop().getDurationInSecs() + "secs.");

		}
		logWriter.close();

	}

	/**
	 * Create MappingStatuses for all overlapping ProbeSet and Transcript pairs.
	 * 
	 * The relationship stores whether the pair are mapped or unmapped and if
	 * unmapped why not.
	 */
	private void findOverlappingProbeSetAndTranscripts() {

		mappingStatuses = new ArrayList();

		for (Iterator iter = probeSets.values().iterator(); iter.hasNext();) {

			ProbeSet ps = (ProbeSet) iter.next();
			// probesetSize * threshold
			int exonFlankThreshold = (int) Math.ceil((((AffyArray) ps.arrays
					.get(0)).getProbeSetSize() * threshold));

			// create a relationship for each probeset2transcript relationship
			for (Iterator iterator = ps.getOverlappingTranscripts().iterator(); iterator
					.hasNext();) {

				Object o = iterator.next();
				MappableTranscript t = (MappableTranscript) o;// iterator.next();

				int reverseStrandHitCount = 0;
				int exonFlankHitCount = 0;
				int intronHitCount = 0;

				for (Iterator probeIter = ps.affyFeatures.iterator(); probeIter
						.hasNext();) {

					final MappableAffyFeature af = (MappableAffyFeature) probeIter
							.next();

					final Location afLoc = af.affyFeature.getLocation();

					boolean strandedHit = false;
					boolean exonHit = false;

					// Note Location.overlaps() ignores strand.

					// probe hit either strand of transcript location?
					// Optimisation: call overlaps() before overlapSize()
					if (afLoc.overlaps(t.getLocation())
							&& afLoc.overlapSize(t.getLocation(), false) == 25) {

						// probe hit exon or flank?
						if (afLoc.overlapSize(t.getCDNALocation(), true) == 25)
							exonFlankHitCount++;

						// probe hit intron?
						else if (afLoc.overlapSize(t.getLocation(), true) == 25)
							intronHitCount++;

						else
							reverseStrandHitCount++;
					}

				}

				MappingStatus r = new MappingStatus(ps, exonFlankThreshold, t,
						exonFlankHitCount, intronHitCount,
						reverseStrandHitCount);
				mappingStatuses.add(r);
			}

		}

		if (verbose) {
			System.out.println("Found " + mappingStatuses.size()
					+ " overlapping probe sets and transcript pairs.");
			int t = 0;
			for (int i = 0, n = mappingStatuses.size(); i < n; i++)
				if (((MappingStatus) mappingStatuses.get(i)).isMapped())
					t += 1;
			System.out.println("Mapped " + t + " probe sets to transcripts.");
		}

	}

	/**
	 * Mark each probe set that maps to too many transcripts.
	 * 
	 * We do this because we do not want to map promiscuous probesets to any
	 * transcripts.
	 * 
	 * @param maxTranscriptsPerCompositeThreshold2
	 */
	private void markProbeSetsThatHitTooManyTranscripts(
			int maxTranscriptsPerCompositeThreshold2) {

		for (Iterator iter = probeSets.keySet().iterator(); iter.hasNext();) {
			ProbeSet ps = (ProbeSet) probeSets.get(iter.next());
			if (ps.getOverlappingTranscripts().size() > maxTranscriptsPerCompositeThreshold)
				ps.tooManyTranscripts = true;
		}

	}

	/**
	 * Writes log file containing composite2transcript relationships for all
	 * relationships where at lease one probe from a composite overlaps the
	 * genomic location where the transcript+flank is (ignores strand).
	 * 
	 * @throws IOException
	 */
	private void writeLog(FileWriter f) throws IOException {

		
		for (Iterator iter = mappingStatuses.iterator(); iter.hasNext();) {

			f.write(iter.next().toString());
			f.write("\n");

		}

		
	}

	/**
	 * Stores mapped probeset-transcripts xrefs and object_xrefs.
	 * 
	 * One xref is created for each probeset-array permutation. They are only
	 * stored if they do not already exist in the database. Each of these xrefs
	 * is linked to the transcipt via an object_xref entry.
	 * 
	 */
	private void store() throws AdaptorException {

		for (int k = 0, n = mappingStatuses.size(); k < n; k++) {

			MappingStatus status = (MappingStatus) mappingStatuses.get(k);

			if (!status.isMapped())
				continue;

			String probeSetName = status.probeSet.probeSetName;

			List arrays = status.probeSet.arrays;
			for (int i = 0, m = arrays.size(); i < m; i++) {

				AffyArray array = (AffyArray) arrays.get(i);
				ExternalDatabase xdb = array.getExternalDatabase();

				String cacheKey = xdb.getName() + "__" + probeSetName;
				ExternalRef xref = (ExternalRef) xrefCache.get(cacheKey);

				if (xref == null) {

					// xref might already be in the database
					// but NOT in the cache
					List xrefs = outputDriver.getExternalRefAdaptor().fetch(
							probeSetName);
					for (int j = 0; xref == null && j < xrefs.size(); ++j) {
						ExternalRef tmp = (ExternalRef) xrefs.get(j);
						if (tmp.getExternalDbId() == xdb.getInternalID())
							xref = tmp;
					}

					if (xref == null) {
						// store this probeset as an xref
						xref = new ExternalRefImpl(outputDriver);
						xref.setExternalDbId(xdb.getInternalID());
						xref.setPrimaryID(probeSetName);
						xref.setDisplayID(probeSetName);
						xref.setVersion("1");
						xref.setDescription(null);
						outputDriver.getExternalRefAdaptor().store(xref);
					}

					xrefCache.put(cacheKey, xref);
				
				}
				
				status.xrefs.add(xref);

				// store one object_xref for each mapped 
				// transcript-probeset-microarray permutation
				outputDriver.getExternalRefAdaptor()
				.storeObjectExternalRefLink(
						status.transcript.getInternalID(),
						ExternalRef.TRANSCRIPT,
						xref.getInternalID());
			}
		}
	}

	/**
	 * 
	 */
	private void mapAffyFeatures2Transcripts() throws IOException {

		int mappedCount = 0;
		int comparisonCount = 0;

		final int nMappableTranscripts = mappableTranscripts.length;
		final int nMappableAffyFeatures = mappableAffyFeatures.length;

		if (true) {

			// We sort the transcripts and affy features by genomic location
			// before we start the search for overlapping pairs as this allows
			// us to avoid many uncessary location comparisons.

			final int nts = mappableTranscripts.length;
			final int nafs = mappableAffyFeatures.length;
			final MappableTranscript[] ts = mappableTranscripts;
			final MappableAffyFeature[] afs = mappableAffyFeatures;

			// sort by location
			Arrays.sort(ts);
			Arrays.sort(afs);

			int afiStart = 0;
			for (int ti = 0; ti < nts; ti++) {

				final MappableTranscript t = ts[ti];
				final Location tLoc = t.getLocation();
				final boolean updateAFiStart = ti + 1 < nts
						&& !tLoc.overlaps(ts[ti + 1].getLocation(), false);

				for (int afi = afiStart; afi < nafs; afi++) {

					final MappableAffyFeature af = afs[afi];
					final Location afLoc = af.location;
					// We ignore strand in this overlap check otherwise
					// we might break out of loop below prematurely.
					final boolean overlap = tLoc.overlaps(afLoc, false);

					comparisonCount++;

					if (overlap) {

						final int overlapSize = t.getCDNALocation()
								.overlapSize(afLoc, true);
						if (overlapSize == afLoc.getLength()) {
							mappedCount += 1;
							af.addTranscript(t);
						}
					} else if (afLoc.compareTo(tLoc) == 1) {
						// don't compare the transcript to the remaining affy
						// features
						// because they are either after it on the seq region or
						// on a different
						// seq region.
						break;
					}

					// Start comparing the next transcript with the next affy
					// feature
					if (updateAFiStart)
						afiStart = afi + 1;
				}
			}

		} else {
			// OLD ALGORITHM m*n

			for (int j = 0; j < nMappableAffyFeatures; j++) {

				final MappableAffyFeature feature = mappableAffyFeatures[j];
				final Location loc = feature.affyFeature.getLocation();

				for (int i = 0; mappableTranscripts != null
						&& i < mappableTranscripts.length; i++) {

					final MappableTranscript t = mappableTranscripts[i];

					// This next statment is an optimization which acts
					// as a fast
					// filter.
					// It compares one location *node* against another
					// *node*in order
					// to
					// decide whether to do the more computationally
					// expensive
					// overlaps() call later.
					comparisonCount++;
					if (t.getLocation().overlaps(loc, true)) {

						final int overlap = loc.overlapSize(
								t.getCDNALocation(), true);
						if (overlap == loc.getLength()) {
							mappedCount += 1;
							feature.addTranscript(t);

						}
					}
				}
			}
		}
		if (verbose)
			System.out
					.println("Found "
							+ mappedCount
							+ " raw mappings between affy_features and transcripts in location."
							+ locationFilter + ". (" + comparisonCount
							+ " comparisons).");

	}

	private void error(String msg) {
		System.err.println("Error: " + msg);
		System.out.println(usage());
		System.exit(0);
	}

	/**
	 * Load the Probes and the Transcripts from the source database(s).
	 */
	public void load() throws AdaptorException {

		// load the misc features first so that the misc feature related
		// validation
		// runs first.
		// We want to do this because a) there is more likely to be problem with
		// them than transcripts (e.g. no corresponding external database
		// entries)
		// and b) they load quicker than transcripts.

		Thread lmf = new Thread() {
			public void run() {
				try {
					loadAffyFeatures();
				} catch (AdaptorException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		};

		Thread lt = new Thread() {
			public void run() {
				try {
					loadTranscripts();
				} catch (AdaptorException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		};

		lmf.start();
		lt.start();

		try {
			lmf.join();
			lt.join();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

	}

	private void loadAffyFeatures() throws AdaptorException {

		List tmp = probeDriver.getAffyFeatureAdaptor()
				.fetchUniqueProbeAndLocation(locationFilter);
		mappableAffyFeatures = new MappableAffyFeature[tmp.size()];
		for (int i = 0; i < mappableAffyFeatures.length; i++) {

			AffyFeature af = (AffyFeature) tmp.get(i);
			String psName = af.getProbeSetName();
			ProbeSet ps = (ProbeSet) probeSets.get(psName);

			// Lazy creation of Probeset using info from one of it's
			// probes.
			if (ps == null) {
				ps = new ProbeSet(psName, probeDriver.getAffyProbeAdaptor()
						.fetch(af.getProbeInternalID())
						.getAffyArraysContainingThisProbe());
				probeSets.put(psName, ps);
			}

			mappableAffyFeatures[i] = new MappableAffyFeature(af, ps);
			ps.addMappableAffyFeature(mappableAffyFeatures[i]);
		}

		if (verbose)
			System.out.println("Loaded " + tmp.size() + " Affy Features and "
					+ probeSets.size() + " Probesets.");
	}

	private void loadTranscripts() throws AdaptorException {

		List buf = new ArrayList();

		for (Iterator iter = transcriptDriver.getGeneAdaptor().fetchIterator(
				locationFilter, true); iter.hasNext();) {
			Gene g = (Gene) iter.next();
			List ts = g.getTranscripts();
			for (int j = 0; j < ts.size(); j++) {
				Transcript t = (Transcript) ts.get(j);
				buf.add(new MappableTranscript(t, downStreamFlank));
			}
		}

		mappableTranscripts = new MappableTranscript[buf.size()];
		buf.toArray(mappableTranscripts);

		if (verbose)
			System.out.println("Loaded " + mappableTranscripts.length
					+ " Transcripts");
	}

	/**
	 * Initialise the instance from command line parameters. Parses the command
	 * line parameters and validates them before initiailising the instance. Run
	 * without parameters or with -h or --help to see usage.
	 * 
	 * @see org.ensembl.driver.plugin.standard.MySQLDriver for db config file
	 *      specification
	 * @param commandLineArgs
	 *            command line arguments, see description for more details.
	 */
	public void parseCommandLineInitParameters(String[] commandLineArgs)
			throws ConfigurationException, ParseException, AdaptorException {

		String outputDriverFilepath = null;
		String transcriptDriverFilepath = null;
		String probeDriverFilepath = null;

		boolean showHelp = commandLineArgs.length == 0;
		verbose = false;

		LongOpt[] longopts = new LongOpt[] {
				new LongOpt("help", LongOpt.NO_ARGUMENT, null, 'h'),
				new LongOpt("verbose", LongOpt.NO_ARGUMENT, null, 'v'),
				new LongOpt("affy-db", LongOpt.REQUIRED_ARGUMENT, null, 'a'),
				new LongOpt("transcript-db", LongOpt.REQUIRED_ARGUMENT, null,
						't'),
				new LongOpt("output-db", LongOpt.REQUIRED_ARGUMENT, null, 'o'),
				new LongOpt("location", LongOpt.REQUIRED_ARGUMENT, null, 'l'),
				new LongOpt("log-file", LongOpt.NO_ARGUMENT, null, 'L'),
				new LongOpt("down-stream-flank", LongOpt.NO_ARGUMENT, null, 'f'),
				new LongOpt("threshold", LongOpt.REQUIRED_ARGUMENT, null, 'T'),
				new LongOpt("max-transcripts-per-composite",
						LongOpt.REQUIRED_ARGUMENT, null, 'n'),
				new LongOpt("dir", LongOpt.REQUIRED_ARGUMENT, null, 'd') };

		Getopt g = new Getopt("ProbeToTranscriptMappingApplication",
				commandLineArgs, "hva:t:o:l:L:f:T:n:d:", longopts, false);
		int c;
		String arg;
		while ((c = g.getopt()) != -1) {
			switch (c) {

			case 'h':
				showHelp = true;
				break;

			case 'v':
				verbose = true;
				break;

			case 'o':
				outputDriverFilepath = g.getOptarg();
				break;

			case 'a':
				probeDriverFilepath = g.getOptarg();
				break;

			case 't':
				transcriptDriverFilepath = g.getOptarg();
				break;

			case 'l':
				locationFilterList.clear();
				locationFilterList.add(new Location(g.getOptarg()));
				break;

			case 'f':
				downStreamFlank = Integer.parseInt(g.getOptarg());
				break;

			case 'T':
				int percentage = Integer.parseInt(g.getOptarg());
				threshold = percentage / 100.0;
				break;

			case 'n':
				maxTranscriptsPerCompositeThreshold = Integer.parseInt(g
						.getOptarg());
				break;

			case 'd':
				workingDirectory = new File(g.getOptarg());
				break;

			case 'L':
				logFilename = g.getOptarg();
				break;
			}
		}

		if (showHelp) {
			System.out.println(usage());
			System.exit(0);
		}

		// connect to separate input and output databases if specified
		// by user
		if (outputDriverFilepath != null)
			outputDriver = loadDriver("output", outputDriverFilepath);
		if (probeDriverFilepath != null)
			probeDriver = loadDriver("probe", probeDriverFilepath);
		if (transcriptDriverFilepath != null)
			transcriptDriver = loadDriver("transcript",
					transcriptDriverFilepath);

		// use "default" driver where necessary
		if (g.getOptind() < commandLineArgs.length) {
			String filepath = commandLineArgs[g.getOptind()];
			if (outputDriverFilepath == null)
				outputDriver = loadDriver("default", filepath);
			if (probeDriverFilepath == null)
				probeDriver = loadDriver("default", filepath);
			if (transcriptDriverFilepath == null)
				transcriptDriver = loadDriver("default", filepath);
		}

		if (outputDriver == null)
			error("Output driver is not set.");
		if (!outputDriver.isConnected())
			error("Cannot connect to Output database");

		if (downStreamFlank < 0)
			error("Down stream flank must be >= 0");

		if (probeDriver == null)
			error("Probe driver is not set.");
		if (!probeDriver.isConnected())
			error("Cannot connect to probe database.");
		if (transcriptDriver == null)
			error("Transcript driver is not set.");
		if (!transcriptDriver.isConnected())
			error("Cannot connect to Transcript database.");

		if (locationFilterList.size()==0) {
			// by default load all data in chromosome coordinate system
			for(Location head = transcriptDriver.getLocationConverter().fetchComplete(new Location("chromosome"));
			head!=null;) {
				Location node = head;
				head = head.next();
				node.setNext(null);
				// bug in fetchComplete() causes location without chromosome to be returned, ignore it.
				if (node.getSeqRegionName()!=null)
					locationFilterList.add(node);
			}
			// put the shortest items at the beginning to help spot bugs quickly during testing
			Collections.sort(locationFilterList, new Comparator() {
				public int compare(Object o1, Object o2) {
					Location l1 = (Location) o1;
					Location l2 = (Location) o2;
					return l1.getLength()-l2.getLength();
				}
			});
			}

		
		if (verbose) {
			System.out.println("ProbeMapper Configuration:");
			System.out.println("==========================");
			System.out.println("Transcript Database: " + transcriptDriver);
			System.out.println("Probe Database: " + probeDriver);
			System.out.println("Output Database: " + outputDriver);
			System.out.println("Location batches: " + locationFilterList);
			System.out.println("Down stream flank: " + downStreamFlank);
			System.out.println("Mapping overlap threshold: "
					+ NumberFormat.getPercentInstance().format(threshold));
			System.out.println("Max transcripts per composite: "
					+ maxTranscriptsPerCompositeThreshold);
			System.out.println("Working directory: "
					+ workingDirectory.getAbsolutePath());
			System.out.println("==========================");
			System.out.flush();
		}

	}

	private Driver loadDriver(String description, String filepath) {
		Driver d = null;
		filepath = new File(workingDirectory, filepath).getAbsolutePath();
		try {
			d = DriverManager.load(filepath);
		} catch (ConfigurationException e) {
			error("Cannot connect to " + description
					+ " database specified in file '" + filepath + "'\n"
					+ e.getMessage());
		}
		if (d == null) {
			error("Failed to initialise " + description
					+ " database specified in file " + filepath);
		}
		return d;
	}

	/**
	 * Usage.
	 * 
	 * @return string containing command line usage for this program.
	 */
	public String usage() {
		String usage = ""
				+ "Maps MicroArray probesets and transcripts from ensembl database(s). "
				+ "\nThe results are stored in an ensembl database as xrefs and "
				+ "\na mapping log file (probeset2transcript.log) is written to the "
				+ " working directory."
				+ "\nUsage: "
				+ "\n  ProbeMapper [OPTIONS] [CONFIG_FILE]"
				+ "\n"
				+ "\nCONFIG_FILE specifies the database to retrieve the probe and transcript data from "
				+ "the database to write results to. Each of these datases can be different and"
				+ "specified using the --affy-probe-db (-a), --transcript-db (-t) and --output_db (-o) options."
				+ "\n"
				+ "\nWhere options are:"
				+ "\n -v                                                        Print verbose information."
				+ "\n     --verbose"
				+ "\n -T THRESHOLD                                              Minimum percentage of composite to transcript hits for mapping. e.g. 60, 75."
				+ "\n     --threshold=THRESHOLD"
				+ "\n -f  FLANK                                                 Down stream flank in bases. Default is "
				+ DEFAULT_DOWN_STREAM_FLANK
				+ "\n     --down-stream-flank FLANK"
				+ "\n -n  MAX_TRANSCRIPTS_PER_COMPOSITE_THRESHOLD               Max number of transcripts allowed per transcript (per sequence region). Default is "
				+ DEFAULT_MAX_TRANSCRIPTS_PER_COMPOSITE
				+ "."
				+ "\n     --max-transcripts-per-composite MAX_TRANSCRIPTS_PER_COMPOSITE_THRESHOLD"
				+ "\n -l LOCATION                                               Only map probes and transcripts in this location. e.g. chromosome:1:20m-21m"
				+ "\n     --location=LOCATION"
				+ "\n -a AFFY_DATABSE_CONFIG_FILE                               File specifying db containing affymetrix probes."
				+ "\n      --affy-probe-db=AFFY_DATABASE_CONFIG_FILE"
				+ "\n -t TRANSCRIPT_DATABASE_CONFIG_FILE                        File specifying db containing transcript."
				+ "\n     --transcript-db=TRANSCRIPT_DATABASE_CONFIG_FILE"
				+ "\n -o OUTPUT_DATABSEB_CONFIG_FILE                            File specifying db to save results to."
				+ "\n     --output-db OUTPUT_DATABSE_CONFIG_FILE"
				+ "\n -d WORKING_DIR                                            Working directory. All file paths are relative to this directory. Default is \".\"."
				+ "\n     --dir=WORKING_DIR"
				+ "\n -L LOG_FILE                                               Log file. Default is "
				+ DEFAULT_LOG_FILENAME
				+ "."
				+ "\n     --log-file LOG_FILE"
				+ "\n -h                                                        Print this help message."
				+ "\n     --help" + "\n";

		return usage;
	}

	/**
	 * Probed transcripts to be mapped.
	 * 
	 * @return zero or more probed transcript.
	 */
	public MappableTranscript[] getMappableTranscripts() {
		return mappableTranscripts;
	}

	/**
	 * @return Returns the threshold.
	 */
	public double getThreshold() {
		return threshold;
	}

	/**
	 * @param threshold
	 *            The threshold to set.
	 */
	public void setThreshold(double threshold) {
		this.threshold = threshold;
	}
}
