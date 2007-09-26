	package gphase.io.gtf;

import gphase.Analyzer;
import gphase.algo.ASAnalyzer;
import gphase.graph.SpliceGraph;
import gphase.io.DefaultIOWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Array;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.IntVector;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.sql.Date;
import java.text.Collator;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//import org.apache.commons.collections.BidiMap;


import com.sun.org.apache.xalan.internal.xsltc.compiler.util.FilterGenerator;

import sun.font.GlyphLayout.GVData;


public class GTFChrReader extends GTFWrapper {
	static boolean sort= false;
	static String usage= "GTFChrReader [options] <inputFile>\n\n"+
		"where options may be\n" +
		"-sort\t sort according to chr, transcriptID, start\n\n" +
		"micha, 2007";
	public static final String NORMED_FILE_TAG= "norman";
	public static final int MAX_GTF_LINE_LENGTH= 1000;
	public static final String[] DEFAULT_CHROMOSOME_FILTER= new String[] {	
			"^CHRM.*$", 
			"^M$",
			"^MT$",
			"^U$",
			"^UN$",
			"^CHRU$",	
			"^CHRUN$",	
			".*UNKNOWN.*",
			".*_RANDOM.*", 
			".*_HAP.*",	//chr5_h2_hap1
			"^\\d{1,}H$",	// droso, hetero?
			"^CHR\\d{1,}H$",
			"^X{1}H$",
			"^CHRXH$",
			"^Y{1}H$",
			"^CHRYH$",
			"^NT_.+$"
	};	// default chromosomes filtered off
	public String transcript_id_tag= GTFObject.TRANSCRIPT_ID_TAG;
	
	
	protected InputStream inputStream= null;
	Species species= null;
	int linesRead= 0;
	long bytesRead= 0;
	int lastPerc= 0;
	Vector readChr= null;
	int readTranscripts= 0;
	int readExons= 0;
	
	int readGenes= 0;
	int readChrs= 0;
	Vector skippedChr= new Vector();
	Vector skippedFeatures= new Vector();
	int skippedObjects= 0;
	Vector skippedTranscripts= new Vector();
	String[] filtChrIDs= null;
	HashMap filtSomeIDs= null;
	boolean[] filtSomeIDSuccess= null;
	String[] noIDs= DEFAULT_CHROMOSOME_FILTER;

	DirectedRegion[] filtRegs= null;
	String[] filtGeneIDs= null;
	String[] filtTrptIDs= null;
	/**
	 * 
	 */
	String[] readFeatures= new String[] {"exon", "CDS"};
	boolean silent= false;
	boolean chromosomeWise= true;
	/**
	 * 
	 */
	boolean readGTF= false;
	/**
	 * 
	 */
	boolean readGene= true;
	Gene[] genes= null;
	
	boolean clustered= false;
	
	// !!! doesnt find the following mutex introns
	// AC051649.8-006	AC051649.8-010	1974632 1975281 1975281 1975505 	1^2- , 3^4-
	// 1bp exonic overlap between the introns
	public Graph getGraph(boolean encode) {
		if (gtfObj== null) {
			try {
				read();
			} catch (Exception e) {
				e.printStackTrace(); 
			}
		}
		return assemble_old(encode);		// <===== check ENCODE here !!!
	}
	
	public static void main(String[] args) {
		if (args.length< 1) {
			System.out.println(usage);
			System.exit(0);
		}
		
		File f= null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-sort")) {
				sort= true;
				continue;
			}
			f= new File(args[i]);
		}
		
		if (f== null|| !f.exists())  {
			System.out.println("input file invalid");
			System.exit(0);
		}
		
		GTFChrReader reader= new GTFChrReader(f.getAbsolutePath());
		if (sort) {
			reader.reformatFile();
		}
	}
	
	public GTFChrReader(String absFName) {
		super(absFName);	
		readChr= new Vector();
		//initSpecies();
	}
	
	void initSpecies() {
		StringTokenizer st= new StringTokenizer(fName, "_.");
		String speName= null, annotationVersion= null, genomeVersion= null;
		int speNb= -1;
		while (st.hasMoreTokens()) {
			String s= st.nextToken();
			int tmpSpNb= Species.getSpeciesNumber(s);
			if (tmpSpNb>= 0) {
				if (speName== null)
					speName= s;
				else {
					if (!speName.equalsIgnoreCase(s))
						System.err.println("Conflicting info on species in file name.");
				}
				if (speNb>= 0&& tmpSpNb!= speNb)
					System.err.println("Conflicting info in file name: "+s+" <> "+Species.SP_NAMES_COMMON[speNb]);
				continue;
			}
			
			tmpSpNb= Species.getSpeciesNumberForGenomeVersion(s);
			if (tmpSpNb>= 0) {
				if (speNb>= 0&& tmpSpNb!= speNb)
					System.err.println("Conflicting info in file name: "+s+" <> "+Species.SP_NAMES_COMMON[speNb]);
				speNb= tmpSpNb;
				speName= Species.SP_NAMES_COMMON[speNb];
				if (genomeVersion== null)
					genomeVersion= s;
				else {					
					if (Species.getGenomeVerNb(genomeVersion)!= Species.getGenomeVerNb(s))
						System.err.println("Conflicting info on genome in file name: "+genomeVersion+" <> "+s);
				}
				continue;
			}

//			tmpSpNb= Species.getSpeciesNumberForAnnotation(s);
//			if (tmpSpNb>= 0) {
//				if (speNb>= 0&& tmpSpNb!= speNb)
//					System.err.println("Conflicting info in file name: "+s+" <> "+Species.SP_NAMES_COMMON[speNb]);
//				speNb= tmpSpNb;
//				speName= Species.SP_NAMES_COMMON[speNb];
//				String tmpGenomeVer= Species.getGenomeVersionForAnnotation(s);
//				if (genomeVersion!= null&& Species.getGenomeVerNb(genomeVersion)!= Species.getGenomeVerNb(s))
//					System.err.println("Conflicting info on genome in file name: "+ tmpGenomeVer+" <> "+genomeVersion);
//				else
//					genomeVersion= tmpGenomeVer;
//				if (annotationVersion== null)
//					annotationVersion= s;
//				else {
//					if (Species.getGenomeVerNbForAnnotation(genomeVersion)!= Species.getGenomeVerNbForAnnotation(s))
//						System.err.println("Conflicting info on annotation in file name: "+annotationVersion+" <> "+s);
//				}
//				continue;
//			}
		}
		
		if (speName== null) {
			if (!silent)
				System.out.println("No valid species name found in file name, guessing \'human\'.");
			species= new Species("human");
		} else 
			species= new Species(speName);
		
		if (genomeVersion== null) {
			genomeVersion= species.getDefaultGenomeVersion();
			if (!silent)
				System.out.println("No valid build found in file name, guessing "+genomeVersion+".");
			species.setGenomeVersion(genomeVersion);
		} else 
			species.setGenomeVersion(genomeVersion);
		
		if (annotationVersion!= null)
			species.setAnnotationVersion(annotationVersion);
	}
	
	public GTFChrReader(InputStream i) {
		System.err.println("init with inputStream not allowed.");
	}
	
	GTFObject createGTFObject(){
		return new France();
	}

	public France[] getFranceObj() {
		return (France[]) getGtfObj();
	}
	
	Vector getClusters(HashMap transHash) {
		Vector trans= new Vector(transHash.values());
		HashMap clusters= new HashMap();	// transcriptID x (vector x vector)
		for (int i = 0; i < trans.size(); i++) {
			for (int j = i+ 1; j < trans.size(); j++) {
				if (overlap((Vector) trans.elementAt(i), (Vector) trans.elementAt(j))) {
					Vector v1= (Vector) trans.elementAt(i);
					Vector v2= (Vector) trans.elementAt(j);
					Vector vc1= (Vector) clusters.remove(((France) v1.elementAt(0)).getTranscriptID());
					// remove from all other transcript ids in the cluster in hash
					if (vc1== null) {
						vc1= new Vector();
						vc1.add(v1);
					}
					Vector vc2= (Vector) clusters.remove(((France) v2.elementAt(0)).getTranscriptID());
					if (vc2== null) {
						vc2= new Vector();
						vc2.add(v2);
					}
					Vector result= (Vector) Arrays.merge(vc1, vc2);	// doppelte eintraege moeglich??
					for (int k = 0; k < result.size(); k++) 
						clusters.put(((France) ((Vector) result.elementAt(k)).elementAt(0)).getTranscriptID(),
								result);
				}
					
			}
		}
		return trans;
	}

	boolean overlap (Vector trans1, Vector trans2) {
//		for (int i = 0; i < trans1.size(); i++) {
//			France f1= (France) trans1.elementAt(i);
//			if (!f1.isExon())
//				continue;
//			for (int j = 0; j < trans2.size(); j++) {
//				France f2= (France) trans2.elementAt(j);
//				if (!f2.isExon())
//					continue;
//				DefaultDirectedRegion d1= new DefaultDirectedRegion(f1.getStrand(), f1.getStart(), f1.getEnd());
//				DefaultDirectedRegion d2= new DefaultDirectedRegion(f1.getStrand(), f2.getStart(), f2.getEnd());
//				if (d1.overlaps(d2))
//					return true;
//			}
//		}
		return false;
	}
	
	static HashMap getChromosomes(HashMap transGTF) {
		Iterator iter= transGTF.values().iterator();
		HashMap chrHash= new HashMap();
		while (iter.hasNext()) {
			Vector gtfsVec= (Vector) iter.next();
			GTFObject o= (GTFObject) (gtfsVec).elementAt(0);
			HashMap tHash= (HashMap) chrHash.remove(o.getChromosome());
			if (tHash== null)
				tHash= new HashMap();
			Vector v= (Vector) tHash.remove(o.getTranscriptID());
			if (v== null)
				v= new Vector();
			for (int i = 0; i < gtfsVec.size(); i++) 
				v.add(gtfsVec.elementAt(i));
			tHash.put(o.getTranscriptID(), v);
			chrHash.put(o.getChromosome(), tHash);
		}
		return chrHash;
	}
	
	private void checkMegaClusters(Transcript[] trans) {
		java.util.Arrays.sort(trans, new DirectedRegion.PositionComparator());	// sort ascending
		for (int i = trans.length- 1; i >= 0; --i) {
			System.out.println(trans[i].getTranscriptID()+ "\t"+ trans[i].getLength()+"\t"+trans[i].getExons().length);
		}
	}
	
	public boolean isEOF() {
		File f= new File(fPath+ File.separator+ fName);
		
		return (bytesRead== f.length());  
	}
	

	Graph assemble_old(boolean encode) {
		
		Species spec= new Species("human");
		spec.setAnnotationVersion("17");
	
			// cluster
		HashMap hash= getGroups("transcript_id", getGtfObj());	// cluster for genes?
		HashMap chrHash= getChromosomes(hash);
		
			// construct transcripts
		Collection co= ((Collection) chrHash.keySet());
		String[] keys= new String[co.size()];
		Iterator iter= co.iterator();
		int x= 0;
		while(iter.hasNext()) 
			keys[x++]= (String) iter.next();
		
		HashMap chr2Hash= new HashMap(chrHash.size());
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap tHash= (HashMap) chrHash.get(chrID);
			Collection co2= ((Collection) tHash.keySet());
			String[] tkeys= new String[co2.size()];
			Iterator iter2= co2.iterator();
			x= 0;
			while (iter2.hasNext())					
				tkeys[x++]= (String) iter2.next();
			HashMap t2Hash= new HashMap(tHash.size());
			chr2Hash.put(chrID, t2Hash);
			for (int j = 0; j < tkeys.length; j++) {	// transcripts
				String tID= tkeys[j];
				GTFObject[] gtfs= (GTFObject[]) Arrays.toField(tHash.get(tID));	// gtf entries for 1 transcript
				France ff= (France) gtfs[0];
				if (encode&& !ff.getSource().contains("VEGA"))
					continue;
				Transcript transcript= new Transcript(tID);
				transcript.setStrand(ff.getStrand());
				for (int k = 0; k < gtfs.length; k++) {		// exons 
					France f= (France) gtfs[k];
					if (f.isExon()) { 
						transcript.updateBoundaries(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
						transcript.setHUGO(gtfs[k].getAttribute("gene_id"));
					}
				}
				t2Hash.put(tID, transcript);	// fill tHash with transcripts
			}
			
		}
		
			// check mega clusters
//		HashMap[] maps= (HashMap[]) Arrays.toField(chr2Hash.values());
//		Vector v= new Vector(maps.length* 100);
//		for (int i = 0; i < maps.length; i++) 
//			v.addAll(maps[i].values());
//		checkMegaClusters(((Transcript[]) Arrays.toField(v)));
		
			// cluster
		HashMap gHash= new HashMap();
		Comparator compi= new DirectedRegion.PositionComparator();
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap t2Hash= (HashMap) chr2Hash.get(chrID);
			Object[] transcripts= t2Hash.values().toArray();
			java.util.Arrays.sort(transcripts, compi);
			Transcript[] t= new Transcript[transcripts.length];
			for (int j = 0; j < t.length; j++) 
				t[j]= (Transcript) transcripts[j];
			Transcript[][] loci= clusterTranscripts(t);
			for (int j = 0; j < loci.length; j++) {
				String gID= Gene.getUniqueID();
				Gene locus= new Gene(spec, gID);
				locus.setStrand(loci[j][0].getStrand());
				locus.setChromosome(chrID);
				for (int k = 0; k < loci[j].length; k++) { // transcripts
					loci[j][k].setGene(locus);
					Vector v= (Vector) ((HashMap) chrHash.get(chrID)).get(loci[j][k].getTranscriptID());
					for (int m = 0; m < v.size(); m++) {
						France f= (France) v.elementAt(m);
						if (f.isExon()) {
							loci[j][k].addExon(new Exon(loci[j][k], f.getExonID(), f.getStart(), f.getEnd()));
						} else if (f.isCoding()) {	// for extending translation
							loci[j][k].addCDS(f.getStart(), f.getEnd());
						}
					}
						// iterate again, all exons have to be there, order in input file not reliable
					for (int m = 0; m < v.size(); m++) {
						France f= (France) v.elementAt(m);
						if (encode&& f.isCoding()) {		// in not-encode no exon ids..
							Exon e= loci[j][k].getExon(f.getAttribute("exon_id"));	// 2nd iteration necessary, cds st b4 exon in file
							e.extendStartCDS(f.getStart());
							e.extendEndCDS(f.getEnd());
							e.setFrame(f.getFrame());
						}
					}
					locus.addTranscript(loci[j][k]);
				}
				gHash.put(gID, locus);
			}
		}		
		
			// build graph
		iter= gHash.values().iterator();
		Graph g= new Graph();
		g.addSpecies(spec);
		while (iter.hasNext()) 
			g.addGene((Gene) iter.next());
		return g;
	}

	/**
		 * depends on the order of gtf-objects from the input for creating transcripts
		 * to create from the outside
		 * @param encode
		 * @return
		 */
		public Gene[] assemble(GTFObject[] gtfObs) {
	
			// build transcripts according to order from the input
			String lastTID= "";
			int lastStrand= 0;
			Transcript trpt= null;
			Vector trptV= new Vector();
			Gene lastGene= null;
			Vector geneV= new Vector();
			String lastChrID= null;
			Comparator startCompi= new AbstractRegion.StartComparator();
			Gene[] geneReg= null;
			boolean merged= true;	// for skiping insertion of init null trpt/gene
			for (int x = 0; x < gtfObs.length; x++) {		// exons
				
				// skip mRNA, gene...
				if (!gtfObs[x].getFeature().equals("exon")&& !gtfObs[x].getFeature().equals("CDS"))
					continue;
				
					// get chrom, GID, TID
				String chrID= gtfObs[x].getSeqname();
				String tid= gtfObs[x].getTranscriptID();
				String gid= gtfObs[x].getGeneID();
				int strand= gtfObs[x].getStrand();
				int start= gtfObs[x].getStart();
				int end= gtfObs[x].getEnd();
				if (tid== null) {
					System.err.println("No TID found.");
					continue;
				}
	
					// new transcript?
				if (!tid.equals(lastTID)|| !chrID.equals(lastChrID)) {
						// save gene
					if (!merged) {	// if last transcript's gene couldnt be merged with anything
						if (geneReg== null)
							geneReg= new Gene[] {trpt.getGene()};
						else {
							int p= java.util.Arrays.binarySearch(geneReg, trpt.getGene(), startCompi);
							geneReg= (Gene[]) Arrays.insert(geneReg, trpt.getGene(), p);
						}
					}
					merged= false;
					lastChrID= chrID;
					Gene ge= new Gene(Gene.getUniqueID());
					ge.setStrand(strand);
					ge.setChromosome(chrID);
					ge.setSpecies(species);
					trpt= new Transcript(tid);
					lastTID= tid;
					trpt.setStrand(strand);
					trpt.setChromosome(chrID);
					trpt.setSource(gtfObs[x].getSource());
					ge.addTranscript(trpt);
				}
										
				if (gtfObs[x].getFeature().equalsIgnoreCase("exon")) {
					String exonID= gtfObs[x].getExonID();
					Exon e= new Exon(trpt, exonID, start, end);
					e.setChromosome(chrID);
					e.setStrand(strand);
					
						// check if gene merge, not necessary when transcript starts are sorted ascending
					if (geneReg!= null) {
						int p= java.util.Arrays.binarySearch(geneReg, e, startCompi);
						if (p>= 0) {
							if (geneReg[p]!= trpt.getGene()) {
								if (merged) {
									trpt.getGene().merge(geneReg[p]);
									geneReg= (Gene[]) Arrays.remove(geneReg, p);
								} else {
									geneReg[p].merge(trpt.getGene());
									merged= true;
								}
							}
						} else {	// between two
							p= -(p+1);	// insertion point
							int q= p-1;
							if (q>= 0&& geneReg[q]!= trpt.getGene()&&
									e.overlaps(geneReg[q])) { 
								if (merged) {	// if already merged, one has to be removed from array
									trpt.getGene().merge(geneReg[q]);
									geneReg= (Gene[]) Arrays.remove(geneReg, q);
								} else
									geneReg[q].merge(trpt.getGene());
								merged= true;
							}
							q= p;	//+1; NO, insertion point is the upper neighbor
							if (q< geneReg.length&& geneReg[q]!= trpt.getGene()&&
									e.overlaps(geneReg[q])) { 
								if (merged) {
									trpt.getGene().merge(geneReg[q]);
									geneReg= (Gene[]) Arrays.remove(geneReg, q);
								} else
									geneReg[q].merge(trpt.getGene());
								merged= true;
							}
						}
					}
					
					trpt.addExon(e);
					
				} else if (gtfObs[x].getFeature().equalsIgnoreCase("CDS")) 
					trpt.addCDS(start, end);			
			}	// end for
			
			if (!merged) {
				if (geneReg== null)
					geneReg= new Gene[] {trpt.getGene()};
				else {
					int p= java.util.Arrays.binarySearch(geneReg, trpt.getGene(), startCompi);
					geneReg= (Gene[]) Arrays.insert(geneReg, trpt.getGene(), p);
				}
			}
			
			return geneReg;
		}

	/**
	 * depends on the order of gtf-objects from the input for creating transcripts
	 * to create from the outside
	 * @param encode
	 * @return
	 */
	public static Species assemble(GTFObject[] gtfObs, String speName, String buildVersion) {
	
		Species spe= new Species(speName);
//		Gene[] ge= assemble(gtfObs);
//		spe.setAnnotationVersion(buildVersion);
//		for (int i = 0; i < ge.length; i++) {
//			spe.addGene(ge[i]);
//			ge[i].setSpecies(spe);
//		}
		
		return spe;
	}
	
	public String getLastReadChr() {
		if (readChr== null)
			return null;
		return (String) readChr.elementAt(readChr.size()- 1);
	}

	/**
	 * doesnt take transcript_id order on the input gtfObs into account
	 * problem f.i. with C.elegans -> same tID on different chromosomes
	 * 
	 * @param encode
	 * @return
	 */
	public static Graph assemble_old(boolean encode, GTFObject[] gtfObs) {
		
		Species spec= new Species("human");
		spec.setAnnotationVersion("17");

			// cluster
		HashMap hash= getGroups("transcript_id", gtfObs);	// cluster for genes?
		HashMap chrHash= getChromosomes(hash);
		
			// construct transcripts
		Collection co= ((Collection) chrHash.keySet());
		String[] keys= new String[co.size()];
		Iterator iter= co.iterator();
		int x= 0;
		while(iter.hasNext()) 
			keys[x++]= (String) iter.next();
		
		HashMap chr2Hash= new HashMap(chrHash.size());
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap tHash= (HashMap) chrHash.get(chrID);
			Collection co2= ((Collection) tHash.keySet());
			String[] tkeys= new String[co2.size()];
			Iterator iter2= co2.iterator();
			x= 0;
			while (iter2.hasNext())					
				tkeys[x++]= (String) iter2.next();
			HashMap t2Hash= new HashMap(tHash.size());
			chr2Hash.put(chrID, t2Hash);
			for (int j = 0; j < tkeys.length; j++) {	// transcripts
				String tID= tkeys[j];
				GTFObject[] gtfs= (GTFObject[]) Arrays.toField(tHash.get(tID));	// gtf entries for 1 transcript
				France ff= (France) gtfs[0];
				if (encode&& !ff.getSource().contains("VEGA"))
					continue;
				Transcript transcript= new Transcript(tID);
				transcript.setStrand(ff.getStrand());
				for (int k = 0; k < gtfs.length; k++) {		// exons 
					France f= (France) gtfs[k];
					if (f.isExon()) 
						transcript.updateBoundaries(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
				}
				t2Hash.put(tID, transcript);	// fill tHash with transcripts
			}
			
		}
		
			// cluster
		HashMap gHash= new HashMap();
		Comparator compi= new DirectedRegion.PositionComparator();
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap t2Hash= (HashMap) chr2Hash.get(chrID);
			Object[] transcripts= t2Hash.values().toArray();
			java.util.Arrays.sort(transcripts, compi);
			Transcript[] t= new Transcript[transcripts.length];
			for (int j = 0; j < t.length; j++) 
				t[j]= (Transcript) transcripts[j];
			Transcript[][] loci= clusterTranscripts(t);
			for (int j = 0; j < loci.length; j++) {
				String gID= Gene.getUniqueID();
				Gene locus= new Gene(spec, gID);
				locus.setStrand(loci[j][0].getStrand());
				locus.setChromosome(chrID);
				for (int k = 0; k < loci[j].length; k++) { // transcripts
					loci[j][k].setGene(locus);
					Vector v= (Vector) ((HashMap) chrHash.get(chrID)).get(loci[j][k].getTranscriptID());
					for (int m = 0; m < v.size(); m++) {
						France f= (France) v.elementAt(m);
						if (f.isExon())
							loci[j][k].addExon(new Exon(loci[j][k], f.getExonID(), f.getStart(), f.getEnd()));
						else if (f.isCDS())
							loci[j][k].addCDS(f.getStart(), f.getEnd());
					}
					locus.addTranscript(loci[j][k]);
				}
				gHash.put(gID, locus);
			}
		}		
		
			// build graph
		iter= gHash.values().iterator();
		Graph g= new Graph();
		g.addSpecies(spec);
		while (iter.hasNext()) 
			g.addGene((Gene) iter.next());
		return g;
	}
	
	

	static protected Transcript[][] clusterTranscripts(DirectedRegion[] regions) {
		
		int maxPlus= Integer.MIN_VALUE, maxMin= Integer.MIN_VALUE;
		Vector clusters= new Vector();
		Vector vPlus= null, vMinus= null;
		for (int i = 0; i < regions.length; i++) {

			DirectedRegion r= regions[i];
			if (regions[i].isForward()) {
				if (regions[i].getStart()> maxPlus) {
					if (vPlus!= null)
						clusters.add(vPlus);
					vPlus= new Vector();
				} 
				vPlus.add(regions[i]);
				if (regions[i].getEnd()> maxPlus)
					maxPlus= regions[i].getEnd();
			} else {
				if (Math.abs(regions[i].getStart())> maxMin) {
					if (vMinus!= null)
						clusters.add(vMinus);
					vMinus= new Vector();
				} 
				vMinus.add(regions[i]);
				if (Math.abs(regions[i].getEnd())> maxMin)
					maxMin= Math.abs(regions[i].getEnd());
			}
		}

		if (vPlus!= null)
			clusters.add(vPlus);
		if (vMinus!= null)
			clusters.add(vMinus);
		
		System.gc();
		
		return (Transcript[][]) Arrays.toField(clusters);
	}

	void assemble_old() {
		
//		Species spec= new Species("human");
//		Vector clusters= getClusters(transHash);	// Vector x Vector x Vector
//		
//		Collection c= transHash.values();
//		HashMap tHash= new HashMap(c.size()); 
//		HashMap gHash= new HashMap(c.size()/ 2);
//		Iterator it= clusters.iterator();
//		while (it.hasNext()) {	// iterate all clusters
//			Vector v= (Vector) it.next();
//			Cluster currCluster= new Cluster(((France) ((Vector) v.elementAt(0)).elementAt(0)).getChromosome(), spec);
//			for (int i = 0; i < v.size(); i++) {
//				Vector vTrans= (Vector) v.elementAt(i);
//				for (int j = 0; j < vTrans.size(); j++) {
//					France f= (France) vTrans.elementAt(j);	// exon, cds, ...
//					Transcript t= (Transcript) tHash.get(f.getTranscriptID());	// get corresponding transcript
//					if (t== null) {
//						t= new Transcript(f.getTranscriptID(),f.getStrand(),currCluster);
//						t.setSource(f.getSource());
//						Gene g= (Gene) gHash.get(f.getGeneID());
//						if (g== null) {
//							g= new Gene(f.getGeneID());
//							g.setAlias(f.getGeneAlias());
//							gHash.put(f.getGeneID(), g);
//						}
//						t.addGene(g);
//						tHash.put(f.getTranscriptID(), t);
//					}
////					if (f.isCDS())
////						t.addCDS(f.getExonID(), f.getStart(), f.getEnd());
//				}
//			}
//		}
	}
	
	static protected HashMap getGroups(String id, GTFObject[] obj) {
		
		HashMap hash= new HashMap();
		for (int i = 0; i < obj.length; i++) {
			String s= obj[i].getAttribute(id);
			Vector tAttrib= (Vector) hash.get(s);
			if (tAttrib== null) {
				tAttrib= new Vector();
				hash.put(s, tAttrib);
			}
			tAttrib.add(obj[i]);
		}
		
		return hash;
	}

	private HashMap getTranscriptInfo() {
		
		HashMap tHash= new HashMap();
		France[] france= getFranceObj();
		for (int i = 0; i < france.length; i++) {
			if (!france[i].isExon()&& !france[i].isCDS())	// collect exons and cds
				continue;
			Vector tAttrib= (Vector) tHash.get(france[i].getTranscriptID());
			if (tAttrib== null) {
				tAttrib= new Vector();
				tHash.put(france[i].getTranscriptID(), tAttrib);
			}
			tAttrib.add(france[i]);
		}
		
		return tHash;
	}

		public String[] getGeneIDs(String[] protIDs) {
			Vector v= new Vector();
			Pattern[] patterns= new Pattern[protIDs.length];
			for (int i = 0; i < protIDs.length; i++) 
				patterns[i]= Pattern.compile(protIDs[i]);	// case sensitive .toUpperCase()
			try {
				BufferedReader buffy;
				if (fPath!= null&& fName!= null)
					buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
				else 
					buffy= new BufferedReader(new InputStreamReader(inputStream));
				String line;
				int lineCtr= 0;
				Vector gtfVec= new Vector();
				System.out.println("Scanning for geneIDs ");
				long t0= System.currentTimeMillis();
				while (buffy.ready()) {
					lineCtr++;
					if (lineCtr%1000== 0) {
						System.out.print("*");
						System.out.flush();
					}
					line= buffy.readLine();
					int x= 0;
					for (int i = 0; i < patterns.length; i++) {
						//Matcher m= patterns[i].matcher(line);	// case sensitive .toUpperCase()
						//if (m.find())
						if (line.contains(protIDs[i]))
							break;
					}
					if (x== protIDs.length)
						continue;
					
					String[] tokens= line.split("\t");
					for (int i = 8; i < tokens.length; i++) {
						if (tokens[i].equals(GTFObject.GENE_ID_TAG)) {
							String id= tokens[i+1];
							if (id.startsWith("\"")|| id.startsWith("\'"))
								id= id.substring(1, id.length()- 1);
							v.add(id);
							break;
						}
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			
			System.out.println("\nPrimary list "+v.size());
			Comparator compi= Collator.getInstance();
			Vector vv= new Vector();
			for (int i = 0; i < v.size(); i++) {
				Arrays.addUnique(vv, v.elementAt(i), compi);
			}
			System.out.println("Found "+vv.size()+" gene ids for "+protIDs.length+" proteins.");
			return (String[]) Arrays.toField(vv);
		}
	public void read_what_is_that(String[] geneIDs) throws Exception {
			this.filtGeneIDs= geneIDs;
			read();
	}
		
		public void read_old() throws Exception {
			reset();
			File f= new File(fPath+ File.separator+ fName);
			long fSize= f.length();
			BufferedReader buffy;			
			if (fPath!= null&& fName!= null)
				buffy= new BufferedReader(new FileReader(f));
			else 
				buffy= new BufferedReader(new InputStreamReader(inputStream));
			String line;
			Vector gtfVec= new Vector();
			int lastPerc= 0;
			System.out.print("\treading annotation ");
			System.out.flush();
			while (buffy.ready()) {
				line= buffy.readLine();
				bytesRead+= line.length()+ 1;
				++linesRead;
				int perc= (int) ((bytesRead* 10d)/ fSize);
				if (perc> lastPerc) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				
					// filter gene IDs
				if (filtGeneIDs!= null) {
					String[] cols= line.split("\t");
					if (linesRead== 1) {
						if (!silent&& cols.length< 9&& !cols[8].startsWith(GTFObject.GENE_ID_TAG))
							System.err.println("gene_id not in col 8");
					}
					cols= cols[8].split(" ");
				
					int x= 0;
					for (x = 0; x < filtGeneIDs.length; x++) {
						String chk= cols[1];
						if (cols[1].contains(";"))
							chk= cols[1].substring(0, cols[1].indexOf(';'));
						if (cols[1].contains("\""))
							chk= cols[1].substring(cols[1].indexOf('\"')+ 1, 
									cols[1].lastIndexOf('\"'));
						if (chk.equalsIgnoreCase(filtGeneIDs[x]))
							break;
					}
					if (x== filtGeneIDs.length)
						continue;
				}
				
				StringTokenizer toki= new StringTokenizer(line, " \t");	// must be tab, see specification
				if (!silent&& toki.countTokens()< 8)
					System.err.println("line "+ linesRead+ ": skipped (<8 token)!\n\t"+ line);
				// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
				GTFObject newObj= createGTFObject();		
				boolean addFlag= true;
				try {				
					newObj.seqname= toki.nextToken();
						// neg chrom filtering
					if (noIDs!= null) {
						int u;
						String chrIDu= newObj.seqname.toUpperCase();
						for (u = 0; u < noIDs.length; u++) {
							if (Pattern.matches(noIDs[u], chrIDu)) 
								break;
						}
						if (u< noIDs.length) {
							if (skippedChr== null)
								skippedChr= new Vector();
							skippedChr= Arrays.addUnique(skippedChr, newObj.seqname);
							continue;		// next line
						}
					}
					
						// pos chromosome filtering
					if (filtChrIDs!= null) {
						int u;
						for (u = 0; u < filtChrIDs.length; u++) {
							if (filtChrIDs[u].equalsIgnoreCase(newObj.seqname))
								break;
						}
						if (u== filtChrIDs.length) {
							if (skippedChr== null)
								skippedChr= new Vector();
							skippedChr= Arrays.addUnique(skippedChr, newObj.seqname);
							continue;		// next
						}
					}

					newObj.source= toki.nextToken();
					newObj.setFeature(toki.nextToken());
					newObj.start= Integer.parseInt(toki.nextToken());
					newObj.end= Integer.parseInt(toki.nextToken());
					newObj.setScore(toki.nextToken());
					newObj.setStrand(toki.nextToken());
					newObj.setFrame(toki.nextToken());
				} catch (Exception e) {
					if (!silent)
						System.err.println("Invalid GTF format: line "+ linesRead+" -- "+e.getMessage());
					//e.printStackTrace();
					//continue;
				}
				
					// optional attributes
				int smc= line.indexOf(';');		// GTF2
				if (smc>= 0) {
					String ss= toki.nextToken();
	//				toki= new StringTokenizer(line, " \t");	// must be tab, see specification
	//				for (int i = 0; i < 8; i++) 
	//					ss= toki.nextToken();
	//				String h= line.substring(0, smc);			// last ';'
	//				h= line.substring(0, h.lastIndexOf(' '));	// two ' ' tokens before
	//				h= line.substring(0, h.lastIndexOf(' '));
					String h= line.substring(line.indexOf(ss), line.length()).trim();	// skip that part
					
					toki= new StringTokenizer(h, ";");		// attributes
					while (toki.hasMoreTokens()) {
						h= toki.nextToken().trim();
						int sep= h.indexOf(' ');
						if (sep < 0) {						// comments
							String s= h;
							while (toki.hasMoreTokens())
								s+= " "+ toki.nextToken();
							newObj.setComments(s);
						}
						if (sep>= 0) {
							String id= h.substring(0, sep);
							String val= h.substring(sep+ 1, h.length());
							
								// filter gene IDs
//							if (id.equals(GTFObject.GENE_ID_TAG)&& filtGeneIDs!= null) {
//								int i;
//								for (i = 0; i < filtGeneIDs.length; i++) {
//									if (filtGeneIDs[i].equals(val))
//										break;
//								}
//								if (i== filtGeneIDs.length) {
//									addFlag= false;
//									break;
//								}
//									
//							}
//								// filter transcript IDs
//							if (id.equals(GTFObject.TRANSCRIPT_ID_TAG)&& filtTrptIDs!= null) {
//								int i;
//								for (i = 0; i < filtTrptIDs.length; i++) {
//									if (filtTrptIDs[i].equals(val))
//										break;
//								}
//								if (i== filtTrptIDs.length) {
//									addFlag= false;
//									break;
//								}
//									
//							}
							newObj.addAttribute(id, val);
						}
					}
				}
				
				if (addFlag)
					gtfVec.add(newObj);
				//System.out.println(gtfVec.size());
			}
			System.out.println();
			
			if (skippedChr!= null) {
				String s= "";
				for (int i = 0; i < skippedChr.size(); i++) 
					s+= skippedChr.elementAt(i)+ ", ";
				System.out.println("\tSkipped "+ skippedChr.size()+ " chroms: "+s.substring(0,s.length()- 2));
			}
			gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
		}

		public Gene readNextGene_old() {
						
				RandomAccessFile raf= null;
				long size= 0l;
				if (fPath!= null&& fName!= null) {
					try {
						raf= new RandomAccessFile(new File(fPath+ File.separator+ fName),"r");
						size= raf.length();
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if(inputStream!= null)
					System.err.println(this.getClass()+" only supports input from Files.");
				
				Gene gene= null;
				try {
					raf.seek(bytesRead);
					
					String line;
					int lineCtr= 0;
					
					String lastTID= null;			
					String lastGID= null;			
					Transcript trpt= null;
					int lastStrand= 0;
					while (true) {
						
							// read line
						line= raf.readLine();
						if (line== null)
							return gene;
							
						bytesRead+= line.length()+ 1;
						int perc= (int) ((bytesRead* 10d)/ size);
						if (perc> lastPerc) {
							++lastPerc;
							System.out.print("*");
							System.out.flush();
						}
						++linesRead;
						
						String[] tokens= line.split("\t");	// must be tab, see specification
						if (tokens.length< 8)
							System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
						int start= Integer.parseInt(tokens[3]);
						int end= Integer.parseInt(tokens[4]);
						int strand= 0;
						if (tokens[6].trim().equals("+"))
							strand= 1;
						else if (tokens[6].trim().equals("-"))
							strand= -1;
						else {
							System.err.println("No strand assignment, line "+linesRead);
							continue;
						}						
						if (lastStrand== 0)
							lastStrand= strand;
						if (strand!= lastStrand)
							break;	// end of gene
						
						// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
						String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
						String[] attrTokens= h.split(";");
						
							// get TID
						String tid= null;
						String gid= null;
						for (int i = 0; i < attrTokens.length; i++) {
							h= attrTokens[i].trim();
							int sep= h.indexOf(' ');
							if (sep < 0) 						// comments
								break;
							if (sep>= 0) {
								String id= h.substring(0, sep);
								if (id.equals(GTFObject.TRANSCRIPT_ID_TAG))
									tid= h.substring(sep+ 1, h.length());
								if (id.equals(GTFObject.GENE_ID_TAG))
									gid= h.substring(sep+ 1, h.length());
							}
						}
						if (tid== null) {
							System.err.println("No TID found.");
							continue;
						}
			
						// new gene?
						if (!gid.equals(lastGID)) {
							if (lastGID!= null)
								break;	// exit
							gene= new Gene(gid);
							gene.setStrand(strand);
							gene.setChromosome(tokens[0]);
							lastGID= gid;
						}
			
							// new transcript?
						if (!tid.equals(lastTID)) {
							trpt= new Transcript(tid);
							trpt.setStrand(strand);
							trpt.setChromosome(tokens[0]);
							trpt.setSource(tokens[1]);
							gene.addTranscript(trpt);
							trpt.setGene(gene);
							lastTID= tid;
							lastStrand= strand;
						}
												
						if (tokens[2].equalsIgnoreCase("exon")) {
							String exonID= null;
							for (int i = 0; i < attrTokens.length; i++) {
								h= attrTokens[i];
								int sep= h.indexOf(' ');
								if (sep < 0) 						// comments
									break;
								if (sep>= 0) {
									String id= h.substring(0, sep);
									if (id.equals(GTFObject.EXON_ID_TAG))
										exonID= h.substring(sep+ 1, h.length());
									break;
								}
							}
							Exon e= new Exon(trpt, exonID, start, end);
							trpt.addExon(e);				
						} else if (tokens[2].equalsIgnoreCase("CDS")) 
							trpt.addCDS(start, end);
			
					}
					raf.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
				
				return gene;
			}

	public void reformatFile_diff() {
			try {
				System.out.println("reformatting..");
				HashMap chrMap= new HashMap();
				File f= new File(fPath+ File.separator+ fName);
				int p= fName.lastIndexOf('.');
				if (p< 0)
					p= fName.length();
				//File fold= new File(fPath+ File.separator+ fName.substring(0, p)+"_original");
				File fnew= new File(fPath+ File.separator+ fName.substring(0, p)+"_"+NORMED_FILE_TAG+"_");
				long size= f.length(); 
				BufferedWriter writer= null;
				
				BufferedReader reader= new BufferedReader(new FileReader(f));
				String line= null;
				long offset= 0l;
				String lastWriterID= null;
				System.out.print("\treformatting ");
				System.out.flush();
				HashMap readChrMap= new HashMap();
				while (reader.ready()) {
					
						// read line
					line= reader.readLine();
					if (line== null)
						break;
						
					offset+= line.length()+ 1;
					int perc= (int) ((offset* 10d)/ size);
					if (perc> lastPerc) {
						++lastPerc;
						System.out.print("*");
						System.out.flush();
					}
					
					String[] tokens= line.split("\t");	// must be tab, see specification
						
					String chrID= tokens[0];		// check chromosome
					if (readChrMap.get(chrID)!= null)
						continue;	// already handled
					
						// get all obj for chr
					BufferedReader buffy= new BufferedReader(new FileReader(f));
					Vector gtfObjV= new Vector();
					while (buffy.ready()) {
						line= buffy.readLine();
						tokens= line.split("\t");
						if (tokens[0].equals(chrID))
							continue;
						
						GTFObject o= new GTFObject(line);
						gtfObjV.add(o);
					}
					
					String writerID= chrID+"_"+tID;
					if (!writerID.equals(lastWriterID)) {	// change writer
						if (writer!= null) {
							writer.flush();
							writer.close();
						}
						File fx= new File(f.getAbsolutePath()+"_"+writerID);
						if (chrMap.get(writerID)== null) {
							if (fx.exists())
								fx.delete();
							chrMap.put(writerID, writerID);
						} 
						writer= new BufferedWriter(new FileWriter(fx, true));
						fx.deleteOnExit();	// doesnt work
						lastWriterID= writerID;
					}
					
					writer.write(line+"\n");
					//writer.newLine();
					writer.flush();
				}
				System.out.println();
				reader.close();
				//f.renameTo(fold);
				writer= new BufferedWriter(new FileWriter(fnew.getAbsolutePath()));
				lastPerc= 0;
				
					// concatenate
				System.out.print("\tfuse ");
				offset= 0l;
				Object[] keys= chrMap.keySet().toArray();
				java.util.Arrays.sort(keys);
				Comparator gtfComp= new GTFObject.PositionComparator();
				for (int i = 0; i < keys.length; i++) {
					File fx= new File(f.getAbsolutePath()+"_"+keys[i]);
					reader= new BufferedReader(new FileReader(fx));
					HashMap trptHash= new HashMap(), startHash= new HashMap(), startHashRev= new HashMap();
					while (reader.ready()) {
						line= reader.readLine();
						GTFObject o= new GTFObject(line);
						if (o.getTranscriptID()== null|| o.getChromosome()== null|| o.getStart()<= 0)
							System.out.println("Invalid GTF line "+line);
						Vector v= (Vector) trptHash.get(o.getTranscriptID());
						if (v== null)
							v= new Vector();
						Object[] a= v.toArray();
						int pp= java.util.Arrays.binarySearch(a, o, gtfComp);
						if (pp< 0)
							pp= (pp+1)* (-1);					
						v.insertElementAt(o, pp);
						trptHash.put(o.getTranscriptID(), v);
						
						Integer start= (Integer) startHash.get(o.getTranscriptID());
						int newStart= o.getStart();
						if ((start== null)|| (newStart< start.intValue()))  
							start= new Integer(newStart);
						
						v= (Vector) startHashRev.get(o.getTranscriptID());
						if (v== null)
							v= new Vector();
						v= Arrays.addUnique(v, o.getTranscriptID());
						startHashRev.put(start, v);
		
						
						offset+= line.length()+ 1;
						int perc= (int) ((offset* 10d)/ size);
						if (perc> lastPerc) {
							++lastPerc;
							System.out.print("*");
							System.out.flush();
						}
					}
					
					keys= startHashRev.keySet().toArray();
					java.util.Arrays.sort(keys);
					for (int j = 0; j < keys.length; j++) {
						Vector v= (Vector) startHashRev.get(keys[i]);
						for (int k = 0; k < v.size(); k++) {
							Vector vv= (Vector) trptHash.get(v.elementAt(k));
							for (int m = 0; m < vv.size(); m++) {
								writer.write(vv.elementAt(m)+"\n");
								//writer.newLine();
							}
						}
					}
					
					
					reader.close();
					fx.delete();
				}
				System.out.println();
				writer.flush(); writer.close();
				
				//System.out.println("\tsuccess, old file in "+fold.getAbsolutePath());
				System.out.println("\tsuccess, new file in "+fnew.getAbsolutePath());
				reset();
			} catch (Exception e) {
				System.err.println("Reformatting failed.");
				e.printStackTrace();
			}
		}

	public void reformatFile() {
		try {
			System.out.println("reformatting..");
			HashMap chrMap= new HashMap();
			File f= new File(fPath+ File.separator+ fName);
			int p= fName.lastIndexOf('.');
			if (p< 0)
				p= fName.length();
			//File fold= new File(fPath+ File.separator+ fName.substring(0, p)+"_original");
			File fnew= new File(fPath+ File.separator+ fName.substring(0, p)+"_"+NORMED_FILE_TAG+"_");
			long size= f.length(); 
			BufferedWriter writer= null;
			
			BufferedReader reader= new BufferedReader(new FileReader(f));
			String line= null;
			long offset= 0l;
			String lastWriterID= null;
			System.out.print("\tsplit ");
			System.out.flush();
			while (reader.ready()) {
				
					// read line
				line= reader.readLine();
				if (line== null)
					break;
					
				offset+= line.length()+ 1;
				int perc= (int) ((offset* 10d)/ size);
				if (perc> lastPerc) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				
				String[] tokens= line.split("\t");	// must be tab, see specification
					
				String chrID= tokens[0];		// check chromosome
				String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
				String[] attrTokens= h.split(";");
				String tID= null;
				for (int i = 0; i < attrTokens.length; i++) {
					h= attrTokens[i].trim();
					int sep= Math.max(h.indexOf(' '), h.indexOf('='));	// = in ASD gff
					if (sep>= 0) {
						String id= h.substring(0, sep);
						if (id.equals(GTFObject.TRANSCRIPT_ID_TAG)) {
							tID= h.substring(sep+ 1, h.length());
							if (tID.charAt(0)== '\"')
								tID= tID.substring(1, tID.length()- 1);
						}
					}
				}
				
				String writerID= chrID+"_"+tID;
				if (!writerID.equals(lastWriterID)) {	// change writer
					if (writer!= null) {
						writer.flush();
						writer.close();
					}
					File fx= new File(f.getAbsolutePath()+"_"+writerID);
					if (chrMap.get(writerID)== null) {
						if (fx.exists())
							fx.delete();
						chrMap.put(writerID, writerID);
					} 
					writer= new BufferedWriter(new FileWriter(fx, true));
					fx.deleteOnExit();	// doesnt work
					lastWriterID= writerID;
				}
				
				writer.write(line+"\n");
				//writer.newLine();
				writer.flush();
			}
			System.out.println();
			reader.close();
			//f.renameTo(fold);
			writer= new BufferedWriter(new FileWriter(fnew.getAbsolutePath()));
			lastPerc= 0;
			
				// concatenate
			System.out.print("\tfuse ");
			offset= 0l;
			Object[] keys= chrMap.keySet().toArray();
			java.util.Arrays.sort(keys);
			Comparator gtfComp= new GTFObject.PositionComparator();
			for (int i = 0; i < keys.length; i++) {
				File fx= new File(f.getAbsolutePath()+"_"+keys[i]);
				reader= new BufferedReader(new FileReader(fx));
				HashMap trptHash= new HashMap(), startHash= new HashMap(), startHashRev= new HashMap();
				while (reader.ready()) {
					line= reader.readLine();
					GTFObject o= new GTFObject(line);
					if (o.getTranscriptID()== null|| o.getChromosome()== null|| o.getStart()<= 0)
						System.out.println("Invalid GTF line "+line);
					Vector v= (Vector) trptHash.get(o.getTranscriptID());
					if (v== null)
						v= new Vector();
					Object[] a= v.toArray();
					int pp= java.util.Arrays.binarySearch(a, o, gtfComp);
					if (pp< 0)
						pp= (pp+1)* (-1);					
					v.insertElementAt(o, pp);
					trptHash.put(o.getTranscriptID(), v);
					
					Integer start= (Integer) startHash.get(o.getTranscriptID());
					int newStart= o.getStart();
					if ((start== null)|| (newStart< start.intValue()))  
						start= new Integer(newStart);
					
					v= (Vector) startHashRev.get(o.getTranscriptID());
					if (v== null)
						v= new Vector();
					v= Arrays.addUnique(v, o.getTranscriptID());
					startHashRev.put(start, v);

					
					offset+= line.length()+ 1;
					int perc= (int) ((offset* 10d)/ size);
					if (perc> lastPerc) {
						++lastPerc;
						System.out.print("*");
						System.out.flush();
					}
				}
				
				keys= startHashRev.keySet().toArray();
				java.util.Arrays.sort(keys);
				for (int j = 0; j < keys.length; j++) {
					Vector v= (Vector) startHashRev.get(keys[i]);
					for (int k = 0; k < v.size(); k++) {
						Vector vv= (Vector) trptHash.get(v.elementAt(k));
						for (int m = 0; m < vv.size(); m++) {
							writer.write(vv.elementAt(m)+"\n");
							//writer.newLine();
						}
					}
				}
				
				
				reader.close();
				fx.delete();
			}
			System.out.println();
			writer.flush(); writer.close();
			
			//System.out.println("\tsuccess, old file in "+fold.getAbsolutePath());
			System.out.println("\tsuccess, new file in "+fnew.getAbsolutePath());
			reset();
		} catch (Exception e) {
			System.err.println("Reformatting failed.");
			e.printStackTrace();
		}
	}
	
	public void reset() {
		bytesRead= 0l;
		linesRead= 0;
		readChr= new Vector();
		skippedChr= new Vector();
		skippedTranscripts= new Vector();
		skippedTranscripts= new Vector();
		// TODO: region and ID-filtering?? or do it method-level.. 
	}

	
	/*
	 * RAF cannot insert :(
	 */
	public void reformatFile_nice_but_not_work() {
		try {
			DualHashBidiMap chrOffsetMap= new DualHashBidiMap();
			
			File f= new File(fPath+ File.separator+ fName);
			long size= f.length(); 
			File fout= new File(fPath+ File.separator+ fName+"_reformatted");
			RandomAccessFile writer= new RandomAccessFile(fout, "rw");
			
			BufferedReader reader= new BufferedReader(new FileReader(f));
			String line= null;
			long offset= 0l;
			String lastChrID= null;
			Vector chrRead= new Vector();
			while (reader.ready()) {
				
					// read line
				line= reader.readLine();
				if (line== null)
					break;
					
				offset+= line.length()+ 1;
				int perc= (int) ((offset* 10d)/ size);
				if (perc> lastPerc) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				
				String[] tokens= line.split("\t");	// must be tab, see specification
					
				String chrID= tokens[0];		// check chromosome
				if (!chrID.equals(lastChrID)) {
					if (chrOffsetMap.get(chrID)== null) 
						chrOffsetMap.put(chrID, new Long(offset));
					else 
						writer.seek(((Long) chrOffsetMap.get(chrID)).longValue());
					lastChrID= chrID;
				}
				
					// update fPointers
				long fPointer= writer.getFilePointer();
				int lineLen= line.length()+ 1;
				Object[] oOff= chrOffsetMap.values().toArray();
				BidiMap invMap= chrOffsetMap.inverseBidiMap();
				for (int i = 0; i < oOff.length; i++) { 
					Long offL= (Long) oOff[i];
					if (offL.longValue()> fPointer) {
						Object key= invMap.get(offL);
						chrOffsetMap.remove(key);
						chrOffsetMap.put(key, new Long(offL.longValue()+ lineLen));
					}
				}
				writer.writeChars(line+"\n");
			}
			
			writer.close();
			reader.close();
		} catch (Exception e) {
			System.err.println("Reformatting failed.");
			e.printStackTrace();
		}
	}
	
	
	public void sweepToChromosome(String chrom) {
		BufferedReader buffy= null;
		long size= 0l;
		
		if (fPath!= null&& fName!= null) {
			try {
				File file= new File(fPath+ File.separator+ fName);
				size= file.length();
				inputStream= new FileInputStream(file);
				inputStream.skip(bytesRead);
				buffy= new BufferedReader(new InputStreamReader(inputStream));
				//buffy= new BufferedReader(new FileReader(file));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
			
		
//		reset();	// start at the beginning
//		int lastPerc= (int) ((bytesRead* 10d)/ size);
		try {
			//raf.seek(bytesRead);
			boolean reset= false;
			while (true) {
				String line= buffy.readLine();			
				if (line== null) {
					if (reset) {
						System.out.println("WARNING: chromosome "+chrom+" not found for sweeping");
						break;
					} else {
						File file= new File(fPath+ File.separator+ fName);
						buffy= new BufferedReader(new FileReader(file));
						reset();
						line= buffy.readLine();
						reset= true;	// only once
					}

				}
				
				String[] tokens= line.split("\t");
				String chr= null;
				if (tokens[0].length()<= 3)
					chr= "chr"+ tokens[0].trim();
				else
					chr= tokens[0].trim();
				if (chr.equalsIgnoreCase(chrom))
					break;
				
				
				bytesRead+= (line.length()+ 1);		// else
				++linesRead;
				int perc= (int) ((bytesRead* 10d)/ size);
//				if (perc> lastPerc&& !silent) {
//					++lastPerc;
//					System.out.print("*");
//					System.out.flush();
//				}
			}
			buffy.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return;
	
	
	}

	void skipToNextChromosome_raf(RandomAccessFile raf) {
			try {
				//raf.seek(bytesRead);
				String lastChrID= null;
				int lastPerc= 0;
				long t0= System.currentTimeMillis();
				while (true) {
					int perc= (int) ((bytesRead* 10d)/ raf.length());
					if (perc> lastPerc&& !silent) {
						++lastPerc;
						System.out.print("*");
						System.out.flush();
					}
					String line= raf.readLine();			
					if (line== null) {
						break;
					}
	//				String[] tokens= line.split("\t");
	//				String chr= tokens[0];
					StringTokenizer toki= new StringTokenizer(line, "\t");
					String chr= toki.nextToken();
					if (chr.length()<= 3)
						chr= "chr"+ chr.trim();
					else
						chr= chr.trim();
					
					
					
					if (lastChrID== null) {
						lastChrID= chr;
						System.out.print("skip "+lastChrID+" -> ");
						System.out.flush();
					}
					if (!chr.equals(lastChrID)) {
						System.out.println(chr+ " "+
								((System.currentTimeMillis()- t0)/ 1000)+ " sec");
						break;
					}
					bytesRead+= (line.length()+ 1);		// else
					++linesRead;
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
			
		}

	void skipToNextChromosome(BufferedReader buffy, long size) {
		try {
			//raf.seek(bytesRead);
			String lastChrID= null;
			int lastPerc= (int) ((bytesRead* 10d)/ size);
//			long t0= System.currentTimeMillis();
			while (true) {
				int perc= (int) ((bytesRead* 10d)/ size);
				if (perc> lastPerc&& !silent) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				
				buffy.mark(MAX_GTF_LINE_LENGTH);
				String line= buffy.readLine();			
				if (line== null) {
					break;
				}
				String[] tokens= line.split("\t");
				String chr= tokens[0];
//				StringTokenizer toki= new StringTokenizer(line, "\t");
//				String chr= toki.nextToken();
				if (chr.length()<= 3)
					chr= "chr"+ chr.trim();
				else
					chr= chr.trim();
				
				
				
				if (lastChrID== null) {
					lastChrID= chr;
//					System.out.print("skip "+lastChrID+" -> ");
//					System.out.flush();
				}
				if (!chr.equals(lastChrID)) {
//					System.out.println(chr+ " "+
//							((System.currentTimeMillis()- t0)/ 1000)+ " sec");
					buffy.reset();
					break;
				}
				bytesRead+= (line.length()+ 1);		// else
				++linesRead;
			}
 			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	
	public String getNextChromosome() {
		RandomAccessFile raf= null;
		long size= 0l;
		if (fPath!= null&& fName!= null) {
			try {
				raf= new RandomAccessFile(new File(fPath+ File.separator+ fName),"r");
				size= raf.length();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else if(inputStream!= null)
			System.err.println(this.getClass()+" only supports input from Files.");
		
		Gene[] geneReg= null;
		Comparator startCompi= new AbstractRegion.StartComparator();
		try {
			raf.seek(bytesRead);
			String line= raf.readLine();
			if (line== null)
				return null;
			raf.close();
			
			String[] tokens= line.split("\t");
			if (tokens[0].length()<= 3)
				return "chr"+ tokens[0];
			else
				return tokens[0];			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;

	}
	
	public boolean checkChr(Gene[] ge) {
		int posStrand= 0;
		int negStrand= 0;
		boolean ok= true;
		for (int i = 0; i < ge.length; i++) {
			if (ge[i].getStrand()> 0)
				++posStrand;
			else if (ge[i].getStrand()< 0)
				++negStrand;
			if (ge[i].getStart()== 0|| ge[i].getEnd()== 0||
					ge[i].getChromosome()== null|| ge[i].getStrand()== 0) {
				System.err.println("Gene not inited "+ge[i].getGeneID());
				ok= false;
			}
			Transcript[] trpt= ge[i].getTranscripts();
			if (trpt.length== 0|| trpt.length> 10) {
				System.err.println("Too many transcripts per gene "+ge[i].getGeneID());
				ok= false;
			}
			Exon[] ex= ge[i].getExons();
			if (ex.length== 0|| ex.length> 100) {
				System.err.println("Too many exons per gene "+ge[i].getGeneID());
				ok= false;
			}
				
			for (int j = 0; j < trpt.length; j++) {
				if (trpt[j].getStart()== 0|| trpt[j].getEnd()== 0||
						trpt[j].getChromosome()== null|| trpt[j].getStrand()== 0) {
					System.err.println("Transcript not inited "+trpt[j].getTranscriptID());
					ok= false;
				}

			}
			for (int j = 0; j < ex.length; j++) {
				if (ex[j].getStart()== 0|| ex[j].getEnd()== 0||
						ex[j].getChromosome()== null|| ex[j].getStrand()== 0) {
					System.err.println("Exon not inited "+ex[j].getTranscripts()[0].getTranscriptID());
					ok= false;
				}

			}
		}
		
		if (Math.min(((double) posStrand)/ negStrand, ((double) negStrand)/ posStrand)> 0.75d) {
			System.err.println("Unequal strand distribution, pos "+posStrand+", neg "+negStrand);
			ok= false;
		}
		return ok;
	}
	/**
				 * Requires the gtf file to be ordered by chromosomes AND that
				 * data (CDS, exon) from the same transcript is in adjacent lines
				 * (which is good for RefSeq in worm, where the same tID/gID exists
				 * on 4 diferent chromosomes). 
				 * It does NOT require any order of the transcript starts within 
				 * the chromosome, clustering is performed exhaustively.
				 * 
				 * retrieves every transcript that is overlapping a certain region of
				 * <code>ovlReg</code>
				 * @return
				 */
				/*
				 * Note that the clustering cannot be separated from the dynamic
				 * gene building, since reduncancy checks for splice sites and exons
				 * within the genes depend on it.
				 */
				public void read() throws Exception { 
								
						BufferedReader buffy= null;
						long size= 0l;
						clustered= false;
						
						if (fPath!= null&& fName!= null) {
							try {
								File file= new File(fPath+ File.separator+ fName);
								size= file.length();
								inputStream= new FileInputStream(file);
								inputStream.skip(bytesRead);
								buffy= new BufferedReader(new InputStreamReader(inputStream));
							} catch (Exception e) {
								e.printStackTrace();
							}
						}
							
						
						Vector geneV= new Vector();
						Comparator startCompi= new DirectedRegion.StartComparator();
							
						Vector gtfV= null;
						if (readGTF)
							gtfV = new Vector();
						try {
							
							String line= null;
							
							String lastTID= null;			
							String lastChrID= null;
							Transcript trpt= null;
							//boolean merged= true;	// for skiping insertion of init null trpt/gene
							while (true) {
								
									// read line
								buffy.mark(MAX_GTF_LINE_LENGTH);	// for skipping chromosomes
								line= buffy.readLine();
								if (line== null) {
									if (readGene&& trpt!= null) 
										geneV.add(trpt.getGene());
									
									
									
										// print statistics
									if (!silent) {
										System.out.println("\nRead:\t"+readChrs+" chromosomes, "+readGenes+" genes, "+readTranscripts+" transcripts, "+
												readExons+" exons.");
										System.out.println("Skip:\t"+skippedChr.size()+" chromosomes, "+skippedTranscripts.size()+" transcripts, "+
												skippedFeatures.size()+" features, "+skippedObjects+" objects.");
										System.out.print("Chromosomes:\t");
										for (int i = 0; i < skippedChr.size(); i++) 
											System.out.print(skippedChr.elementAt(i)+" ");
										System.out.print("\nTranscripts:\t");
										for (int i = 0; i < skippedTranscripts.size(); i++) 
											System.out.print(skippedTranscripts.elementAt(i)+" ");
										System.out.print("\nFeatures:\t");
										for (int i = 0; i < skippedFeatures.size(); i++) 
											System.out.print(skippedFeatures.elementAt(i)+" ");
										System.out.println();
									}
									break;
								}
								++linesRead;
								bytesRead+= line.length()+ 1;
									
								int perc= (int) ((bytesRead* 10d)/ size);
								if (perc> lastPerc&& !silent) {
									++lastPerc;
									System.out.print("*");
									System.out.flush();
								}
															
								String[] tokens= line.split("\t");	// must be tab, see specification
								if (!silent&& tokens.length< 8)
									System.err.println("line "+ linesRead+ ": skipped (<8 token)!\n\t"+ line);
								
									
								int x;		// skip mRNA, gene... 
								for (x = 0; readFeatures!= null&& x < readFeatures.length; x++) {
									if (readFeatures[x].equalsIgnoreCase(tokens[2]))
										break;
								}
								if (readFeatures!= null&& x== readFeatures.length) {
									skippedFeatures= Arrays.addUnique(skippedFeatures, tokens[2]);
									continue;
								}
								
								int start= Integer.parseInt(tokens[3]);
								int end= Integer.parseInt(tokens[4]);
								int strand= 0;
								strand= GTFObject.parseStrand(tokens[6].trim());
								if (strand== 0) {
									if (!silent)
										System.err.println("No strand assignment, line "+linesRead);
									continue;
								}						
								
								// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
								GTFObject obj= new GTFObject();	// build gtf-obj, handling of attributes easier
								obj.setStart(start);
								obj.setEnd(end);
								obj.setStrand(strand);
								obj.setFeature(tokens[2]);
								obj.setSource(tokens[1]);
								obj.setSeqname(tokens[0]);
								obj.setScore(tokens[5]);
								obj.setFrame(tokens[7]);
								String[] attrTokens= new String[0];
								if (tokens.length> 8) {
									String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
									attrTokens= h.split(";");
									for (int i = 0; i < attrTokens.length; i++) {
										h= attrTokens[i].trim();
										h= h.replaceAll("\\s+", " ");
										int sep= Math.max(h.indexOf(' '), h.indexOf('='));	// = in ASD gff
										if (sep < 0) 						// comments
											break;
										if (sep>= 0) {
											String id= h.substring(0, sep);
											String val= h.substring(sep+ 1, h.length());
											if (val.charAt(0)== '\"')
												val= val.substring(1, val.length()- 1);
											obj.addAttribute(id, val);
										}
									}
								}
	
									// filter gtfs that have ID...
								if (filtSomeIDs!= null) {
									boolean found= false;
									for (int i = 0; i < attrTokens.length; i++) {
										String h= attrTokens[i].trim(); 
										int sep= h.indexOf(' ');
										if (sep < 0) 						// comments
											break;
										if (sep>= 0) {
											String val= h.substring(sep+1, h.length()).trim();
											if (val.startsWith("\""))
												val= val.substring(1, val.length()-1);
											String[] t= val.split(" ");
											for (int j = 0; j < t.length; j++) {
												Integer in= (Integer) filtSomeIDs.get(t[j]);
												if (in== null)
													continue;	// ID not in filter list
												filtSomeIDs.put(t[j], new Integer(in.intValue()+ 1));
												found= true;
											}
										}
									}
									if (!found) {
										skippedObjects++;
										continue;
									}
								}
								
									// check chromosome
								String chrID= tokens[0];
								if (fName.contains("cow")&& chrID.equals("30"))
									chrID= "X";
								if (chromosomeWise&& !chrID.equals(lastChrID)) {
									if (lastChrID== null) {
										// neg chromosome filtering
										if (noIDs!= null) {
											int u;
											String chrIDu= chrID.toUpperCase();
											for (u = 0; u < noIDs.length; u++) {
												if (Pattern.matches(noIDs[u], chrIDu)) 
													break;
											}
											if (u< noIDs.length) {
												--linesRead;
												bytesRead-= line.length()+ 1;
												buffy.reset();
												if (skippedChr== null)
													skippedChr= new Vector();
												skippedChr= Arrays.addUnique(skippedChr, chrID);
												skipToNextChromosome(buffy, size);
												continue;
											}
										}
											// pos chromosome filtering
										if (filtChrIDs!= null) {
											int u;
											for (u = 0; u < filtChrIDs.length; u++) {
												if (filtChrIDs[u].equalsIgnoreCase(chrID))
													break;
											}
											if (u== filtChrIDs.length) {
												--linesRead;
												bytesRead-= line.length()+ 1;
												buffy.reset();
												if (skippedChr== null)
													skippedChr= new Vector();
												skippedChr= Arrays.addUnique(skippedChr, chrID);
												skipToNextChromosome(buffy, size);
												continue;
											}
										}
										
										lastChrID= chrID;
										for (int i = 0; readChr!= null&& i < readChr.size(); i++) 
											if (readChr.elementAt(i).equals(chrID)) {
												if (!silent) {
													System.err.println("Chromosome "+chrID+" already read! ("+linesRead+")");
													System.err.println(line);
												}
												return;
											}
										
										
									} else {
										if (readGene)
											geneV.add(trpt.getGene());
										if (readChr== null)
											readChr= new Vector();
										readChr.add(lastChrID);
										--linesRead;
										bytesRead-= line.length()+ 1;										
										break;
									}
									
									readChrs++;
								}
	
								/*	--- END OF FILTERING --- */
								
								if (readGTF) 
									gtfV.add(obj);
								if (!readGene)
									continue;
								
									// get GID, TID
								String tid= obj.getAttribute(transcript_id_tag);
								if (tid== null) {
									if (!silent)
										System.err.println("No TID found.");
									continue;
								} 
					
									// new transcript?
								if (!tid.equals(lastTID)) {
									boolean overlap= true;
									if (filtRegs!= null) {	// should be ok for the moment. maybe also check in merging step
										int i;
										for (i = 0; i < filtRegs.length; i++) 
											if (filtRegs[i].overlaps(trpt))
												break;
										if (i== filtRegs.length)
											overlap= false;
									}
										// pull identical exon attributes up to transcript level
									if (trpt!= null&& trpt.getExons()[0].getAttributes()!= null) {
										Object[] keys= trpt.getExons()[0].getAttributes().keySet().toArray();
										for (int i = 0; i < keys.length; i++) {
											int j;
											for (j = 1; j < trpt.getExons().length; j++) {
												if (!trpt.getExons()[j].getAttribute(keys[i]).equals(trpt.getExons()[0].getAttribute(keys[i])))
													break;
											}
											if (j== trpt.getExons().length) {
												Object val= trpt.getExons()[0].getAttribute(keys[i]);
												for (int k = 0; k < trpt.getExons().length; k++) 
													trpt.getExons()[k].getAttributes().remove(keys[i]);
												trpt.addAttribute(keys[i], val);
											}
										}
									}
									
		
									// save gene
									if (overlap) {	
										//if (!merged) 	// if last transcript's gene couldnt be merged with anything
										if (trpt!= null)
											geneV.add(trpt.getGene());
									} else {
										skippedTranscripts.add(lastTID);
									}
									//merged= false;
									Gene ge= new Gene(Gene.getUniqueID());
									ge.setStrand(strand);
									ge.setChromosome(tokens[0]);
									ge.setSpecies(species);
									trpt= new Transcript(tid);
									lastTID= tid;
									trpt.setStrand(strand);
									trpt.setChromosome(tokens[0]);
									trpt.setSource(tokens[1]);
									ge.addTranscript(trpt);
									readTranscripts++;
								}
														
								if (obj.getFeature().equalsIgnoreCase("exon")) {
									String exonID= obj.getAttribute(GTFObject.EXON_ID_TAG);
									Exon e= new Exon(trpt, exonID, start, end);
									e.setChromosome(chrID);
									e.setStrand(strand);
									++readExons;
									
										// check if gene merge, not necessary when transcript starts are sorted ascending
									//trpt.addExon(e);		// automatically now in Exon constructor	
	//								if (geneReg!= null) {
	//									int p= java.util.Arrays.binarySearch(geneReg, e, startCompi);	// 
	//									if (p>= 0) {
	//										if (geneReg[p]!= trpt.getGene()) {
	//											if (merged) {
	//												trpt.getGene().merge(geneReg[p]);
	//												geneReg= (Gene[]) Arrays.remove(geneReg, p);
	//											} else {
	//												geneReg[p].merge(trpt.getGene());
	//												merged= true;
	//											}
	//										}
	//									} else {	// between two
	//										p= -(p+1);	// insertion point
	//										int q= p-1;
	//										if (q>= 0&& geneReg[q]!= trpt.getGene()&&
	//												e.overlaps(geneReg[q])) { 
	//											if (merged) {	// if already merged, one has to be removed from array
	//												trpt.getGene().merge(geneReg[q]);
	//												geneReg= (Gene[]) Arrays.remove(geneReg, q);
	//											} else
	//												geneReg[q].merge(trpt.getGene());
	//											merged= true;
	//										}
	//										q= p;	//+1; NO, insertion point is the upper neighbor
	//										if (q< geneReg.length&& geneReg[q]!= trpt.getGene()&&
	//												e.overlaps(geneReg[q])) { 
	//											if (merged) {
	//												trpt.getGene().merge(geneReg[q]);
	//												geneReg= (Gene[]) Arrays.remove(geneReg, q);
	//											} else
	//												geneReg[q].merge(trpt.getGene());
	//											merged= true;
	//										}
	//									}
	//								}
									
									
									
								} else if (tokens[2].equalsIgnoreCase("CDS")) {
									trpt.addCDS(start, end);
									Object[] keys= obj.getAttributes().keySet().toArray();
									for (int i = 0; i < keys.length; i++) {
										if ((!GTFObject.GENE_ID_TAG.equals(keys[i]))&&
												(!GTFObject.TRANSCRIPT_ID_TAG.equals(keys[i]))&&
												(!GTFObject.EXON_ID_TAG.equals(keys[i])))	{		// TODO make it nice: protein_id "...", other_id "..." 
												String putProtID= obj.getAttribute((String) keys[i]);
												String[] idTokens= putProtID.split(" ");
												for (int j = 0; j < idTokens.length; j++) 
													trpt.getTranslations()[0].addProteinID(idTokens[j]);
										}
									}
									
								} else {
									int y;
									for (y = 0; readFeatures!= null&& y < readFeatures.length; y++) 
										if (readFeatures[y].equalsIgnoreCase(tokens[2]))
											break;
									if (readFeatures== null|| y< readFeatures.length) {
										String regID= tokens[2];
										DirectedRegion newReg= new DirectedRegion(
												obj.getStart(),
												obj.getEnd(),
												trpt.getStrand()
										);
										newReg.setID(obj.getAttribute(tokens[2]+"_id"));
										newReg.setScore(obj.getScore());
										newReg.setChromosome(trpt.getChromosome());
										DirectedRegion[] reg= null;
										if (trpt.getAttributes()!= null)
											reg= (DirectedRegion[]) trpt.getAttributes().remove(regID);
										if (reg== null)
											reg= new DirectedRegion[] {newReg};
										else {
											int i;
											for (i = 0; i < reg.length; i++) {
												if (reg[i].getID().equals(newReg.getID())) {
													if (reg[i].getScore()!= newReg.getScore())
														System.out.println("WARNING: regions with non-matching scores "+newReg+", "+reg[i]);
													reg[i].set5PrimeEdge(Math.min(newReg.get5PrimeEdge(), reg[i].get5PrimeEdge()));
													reg[i].set3PrimeEdge(Math.max(newReg.get3PrimeEdge(), reg[i].get3PrimeEdge()));
													break;
												}
											}
											if (i== reg.length)
												reg= (DirectedRegion[]) Arrays.add(reg, newReg);
										}
										
										trpt.addAttribute(regID, reg);	// "domain", 
									}
										
									
								}
					
							}
							buffy.close();
						} catch (IOException e) {
							if (!silent)
								System.err.println("line "+ linesRead);
							e.printStackTrace();
						}
					
						
						if (readGene) {
							this.genes= (Gene[]) Arrays.toField(geneV);	// doesnt matter which sorting, no?!
							clustered= false;
						}
						System.gc();
						if (readGTF)
							this.gtfObj= (GTFObject[]) Arrays.toField(gtfV);
						
						if (bytesRead== size&& !silent)
							System.out.println();
						
					}

	private Gene[] clusterLoci(Gene[] g) {
		if (g== null)
			return null;
		
		Comparator compi= new DirectedRegion.StartComparator();
		java.util.Arrays.sort(g, compi);
		
		Vector v= new Vector();
		for (int i = 0; i < g.length; i++) {
			v.add(g[i]);
			int j;
			for (j = i+1; j < g.length; j++) {
				if (g[i].overlaps(g[j]))
					g[i].merge(g[j]);
				else 
					break;
			}		
			i= j-1;
		}
		
		return (Gene[]) Arrays.toField(v);
	}

	/**
			 * Requires the gtf file to be ordered by chromosomes AND that
			 * data (CDS, exon) from the same transcript is in adjacent lines
			 * (which is good for RefSeq in worm, where the same tID/gID exists
			 * on 4 diferent chromosomes). 
			 * It does NOT require any order of the transcript starts within 
			 * the chromosome, clustering is performed exhaustively.
			 * 
			 * retrieves every transcript that is overlapping a certain region of
			 * <code>ovlReg</code>
			 * @return
			 */
			/*
			 * Note that the clustering cannot be separated from the dynamic
			 * gene building, since reduncancy checks for splice sites and exons
			 * within the genes depend on it.
			 */
			public void read_transcript_clustering_while_reading_iefficient() throws Exception { 
							
			
					BufferedReader buffy= null;
					long size= 0l;
					
					if (fPath!= null&& fName!= null) {
						try {
							File file= new File(fPath+ File.separator+ fName);
							size= file.length();
							inputStream= new FileInputStream(file);
							inputStream.skip(bytesRead);
							buffy= new BufferedReader(new InputStreamReader(inputStream));
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
						
					
					Gene[][] geneReg= new Gene[2][];
					Comparator startCompi= new AbstractRegion.StartComparator();
					Comparator endCompi= new AbstractRegion.StartComparator();
						
					Vector gtfV= null;
					if (readGTF)
						gtfV = new Vector();
					try {
						
						String line= null;
						
						String lastTID= null;			
						String lastChrID= null;
						Transcript trpt= null;
						//boolean merged= true;	// for skiping insertion of init null trpt/gene
						while (true) {
							
								// read line
							buffy.mark(MAX_GTF_LINE_LENGTH);	// for skipping chromosomes
							line= buffy.readLine();
							if (line== null) {
								if (readGene&& trpt!= null) 
									mergeGene(geneReg, trpt.getGene(), startCompi, endCompi);
								
								
								
									// print statistics
								if (!silent) {
									System.out.println("\nRead:\t"+readChrs+" chromosomes, "+readGenes+" genes, "+readTranscripts+" transcripts, "+
											readExons+" exons.");
									System.out.println("Skip:\t"+skippedChr.size()+" chromosomes, "+skippedTranscripts.size()+" transcripts, "+
											skippedFeatures.size()+" features, "+skippedObjects+" objects.");
									System.out.print("Chromosomes:\t");
									for (int i = 0; i < skippedChr.size(); i++) 
										System.out.print(skippedChr.elementAt(i)+" ");
									System.out.print("\nTranscripts:\t");
									for (int i = 0; i < skippedTranscripts.size(); i++) 
										System.out.print(skippedTranscripts.elementAt(i)+" ");
									System.out.print("\nFeatures:\t");
									for (int i = 0; i < skippedFeatures.size(); i++) 
										System.out.print(skippedFeatures.elementAt(i)+" ");
									System.out.println();
								}
								break;
							}
							++linesRead;
							bytesRead+= line.length()+ 1;
								
							int perc= (int) ((bytesRead* 10d)/ size);
							if (perc> lastPerc&& !silent) {
								++lastPerc;
								System.out.print("*");
								System.out.flush();
							}
														
							String[] tokens= line.split("\t");	// must be tab, see specification
							if (!silent&& tokens.length< 8)
								System.err.println("line "+ linesRead+ ": skipped (<8 token)!\n\t"+ line);
							
								
							int x;		// skip mRNA, gene... 
							for (x = 0; readFeatures!= null&& x < readFeatures.length; x++) {
								if (readFeatures[x].equalsIgnoreCase(tokens[2]))
									break;
							}
							if (readFeatures!= null&& x== readFeatures.length) {
								skippedFeatures= Arrays.addUnique(skippedFeatures, tokens[2]);
								continue;
							}
							
							int start= Integer.parseInt(tokens[3]);
							int end= Integer.parseInt(tokens[4]);
							int strand= 0;
							strand= GTFObject.parseStrand(tokens[6].trim());
							if (strand== 0) {
								if (!silent)
									System.err.println("No strand assignment, line "+linesRead);
								continue;
							}						
							
							// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
							GTFObject obj= new GTFObject();	// build gtf-obj, handling of attributes easier
							obj.setStart(start);
							obj.setEnd(end);
							obj.setStrand(strand);
							obj.setFeature(tokens[2]);
							obj.setSource(tokens[1]);
							obj.setSeqname(tokens[0]);
							obj.setScore(tokens[5]);
							obj.setFrame(tokens[7]);
							String[] attrTokens= new String[0];
							if (tokens.length> 8) {
								String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
								attrTokens= h.split(";");
								for (int i = 0; i < attrTokens.length; i++) {
									h= attrTokens[i].trim();
									h= h.replaceAll("\\s+", " ");
									int sep= Math.max(h.indexOf(' '), h.indexOf('='));	// = in ASD gff
									if (sep < 0) 						// comments
										break;
									if (sep>= 0) {
										String id= h.substring(0, sep);
										String val= h.substring(sep+ 1, h.length());
										if (val.charAt(0)== '\"')
											val= val.substring(1, val.length()- 1);
										obj.addAttribute(id, val);
									}
								}
							}

								// filter gtfs that have ID...
							if (filtSomeIDs!= null) {
								boolean found= false;
								for (int i = 0; i < attrTokens.length; i++) {
									String h= attrTokens[i].trim(); 
									int sep= h.indexOf(' ');
									if (sep < 0) 						// comments
										break;
									if (sep>= 0) {
										String val= h.substring(sep+1, h.length()).trim();
										if (val.startsWith("\""))
											val= val.substring(1, val.length()-1);
										String[] t= val.split(" ");
										for (int j = 0; j < t.length; j++) {
											Integer in= (Integer) filtSomeIDs.get(t[j]);
											if (in== null)
												continue;	// ID not in filter list
											filtSomeIDs.put(t[j], new Integer(in.intValue()+ 1));
											found= true;
										}
									}
								}
								if (!found) {
									skippedObjects++;
									continue;
								}
							}
							
								// check chromosome
							String chrID= tokens[0];
							if (fName.contains("cow")&& chrID.equals("30"))
								chrID= "X";
							if (chromosomeWise&& !chrID.equals(lastChrID)) {
								if (lastChrID== null) {
									// neg chromosome filtering
									if (noIDs!= null) {
										int u;
										String chrIDu= chrID.toUpperCase();
										for (u = 0; u < noIDs.length; u++) {
											if (Pattern.matches(noIDs[u], chrIDu)) 
												break;
										}
										if (u< noIDs.length) {
											--linesRead;
											bytesRead-= line.length()+ 1;
											buffy.reset();
											if (skippedChr== null)
												skippedChr= new Vector();
											skippedChr= Arrays.addUnique(skippedChr, chrID);
											skipToNextChromosome(buffy, size);
											continue;
										}
									}
										// pos chromosome filtering
									if (filtChrIDs!= null) {
										int u;
										for (u = 0; u < filtChrIDs.length; u++) {
											if (filtChrIDs[u].equalsIgnoreCase(chrID))
												break;
										}
										if (u== filtChrIDs.length) {
											--linesRead;
											bytesRead-= line.length()+ 1;
											buffy.reset();
											if (skippedChr== null)
												skippedChr= new Vector();
											skippedChr= Arrays.addUnique(skippedChr, chrID);
											skipToNextChromosome(buffy, size);
											continue;
										}
									}
									
									lastChrID= chrID;
									for (int i = 0; readChr!= null&& i < readChr.size(); i++) 
										if (readChr.elementAt(i).equals(chrID)) {
											if (!silent) {
												System.err.println("Chromosome "+chrID+" already read! ("+linesRead+")");
												System.err.println(line);
											}
											return;
										}
									
									
								} else {
									//if (!merged) 
										mergeGene(geneReg, trpt.getGene(), startCompi, endCompi);
									if (readChr== null)
										readChr= new Vector();
									readChr.add(lastChrID);
									--linesRead;
									bytesRead-= line.length()+ 1;
									break;
								}
								
								readChrs++;
							}

							/*	--- END OF FILTERING --- */
							
							if (readGTF) 
								gtfV.add(obj);
							if (!readGene)
								continue;
							
								// get GID, TID
							String tid= obj.getAttribute(transcript_id_tag);
							if (tid== null) {
								if (!silent)
									System.err.println("No TID found.");
								continue;
							} 
				
								// new transcript?
							if (!tid.equals(lastTID)) {
								boolean overlap= true;
								if (filtRegs!= null) {	// should be ok for the moment. maybe also check in merging step
									int i;
									for (i = 0; i < filtRegs.length; i++) 
										if (filtRegs[i].overlaps(trpt))
											break;
									if (i== filtRegs.length)
										overlap= false;
								}
									// pull identical exon attributes up to transcript level
								if (trpt!= null&& trpt.getExons()[0].getAttributes()!= null) {
									Object[] keys= trpt.getExons()[0].getAttributes().keySet().toArray();
									for (int i = 0; i < keys.length; i++) {
										int j;
										for (j = 1; j < trpt.getExons().length; j++) {
											if (!trpt.getExons()[j].getAttribute(keys[i]).equals(trpt.getExons()[0].getAttribute(keys[i])))
												break;
										}
										if (j== trpt.getExons().length) {
											Object val= trpt.getExons()[0].getAttribute(keys[i]);
											for (int k = 0; k < trpt.getExons().length; k++) 
												trpt.getExons()[k].getAttributes().remove(keys[i]);
											trpt.addAttribute(keys[i], val);
										}
									}
								}
								
	
								// save gene
								if (overlap) {	
									//if (!merged) 	// if last transcript's gene couldnt be merged with anything
									if (trpt!= null)
										mergeGene(geneReg, trpt.getGene(), startCompi, endCompi);
									readGenes++;
								} else {
									skippedTranscripts.add(lastTID);
								}
								//merged= false;
								Gene ge= new Gene(Gene.getUniqueID());
								ge.setStrand(strand);
								ge.setChromosome(tokens[0]);
								ge.setSpecies(species);
								trpt= new Transcript(tid);
								lastTID= tid;
								trpt.setStrand(strand);
								trpt.setChromosome(tokens[0]);
								trpt.setSource(tokens[1]);
								ge.addTranscript(trpt);
								readTranscripts++;
							}
													
							if (obj.getFeature().equalsIgnoreCase("exon")) {
								String exonID= obj.getAttribute(GTFObject.EXON_ID_TAG);
								Exon e= new Exon(trpt, exonID, start, end);
								e.setChromosome(chrID);
								e.setStrand(strand);
								++readExons;
								
									// check if gene merge, not necessary when transcript starts are sorted ascending
								//trpt.addExon(e);		// automatically now in Exon constructor	
//								if (geneReg!= null) {
//									int p= java.util.Arrays.binarySearch(geneReg, e, startCompi);	// 
//									if (p>= 0) {
//										if (geneReg[p]!= trpt.getGene()) {
//											if (merged) {
//												trpt.getGene().merge(geneReg[p]);
//												geneReg= (Gene[]) Arrays.remove(geneReg, p);
//											} else {
//												geneReg[p].merge(trpt.getGene());
//												merged= true;
//											}
//										}
//									} else {	// between two
//										p= -(p+1);	// insertion point
//										int q= p-1;
//										if (q>= 0&& geneReg[q]!= trpt.getGene()&&
//												e.overlaps(geneReg[q])) { 
//											if (merged) {	// if already merged, one has to be removed from array
//												trpt.getGene().merge(geneReg[q]);
//												geneReg= (Gene[]) Arrays.remove(geneReg, q);
//											} else
//												geneReg[q].merge(trpt.getGene());
//											merged= true;
//										}
//										q= p;	//+1; NO, insertion point is the upper neighbor
//										if (q< geneReg.length&& geneReg[q]!= trpt.getGene()&&
//												e.overlaps(geneReg[q])) { 
//											if (merged) {
//												trpt.getGene().merge(geneReg[q]);
//												geneReg= (Gene[]) Arrays.remove(geneReg, q);
//											} else
//												geneReg[q].merge(trpt.getGene());
//											merged= true;
//										}
//									}
//								}
								
								
								
							} else if (tokens[2].equalsIgnoreCase("CDS")) {
								trpt.addCDS(start, end);
								Object[] keys= obj.getAttributes().keySet().toArray();
								for (int i = 0; i < keys.length; i++) {
									if ((!GTFObject.GENE_ID_TAG.equals(keys[i]))&&
											(!GTFObject.TRANSCRIPT_ID_TAG.equals(keys[i]))&&
											(!GTFObject.EXON_ID_TAG.equals(keys[i])))	{		// TODO make it nice: protein_id "...", other_id "..." 
											String putProtID= obj.getAttribute((String) keys[i]);
											String[] idTokens= putProtID.split(" ");
											for (int j = 0; j < idTokens.length; j++) 
												trpt.getTranslations()[0].addProteinID(idTokens[j]);
									}
								}
								
							} else {
								int y;
								for (y = 0; readFeatures!= null&& y < readFeatures.length; y++) 
									if (readFeatures[y].equalsIgnoreCase(tokens[2]))
										break;
								if (readFeatures== null|| y< readFeatures.length) {
									String regID= tokens[2];
									DirectedRegion newReg= new DirectedRegion(
											obj.getStart(),
											obj.getEnd(),
											trpt.getStrand()
									);
									newReg.setID(obj.getAttribute(tokens[2]+"_id"));
									newReg.setScore(obj.getScore());
									newReg.setChromosome(trpt.getChromosome());
									DirectedRegion[] reg= null;
									if (trpt.getAttributes()!= null)
										reg= (DirectedRegion[]) trpt.getAttributes().remove(regID);
									if (reg== null)
										reg= new DirectedRegion[] {newReg};
									else {
										int i;
										for (i = 0; i < reg.length; i++) {
											if (reg[i].getID().equals(newReg.getID())) {
												if (reg[i].getScore()!= newReg.getScore())
													System.out.println("WARNING: regions with non-matching scores "+newReg+", "+reg[i]);
												reg[i].set5PrimeEdge(Math.min(newReg.get5PrimeEdge(), reg[i].get5PrimeEdge()));
												reg[i].set3PrimeEdge(Math.max(newReg.get3PrimeEdge(), reg[i].get3PrimeEdge()));
												break;
											}
										}
										if (i== reg.length)
											reg= (DirectedRegion[]) Arrays.add(reg, newReg);
									}
									
									trpt.addAttribute(regID, reg);	// "domain", 
								}
									
								
							}
				
						}
						buffy.close();
					} catch (IOException e) {
						if (!silent)
							System.err.println("line "+ linesRead);
						e.printStackTrace();
					}
				
					
					if (bytesRead== size&& !silent)
						System.out.println();
					
					if (readGene)
						this.genes= geneReg[0];	// doesnt matter which sorting, no?!
					if (readGTF)
						this.gtfObj= (GTFObject[]) Arrays.toField(gtfV);
				}

	/**
	 * per definition there cannot 2 genes with same start in the arrays (merged!) 
	 * therefore binary search works
	 * 
	 * not so sure, during the multiple merges there can happen the case that two have
	 * the same 
	 * 
	 * @param geneReg
	 * @param gene
	 * @param startCompi
	 * @param endCompi
	 * @return
	 */
	private boolean mergeGene(Gene[][] geneReg, Gene gene, Comparator startCompi, Comparator endCompi) {
		boolean merged= false;
		if (geneReg[0]== null) { 
			geneReg[0]= new Gene[] {gene};
			geneReg[1]= geneReg[0];	// not needed here, see efficient version..
		} else {
			Vector vRest= new Vector(geneReg[0].length);
			for (int i = 0; i < geneReg[0].length; i++) {
				if (geneReg[0][i].overlaps(gene)) 
					gene.merge(geneReg[0][i]);
				else
					vRest.add(geneReg[0][i]);
			}
			vRest.add(gene);
			geneReg[0]= (Gene[]) Arrays.toField(vRest);
			geneReg[1]= geneReg[0];	// not needed here, see efficient version..
		}
	
		return merged;
	}

	/**
	 * per definition there cannot 2 genes with same start in the arrays (merged!) 
	 * therefore binary search works
	 * 
	 * not so sure, during the multiple merges there can happen the case that two have
	 * the same 
	 * 
	 * @param geneReg
	 * @param gene
	 * @param startCompi
	 * @param endCompi
	 * @return
	 * @deprecated http://genome.ucsc.edu/cgi-bin/hgTracks?position=chr5:140207541-140372113&hgsid=95135519&refGene=pack&hgFind.matches=NM_031857,
	 */
	private boolean mergeGene_efficient_notwork(Gene[][] geneReg, Gene gene, Comparator startCompi, Comparator endCompi) {
		boolean merged= false;
		if (geneReg[0]== null) { 
			geneReg[0]= new Gene[] {gene};
			geneReg[1]= new Gene[] {gene};
		} else {
			
				// check forward in start array for overlaps
			int p= java.util.Arrays.binarySearch(geneReg[0], gene, startCompi);			
			if (p< 0)
				p= -(p+1);
			
			for (int i = p; i < geneReg[0].length; i++) {
				if (geneReg[0][i]== gene)
					continue;
				if (!geneReg[0][i].overlaps(gene))
					break;
				if (merged) {	// both in array, merge step doesnt affect ordering
					gene.merge(geneReg[0][i]);
					
					int q= java.util.Arrays.binarySearch(geneReg[1], geneReg[0][i], endCompi);
					geneReg[1]= (Gene[]) Arrays.remove(geneReg[1], q);
					geneReg[0]= (Gene[]) Arrays.remove(geneReg[0], i--);
					// gene has not to be inserted, because already in array
					
				} else {	// new start/end in array
					geneReg[0][i].merge(gene);	//.. for start it is the next one in the array, nothing happend
					java.util.Arrays.sort(geneReg[1],endCompi);	// .. but maybe end order out of control now..
					gene= geneReg[0][i];
					merged= true;
				}
			}
			
			// .. and backward in reverse
			int r= java.util.Arrays.binarySearch(geneReg[1], gene, endCompi);
			if (r< 0)
				r= -(r+1);
			
			int s= r;
			if (r== geneReg[1].length)
				--s;
			for (int i = s; i >= 0; --i) {
				if (geneReg[1][i]== gene)
					continue;
				if (!geneReg[1][i].overlaps(gene))
					break;
				if (merged) {					
					gene.merge(geneReg[1][i]);
					
					int q= java.util.Arrays.binarySearch(geneReg[0], geneReg[1][i], startCompi);					
					geneReg[0]= (Gene[]) Arrays.remove(geneReg[0], q);					
					geneReg[1]= (Gene[]) Arrays.remove(geneReg[1], i++);	// re-iterate i
					
				} else {
					geneReg[1][i].merge(gene);
					java.util.Arrays.sort(geneReg[0],startCompi);	
					gene= geneReg[1][i];
					merged= true;
				}
			}
			
				// nowhere merged
			if (!merged) {
				//p= java.util.Arrays.binarySearch(geneReg[0], gene, startCompi);	// can have changed in size, NO: no merge ops
				geneReg[0]= (Gene[]) Arrays.insert(geneReg[0], gene, p);
				//r= java.util.Arrays.binarySearch(geneReg[1], gene, endCompi);
				geneReg[1]= (Gene[]) Arrays.insert(geneReg[1], gene, r);
			}

		}

		return merged;
	}

	/**
		 * Requires the gtf file to be ordered by chromosomes AND that
		 * data (CDS, exon) from the same transcript is in adjacent lines
		 * (which is good for RefSeq in worm, where the same tID/gID exists
		 * on 4 diferent chromosomes). 
		 * It does NOT require any order of the transcript starts within 
		 * the chromosome, clustering is performed exhaustively.
		 * 
		 * retrieves every transcript that is overlapping a certain region of
		 * <code>ovlReg</code>
		 * @return
		 */
		/*
		 * Note that the clustering cannot be separated from the dynamic
		 * gene building, since reduncancy checks for splice sites and exons
		 * within the genes depend on it.
		 */
		public Gene[] readNextChromosome_raf() { 
						
		
				RandomAccessFile raf= null;
				long size= 0l;
				if (fPath!= null&& fName!= null) {
					try {
						raf= new RandomAccessFile(new File(fPath+ File.separator+ fName),"r");
						size= raf.length();
					} catch (Exception e) {
						e.printStackTrace();
					}
				} else if(inputStream!= null)
					System.err.println(this.getClass()+" only supports input from Files.");
				
				Gene[] geneReg= null;
				Comparator startCompi= new AbstractRegion.StartComparator();
					
				try {
					raf.seek(bytesRead);
					
					String line;
					
					String lastTID= null;			
					String lastChrID= null;
					Transcript trpt= null;
					boolean merged= true;	// for skiping insertion of init null trpt/gene
					while (true) {
						
							// read line
						line= raf.readLine();
						if (line== null) 
							break;
						++linesRead;
						bytesRead+= line.length()+ 1;
							
						int perc= (int) ((bytesRead* 10d)/ size);
						if (perc> lastPerc&& !silent) {
							++lastPerc;
							System.out.print("*");
							System.out.flush();
						}
						
						String[] tokens= line.split("\t");	// must be tab, see specification
						if (tokens.length< 8)
							System.err.println("line "+ linesRead+ ": skipped (<8 token)!\n\t"+ line);
						int start= Integer.parseInt(tokens[3]);
						int end= Integer.parseInt(tokens[4]);
						int strand= 0;
						strand= GTFObject.parseStrand(tokens[6].trim());
						if (strand== 0) {
							System.err.println("No strand assignment, line "+linesRead);
							continue;
						}						
						
						// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
						String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
						String[] attrTokens= h.split(";");
	
							// skip mRNA, gene...
						if (!tokens[2].equals("exon")&& !tokens[2].equals("CDS"))
							continue;
						
							// check chromosome
						String chrID= tokens[0];
						if (fName.contains("cow")&& chrID.equals("30"))
							chrID= "X";
						if (!chrID.equals(lastChrID)) {
							String[] noIDs= new String[] {
									"M", "U", "RANDOM", "HAPLO"
							};
							int u;
							for (u = 0; u < noIDs.length; u++) {
								if (chrID.toUpperCase().contains(noIDs[u])) 
									break;
							}
							if (u< noIDs.length) {
								//skipToNextChromosome(raf);
								continue;
							}
							
								// pos chromosome filtering
							if (filtChrIDs!= null) {
								for (u = 0; u < filtChrIDs.length; u++) {
									if (filtChrIDs[u].equalsIgnoreCase(chrID))
										break;
								}
								if (u== filtChrIDs.length) {
									//skipToNextChromosome(raf);
									continue;
								}
							}
							 
							if (lastChrID== null) {
								lastChrID= chrID;
								for (int i = 0; readChr!= null&& i < readChr.size(); i++) 
									if (readChr.elementAt(i).equals(chrID)) {
										System.err.println("Chromosome "+chrID+" already read! ("+linesRead+")");
										System.err.println(line);
										return null;
									}
							} else {
								if (readChr== null)
									readChr= new Vector();
								readChr.add(lastChrID);
								--linesRead;
								bytesRead-= line.length()+ 1;
								break;
							}
						}
						
							// get GID, TID
						String tid= null;
						String gid= null;
						for (int i = 0; i < attrTokens.length; i++) {
							h= attrTokens[i].trim();
							int sep= Math.max(h.indexOf(' '), h.indexOf('='));	// = in ASD gff
							if (sep < 0) 						// comments
								break;
							if (sep>= 0) {
								String id= h.substring(0, sep);
								if (id.equals(transcript_id_tag)) {
									tid= h.substring(sep+ 1, h.length());
									if (tid.charAt(0)== '\"')
										tid= tid.substring(1, tid.length()- 1);
								}
							}
						}
						if (tid== null) {
							System.err.println("No TID found.");
							continue;
						}
			
							// new transcript?
						if (!tid.equals(lastTID)) {
							boolean overlap= true;
							if (filtRegs!= null) {	// should be ok for the moment. maybe also check in merging step
								int i;
								for (i = 0; i < filtRegs.length; i++) 
									if (filtRegs[i].overlaps(trpt))
										break;
								if (i== filtRegs.length)
									overlap= false;
							}

							// save gene
							if (overlap&& (!merged)) {	// if last transcript's gene couldnt be merged with anything
								if (geneReg== null)
									geneReg= new Gene[] {trpt.getGene()};
								else {
									int p= java.util.Arrays.binarySearch(geneReg, trpt.getGene(), startCompi);
									geneReg= (Gene[]) Arrays.insert(geneReg, trpt.getGene(), p);
								}
							}
							merged= false;
							Gene ge= new Gene(Gene.getUniqueID());
							ge.setStrand(strand);
							ge.setChromosome(tokens[0]);
							trpt= new Transcript(tid);
							lastTID= tid;
							trpt.setStrand(strand);
							trpt.setChromosome(tokens[0]);
							trpt.setSource(tokens[1]);
							ge.addTranscript(trpt);
						}
												
						if (tokens[2].equalsIgnoreCase("exon")) {
							String exonID= null;
							for (int i = 0; i < attrTokens.length; i++) {
								h= attrTokens[i];
								int sep= h.indexOf(' ');
								if (sep < 0) 						// comments
									break;
								if (sep>= 0) {
									String id= h.substring(0, sep);
									if (id.equals(GTFObject.EXON_ID_TAG))
										exonID= h.substring(sep+ 1, h.length());
									break;
								}
							}
							Exon e= new Exon(trpt, exonID, start, end);
							e.setChromosome(chrID);
							e.setStrand(strand);
							
								// check if gene merge, not necessary when transcript starts are sorted ascending
							if (geneReg!= null) {
								int p= java.util.Arrays.binarySearch(geneReg, e, startCompi);
								if (p>= 0) {
									if (geneReg[p]!= trpt.getGene()) {
										if (merged) {
											trpt.getGene().merge(geneReg[p]);
											geneReg= (Gene[]) Arrays.remove(geneReg, p);
										} else {
											geneReg[p].merge(trpt.getGene());
											merged= true;
										}
									}
								} else {	// between two
									p= -(p+1);	// insertion point
									int q= p-1;
									if (q>= 0&& geneReg[q]!= trpt.getGene()&&
											e.overlaps(geneReg[q])) { 
										if (merged) {	// if already merged, one has to be removed from array
											trpt.getGene().merge(geneReg[q]);
											geneReg= (Gene[]) Arrays.remove(geneReg, q);
										} else
											geneReg[q].merge(trpt.getGene());
										merged= true;
									}
									q= p;	//+1; NO, insertion point is the upper neighbor
									if (q< geneReg.length&& geneReg[q]!= trpt.getGene()&&
											e.overlaps(geneReg[q])) { 
										if (merged) {
											trpt.getGene().merge(geneReg[q]);
											geneReg= (Gene[]) Arrays.remove(geneReg, q);
										} else
											geneReg[q].merge(trpt.getGene());
										merged= true;
									}
								}
							}
							
							trpt.addExon(e);
							
						} else if (tokens[2].equalsIgnoreCase("CDS")) 
							trpt.addCDS(start, end);
			
					}
					raf.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			
				
				if (geneReg== null&& !silent)
					System.out.println();
//				else {
//					System.out.print("*");
//					System.out.flush();
//				}
				return geneReg;
			}

	/**
	 * Requires the gtf file to be ordered by chromosomes
	 * AND ascending transcript starts (both strands mixed or
	 * separated).
	 * @deprecated
	 */
	public Gene readNextGene() {
				
		RandomAccessFile raf= null;
		long size= 0l;
		if (fPath!= null&& fName!= null) {
			try {
				raf= new RandomAccessFile(new File(fPath+ File.separator+ fName),"r");
				size= raf.length();
			} catch (Exception e) {
				e.printStackTrace();
			}
		} else if(inputStream!= null)
			System.err.println(this.getClass()+" only supports input from Files.");
		
		readChr= new Vector();
		Gene[] geneReg= null;
		Comparator startCompi= new AbstractRegion.StartComparator();
		try {
			raf.seek(bytesRead);
			
			String line;
			int lineCtr= 0;
			
			String lastTID= null;			
			String lastChrID= null;
			Transcript trpt= null;
			while (true) {
				
					// read line
				line= raf.readLine();
				if (line== null)
					break;
					
				bytesRead+= line.length()+ 1;
				int perc= (int) ((bytesRead* 10d)/ size);
				if (perc> lastPerc) {
					++lastPerc;
					System.out.print("*");
					System.out.flush();
				}
				++linesRead;
				
				String[] tokens= line.split("\t");	// must be tab, see specification
				if (tokens.length< 8)
					System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
				int start= Integer.parseInt(tokens[3]);
				int end= Integer.parseInt(tokens[4]);
				int strand= 0;
				if (tokens[6].trim().equals("+"))
					strand= 1;
				else if (tokens[6].trim().equals("-"))
					strand= -1;
				else {
					System.err.println("No strand assignment, line "+linesRead);
					continue;
				}						
				
				// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
				String h= line.substring(line.indexOf(tokens[8]), line.length()).trim();	// attributes, comments
				String[] attrTokens= h.split(";");
				
					// check chromosome
				String chrID= tokens[0];
				if (!chrID.equals(lastChrID)) {
					if (lastChrID== null) {
						lastChrID= chrID;
						for (int i = 0; i < readChr.size(); i++) 
							if (readChr.elementAt(i).equals(chrID)) {
								System.err.println("Chromosome "+chrID+" already read! ("+lineCtr+")");
								return null;
							}
						readChr.add(chrID);
					} else
						break;
				}
				
					// get GID, TID
				String tid= null;
				String gid= null;
				for (int i = 0; i < attrTokens.length; i++) {
					h= attrTokens[i].trim();
					int sep= h.indexOf(' ');
					if (sep < 0) 						// comments
						break;
					if (sep>= 0) {
						String id= h.substring(0, sep);
						if (id.equals(GTFObject.TRANSCRIPT_ID_TAG)) {
							tid= h.substring(sep+ 1, h.length());
							if (tid.charAt(0)== '\"')
								tid= tid.substring(1, tid.length()- 1);
						}
					}
				}
				if (tid== null) {
					System.err.println("No TID found.");
					continue;
				}
	
					// new transcript?
				if (!tid.equals(lastTID)) {
					Gene ge= new Gene(Gene.getUniqueID());
					ge.setStrand(strand);
					ge.setChromosome(tokens[0]);
					trpt= new Transcript(tid);
					lastTID= tid;
					trpt.setStrand(strand);
					trpt.setChromosome(tokens[0]);
					trpt.setSource(tokens[1]);
					ge.addTranscript(trpt);
				}
										
				if (tokens[2].equalsIgnoreCase("exon")) {
					String exonID= null;
					for (int i = 0; i < attrTokens.length; i++) {
						h= attrTokens[i];
						int sep= h.indexOf(' ');
						if (sep < 0) 						// comments
							break;
						if (sep>= 0) {
							String id= h.substring(0, sep);
							if (id.equals(GTFObject.EXON_ID_TAG))
								exonID= h.substring(sep+ 1, h.length());
							break;
						}
					}
					Exon e= new Exon(trpt, exonID, start, end);
						// check if gene merge, not necessary when transcript starts are sorted ascending
//					if (trpt.getGene().getExons()!= null) {
//						int p= java.util.Arrays.binarySearch(geneReg, e, startCompi);
//						if (p>= 0) {
//							geneReg[p].merge(trpt.getGene());
//							geneReg= (Gene[]) Arrays.remove(geneReg, geneReg.length- 1);
//						} else {
//							p= -(p+1);	// insertion point
//							if (trpt.isForward())
//								--p;
//							if (p>= 0&& Math.abs(geneReg[p].getEnd())>= Math.abs(e.getEnd())) {
//								geneReg[p].merge(trpt.getGene());
//								geneReg= (Gene[]) Arrays.remove(geneReg, geneReg.length- 1);
//							}
//						}
//					}
					int chkStart= Math.abs(trpt.getGene().getStart());
					trpt.addExon(e);
					if (Math.abs(trpt.getGene().getStart())< chkStart)
						System.err.println("Transcripts/Gene not sorted ascending.");
						// insert gene when first exon inserted (region inited!) 
					if (trpt.getGene().getExons().length== 1) {	
						if (geneReg== null)
							geneReg= new Gene[] {trpt.getGene()};
						else {
							int p= java.util.Arrays.binarySearch(geneReg, trpt.getGene(), startCompi);
							if (p>= 0) 
								geneReg[p].merge(trpt.getGene());	// not in, not remove
							else {
								int q= -(p+1);
								if (trpt.isForward())
									--q;
								if (q>= 0&& Math.abs(geneReg[q].getEnd())>= Math.abs(trpt.getGene().getStart())) 
									geneReg[q].merge(trpt.getGene());
								else 
									geneReg= (Gene[]) Arrays.insert(geneReg, trpt.getGene(), p);
							}
						}
					}
				} else if (tokens[2].equalsIgnoreCase("CDS")) 
					trpt.addCDS(start, end);
	
			}
			raf.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		return null;	// geneReg;
	}

	public boolean isApplicable() {
				
		BufferedReader buffy= null;
		File f= new File(fPath+ File.separator+ fName);
		if (fPath!= null&& fName!= null)
			try {
				buffy= new BufferedReader(new FileReader(f));
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
		else 
			buffy= new BufferedReader(new InputStreamReader(inputStream));
		
		System.out.print("Checking format of file "+f.getAbsolutePath()+" ");
		System.out.flush();
		long bytesRead= 0l;
		long size= f.length();
		int perc= 0;
		String lastChrID= null;
		try {
			String line;
			int lineCtr= 0;
			Vector v= new Vector();
			while (buffy.ready()) {
				
				line= buffy.readLine();
				if (line== null)
					break;
				
				++lineCtr;
				bytesRead+= line.length()+ 1;
				if (bytesRead*10/size> perc) {
					System.out.print("*");
					System.out.flush();
					++perc;
				}
				String[] tokens= line.split("\t");	// must be tab, see specification
				if (tokens.length< 8) {
					System.err.println("Warning: line "+ lineCtr+ ": has <8 token !\n\t"+ line);
					return false;
				}
				
				if (!tokens[0].equals(lastChrID)) {
					for (int i = 0; i < v.size(); i++) 
						if (tokens[0].equals(v.elementAt(i))) {
							System.err.println("Error: "+tokens[0]+" already read (line "+lineCtr+")");
							return false;
						} 
							
					v.add(tokens[0]);
					lastChrID= tokens[0];
				}
			}
			buffy.close();
			System.out.println("\ttrue");
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		return true;
	}

	public void setFiltTrptIDs(String[] filtTrptIDs) {
		this.filtTrptIDs = filtTrptIDs;
	}

	public void setSilent(boolean silent) {
		this.silent = silent;
	}

	public String[] getFiltSomeIDs() {
		if (filtSomeIDs== null)
			return null;
		
		String[] result= new String[filtSomeIDs.size()];
		Object[] keys= filtSomeIDs.keySet().toArray();
		for (int i = 0; i < result.length; i++) 
			result[i]= (String) keys[i];
		return result;
	}

	public void setFiltSomeIDs(String[] filtSomeIDs) {
		this.filtSomeIDs= new HashMap(filtSomeIDs.length);
		for (int i = 0; i < filtSomeIDs.length; i++)  
			this.filtSomeIDs.put(filtSomeIDs[i].toUpperCase(), new Integer(0));
		
	}
	
	public void setReadAllLines() {
		filtChrIDs= null;
		filtGeneIDs= null;
		filtRegs= null;
		filtSomeIDs= null;
		filtTrptIDs= null;
		readFeatures= null;
		noIDs= null;
	}
	
	public String[] getFiltSomeIDsNotFound() {
		if (filtSomeIDs== null)
			return null;
		
		Vector v= new Vector();
		Object[] keys= filtSomeIDs.keySet().toArray();
		for (int i = 0; i < keys.length; i++) { 
			Integer cnt= (Integer) filtSomeIDs.get(keys[i]);
			if (cnt.intValue()== 0)
				v.add(keys[i]);
		}
		return (String[]) Arrays.toField(v);
	}

	public boolean isChromosomeWise() {
		return chromosomeWise;
	}

	public void setChromosomeWise(boolean chromosomeWise) {
		this.chromosomeWise = chromosomeWise;
	}

	public String[] getFiltChrIDs() {
		return filtChrIDs;
	}

	public void setFiltChrIDs(String[] filtChrIDs) {
		this.filtChrIDs = filtChrIDs;
	}

	public Gene[] getGenes() {
		if (genes== null)
			return null;
		if (!clustered) {
			for (int i = 0; i < genes.length; i++) 
				genes[i].setConstruct(false);
			genes= clusterLoci(genes);
			clustered= true;
			readGenes+= this.genes.length;
		}
		
		return genes;
	}

	public void setGenes(Gene[] genes) {
		this.genes = genes;
	}

	public DirectedRegion[] getFiltRegs() {
		return filtRegs;
	}

	public void setFiltRegs(DirectedRegion[] filtRegs) {
		this.filtRegs = filtRegs;
	}

	public boolean isReadGene() {
		return readGene;
	}

	public void setReadGene(boolean readGene) {
		this.readGene = readGene;
	}

	public boolean isReadGTF() {
		return readGTF;
	}

	public void setReadGTF(boolean readGTF) {
		this.readGTF = readGTF;
	}

	public String[] getReadFeatures() {
		return readFeatures;
	}

	public void setReadFeatures(String[] readFeatures) {
		this.readFeatures = readFeatures;
	}
	
	public void addReadFeature(String newFeature) {
		String[] r= new String[readFeatures.length+ 1];
		for (int i = 0; i < readFeatures.length; i++) {
			r[i]= readFeatures[i];
		}
		r[r.length- 1]= newFeature;
		readFeatures= r;
	}

	public int getReadGenes() {
		return readGenes;
	}

	public long getBytesRead() {
		return bytesRead;
	}
	
	public int getUnclusteredGeneNb() {
		if (this.genes== null)
			return 0;
		return this.genes.length;
	}

}
