	package gphase.io.gtf;

import gphase.algo.ASAnalyzer;
import gphase.graph.SpliceGraph;
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
import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
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


public class EncodeWrapper extends GTFWrapper {
	
	protected InputStream inputStream= null;
	String speName= null;
	String genomeVer= null;
	String buildVersion= null;
	int linesRead= 0;
	boolean silent= false;
	
	public Graph getGraph(boolean encode) {
		if (gtfObj== null) {
			try {
				read();
			} catch (Exception e) {
				e.printStackTrace(); 
			}
		}
		return assemble(encode);		// <===== check ENCODE here !!!
	}
	
	public static void main(String[] args) {
		//"encode/44regions_genes_CHR_coord.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		//"encode/EnsemblGenes_fromUCSC.gtf"
		
		//"encode/gencode_races.gtf"
		//"encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		// "encode/Sequences_mapped_HAVANA_136.gtf"
		String fName= "encode/RefSeqGenes_fromUCSC.gtf";
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		boolean encode= false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		
		Graph g= myWrapper.getGraph(encode);		// <===== check ENCODE here !!!
		//g.filterNonCodingTranscripts();
		//g.filterCodingTranscripts();
		//g.initTU();

//		PrintStream pr= null;
//		try {
//			pr = new PrintStream(new File("deg_ensembl.txt"));
//		} catch (FileNotFoundException e1) {
//			// TODO Auto-generated catch block
//			e1.printStackTrace();
//		}
//		ASAnalyzer.test04_determineVariations(g, pr);
		
		
//		ASVariation[][] classes= g.getASVariations(1);
//		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; j < classes[i].length; j++) {
//				if (classes[i][j].isProteinCoding_old_publ()&& !classes[i][j].isProteinCoding()) {
//					classes[i][j].isProteinCoding_old_publ();
//					classes[i][j].isProteinCoding();
//				}
//			}
//		}
		

			// multi variations
//		ASMultiVariation[][] vars= g.getASMultiVariations();
//		vars= (ASMultiVariation[][]) Arrays.sort2DFieldRev(vars);
//		for (int i = 0; i < vars.length; i++) {
//			System.out.println(vars[i].length+"\t"+vars[i][0]);
//		}
		
//		ASVariation[] vars= ASAnalyzer.getVariation("( // 1=2^)", g.getASVariations(ASMultiVariation.FILTER_NONE));
//		System.out.println(vars.length+ " events");
//		for (int i = 0; i < vars.length; i++) 
//			;//vars[i].outputDetail(System.out);
		
//		ASAnalyzer.check_AA_AD(g, true, false, false);
//		ASAnalyzer.check_AA_AD(g, true, true, false);
//		ASAnalyzer.check_AA_AD(g, false, false, false);
//		ASAnalyzer.check_AA_AD(g, false, true, false);
		
		ASAnalyzer.test01_clusters_coverage_as(g, System.out);
//		if (1== 1)
//			System.exit(0);
		
		
//		try {
//			PrintStream buffy= new PrintStream(new File("atg_aa_analysis_vars_refseq").getAbsolutePath());
//			ASAnalyzer.outputFirstExonIntronAtg2(g, buffy);
//			SpliceSite[] ss1= ASAnalyzer.getFirstExonIntronAtg3(g, null);
//			SpliceSite[] ss2= ASAnalyzer.getFirstExonIntronAtg4(g, null);
//			int dnrEq= 0, accEq= 0, dnrNe= 0, accNe= 0;
//			for (int i = 0; i < ss2.length; i++) {
//				int j;
//				for (j = 0; j < ss1.length; j++) {
//					if (ss2[i]== ss1[j]) {
//						if (ss2[i].isDonor())
//							dnrEq++;
//						else
//							accEq++;
//						break;
//					}						
//				}
//				if (j== ss1.length) {
//					if (ss2[i].isDonor())
//						dnrNe++;
//					else
//						accNe++;
//				}
//			}
//			System.out.println(dnrEq+","+accEq+" : "+dnrNe+","+accNe);
//			buffy.flush(); buffy.close();
//		} catch (Exception e) {
//			e.printStackTrace();
//		}
		

//		ASAnalyzer.outputVariations(new ASVariation[][] {ASAnalyzer.getVariation("( 1^3=// 2^4=)", g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY))}, false, false, System.out);
		// "output5UTR_REFSEQ"
		try {
			PrintStream p= new PrintStream("SSout_encode");
//			//ASAnalyzer.output5UTRLengthAnalysis(g, p);
//			//ASAnalyzer.outputInternalIntrons(g, p);
			ASAnalyzer.outputSSOutCdsUtr(g, p);
			p.flush();
			p.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
//		ASAnalyzer.getSylvainsSize(g, System.out);
//		ASAnalyzer.test01_clusters_coverage_as(g, System.out);
//		ASAnalyzer.test(new String[]{"test01_clusters_coverage_as",
//				"test02_ss_statistics",
//				"test03_lengthVariationModulo",
//				"test04_determineVariations"}, 
//				new Class[] {Graph.class, PrintStream.class}, 
//				new boolean[] {true, true, false, false});
//		ASAnalyzer.test(new String[]{"test02_ss_statistics"}, 
//				new Class[] {Graph.class, PrintStream.class}, 
//				new boolean[] {true});
//		ASAnalyzer.test02_ss_statistics(g, System.out);
//		ASAnalyzer.test02_ss_statistics_outCDSalt(g);
//		ASAnalyzer.test02b_ss_statistics_3P(g, System.out);
//		ASAnalyzer.test03_lengthVariationModulo(g, System.out);
//		ASAnalyzer.test03b_lengthVariationFirstExon(g, System.out);
//		ASAnalyzer.test04_determineVariations(g, System.out);
		
//		ASAnalyzer.check_AA_AD(g, true, false, false);
//		ASAnalyzer.check_AA_AD(g, true, true, false);
//		ASAnalyzer.check_AA_AD(g, true, false, true);
//		ASAnalyzer.check_AA_AD(g, true, true, true);
//		ASAnalyzer.check_AA_AD(g, false, false, false);
//		ASAnalyzer.check_AA_AD(g, false, true, false);
//		ASAnalyzer.check_AA_AD(g, false, false, true);
//		ASAnalyzer.check_AA_AD(g, false, true, true);
		
//		long[] r=ASAnalyzer.getCDSUTRSizeComplement(g);
//		System.out.println("5UTR "+r[0]+", CDS "+r[1]+", 3UTR "+r[2]);

		//		ASVariation[] var= ASAnalyzer.getVariation("(1= // 2=)", g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY));
//		for (int i = 0; i < var.length; i++) {
//			var[i].outputDetail(System.out);
//		}
		
//		ASVariation[][] as= ASAnalyzer.determineVariations(g, (PrintStream) null);
//		for (int i = 0; i < as.length; i++) {
//			System.out.println(as[i][0].toString()+"\t"+as[i][0].toBitString());
//		}
		
//		ASAnalyzer.lengthVariationModulo(g);
//		long t1= System.currentTimeMillis();
//		System.out.println("build up"+(t1- t0));
//		g.getASVariations(ASMultiVariation.FILTER_NONE);
//		long t2= System.currentTimeMillis();
//		System.out.println("filter: "+ (t2- t1));
//		g.getASVariations(ASMultiVariation.FILTER_REDUNDANT);
//		long t3= System.currentTimeMillis();
//		System.out.println("filter redundant: "+(t3- t2));
		
			// analyze1
//		System.out.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
//		ASAnalyzer.filterSingleExonClusters(g);
//		System.out.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
//		System.out.println("clusters w as: "+g.getClustersWithAS());
		
		//g.filterNonCodingTranscripts();
		//System.out.println(ASAnalyzer.getPercASSpliceSites(g));
		//ASVariation[][] classes= ASAnalyzer.determineVariations(g, "isTrue");
		//ASAnalyzer.lengthVariation(g);
		//ASAnalyzer.getSpecificVariations(g, System.out);
		//ASAnalyzer.getVincenco(g);
		//ASAnalyzer.getAllASVariation(g);
		//ASAnalyzer.determineVariations(g, System.out);
//		System.out.println(g.getTranscriptCount());
		
		
//		ASAnalyzer.analyzeGene(g, "RNPC2", "isTrue", ASMultiVariation.FILTER_REDUNDANT);
//		int cnt= 0, not= 0;
//		BufferedWriter buffy= new BufferedWriter(new FileWriter)
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; j < classes[i].length; j++) {
//				classes[i][j].outputDetail(System.out);
//			}
//		}
//		System.out.println(cnt+","+not+" ("+(cnt+not)+")");
		
		
//		ASVariation[][] classes= null;
//		ASVariation[] events= ASAnalyzer.getVariation("(1=4^ // 2=3^)", classes);
//		int cnt= 0;
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; events!= null&& j < events.length; j++) {
//				if (events[j].includesFirstSpliceSite()) {
//					events[j].outputDetail(System.out);
//				}
			//}
//		}
//		System.out.println("--");
				
		
//			events= ASAnalyzer.getVariation("(1=3^ // 2=4^)", classes);
////			int cnt= 0;
////			for (int i = 0; i < classes.length; i++) {
//				for (int j = 0; j < events.length; j++) {
////					if (events[j].includesFirstSpliceSite()) {
//						events[j].outputDetail(System.out);
//					}
				//}
//			}
//			System.out.println(cnt);
		//ASVariation[][] classes2= ASAnalyzer.determineVariations(g, "isNotAtAllCoding", true);
//		ASAnalyzer.diff(classes,classes2);
//		classes= ASAnalyzer.getVariation("(1^ // 2^)", classes);
//		classes= ASAnalyzer.getVariation("( // 1^2=)", classes);
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; j < 10; j++) {
//				classes[i][j].outputDetail(System.out);
//			}
//		}
		
		//ASAnalyzer.determineVariations(g, "isPartiallyCoding");
		
		//ASAnalyzer.analyzeGene(g, "")
		//ASAnalyzer.lowRepresented(g, 10);
		//ASAnalyzer.debug(g, 10);
		//ASAnalyzer.lengthPlot(g);
		//ASAnalyzer.getVariation("( // 1=2^3=4^5=6^7=8^9=10^11=12^13=14^15=16^17=18^19=20^21=22^23=24^25=26^27=28^)",
		//		g);

	}
	
	public EncodeWrapper(String absFName) {
		super(absFName);		
	}
	
	public EncodeWrapper(String absFName, String speName) {
		this(absFName);
		this.speName= speName;
	}
	
	public EncodeWrapper(String absFName, String speName, String buildVersion) {
		this(absFName, speName);
		this.buildVersion= buildVersion;
	}
	
	public EncodeWrapper(InputStream i) {
		this.inputStream= i;
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
	
	String getBuildVersion() {
		if (buildVersion== null) {
			String[] tokens= fName.split("_");
			for (int i = 0; buildVersion== null&& i < tokens.length; i++) 
				for (int j = 0; buildVersion== null&& j < Species.ASSEMBLY_PFX.length; j++) 
					if (tokens[i].startsWith(Species.ASSEMBLY_PFX[j])) 
						buildVersion= tokens[i];
			if (buildVersion== null) {		
				System.err.println("Coudn't guess build version.");
				buildVersion= "Unknown";
			}
		}
		return buildVersion;
	}

	String getSpeciesName() {
		if (speName== null) {
			String[] tokens= fName.split("_");
			for (int i = 0; speName== null&& i < tokens.length; i++) {
				for (int j = 0; speName== null&& j < Species.SP_NAMES_COMMON.length; j++) {
					if (Species.SP_NAMES_COMMON[j].equalsIgnoreCase(tokens[i]))
						speName= Species.SP_NAMES_COMMON[j];
					else if (Species.getAbbrevNameForBinomial(Species.SP_NAMES_BINOMIAL[j]).equalsIgnoreCase(tokens[i]))
						speName= Species.getAbbrevNameForPrefix(Species.SP_NAMES_BINOMIAL[j]);
					else if (Species.SP_NAMES_BINOMIAL[j].equalsIgnoreCase(tokens[i]))
						speName= Species.SP_NAMES_BINOMIAL[j];
				}
			}
			if (speName== null) {
				System.err.println("Coudn't guess species name");
				speName= "Unknown";
			}
		}
		return speName;
	}
	
	Graph assemble(boolean encode) {
		
		Species spec;
		if (speName== null|| genomeVer== null) {
			spec= new Species("human");
			spec.setGenomeVersion("hg18");
		} else { 
			spec= new Species(speName);
			spec.setGenomeVersion(genomeVer);
		}
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
				if (!ff.isExon()&& !ff.isCDS())
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
						if (encode&& f.isCoding()&& false) {		// in not-encode no exon ids..
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
	 * doesnt take transcript_id order on the input gtfObs into account
	 * problem f.i. with C.elegans -> same tID on different chromosomes
	 * 
	 * @param encode
	 * @return
	 */
	public static Graph assemble(boolean encode, GTFObject[] gtfObs) {
		
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
		
		return (Transcript[][]) Arrays.toField(clusters);
	}

	Graph assemble_stable(boolean encode) {
		
		Species spec= new Species("human");
		spec.setAnnotationVersion("17");
		
			// cluster
		HashMap gHash= getGroups("gene_id", getGtfObj());	// cluster for genes
		
			// infer objects
		Collection co= ((Collection) gHash.keySet());
		String[] keys= new String[co.size()];
		Iterator iter= co.iterator();
		int x= 0;
		while(iter.hasNext()) 
			keys[x++]= (String) iter.next();
		
		for (int i = 0; i < keys.length; i++) {	// genes
			String key= keys[i];
			Gene gene= new Gene(spec, key);
			GTFObject[] gtfs= (GTFObject[]) Arrays.toField(gHash.remove(key));	// transcripts in gene
			France ff= (France) gtfs[0];
			gene.setStrand(ff.getStrand());
			gene.setChromosome(ff.getChromosome());
			gene.setConfidence(ff.getSource());
			HashMap tHash= getGroups("transcript_id", gtfs);
			Collection co2= ((Collection) tHash.keySet());
			String[] tkeys= new String[co2.size()];
			Iterator iter2= co2.iterator();
			x= 0;
			while (iter2.hasNext())
				tkeys[x++]= (String) iter2.next(); 
			for (int j = 0; j < tkeys.length; j++) {	// transcripts
				String tID= tkeys[j];
				Transcript transcript= new Transcript(gene, tID);
				Vector v= (Vector) tHash.remove(tID);
				for (int k = 0; k < v.size(); k++) {		// exons 
					France f= (France) v.elementAt(k);
					if (f.isExon()) 
						transcript.addExon(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
					else if (f.isCDS())
						transcript.addCDS(f.getStart(), f.getEnd());
				}
				gene.addTranscript(transcript);
			}
			if (!encode|| gene.getConfidence().contains("VEGA"))
				gHash.put(key, gene);	// re-fill hash
		}
		
			// build graph
		iter= gHash.values().iterator();
		Graph g= new Graph();
		g.addSpecies(spec);
		while (iter.hasNext()) 
			g.addGene((Gene) iter.next());
		return g;
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
	/* (non-Javadoc)
		 * @see gphase.io.IOWrapper#read()
		 */
		public void read(String[] ids) throws Exception {
			
			if (ids== null|| ids.length< 1)
				return;
			
			BufferedReader buffy;
			File file= new File(fPath+ File.separator+ fName);
			if (fPath!= null&& fName!= null)
				buffy= new BufferedReader(new FileReader(file));
			else 
				buffy= new BufferedReader(new InputStreamReader(inputStream));
			String line;
			int lineCtr= 0;
			Vector gtfVec= new Vector();
			int lastPerc= 0;
			long bytesRead= 0;
			while (buffy.ready()) {
				lineCtr++;
				line= buffy.readLine();
				bytesRead+= (line.length()+ 1);		// else
				int perc= (int) ((bytesRead* 10d)/ file.length());
				if (perc> lastPerc) {
					System.out.print("*");
					System.out.flush();
					++lastPerc;
				}
				String[] cols= line.split("\t");
				if (lineCtr== 1) {
					if (cols.length< 9&& !cols[8].startsWith(GTFObject.GENE_ID_TAG))
						System.err.println("gene_id not in col 8");
				}
				cols= cols[8].split(" ");
				int x= 0;
				for (x = 0; x < ids.length; x++) 
					if (cols[1].substring(1).startsWith(ids[x])||
							cols[1].startsWith(ids[x]))
						break;
				if (x== ids.length)
					continue;
				StringTokenizer toki= new StringTokenizer(line, " \t");	// must be tab, see specification
				if (toki== null) {
					System.err.println("Invalid GTF format: line "+ lineCtr+": "+line);
					continue;
				}
				if (toki.countTokens()< 8)
					System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
				// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
				GTFObject newObj= createGTFObject();
				try {				
					newObj.setSeqname(toki.nextToken());
					newObj.source= toki.nextToken();
					newObj.setFeature(toki.nextToken());
					newObj.start= Integer.parseInt(toki.nextToken());
					newObj.end= Integer.parseInt(toki.nextToken());
					newObj.setScore(toki.nextToken());
					newObj.setStrand(toki.nextToken());
					newObj.setFrame(toki.nextToken());
				} catch (Exception e) {
					System.err.println("Invalid GTF format: line "+ lineCtr+" -- "+e.getMessage());
					//e.printStackTrace();
					continue;
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
							newObj.addAttribute(id, val);
						}
					}
				}
					
				gtfVec.add(newObj);
				//System.out.println(gtfVec.size());
			}
			System.out.println();
			
			gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
		}

		/* (non-Javadoc)
			 * @see gphase.io.IOWrapper#read()
			 */
			public void read() throws Exception {
				
				
				BufferedReader buffy;
				if (fPath!= null&& fName!= null)
					buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
				else 
					buffy= new BufferedReader(new InputStreamReader(inputStream));
				String line;
				int lineCtr= 0;
				Vector gtfVec= new Vector();
				while (buffy.ready()) {
					lineCtr++;
					line= buffy.readLine();
					StringTokenizer toki= new StringTokenizer(line, " \t");	// must be tab, see specification
					if (toki.countTokens()< 8)
						System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
					// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
					GTFObject newObj= createGTFObject();
					try {				
						newObj.seqname= toki.nextToken();
						newObj.source= toki.nextToken();
						newObj.setFeature(toki.nextToken());
						newObj.start= Integer.parseInt(toki.nextToken());
						newObj.end= Integer.parseInt(toki.nextToken());
						newObj.setScore(toki.nextToken());
						newObj.setStrand(toki.nextToken());
						newObj.setFrame(toki.nextToken());
					} catch (Exception e) {
						if (!silent)
							System.err.println("Invalid GTF format: line "+ lineCtr);
						//e.printStackTrace();
						continue;
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
								newObj.addAttribute(id, val);
							}
						}
					}
						
					gtfVec.add(newObj);
					//System.out.println(gtfVec.size());
				}
				
				gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
			}

			public String getGenomeVer() {
				return genomeVer;
			}

			public void setGenomeVer(String genomeVer) {
				this.genomeVer = genomeVer;
			}

			public String getSpeName() {
				return speName;
			}

			public void setSpeName(String speName) {
				this.speName = speName;
			}

			public boolean isSilent() {
				return silent;
			}

			public void setSilent(boolean silent) {
				this.silent = silent;
			}
}
