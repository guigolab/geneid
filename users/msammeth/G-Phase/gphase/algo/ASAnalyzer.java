/*
 * Created on Mar 1, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.algo;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Calendar;
import java.util.Collection;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JWindow;
import javax.swing.plaf.multi.MultiViewportUI;

import com.p6spy.engine.logging.appender.StdoutLogger;
import com.sun.org.apache.xalan.internal.xsltc.compiler.util.FilterGenerator;
import com.sun.org.apache.xml.internal.utils.StringToStringTable;

import gphase.Toolbox;
import gphase.db.EnsemblDBAdaptor;
import gphase.gui.Paparazzi;
import gphase.gui.SpliceOSigner;
import gphase.io.gtf.ATDWrapper;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASEvent;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Region;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.model.SpliceSite.PositionComparator;
import gphase.regex.Automaton;
import gphase.regex.ChessMaster;
import gphase.regex.RegExp;
import gphase.tools.Arrays;
import gphase.tools.ENCODE;
import gphase.tools.Time;

/**
 * 
 * 
 * @author msammeth
 */
public class ASAnalyzer {

	Graph graph;
	
	/**
	 * Generates a (highly redundant) file with all 5UTR exons/exon parts.
	 * Iterated by transcripts.
	 */
	
	public static void output5UTRExonFragments() {
		
		GTFWrapper gtf= null;
		String fName= "5UTRexons.gff";
		fName= Toolbox.checkFileExists(fName);
		gtf= new GTFWrapper(fName);
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		Comparator compi= new DirectedRegion.OrderComparator();
		Vector regV= new Vector();
		HashMap attMap= new HashMap();
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				Exon[] e= t[j].getExons();
				java.util.Arrays.sort(e, compi);
				for (int k = 0; k < e.length; k++) {
					String[] att= new String[] {t[j].getTranscriptID(), e[k].getExonID(), Gene.REGION_ID[Gene.REGION_TRANSCRIPT_5UTR]};
					if (e[k].isCoding()) {
						if (!e[k].isCoding5Prime()) {
							int p1= Math.abs(e[k].get5PrimeEdge());
							int p2= Math.abs(e[k].get5PrimeCDS()- 1);
							if (p1> p2) {
								int h= p1;
								p1= p2;
								p2= h;
							}
							DirectedRegion reg= new DirectedRegion(p1, p2, e[k].getGene().getStrand());
							reg.setChromosome(e[k].getChromosome());
							reg.setStrand(e[k].getGene().getStrand());
							reg.setID("part_exon");
							regV.add(reg);
							attMap.put(reg, att);
						}
						break;
					}
					regV.add(e[k]);
					e[k].setStrand(e[k].getGene().getStrand());
					attMap.put(e[k], att);
				}
			}
		}
		
		Vector v= new Vector(regV.size());		
		for (int i = 0; i < regV.size(); i++) 
			v.add(GTFObject.createGFFObject((DirectedRegion) regV.elementAt(i), "gencode", (String[]) attMap.get(regV.elementAt(i))));
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v)); 
		try {
			gtf.writeGFF();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Generates a (highly redundant) file with all 5UTR exons/exon parts.
	 * Iterated by transcripts.
	 */
	
	public static void output5UTRIntrons() {
		
		GTFWrapper gtf= null;
		String fName= "5UTRintrons.gff";
		fName= Toolbox.checkFileExists(fName);
		gtf= new GTFWrapper(fName);
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		Comparator compi= new DirectedRegion.OrderComparator();
		Vector regV= new Vector();
		HashMap attMap= new HashMap();
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				Exon[] e= t[j].getExons();
				java.util.Arrays.sort(e, compi);
				for (int k = 1; k < e.length; k++) {	// 1st exon doesnt have an intron before..
					String[] att= new String[] {t[j].getTranscriptID(), "before "+ e[k].getExonID(), Gene.REGION_ID[Gene.REGION_TRANSCRIPT_5UTR]};
					int p1= Math.abs(e[k-1].get3PrimeEdge()+ 1);
					int p2= Math.abs(e[k].get5PrimeEdge()- 1);
					if (p1> p2) {
						int h= p1;
						p1= p2;
						p2= h;
					}					
					DirectedRegion reg= new DirectedRegion(p1, p2, e[k].getStrand());
					reg.setChromosome(e[k].getChromosome());
					reg.setStrand(e[k].getGene().getStrand());
					reg.setID("intron");
					regV.add(reg);
					attMap.put(reg, att);
					if (e[k].isCoding())
						break;
				}
			}
		}
		
		Vector v= new Vector(regV.size());		
		for (int i = 0; i < regV.size(); i++) 
			v.add(GTFObject.createGFFObject((DirectedRegion) regV.elementAt(i), "gencode", (String[]) attMap.get(regV.elementAt(i))));
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v)); 
		try {
			gtf.writeGFF();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Generates a (highly redundant) file with all CDS introns.
	 * Iterated by transcripts.
	 */
	
	public static void outputCDSIntrons() {
		
		GTFWrapper gtf= null;
		String fName= "CDSintrons.gff";
		fName= Toolbox.checkFileExists(fName);
		gtf= new GTFWrapper(fName);
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		Comparator compi= new DirectedRegion.OrderComparator();
		Vector regV= new Vector();
		HashMap attMap= new HashMap();
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				Exon[] e= t[j].getExons();
				java.util.Arrays.sort(e, compi);
				boolean start= false;
				for (int k = 0; k < e.length; k++) {
					if (!start) {
						if (e[k].isCoding()) 
							start= true;
						continue;
					}
					String[] att= new String[] {t[j].getTranscriptID(), "before "+ e[k].getExonID(), Gene.REGION_ID[Gene.REGION_TRANSCRIPT_CDS]};
					int p1= Math.abs(e[k-1].get3PrimeEdge()+ 1);
					int p2= Math.abs(e[k].get5PrimeEdge()- 1);
					if (p1> p2) {
						int h= p1;
						p1= p2;
						p2= h;
					}					
					DirectedRegion reg= new DirectedRegion(p1, p2, e[k].getStrand());
					reg.setChromosome(e[k].getChromosome());
					reg.setStrand(e[k].getGene().getStrand());
					reg.setID("intron");
					regV.add(reg);
					attMap.put(reg, att);

					if (start&& !e[k].isCoding3Prime())
						break;
					
				}
			}
		}
		
		Vector v= new Vector(regV.size());		
		for (int i = 0; i < regV.size(); i++) 
			v.add(GTFObject.createGFFObject((DirectedRegion) regV.elementAt(i), "gencode", (String[]) attMap.get(regV.elementAt(i))));
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v)); 
		try {
			gtf.writeGFF();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	/**
	 * Generates a (highly redundant) file with all 5UTR exons/exon parts.
	 * Iterated by transcripts.
	 */
	
	public static void outputSpliceSites() {
		
		GTFWrapper gtf= null;
		String fName= "splicesites_all_cod_transcripts.gtf";
		Toolbox.checkFileExists(fName);
		gtf= new GTFWrapper(fName);
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		g.getASVariations(ASMultiVariation.FILTER_NONE);	// for init const/alt
		Comparator compi= new DirectedRegion.OrderComparator();
		Vector v= new Vector();
		String[] attMap= null;
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				SpliceSite[] ss= t[j].getSpliceChain();
				for (int k = 0; k < ss.length; k++) {
					if (ss[k].getPos()< t[j].getTranslations()[0].get5PrimeEdge()) 
						attMap= new String[]{t[j].getTranscriptID(), "5UTR_transcript"};
					else if (ss[k].getPos()> t[j].getTranslations()[0].get3PrimeEdge()) 
						attMap= new String[]{t[j].getTranscriptID(), "3UTR_transcript"};
					else
						attMap= new String[]{t[j].getTranscriptID(), "CDS_transcript"};
					v.add(GTFObject.createGFFObject(ss[k], "gencode", attMap));
				}
			}
		}
		
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v)); 
		try {
			gtf.writeGFF();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void output5UTRLengthAnalysis(Graph g, PrintStream pr) {
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
		for (int i = 0; i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				if (vars[i][j].isCompletelyIn5UTR()) {
					int len1before= vars[i][j].getTranscript1().get5UTR(true);
					int len1after= vars[i][j].getTranscript1().get5UTR(false);
					int len2before= vars[i][j].getTranscript2().get5UTR(true);
					int len2after= vars[i][j].getTranscript2().get5UTR(false);
					pr.println(
						vars[i][j].getTranscript1().getTranscriptID()+ " "+
						len1before+ " "+
						len1after+ " "+
						(len1before- len1after)+ "\t"+
						
						vars[i][j].getTranscript2().getTranscriptID()+ " "+
						len2before+ " "+
						len2after+ " "+
						(len2before- len2after)+"\t"+
						
						(len1before- len2before)+" "+	// TSS
						(len1after- len2after)+"\t"+	//TSS+ AS
						
						(vars[i][j].getLengthDiff_relative12(true))+"\t"+	// AS
						vars[i][j]);
				}
			}
		}
	}
	
	public static void outputConstitutiveExons(Graph g) {
		Exon[] ex= g.getExons();
		for (int i = 0; i < ex.length; i++) {
			if (ex[i].getTranscripts().length== 1)
				System.out.println(ex[i].toPosString());
		}
	}
	
	public static long[] getCDSUTRSizeComplement(Graph g) {
		Species[] spe= g.getSpecies();
		long[] regions= new long[3];
		for (int i = 0; i < regions.length; i++) 
			regions[i]= 0;
		for (int i = 0; i < spe.length; i++) {
			Gene[] ge= spe[i].getGenes();
			for (int j = 0; j < ge.length; j++) {
				try {
					regions[0]+= Math.abs(ge[j].getReal5UTR().getEnd()- ge[j].getReal5UTR().getStart())+ 1;
				} catch (NullPointerException e) {;}
				try {
					regions[1]+= Math.abs(ge[j].getCDSRegion().getEnd()- ge[j].getCDSRegion().getStart())+ 1;
				} catch (NullPointerException e) {;}
				try {
					regions[2]+= Math.abs(ge[j].getReal3UTR().getEnd()- ge[j].getReal3UTR().getStart())+ 1;
				} catch (NullPointerException e) {;}
			}
		}
		
		return regions;
	}
	
	public static void check_AA_AD_old(Graph g, boolean pure, boolean filterLow, boolean filterFirstExon) {
		g.filterNonCodingTranscripts();
		ASVariation[][] as= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		
			// get all events containing AA, AD
		Vector adVec= new Vector();
		Vector aaVec= new Vector();
		for (int i = 0; i < as.length; i++) {
			String s= as[i][0].toBitString();
			int x= 0;
			for (int j = 0; j < as[i].length; j++) {
				Vector va= new Vector();
				Vector vd= new Vector();
				//if ((x= s.indexOf("DBA"))>= 0|| (x= s.indexOf("DCA"))>= 0) {
				//if (s.startsWith("DBA")|| s.startsWith("DCA")) {
				if (s.equals("DBA")|| s.equals("DCA")) {
					SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
					//if (Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10)	// filter short
						vd.add(new ASVariation(
							as[i][j].getTranscript1(),
							as[i][j].getTranscript2(),
							ss[0],
							ss[1]
						));
				}
				//if ((x= s.indexOf("ABD"))>= 0|| (x= s.indexOf("ACD"))>= 0) {
				if (s.equals("ABD")|| s.equals("ACD")) {
				//if (s.startsWith("ABD")|| s.startsWith("ACD")) {
					SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
					//if (Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10)	// filter short
						va.add(new ASVariation(
							as[i][j].getTranscript1(),
							as[i][j].getTranscript2(),
							ss[0],
							ss[1]
						));
				}
				if (va.size()> 0)
					aaVec.add(va);
				if (vd.size()> 0)
					adVec.add(vd);
			}
		}
		ASVariation[][] aa= (ASVariation[][]) Arrays.toField(aaVec);
		ASVariation[][] ad= (ASVariation[][]) Arrays.toField(adVec);

			// test out
//		System.out.println("-- donors:");
//		for (int i = 0; i < ad.length; i++) 
//			System.out.println(ad[i][0]);		
//		System.out.println("\n-- acceptors:");
//		for (int i = 0; i < aa.length; i++) 
//			System.out.println(aa[i][0]);

			// count coding/non-coding
		int cds= 0, utr= 0;
		for (int i = 0; i < aa.length; i++) 
			for (int j = 0; j < aa[i].length; j++) 
				if (aa[i][j].isProteinCoding())
					aa[i][j].outputDetail(System.out); //++cds;
				else if (aa[i][j].isNotAtAllCoding())
					aa[i][j].outputDetail(System.out); //++utr;
			
		System.out.println("acceptors: "+cds+" cds, "+utr+"utr");
		cds= 0; utr= 0;
		for (int i = 0; i < ad.length; i++) 
			for (int j = 0; j < ad[i].length; j++) {
				if (ad[i][j].isProteinCoding()) 
					ad[i][j].outputDetail(System.out); //++cds;
				else if (ad[i][j].isNotAtAllCoding()) {
					ad[i][j].outputDetail(System.out); //++utr;
					AbstractRegion rcds= ad[i][j].getGene().getRealCDS();
					SpliceSite[] su= ad[i][j].getSpliceUniverse();
					// TODO stopped here
				}
			}
		System.out.println("donors: "+cds+" cds, "+utr+"utr");
		
	}


	public static ASVariation[][] get_AA_AD(Graph g, boolean pure, boolean filterLow, boolean filterFirstExon) {
		
		//g.filterNonCodingTranscripts();
		ASVariation[][] as= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		
			// get all events containing AA, AD
		Vector adVec= new Vector();
		Vector aaVec= new Vector();
		for (int i = 0; i < as.length; i++) {
			String s= as[i][0].toBitString();
			int x= 0;
			for (int j = 0; j < as[i].length; j++) {
				Vector va= new Vector();
				Vector vd= new Vector();
				x= 0;
				if ((pure&& ((s.equals("ABD")|| s.equals("ACD"))))||
						(!pure&& ((x= s.indexOf("ABD"))>= 0|| (x= s.indexOf("ACD"))>= 0))) {
					//if (s.startsWith("ABD")|| s.startsWith("ACD")) {
						SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
						SpliceSite[] flanks= as[i][j].getFlankingSpliceSites();		// skip terminal exon
						if (((filterFirstExon&& flanks[1]!= null)|| !filterFirstExon)&&
								(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10))	// filter short
							va.add(new ASVariation(
								as[i][j].getTranscript1(),
								as[i][j].getTranscript2(),
								ss[0],
								ss[1]
							));
					}
				
				
				x= 0;
				if ((pure&& ((s.equals("DBA")|| s.equals("DCA"))))||
					(!pure&& ((x= s.indexOf("DBA"))>= 0|| (x= s.indexOf("DCA"))>= 0))) {
				//if (s.startsWith("DBA")|| s.startsWith("DCA")) {
					SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
					SpliceSite[] flanks= as[i][j].getFlankingSpliceSites();		// skip initial exons
					if (((filterFirstExon&& flanks[0]!= null)|| !filterFirstExon)&&
							(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10)) {	// filter short
						vd.add(new ASVariation(
							as[i][j].getTranscript1(),
							as[i][j].getTranscript2(),
							ss[0],
							ss[1]
						));
//					} else {
//						if (((filterFirstExon&& flanks[0]!= null)|| !filterFirstExon)&&
//								(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())< 10)) {	// filter short
//							if (as[i][j].isNotAtAllCoding())
//								as[i][j].outputDetail(System.out);
//						}
					}
				}
				if (va.size()> 0)
					aaVec.add(va);
				if (vd.size()> 0)
					adVec.add(vd);
			}
		}
		ASVariation[][] aa= (ASVariation[][]) Arrays.toField(aaVec);
		ASVariation[][] ad= (ASVariation[][]) Arrays.toField(adVec);
		
		return (ASVariation[][]) Arrays.merge(aa, ad).toArray();		
	}
		
	public static void check_AA_AD(Graph g, boolean pure, boolean filterLow, boolean filterFirstIntron) {
		
		g.filterNonCodingTranscripts();
		ASVariation[][] as= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		
			// get all events containing AA, AD
		Vector adVec= new Vector();
		Vector aaVec= new Vector();
		for (int i = 0; i < as.length; i++) {
			String s= as[i][0].toBitString();
			int x= 0;
			for (int j = 0; j < as[i].length; j++) {
				Vector va= new Vector();
				Vector vd= new Vector();
				x= 0;
				if ((pure&& ((s.equals("ABD")|| s.equals("ACD"))))||
						(!pure&& ((x= s.indexOf("ABD"))>= 0|| (x= s.indexOf("ACD"))>= 0))) {
					//if (s.startsWith("ABD")|| s.startsWith("ACD")) {
						SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
						SpliceSite[] flanks= as[i][j].getFlankingSpliceSites();		// skip terminal exon
//						if (((filterFirstExon&& flanks[1]!= null)|| !filterFirstExon)&&
//								(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10))	// filter short
						SpliceSite s1= as[i][j].getTranscript1().getSpliceChain()[0];	// skip first intron
						if (s1.getPos()> as[i][j].getTranscript2().getSpliceChain()[0].getPos())
							s1= as[i][j].getTranscript1().getSpliceChain()[0];
						PositionComparator compi= new SpliceSite.PositionComparator();
						if (((filterFirstIntron&& flanks[0]!= null&& compi.compare(flanks[0],s1)!= 0)|| !filterFirstIntron)&&
								(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10))	// filter short
							va.add(new ASVariation(
								as[i][j].getTranscript1(),
								as[i][j].getTranscript2(),
								ss[0],
								ss[1]
							));
					}
				
				
				x= 0;
				if ((pure&& ((s.equals("DBA")|| s.equals("DCA"))))||
					(!pure&& ((x= s.indexOf("DBA"))>= 0|| (x= s.indexOf("DCA"))>= 0))) {
				//if (s.startsWith("DBA")|| s.startsWith("DCA")) {
					SpliceSite[][] ss= as[i][j].getSpliceSites(x, x+2);
					SpliceSite[] flanks= as[i][j].getFlankingSpliceSites();		// skip initial exons
					if (((filterFirstIntron&& flanks[0]!= null)|| !filterFirstIntron)&&
							(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())>= 10)) {	// filter short
						vd.add(new ASVariation(
							as[i][j].getTranscript1(),
							as[i][j].getTranscript2(),
							ss[0],
							ss[1]
						));
//					} else {
//						if (((filterFirstExon&& flanks[0]!= null)|| !filterFirstExon)&&
//								(!filterLow|| Math.abs(ss[0][0].getPos()- ss[1][0].getPos())< 10)) {	// filter short
//							if (as[i][j].isNotAtAllCoding())
//								as[i][j].outputDetail(System.out);
//						}
					}
				}
				if (va.size()> 0)
					aaVec.add(va);
				if (vd.size()> 0)
					adVec.add(vd);
			}
		}
		ASVariation[][] aa= (ASVariation[][]) Arrays.toField(aaVec);
		ASVariation[][] ad= (ASVariation[][]) Arrays.toField(adVec);
		
		
			// test out
//		System.out.println("-- donors:");
//		for (int i = 0; i < ad.length; i++) 
//			System.out.println(ad[i][0]);		
//		System.out.println("\n-- acceptors:");
//		for (int i = 0; i < aa.length; i++) 
//			System.out.println(aa[i][0]);

		
		System.out.print("\nAA/AD");
		if (pure)
			System.out.print(" pure");
		else 
			System.out.print(" all");
		if (filterLow)
			System.out.print(" filter <10 nt");
		if (filterFirstIntron)
			System.out.print(" filter 1st intron");
		System.out.println(" (cds, utr)");
		
			// count coding/non-coding
		int cds= 0, utr= 0;
		for (int i = 0; i < aa.length; i++) 
			for (int j = 0; j < aa[i].length; j++) 
				if (aa[i][j].isProteinCoding())
					++cds;	//aa[i][j].outputDetail(System.out); // aa[i][j].outputDetail(System.out); //
				else if (aa[i][j].isNotAtAllCoding())
					++utr; //aa[i][j].outputDetail(System.out); //
			
		System.out.println("AA: "+cds+" cds,\t"+utr+"utr");
		
		cds= 0; utr= 0;
		for (int i = 0; i < ad.length; i++) 
			for (int j = 0; j < ad[i].length; j++) {
				if (ad[i][j].isProteinCoding()) 
					++cds;	//ad[i][j].outputDetail(System.out); //
				else if (ad[i][j].isNotAtAllCoding()) {
					++utr;	// ad[i][j].outputDetail(System.out); //
					AbstractRegion rcds= ad[i][j].getGene().getRealCDS();
					SpliceSite[] su= ad[i][j].getSpliceUniverse();
					// TODO stopped here
				}
			}
		System.out.println("AD: "+cds+" cds,\t"+utr+"utr");
		

	}
	
	public ASAnalyzer(Graph newGraph) {
		this.graph= newGraph;
	}

	public static void test01_clusters_coverage_as() {
		
		String fName= "test01_clusters_coverage_encode_cds.txt";			
		PrintStream p= null;
		try {
			fName= Toolbox.checkFileExists(fName);
			p= new PrintStream(fName);
			//p= System.out;
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		//g.filterNMDTranscripts();
		
		p.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
		ASAnalyzer.filterSingleExonClusters(g);
		int x= g.getGenes().length;
		int y= g.getTranscriptCount();
	
			// median
		Gene[] ge= g.getGenes();
		Vector distrV= new Vector();
		for (int i = 0; i < ge.length; i++) {
			Integer a= new Integer(ge[i].getTranscriptCount());
			distrV.add(a);
		}
		Integer[] distr= (Integer[]) Arrays.toField(distrV);
		int[] dist= new int[distr.length];
		for (int i = 0; i < dist.length; i++) 
			dist[i]= distr[i].intValue();
		java.util.Arrays.sort(dist);
		float med;
		if (dist.length% 2== 0)
			med=  dist[dist.length/2];
		else
			med= (dist[dist.length/2]+ dist[dist.length/2+1])/ 2f;
		p.println("clusters: "+x+" transcripts: "+y+"\t"+ ((float) y/ (float) x)+"\t"+ med);
		
			// median
		ge= g.getClustersWithAS();
		if (ge!= null&& ge.length> 0) {
			distrV= new Vector();
			for (int i = 0; i < ge.length; i++) {
				Integer a= new Integer(ge[i].getTranscriptCount());
				distrV.add(a);
			}
			distr= (Integer[]) Arrays.toField(distrV);
			dist= new int[distr.length];
			for (int i = 0; i < dist.length; i++) 
				dist[i]= distr[i].intValue();
			java.util.Arrays.sort(dist);
			if (dist.length% 2== 0)
				med=  dist[dist.length/2];
			else
				med= (dist[dist.length/2]+ dist[dist.length/2+1])/ 2f;
			int z= ge.length;
			p.println("clusters w as: "+z+"\t"+((float) z/ (float) x)+ "\t"+ med);
		} else
			p.println("clusters w as: "+0+"\t"+0+ "\t"+ 0);
	}

	public static void test01_clusters_coverage_as(Graph g, PrintStream p) {
		p.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
		ASAnalyzer.filterSingleExonClusters(g);
		int x= g.getGenes().length;
		int y= g.getTranscriptCount();

			// median
		Gene[] ge= g.getGenes();
		Vector distrV= new Vector();
		for (int i = 0; i < ge.length; i++) {
			Integer a= new Integer(ge[i].getTranscriptCount());
			distrV.add(a);
		}
		Integer[] distr= (Integer[]) Arrays.toField(distrV);
		int[] dist= new int[distr.length];
		for (int i = 0; i < dist.length; i++) 
			dist[i]= distr[i].intValue();
		java.util.Arrays.sort(dist);
		float med;
		if (dist.length% 2== 0)
			med=  dist[dist.length/2];
		else
			med= (dist[dist.length/2]+ dist[dist.length/2+1])/ 2f;
		p.println("clusters: "+x+" transcripts: "+y+"\t"+ ((float) y/ (float) x)+"\t"+ med);
		
			// median
		ge= g.getClustersWithAS();
		if (ge!= null&& ge.length> 0) {
			distrV= new Vector();
			for (int i = 0; i < ge.length; i++) {
				Integer a= new Integer(ge[i].getTranscriptCount());
				distrV.add(a);
			}
			distr= (Integer[]) Arrays.toField(distrV);
			dist= new int[distr.length];
			for (int i = 0; i < dist.length; i++) 
				dist[i]= distr[i].intValue();
			java.util.Arrays.sort(dist);
			if (dist.length% 2== 0)
				med=  dist[dist.length/2];
			else
				med= (dist[dist.length/2]+ dist[dist.length/2+1])/ 2f;
			int z= ge.length;
			p.println("clusters w as: "+z+"\t"+((float) z/ (float) x)+ "\t"+ med);
		} else
			p.println("clusters w as: "+0+"\t"+0+ "\t"+ 0);
	}

	public static void test01b_clusters_coverage_as(Graph g, PrintStream p) {
		filterSingleExonClusters(g);
		Gene[] ge= g.getClustersWithAS();
		System.out.println(g.getGenes().length+ " clusters, "+ ge.length+ " with AS ("+((float) ge.length/ g.getGenes().length)+")");
		int ctr= 0;
		for (int i = 0; i < ge.length; i++) {
			System.out.println(ge[i].getTranscripts().length);
			ctr+= ge[i].getTranscripts().length;
		}
		System.out.println(ctr+ " transcripts over "+ ge.length+ " clusters ("+((float) ctr/ ge.length)+")");
	}

	static int getOutfileNumber() {
		File dir= new File("./");
		File[] files= dir.listFiles();
		int max= 0;
		for (int i = 0; i < files.length; i++) {
			String s= files[i].getName();
			if (!s.startsWith("out"))
				continue;
			s= s.substring(3, s.indexOf("_"));
			int x= 0;
			try {
				x= Integer.parseInt(s);
			} catch (NumberFormatException e) {
				;// :)
			}
			max= Math.max(x, max);
		}
		return max+1;
	}
	
	public static void test(String[] methodName, Class[] args, boolean[] concat) {
		
		System.out.println(new Date(System.currentTimeMillis()));
		try {
			String[] outF= new String[methodName.length];
			Method[] m= new Method[methodName.length];
			int base= getOutfileNumber();
			for (int i = 0; i < outF.length; i++) {
				outF[i]= ""+ base++;
				for (int j = 0; outF[i].length() < 3; j++) 
					outF[i]= "0"+ outF[i];
				outF[i]= "out"+ outF[i]+ "_"+ methodName[i];
				
				m[i]= null;
				try {
					m[i] = ASAnalyzer.class.getMethod(methodName[i], args);
				} catch (Exception e) {
					e.printStackTrace();
				} 
				
			}
		
			PrintStream pr;
			EncodeWrapper myEncode= new EncodeWrapper(new File("encode/44regions_genes_CHR_coord.gtf").getAbsolutePath());
			Graph g= myEncode.getGraph(true);
			for (int i = 0; i < m.length; i++) {
				if (concat[i]) 
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_summary", true)));
				else 
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_encode_total")));
				pr.println("ENCODE total");
				m[i].invoke(null, new Object[]{g,pr});
				pr.flush(); pr.close();
				System.out.println("ENCODE total finished with analysis "+m[i]);
			}
			
//				// rebuild graph: otherwise splice sites are not removed, test02 fails..
//			myEncode= new EncodeWrapper(new File("encode/44regions_genes_CHR_coord.gtf").getAbsolutePath());
//			g= myEncode.getGraph(true);
//			g.filterNonCodingTranscripts();
//			for (int i = 0; i < m.length; i++) {
//				if (concat[i]) {
//					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_summary", true)));
//					pr.println("\n\nENCODE coding");
//				} else
//					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_encode_coding")));
//				m[i].invoke(null, new Object[]{g,pr}); 	
//				pr.flush(); pr.close();
//				System.out.println("ENCODE coding finished with analysis "+m[i]);
//			}
//			g= null;
//			
			
			myEncode= new EncodeWrapper(new File("encode/EnsemblGenes_fromUCSC.gtf").getAbsolutePath());
			g= myEncode.getGraph(false);
			g.filterNonCodingTranscripts();
			for (int i = 0; i < m.length; i++) {
				if (concat[i]) {
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_summary", true)));
					pr.println("\n\nENSEMBL coding");
				} else
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_ensembl_coding")));
				m[i].invoke(null, new Object[]{g,pr});
				pr.flush(); pr.close();
				System.out.println("ENSEMBL coding finished with analysis "+m[i]);
			}
			g= null;
			
			
			myEncode= new EncodeWrapper(new File("encode/RefSeqGenes_fromUCSC.gtf").getAbsolutePath());
			g= myEncode.getGraph(false);
			g.filterNonCodingTranscripts();
			for (int i = 0; i < m.length; i++) {
				if (concat[i]) {
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_summary", true)));
					pr.println("\n\nREFSEQ coding");
				} else
					pr= new PrintStream(new BufferedOutputStream(new FileOutputStream(outF[i]+"_refseq_coding")));
				m[i].invoke(null, new Object[]{g,pr}); 	
				pr.flush(); pr.close();
				System.out.println("REFSEQ coding finished with analysis "+m[i]);
			}

		
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(new Date(System.currentTimeMillis()));
		
	}

	public static void filterSingleExonClusters(Graph g) {
		Species[] spe= g.getSpecies();
		for (int i = 0; i < spe.length; i++) {
			Gene[] ge= spe[i].getGenes();
			for (int j = 0; j < ge.length; j++) {
				if (ge[j].getExons().length== 1)
					g.removeKill(ge[j]);
			}
		}
	}
	
	public static void getSylvainsSize(Graph g, PrintStream p) {
		
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
		
		// skipped
		p.println("--- skipped exon ---");
		ASVariation[] sk= getVariation("( // 1=2^)", vars);
		for (int i = 0; i < sk.length; i++) {
			SpliceSite[] sc= (sk[i].getSpliceChain1()== null|| sk[i].getSpliceChain1().length< 1)?
					sk[i].getSpliceChain2():sk[i].getSpliceChain1();
			p.println(sk[i].getGene().getChromosome()+ " "
					+ sk[i].getGene().getStrand()+ " "
					+ sc[0].getPos()+ " "
					+ sc[1].getPos()+ " "
					+ sk[i].toCodingInfoString());
		}

		// donors
		p.println();
		p.println("--- alt donor ---");
		ASVariation[] don= getVariation("(1^ // 2^)", vars);
		for (int i = 0; i < don.length; i++) {
			SpliceSite[] flanks= don[i].getFlankingSpliceSites();
			if (flanks[0]== null)
				continue;	// skip initial exons
			p.println(don[i].getGene().getChromosome()+ " "
					+ don[i].getGene().getStrand()+ " "
					+ flanks[0].getPos()+ " "
					+ don[i].getSpliceChain1()[0].getPos()+ " "
					+ don[i].toCodingInfoString());
			p.println(don[i].getGene().getChromosome()+ " "
					+ don[i].getGene().getStrand()+ " "
					+ flanks[0].getPos()+ " "
					+ don[i].getSpliceChain2()[0].getPos()+ " "
					+ don[i].toCodingInfoString());
		}

		// acceptors
		p.println();
		p.println("--- alt acceptor ---");
		ASVariation[] acc= getVariation("(1= // 2=)", vars);
		for (int i = 0; i < acc.length; i++) {
			SpliceSite[] flanks= acc[i].getFlankingSpliceSites();
			if (flanks[1]== null)
				continue;	// skip terminal exons
			p.println(acc[i].getGene().getChromosome()+ " "
					+ acc[i].getGene().getStrand()+ " "
					+ acc[i].getSpliceChain1()[0].getPos()+ " "
					+ flanks[1].getPos()+ " "
					+ acc[i].toCodingInfoString());
			p.println(acc[i].getGene().getChromosome()+ " "
					+ acc[i].getGene().getStrand()+ " "
					+ acc[i].getSpliceChain2()[0].getPos()+ " "
					+ flanks[1].getPos()+ " "
					+ acc[i].toCodingInfoString());
		}
}
	
	public static void getAllASVariation(Graph g) {
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		outputVariations(classes, false, false, System.out);
		int cnt= 0;
		try {
			PrintStream printer= new PrintStream("out005_encode_vars_all_non-redundant");
			for (int i = 0; i < classes.length; i++) {
				cnt+= classes[i].length;
				for (int j = 0; j < classes[i].length; j++) {
					classes[i][j].outputDetail(printer);
				}
				printer.print("\n\n");
			}
			printer.flush();
			printer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(cnt);
	}
	public static void getVincenco(Graph g) {
	
		ASVariation[][] classes= g.getASVariations(1);
		ASVariation[] events= getVariation("( // 1^2=)", classes);
		int cnt1= 0, cnt2= 0, cnt3= 0, cnt4= 0;
		Comparator compi= new SpliceSite.PositionComparator();
		for (int i = 0; i < events.length; i++) {
			if ((events[i].getTranscript1().getSpliceChain().length> 0&&
					compi.compare(events[i].getBorderingSpliceSites()[0],
					events[i].getTranscript1().getSpliceChain()[0])== 0)||
				(events[i].getTranscript2().getSpliceChain().length> 0&& 
					compi.compare(events[i].getBorderingSpliceSites()[0],
						events[i].getTranscript2().getSpliceChain()[0])== 0)) {
				cnt1++;
				if (events[i].getLengthDiff(false)< 150) {
					cnt3++;
					if (events[i].isSemiProteinCoding()) {
						cnt4++;
						events[i].outputDetail(System.out);
					}
				}
			}
			if (events[i].getLengthDiff(false)< 150) {
				cnt2++;
			}
		} 
		
		System.out.println(cnt1+","+cnt2+","+cnt3+","+cnt4);
	}
	
	
	public static void lengthVariation(Graph g) {
			ASVariation[][] classes= g.getASVariations(1);
			classes= (ASVariation[][]) filter(classes, "isNotAtAllCoding");
			classes= (AS.ArraysVariation[][]) Arrays.sort2DFieldRev(classes);
			//writeLabels(classes);
			classes= sortLabels(classes);
			
			
			
				// for plots
//			ASVariation[][] vRev= new ASVariation[classes[0].length][];
//			for (int i = 0; i < vRev.length; i++) {
//				vRev[i]= new ASVariation[classes.length];
//				for (int j = 0; j < vRev[i].length; j++) {
//					if (classes[j]== null|| classes[j].length<= i) 
//						vRev[i][j]= null;
//					else
//						vRev[i][j]= classes[j][i];
//				}
//			}
//			classes= vRev;
//			
//			//ASVariation[] events= getVariation("( // 1^2=)", classes);
//			for (int x = 0; x < classes[0].length; x++) 
//				System.out.print(classes[0][x]+ ",");
//			System.out.println();
//			for (int x = 0; x < classes.length; x++) {
//				ASVariation[] events= classes[x];
//				for (int i = 0; i < events.length; i++) {
//					if (events[i]== null) {
//						System.out.print(",");
//						continue;
//					}
//					int[] a= events[i].getLength(true);
//					int diffA= Math.abs(a[0]- a[1]);
//	//				events[i].outputDetail(System.out);
//					System.out.print(diffA+ ",");
//				}
//				System.out.println();
//			}
		}

	public static void test03b_lengthVariationModuloAA_AD(Graph g, PrintStream p) {
			ASVariation[][] classes= g.getASVariations(1);
			classes= (ASVariation[][]) Arrays.filter(classes, "isNotAtAllCoding");
			classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
			
			//ASVariation[] events= getVariation("( // 1^2=)", classes);
			int cnt0= 0, cnt1= 0, cnt2= 0;
			int scnt0= 0, scnt1= 0, scnt2= 0, stot= 0;
			int tot= 0;
			for (int x = 0; classes!= null&& x < classes.length; x++) {
				cnt1= 0; cnt2= 0; cnt0= 0;
				tot= 0;
				ASVariation[] events= classes[x];
				p.println(events[0].toString());
				boolean chkLen= false;
				if (classes[x][0].toString().equals("(1^ // 2^)")||
						classes[x][0].toString().equals("(1= // 2=)"))
					chkLen= true;
				for (int i = 0; i < events.length; i++) {
					int[] a= events[i].getLength(true);
					int[] b= events[i].getLength(false);
					int diffA= Math.abs(a[0]- a[1]);
					int diffB= a[0]- a[1];
					if (chkLen&& diffA< 11)
						continue;
	//				if (x== 0)
	//					events[i].outputDetail(System.out);
	//				System.out.println(a[0]+","+a[1]+": "+diffA+"("+(diffA%3)+")");
					if (diffA%3== 0)
						cnt0++;
					if (diffA%3== 1)
						cnt1++;
					if (diffA%3== 2)
						cnt2++;
					++tot;
				}
				p.println("0: "+cnt0+"\t1: "+cnt1+"\t2: "+cnt2+"\t("+tot+"): "+((double) cnt0/ tot));
				scnt1+= cnt1;
				scnt0+= cnt0;
				scnt2+= cnt2;
				stot+= tot;
			}
			p.println("===> 0: "+scnt0+"\t1: "+scnt1+"\t2: "+scnt2+"\t("+stot+"): "+((double) scnt0/ stot));
		}

	public static void test03b_lengthVariationFirstExon(Graph g, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.filter(classes, "isNotAtAllCoding");
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		
		//ASVariation[] events= getVariation("( // 1^2=)", classes);
		int bal= 0, bal0= 0;
		for (int x = 0; classes!= null&& x < classes.length; x++) {
			if (!classes[x][0].toString().equals("(1^ // 2^)"))
				continue;
			ASVariation[] events= classes[x];
			p.println(events[0].toString());
			for (int i = 0; i < events.length; i++) {
				if (events[i].getTranscript1().getPredSpliceSite(events[i].getSpliceChain1()[0])!= null||
						events[i].getTranscript2().getPredSpliceSite(events[i].getSpliceChain2()[0])!= null)
					continue;
				int[] a= events[i].getLength(true);
				int diffA= a[0]- a[1];
				bal0+= diffA;
				int diffB;
				diffB= events[i].getTranscript2().get5PrimeEdge()- 
							events[i].getTranscript1().get5PrimeEdge();
				
				bal+= diffB;
				System.out.println(i+":\t"+diffA+"\t"+diffB+"\t"+(diffA+ diffB));
			}
		}
		p.println("===> "+bal0+ ", "+bal);
	}
	
	public static void test03_lengthVariationModulo(Graph g, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.filter(classes, "isProteinCoding.Arrays");
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		
		//ASVariation[] events= getVariation("( // 1^2=)", classes);
		int cnt0= 0, cnt1= 0, cnt2= 0;
		int scnt0= 0, scnt1= 0, scnt2= 0, stot= 0;
		int tot= 0;
		for (int x = 0; classes!= null&& x < classes.length; x++) {
			cnt1= 0; cnt2= 0; cnt0= 0;
			tot= 0;
			ASVariation[] events= classes[x];
			p.println(events[0].toString());
			for (int i = 0; i < events.length; i++) {
				int[] a= events[i].getLength(true);
				int[] b= events[i].getLength(false);
				int diffA= Math.abs(a[0]- a[1]);
//				if (x== 0)
//					events[i].outputDetail(System.out);
//				System.out.println(a[0]+","+a[1]+": "+diffA+"("+(diffA%3)+")");
				if (diffA%3== 0)
					cnt0++;
				if (diffA%3== 1)
					cnt1++;
				if (diffA%3== 2)
					cnt2++;
			}
			tot+= events.length;
			p.println("0: "+cnt0+"\t1: "+cnt1+"\t2: "+cnt2+"\t("+tot+"): "+((double) cnt0/ tot));
			scnt1+= cnt1;
			scnt0+= cnt0;
			scnt2+= cnt2;
			stot+= tot;
		}
		p.println("===> 0: "+scnt0+"\t1: "+scnt1+"\t2: "+scnt2+"\t("+stot+"): "+((double) scnt0/ stot));
	}
	
	
	public static void outputCodingRegions() {
		String fName= "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		boolean encode= false;
		if (fName.equals("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		Graph g= myWrapper.getGraph(encode);
		g.filterNonCodingTranscripts();
		Gene[] ge= g.getGenes();
		try {
			PrintStream buffy= new PrintStream("regions_coding.txt");
			AbstractRegion re;
			for (int i = 0; i < ge.length; i++) {
				re= ge[i].getReal5UTR();
				if (re!= null)
					buffy.println(ge[i].getChromosome()+" "
							+ re.getStart()+ " "
							+ re.getEnd()+ " "
							+ "5utr");
				re= ge[i].getRealCDS();
				if (re!= null)
					buffy.println(ge[i].getChromosome()+" "
							+ re.getStart()+ " "
							+ re.getEnd()+ " "
							+ "cds");
				re= ge[i].getReal3UTR();
				if (re!= null)
					buffy.println(ge[i].getChromosome()+" "
							+ re.getStart()+ " "
							+ re.getEnd()+ " "
							+ "3utr");
			}
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void output3UTR() {
		String fName= "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		boolean encode= false;
		if (fName.equals("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		Graph g= myWrapper.getGraph(encode);
		g.filterNonCodingTranscripts();
		Gene[] ge= g.getGenes();
		try {
			PrintStream buffy= new PrintStream("3utr_ss.txt");
			AbstractRegion re;
			for (int i = 0; i < ge.length; i++) {
				SpliceSite[] ss= 
					ge[i].getSpliceSites(SpliceSite.CONSTITUTIVE_SS, Gene.REGION_REAL_3UTR);
				for (int j = 0; ss!= null&& j < ss.length; j++) 
					buffy.println(ge[i].getChromosome()+ " "+ ss[j].getPos());
			}
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void outputCheckFile() {
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		String fName= "encode/44regions_genes_CHR_coord.gtf";
		System.out.println("Writing consistency file for: "+fName);
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		boolean encode= false;
		if (fName.equals("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		Graph g= myWrapper.getGraph(encode);
		
		PrintStream p= null;
		Calendar c= Calendar.getInstance();
			// check c.MONTH for bug?!
		String outName= "consistency_"+(c.get(c.MONTH)+ 1)+"_"+c.get(c.DAY_OF_MONTH)+"_"+c.get(c.HOUR_OF_DAY)+"_"+c.get(c.MINUTE);
		String sffx= ".chk";
		try {
			p= new PrintStream(outName+"_all"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//determineVariations(g, p);
		determineVariationsDetails(g, "isTrue", p);
		p.flush();
		p.close();
		
		try {
			p= new PrintStream(outName+"_cod"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isProteinCoding", p);
		p.flush();
		p.close();

		try {
			p= new PrintStream(outName+"_part"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isPartiallyCoding", p);
		p.flush();
		p.close();

		try {
			p= new PrintStream(outName+"_non"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isNotAtAllCoding", p);
		p.flush();
		p.close();
}
	
	public static void test02_ss_statistics() {
	
		String iname= "encode/44regions_genes_CHR_coord.gtf.conservation_high";
		Graph g= getGraph(iname);
		
		//String fName= "test02_ss_statistics_refseq_in_encode_nmd.txt";
		String fName= iname+ ".utr";
		PrintStream p= null;
		try {
			fName= Toolbox.checkFileExists(fName);
			p= new PrintStream(fName);
			//p= System.out;
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//g.filterNonCodingTranscripts();
		g.filterNMDTranscripts();
		boolean chk= false;
	
		SpliceSite[][] res= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
		p.println("tot\t"+res[0].length+"\t"+res[1].length);
		int ctrAlt5UTR= 0, ctrAltCDS= 0, ctrAlt3UTR= 0, ctrAltRest= 0;
		Vector chkV= new Vector();
		for (int i = 0; i < res[0].length; i++) {
			if (res[0][i].isConstitutive())
				continue;
			if (res[0][i].is5UTRMaxTranscript()) {
				++ctrAlt5UTR;
				chkV.add(res[0][i]);
			}
			if (res[0][i].isCDSMaxTranscript())
				++ctrAltCDS;
			if (res[0][i].is3UTRMaxTranscript())
				++ctrAlt3UTR;
			if (!res[0][i].is5UTRMaxTranscript()&& !res[0][i].isCDSMaxTranscript()
					&& !res[0][i].is3UTRMaxTranscript())
				++ctrAltRest;
		}
		int ctrCon5UTR= 0, ctrConCDS= 0, ctrCon3UTR= 0, ctrConRest= 0;
		for (int i = 0; i < res[1].length; i++) {
			if (!res[1][i].isConstitutive())
				continue;
			if (res[1][i].is5UTRMaxTranscript()) {
				++ctrCon5UTR;
				chkV.add(res[1][i]);
			}
			if (res[1][i].isCDSMaxTranscript())
				++ctrConCDS;
			if (res[1][i].is3UTRMaxTranscript())
				++ctrCon3UTR;
			if (!res[1][i].is5UTRMaxTranscript()&& !res[1][i].isCDSMaxTranscript()
					&& !res[1][i].is3UTRMaxTranscript())
				++ctrConRest;
		}
	
		p.println("5UTR_alt\t"+ctrAlt5UTR+"("+((float) ctrAlt5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")"+
				"\t5UTR_con "+ctrCon5UTR+"("+((float) ctrCon5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")");
		p.println("CDS_alt\t"+ctrAltCDS+"("+((float) ctrAltCDS/(ctrAltCDS+ ctrConCDS))+")"+
				"\tCDS_con "+ctrConCDS+"("+((float) ctrConCDS/(ctrAltCDS+ ctrConCDS))+")");
		p.println("3UTR_alt\t"+ctrAlt3UTR+"("+((float) ctrAlt3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")"+
				"\t3UTR_con "+ctrCon3UTR+"("+((float) ctrCon3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")");
		p.println("REST_alt\t"+ctrAltRest+"("+((float) ctrAltRest/(ctrAltRest+ ctrConRest))+")"+
				"\tREST_con "+ctrConRest+"("+((float) ctrConRest/(ctrAltRest+ ctrConRest))+")");
		
		if (chk) {
			try {
				PrintStream pr= new PrintStream("test02_ss_statistics_SS_5UTR_encode_coding_nmd.txt");
				for (int i = 0; i < chkV.size(); i++) {
					SpliceSite ss= (SpliceSite) chkV.elementAt(i);
					pr.println(ss.getGene().getChromosome()+" "+ss);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

		}
	}

	public static void test02_ss_statistics(Graph g, PrintStream p, boolean codTranscripts) {
		
			if (codTranscripts)
				g.filterNonCodingTranscripts();
		
			g.getASVariations(ASMultiVariation.FILTER_NONE);	// init
			
			SpliceSite[][] res= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
			p.println("tot\t"+res[0].length+"\t"+res[1].length);
			int ctrAlt5UTR= 0, ctrAltCDS= 0, ctrAlt3UTR= 0, ctrAltRest= 0;
			for (int i = 0; i < res[0].length; i++) {
				if (res[0][i].isConstitutive())
					continue;
				if (res[0][i].is5UTRMaxTranscript())
					++ctrAlt5UTR;
				if (res[0][i].isCDSMaxTranscript())
					++ctrAltCDS;
				if (res[0][i].is3UTRMaxTranscript())
					++ctrAlt3UTR;
				if (!res[0][i].is5UTRMaxTranscript()&& !res[0][i].isCDSMaxTranscript()
						&& !res[0][i].is3UTRMaxTranscript())
					++ctrAltRest;
			}
			int ctrCon5UTR= 0, ctrConCDS= 0, ctrCon3UTR= 0, ctrConRest= 0;
			for (int i = 0; i < res[1].length; i++) {
				if (!res[1][i].isConstitutive())
					continue;
				if (res[1][i].is5UTRMaxTranscript())
					++ctrCon5UTR;
				if (res[1][i].isCDSMaxTranscript())
					++ctrConCDS;
				if (res[1][i].is3UTRMaxTranscript())
					++ctrCon3UTR;
				if (!res[1][i].is5UTRMaxTranscript()&& !res[1][i].isCDSMaxTranscript()
						&& !res[1][i].is3UTRMaxTranscript())
					++ctrConRest;
			}

			p.println("5UTR_alt\t"+ctrAlt5UTR+"("+((float) ctrAlt5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")"+
					"\t5UTR_con "+ctrCon5UTR+"("+((float) ctrCon5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")");
			p.println("CDS_alt\t"+ctrAltCDS+"("+((float) ctrAltCDS/(ctrAltCDS+ ctrConCDS))+")"+
					"\tCDS_con "+ctrConCDS+"("+((float) ctrConCDS/(ctrAltCDS+ ctrConCDS))+")");
			p.println("3UTR_alt\t"+ctrAlt3UTR+"("+((float) ctrAlt3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")"+
					"\t3UTR_con "+ctrCon3UTR+"("+((float) ctrCon3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")");
			p.println("REST_alt\t"+ctrAltRest+"("+((float) ctrAltRest/(ctrAltRest+ ctrConRest))+")"+
					"\tREST_con "+ctrConRest+"("+((float) ctrConRest/(ctrAltRest+ ctrConRest))+")");
		}

	public static void test02_ss_statistics_1st_submission(Graph g, PrintStream p, boolean codTranscripts) {
		
			if (codTranscripts)
				g.filterNonCodingTranscripts();
		
			g.getASVariations(ASMultiVariation.FILTER_NONE);	// init
			
			SpliceSite[][] res= g.getSpliceSites(Gene.REGION_REAL_5UTR);
			p.println("\t#AS-SS\t#CONST-SS\tpercAS");
			p.println("5UTR: "+res[0].length+" "+res[1].length+"\t"+((float) res[0].length/(float) (res[0].length+ res[1].length)));
			int utr5= res[1].length;
	//		for (int i = 0; i < res[1].length; i++) 
	//			System.out.println(res[1][i]);
			res= g.getSpliceSites(Gene.REGION_REAL_CDS);
			p.println("CDS: "+res[0].length+" "+res[1].length+"\t"+((float) res[0].length/(float) (res[0].length+ res[1].length)));
			int cds= res[1].length;
	//		for (int i = 0; i < res[1].length; i++) 
	//			System.out.println(res[1][i]);
			res= g.getSpliceSites(Gene.REGION_REAL_3UTR);
			p.println("3UTR: "+res[0].length+" "+res[1].length+"\t"+((float) res[0].length/(float) (res[0].length+ res[1].length)));
			int utr3= res[1].length;
	//		for (int i = 0; i < res[1].length; i++) 
	//			System.out.println(res[1][i]);
			res= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
			p.println("tot: "+res[0].length+"\t"+res[1].length+"\t"+((float) res[0].length/(float) (res[0].length+ res[1].length)));
			float f1= ((float) utr5/(float) res[1].length);
			float f2= ((float) cds/(float) res[1].length);
			float f3= ((float) utr3/(float) res[1].length);
			p.println("%5UTR\t%CDS\t%3UTR\t%Rest");
			p.println(f1+"\t"+f2+"\t"+f3+"\t"+(1f- f1- f2- f3));
		}

	public static void outputFirstExonIntronAtg(Graph g, PrintStream p) {
			
			g.filterNonCodingTranscripts();
			g.getASVariations(ASMultiVariation.FILTER_NONE);	// to init alt/const SSs
			
			Gene[] ge= g.getGenes();
			for (int i = 0; i < ge.length; i++) {
				Transcript[] tr= ge[i].getTranscripts();
				for (int j = 0; j < tr.length; j++) {
					
					if (tr[j].getExons().length< 2)
						continue;
					Exon first= tr[j].getExons()[0];
					Exon second= tr[j].getExons()[1];
					int atg= tr[j].getTranslations()[0].get5PrimeEdge();
	
					p.println(
						ge[i].getChromosome()+ " "+
						ge[i].getStrand()+ " "+
						tr[j].getTranscriptID()+ " "+
						"1st_exon "+
	//					Math.abs(first.getStart())+ " "+
	//					Math.abs(first.getEnd())+ " "+
						Math.abs(first.get3PrimeEdge())+ " "+
						(first.getDonor().isConstitutive()?"const":"alt")+" "+
						(first.contains(atg)?"atg ":"non ")+
						(first.contains(atg)?first.get3PrimeEdge()- atg+" ":
							(first.get3PrimeEdge()- first.get5PrimeEdge()+ 1)+ " ")+
						(second.get5PrimeEdge()- first.get3PrimeEdge())+ " "
					);
					
					p.println(
						ge[i].getChromosome()+ " "+
						ge[i].getStrand()+ " "+
						tr[j].getTranscriptID()+ " "+
						"2nd_exon "+
	//					Math.abs(second.getStart())+ " "+
	//					Math.abs(second.getEnd())+ " "+
						Math.abs(second.get5PrimeEdge())+ " "+
						(second.getAcceptor().isConstitutive()?"const":"alt")+" "+
						(second.contains(atg)?"atg ":"non ")+
						(second.contains(atg)?Math.abs(atg- second.get5PrimeEdge())+" ":
							(second.get3PrimeEdge()- second.get5PrimeEdge()+ 1)+ " ")+
						(second.get5PrimeEdge()- first.get3PrimeEdge())+ " "
					);
					
				}
			}
		}

	public static SpliceSite[] getFirstExonIntronAtg4(Graph g, PrintStream p) {
			
			g.filterNonCodingTranscripts();
			g.getASVariations(ASMultiVariation.FILTER_NONE);	// to init alt/const SSs
			
			Gene[] ge= g.getGenes();
			Vector ssV= new Vector();
			for (int i = 0; i < ge.length; i++) {
				SpliceSite[] ss= ge[i].getSpliceSites();
				for (int j = 0; ss!= null&& ss.length>= 2&& j < 2; j++) {
					ssV.add(ss[j]);
				}
			}
			
			return (SpliceSite[]) Arrays.toField(ssV);
		}

	public static void outputFirstExonIntronAtg2(Graph g, PrintStream p) {
		
		g.filterNonCodingTranscripts(); 
		//g.filterCodingTranscripts();
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);	// to init alt/const SSs
	
		
		int ctrDon= 0, ctrAcc= 0;
		for (int i = 0; i < vars.length; i++) {
			if (vars[i][0].toString().equals("(1^ // 2^)")||
					vars[i][0].toString().equals("(1= // 2=)")) {
				for (int j = 0; j < vars[i].length; j++) {
					if (vars[i][0].toString().equals("(1^ // 2^)")
							&& vars[i][j].getFlankingSpliceSites()[0]== null) {
						vars[i][j].outputDetail(System.out);
						String out= "";
						out+= "AD "+ vars[i][j].getLengthDiff(true)+ " ";
						int atg1= vars[i][j].getTranscript1().getTranslations()[0].get5PrimeEdge();
						if (vars[i][j].getTranscript1().getExons()[0].contains(atg1))
							out+= "atg "+ (vars[i][j].getTranscript1().getExons()[0].get3PrimeEdge()- atg1);
						else if (vars[i][j].getTranscript1().getExons()[1].contains(atg1))
							out+= "2cod "+ (atg1- vars[i][j].getTranscript2().getExons()[1].get5PrimeEdge());
						else
							out+= "ncod "+ (vars[i][j].getTranscript1().getExons()[0].getLength());
						out+= " ";
						int atg2= vars[i][j].getTranscript2().getTranslations()[0].get5PrimeEdge();
						if (vars[i][j].getTranscript2().getExons()[0].contains(atg2))
							out+= "atg "+ (vars[i][j].getTranscript2().getExons()[0].get3PrimeEdge()- atg2);
						else if (vars[i][j].getTranscript2().getExons()[1].contains(atg2))
							out+= "2cod "+ (atg2- vars[i][j].getTranscript2().getExons()[1].get5PrimeEdge());
						else
							out+= "ncod "+ (vars[i][j].getTranscript2().getExons()[0].getLength());
						out+= " ";
						out+= vars[i][j].getTranscript1().getExons()[0].getLength()+ " ";
						out+= vars[i][j].getTranscript2().getExons()[0].getLength()+ " ";
						p.println(out);
						ctrDon++;
					} else if (vars[i][0].toString().equals("(1= // 2=)")
							&& vars[i][j].getFlankingSpliceSites()[0]== vars[i][j].getTranscript1().getSpliceChain()[0]&&
							vars[i][j].getFlankingSpliceSites()[0]== vars[i][j].getTranscript2().getSpliceChain()[0]) {
						vars[i][j].outputDetail(System.out);
						String out= "";
						out+= "AA "+ vars[i][j].getLengthDiff(true)+ " ";
						int atg1= vars[i][j].getTranscript1().getTranslations()[0].get5PrimeEdge();
						if (vars[i][j].getTranscript1().getExons()[1].contains(atg1))
							out+= "atg "+ (atg1- vars[i][j].getTranscript1().getExons()[1].get5PrimeEdge());
						else if(vars[i][j].getTranscript1().getExons()[0].contains(atg1))
							out+= "cod "+ (vars[i][j].getTranscript2().getExons()[1].getLength());
						else
							out+= "ncod "+ (vars[i][j].getTranscript1().getExons()[1].getLength());
						out+= " ";
						int atg2= vars[i][j].getTranscript2().getTranslations()[0].get5PrimeEdge();
						if (vars[i][j].getTranscript2().getExons()[1].contains(atg2))
							out+= "atg "+ (atg2- vars[i][j].getTranscript2().getExons()[1].get5PrimeEdge());
						else if(vars[i][j].getTranscript2().getExons()[0].contains(atg2))
							out+= "cod "+ (vars[i][j].getTranscript2().getExons()[1].getLength());
						else
							out+= "ncod "+ (vars[i][j].getTranscript2().getExons()[1].getLength());
						out+= " ";
						out+= vars[i][j].getTranscript1().getExons()[1].getLength()+ " ";
						out+= vars[i][j].getTranscript2().getExons()[1].getLength()+ " ";
						p.println(out);
						ctrAcc++;
					}
				}
			}
		}
		
		System.out.println(ctrDon+", "+ctrAcc);
	}

	public static SpliceSite[] getFirstExonIntronAtg3(Graph g, PrintStream p) {
		
		g.filterNonCodingTranscripts();
		//g.filterCodingTranscripts();
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);	// to init alt/const SSs

		Vector ssV= new Vector();
		int ctrDon= 0, ctrAcc= 0;
		for (int i = 0; i < vars.length; i++) {
			if (vars[i][0].toString().equals("(1^ // 2^)")||
					vars[i][0].toString().equals("(1= // 2=)")) {
				for (int j = 0; j < vars[i].length; j++) {
					if (vars[i][0].toString().equals("(1^ // 2^)")
							&& vars[i][j].getFlankingSpliceSites()[0]== null) {
						int k;
						for (k = 0; k < ssV.size(); k++) 
							if (ssV.elementAt(k)== vars[i][j].getSpliceChain1()[0])
								break;
						if (k== ssV.size())
							ssV.add(vars[i][j].getSpliceChain1()[0]);
						for (k = 0; k < ssV.size(); k++) 
							if (ssV.elementAt(k)== vars[i][j].getSpliceChain2()[0])
								break;
						if (k== ssV.size())
							ssV.add(vars[i][j].getSpliceChain2()[0]);
						ctrDon++;
					} else if (vars[i][0].toString().equals("(1= // 2=)")
							&& vars[i][j].getFlankingSpliceSites()[0]== vars[i][j].getTranscript1().getSpliceChain()[0]&&
							vars[i][j].getFlankingSpliceSites()[0]== vars[i][j].getTranscript2().getSpliceChain()[0]) {
						int k;
						for (k = 0; k < ssV.size(); k++) 
							if (ssV.elementAt(k)== vars[i][j].getSpliceChain1()[0])
								break;
						if (k== ssV.size())
							ssV.add(vars[i][j].getSpliceChain1()[0]);
						for (k = 0; k < ssV.size(); k++) 
							if (ssV.elementAt(k)== vars[i][j].getSpliceChain2()[0])
								break;
						if (k== ssV.size())
							ssV.add(vars[i][j].getSpliceChain2()[0]);
						ctrAcc++;
					}
				}
			}
		}
		
		System.out.println(ctrDon+", "+ctrAcc+": "+ssV.size());
		return (SpliceSite[]) Arrays.toField(ssV);		
	}

	public static void outputFirstExons(Graph g, PrintStream p) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr= ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {
				Exon first= tr[j].getExons()[0];
				if (first.getDonor()== null)
					continue;
				p.println(
					ge[i].getChromosome()+ " "+
					Math.abs(first.getStart())+ " "+
					Math.abs(first.getEnd())+ " "+
					(first.getDonor().isConstitutive()?"const":"alt")+" "+
					ge[i].getStrand()
				);
			}
		}
	}

	public static void outputFirstExons1(Graph g, PrintStream p) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr= ge[i].getTranscripts();
			//ge[i].getASVariations();
			for (int j = 0; j < tr.length; j++) {
				Exon first= tr[j].getExons()[0];
				if (first.getDonor()== null)
					continue;
				p.println(
					ge[i].getChromosome()+ " "+
					Math.abs(first.getStart())+ " "+
					Math.abs(first.getEnd())+ " "+
					(first.getDonor().isConstitutive()?"const":"alt")+" "+
					ge[i].getStrand()
				);
			}
		}
	}

	public static void outputInternalExons(Graph g, PrintStream p) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr= ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {
				if (tr[j].getExons().length< 2)
					continue;
				for (int k = 0; k < tr[j].getExons().length; k++) {
					p.println(
							ge[i].getChromosome()+ " "+
							Math.abs(tr[j].getExons()[k].getStart())+ " "+
							Math.abs(tr[j].getExons()[k].getEnd())+ " "+
							(tr[j].getExons()[k].getDonor()== null?"null":(tr[j].getExons()[k].getDonor().isConstitutive()?"don_const":"don_alt"))+" "+
							(tr[j].getExons()[k].getAcceptor()== null?"null":(tr[j].getExons()[k].getAcceptor().isConstitutive()?"acc_const":"acc_alt"))+" "+
							ge[i].getStrand()
						);
				}
			}
		}
	}

	public static void outputFirstIntrons(Graph g, PrintStream p) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr= ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {
				if (tr[j].getExons().length< 2)	// no intron
					continue;
				int start= tr[j].getExons()[0].get3PrimeEdge()+ 1;	// for pos and neg
				SpliceSite donor= tr[j].getExons()[0].getDonor();
				int end= tr[j].getExons()[1].get5PrimeEdge()- 1;
				SpliceSite acceptor= tr[j].getExons()[1].getAcceptor();
				if (!tr[j].isLeadingStrand()) {	// swap start/end for neg strand
					int h= start;
					start= end;
					end= h;
				}
					// check for (real)5UTR?
				p.println(
					ge[i].getChromosome()+ " "+
					Math.abs(start)+ " "+
					Math.abs(end)+ " "+
					(donor.isConstitutive()?"don_const":"don_alt")+" "+
					(acceptor.isConstitutive()?"acc_const":"acc_alt")+" "+
					ge[i].getStrand()
				);
			}
		}
	}

	public static void outputInternalIntrons(Graph g, PrintStream p) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] tr= ge[i].getTranscripts();
			for (int j = 0; j < tr.length; j++) {
				
				DirectedRegion[] introns= tr[j].getIntrons();
				Exon[] exons= tr[j].getExons();
				for (int k = 1; introns!= null&& k < introns.length; k++) {
					p.println(
							ge[i].getChromosome()+ " "+
							Math.abs(introns[k].getStart())+ " "+
							Math.abs(introns[k].getEnd())+ " "+
							(exons[k].getDonor().isConstitutive()?"don_const":"don_alt")+" "+
							(exons[k+1].getAcceptor().isConstitutive()?"acc_const":"acc_alt")+" "+
							ge[i].getStrand()
						);
				}
			}
		}
	}

	public static void test02_ss_statistics_outCDSalt(Graph g) {
			g.getASVariations(ASMultiVariation.FILTER_NONE);
			SpliceSite[][] res= g.getSpliceSites(Gene.REGION_REAL_CDS);
			try {
				PrintStream p= new PrintStream("check_as_cds");
				for (int i = 0; i < res[0].length; i++) 
					p.println(res[0][i].toOutString());
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			

		}

	public static void test02b_ss_statistics_3P(Graph g, PrintStream p) {
			g.getASVariations(ASMultiVariation.FILTER_NONE);
			SpliceSite[][] res= g.getSpliceSites(Gene.REGION_REAL_5UTR);

			res= g.getSpliceSites(Gene.REGION_REAL_3UTR);
			int ir= 0;
			for (int i = 0; i < res[0].length; i++) {
				ASVariation[] as= res[0][i].getAsVars();
				for (int j = 0; j < as.length; j++) 
					if (as[j].toString().equals("( // 1^2=)")) { 
						++ir;
						break;
					}
				
			}
			
			p.println("SSir "+ir+", SSas "+res[0].length+", SStot "+res[1].length);
		}

	/** @deprecated ?
	 * 
	 * 
	 * @param g
	 * @return
	 */
	public static double getPercASSpliceSites(Graph g) {
		Gene[] ge= g.getGenes();
		int altPS= 0;
		int totPS= 0;
		Vector v= new Vector();
		for (int i = 0; i < ge.length; i++) {
			SpliceSite[] ps= ge[i].getSpliceSites();
			if (ps== null)
				continue;
			totPS+= ps.length;
			for (int j = 0; j < ps.length; j++) {
				v.add(ps[j]);
				if (!ps[j].isConstitutive()) {
					++altPS;
				}
			}
		}
		BufferedWriter buffy= null;
		try {
			buffy = new BufferedWriter(new FileWriter("ss.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		for (int i = 0; i < v.size(); i++) {
			SpliceSite ss= (SpliceSite) v.elementAt(i);
			System.out.println(ss.getGene().getChromosome()+" "+ss);
			try {
				buffy.write(ss.getGene().getChromosome()+" "+ss+"\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		try {
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println(altPS+","+totPS);
		return (double) altPS/ totPS;
	}
	
	public static void diff(ASVariation[][] classes1, ASVariation[][] classes2) {
		int cnt= 0;
		Comparator compi= null; //new ASVariation.IdentityComparator();
		for (int i = 0; i < classes1.length; i++) {
			String id= classes1[i][0].toString();
				// check missing in 2nd file
			int j;
			for (j = 0; j < classes2.length; j++) 
				if (classes2[j][0].toString().equals(id))
					break;
			if (j== classes2.length) {
				System.err.println("Missing in 2nd file "+id);
				for (int k = 0; k < classes2.length; k++) 
					classes1[i][k].outputDetail(System.out);
				cnt+= classes1[i].length;
				continue;
			}
			
				// check in 1st file
			for (int k = 0; k < classes1[i].length; k++) {
				int m;
				for (m = 0; m < classes2[j].length; m++) 
					if (compi.compare(classes1[i][k], classes2[j][m])== 0)
						break;
				if (m== classes2.length) {
					System.err.println("Missing in 2nd file ");
					classes1[i][k].outputDetail(System.out);
					++cnt;
				}
			}
			// check in 2nd file
			for (int k = 0; k < classes2[j].length; k++) {
				int m;
				for (m = 0; m < classes1[i].length; m++) 
					if (compi.compare(classes1[i][k], classes2[j][m])== 0)
						break;
				if (m== classes2.length) {
					System.err.println("Missing in 1st file ");
					classes1[i][k].outputDetail(System.out);
					++cnt;
				}
			}
		}
		
		//check missing in 1st file 
		for (int i = 0; i < classes2.length; i++) {
			String id= classes2[i][0].toString();
			int j;
			for (j = 0; j < classes1.length; j++) 
				if (classes1[j][0].toString().equals(id))
					break;
			if (j== classes1.length) {
				System.err.println("Missing in 1st file "+id);
				for (int k = 0; k < classes1[j].length; k++) 
					classes1[i][k].outputDetail(System.out);
				cnt+= classes1[i].length;
				continue;
			}
		}
			

	}
	
	public static void analyzeGene(Graph g, String geneID, String methodName, int filter) {
		Gene ge= g.getGene(geneID, "human");
		System.out.println("--- Report for Gene "+geneID+" ("+methodName+"): ---");
		ASVariation[] as= ge.getASVariations(filter);
		Method m= null;
		try {
			m = as[0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		as= ((ASVariation[][]) Arrays.filter( new ASVariation[][] {as}, m))[0];
		
		
		for (int i = 0; i < as.length; i++) {
			as[i].outputDetail(System.out);
			System.out.println();
		}
	}

	public final static String INPUT_ENCODE= "encode/44regions_genes_CHR_coord.gtf";
	public final static String INPUT_ENCODE_RACES= "encode/gencode_races.gtf";
	public final static String INPUT_ENSEMBL_CODING_FROM_UCSC= "encode/EnsemblGenes_fromUCSC.gtf";
	public final static String INPUT_ENSEMBL_CODING_FROM_UCSC_ENCODE_ONLY= "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf";
	public final static String INPUT_REFSEQ_CODING_FROM_UCSC_ENCODE_ONLY= "encode/RefSeqGenes_fromUCSC.inENCODE.gtf";
	public final static String INPUT_REFSEQ_CODING_FROM_UCSC= "encode/RefSeqGenes_fromUCSC.gtf";
	public final static String INPUT_ENSEMBL_FROM_ENSEMBL= "encode/EnsemblGenes_all_fromENSEMBL.gtf";
	public final static String INPUT_HAVANA= "encode/Sequences_mapped_HAVANA_136.gtf";
	public final static String INPUT_ASD= "Sequences_mapped_HAVANA_136.gtf";
	public final static String[] INPUT_GO_LEVEL1= new String[] {
		   "0008150",	// biological_process [137651]  
           "0005575", 	// cellular_component [115920]  
           "0003674" 	// molecular_function [136131]  
	};
	public final static String[] INPUT_GO_LEVEL2= new String[] {
		
		  "0009987",	// : cellular process [79773] 
		  "0007275",	// : development [14733] 
		  "0021700",	// : developmental maturation [219] 
		  "0040007",	// : growth [3472] 
		  "0051704",	// : interaction between organisms [1562] 
		  "0007582",	// : physiological process [82865] 
		  "0043473",	// : pigmentation [126] 
		  "0050789",	// : regulation of biological process [17317] 
		  "0000003",	// : reproduction [4588] 
		  "0050896",	// : response to stimulus [16684] 
		  "0016032",	// : viral life cycle [332] 
		  "0005623",	// : cell [73091] 
		  "0044464",	// : cell part [73053] 
		  "0031975",	// : envelope [2954] 
		  "0031012",	// : extracellular matrix [584] 
		  "0044420",	// : extracellular matrix part [284] 
		  "0005576",	// : extracellular region [6259] 
		  "0044421",	// : extracellular region part [4726] 
		  "0031974",	// : membrane-enclosed lumen [4477] 
		  "0043226",	// : organelle [47321] 
		  "0044422",	// : organelle part [15544] 
		  "0043234",	// : protein complex [13680] 
		  "0045202",	// : synapse [270] 
		  "0044456",	// : synapse part [112] 
		  "0019012",	// : virion [154] 
		  "0044423",	// : virion part [121] 
		  "0016209",	// : antioxidant activity [535] 
		  "0015457",	// : auxiliary transport protein activity [169] 
		  "0005488",	// : binding [38199] 
		  "0003824",	// : catalytic activity [44092] 
		  "0030188",	// : chaperone regulator activity [54] 
		  "0042056",	// : chemoattractant activity [13] 
		  "0045499",	// : chemorepellant activity [8] 
		  "0031992",	// : energy transducer activity [0] 
		  "0030234",	// : enzyme regulator activity [2420] 
		  "0003774",	// : motor activity [563] 
		  "0045735",	// : nutrient reservoir activity [51] 
		  "0031386",	// : protein tag [18] 
		  "0004871",	// : signal transducer activity [7817] 
		  "0005198",	// : structural molecule activity [4022] 
		  "0030528",	// : transcription regulator activity [8944] 
		  "0045182",	// : translation regulator activity [844] 
		  "0005215",	// : transporter activity [9657] 
		  "0030533"	// : triplet codon-amino acid adaptor activity [1269]
	};
	
	static void addVirtualGene1(Graph g) {
		Gene newGene= new Gene(g.getSpeciesByEnsemblPrefix("ENS"), "ENSG00001000001");
		Transcript trans1= new Transcript(newGene, "ENST00000100001");
		newGene.addTranscript(trans1);
		Transcript trans2= new Transcript(newGene, "ENST00000100002");
		newGene.addTranscript(trans2);
		Exon e= new Exon(trans1, "ENSE00000100001", 10,40);
		trans1.addExon(e);
		trans2.addExon(e);
		e= new Exon(trans1, "ENSE00000100001", 50,150);
		trans1.addExon(e);
		e= new Exon(trans1, "ENSE00000100002", 100,250);
		trans2.addExon(e);
		e= new Exon(trans1, "ENSE00000100003", 200,250);
		trans1.addExon(e);
		
		g.addGene(newGene);
	}
	
	public static Graph getGraph(String fName)  {
		EncodeWrapper myWrapper;
		if (fName.equals(INPUT_ASD))
			myWrapper= new ATDWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		else 
			myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		boolean encode= false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode= true;
		
		Graph g= myWrapper.getGraph(encode);		// <===== check ENCODE here !!!
		return g;
	}
	public static Graph getGraph_now(String fName)  {
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		boolean encode= false;
		if (fName.startsWith("encode/44regions_genes_CHR_coord"))
			encode= true;
		
		Graph g= myWrapper.getGraph(encode);		// <===== check ENCODE here !!!
		return g;
	}	
	public static ASVariation[] getVariation(String s, ASVariation[][] classes) {
		for (int i = 0; i < classes.length; i++) {
			if (classes[i][0].toString().equals(s))
				return classes[i];
		}
		return null;
	}
	
	public static ASVariation[] getVariation(String s, ASVariation[] classes) {
		if (classes== null)
			return null;
		
		Vector v= new Vector();
		for (int i = 0; i < classes.length; i++) {
			if (classes[i].toString().equals(s))
				v.add(classes[i]);
		}
		return (ASVariation[]) Arrays.toField(v);
	}
	
	public static ASVariation[] getAlternateDonors(ASVariation[][] classes) {
		Vector v= new Vector();
		for (int i = 0; i < classes.length; i++) {
			if ((classes[i][0].toString().startsWith("(1^")&& 
					classes[i][0].toString().contains("// 2^"))
					|| (classes[i][0].toString().startsWith("(2^")&& 
							classes[i][0].toString().contains("// 1^")))
				v.add(classes[i]);
		}
		return (ASVariation[]) Arrays.toField(v);
	}
	

	public static ASVariation[][] getSpecificVariations(Graph g, PrintStream p) {
		
		ASVariation[][] classes= g.getASVariations(0);
		Comparator compi= new ASVariation.UniqueCodingComparator();
		Comparator c= new ASVariation.CodingHierarchyFilter();
		ASVariation[][] newClasses= new ASVariation[classes.length][];
		for (int i = 0; i < classes.length; i++) {
			Vector v= new Vector();
			for (int j = 0; j < classes[i].length; j++) {
				for (int k = (j+1); k < classes[i].length; k++) {
					if (compi.compare(classes[i][j], classes[i][k])== 0) {	 
						v.add(classes[i][j]);	// keep one
					} else {
						for (int m = 0; m < v.size(); m++) {				
							if ((c.compare(v.elementAt(m), classes[i][k])== 0)
									|| (c.compare(v.elementAt(m), classes[i][j])== 0))
								v.remove(m);	// remove both
						}
					}
				}
			}
			newClasses[i]= (ASVariation[]) Arrays.toField(v);
		}
		
		return outputVariations(newClasses, true, false, System.out);
	}
	/**
		 * @param g
		 * @param respectCoding
		 * @return
		 */
		public static void test04_determineVariations(Graph g, PrintStream p) {
			
			ASVariation[][] classes= g.getASVariations(1);
			classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
			
			try {
				Method m = classes[0][0].getClass().getMethod("isTrue", null);
				ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
				outputVariationsDegree(filtClasses, p);	//false, true, 
				
//				m = classes[0][0].getClass().getMethod("isProteinCoding_old_publ", null);
//				filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
//				outputVariations(filtClasses, true, false, p);
//		
//				m = classes[0][0].getClass().getMethod("isPartiallyCoding", null);
//				filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
//				outputVariations(filtClasses, true, false, p);
//		
//				m = classes[0][0].getClass().getMethod("isNotAtAllCoding", null);
//				filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
//				outputVariations(filtClasses, true, false, p);
				
	//			Comparator compi= new Comparator()
	//			for (int i = 0; i < filtClasses.length; i++) {
	//				for (int j = 0; j < filtClasses[i].length; j++) {
	//					for (int k = 0; k < filtClasses1.length; k++) {
	//						if (!filtClasses[i][0].toString().equals(filtClasses1[k][0].toString()))
	//							continue;
	//						for (int m = 0; m < filtClasses1[k].length; m++) {
	//							if (filtClasses[i][j])
	//						}
	//					}
	//				}
	//			}
				
			} catch (Exception e) {
				e.printStackTrace();
			} 
		}


		/**
			 * @param g
			 * @param respectCoding
			 * @return
			 */
			public static void test04_determineVariations_rev(Graph g, PrintStream p) {
				
				ASVariation[][] classes= g.getASVariations(1);
				classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
				
				try {
					
					Method m = classes[0][0].getClass().getMethod("isTrue", null);
		
					ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, false, true, p);
					
					m = classes[0][0].getClass().getMethod("isProteinCoding", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
			
					m = classes[0][0].getClass().getMethod("isNotProteinCoding", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
			
					m = classes[0][0].getClass().getMethod("isCompletelyIn5UTR", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
					m = classes[0][0].getClass().getMethod("isCompletelyIn3UTR", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
					m = classes[0][0].getClass().getMethod("isTwilightZone", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
					m = classes[0][0].getClass().getMethod("isCodingFunc", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
			
					m = classes[0][0].getClass().getMethod("is5UTRFunc", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
			
					m = classes[0][0].getClass().getMethod("is3UTRFunc", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
			
					m = classes[0][0].getClass().getMethod("isTwilightFunc", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
	//				for (int i = 0; i < filtClasses.length; i++) {
	//					if (filtClasses[i]!= null&& filtClasses[i].length> 0) 
	//						for (int j = 0; j < filtClasses[i].length; j++) {
	//							System.out.println(filtClasses[i][j].getTranscript1()+"\t"+filtClasses[i][j].getTranscript2()+"\t"+filtClasses[i][j].getSpliceUniverse()[0].getPos());
	//						}
	//				}
					
			
					m = classes[0][0].getClass().getMethod("isTwilightCDS", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
					m = classes[0][0].getClass().getMethod("isTwilight5UTR", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
		
					m = classes[0][0].getClass().getMethod("isTwilight3UTR", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
					m = classes[0][0].getClass().getMethod("isTwilightSpooky", null);
					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
					p.println(m.getName());
					outputVariations(filtClasses, true, false, p);
					
		//			Comparator compi= new Comparator()
		//			for (int i = 0; i < filtClasses.length; i++) {
		//				for (int j = 0; j < filtClasses[i].length; j++) {
		//					for (int k = 0; k < filtClasses1.length; k++) {
		//						if (!filtClasses[i][0].toString().equals(filtClasses1[k][0].toString()))
		//							continue;
		//						for (int m = 0; m < filtClasses1[k].length; m++) {
		//							if (filtClasses[i][j])
		//						}
		//					}
		//				}
		//			}
					
				} catch (Exception e) {
					e.printStackTrace();
				} 
			}

		/**
				 * @param g
				 * @param respectCoding
				 * @return
				 */
				public static void test04_determineVariations_nmd() {
					
					//"encode/44regions_genes_CHR_coord.gtf.TissueSpecificity_high";
					//"encode/44regions_genes_CHR_coord.gtf.conservation_high";
					String iname= INPUT_ENCODE;
					String sfx= "_GENCODE.landscape";			
					Graph g= getGraph(iname);
					
					PrintStream p= null;
					
					//g.filterNonCodingTranscripts();
					g.filterNMDTranscripts();
					ASVariation[][] classes= g.getASVariations(ASMultiVariation.FILTER_NONE);
					classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
					ASVariation[][] filtClasses= new ASVariation[classes.length][];
					Comparator compx= new ASVariation.HierarchyMergeFlags();
					for (int i = 0; i < filtClasses.length; i++) 
						filtClasses[i]= ASMultiVariation.removeRedundancyHierachically(classes[i], compx);
					classes= filtClasses;
					
					try {
						Method m = classes[0][0].getClass().getMethod("is_all", null);	// warum eigentlich?
						try {
							String fName= Toolbox.checkFileExists(iname+"."+m.getName()+sfx);
							p= new PrintStream(fName);
							//p= System.out;
						} catch (Exception e) {
							e.printStackTrace();
						}
						Comparator compi= new ASVariation.StructureComparator();
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(classes[i], compi);
						filtClasses= (ASVariation[][]) Arrays.filter(filtClasses, m);
						p.println(m.getName());
						outputVariations(filtClasses, false, true, p);
						p.flush(); p.close();
						
//						m = classes[0][0].getClass().getMethod("isUTRnotLastIntron", null);
//						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
//						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
//							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
//						p.println(m.getName());
//						outputVariations(filtClasses, true, false, p);

						
						
						// "isProteinCoding_1cover"
						m = classes[0][0].getClass().getMethod("is_affecting_CDS", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						p.println(m.getName());
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
				
						m = classes[0][0].getClass().getMethod("is_contained_in_CDS", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
				
						// "isNotProteinCoding_1cover"
						m = classes[0][0].getClass().getMethod("is_affecting_5UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//				System.out.println(m.getName());
		//				for (int i = 0; i < filtClasses.length; i++) {
		//					if (filtClasses[i]!= null&& filtClasses[i][0].toString().equals("(1^ // 2^)"))
		//						for (int j = 0; j < filtClasses[i].length; j++) {
		//							System.out.println(filtClasses[i][j].getTranscript1()+","+filtClasses[i][j].getTranscript2());
		//						}
		//				}
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
				
						m = classes[0][0].getClass().getMethod("is_contained_in_5UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
				
						m = classes[0][0].getClass().getMethod("is_affecting_3UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
						m = classes[0][0].getClass().getMethod("is_contained_in_3UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
				
						m = classes[0][0].getClass().getMethod("is_twilight_5UTR_CDS", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
						m = classes[0][0].getClass().getMethod("is_twilight_CDS_3UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
						m = classes[0][0].getClass().getMethod("is_twilight_5UTR_3UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
						m = classes[0][0].getClass().getMethod("is_twilight_5UTR_CDS_3UTR", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
						m = classes[0][0].getClass().getMethod("is_non_classifiable", null);
						try {p= new PrintStream(iname+"."+m.getName()+sfx);
						} catch (Exception e) {e.printStackTrace();}
						filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//				for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//					filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
						p.println(m.getName());
						outputVariations(filtClasses, true, false, p);
						p.flush(); p.close();
						
					} catch (Exception e) {
						e.printStackTrace();
					} 
				}

		/**
				 * @param g
				 * @param respectCoding
				 * @return
				 */
				public static void test04_determineVariations_go() {
		
//					Graph g= getGraph( "encode/RefSeqGenes_fromUCSC_geneid.gtf");
//					g.filterNMDTranscripts();
//					ASVariation[][] classes= g.getASVariations(ASMultiVariation.FILTER_NONE);
//					classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
					Comparator compi= new ASVariation.StructureComparator();
//					for (int i = 0; i < classes.length; i++) 
//						classes[i]= ASMultiVariation.removeRedundancy(classes[i], compi);
//					outputVariations(classes, false, true, System.out);
					
					//String[] extensions= {"Brain", "Colon", "Intest", "Liver", "Kidney", "Stomach"};
					int[] extensions= GOFragmenter.GO_ALL;
					String bName= "encode/44regions_genes_CHR_coord.gtf";	// INPUT_ENCODE;
					for (int x = 0; x < extensions.length; x++) {
						String fName= bName+"."+extensions[x];	// *
						File f= new File(fName);
						System.out.print("parsing "+fName+"..");
						if (f.length()< 1) {
							System.out.println("skipped.");
							continue;
						}
						EncodeWrapper myWrapper= new EncodeWrapper(f.getAbsolutePath()); // testGTF.gtf
						try {
							myWrapper.read();
						} catch (Exception e) {
							e.printStackTrace(); 
						}
						boolean encode= false;
						if (fName.contains("44regions_genes_CHR_coord.gtf"))
							encode= true;
						Graph g= myWrapper.getGraph(encode);		// <===== check ENCODE here !!!
						//g.filterNonCodingTranscripts();
						g.filterNMDTranscripts();
						
						PrintStream p= null;
						try {
							//fName= Toolbox.checkFileExists(fName);
							p= new PrintStream(fName+ ".landscape");
							//p= System.out;
						} catch (Exception e) {
							e.printStackTrace();
						}
						
						ASVariation[][] classes= g.getASVariations(ASMultiVariation.FILTER_NONE);
						classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
						if (classes== null)
							continue;
						ASVariation[][] filtClasses= new ASVariation[classes.length][];
						Comparator compx= new ASVariation.HierarchyMergeFlags();
						for (int i = 0; i < filtClasses.length; i++) 
							filtClasses[i]= ASMultiVariation.removeRedundancyHierachically(classes[i], compx);
						classes= filtClasses;
						
						try {
							Method m = classes[0][0].getClass().getMethod("isTrue", null);	// warum eigentlich?
							for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
								filtClasses[i]= ASMultiVariation.removeRedundancy(classes[i], compi);
							filtClasses= (ASVariation[][]) Arrays.filter(filtClasses, m);
							filtClasses= (ASVariation[][]) Arrays.sort2DFieldRev(filtClasses);
							p.println(m.getName());
							outputVariations(filtClasses, false, false, p);
							
		//					// "isProteinCoding_1cover"
		//					m = classes[0][0].getClass().getMethod("isCDSRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//			
		//					m = classes[0][0].getClass().getMethod("isCDSNonRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//			
		//					// "isNotProteinCoding_1cover"
		//					m = classes[0][0].getClass().getMethod("is5UTRRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//			
		//					m = classes[0][0].getClass().getMethod("is5UTRNonRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//			
		//					m = classes[0][0].getClass().getMethod("is3UTRRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
		//					m = classes[0][0].getClass().getMethod("is3UTRNonRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//			
		//					m = classes[0][0].getClass().getMethod("isTwilightRedundant5CDS", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
		//					m = classes[0][0].getClass().getMethod("isTwilightRedundantCDS3", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
		//					m = classes[0][0].getClass().getMethod("isTwilightRedundant53", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
		//					m = classes[0][0].getClass().getMethod("isTwilightRedundantAll", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					for (int i = 0; filtClasses!= null&& i < filtClasses.length; i++) 
		//						filtClasses[i]= ASMultiVariation.removeRedundancy(filtClasses[i], compi);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
		//					m = classes[0][0].getClass().getMethod("isNothingRedundant", null);
		//					filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
		//					p.println(m.getName());
		//					outputVariations(filtClasses, true, false, p);
		//					
						} catch (Exception e) {
							e.printStackTrace();
						}
						System.out.println("done.");
					}
				}

		/**
		 * @param g
		 * @param respectCoding
		 * @return
		 */
		public static void test02_UTRvsCDS_go() {

			//String[] extensions= {"Brain", "Colon", "Intest", "Liver", "Kidney", "Stomach"};
			int[] extensions= GOFragmenter.GO_ALL;
			// "encode/44regions_genes_CHR_coord.gtf.expres_high";	
			// INPUT_ENCODE;
			String bName= "encode/RefSeqGenes_fromUCSC_geneid.gtf";	
			for (int x = 0; x < extensions.length; x++) {
				String fName= bName+"."+extensions[x];	// *
				System.out.print("parsing "+fName);
				File f= new File(fName);
				if (f.length()< 1) {
					System.out.println("skipped.");
					continue;
				}
				EncodeWrapper myWrapper= new EncodeWrapper(f.getAbsolutePath()); // testGTF.gtf
				try {
					myWrapper.read();
				} catch (Exception e) {
					e.printStackTrace(); 
				}
				boolean encode= false;
				if (fName.contains("44regions_genes_CHR_coord.gtf"))
					encode= true;
				Graph g= myWrapper.getGraph(encode);		// <===== check ENCODE here !!!
				//g.filterNonCodingTranscripts();
				g.filterNMDTranscripts();
				
				PrintStream p= null;
				try {
					//fName= Toolbox.checkFileExists(fName);
					p= new PrintStream(fName+ ".utr");
					//p= System.out;
				} catch (Exception e) {
					e.printStackTrace();
				}
				
				SpliceSite[][] res= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
				if (res== null)
					continue;
				p.print("tot\t");
				if (res[0]== null)
					p.print("0\t");
				else
					p.print(res[0].length+"\t");
				if (res[1]== null)
					p.println("0");
				else
					p.println(res[1].length);
				int ctrAlt5UTR= 0, ctrAltCDS= 0, ctrAlt3UTR= 0, ctrAltRest= 0;
				Vector chkV= new Vector();
				for (int i = 0; res[0]!= null&& i < res[0].length; i++) {
					if (res[0][i].isConstitutive())
						continue;
					if (res[0][i].is5UTRMaxTranscript()) {
						++ctrAlt5UTR;
						chkV.add(res[0][i]);
					}
					if (res[0][i].isCDSMaxTranscript())
						++ctrAltCDS;
					if (res[0][i].is3UTRMaxTranscript())
						++ctrAlt3UTR;
					if (!res[0][i].is5UTRMaxTranscript()&& !res[0][i].isCDSMaxTranscript()
							&& !res[0][i].is3UTRMaxTranscript())
						++ctrAltRest;
				}
				int ctrCon5UTR= 0, ctrConCDS= 0, ctrCon3UTR= 0, ctrConRest= 0;
				for (int i = 0; res[1]!= null&& i < res[1].length; i++) {
					if (!res[1][i].isConstitutive())
						continue;
					if (res[1][i].is5UTRMaxTranscript()) {
						++ctrCon5UTR;
						chkV.add(res[1][i]);
					}
					if (res[1][i].isCDSMaxTranscript())
						++ctrConCDS;
					if (res[1][i].is3UTRMaxTranscript())
						++ctrCon3UTR;
					if (!res[1][i].is5UTRMaxTranscript()&& !res[1][i].isCDSMaxTranscript()
							&& !res[1][i].is3UTRMaxTranscript())
						++ctrConRest;
				}
			
				p.println("5UTR_alt\t"+ctrAlt5UTR+"("+((float) ctrAlt5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")"+
						"\t5UTR_con "+ctrCon5UTR+"("+((float) ctrCon5UTR/(ctrAlt5UTR+ ctrCon5UTR))+")");
				p.println("CDS_alt\t"+ctrAltCDS+"("+((float) ctrAltCDS/(ctrAltCDS+ ctrConCDS))+")"+
						"\tCDS_con "+ctrConCDS+"("+((float) ctrConCDS/(ctrAltCDS+ ctrConCDS))+")");
				p.println("3UTR_alt\t"+ctrAlt3UTR+"("+((float) ctrAlt3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")"+
						"\t3UTR_con "+ctrCon3UTR+"("+((float) ctrCon3UTR/(ctrAlt3UTR+ ctrCon3UTR))+")");
				p.println("REST_alt\t"+ctrAltRest+"("+((float) ctrAltRest/(ctrAltRest+ ctrConRest))+")"+
						"\tREST_con "+ctrConRest+"("+((float) ctrConRest/(ctrAltRest+ ctrConRest))+")");
				
			}
		}

		/**
		 * @param g
		 * @param respectCoding
		 * @return
		 */
		public static void test04_checkOutsideEncode(Graph g, PrintStream p) {
			
			ASVariation[][] classes= g.getASVariations(ASMultiVariation.FILTER_CODING_REDUNDANT);
			classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
			
			try {
				ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, "is5UTRFunc");
				ENCODE gp= new ENCODE();
				for (int i = 0; i < filtClasses.length; i++) {
					if (filtClasses[i][0].toString().equals(ASVariation.ID_PURE_AD)||
							filtClasses[i][0].toString().equals(ASVariation.ID_PURE_AA)) {
						int ctr= 0;
						for (int j = 0; j < filtClasses[i].length; j++) {
							if (gp.isInEncode(filtClasses[i][j].getVariationArea().getAbsoluteRegion()))
								++ctr;
						}
						p.println(filtClasses[i][0].toString()+ "\t"+ ctr+"\t"+ filtClasses[i].length);
					}
				}
				
				filtClasses= (ASVariation[][]) Arrays.filter(classes, "isCodingFunc");
				for (int i = 0; i < filtClasses.length; i++) {
					if (filtClasses[i][0].toString().equals(ASVariation.ID_PURE_AD)||
							filtClasses[i][0].toString().equals(ASVariation.ID_PURE_AA)) {
						int ctr= 0;
						for (int j = 0; j < filtClasses[i].length; j++) {
							if (gp.isInEncode(filtClasses[i][j].getVariationArea().getAbsoluteRegion()))
								++ctr;
						}
						p.println(filtClasses[i][0].toString()+ "\t"+ ctr+"\t"+ filtClasses[i].length);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
 		}

	/**
		 * @param g
		 * @param respectCoding
		 * @return
		 */
		public static void determineMultiVariations(Graph g, PrintStream p) {
			
			ASMultiVariation[][] mVar= g.getASMultiVariations2();
			mVar= (ASMultiVariation[][]) Arrays.sort2DFieldRev(mVar);
			
			try {
				outputVariations(mVar, p);
				
			} catch (Exception e) {
				e.printStackTrace();
			} 
		}

	/**
	 * @param g
	 * @param respectCoding
	 * @return
	 */
	public static void test04a_determineVarDegree(Graph g, PrintStream p) {
		
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		
		try {
			Method m = classes[0][0].getClass().getMethod("isTrue", null);
			ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
			p.println("all");
			outputDegree(filtClasses, p);
			
			m = classes[0][0].getClass().getMethod("isProteinCoding", null);
			filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
			p.println("cds");
			outputDegree(filtClasses, p);
	
			m = classes[0][0].getClass().getMethod("isPartiallyCoding", null);
			filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
			p.println("part");
			outputDegree(filtClasses, p);
	
			m = classes[0][0].getClass().getMethod("isNotAtAllCoding", null);
			filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
			p.println("utr");
			outputDegree(filtClasses, p);
			
		} catch (Exception e) {
			e.printStackTrace();
		} 
	}

	/**
	 * @param g
	 * @param respectCoding
	 * @return
	 */
	public static void test05_GO_splice_utr() {
		String fName= "test05_go_splice_utr.txt";			
		PrintStream p= null;
		try {
			fName= Toolbox.checkFileExists(fName);
			p= new PrintStream(fName);
			//p= System.out;
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		Graph g= getGraph(INPUT_ENCODE);
		g.filterNonCodingTranscripts();
		g.filterNMDTranscripts();
		Gene[] ge= g.getGenes();
		p.println("transcriptID // flag spliced 5UTR // flag AS 5UTR");
		for (int i = 0; i < ge.length; i++) {
			Transcript[] trpts= ge[i].getTranscripts();
			ASMultiVariation[] asMVars= ge[i].getASMultiVariations();
			for (int j = 0; j < trpts.length; j++) {
				int sflag= 0;
				int asflag= 0;
				if (trpts[j].getSpliceChain()== null|| trpts[j].getSpliceChain().length> 0) {
					if (trpts[j].getSpliceChain()[0].getPos()< trpts[j].getTranslations()[0].get5PrimeEdge())
						sflag= 1;
					for (int k = 0; asMVars!= null&& k < asMVars.length; k++) {
						ASVariation[] asVars= asMVars[k].getASVariationsAll();
						int m;
						for (m = 0; m < asVars.length; m++) {
							if (asVars[m].getTranscript1()!= trpts[j]&& asVars[m].getTranscript2()!= trpts[j])
								continue;
							SpliceSite[] su= asVars[m].getSpliceUniverse();
							if (su[su.length- 1].getPos()< trpts[j].getTranslations()[0].get5PrimeEdge()) {
								asflag= 1;
								break;
							}
						}
						if (m< asVars.length)
							break;
					}
				}
				p.println(trpts[j].getTranscriptID()+ "\t"+ sflag+"\t"+asflag);
			}
		}
	
		p.flush(); p.close();
	}

	/**
		 * @param g
		 * @param respectCoding
		 * @return
		 */
		public static void test05_GO_splice_utr_sigTest() {
			
				// "encode/RefSeqGenes_fromUCSC_geneid.gtf"
			Graph g= getGraph(INPUT_ENCODE);
			
			//g.filterNonCodingTranscripts();
			g.filterNMDTranscripts();
			HashMap mapUTRspliced= new HashMap();
			HashMap mapUTRas= new HashMap();
			HashMap mapAll= new HashMap();
			Gene[] ge= g.getGenes();
			for (int i = 0; i < ge.length; i++) {
				Transcript[] trpts= ge[i].getTranscripts();
				ASMultiVariation[] asMVars= ge[i].getASMultiVariations();
				for (int j = 0; j < trpts.length; j++) {
					if (mapAll.get(trpts[j].getHUGO())== null)
						mapAll.put(trpts[j].getHUGO().toUpperCase(), trpts[j].getHUGO().toUpperCase());
					int sflag= 0;
					int asflag= 0;
					if (trpts[j].getSpliceChain()== null|| trpts[j].getSpliceChain().length> 0) {
						if (trpts[j].getTranslations()!= null&&
								trpts[j].getSpliceChain()[0].getPos()< trpts[j].getTranslations()[0].get5PrimeEdge())
							sflag= 1;
						for (int k = 0; asMVars!= null&& k < asMVars.length; k++) {
							ASVariation[] asVars= asMVars[k].getASVariationsAll();
							int m;
							for (m = 0; m < asVars.length; m++) {
								if (asVars[m].getTranscript1()!= trpts[j]&& asVars[m].getTranscript2()!= trpts[j])
									continue;
								SpliceSite[] su= asVars[m].getSpliceUniverse();
								if (trpts[j].getTranslations()!= null&&
										su[su.length- 1].getPos()< trpts[j].getTranslations()[0].get5PrimeEdge()) {
									asflag= 1;
									break;
								}
							}
							if (m< asVars.length)
								break;
						}
					}
					if ((sflag== 1)&& (mapUTRspliced.get(trpts[j].getHUGO())== null))
						mapUTRspliced.put(trpts[j], trpts[j].getHUGO().toUpperCase());
					if ((asflag== 1)&& (mapUTRas.get(trpts[j].getHUGO())== null))
						mapUTRas.put(trpts[j], trpts[j].getHUGO().toUpperCase());
					
					if (trpts[j].getHUGO()== null)
						System.out.println(trpts[j].getTranscriptID());
					
				}
			}

			Vector testNodes= new Vector(mapUTRas.keySet());
			System.out.println(testNodes.size());
//			StringTokenizer st= new StringTokenizer("AC000123.1	AC004500.2	AC009955.4	AC019011.2	AC053503.10	AC068580.3	AC113617.1	ACSL6	APOA1	APOC3	APXL2	C16ORF9	C20ORF128	C20ORF173	C20ORF52	C21ORF13	C21ORF49	C21ORF59	C21ORF62	CATSPER2	CEP2	CKMT1	CLDN12	CPNE1	DNASE1L1	EEF1A1	EHD1	ELL3	FAM3A	G6PD	GART	IFNAR2	IGF2	IKBKG	INS	ITGB4BP	KATNAL1	LENG5	LENG8	LILRA2	LILRB1	LSP1	MEN1	MRPL28	MTO1	MYADM	OLIG2	P4HA2	PIK4CB	PRPF31	RAD50	RASGRP2	RBM12	RFX5	RGS11	RHBDF1	RP11-126K1.3	RPL10	RPS9	SERPINB8	SH3BGR	SLC10A3	STAG2	STEAP2	TFEB	TFPT	TRIM22	TRIM5	UBQLN3");
//			while (st.hasMoreTokens())
//				testNodes.add(st.nextElement());
			
			for (int i = 0; i < testNodes.size(); i++) 
				System.out.print(testNodes.elementAt(i)+" ");
			
			testNodes= new Vector(mapUTRspliced.keySet());
			System.out.println(testNodes.size()+"*");
			for (int i = 0; i < testNodes.size(); i++) 
				System.out.print(testNodes.elementAt(i)+" ");
			System.out.println();
			Vector allNodes= new Vector(mapAll.keySet());
//			st= new StringTokenizer(new String("ACCN4 ACSL6 APOA1 APOA4 APOA5 APOC3 ARF5 ARHGAP26 ARHGAP4 ARHGDIG ASCL2 ASZ1 ATP11A ATP5O AVPR2 AXIN1 BCL11B BIRC4 BPIL2 C16orf33 C16orf35 C20orf52 C21orf108 C21orf55 C21orf59 C21orf62 C21orf63 C21orf66 C22orf24 CACNG6 CACNG7 CACNG8 CAPZA2 CATSPER2 CAV1 CAV2 CDC42BPG CDC42EP5 CDH2 CFTR CGN CLDN12 CNOT3 CPNE1 CRAT CRYZL1 CSF2 CTAG1B CTAG2 CTGF CTSD CXorf12 CXorf2 DDX18 DDX43 DEPDC5 DES DKC1 DNASE1L1 DOLPP1 DONSON DRG1 DSCR2 EEF1A1 EHD1 EIF4ENIF1 ELL3 EMD ENPP1 EVX1 EXT1 F10 F7 F8 F8A1 FAM3A FAM50A FBXO7 FGF1 FLNA FOXP2 FOXP4 FRS3 FSCN3 FZD1 G6PD GART GCC1 GDF5 GDF9 GRM8 HBA2 HBB HBD HBE1 HBG1 HBG2 HBQ1 HBZ HCFC1 HMGN1 HNT HOXA1 HOXA10 HOXA11 HOXA13 HOXA2 HOXA3 HOXA4 HOXA5 HOXA6 HOXA7 HOXA9 HS3ST4 HYPK IFNAR1 IFNAR2 IFNGR2 IGF2 IKBKG IL10RB IL13 IL3 IL4 IL5 INHA INS IRAK1 IRF1 ITGB4BP ITSN1 KATNAL1 KCNQ5 KIF3A KIR2DL1 KIR2DL3 KIR2DL4 KIR3DL1 L1CAM LAIR1 LAIR2 LILRA1 LILRA2 LILRA3 LILRB1 LILRB2 LILRB3 LILRB4 LILRB5 LRRK2 LSP1 LUC7L MAP1A MAP3K1 MAP4K2 MCF2L MDFI MECP2 MEN1 MET MFAP1 MMP24 MMP26 MOXD1 MPG MPP1 MRPL23 MRPL28 MTCP1 MTO1 MYADM NCR2 NDUFA3 NFS1 NME4 NR2E1 NRXN2 NUP188 OLIG1 OLIG2 OPN1LW OPN1MW OR51A2 OR51A4 OR51A7 OR51B2 OR51B4 OR51B5 OR51B6 OR51F2 OR51G1 OR51G2 OR51I1 OR51I2 OR51L1 OR51M1 OR51Q1 OR51S1 OR51T1 OR51V1 OR52A1 OR52A5 OR52B6 OR52D1 OR52E2 OR52H1 OR52J3 OR52N4 OR52R1 OR56B1 OSCAR OSTM1 P4HA2 PCDH15 PDLIM4 PFTK1 PGC PHYHD1 PIK4CB PIP5K1A PISD PLXNA3 POGZ POLR3K PPP2R4 PRKCG PRPF31 PSMB4 PSMD4 PYGM RAB11FIP3 RAD50 RBM12 RENBP RFPL2 RFPL3 RFX5 RGS11 RNPC2 RPL10 RPS9 SEC63 SELENBP1 SEPT8 SERPINB10 SERPINB11 SERPINB13 SERPINB2 SERPINB3 SERPINB4 SERPINB7 SERPINB8 SF1 SH3BGR SH3GLB2 SLC10A3 SLC22A11 SLC22A12 SLC22A4 SLC22A5 SLC2A13 SLC4A3 SLC5A1 SLC5A4 SNX27 SNX3 SON SPAG4 SPP2 ST7 STAG2 STEAP2 STRC SYN3 SYNJ1 SYT8 TAZ TES TFEB TH THOC2 TIMP3 TKTL1 TMEM15 TMEM8 TNNI2 TNNT3 TP53BP1 TRIM22 TRIM34 TRIM5 TRIM6 TRPM8 TUFT1 UBQLN3 UGT1A1 UGT1A10 UGT1A3 UGT1A4 UGT1A5 UGT1A6 UGT1A7 UGT1A8 UGT1A9 USP49 WNT2 WRB XX-FW81657B9.4 YWHAH ZNF259").toUpperCase());
//			while (st.hasMoreTokens())
//				allNodes.add(st.nextElement());

//			for (int i = 0; i < allNodes.size(); i++) 
//				System.out.print(allNodes.elementAt(i)+"|");			
			System.out.println(allNodes.size());
			SignificanceTester tester= new SignificanceTester();
			tester.initGO();
			
			System.out.println("=== Hyper, Bonferroni ===");
			tester.significanceTestHyperGeometric(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BONFERRONI);
			System.out.println("=== Hyper, Benjamini-Hochberg ===");
			tester.significanceTestHyperGeometric(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BENJAMINI_HOCHBERG);
			System.out.println("=== Binomial, Bonferroni ===");
			tester.significanceTestBinomial(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BONFERRONI);
			System.out.println("=== Binomial, Benjamini-Hochberg ===");
			tester.significanceTestBinomial(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BENJAMINI_HOCHBERG);
			
			System.out.println("=== UNDER: Hyper, Bonferroni ===");
			tester.significanceTestHyperGeometricUnder(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BONFERRONI);
			System.out.println("=== UNDER: Hyper, Benjamini-Hochberg ===");
			tester.significanceTestHyperGeometricUnder(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BENJAMINI_HOCHBERG);
			System.out.println("=== UNDER: Binomial, Bonferroni ===");
			tester.significanceTestBinomialUnder(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BONFERRONI);
			System.out.println("=== UNDER: Binomial, Benjamini-Hochberg ===");
			tester.significanceTestBinomialUnder(testNodes, allNodes, "0.05", SignificanceTester.SIGNIFICANT_AFTER_CORRECTION, SignificanceTester.CORRECTION_BENJAMINI_HOCHBERG);
		}

	private static void outputDegree(ASVariation[][] filtClasses, PrintStream p) {
		
			// cluster complexity classes
		HashMap h= new HashMap();
		for (int i = 0; i < filtClasses.length; i++) {
			int deg= filtClasses[i][0].getDegree();
			Integer key= new Integer(deg);
			Integer val= (Integer) h.remove(key);
			if (val!= null) 
				val= new Integer(filtClasses[i].length+ val.intValue());
			else 
				val= new Integer(filtClasses[i].length);
			h.put(key, val);
		}
		
			// sort ascending
		Object[] o1= h.values().toArray();
		Object[] o2= h.keySet().toArray();
		Integer[] vals= new Integer[o1.length]; 
		Integer[] keys= new Integer[o2.length];
		for (int i = 0; i < keys.length; i++) {
			vals[i]= (Integer) o1[i];
			keys[i]= (Integer) o2[i];
		}
		// key= complexity
		// val= occurrences
		for (int i = 0; i < keys.length; i++) {
			for (int j = 0; j < keys.length; j++) {
//				if (keys[j].intValue()< keys[i].intValue()) {
				if (vals[j].intValue()< vals[i].intValue()) {
					Integer m= vals[i];
					vals[i]= vals[j];
					vals[j]= m;
					m= keys[i];
					keys[i]= keys[j];
					keys[j]= m;
				}
			}
		}
		
			// output
		for (int i = 0; i < keys.length; i++) {
			p.println(keys[i]+"\t"+vals[i]);
		}
	}

	public static void determineVariationsDetails(Graph g, String methodName, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= (ASVariation[][]) Arrays.filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		int cnt= 0;
		for (int i = 0; i < classes.length; i++) {
			cnt+= classes[i].length;
			for (int j = 0; j < classes[i].length; j++) {
				classes[i][j].outputDetail(p);
			}
		}
	}

	public static ASVariation[][] determineVariations(Graph g, String methodName, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= (ASVariation[][]) Arrays.filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		return outputVariations(classes, false, false, p);
	}

	public static ASVariation[][] determineVariations(Graph g, String methodName) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= (ASVariation[][]) Arrays.filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		return outputVariations(classes, true, false, System.out);
	}

	static String[] readLabels() {
		String fName= "labelSort.txt";
		Vector v= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			while (buffy.ready()) {
				String s= buffy.readLine().trim();
				if (s.length()> 0)
					v.add(s);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return (String[]) Arrays.toField(v);
	}
	
	
	static void writeLabels(ASVariation[][] classes) {
		String fName= "labelSort.txt";
		Vector v= new Vector();
		try {
			BufferedWriter buffy= new BufferedWriter(new FileWriter(fName));
			for (int i = 0; i < classes.length; i++) {
				buffy.write(classes[i][0].toString()+"\n");
			}
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public static ASVariation[][] sortLabels(ASVariation[][] classes) {
		String[] labels= readLabels();
		Vector v= new Vector(classes.length);
		for (int i = 0; i < labels.length; i++) {
			int j;
			for (j = 0; j < classes.length; j++) { 
				if (classes[j][0].toString().equals(labels[i])) {
					v.add(classes[j]);
					break;
				}
			}
			if (j== classes.length) {
				v.add(null);
			}
			
		}
			
		return (ASVariation[][]) Arrays.toField(v);

	}
	
	/**
			 * @param classes
			 * @param sortConsistently
			 * @param writeFile
			 * @param p
			 * @return
			 */
			public static ASVariation[][] outputVariations(ASVariation[][] classes, boolean sortConsistently, boolean writeFile, PrintStream p) {
				
				if (classes== null|| p== null)
					return classes;
				
					// sort
				int nullObj= 0;
				if (sortConsistently) {
					String[] labels= readLabels();
					Vector v= new Vector(classes.length);
					for (int i = 0; i < labels.length; i++) {
						int j;
						for (j = 0; j < classes.length; j++) { 
							if (classes[j]== null)
								continue;
							if (classes[j][0].toString().equals(labels[i])) {
								v.add(classes[j]);
								break;
							}
						}
						if (j== classes.length) {	// empty lines
							v.add(null);
							nullObj++;
						}
					}
					for (int i = 0; i < classes.length; i++) {
						int j;
						for (j = 0; j < labels.length; j++) 
							if (classes[i][0].toString().equals(labels[j]))
								break;
						if (j== labels.length)		// add missing ones
							v.add(classes[i]);
					}
					classes= (ASVariation[][]) Arrays.toField(v);
				}
				
				int threshold= 0;
				int totNum= 0, thrNum= 0;
				for (int i = 0; i < classes.length; i++) { 
					if (classes[i]== null)
						continue;
					totNum+= classes[i].length;
					if (classes[i]!= null&& classes[i].length>= threshold)
						thrNum+= classes[i].length;
				}
				p.println("ASVariations: "+ totNum+ ", structural distinct classes: "+ (classes.length- nullObj));
				int cnt1= 0, cnt2= 0;
				// for excluding low length population
		//		int tot= 0;
		//		for (int i = 0; i < classes.length; i++) {
		//			if (classes[i]== null) 
		//				continue;
		//			int rel= 0;
		//			for (int j = 0; j < classes[i].length; j++) 
		//				if (classes[i][j].getLengthDiff(true)> 10)
		//					tot++;
		//		}
				for (int i = 0; i < classes.length; i++) {
					if (classes[i]== null) {
						p.println();
						continue;
					}
					// exlcuding short length population
		//			int rel= 0;
		//			for (int j = 0; j < classes[i].length; j++) {
		//				if (classes[i][j].getLengthDiff(true)> 10)
		//					rel++;
		//			}
		//			p.println(rel+ " "+ ((float) rel/ tot)+" "+classes[i][0]);
					p.println(classes[i].length+ " "+ ((float) classes[i].length/ totNum)+ " "+classes[i][0]);
						// output all events + coordinates
	//				for (int j = 0; j < classes[i].length; j++) {
	//					SpliceSite[] ss= classes[i][j].getSpliceUniverse();
	//					p.println("\t"+ classes[i][j].getTranscript1().getTranscriptID()+ 
	//											"\t"+ classes[i][j].getTranscript2().getTranscriptID()+
	//											"\t"+ ss[0].getPos()+
	//											"\t"+ ss[ss.length- 1].getPos());
	//				}
		
					if (classes[i].length< threshold) {
						++cnt1;
						cnt2+= classes[i].length;
					}
				}
				//p.println("low represented ("+threshold+") events: "+ cnt1+ " ("+cnt2+")\n");
				p.println("\n\n");
				
				//GraphHandler.writeOut(g);
				if (writeFile)
					writeLabels(classes);
				
				return classes;
			}

	/**
		 * @param classes
		 * @param sortConsistently
		 * @param writeFile
		 * @param p
		 * @return
		 */
		public static ASVariation[][] outputVariationsDegree(ASVariation[][] classes, PrintStream p) {
			
			if (classes== null|| p== null)
				return classes;

				// sort according to degrees
			HashMap mp= new HashMap();
			for (int i = 0; i < classes.length; i++) {
				Integer deg= new Integer(classes[i][0].getDegree());
				Vector v= (Vector) mp.get(deg);
				if (v== null)
					v= new Vector();
				for (int j = 0; j < classes[i].length; j++) 
					v.add(classes[i][j]);
				mp.put(deg, v);
			}
			classes= (ASVariation[][]) Arrays.toField(mp.values().toArray());
			classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);

				// count
			int nullObj= 0;
			int totNum= 0;
			for (int i = 0; i < classes.length; i++) { 
				if (classes[i]== null)
					continue;
				totNum+= classes[i].length;
			}
			p.println("ASVariations: "+ totNum+ ", classes of different degree: "+ (classes.length- nullObj));
			for (int i = 0; i < classes.length; i++) {
				if (classes[i]== null) {
					p.println();
					continue;
				} else
					p.println(classes[i][0].getDegree()+"\t"+classes[i].length);
			}
			p.println("\n\n");
			
			return classes;
		}

	/**
	 * @param classes
	 * @param sortConsistently
	 * @param writeFile
	 * @param p
	 * @return
	 */
	public static void outputVariations(ASMultiVariation[][] classes, PrintStream p) {
		
		if (classes== null|| p== null)
			return;
		
		
		int threshold= 0;
		int totNum= 0, thrNum= 0;
		for (int i = 0; i < classes.length; i++) { 
			if (classes[i]== null)
				continue;
			totNum+= classes[i].length;
			if (classes[i]!= null&& classes[i].length>= threshold)
				thrNum+= classes[i].length;
		}
		p.println("ASVariations: "+ totNum+ ", structural distinct classes: "+ classes.length);
		int cnt1= 0, cnt2= 0;
		for (int i = 0; i < classes.length; i++) {
			if (classes[i]== null) {
				p.println();
				continue;
			}
			p.println(classes[i].length+ " "+ ((float) classes[i].length/ totNum)+ " "+classes[i][0]);

			if (classes[i].length< threshold) {
				++cnt1;
				cnt2+= classes[i].length;
			}
		}
		p.println("\n\n");
		
	}

	public static ASVariation[][] debug(Graph g, int howMany, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;
		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		for (int i = classes.length-1; i > classes.length-(10+1); i--) {
			System.out.println(classes[i][0]);
			
			SpliceSite[] sc= classes[i][0].getSpliceChain1();
			System.out.print(classes[i][0].getGene().getChromosome()+"\t");
			for (int j = 0; j < sc.length; j++) {
				System.out.print(sc[j]+" ");
			}
			System.out.println();
			
			sc= classes[i][0].getSpliceChain2();
			System.out.print(classes[i][0].getGene().getChromosome()+"\t");
			for (int j = 0; j < sc.length; j++) {
				System.out.print(sc[j]+" ");
			}
			System.out.println("\n\n");
		}
		return classes;
	}

	public static ASVariation[][] lowRepresented(Graph g, int howMany, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;
		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		for (int i = 0; i < classes.length; i++) {
			
			if (classes[i].length> howMany)
				continue;
			
			System.out.println("#"+classes[i].length+"\t"+classes[i][0]);
			for (int k = 0; k < classes[i].length; k++) {
				
				SpliceSite[] sc= classes[i][k].getSpliceChain1();
				System.out.print(classes[i][k].getGene().getChromosome()+"\t");
				for (int j = 0; j < sc.length; j++) {
					System.out.print(sc[j]+" ");
				}
				System.out.println();
				
				sc= classes[i][k].getSpliceChain2();
				System.out.print(classes[i][k].getGene().getChromosome()+"\t");
				for (int j = 0; j < sc.length; j++) {
					System.out.print(sc[j]+" ");
				}
				System.out.println();
			}
			System.out.println("\n\n");
		}
		return classes;
	}
	
	static Graph identifyOrthologSpliceSites(Graph g) {
		
		for (int i = 0; i < g.getSpecies().length; i++) {
			Gene[] ge= g.getSpecies()[i].getGenes();
			for (int j = 0; j < ge.length; j++) {
				GeneHomology[] h= ge[j].getHomologies();	// work on homologies
				for (int k = 0; k < h.length; k++) {
					if (!h[k].orthologSpliceSitesIdentified())
						h[k].identifyOrthologSpliceSites();	// only identify if necessary
				}
			}
		}
		
		return g;
	}
	
	static ASVariation[] findClass(Graph g, String s, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;

		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		int totNum= 0;
		for (int i = 0; i < classes.length; i++) {
			if (classes[i][0].toString().equals(s))
				return classes[i];
		}
		return null;
	}
	
	static void filterEvents(File rawFile) {
		try {
		
			BufferedReader buffy= new BufferedReader(new FileReader(rawFile));
			BufferedWriter buffy2= new BufferedWriter(new FileWriter(rawFile.getPath()+ "_filtered"));
			while (buffy.ready()) {
				String s= buffy.readLine();
				if (s.length()< 80)
					buffy2.write(s+ "\n");				
			}
			buffy2.flush();
			buffy2.close();
		} catch (Exception e) {
			System.err.println("Problems:");
			e.printStackTrace();
		}
	}
	
	static ASVariation[][] match(Graph g, String r1, String r2, boolean respectCoding) {
		
		RegExp reg1= new RegExp(r1);
		RegExp reg2= new RegExp(r2);
		Automaton a1= reg1.toAutomaton();
		Automaton a2= reg2.toAutomaton();
		
		ASVariation[][] vars= g.getASVariations(1);
		int hits= 0;
		Vector hitClasses= new Vector();
		for (int i = 0; i < vars.length; i++) {
			a1.init(); a2.init();
			ChessMaster c= new ChessMaster(
					new Automaton[] {a1, a2},
					new SpliceSite[][] {vars[i][0].getSpliceChain1(), vars[i][0].getSpliceChain2()}
			);
			for (int j = 0; j < vars[i].length; j++) {
				if (vars[i][j].getGene().getStableID().equals("ENSG00000186777"))
					System.currentTimeMillis();
			}
			if (c.match()) {
				hits+= vars[i].length;
				hitClasses.add(vars[i]);
			}
		}
		
		return (ASVariation[][]) Arrays.toField(hitClasses);
	}
	
	/**
	 * 
	 * @param vars
	 * @param regions relative positions (1-based) of the splice sites including
	 * the regions for the length plot. Negative values are counted from the end. 
	 */
	public static void lengthPlot(Graph g, boolean respectCoding) {
		
		ASVariation[][] vars= g.getASVariations(1);
		int base= 0;
		for (int i = 0; i < vars.length; i++) 
			base+= vars[i].length;
		long t0= System.currentTimeMillis();
		ASVariation[][] hits= null;
		int filterVal= -1;
		
			// donors
		hits= match(g, "^D1$", "^D2(3,4)*$", true);
		//hits= filter5UTR(hits);
		int cnt= countHits(hits);
		System.out.println("\nDonor type 1 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^D2$", "^D1(3,4)*$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDonor type 2 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^,A,4$", "^,A,(1,2)*3$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nAcceptor type 1 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^,A,3$", "^,A,(1,2)*4$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nAcceptor type 2 ("+cnt+")");
		outputLength(hits);

		hits= match(g, "^,D,1,6$", "^,D,2(3,4)*5$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDon/Acc type 1 ("+cnt+")");
		outputLength(hits);

		hits= match(g, "^,D,2,5$", "^,D,1(3,4)*6$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDon/Acc type 1 ("+cnt+")");
		outputLength(hits);
	}
	private static void outputLength(ASVariation[][] hits) {
		for (int i = 0; i < hits.length; i++) {
			for (int j = 0; j < hits[i].length; j++) {
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= (sUniverse[k+1].getPos()- sUniverse[k].getPos())+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)% 3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
			}
		}
	}
	
	public static void outputSSOut(Graph g, PrintStream pr) {
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_NONE);		
		System.out.println(g.countGenesTranscriptsExons());
		SpliceSite[][] res= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
		for (int j = 0; j < res[0].length; j++) {
			pr.println(res[0][j].getGene().getChromosome()+ " "
					+ Math.abs(res[0][j].getPos())+ " "+ res[0][j].getGene().getStrand()+ " "
					+ (res[0][j].isConstitutive()?"const":"alt")+ " "+ (res[0][j].isDonor()?"don":"acc"));
		}
		for (int j = 0; j < res[1].length; j++) {
			pr.println(res[1][j].getGene().getChromosome()+ " "
					+ Math.abs(res[1][j].getPos())+ " "+ res[1][j].getGene().getStrand()+ " "
					+ (res[1][j].isConstitutive()?"const":"alt")+ " "+ (res[1][j].isDonor()?"don":"acc"));
		}
	}

	public static void outputSSOutCdsUtr(Graph g, PrintStream pr) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		
			// map for x-check
		SpliceSite[][] ss= g.getSpliceSites(Gene.REGION_COMPLETE_GENE);
		HashMap map= new HashMap(ss[0].length+ ss[1].length);
		for (int i = 0; i < ss[0].length; i++) 
			map.put(ss[0][i], ss[0][i]);
		for (int i = 0; i < ss[1].length; i++) 
			map.put(ss[1][i], ss[1][i]);
		
		int[] types= new int[] {Gene.REGION_REAL_5UTR, Gene.REGION_REAL_CDS};
		for (int i = 0; i < types.length; i++) {
			SpliceSite[][] res= g.getSpliceSites(types[i]);
			for (int j = 0; j < res[0].length; j++) {
				pr.println(res[0][j].getGene().getChromosome()+ " "
						+ Math.abs(res[0][j].getPos())+ " "+ res[0][j].getGene().getStrand()+ " "
						+ (res[0][j].isConstitutive()?"const":"alt")+ " "+ (res[0][j].isDonor()?"don":"acc")+ " "+
						(i== 0?"utr":"cds"));
				Object o= map.remove(res[0][j]);
				if (o== null)
					System.err.println(res[0][j]+ " not found");
			}
			for (int j = 0; j < res[1].length; j++) {
				if (res[1][j].isConstitutive())
					pr.println(res[1][j].getGene().getChromosome()+ " "
							+ Math.abs(res[1][j].getPos())+ " "+ res[1][j].getGene().getStrand()+ " "
							+ (res[1][j].isConstitutive()?"const":"alt")+ " "+ (res[1][j].isDonor()?"don":"acc")+ " "+
							(i== 0?"utr":"cds"));
				Object o= map.remove(res[1][j]);
				if (o== null)
					System.err.println(res[1][j]+ " not found");
			}
		}
		
		Object[] oo= map.values().toArray();
		for (int i = 0; i < oo.length; i++) {
			SpliceSite s= (SpliceSite) oo[i];
			pr.println(s.getGene().getChromosome()+ " "
					+ Math.abs(s.getPos())+ " "+ s.getGene().getStrand()+ " "
					+ (s.isConstitutive()?"const":"alt")+ " "+ (s.isDonor()?"don":"acc")+ " "+
					"twilight");
		}
	}
	private static int countHits(ASVariation[][] vars) {
		int sum= 0;
		for (int i = 0; i < vars.length; i++) 
			sum+= vars[i].length;
		return sum;
	}
	
	/**
	 * 
	 * @param vars
	 * @param regions relative positions (1-based) of the splice sites including
	 * the regions for the length plot. Negative values are counted from the end. 
	 */
	static void lengthPlot_old(Graph g) {
		
		ASVariation[][] vars= g.getASVariations(1);
		int base= 0;
		for (int i = 0; i < vars.length; i++) 
			base+= vars[i].length;
		long t0= System.currentTimeMillis();
		ASVariation[][] hits= null;
		int filterVal= -1;
		
			// donors
		hits= match(g, "^D1$", "^D2(3,4)*$", true);
		int hitcount= 0;
		int filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("Donor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		int nul= 0, eins= 0, zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= (sUniverse[k+1].getPos()- sUniverse[k].getPos())+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)% 3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");
		
		hits= match(g, "^D2$", "^D1(3,4)*$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		int eq= 0, diff11= 0, diff12= 0, diff21= 0, diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+"\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
		
			// acceptors
		hits= match(g, "^,A,4$", "^,A,(1,2)*3$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nAcceptor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		nul= 0; eins= 0; zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos();
				int skip= 0;
				for (int k = 0; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos() + 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getPredSpliceSite(sUniverse[0]);
				if (s== null)
					s= hits[i][j].getTranscript2().getPredSpliceSite(sUniverse[0]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)%3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+ alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");

		hits= match(g, "^,A,3$", "^,A,(1,2)*4$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nAcceptor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		eq= 0; diff11= 0; diff12= 0; diff21= 0; diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos();
				int skip= 0;
				for (int k = 0; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getPredSpliceSite(sUniverse[0]);
				if (s== null)
					s= hits[i][j].getTranscript2().getPredSpliceSite(sUniverse[0]);
				int b= (sUniverse[0].getPos()- s.getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+ "\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
		
		hits= match(g, "^,D,1,6$", "^,D,2(3,4)*5$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor/Acceptor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		nul= 0; eins= 0; zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= (sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos())+
						(sUniverse[1].getPos()- sUniverse[0].getPos());
				
				int skip= 0;	
				for (int k = 2; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				
				int a= (sUniverse[2].getPos()- sUniverse[1].getPos())+ 1;
				int b= (sUniverse[sUniverse.length- 2].getPos()- sUniverse[sUniverse.length- 3].getPos())+ 1;
				int x= (alt+skip)%3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+"\t"+
						+x+"\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");

		
		hits= match(g, "^,D,2,5$", "^,D,1(3,4)*6$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor/Acceptor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		eq= 0; diff11= 0; diff12= 0; diff21= 0; diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= (sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos())+
						(sUniverse[1].getPos()- sUniverse[0].getPos());
				
				int skip= 0;	
				for (int k = 2; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				
				int a= (sUniverse[2].getPos()- sUniverse[1].getPos())+ 1;
				int b= (sUniverse[sUniverse.length- 2].getPos()- sUniverse[sUniverse.length- 3].getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+"\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
	}
	
	public static void getSpliceSites(Graph g, String fName, String source) {
		g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);	// init SSs
		Gene[] ge= g.getGenes();
		Vector gffV= new Vector(ge.length* 10);
		SpliceSite[] ss;
		String[] attributes;
		for (int i = 0; i < ge.length; i++) {
			ss= ge[i].getSpliceSites(SpliceSite.CONSTITUTIVE_SS, Gene.REGION_REAL_5UTR);
			attributes= new String[] {"CONSTITUTIVE_SS", "REGION_REAL_5UTR"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));
			
			ss= ge[i].getSpliceSites(SpliceSite.ALTERNATE_SS, Gene.REGION_REAL_5UTR);
			attributes= new String[] {"ALTERNATE_SS", "REGION_REAL_5UTR"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));

			ss= ge[i].getSpliceSites(SpliceSite.CONSTITUTIVE_SS, Gene.REGION_REAL_CDS);
			attributes= new String[] {"CONSTITUTIVE_SS", "REGION_REAL_CDS"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));

			ss= ge[i].getSpliceSites(SpliceSite.ALTERNATE_SS, Gene.REGION_REAL_CDS);
			attributes= new String[] {"ALTERNATE_SS", "REGION_REAL_CDS"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));
			
			ss= ge[i].getSpliceSites(SpliceSite.CONSTITUTIVE_SS, Gene.REGION_REAL_3UTR);
			attributes= new String[] {"CONSTITUTIVE_SS", "REGION_REAL_3UTR"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));
			
			ss= ge[i].getSpliceSites(SpliceSite.ALTERNATE_SS, Gene.REGION_REAL_3UTR);
			attributes= new String[] {"ALTERNATE_SS", "REGION_REAL_3UTR"};
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				gffV.add(GTFObject.createGFFObject(ss[j], source, attributes));
		}
		
		GTFWrapper wrapper= new GTFWrapper(fName);
		wrapper.setGtfObj((GTFObject[]) Arrays.toField(gffV));
		try {
			wrapper.writeGFF();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void getUTRExons(Graph g, String fName) {
		Exon[] ex5= g.getExons(Gene.REGION_REAL_5UTR);
		Exon[] ex3= g.getExons(Gene.REGION_REAL_3UTR);
		
		Vector v= new Vector(ex5.length+ ex3.length);
		for (int i = 0; i < ex5.length; i++) {
			GTFObject gtf= new GTFObject();
			try {
				gtf.setSeqname(ex5[i].getChromosome());
				gtf.setSource("REAL_5UTR");
				gtf.setFeature("exon");
				gtf.setStart(Math.abs(ex5[i].getStart()));
				gtf.setEnd(Math.abs(ex5[i].getEnd()));
				gtf.setStrand(ex5[i].getStrand());
				Transcript[] tr= ex5[i].getTranscripts();
				for (int j = 0; j < tr.length; j++)
					gtf.addAttribute("trascript_id", tr[j].getTranscriptID());
				v.add(gtf);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		for (int i = 0; i < ex3.length; i++) {
			try {
				GTFObject gtf= new GTFObject();
				gtf.setSeqname(ex3[i].getChromosome());
				gtf.setSource("REAL_5UTR");
				gtf.setFeature("exon");
				gtf.setStart(Math.abs(ex3[i].getStart()));
				gtf.setEnd(Math.abs(ex3[i].getEnd()));
				gtf.setStrand(ex3[i].getStrand());
				Transcript[] tr= ex3[i].getTranscripts();
				for (int j = 0; j < tr.length; j++)
					gtf.addAttribute("trascript_id", tr[j].getTranscriptID());
				v.add(gtf);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		GTFObject[] obj= (GTFObject[]) Arrays.toField(v);
		GTFWrapper wrapper= new GTFWrapper(fName);
		wrapper.setGtfObj(obj);
		try {
			wrapper.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void outputMultiVars(ASMultiVariation[] mVars, PrintStream p) {
		for (int j = 0; j < mVars.length; j++) {
			HashMap map= mVars[j].getTransHash();
			Object[] sc= map.keySet().toArray();
			for (int k = 0; k < sc.length; k++) {
				SpliceSite[] sc1= (SpliceSite[]) sc[k];
				for (int m = 0; m < sc1.length; m++) 
					p.print("\t"+ sc1[m]);
				p.print("\t[");
				Transcript[] trpt= (Transcript[]) map.get(sc1);
				for (int m = 0; m < trpt.length; m++) 
					System.out.print(trpt[m]+",");
				p.println("]");
			}
			//analyzeMutex(mVars[j], System.out);			
			p.println("\n");
		}

	}
	
	public static void analyzeMutex(ASMultiVariation var, PrintStream p) {
			// isolate all sc with 2 SSs
		SpliceSite[][] sc= var.getSpliceChains();
		Vector v= new Vector();
		DirectedRegion reg;
		for (int i = 0; i < sc.length; i++) {
			if (sc[i]== null|| sc[i].length!= 2)
				continue;
			Transcript t[]= (Transcript[]) var.getTransHash().get(sc[i]);
			int j;
			for (j = 0; j < t.length; j++) {
				Translation trln= t[j].getTranslations()[0];
				if (trln.contains(sc[i][0])|| trln.contains(sc[i][1])) {
					p.print("cod, ");
					break;
				}
			}
			if (j== t.length)
				p.print("ncod, ");
			
			reg= new DirectedRegion(sc[i][0].getPos(), sc[i][1].getPos(), 
					var.getGene().getStrand());
			reg.setChromosome(var.getGene().getChromosome());
			reg.setSpecies(var.getGene().getSpecies());
			v.add(reg);
		}
		p.println();

			// length
		p.println("lengthes:");
		DirectedRegion[] regs= (DirectedRegion[]) Arrays.toField(v);
		for (int i = 0; i < regs.length; i++) {
			for (int j = (i+1); j < regs.length; j++) {
				p.println(regs[i].getLength()+" ("+(regs[i].getLength()%3)+"), "+
						regs[j].getLength()+" ("+(regs[j].getLength()%3)+"), "+
						(regs[i].getLength()+ regs[j].getLength())+ " ("+ ((regs[i].getLength()+ regs[j].getLength())% 3)+ ")"
						);
			}
		}
				
			// ntID
		for (int i = 0; i < regs.length; i++) {
			for (int j = (i+1); j < regs.length; j++) {
				String seq1= Graph.readSequence(regs[i]);
				String seq2= Graph.readSequence(regs[j]);
				int id= Toolbox.seqIdentity(seq1, seq2);
				int percID= id* 100/ Math.min(seq1.length(), seq2.length());
				p.print(percID+ ", ");
			}
		}
	}
	
	public static void main(String[] args) {

		//getSpliceSites(g, "sylvain_spicesites.gff", "gencode");
		//output5UTRExonFragments();
		//outputSpliceSites();
		//output5UTRIntrons();
		//outputCDSIntrons();
//		if (1== 1)
//			System.exit(0);
		
		
		Graph g= null;
		//g= getGraph(INPUT_ENCODE);
		//g.filterNMDTranscripts();
		//test02_ss_statistics(g, System.out, false);
		//test04_determineVariations_nmd();
		//test02_ss_statistics();
		//test01_clusters_coverage_as();
		//mvar_output();
		//test05_GO_splice_utr();
		//test02_UTRvsCDS_go();
		test04_determineVariations_nmd();
		//test02_ss_statistics();
		//test05_GO_splice_utr_sigTest();
		if (1== 1)
			System.exit(0);
		
		
		g.filterNonCodingTranscripts();
		//g.filterCodingTranscripts();
		//test04a_determineVarDegree(g, System.out);
		//test04_checkOutsideEncode(g, System.out);

//		ASMultiVariation[][] vars= g.getASMultiVariations(-1);	
//		ASMultiVariation[][] vars= g.getASMultiVariations(2);	// 2
		ASMultiVariation[][] vars= g.getASMultiVariations(2);	// 2
		vars= (ASMultiVariation[][]) Arrays.sort2DFieldRev(vars);
		for (int i = 0; i < vars.length; i++) {
			System.out.println(vars[i].length+ "\t"+ vars[i][0]);
//			if (vars[i][0].toString().startsWith("(1-2^ // 3-4^)")||
//					vars[i][0].toString().startsWith("(1-2^3-4^ // 1-2^ // 3-4^)")||
//					vars[i][0].toString().startsWith("(1-2^ // 3-4^ // )")) {
			if ((vars[i][0].toString().startsWith("(1^ // 2^ // 3^)"))||
					(vars[i][0].toString().startsWith("(1- // 2- // 3-)"))){
				outputMultiVars(vars[i], System.out);
			}
		}
		
//		ASVariation[][] vars2= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
//		vars2= (ASVariation[][]) Arrays.sort2DFieldRev(vars2);
//		for (int i = 0; i < vars2.length; i++) {
//			System.out.println(vars2[i].length+ "\t"+ vars2[i][0]);
//			if (vars2[i][0].toString().equals("(1=2^ // 3=4^)")) {
//				System.out.println(vars2[i].length+ "\t"+ vars2[i][0]);
//				for (int j = 0; j < vars2[i].length; j++) {
//					System.out.print("\t"+ vars2[i][j].getGene().getChromosome()+" ");
//					SpliceSite[] sc1= vars2[i][j].getSpliceChain1();
//					SpliceSite[] sc2= vars2[i][j].getSpliceChain2();
//					SpliceSite[][] sc= new SpliceSite[][] {sc1, sc2};
//					if (sc2.length> sc1.length)
//						Arrays.swap(sc);
//					for (int k = 0; k < sc.length; k++) {
//						for (int m = 0; m < sc[k].length; m++) 
//							System.out.print(sc[k][m]+ " ");
//						System.out.print("\t// ");
//					}
//					System.out.println();
//				}
//			}
//		}
		
		
		
//		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
//		for (int i = 0; i < vars.length; i++) {
//			if (vars[i][0].toString().equals("(2^ // 1^3=4^)"))
//				for (int j = 0; j < vars[i].length; j++) {
//					System.out.println(vars[i][j].getGene().getGeneID()+": "+ 
//							vars[i][j].getTranscript1().getTranscriptID()+
//							vars[i][j].getTranscript2().getTranscriptID());
//				}
//		}
		
		//		filterEvents(new File("human_graph_as_variations_non-redundant.log"));

		//outputCheckFile();
		//outputCodingRegions();
		//output3UTR();
		
//		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
//		//Graph g= null;
//		Graph g= adaptor.getGraphAllGenes(new Species("chimp"));
//		
//		lengthPlot(g);
//		//determineVariations(g);
//		if (1== 1)
//			System.exit(0);
//		
//		determineVariations(g);
//		
//		ASVariation[] var= findClass(g, "(1^3= // 2^4=)");
//
//		for (int i = 0;  var!= null&& i < var.length; i++) {
//			Gene ge= var[i].getGene();
//			SpliceOSigner painter= new SpliceOSigner(ge);
//			JFrame frame= new JFrame();
//			frame.getContentPane().add(painter);
//			frame.pack();
//			frame.setVisible(true);
//		}
//		
//		if (1== 1)
//			System.exit(0);
//		
//		ASVariation[][][] clazzes= g.getASEvents();
//		try {
//			Paparazzi.writeGIF(null, new File("test.gif"));
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		if (1== 1)
//			System.exit(0);
//		clazzes= (ASVariation[][][]) Arrays.sort2DField(clazzes);
//		System.out.println("asSuperClasses "+ clazzes.length);		
//		for (int i = 0; i < clazzes.length; i++) {
//			System.out.println(clazzes[i].length+ " "+clazzes[i][0][0].toASEventString());
//
//			for (int j = 0; j < clazzes[i].length; j++) 
//				System.out.println("\t"+ clazzes[i][j].length+ " "+clazzes[i][j][0]);
//		}
//
//		//GraphHandler.writeOut(g);
//		
//		
//		
//		
//			//System.out.print(ge[i].getStableID()+"\t");
////			ASMultiVariation[] as= ge[i].getASComplexes();
////			boolean breakk= false;
////			for (int j = 0; as!= null&& !breakk&& j < as.length; j++) 
////				for (int k = 0; !breakk&& k < as[j].getASVariations().length; k++) 
////					for (int l = 0; !breakk&& l < as[j].getASVariations()[k].getASEvents().length; l++) 
////						if (as[j].getASVariations()[k].getASEvents()[l].isIntronRetention())
////							breakk= true;
////			if (ge[i].getStableID().equals("ENSG00000107643")) {
////				JFrame myFrame= new JFrame(); 
////				myFrame.getContentPane().add(new TranscriptPainter(ge[i]));
////				myFrame.pack();
////				myFrame.setVisible(true);
////			}				
////		}
////		System.out.println(asevents+ "AS Variation");
//		
	}

	public static void mvar_output() {
	
			Graph g= null;
			g= getGraph(INPUT_ENCODE);
			g.filterNonCodingTranscripts();
			//g.filterCodingTranscripts();
			
	//		ASMultiVariation[][] vars= g.getASMultiVariations(-1);	
			ASMultiVariation[][] vars= g.getASMultiVariations(-1);	// 2
	//		ASMultiVariation[][] vars= g.getASMultiVariations2();	// 2
			vars= (ASMultiVariation[][]) Arrays.sort2DFieldRev(vars);
			for (int i = 0; i < vars.length; i++) {
				System.out.println(vars[i].length+ "\t"+ vars[i][0]);
	//			if (vars[i][0].toString().startsWith("(1-2^ // 3-4^)")||
	//					vars[i][0].toString().startsWith("(1-2^3-4^ // 1-2^ // 3-4^)")||
	//					vars[i][0].toString().startsWith("(1-2^ // 3-4^ // )")) {
//				if ((vars[i][0].toString().startsWith("(1^ // 2^ // 3^)"))||
//						(vars[i][0].toString().startsWith("(1- // 2- // 3-)"))){
//					outputMultiVars(vars[i], System.out);
//				}
			}
			
			ASVariation[][] vars2= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			vars2= (ASVariation[][]) Arrays.sort2DFieldRev(vars2);
			for (int i = 0; i < vars2.length; i++) {
				System.out.println(vars2[i].length+ "\t"+ vars2[i][0]);
				if (vars2[i][0].toString().equals("(1=2^ // 3=4^)")) {
					System.out.println(vars2[i].length+ "\t"+ vars2[i][0]);
					for (int j = 0; j < vars2[i].length; j++) {
						System.out.print("\t"+ vars2[i][j].getGene().getChromosome()+" ");
						SpliceSite[] sc1= vars2[i][j].getSpliceChain1();
						SpliceSite[] sc2= vars2[i][j].getSpliceChain2();
						SpliceSite[][] sc= new SpliceSite[][] {sc1, sc2};
						if (sc2.length> sc1.length)
							Arrays.swap(sc);
						for (int k = 0; k < sc.length; k++) {
							for (int m = 0; m < sc[k].length; m++) 
								System.out.print(sc[k][m]+ " ");
							System.out.print("\t// ");
						}
						System.out.println();
					}
				}
			}
			
			
			
	//		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
	//		for (int i = 0; i < vars.length; i++) {
	//			if (vars[i][0].toString().equals("(2^ // 1^3=4^)"))
	//				for (int j = 0; j < vars[i].length; j++) {
	//					System.out.println(vars[i][j].getGene().getGeneID()+": "+ 
	//							vars[i][j].getTranscript1().getTranscriptID()+
	//							vars[i][j].getTranscript2().getTranscriptID());
	//				}
	//		}
			
			//		filterEvents(new File("human_graph_as_variations_non-redundant.log"));
	
			//outputCheckFile();
			//outputCodingRegions();
			//output3UTR();
			
	//		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
	//		//Graph g= null;
	//		Graph g= adaptor.getGraphAllGenes(new Species("chimp"));
	//		
	//		lengthPlot(g);
	//		//determineVariations(g);
	//		if (1== 1)
	//			System.exit(0);
	//		
	//		determineVariations(g);
	//		
	//		ASVariation[] var= findClass(g, "(1^3= // 2^4=)");
	//
	//		for (int i = 0;  var!= null&& i < var.length; i++) {
	//			Gene ge= var[i].getGene();
	//			SpliceOSigner painter= new SpliceOSigner(ge);
	//			JFrame frame= new JFrame();
	//			frame.getContentPane().add(painter);
	//			frame.pack();
	//			frame.setVisible(true);
	//		}
	//		
	//		if (1== 1)
	//			System.exit(0);
	//		
	//		ASVariation[][][] clazzes= g.getASEvents();
	//		try {
	//			Paparazzi.writeGIF(null, new File("test.gif"));
	//		} catch (IOException e) {
	//			// TODO Auto-generated catch block
	//			e.printStackTrace();
	//		}
	//		if (1== 1)
	//			System.exit(0);
	//		clazzes= (ASVariation[][][]) Arrays.sort2DField(clazzes);
	//		System.out.println("asSuperClasses "+ clazzes.length);		
	//		for (int i = 0; i < clazzes.length; i++) {
	//			System.out.println(clazzes[i].length+ " "+clazzes[i][0][0].toASEventString());
	//
	//			for (int j = 0; j < clazzes[i].length; j++) 
	//				System.out.println("\t"+ clazzes[i][j].length+ " "+clazzes[i][j][0]);
	//		}
	//
	//		//GraphHandler.writeOut(g);
	//		
	//		
	//		
	//		
	//			//System.out.print(ge[i].getStableID()+"\t");
	////			ASMultiVariation[] as= ge[i].getASComplexes();
	////			boolean breakk= false;
	////			for (int j = 0; as!= null&& !breakk&& j < as.length; j++) 
	////				for (int k = 0; !breakk&& k < as[j].getASVariations().length; k++) 
	////					for (int l = 0; !breakk&& l < as[j].getASVariations()[k].getASEvents().length; l++) 
	////						if (as[j].getASVariations()[k].getASEvents()[l].isIntronRetention())
	////							breakk= true;
	////			if (ge[i].getStableID().equals("ENSG00000107643")) {
	////				JFrame myFrame= new JFrame(); 
	////				myFrame.getContentPane().add(new TranscriptPainter(ge[i]));
	////				myFrame.pack();
	////				myFrame.setVisible(true);
	////			}				
	////		}
	////		System.out.println(asevents+ "AS Variation");
	//		
		}
}
