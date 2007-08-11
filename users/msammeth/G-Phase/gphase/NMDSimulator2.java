package gphase;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Arrays;

public class NMDSimulator2 {
	public static int MIN_DIST_NC_EJC= 50;
	public static int MAX_DIST_NC_EJC= 55;
	public static int MIN_ORF_LENGTH_AA= 35;

	public static void checkNMD() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		g.filterNonCodingTranscripts();
		
		NMDSimulator2 sim;
		int cntNMD= 0, cntAll= 0;
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				sim= new NMDSimulator2(t[j]);
				if (sim.isNMD())
					++cntNMD;
				else
					++cntAll;
			}
		}
	}
	
	public static void testHAVANAAnnotation() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		g.filterNonCodingTranscripts();
		
		int cntWrong= 0, cntCorr= 0;
		HashMap wrongHash= new HashMap();
		Gene[] ge= g.getGenes();
		Translation predORF, annORF;
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();			
			for (int j = 0; j < t.length; j++) {
				predORF= t[j].findLongestORF();
				annORF= t[j].getTranslations()[0];
				if (predORF.getStart()!= annORF.getStart()||
						predORF.getEnd()!= annORF.getEnd()) {
					++cntWrong;
					wrongHash.put(annORF, predORF);
				} else
					++cntCorr;
			}
		}
		
		System.out.println("correct "+cntCorr+", wrong "+cntWrong);
		String fName= "tstHAVANAAnnotation.txt";
		fName= Toolbox.checkFileExists(fName);
		PrintStream p= null;
		try {
			p= new PrintStream(fName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		Object[] wrongKeys= wrongHash.keySet().toArray();
		for (int i = 0; i < wrongKeys.length; i++) {
			annORF= (Translation) wrongKeys[i];
			predORF= (Translation) wrongHash.get(wrongKeys[i]);
			p.println(">"+ annORF.getTranscript());
			p.println(annORF.getStart()+ " \t"+ annORF.getEnd());
			p.println(predORF.getStart()+ " \t"+ predORF.getEnd());
			p.println();
		}
		p.flush(); p.close();
	}

	public static void outputNMDEvents() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		

			// mark transcripts
		int cntFP= 0, cntFN= 0;
		Gene[] ge= g.getGenes();
		Translation predORF;
		NMDSimulator2 sim;
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();			
			for (int j = 0; j < t.length; j++) {
				predORF= t[j].findLongestORF();
				sim= new NMDSimulator2(t[j]);
				if (predORF== null) {
					if (t[j].isCoding())
						++cntFP;
					t[j].setNmd((byte) 0);
					continue;
				}
				t[j].setPredORF(predORF);
				sim= new NMDSimulator2(t[j]);
				if (sim.isTerminatingUpstreamOfEJC(predORF, MIN_DIST_NC_EJC))
					t[j].setNmd((byte) 3);
				else if (sim.hasUsORF(predORF, MIN_ORF_LENGTH_AA))
					t[j].setNmd((byte) 5);
				if ((!t[j].isCoding())&& (t[j].getNmd()> 0))
					++cntFN;
			}
		}
		
			// extract events
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_NONE);
		vars= (ASVariation[][]) Arrays.sort2DFieldRev(vars);
		
		
			// output 
		System.out.println("HAVANA protein coding FP "+cntFP+", FN"+cntFN);
		String fName= "outEventsNMDSimulator.txt";
		fName= Toolbox.checkFileExists(fName);
		PrintStream p= null;
		try {
			p= new PrintStream(fName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		p.println("AS Event // coding // trpt1 ID // annORF start // end // predORF start // end // NMD(3'/5'/no, hierachically) // schain ///" +
				"trpt2 ID // annORF start // end // predORF start // end // NMD(3'/5'/no, hierachically) // schain");
		Comparator compi= new ASVariation.SpliceSiteOrderComparator();
		for (int i = 0; i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				String s= vars[i][j].toString();
				String codS= "";
				if (vars[i][j].isTouching5UTR())
					codS+= "5UTR-";
				if (vars[i][j].isProteinCoding())
					codS+= "cod-";
				if (vars[i][j].isTouching3UTR())
					codS+= "3UTR-";
				if (codS.length()> 0)
					codS= codS.substring(0, codS.length()- 1);
				else
					codS= ".";	// for 2 nc transcripts
				s+= "\t"+ codS;
					
				Transcript t0= vars[i][j].getTranscript1();
				Transcript t1= vars[i][j].getTranscript2();
				SpliceSite[] sc0= vars[i][j].getSpliceChain1(); 
				SpliceSite[] sc1= vars[i][j].getSpliceChain2();
				if (compi.compare(sc0, sc1)> 0) {
					SpliceSite[] sch= sc0; sc0= sc1; sc1= sch;
					Transcript h= t0; t0= t1; t1= h;
				}
				s+= "\t"+ t0.toStringNMD()+ "\t"+ Transcript.toStringSChain(sc0);
				s+= "\t"+ t1.toStringNMD()+ "\t"+ Transcript.toStringSChain(sc1);
				p.println(s);
			}
		}
		p.flush(); p.close();
	}

	public static void testNMD() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		g.filterNonCodingTranscripts();
		// g.getASVariations(ASMultiVariation.FILTER_NONE);		// for AS flags..
		
		int cntWrong= 0, cntCorr= 0, cntWrongNMD= 0, cntCorrNMD= 0, cntFilt1= 0, cntFilt2= 0;
		HashMap wrongHash= new HashMap();
		HashMap nmdHash= new HashMap();
		Gene[] ge= g.getGenes();
		Translation predORF, annORF;
		NMDSimulator2 sim;
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();			
			for (int j = 0; j < t.length; j++) {
				annORF= t[j].getTranslations()[0];
				if (annORF.get5PrimeEdge()== t[j].get5PrimeEdge()||
						annORF.get3PrimeEdge()== t[j].get3PrimeEdge()) {
					//System.out.println("skipping (open ended HAVANA frame) "+ t[j]);
					++cntFilt1;
					continue;
				}
				if (annORF.getSplicedLength()< NMDSimulator2.MIN_ORF_LENGTH_AA* 3) {
					//System.out.println("skipping (HAVANA ORF too small) "+ t[j]);
					++cntFilt2;
					continue;
				}
				predORF= t[j].findHavanaORF(); 	// findLongestORF();
				sim= new NMDSimulator2(t[j]);
				if (predORF== null) {
					++cntWrong;
					wrongHash.put(annORF, predORF);
					continue;
				}
				if (predORF.getStart()!= annORF.getStart()||
						predORF.getEnd()!= annORF.getEnd()) {
					++cntWrong;
					if (sim.isNMD(predORF))
						++cntWrongNMD;
					else
						wrongHash.put(annORF, predORF);
				} else {
					++cntCorr;
					if (sim.isNMD(annORF))
						++cntCorrNMD;
				}
			}
		}
		
			// output wrong list
		System.out.println("skipped "+(cntFilt1+ cntFilt2)+ " trpts ("+cntFilt1+" open ended, "+cntFilt2+" too short)");
		Object[] wrongKeys= wrongHash.keySet().toArray();
		System.out.println("correct "+cntCorr+"("+cntCorrNMD+
				"), wrong "+cntWrong+"("+cntWrongNMD+")");
		String fName= "tstHAVANAAnnotation.txt";
		fName= Toolbox.checkFileExists(fName);
		PrintStream p= null;
		try {
			p= new PrintStream(fName);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < wrongKeys.length; i++) {
			annORF= (Translation) wrongKeys[i];
			predORF= (Translation) wrongHash.get(wrongKeys[i]);
			p.println(">"+ annORF.getTranscript());
			p.println(annORF.getStart()+ " \t"+ annORF.getEnd());
			if (predORF!= null)
				p.println(predORF.getStart()+ " \t"+ predORF.getEnd());			
			p.println();
		}
		p.flush(); p.close();
	}

	public static void testNMD_noncoding() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		
		int cntNMD3= 0, cntNMD0= 0, cntNoORF= 0;
		HashMap wrongHash= new HashMap();
		Gene[] ge= g.getGenes();
		Translation predORF, annORF;
		NMDSimulator2 sim;
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				if (t[j].isCoding())
					continue;
				predORF= t[j].findLongestORF();
				sim= new NMDSimulator2(t[j]);
				if (predORF== null) { 
					++cntNoORF;
					continue;
				}
				if (sim.isNMD(predORF))
					++cntNMD3;
				else
					++cntNMD0;
			}
		}
		
			// output wrong list
		int total= (cntNMD3+ cntNMD0+ cntNoORF);
		int total2= (cntNMD3+ cntNMD0);
		System.out.println("total "+total+": no ORF "+cntNoORF+ "("+ ((float) cntNoORF/ total)+
				"), NMD3 "+cntNMD3+"("+((float) cntNMD3/ total2)+"), NMD0 "+cntNMD0+"("+((float) cntNMD0/ total2)+")");
	}
	
	public static void main(String[] args) {
		//testNMD();
		//testNMD_noncoding();
		outputNMDEvents();
	}
	
	Transcript trpt= null;
	
	public NMDSimulator2(Transcript t) {
		this.trpt= t;
	}
	
	/**
	 * The stop codon must be in the last exon or no further than 50bp`
	 * from the end of the penultimate exon. [HAVANA]
	 * 
	 * @param minDistNts
	 * @return
	 */
	public boolean isTerminatingUpstreamOfEJC(Translation tln, int minDistNt) {
		
		int critPoint= tln.get3PrimeEdge(); 
		int[] prematStops= tln.getPrematureStops();	// genomic positions
		if (prematStops.length> 0&& prematStops[0]< critPoint)
			critPoint= prematStops[0];
		critPoint= trpt.getGenomicPosition(trpt.getExonicPosition(critPoint)+ minDistNt);
		
		int i;
		for (i = 0; i < trpt.getExons().length; i++) 
			if (trpt.getExons()[i].get5PrimeEdge()> critPoint)	// compare 5' !!
				break;
				
		if (i< trpt.getExons().length)
			return true;
		return false;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS > 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans
	 * @param maxDistAA
	 * @deprecated HAVANA seems to mean something different by this phrase,
	 * see <code>hasUsORF</code>
	 * @return
	 */
	public boolean isInternalATG(Translation trans, int maxDistAA) {
		int maxDistNt= maxDistAA* 3;
		int startPos= trpt.getExonicPosition(trans.get5PrimeEdge());
		if (startPos< maxDistNt)
			return false;
		int frame= startPos% 3;
		String seq= trpt.getSplicedSequence();
		seq= seq.substring(frame);	// to search in the correct frame
		int limit= startPos- frame- maxDistNt;
		int i;	// find last inframe start before ATG (outside the safe distance)
		int usATGPos= 0;
		for (i = 0; (i+3) < limit; i+=3) 
			if (seq.substring(i, i+3).equalsIgnoreCase(Translation.START_CODON))
				usATGPos= i;
		if (usATGPos== 0)
			return false;	// no upstream atg >35 aa found
		
		seq= seq.substring(usATGPos, startPos- frame);
		for (i = 0; (i+3) < seq.length(); i+=3) 
			if (seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[0])||
				seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[1])||
				seq.substring(i, i+3).equalsIgnoreCase(Translation.STOP_CODONS[2]))
					return false;	// inframe stop found, ATG is not internal to other ORF
		return true;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS > 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans
	 * @param maxDistAA
	 * @return
	 */
	public boolean hasUsORF(Translation trans, int minSizeAA) {
		int minSizeNt= minSizeAA* 3;
		int startPos= trans.getTranscript().getExonicPosition(trans.get5PrimeEdge());
		if (startPos< minSizeNt)
			return false;	// cannot have a us ORF of minSize
	
		Translation[] uORFs= trans.getUsORF();
		int cntORF= 0;
		Vector usORFV= new Vector();
		for (int i = 0; uORFs!= null&& i < uORFs.length; i++) {
			if (uORFs[i].getSplicedLength()< minSizeNt)	// bit redundant since already not predicted, but important when predMinSize<> minSize here
				continue;
			++cntORF;
			usORFV.add(uORFs[i]);
		}
		
		return (cntORF> 0);
	}
	
	void checkATG() {
		Translation trans= trpt.getTranslations()[0];
//		int[] ncPos= trans.getNCPosAA();
//		
//		if (ncPos== null|| ncPos.length< 1)
//			return false;
		
		PrintStream p= null;
		try {
			FileOutputStream f= new FileOutputStream("_nonATGseqs.txt", true);
			p= new PrintStream(f);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String seq= trpt.getCDSSequenceNt().toUpperCase();
		if (!seq.startsWith("ATG")) {
			p.println(">"+trans.getTranscript().getTranscriptID()+
					" "+ trans.getTranscript().getChromosome()+
					" "+ trans.getStrand());
			p.println(seq);
		} else {
			if (trans.getStrand()< 0)
				System.out.println(trans.getTranscript().getTranscriptID());
		}
		p.flush();
		p.close();

	}
	
	public boolean isNMD(Translation trln) {
		boolean ejcDS= isTerminatingUpstreamOfEJC(trln, MIN_DIST_NC_EJC);
		boolean atgUS= false;	//hasUsORF(trln, MIN_ORF_LENGTH_AA);
		return (ejcDS|| atgUS);
	}
	public boolean isNMD() {
		Translation trln= trpt.findLongestORF();
		return isNMD(trln);
	}
}
