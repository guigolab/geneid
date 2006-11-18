package gphase;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Transcript;
import gphase.model.Translation;

public class NMDSimulator {
	public static int MIN_DIST_NC_EJC= 50;
	public static int MAX_DIST_NC_EJC= 55;
	public static int MIN_ORF_LENGTH_AA= 35;

	public static void checkNMD() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		g.filterNonCodingTranscripts();
		
		NMDSimulator sim;
		int cntNMD= 0, cntAll= 0;
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();
			for (int j = 0; j < t.length; j++) {
				sim= new NMDSimulator(t[j]);
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

	public static void testNMD() {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);		
		g.filterNonCodingTranscripts();
		
		int cntWrong= 0, cntCorr= 0, cntWrongNMD= 0, cntCorrNMD= 0;
		HashMap wrongHash= new HashMap();
		Gene[] ge= g.getGenes();
		Translation predORF, annORF;
		NMDSimulator sim;
		for (int i = 0; i < ge.length; i++) {
			Transcript[] t= ge[i].getTranscripts();			
			for (int j = 0; j < t.length; j++) {
				annORF= t[j].getTranslations()[0];
				if (annORF.get5PrimeEdge()== t[j].get5PrimeEdge()||
						annORF.get3PrimeEdge()== t[j].get3PrimeEdge()) {
					System.out.println("skipping (open ended HAVANA frame) "+ t[j]);
					continue;
				}
				if (annORF.getSplicedLength()< NMDSimulator.MIN_ORF_LENGTH_AA* 3) {
					System.out.println("skipping (HAVANA ORF too small) "+ t[j]);
					continue;
				}
				predORF= t[j].findLongestORF();
				sim= new NMDSimulator(t[j]);
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
	
	public static void main(String[] args) {
		testNMD();
	}
	
	Transcript trpt= null;
	
	public NMDSimulator(Transcript t) {
		this.trpt= t;
	}
	
	/**
	 * The stop codon must be in the last exon or no further than 50bp`
	 * from the end of the penultimate exon. [HAVANA]
	 * 
	 * @param minDistNts
	 * @return
	 */
	public boolean isTerminatingUpstreamOfEJC(Translation trans, int minDistNt) {
		int critPoint= trpt.getGenomicPosition(trpt.getExonicPosition(trans.get3PrimeEdge())+ minDistNt);
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
		for (int i = 0; uORFs!= null&& i < uORFs.length; i++) {
			if (uORFs[i].getSplicedLength()< minSizeNt)	// bit redundant since already not predicted, but important when predMinSize<> minSize here
				continue;
			++cntORF;
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
		boolean atgUS= false; //hasUsORF(trln, MIN_ORF_LENGTH_AA);
		return (ejcDS|| atgUS);
	}
	public boolean isNMD() {
		Translation trln= trpt.findLongestORF();
		return isNMD(trln);
	}
}
