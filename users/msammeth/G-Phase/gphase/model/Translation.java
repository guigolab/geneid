/*
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.NMDSimulator;
import gphase.StopCodons;
import gphase.tools.Arrays;
import gphase.tools.IntVector;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Vector;

/**
 * 
 * 
 * @author micha
 */
public class Translation extends DirectedRegion {

	static final long serialVersionUID= 8996021902187779155L;
	
	public static final String START_CODON= "ATG";
	public static final String[] STOP_CODONS= new String[] {"TAA", "TGA", "TAG"};
	
	String translationID= null;
	Transcript transcript= null;
	int splicedLength= -1;
	
	public static int[] getCodonCount(String[] codons, String seq) {
		int[] res= new int[3];
		for (int i = 0; i < res.length; i++) 
			res[i]= getCodonCount(codons, seq, i);
		return res;
	}
	
	public static int getCodonCount(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		int cnt= 0;
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				++cnt;
		}
		return cnt;
	}

	/**
	 * condenses all 3 frames to one position array
	 * @param codons
	 * @param seq
	 * @return
	 */
	public static int[] getCodonPositions(String[] codons, String seq) {
		IntVector res= new IntVector();
		for (int i = 0; i < 3; i++) { 
			int[] pos= getCodonPositions(codons, seq, i);
			for (int j = 0; j < pos.length; j++) {
				int ins= java.util.Arrays.binarySearch(res.toIntArray(), pos[j]);
				if (ins< 0) 
					res.insert(pos[j], ins);
			}
		}
		return res.toIntArray();
	}
	
	public static int[] getCodonPositions(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		IntVector pos= new IntVector();
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				pos.add(i);
		}
		return pos.toIntArray();
	}
	
	public static int[] getStartCount(String seq) {
		return getCodonCount(new String[] {START_CODON}, seq);
	}
	
	public static int[] getStopCount(String seq) {
		return getCodonCount(STOP_CODONS, seq);
	}
	
	public Translation(Transcript newTranscript, int newStart, int newEnd, int newStrand) {
		super(newStart, newEnd, newStrand);
		this.transcript= newTranscript;
	}
	
	public Translation(Transcript newTranscript) {
		this.transcript= newTranscript;
		this.strand= getTranscript().getStrand();
		setID("translation");
		setStrand(getTranscript().getGene().getStrand());
	}
	
	public int getSplicedLength() {
		if (splicedLength< 0) {
			DirectedRegion[] regs= transcript.getCDSRegions();
			splicedLength= 0;
			for (int i = 0; i < regs.length; i++) 
				splicedLength+= regs[i].getLength();
		}
		return splicedLength;
	}
	
	public Translation(Transcript newTranscript, String stableTranslationID) {

		this(newTranscript);
		this.translationID= stableTranslationID;
	}

	public String getChromosome() {
		return getTranscript().getChromosome();
	}
	/**
	 * @return Returns the transcript.
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	
	public Species getSpecies() {
		return getTranscript().getSpecies();
	}
	/**
	 * @return Returns the translationID.
	 */
	public String getTranslationID() {
		return translationID;
	}
	/**
	 * @param translationID The translationID to set.
	 */
	public void setTranslationID(String newTranslationID) {
		this.translationID = newTranslationID;
	}

	public void setSplicedLength(int splicedLength) {
		this.splicedLength = splicedLength;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS > 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans
	 * @param maxDistAA
	 * @return
	 */
	public Translation[] getUsORF() {
		
		Translation[] trns= getTranscript().getAllORFs();
		Vector uOrfV= new Vector();
		for (int i = 0; trns!= null&& i < trns.length; i++) {
			if (trns[i].get3PrimeEdge()< get5PrimeEdge())
				uOrfV.add(trns[i]);
		}
		
		Translation[] uOrf= (Translation[]) Arrays.toField(uOrfV);
		return uOrf;
	}
	
	public boolean isOpenEnded5() {
		return (get5PrimeEdge()== transcript.get5PrimeEdge());
	}

	public boolean isOpenEnded3() {
		return (get3PrimeEdge()== transcript.get3PrimeEdge());
	}
	
	public boolean isOpenEnded() {
		// TODO: check seq for ATG/stop
		if (isOpenEnded5()|| isOpenEnded3())
			return true;
		return false;
	}
	
	

	
}
