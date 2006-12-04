/*
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.NMDSimulator;
import gphase.tools.Arrays;

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
		for (int i = 0; i < trns.length; i++) {
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
