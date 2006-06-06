/*
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.Serializable;
import java.util.HashMap;

/**
 * 
 * 
 * @author micha
 */
public class Translation extends DirectedRegion {

	static final long serialVersionUID= 8996021902187779155L;
	int translationID= -1;
	String translationID_triv= null;
	Transcript transcript= null;
	public Translation(Transcript newTranscript) {
		this.transcript= newTranscript;
		this.strand= getTranscript().getStrand();
	}
		
	public Translation(Transcript newTranscript, String stableTranslationID) {

		this(newTranscript);
			
		if (stableTranslationID.startsWith(getSpecies().getEnsemblPrefix())) 
			translationID= Integer.parseInt(Graph.decodeStableID(stableTranslationID, true));
		else	// HOX-genes in tetraodon: eg HOXA1, HOXAb10, HOXB-EVX, HOXB1, ...
			translationID_triv= stableTranslationID;
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
	public int getTranslationID() {
		return translationID;
	}
	/**
	 * @param translationID The translationID to set.
	 */
	public void setTranslationID(int translationID) {
		this.translationID = translationID;
	}

}
