/*
 * Created on Mar 12, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model_heavy;

import gphase.tools.Arrays;

/**
 * 
 * 
 * @author msammeth
 */
public class TSSite extends AbstractSite {
	boolean tsStart= false;
	Transcript[] transcripts= null;
	
	public TSSite(Transcript newTrans, int newPos, boolean isTSS) {
		super(newPos);
		this.tsStart= isTSS;
		this.pos= newPos;
		addTranscript(newTrans);
	}
	
	public Transcript[] getTranscripts() {
		return transcripts;
	}
	public boolean addTranscript(Transcript transcript) {
		if (transcripts== null) { 
			transcripts= new Transcript[]{transcript};
			return true;
		}
		
		for (int i = 0; i < transcripts.length; i++) 
			if (transcripts[i].equals(transcript))
				return false;
		
		transcripts= (Transcript[]) Arrays.add(transcripts, transcript);
		return true;
	}
	public boolean isTsStart() {
		return tsStart;
	}
	
}
