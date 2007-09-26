/*
 * Created on Apr 27, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.util.Vector;

import qalign.algo.dialign.MultiFrag;
import qalign.algo.dialign.MultiFragExt;
import qalign.tools.MSFWrapper;

/**
 * 
 * 
 * @author micha
 */
public class RepeatMSFWrapper extends MSFWrapper {

	protected MultiFragExt[] frags= null;
	boolean highlight= true;
	
	public static void main(String[] args) {
	}

	public RepeatMSFWrapper(String absFName) {
		super(absFName);
		highlight= true;
	}
	
	public RepeatMSFWrapper(String newFName, String newFPath) {
		super(newFName, newFPath);
		highlight= true;
	}
	
	public void setRepeatFragments(MultiFrag[] newFrags) {
		
			// catch
		if (newFrags== null)
			return;
			
			// collect in vector
		Vector fragVec= new Vector(newFrags.length);
		for (int i= 0; i < newFrags.length; i++) 
			try {
				fragVec.add((MultiFragExt) newFrags[i]);
			} catch (ClassCastException e) {
				fragVec.add(new MultiFragExt(newFrags[i])); // continue
			}
		
			// convert to array
		MultiFragExt[] castFrags= new MultiFragExt[fragVec.size()];
		for (int i= 0; i < castFrags.length; i++) 
			castFrags[i]= (MultiFragExt) fragVec.elementAt(i);
			
			// call setter
		setRepeatFragments(castFrags);
	}
	
	public void setRepeatFragments(MultiFragExt[] newFrags) {
		
		this.frags= newFrags;
		if (frags!= null&& sequences!= null)
			formatSequences();
	}
	
	public void setSequences(String[] newSeqs) {
		
		super.setSequences(newSeqs);
		if (frags!= null&& sequences!= null)
			formatSequences();
	}
	
	public void formatSequences() {
		
		if (!highlight)
			return;
			
			// init w lowercase
		for (int i= 0; i < sequences.length; i++) 
			sequences[i]= sequences[i].toLowerCase();
		
			// rise all letters contained in fragments to uppercase
		for (int i= 0; i < frags.length; i++) {
			
			MultiFragExt frag= frags[i];
			
			if (!frag.isRepeat())
				continue;
			
				// first seq of fragment
			int start= frag.getSequenceStart(true);
			int pstart;
			for (pstart= 0; start>= 0; pstart++) 
				if (Character.isLetter(sequences[frag.getSequenceNo(true)].charAt(pstart)))
					--start;
			--pstart;
			sequences[frag.getSequenceNo(true)]= 
				sequences[frag.getSequenceNo(true)].substring(0, pstart)+
				sequences[frag.getSequenceNo(true)].substring(pstart, pstart+frag.getLength()).toUpperCase()+
				sequences[frag.getSequenceNo(true)].substring(pstart+ frag.getLength());
				
				// second seq of fragment
			start= frag.getSequenceStart(false);
			for (pstart= 0; start>= 0; pstart++)
				if (Character.isLetter(sequences[frag.getSequenceNo(false)].charAt(pstart)))
					--start;
			--pstart;
			sequences[frag.getSequenceNo(false)]= 
				sequences[frag.getSequenceNo(false)].substring(0, pstart)+
				sequences[frag.getSequenceNo(false)].substring(pstart, pstart+frag.getLength()).toUpperCase()+
				sequences[frag.getSequenceNo(false)].substring(pstart+ frag.getLength());
		}
	}
	/**
	 * @return
	 */
	public boolean isHighlight() {
		return highlight;
	}

	/**
	 * @param b
	 */
	public void setHighlight(boolean b) {
		highlight = b;
	}

}
