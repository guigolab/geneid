/*
 * Created on Dec 8, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import repeat.data.Repeat;

/**
 * 
 * 
 * @author micha
 */
public class BAliRepeatAlignment extends BAliManager {

	/**
	 * @param newFileBase
	 */
	public BAliRepeatAlignment(String newFileBase) {
		super(newFileBase);
		// TODO Auto-generated constructor stub
	}
	
	public void setRepeats(Repeat[] newRepeats) {
		this.repeats= newRepeats;
	}
	
	public void setSeqNames(String[] newSeqNames) {
		this.seqNames= newSeqNames;
	}
}
