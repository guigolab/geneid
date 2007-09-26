/*
 * Created on Mar 3, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model_heavy;

import java.io.Serializable;

/**
 * 
 * 
 * @author msammeth
 */
public class Phase implements Serializable {

	int startPhase= -1;
	int endPhase= -1;

	public Phase(int newStartPhase, int newEndPhase) {
		this.startPhase= newStartPhase;
		this.endPhase= newEndPhase;
	}

}
