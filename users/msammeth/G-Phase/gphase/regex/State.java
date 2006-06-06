/*
 * Created on May 5, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.regex;

import gphase.tools.Arrays;

/**
 * 
 * 
 * @author msammeth
 */
public class State {

	static int idCounter= 0;
	
	int id;
	boolean start= false;
	boolean end= false;
	Transition[] transitions;
	
	public State()  {
		this.id= idCounter++;
	}
	
	public void addTransition(int label, boolean newDonor, State toState) {
		Transition t= new Transition(label, newDonor, toState);
		if (transitions== null) {
			transitions= new Transition[] {t};
			return;
		}
		transitions= (Transition[]) Arrays.extendField(transitions, t); 
	}
	
	public String toString() {
		String result= "";
		for (int i = 0; i < transitions.length; i++) 
			result+= transitions[i]+ ",";
		return result;
	}

	public boolean isEnd() {
		return end;
	}
	public void setEnd(boolean end) {
		this.end = end;
	}
	public boolean isStart() {
		return start;
	}
	public void setStart(boolean start) {
		this.start = start;
	}

	public Transition[] getTransitions() {
		return transitions;
	}
	
	public boolean hasHigherTransition(int val) {
		for (int i = 0; i < transitions.length; i++) 
			if (transitions[i].getLabel()> val)
				return true;
		return false;
	}
}
