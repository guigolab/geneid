/*
 * Created on May 5, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.regex;

import java.util.Arrays;
import java.util.Comparator;

/**
 * 
 * 
 * @author msammeth
 */
public class Transition {

	public static class TransitionComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			Transition t1= (Transition) arg0;
			Transition t2= (Transition) arg1;
			if (t1.getLabel()< t2.getLabel())
				return -1;
			if (t1.getLabel()> t2.getLabel())
				return 1;
			return 0;
		}
	}
	int label;
	State toState;
	boolean donor= false;
	public boolean isDonor() {
		return donor;
	}

	public Transition(int newLabel, boolean newDonor, State target) {
		this.label= newLabel;
		this.donor= newDonor;
		this.toState= target;
	}
	
	public String toString() {
		String result= Integer.toString(label);
		if (donor)
			result+= "^";
		else
			result+= "=";
		return result;
	}

	public int getLabel() {
		return label;
	}

	public State getToState() {
		return toState;
	}
}
