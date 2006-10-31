/*
 * Created on May 5, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.regex;

import gphase.tools.Arrays;

import java.util.Comparator;
import java.util.Stack;
import java.util.Vector;

/**
 * 
 * 
 * @author msammeth
 */
public class Automaton {
	
	State initialState;
	State endState;
	State currentState;
	Stack backtrackStack;
	
	public Automaton() {
		initialState= new State();
		currentState= initialState;
		endState= initialState;
	}
	
	public void sortTransitions() {
		sortTransitions(initialState, new Transition.TransitionComparator(), -1);
	}
	
	private void sortTransitions(State s, Comparator c, int val) {
		
		Transition[] trans= s.getTransitions();
		if (trans== null)
			return;
		
		java.util.Arrays.sort(trans, c);
		for (int i = 0; i < trans.length; i++) 
			if (trans[i].getLabel()> val)
				sortTransitions(trans[i].getToState(), c, trans[i].getLabel());
			
	}
	
	/**
	 * 
	 * @deprecated ambigous, no longer in use
	 */
	public void move(int val) {
		Transition[] trans= currentState.getTransitions();
		for (int i = 0; i < trans.length; i++) {
			if (trans[i].getLabel()== val) {
				currentState= trans[i].getToState();
				return;
			}
		}		
	}
	
	public void move(Transition t) {
		Transition[] trans= currentState.getTransitions();
		for (int i = 0; i < trans.length; i++) {
			if (trans[i]== t) {
				currentState= trans[i].getToState();
				
				Transition[] trans2= currentState.getTransitions();	// skip single epsilon edges
				if (trans2!= null&& trans2.length== 1&& trans2[0].getLabel()== -1)
					currentState= trans2[0].toState;
				
				return;
			}
		}		
		System.err.println("Transition "+t+" not found!");
	}
	
	
	public Transition[] getChunk(int val) {
		
		Vector result= new Vector();
		Transition[] trans= currentState.getTransitions();
		for (int i = 0; i < trans.length; i++) 
			if (trans[i].getLabel()== val)
				result.add(trans[i]);
		
		return (Transition[]) Arrays.toField(result);
	}
	
	public Transition createTransition(int label, boolean donor) {
		State newState= new State();
		Transition t= endState.addTransition(label, donor, newState);
		endState= newState;
		return t;
	}
	public State getInitialState() {
		return initialState;
	}
	public State getEndState() {
		return endState;
	}
	
	public void step(Transition t) {
		Transition[] trans= getCurrentState().getTransitions();
		for (int i = 0; i < trans.length; i++) 
			if (trans[i].equals(t)) {
				currentState= trans[i].toState;
				return;
			}
		System.err.println("Transition "+ t+ " not found!");
	}
	
	public void init() {
		currentState= initialState;
	}
	
	public String toString() {
		return initialState.toString();
	}
	public State getCurrentState() {
		return currentState;
	}
	
	public Transition getCurrentMinTransition() {
		Transition[] trans= currentState.getTransitions();
		if (trans== null)
			return null;
		int min= Integer.MAX_VALUE;
		Transition result= null;
		for (int i = 0; i < trans.length; i++) {
			if (trans[i].getLabel()< min) {
				min= trans[i].getLabel();
				result= trans[i];
			} else if (trans[i].getLabel()== min)
				System.err.println("Two outgoing edges with same label "+ trans[i].getLabel());
		}
		
		return result;
	}

	public Transition getCurrentMaxTransition() {
		Transition[] trans= currentState.getTransitions();
		if (trans== null)
			return null;
		int max= Integer.MIN_VALUE;
		Transition result= null;
		for (int i = 0; i < trans.length; i++) {
			if (trans[i].getLabel()> max) {
				max= trans[i].getLabel();
				result= trans[i];
			} else if (trans[i].getLabel()== max)
				System.err.println("Two outgoing edges with same label "+ trans[i].getLabel());
		}
		
		return result;
	}
	
	public boolean isEndState() {
		return currentState== endState;
	}

	public void setCurrentState(State currentState) {
		this.currentState = currentState;
	}
}
