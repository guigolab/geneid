/*
 * dk.brics.automaton
 * Copyright (C) 2001-2004 Anders Moeller
 * All rights reserved.
 */

package gphase.regregex;

import java.util.*;
import java.io.*;

/**
 * <tt>Automaton</tt> state.
 * 
 * @author Anders M&oslash;ller &lt; <a
 *         href="mailto:amoeller@brics.dk">amoeller@brics.dk </a>&gt;
 */
public class State implements Serializable, Comparable {
	static final long serialVersionUID = 30001;

	boolean accept;

	HashSet transitions;

	int number;

	int id;

	static int next_id;

	/** Constructs new state. Initially, the new state is a reject state. */
	public State() {
		resetTransitions();
		id = next_id++;
	}

	/** Resets transition set. */
	void resetTransitions() {
		transitions = new HashSet();
	}

	/**
	 * Returns set of outgoing transitions. Subsequent changes are reflected in
	 * the automaton.
	 * 
	 * @return transition set
	 */
	public Set getTransitions() {
		return transitions;
	}

	/**
	 * Adds outgoing transition.
	 * 
	 * @param t
	 *            transition
	 */
	public void addTransition(Transition t) {
		transitions.add(t);
	}

	/**
	 * Sets acceptance for this state.
	 * 
	 * @param accept
	 *            if true, this state is an accept state
	 */
	public void setAccept(boolean accept) {
		this.accept = accept;
	}

	/**
	 * Returns acceptance status.
	 * 
	 * @return true is this is an accept state
	 */
	public boolean isAccept() {
		return accept;
	}

	/**
	 * Performs lookup in transitions.
	 * 
	 * @param c
	 *            character to look up
	 * @return destination state, null if no matching outgoing transition
	 */
	public State step(char c) {
		Iterator i = transitions.iterator();
		while (i.hasNext()) {
			Transition t = (Transition) i.next();
			if (t.min <= c && c <= t.max)
				return t.to;
		}
		return null;
	}

	void addEpsilon(State to) {
		if (to.accept)
			accept = true;
		Iterator i = to.transitions.iterator();
		while (i.hasNext()) {
			Transition t = (Transition) i.next();
			transitions.add(t);
		}
	}

	/**
	 * Returns transitions sorted by (min, reverse max, to) or (to, min, reverse
	 * max)
	 */
	Transition[] getSortedTransitionArray(boolean to_first) {
		Transition[] e = (Transition[]) transitions.toArray(new Transition[0]);
		TransitionComparator c = new TransitionComparator();
		c.to_first = to_first;
		Arrays.sort(e, c);
		return e;
	}

	List getSortedTransitions(boolean to_first) {
		return Arrays.asList(getSortedTransitionArray(to_first));
	}

	/**
	 * Returns string describing this state. Normally invoked via
	 * {@link Automaton#toString()}.
	 */
	public String toString() {
		StringBuffer b = new StringBuffer();
		b.append("state ").append(number);
		if (accept)
			b.append(" [accept]");
		else
			b.append(" [reject]");
		b.append(":\n");
		Iterator i = transitions.iterator();
		while (i.hasNext()) {
			Transition t = (Transition) i.next();
			b.append("  ").append(t.toString()).append("\n");
		}
		return b.toString();
	}

	/**
	 * Compares this object with the specified object for order. States are
	 * ordered by the time of construction.
	 */
	public int compareTo(Object o) {
		return ((State) o).id - id;
	}
}