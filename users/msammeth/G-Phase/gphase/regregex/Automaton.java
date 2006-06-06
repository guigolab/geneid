package gphase.regregex;

/*
 * dk.brics.automaton
 * Copyright (C) 2001-2004 Anders Moeller
 * All rights reserved.
 */

import java.util.*;
import java.io.*;
import java.net.*;

/* Class invariants:
 *
 * - An automaton is either represented explicitly (with State and Transition objects)
 *   or with a singleton string in case the automaton accepts exactly one string.
 * - Automata are always reduced (see reduce()) 
 *   and have no transitions to dead states (see removeDeadTransitions()).
 * - If an automaton is nondeterministic, then isDeterministic() returns false (but
 *   the converse is not required).
 */

/**
 * Finite-state automaton with regular expression operations.
 * <p>
 * Automata are represented using {@link State} and {@link Transition} objects.
 * Implicitly, all states and transitions of an automaton are reachable from its initial state.
 * If the states or transitions are manipulated manually, the {@link #restoreInvariant()}
 * and {@link #setDeterministic(boolean)} methods should be used afterwards to restore 
 * certain representation invariants that are assumed by the built-in automata operations.
 * @author Anders M&oslash;ller &lt;<a href="mailto:amoeller@brics.dk">amoeller@brics.dk</a>&gt;
 */
public class Automaton implements Serializable
{
    static final long serialVersionUID = 10001;

    /**
     * Minimize using Huffman's O(n<sup>2</sup>) algorithm. 
     * This is the standard text-book algorithm.
     * @see #setMinimization(int)
     */
    public static final int MINIMIZE_HUFFMAN = 0;

    /**
     * Minimize using Brzozowski's O(2<sup>n</sup>) algorithm. 
     * This algorithm uses the reverse-determinize-reverse-determinize trick, which has a bad
     * worst-case behavior but often works very well in practice 
     * (even better than Hopcroft's!).
     * @see #setMinimization(int)
     */
    public static final int MINIMIZE_BRZOZOWSKI = 1;

    /**
     * Minimize using Hopcroft's O(n log n) algorithm.
     * This is regarded as one of the most generally efficient algorithms that exist.
     * @see #setMinimization(int)
     */
    public static final int MINIMIZE_HOPCROFT = 2;

    /** Selects minimization algorithm (default: <code>MINIMIZE_HOPCROFT</code>). */
    static int minimization = MINIMIZE_HOPCROFT;

    /** Initial state of this automaton. */
    State initial;

    /** If true, then this automaton is definitely deterministic 
	(i.e., there are no choices for any run, but a run may crash). */
    boolean deterministic;

    /** Extra data associated with this automaton. */
    Object info;

    /** Hash code. Recomputed by {@link #minimize()}. */
    int hash_code;
    
    /** Singleton string. Null if not applicable. */
    String singleton;

    /** Minimize always flag. */
    static boolean minimize_always;

    /** 
     * Constructs new automaton that accepts the empty language.
     * Using this constructor, automata can be constructed manually from
     * {@link State} and {@link Transition} objects.
     * @see #setInitialState(State)
     * @see State
     * @see Transition
     */
    public Automaton() 
    {
	initial = new State();
	deterministic = true;
    }

    /** 
     * Selects minimization algorithm (default: <code>MINIMIZE_HOPCROFT</code>). 
     * @param algorithm minimization algorithm
     */
    static public void setMinimization(int algorithm)
    {
	minimization = algorithm;
    }

    /**
     * Sets or resets minimize always flag.
     * If this flag is set, then {@link #minimize()} will automatically
     * be invoked after all operations that otherwise may produce non-minimal automata.
     * By default, the flag is not set.
     * @param flag if true, the flag is set
     */
    static public void setMinimizeAlways(boolean flag)
    {
	minimize_always = flag;
    }

    void checkMinimizeAlways()
    {
	if (minimize_always)
	    minimize();
    }

    boolean isSingleton()
    {
	return singleton!=null;
    }

    /** 
     * Sets initial state. 
     * @param s state
     */
    public void setInitialState(State s)
    {
	initial = s;
	singleton = null;
    }

    /** 
     * Gets initial state. 
     * @return state
     */
    public State getInitialState()
    {
	expandSingleton();
	return initial;
    }

    /**
     * Returns deterministic flag for this automaton.
     * @return true if the automaton is definitely deterministic, false if the automaton
     *         may be nondeterministic
     */
    public boolean isDeterministic()
    {
	return deterministic;
    }

    /**
     * Sets deterministic flag for this automaton.
     * This method should (only) be used if automata are constructed manually.
     * @param deterministic true if the automaton is definitely deterministic, false if the automaton
     *                      may be nondeterministic
     */
    public void setDeterministic(boolean deterministic)
    {
	this.deterministic = deterministic;
    }

    /**
     * Associates extra information with this automaton. 
     * @param info extra information
     */
    public void setInfo(Object info)
    {
	this.info = info;
    }

    /**
     * Returns extra information associated with this automaton. 
     * @return extra information
     * @see #setInfo(Object)
     */
    public Object getInfo()
    {
	return info;
    }

    /** 
     * Returns the set states that are reachable from the initial state.
     * @return set of {@link State} objects
     */
    public Set getStates()
    {
	expandSingleton();
	HashSet visited = new HashSet();
	LinkedList worklist = new LinkedList();
	worklist.add(initial);
	visited.add(initial);
	while (worklist.size()>0) {
	    State s = (State) worklist.removeFirst();
	    Iterator i = s.transitions.iterator();
	    while (i.hasNext()) {
		Transition t = (Transition) i.next();
		if (!visited.contains(t.to)) {
		    visited.add(t.to);
		    worklist.add(t.to);
		}
	    }
	}
	return visited;
    }

    /** 
     * Returns the set of reachable accept states. 
     * @return set of {@link State} objects
     */
    public Set getAcceptStates()
    {
	expandSingleton();
	HashSet accepts = new HashSet();
	HashSet visited = new HashSet();
	LinkedList worklist = new LinkedList();
	worklist.add(initial);
	visited.add(initial);
	while (worklist.size()>0) {
	    State s = (State) worklist.removeFirst();
	    if (s.accept)
		accepts.add(s);
	    Iterator i = s.transitions.iterator();
	    while (i.hasNext()) {
		Transition t = (Transition) i.next();
		if (!visited.contains(t.to)) {
		    visited.add(t.to);
		    worklist.add(t.to);
		}
	    }
	}
	return accepts;
    }

    /** Assigns consecutive numbers to the given states. */
    void setStateNumbers(Set states)
    {
	Iterator i = states.iterator();
	int number = 0;
	while (i.hasNext()) {
	    State s = (State) i.next();
	    s.number = number++;
	}
    }

    /** 
     * Checks whether there is a loop containing s. (This is sufficient since 
     * there are never transitions to dead states.) 
     */
    boolean isFinite(State s, HashSet path)
    {
	path.add(s);
	Iterator i = s.transitions.iterator();
	while (i.hasNext()) {
	    Transition t = (Transition) i.next();
	    if (path.contains(t.to) || !isFinite(t.to, path))
		return false;
	}
	path.remove(s);
	return true;
    }

    /** Returns the strings that can be produced from s, returns false if more than 
	<code>limit</code> strings are found. <code>limit</code>&lt;0 means "infinite". */
    boolean getFiniteStrings(State s, HashSet pathstates, HashSet strings, StringBuffer path, int limit)
    {
	pathstates.add(s);
	Iterator i = s.transitions.iterator();
	while (i.hasNext()) {
	    Transition t = (Transition) i.next();
	    if (pathstates.contains(t.to))
		return false;
	    for (int n = t.min; n<=t.max; n++) {
		path.append((char) n);
		if (t.to.accept) {
		    strings.add(path.toString());
		    if (limit>=0 && strings.size()>limit)
			return false;
		}
		if (!getFiniteStrings(t.to, pathstates, strings, path, limit))
		    return false;
		path.deleteCharAt(path.length()-1);
	    }
	}
	pathstates.remove(s);
	return true;
    }

    /** Adds transitions to explicit crash state to ensure that
	transition function is total. */
    void totalize()
    {
	State s = new State();
	s.transitions.add(new Transition(Character.MIN_VALUE, Character.MAX_VALUE, s));
	Iterator i = getStates().iterator();
	while (i.hasNext()) {
	    State p = (State) i.next();
	    Iterator j = p.getSortedTransitions(false).iterator();
	    int maxi = Character.MIN_VALUE;
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		if (t.min>maxi)
		    p.transitions.add(new Transition((char) maxi, (char) (t.min-1), s));
		if (t.max+1>maxi)
		    maxi = t.max+1;
	    }
	    if (maxi<=Character.MAX_VALUE)
		p.transitions.add(new Transition((char) maxi, Character.MAX_VALUE, s));
	}
    }

    /**
     * Restores representation invariant.
     * This method must be invoked before any built-in automata operation is performed 
     * if automaton states or transitions are manipulated manually.
     * @see #setDeterministic(boolean)
     */
    public void restoreInvariant()
    {
	removeDeadTransitions();
	hash_code = 0;
    }
    
    /** 
     * Reduces this automaton.
     * An automaton is "reduced" by combining overlapping and adjacent edge intervals with same destination. 
     */
    public void reduce()
    {
	if (isSingleton())
	    return;
	Set states = getStates();
	setStateNumbers(states);
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    Iterator j = s.getSortedTransitions(true).iterator();
	    s.resetTransitions();
	    State p = null;
	    int min = -1, max = -1;
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		if (p==t.to) {
		    if (t.min<=max+1) {
			if (t.max>max)
			    max = t.max;
		    } else {
			if (p!=null)
			    s.transitions.add(new Transition((char) min, (char) max, p));
			min = t.min;
			max = t.max;
		    }
		} else {
		    if (p!=null)
			s.transitions.add(new Transition((char) min, (char) max, p));
		    p = t.to;
		    min = t.min;
		    max = t.max;
		}
	    }
	    if (p!=null)
		s.transitions.add(new Transition((char) min, (char) max, p));
	}
    }

    /** Gets sorted array of all interval start points. */
    char[] getStartPoints()
    {
	Set pointset = new HashSet();
	Iterator i = getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    pointset.add(new Character(Character.MIN_VALUE));
	    Iterator j = s.transitions.iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		pointset.add(new Character(t.min));
		if (t.max<Character.MAX_VALUE)
		    pointset.add(new Character((char) (t.max+1)));
	    }
	}
	char[] points = new char[pointset.size()];
	Iterator k = pointset.iterator();
	int n = 0;
	while (k.hasNext()) {
	    Character m = (Character) k.next();
	    points[n++] = m.charValue();
	}
	Arrays.sort(points);
	return points;
    }

    /** 
     * Returns set of live states. A state is "live" if an accept state is reachable from it. 
     * @return set of {@link State} objects
     */
    public Set getLiveStates()
    {
	expandSingleton();
	return getLiveStates(getStates());
    }

    Set getLiveStates(Set states)
    {
	HashMap map = new HashMap();
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    map.put(s, new HashSet());
	}
	i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    Iterator j = s.transitions.iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		((HashSet) map.get(t.to)).add(s);
	    }
	}
	Set live = new HashSet(getAcceptStates());
	LinkedList worklist = new LinkedList(live);
	while (worklist.size()>0) {
	    State s = (State) worklist.removeFirst();
	    Iterator j = ((HashSet) map.get(s)).iterator();
	    while (j.hasNext()) {
		State p = (State) j.next();
		if (!live.contains(p)) {
		    live.add(p);
		    worklist.add(p);
		}
	    }
	}	
	return live;
    }

    /** 
     * Removes transitions to dead states and calls {@link #reduce()}
     * (a state is "dead" if no accept state is reachable from it). 
     */
    public void removeDeadTransitions()
    {
	if (isSingleton())
	    return;
	Set states = getStates();
	Set live = getLiveStates(states);
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    Iterator j = s.transitions.iterator();
	    s.resetTransitions();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		if (live.contains(t.to))
		    s.transitions.add(t);
	    }
	}
	reduce();
    }

    /** Returns sorted array of transitions for each state (and sets state numbers). */
    Transition[][] getSortedTransitions(Set states)
    {
	setStateNumbers(states);
	Transition[][] transitions = new Transition[states.size()][];
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    transitions[s.number] = s.getSortedTransitionArray(false);
	}
	return transitions;
    }
    
    /** Returns new (deterministic) automaton with the empty language. */
    public static Automaton makeEmpty()
    {
	Automaton a = new Automaton();
	State s = new State();
	a.initial = s;
	a.deterministic = true;
	return a;
    }

    /** Returns new (deterministic) automaton that accepts only the empty string. */
    public static Automaton makeEmptyString()
    {
	Automaton a = new Automaton();
	a.singleton = "";
	a.deterministic = true;
	return a;
    }

    /** Returns new (deterministic) automaton that accepts all strings. */
    public static Automaton makeAnyString()
    {
	Automaton a = new Automaton();
	State s = new State();
	a.initial = s;
	s.accept = true;
	s.transitions.add(new Transition(Character.MIN_VALUE, Character.MAX_VALUE, s));
	a.deterministic = true;
	return a;
    }

    /** Returns new (deterministic) automaton that accepts any single character. */
    public static Automaton makeAnyChar()
    {
	return makeCharRange(Character.MIN_VALUE, Character.MAX_VALUE);
    }

    /** Returns new (deterministic) automaton that accepts a single 
	character of the given value. */
    public static Automaton makeChar(char c)
    {
	return makeCharRange(c, c);
    }

    /** Returns new (deterministic) automaton that accepts a single 
	char whose value is in the given interval (including both end points). */
    public static Automaton makeCharRange(char min, char max)
    {
	Automaton a = new Automaton();
	State s1 = new State();
	State s2 = new State();
	a.initial = s1;
	s2.accept = true;
	if (min<=max)
	    s1.transitions.add(new Transition(min, max, s2));
	a.deterministic = true;
	return a;
    }

    /** Returns new (deterministic) automaton that accepts 
	a single character in the given set. */
    public static Automaton makeCharSet(String set)
    {
	Automaton a = new Automaton();
	State s1 = new State();
	State s2 = new State();
	a.initial = s1;
	s2.accept = true;
	for (int i = 0; i < set.length(); i++)
	    s1.transitions.add(new Transition(set.charAt(i), s2));
	a.deterministic = true;
	a.reduce();
	return a;
    }

    /**
     * Constructs sub-automaton corresponding to decimal numbers of 
     * length x.substring(n).length().
     */
    private static State anyOfRightLength(String x, int n)
    {
	State s = new State();
	if (x.length()==n)
	    s.setAccept(true);
	else
	    s.addTransition(new Transition('0', '9', anyOfRightLength(x, n+1)));
	return s;
    }

    /**
     * Constructs sub-automaton corresponding to decimal numbers of value 
     * at least x.substring(n) and length x.substring(n).length().
     */
    private static State atLeast(String x, int n, Collection initials, boolean zeros)
    {
	State s = new State();
	if (x.length()==n)
	    s.setAccept(true);
	else {
	    if (zeros)
		initials.add(s);
	    char c = x.charAt(n);
	    s.addTransition(new Transition(c, atLeast(x, n+1, initials, zeros && c=='0')));
	    if (c<'9')
		s.addTransition(new Transition((char) (c+1), '9', anyOfRightLength(x, n+1)));
	}
	return s;
    }

    /**
     * Constructs sub-automaton corresponding to decimal numbers of value 
     * at most x.substring(n) and length x.substring(n).length().
     */
    private static State atMost(String x, int n)
    {
	State s = new State();
	if (x.length()==n)
	    s.setAccept(true);
	else {
	    char c = x.charAt(n);
	    s.addTransition(new Transition(c, atMost(x, (char) n+1)));
	    if (c>'0')
		s.addTransition(new Transition('0', (char) (c-1), anyOfRightLength(x, n+1)));
	}
	return s;
    }

    /**
     * Constructs sub-automaton corresponding to decimal numbers of value 
     * between x.substring(n) and y.substring(n) and of
     * length x.substring(n).length() (which must be equal to y.substring(n).length()).
     */
    private static State between(String x, String y, int n, Collection initials, boolean zeros)
    {
	State s = new State();
	if (x.length()==n)
	    s.setAccept(true);
	else {
	    if (zeros)
		initials.add(s);
	    char cx = x.charAt(n);
	    char cy = y.charAt(n);
	    if (cx==cy) 
		s.addTransition(new Transition(cx, between(x, y, n+1, initials, zeros && cx=='0')));
	    else { // cx<cy
		s.addTransition(new Transition(cx, atLeast(x, n+1, initials, zeros && cx=='0')));
		s.addTransition(new Transition(cy, atMost(y, n+1)));
		if (cx+1<cy)
		    s.addTransition(new Transition((char) (cx+1), (char) (cy-1), anyOfRightLength(x, n+1)));
	    }
	}
	return s;
    }

    /** 
     * Returns new automaton that accepts strings representing 
     * decimal non-negative integers in the given interval.
     * @param min minimal value of interval
     * @param max maximal value of inverval (both end points are included in the interval)
     * @param digits if >0, use fixed number of digits (strings must be prefixed 
     *               by 0's to obtain the right length) -
     *               otherwise, the number of digits is not fixed
     * @exception IllegalArgumentException if min>max or if numbers in the interval cannot be expressed
     *                                     with the given fixed number of digits
     */
    public static Automaton makeInterval(int min, int max, int digits)
	throws IllegalArgumentException
    {
	Automaton a = new Automaton();
	String x = (new Integer(min)).toString();
	String y = (new Integer(max)).toString();
	if (min>max || (digits>0 && y.length()>digits))
	    throw new IllegalArgumentException();
	int d;
	if (digits>0)
	    d = digits;
	else
	    d = y.length();
	StringBuffer bx = new StringBuffer();
	for (int i = x.length(); i<d; i++)
	    bx.append('0');
	bx.append(x);
	x = bx.toString();
	StringBuffer by = new StringBuffer();
	for (int i = y.length(); i<d; i++)
	    by.append('0');
	by.append(y);
	y = by.toString();
	Collection initials = new ArrayList();
	a.initial = between(x, y, 0, initials, digits<=0);
	if (digits<=0) {
	    ArrayList pairs = new ArrayList();
	    Iterator j = initials.iterator();
	    while (j.hasNext()) {
		State p = (State) j.next();
		if (a.initial!=p)
		    pairs.add(new StatePair(a.initial, p));
	    }
	    a.addEpsilons(pairs);
	    a.initial.addTransition(new Transition('0', a.initial));
	} else
	    a.deterministic = true;
	a.checkMinimizeAlways();
	return a;
    }

    /** Expands singleton representation to normal representation. */
    void expandSingleton()
    {
	if (isSingleton()) {
	    State p = new State();
	    initial = p;
	    for (int i = 0; i < singleton.length(); i++) {
		State q = new State();
		p.transitions.add(new Transition(singleton.charAt(i), q));
		p = q;
	    }
	    p.accept = true;
	    deterministic = true;
	    singleton = null;
	}
    }
    
    /** Returns new (deterministic) automaton that accepts
	the single given string. 
	<p>
	Complexity: constant. */
    public static Automaton makeString(String s)
    {
	Automaton a = new Automaton();
	a.singleton = s;
	a.deterministic = true;
	return a;
    }

    /** Returns new automaton that accepts the concatenation of the languages
	of this	and the given automaton. 
	<p>
	Complexity: linear in number of states. */
    public Automaton concatenate(Automaton a)
    {
	if (isSingleton() && a.isSingleton()) 
	    return makeString(singleton + a.singleton);
	a = a.cloneExpanded();
	Automaton b = cloneExpanded();
	Iterator i = b.getAcceptStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    s.accept = false;
	    s.addEpsilon(a.initial);
	}
	b.deterministic = false;
	b.checkMinimizeAlways();
	return b;
    }

    /** Returns new automaton that accepts the concatenation of the languages
	of the given automata.
	<p>
	Complexity: linear in total number of states. */
    static public Automaton concatenate(List l)
    {
	if (l.isEmpty())
	    return makeEmptyString();
	Iterator i = l.iterator();
	Automaton b = ((Automaton) i.next()).cloneExpanded();
	Iterator j = b.getAcceptStates().iterator();
	while (i.hasNext()) {
	    Automaton a = ((Automaton) i.next()).cloneExpanded();
	    Set ns = a.getAcceptStates();
	    while (j.hasNext()) {
		State s = (State) j.next();
		s.accept = false;
		s.addEpsilon(a.initial);
		if (s.accept)
		    ns.add(s);
	    }
	    j = ns.iterator();
	}
	b.deterministic = false;
	b.checkMinimizeAlways();
	return b;
    }

    /** Returns new automaton that accepts the union of the empty string and
	the language of this automaton. 
	<p>
	Complexity: linear in number of states. */
    public Automaton optional()
    {
	Automaton a = cloneExpanded();
	State s = new State();
	s.addEpsilon(a.initial);
        s.accept = true;
	a.initial = s;
	a.deterministic = false;
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns new automaton that accepts the Kleene star (zero or
	more concatenated repetitions) of the language of this automaton. 
	<p>
	Complexity: linear in number of states. */
    public Automaton repeat()
    {
	Automaton a = cloneExpanded();
	State s = new State();
        s.accept = true;
	s.addEpsilon(a.initial);
	Iterator i = a.getAcceptStates().iterator();
	while (i.hasNext()) {
	    State p = (State) i.next();
	    p.addEpsilon(s);
	}
	a.initial = s;
	a.deterministic = false;
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns new automaton that accepts <code>min</code> or more
	concatenated repetitions of the language of this automaton.  
	<p>
	Complexity: linear in number of states and in <code>min</code>. */
    public Automaton repeat(int min)
    {
	Automaton a = repeat();
	while (min-->0)
	    a = concatenate(a);
	return a;
    }

    /** Returns new automaton that accepts between <code>min</code> and <code>max</code>
	(including both) concatenated repetitions of the language of this automaton.  
	<p>
	Complexity: linear in number of states and in <code>min</code>
	and <code>max</code>.  */
    public Automaton repeat(int min, int max)
    {
	expandSingleton();
	if (min>max)
	    return makeEmpty();
	max -= min;
	Automaton a;
	if (min==0)
	    a = makeEmptyString();
	else if (min==1)
	    a = (Automaton) clone();
	else {
	    a = this;
	    while (--min>0)
		a = concatenate(a);
	}
	if (max==0)
	    return a;
	Automaton d = (Automaton) clone();
	while (--max>0) {
	    Automaton c = (Automaton) clone();
	    Iterator i = c.getAcceptStates().iterator();
	    while (i.hasNext()) {
		State p = (State) i.next();
		p.addEpsilon(d.initial);
	    }
	    d = c;
	}
	Iterator i = a.getAcceptStates().iterator();
	while (i.hasNext()) {
	    State p = (State) i.next();
	    p.addEpsilon(d.initial);
	}
	a.deterministic = false;
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns new (deterministic) automaton that accepts the 
	complement of the language of this automaton. 
	<p>
	Complexity: linear in number of states (if already deterministic). */
    public Automaton complement()
    {
	Automaton a = cloneExpanded();
	a.determinize();
	a.totalize();
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State p = (State) i.next();
	    p.accept = !p.accept;
	}
	a.removeDeadTransitions();
	return a;
    }

    /** Returns new (deterministic) automaton that 
	accepts the intersection of the languages of this and the given automaton.
	As a side-effect, both this and the given automaton are determinized, 
	if not already deterministic. 
	<p>
	Complexity: quadratic in number of states (if already deterministic). */
    public Automaton intersection(Automaton a)
    {
	if (isSingleton() && a.isSingleton()) {
	    if (singleton.equals(a.singleton))
		return makeString(singleton);
	    else
		return makeEmpty();
	}
	determinize();
	a.determinize();
	Transition[][] transitions1 = getSortedTransitions(getStates());
	Transition[][] transitions2 = getSortedTransitions(a.getStates());
	Automaton c = new Automaton();
	LinkedList worklist = new LinkedList();
	HashMap newstates = new HashMap();
	State s = new State();
	c.initial = s;
	StatePair p = new StatePair(s, initial, a.initial);
	worklist.add(p);
	newstates.put(p, p);
	while (worklist.size()>0) {
	    p = (StatePair) worklist.removeFirst();
	    p.s.accept = p.s1.accept && p.s2.accept;
	    Transition[] t1 = transitions1[p.s1.number];
	    Transition[] t2 = transitions2[p.s2.number];
	    for (int n1 = 0, n2 = 0; n1<t1.length && n2<t2.length;) {
		if (t1[n1].max<t2[n2].min)
		    n1++;
		else if (t2[n2].max<t1[n1].min)
		    n2++;
		else {
		    StatePair q = new StatePair(t1[n1].to, t2[n2].to);
		    StatePair r = (StatePair) newstates.get(q);
		    if (r==null) {
		        q.s = new State();
			worklist.add(q);
			newstates.put(q, q);
			r = q;
		    }
		    char min = t1[n1].min>t2[n2].min ? t1[n1].min : t2[n2].min;
		    char max = t1[n1].max<t2[n2].max ? t1[n1].max : t2[n2].max;
		    p.s.transitions.add(new Transition(min, max, r.s));
		    if (t1[n1].max<t2[n2].max)
			n1++;
		    else
			n2++;
		}
	    }
	}
	c.deterministic = true;
	c.removeDeadTransitions();
	c.checkMinimizeAlways();
	return c;
    }

    /** Returns new automaton that accepts the union of the languages
	of this and the given automaton. 
	<p>
	Complexity: linear in number of states. */
    public Automaton union(Automaton a)
    {
	if (isSingleton() && a.isSingleton() && singleton.equals(a.singleton))
	    return (Automaton) clone();
	a = a.cloneExpanded();
	Automaton b = cloneExpanded();
	State s = new State();
	s.addEpsilon(a.initial);
	s.addEpsilon(b.initial);
	a.initial = s;
	a.deterministic = false;
	a.checkMinimizeAlways();
	return a;
    }
    
    /** Returns new automaton that accepts the union of the languages
	of the given automata.
	<p>
	Complexity: linear in number of states. */
    static public Automaton union(List l)
    {
	State s = new State();
	Iterator i = l.iterator();
	while (i.hasNext()) {
	    Automaton b = ((Automaton) i.next()).cloneExpanded();
	    s.addEpsilon(b.initial);
	}
	Automaton a = new Automaton();
	a.initial = s;
	a.deterministic = false;
	a.checkMinimizeAlways();
	return a;
    }

    /** Determinizes this automaton. 
	<p>
	Complexity: exponential in number of states. */
    public void determinize()
    {
	if (deterministic || isSingleton())
	    return;
	Set initialset = new HashSet();
	initialset.add(initial);
	determinize(initialset);
    }

    /** Determinizes this automaton using the given set of initial states. */
    private void determinize(Set initialset)
    {
	char[] points = getStartPoints();
	// subset construction
	Map sets = new HashMap();
	LinkedList worklist = new LinkedList();
	Map newstate = new HashMap();
	sets.put(initialset, initialset);
	worklist.add(initialset);
	initial = new State();
	newstate.put(initialset, initial);
	while (worklist.size()>0) {
	    Set s = (Set) worklist.removeFirst();
	    State r = (State) newstate.get(s);
	    Iterator i = s.iterator();
	    while (i.hasNext()) {
		State q = (State) i.next();
		if (q.accept) {
		    r.accept = true;
		    break;
		}
	    }
	    for (int n = 0; n<points.length; n++) {
		Set p = new HashSet();
		Iterator j = s.iterator();
		while (j.hasNext()) {
		    State q = (State) j.next();
		    Iterator k = q.transitions.iterator();
		    while (k.hasNext()) {
			Transition t = (Transition) k.next();
			if (t.min<=points[n] && points[n]<=t.max) 
			    p.add(t.to);
		    }
		}
		if (!sets.containsKey(p)) {
		    sets.put(p, p);
		    worklist.add(p);
		    newstate.put(p, new State());
		}
		State q = (State) newstate.get(p);
		char min = points[n];
		char max;
		if (n+1<points.length)
		    max = (char) (points[n+1]-1);
		else
		    max = Character.MAX_VALUE;
		r.transitions.add(new Transition(min, max, q));
	    }
	}
	deterministic = true;
	removeDeadTransitions();
    }

    /**
     * Minimizes (and determinizes if not already deterministic) this automaton.
     * @see #setMinimization(int)
     */
    public void minimize()
    {
	if (!isSingleton()) {
	    switch (minimization) {
	    case MINIMIZE_HUFFMAN:
		minimizeHuffman();
		break;
	    case MINIMIZE_BRZOZOWSKI:
		minimizeBrzozowski();
		break;
	    default:
		minimizeHopcroft();
	    }
	}
	// recompute hash code
	hash_code = getNumberOfStates()*3+getNumberOfTransitions()*2;
	if (hash_code==0)
	    hash_code = 1;
    }

    private boolean statesAgree(Transition[][] transitions, boolean[][] mark, int n1, int n2)
    {
	Transition[] t1 = transitions[n1];
	Transition[] t2 = transitions[n2];
	for (int k1 = 0, k2 = 0; k1<t1.length && k2<t2.length;) {
	    if (t1[k1].max<t2[k2].min)
		k1++;
	    else if (t2[k2].max<t1[k1].min)
		k2++;
	    else {
		int m1 = t1[k1].to.number;
		int m2 = t2[k2].to.number;
		if (m1>m2) {
		    int t = m1;
		    m1 = m2;
		    m2 = t;
		}
		if (mark[m1][m2])
		    return false;
		if (t1[k1].max<t2[k2].max)
		    k1++;
		else
		    k2++;
	    }
	}
	return true;
    }

    private void addTriggers(Transition[][] transitions, boolean[][] mark, 
			     HashSet[][] triggers, int n1, int n2)
    {
	Transition[] t1 = transitions[n1];
	Transition[] t2 = transitions[n2];
	for (int k1 = 0, k2 = 0; k1<t1.length && k2<t2.length;) {
	    if (t1[k1].max<t2[k2].min)
		k1++;
	    else if (t2[k2].max<t1[k1].min)
		k2++;
	    else {
		if (t1[k1].to!=t2[k2].to) {
		    int m1 = t1[k1].to.number;
		    int m2 = t2[k2].to.number;
		    if (m1>m2) {
			int t = m1;
			m1 = m2;
			m2 = t;
		    }
		    if (triggers[m1][m2]==null)
			triggers[m1][m2] = new HashSet();
		    triggers[m1][m2].add(new IntPair(n1, n2));
		}
		if (t1[k1].max<t2[k2].max)
		    k1++;
		else
		    k2++;
	    }
	}
    }
	    
    private void markPair(boolean[][] mark, HashSet[][] triggers, int n1, int n2)
    {
	mark[n1][n2] = true;
	if (triggers[n1][n2]!=null) {
	    Iterator i = triggers[n1][n2].iterator();
	    while (i.hasNext()) {
		IntPair p = (IntPair) i.next();
		int m1 = p.n1;
		int m2 = p.n2;
		if (m1>m2) {
		    int t = m1;
		    m1 = m2;
		    m2 = t;
		}
		if (!mark[m1][m2])
		    markPair(mark, triggers, m1, m2);
	    }
	}
    }

    /** Minimize using Huffman's algorithm. */
    private void minimizeHuffman()
    {
	determinize();
	totalize();
	Set ss = getStates();
	Transition[][] transitions = new Transition[ss.size()][];
	State[] states = (State[]) ss.toArray(new State[0]);
	boolean[][] mark = new boolean[states.length][states.length];
	HashSet[][] triggers = new HashSet[states.length][states.length];
	// initialize marks based on acceptance status and find transition arrays
	for (int n1 = 0; n1<states.length; n1++) {
	    states[n1].number = n1;
	    transitions[n1] = states[n1].getSortedTransitionArray(false);
	    for (int n2 = n1+1; n2<states.length; n2++)
		if (states[n1].accept!=states[n2].accept)
		    mark[n1][n2] = true;
	}
	// for all pairs, see if states agree
	for (int n1 = 0; n1<states.length; n1++)
	    for (int n2 = n1+1; n2<states.length; n2++)
		if (!mark[n1][n2]) {
		    if (statesAgree(transitions, mark, n1, n2))
			addTriggers(transitions, mark, triggers, n1, n2);
		    else
			markPair(mark, triggers, n1, n2);
		}
	// assign equivalence class numbers to states
	int numclasses = 0;
	for (int n = 0; n<states.length; n++)
	    states[n].number = -1;
	for (int n1 = 0; n1<states.length; n1++)
	    if (states[n1].number==-1) {
		states[n1].number = numclasses;
		for (int n2 = n1+1; n2<states.length; n2++)
		    if (!mark[n1][n2])
			states[n2].number = numclasses;
		numclasses++;
	    }
	// make a new state for each equivalence class
	State[] newstates = new State[numclasses];
	for (int n = 0; n<numclasses; n++)
	    newstates[n] = new State();
	// select a class representative for each class and find the new initial state
	for (int n = 0; n<states.length; n++) {
	    newstates[states[n].number].number = n;
	    if (states[n]==initial)
		initial = newstates[states[n].number];
	}
	// build transitions and set acceptance
	for (int n = 0; n<numclasses; n++) {
	    State s = newstates[n];
	    s.accept = states[s.number].accept;
	    Iterator i = states[s.number].transitions.iterator();
	    while (i.hasNext()) {
		Transition t = (Transition) i.next();
		s.transitions.add(new Transition(t.min, t.max, newstates[t.to.number]));
	    }
	}
	removeDeadTransitions();
    }

    /** Minimize using Brzozowski's algorithm. */
    private void minimizeBrzozowski()
    {
	if (isSingleton())
	    return;
	determinize(reverse());
	determinize(reverse());
    }

    /** Minimize using Hopcroft's algorithm. */
    private void minimizeHopcroft()
    {
	determinize();
	Set tr = initial.getTransitions();
	if (tr.size()==1) {
	    Transition t = (Transition) tr.iterator().next();
	    if (t.to==initial && t.min==Character.MIN_VALUE && t.max==Character.MAX_VALUE)
		return;
	}
       	totalize();
	// make arrays for numbered states and effective alphabet
	Set ss = getStates();
	State[] states = new State[ss.size()];
	Iterator i = ss.iterator();
	int number = 0;
	while (i.hasNext()) {
	    State q = (State) i.next();
	    states[number] = q;
	    q.number = number++;
	}
	char[] sigma = getStartPoints();
	// initialize data structures
	LinkedList[][] reverse = new LinkedList[states.length][sigma.length];
	boolean[][] reverse_nonempty = new boolean[states.length][sigma.length];
	LinkedList[] partition = new LinkedList[states.length];
	int[] block = new int[states.length];
	StateList[][] active = new StateList[states.length][sigma.length];
	StateListNode[][] active2 = new StateListNode[states.length][sigma.length];
	LinkedList pending = new LinkedList();
	boolean[][] pending2 = new boolean[sigma.length][states.length];
	ArrayList split = new ArrayList();
	boolean[] split2 = new boolean[states.length];
	ArrayList refine = new ArrayList();
	boolean[] refine2 = new boolean[states.length];
	ArrayList[] splitblock = new ArrayList[states.length];
	for (int q = 0; q<states.length; q++) {
	    splitblock[q] = new ArrayList();
	    partition[q] = new LinkedList();
	    for (int a = 0; a<sigma.length; a++) {
		reverse[q][a] = new LinkedList();
		active[q][a] = new StateList();
	    }
	}
	// find initial partition and reverse edges
	for (int q = 0; q<states.length; q++) {
	    State qq = states[q];
	    int j;
	    if (qq.accept) 
		j = 0;
	    else
		j = 1;
	    partition[j].add(qq);
	    block[qq.number] = j;
	    for (int a = 0; a<sigma.length; a++) {
		char aa = sigma[a];
		State p = qq.step(aa);
		reverse[p.number][a].add(qq);
		reverse_nonempty[p.number][a] = true;
	    }
	}	
	// initialize active sets
	for (int j = 0; j<=1; j++) 
	    for (int a = 0; a<sigma.length; a++) {
		Iterator q = partition[j].iterator();
		while (q.hasNext()) {
		    State qq = (State) q.next();
		    if (reverse_nonempty[qq.number][a])
			active2[qq.number][a] = active[j][a].add(qq);
		}
	    }
	// initialize pending
	for (int a = 0; a<sigma.length; a++) {
	    char aa = sigma[a];
	    int a0 = active[0][a].size;
	    int a1 = active[1][a].size;
	    int j;
	    if (a0 <= a1) 
		j = 0;
	    else
		j = 1;
	    pending.add(new IntPair(j, a));
	    pending2[a][j] = true;
	}
	// process pending until fixed point
	int k = 2;
	while (!pending.isEmpty()) {
	    IntPair ip = (IntPair) pending.removeFirst();
	    int p = ip.n1;
	    int a = ip.n2;
	    pending2[a][p] = false;
	    // find states that need to be split off their blocks
	    for (StateListNode m = active[p][a].first; m!=null; m = m.next) {
		Iterator r = reverse[m.q.number][a].iterator();
		while (r.hasNext()) {
		    State s = (State) r.next();
		    if (!split2[s.number]) {
			split2[s.number] = true;
			split.add(s);
			int j = block[s.number];
			splitblock[j].add(s);
			if (!refine2[j]) {
			    refine2[j] = true;
			    refine.add(new Integer(j));
			}
		    }
		}
	    }  
	    // refine blocks
	    Iterator m = refine.iterator();
	    while (m.hasNext()) {
		int j = ((Integer) m.next()).intValue();
		if (splitblock[j].size()<partition[j].size()) {
		    LinkedList b1 = partition[j];
		    LinkedList b2 = partition[k];
		    Iterator r = splitblock[j].iterator();
		    while (r.hasNext()) {
			State s = (State) r.next();
			b1.remove(s);
			b2.add(s);
			block[s.number] = k;
			for (int c = 0; c<sigma.length; c++) {
			    StateListNode sn = active2[s.number][c];
			    if (sn!=null && sn.sl==active[j][c]) {
				sn.remove();
				active2[s.number][c] = active[k][c].add(s);
			    }
			}
		    }
		    // update pending
		    for (int c = 0; c<sigma.length; c++) {
			int aj = active[j][c].size;
			int ak = active[k][c].size;
			if (!pending2[c][j] && 0<aj && aj<=ak) {
			    pending2[c][j] = true;
			    pending.add(new IntPair(j, c));
			} else {
			    pending2[c][k] = true;
			    pending.add(new IntPair(k, c));
			}
		    }
		    k++;
		}
		Iterator r = splitblock[j].iterator();
		while (r.hasNext()) {
		    State s = (State) r.next();
		    split2[s.number] = false;
		}
		refine2[j] = false;
		splitblock[j].clear();
	    }
	    split.clear();
	    refine.clear();
	}
	// make a new state for each equivalence class, set initial state
	State[] newstates = new State[k];
	for (int n = 0; n<newstates.length; n++) {
	    State s = new State();
	    newstates[n] = s;
	    Iterator j = partition[n].iterator();
	    while (j.hasNext()) {
		State q = (State) j.next();
		if (q==initial)
		    initial = s;
		s.accept = q.accept;
		s.number = q.number; // select representative
		q.number = n;
	    }
	}
	// build transitions and set acceptance
	for (int n = 0; n<newstates.length; n++) {
	    State s = newstates[n];
	    s.accept = states[s.number].accept;
	    Iterator j = states[s.number].transitions.iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		s.transitions.add(new Transition(t.min, t.max, newstates[t.to.number]));
	    }
	}
	removeDeadTransitions();
    }

    /** 
     * Reverses the language of this (non-singleton) automaton 
     * while returning the set of new initial states. 
     * (Used for <code>minimizeBrzozowski()</code>) 
     */
    private Set reverse()
    {
	// reverse all edges
	HashMap m = new HashMap();
	Set states = getStates();
	Set accept = getAcceptStates();
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State r = (State) i.next();
	    m.put(r, new HashSet());
	    r.accept = false;
	}
	i = states.iterator();
	while (i.hasNext()) {
	    State r = (State) i.next();
	    Iterator j = r.getTransitions().iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();		    
		((HashSet) m.get(t.to)).add(new Transition(t.min, t.max, r));
	    }
	}
	i = states.iterator();
	while (i.hasNext()) {
	    State r = (State) i.next();
	    r.transitions = (HashSet) m.get(r);
	}
	// make new initial+final states
	initial.accept = true;
	initial = new State();
	i = accept.iterator();
	while (i.hasNext()) {
	    State r = (State) i.next();
	    initial.addEpsilon(r); // ensures that all initial states are reachable
	}
	deterministic = false;
	return accept;
    }

    /** Returns new automaton that accepts the single chars that occur
	in strings that are accepted by this automaton. */
    public Automaton singleChars()
    {
	Automaton a = new Automaton();
	State s = new State();
	a.initial = s;
	State q = new State();
	q.accept = true; 
	if (isSingleton()) {
	    for (int i = 0; i<singleton.length(); i++)
		s.transitions.add(new Transition(singleton.charAt(i), q));		
	} else {
	    Iterator i = getStates().iterator();
	    while (i.hasNext()) {
		State p = (State) i.next();
		Iterator j = p.transitions.iterator();
		while (j.hasNext()) {
		    Transition t = (Transition) j.next();
		    s.transitions.add(new Transition(t.min, t.max, q));
		}
	    }
	}
	a.deterministic = true;
	a.removeDeadTransitions();
	return a;
    }

    private void addSetTransitions(State s, String set, State p)
    {
	for (int n = 0; n<set.length(); n++)
	    s.transitions.add(new Transition(set.charAt(n), p));
    }

    /** Returns a new automaton that accepts the trimmed language of this automaton.
	The resulting automaton is constructed as follows:
	1) Whenever a <code>c</code> character is allowed in the original automaton,
	one or more <code>set</code> characters are allowed in the new automaton.
	2) The automaton is prefixed and postfixed with any number of 
	<code>set</code> characters.
	@param set set of characters to be trimmed
	@param c canonical trim character (assumed to be in <code>set</code>) */
    public Automaton trim(String set, char c)
    {
	Automaton a = cloneExpanded();
	State f = new State();
	addSetTransitions(f, set, f);
	f.accept = true;
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    State r = s.step(c);
	    if (r!=null) {
		// add inner
		State q = new State();
		addSetTransitions(q, set, q);
		addSetTransitions(s, set, q);
		q.addEpsilon(r);
	    }
	    // add postfix
	    if (s.accept)
		s.addEpsilon(f);
	}
	// add prefix
	State p = new State();
	addSetTransitions(p, set, p);
	p.addEpsilon(a.initial);
	a.initial = p;
	a.deterministic = false;
	a.removeDeadTransitions();
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns a new automaton that accepts the compressed language of this automaton.
	Whenever a <code>c</code> character is allowed in the original automaton,
	one or more <code>set</code> characters are allowed in the new automaton.
	@param set set of characters to be compressed
	@param c canonical compress character (assumed to be in <code>set</code>) */
    public Automaton compress(String set, char c)
    {
	Automaton a = cloneExpanded();
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    State r = s.step(c);
	    if (r!=null) {
		// add inner
		State q = new State();
		addSetTransitions(q, set, q);
		addSetTransitions(s, set, q);
		q.addEpsilon(r);
	    }
	}
	// add prefix
	a.deterministic = false;
	a.removeDeadTransitions();
	a.checkMinimizeAlways();
	return a;
    }

    /** finds the largest entry whose value is less than or equal to c, or 0 if there is no such entry */
    static int findIndex(char c, char[] points)
    {
	int a = 0;
	int b = points.length;
	while (b-a>1) {
	    int d = (a+b)/2;
	    if (points[d]>c)
		b = d;
	    else if (points[d]<c)
		a = d;
	    else
		return d;
	}
	return a;
    }

    /** Returns new automaton where all transition labels have been substituted.
	<p>
	Each transition labeled <code>c</code> is changed to a set of transitions, one for each character
	in <code>map(c)</code>. If <code>map(c)</code> is null, then the transition is unchanged.
	@param map map from characters to sets of characters (where characters are <code>Character</code> objects) */
    public Automaton subst(Map map)
    {
	if (map.isEmpty())
	    return (Automaton) clone();
	Set ckeys = new TreeSet(map.keySet());
	char[] keys = new char[ckeys.size()];
	Iterator k = ckeys.iterator();
	for (int j = 0; j<keys.length; j++)
	    keys[j] = ((Character) k.next()).charValue();
      	Automaton a = cloneExpanded();
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    Iterator j = s.transitions.iterator();
	    s.resetTransitions();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		int index = findIndex(t.min, keys);
		while (t.min<=t.max) {
		    if (keys[index]>t.min) {
			char m = (char) (keys[index]-1);
			s.transitions.add(new Transition(t.min, m, t.to));
			t.min = (char) (m+1);
		    } else if (keys[index]<t.min) {
			char m;
			if (index+1<keys.length) 
			    m = (char) (keys[++index]-1);
			else
			    m = Character.MAX_VALUE;
			if (t.max<m)
			    m = t.max;
			s.transitions.add(new Transition(t.min, m, t.to));
			if (m+1 > Character.MAX_VALUE)
			    break;
			t.min = (char) (m+1);
		    } else { // found t.min in substitution map
			Set ss = (Set) map.get(new Character(t.min));
			k = ss.iterator();
			while (k.hasNext()) {
			    char c = ((Character) k.next()).charValue();
			    s.transitions.add(new Transition(c, t.to));
			}
			if (t.min+1 > Character.MAX_VALUE)
			    break;
			t.min++;
			if (index+1<keys.length && keys[index+1]==t.min)
			    index++;
		    }
		}
	    }
	}
	a.deterministic = false;
	a.removeDeadTransitions();
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns new automaton accepting the homomorphic image of this automaton
	using the given function.
	<p>
	This method maps each transition label to a new value.  
	<code>source</code> and <code>dest</code> are assumed to be
	arrays of same length, and <code>source</code> must be sorted in
	increasing order and contain no duplicates. 
	<code>source</code> defines the starting points of char intervals,
	and the corresponding entries in <code>dest</code> define the 
	starting points of corresponding new intervals.  */
    public Automaton homomorph(char[] source, char[] dest)
    {
	Automaton a = cloneExpanded();
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    Iterator j = s.transitions.iterator();
	    s.resetTransitions();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		int min = t.min;
		while (min<=t.max) {
		    int n = findIndex((char) min, source);
		    char nmin = (char) (dest[n]+min-source[n]);
		    int end = (n+1==source.length) ? Character.MAX_VALUE : source[n+1]-1;
		    int length;
		    if (end<t.max)
			length = end+1 - min;
		    else
			length = t.max+1 - min;
		    s.transitions.add(new Transition(nmin, (char) (nmin+length-1), t.to));
		    min += length;
		}
	    }
	}
	a.deterministic = false;
	a.removeDeadTransitions();
	a.checkMinimizeAlways();
	return a;
    }

    /** Returns new automaton with projected alphabet.  
        The new automaton accepts all strings that are projections of
        strings accepted by this automaton onto the given characters
        (represented by <code>Character</code>). If <code>null</code> is
        in the set, it abbreviates the intervals u0000-uDFFF 
	and uF900-uFFFF (i.e., the non-private code points). It is assumed that
	all other characters from <code>chars</code> are in the interval
        uE000-uF8FF. */
    public Automaton projectChars(Set chars)
    {
	Character[] c = (Character[]) chars.toArray(new Character[0]);
	char[] cc = new char[c.length];
	boolean normalchars = false;
	for (int i = 0; i<c.length; i++)
	    if (c[i]==null)
		normalchars = true;
	    else
		cc[i] = c[i].charValue();
	Arrays.sort(cc);
	HashSet epsilons = new HashSet();
	Automaton a = cloneExpanded();
	Iterator i = a.getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    HashSet new_transitions = new HashSet();
	    Iterator j = s.transitions.iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		boolean addepsilon = false;
		if (t.min<'\uf900' && t.max>'\udfff') {
		    int w1 = Arrays.binarySearch(cc, t.min>'\ue000' ? t.min : '\ue000');
		    if (w1<0) {
			w1 = -w1-1;
			addepsilon = true;
		    }
		    int w2 = Arrays.binarySearch(cc, t.max<'\uf8ff' ? t.max : '\uf8ff');
		    if (w2<0) {
			w2 = -w2-2;
			addepsilon = true;
		    }
		    for (int w = w1; w<=w2; w++) {
			new_transitions.add(new Transition(cc[w], t.to));
			if (w>w1 && cc[w-1]+1!=cc[w])
			    addepsilon = true;
		    }
		}
		if (normalchars) {
		    if (t.min<='\udfff')
			new_transitions.add(new Transition(t.min, 
							   t.max<'\udfff' ? t.max : '\udfff', 
							   t.to));
		    if (t.max>='\uf900')
			new_transitions.add(new Transition(t.min>'\uf900' ? t.min : '\uf900', 
							   t.max,
							   t.to));
		} else if (t.min<='\udfff' || t.max>='\uf900')
		    addepsilon = true;
		if (addepsilon)
		    epsilons.add(new StatePair(s, t.to));
	    }
	    s.transitions = new_transitions;
	}
	a.reduce();
	a.addEpsilons(epsilons);
	a.removeDeadTransitions();
	a.checkMinimizeAlways();
	return a;
    }

    /** 
     * Adds epsilon transitions to this automaton.
     * This method adds extra character interval transitions that are equivalent to the given
     * set of epsilon transitions. 
     * @param pairs collection of {@link StatePair} objects representing pairs of source/destination states 
     *        where epsilon transitions should be added
     */
    public void addEpsilons(Collection pairs)
    {
	expandSingleton();
	HashMap forward = new HashMap();
	HashMap back = new HashMap();
	Iterator i = pairs.iterator();
	while (i.hasNext()) {
	    StatePair p = (StatePair) i.next();
	    HashSet to = (HashSet) forward.get(p.s1);
	    if (to==null) {
		to = new HashSet();
		forward.put(p.s1, to);
	    }
	    to.add(p.s2);
	    HashSet from = (HashSet) back.get(p.s2);
	    if (from==null) {
		from = new HashSet();
		back.put(p.s2, from);
	    }
	    from.add(p.s1);
	}
	// calculate epsilon closure
	LinkedList worklist = new LinkedList(pairs);
	HashSet workset = new HashSet(pairs);
	while (!worklist.isEmpty()) {
	    StatePair p = (StatePair) worklist.removeFirst();
	    workset.remove(p);
	    HashSet to = (HashSet) forward.get(p.s2);
	    HashSet from = (HashSet) back.get(p.s1);
	    if (to!=null) {
		Iterator j = to.iterator();
		while (j.hasNext()) {
		    State s = (State) j.next();
		    StatePair pp = new StatePair(p.s1, s);
		    if (!pairs.contains(pp)) {
			pairs.add(pp);
			((HashSet) forward.get(p.s1)).add(s);
			((HashSet) back.get(s)).add(p.s1);
			worklist.add(pp);
			workset.add(pp);
			if (from!=null) {
			    Iterator k = from.iterator();
			    while (k.hasNext()) {
				State q = (State) k.next();
				StatePair qq = new StatePair(q, p.s1);
				if (!workset.contains(qq)) {
				    worklist.add(qq);
				    workset.add(qq);
				}
			    }
			}
		    }
		}
	    }
	}
	// add transitions
	i = pairs.iterator();
	while (i.hasNext()) {
	    StatePair p = (StatePair) i.next();
	    p.s1.addEpsilon(p.s2);
	}
	deterministic = false;
	checkMinimizeAlways();
    }

    /** Returns true if the given string is accepted by this automaton. 
	As a side-effect, this automaton is determinized
	if not already deterministic. 
	<p>
	Complexity: linear in length of string (if automaton is already deterministic)
	and in number of transitions. 
	<p>
	<b>Note:</b> to obtain maximum speed, use the {@link RunAutomaton} class. */
    public boolean run(String s)
    {
	if (isSingleton())
	    return s.equals(singleton);
	determinize();
	State p = initial;
	for (int i = 0; i < s.length(); i++) {
	    State q = p.step(s.charAt(i));
	    if (q==null)
		return false;
	    p = q;
	}
	return p.accept;
    }

    /** Returns number of states in this automaton. */
    public int getNumberOfStates()
    {
	if (isSingleton())
	    return singleton.length()+1;
	return getStates().size();
    }

    /** Returns number of transitions in this automaton. 
	This number is counted as the total number of edges, where one edge
	may be a character interval. */
    public int getNumberOfTransitions()
    {
	if (isSingleton())
	    return singleton.length();
	int c = 0;
        Iterator i = getStates().iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    c += s.transitions.size();
	}
	return c;
    }

    /** Returns true if this automaton accepts no strings. */
    public boolean isEmpty()
    {
	if (isSingleton())
	    return false;
	return initial.accept==false && initial.transitions.size()==0;
    }

    /** Returns true if this automaton accepts all strings. */
    public boolean isTotal()
    {
	if (isSingleton())
	    return false;
	if (initial.accept==true && initial.transitions.size()==1) {
	    Transition t = (Transition) initial.transitions.iterator().next();
	    return t.to==initial && t.min==Character.MIN_VALUE && t.max==Character.MAX_VALUE;
	}
	return false;
    }

    /** Returns true if the language of this automaton is finite. */
    public boolean isFinite()
    {
	if (isSingleton())
	    return true;
	return isFinite(initial, new HashSet());
    }

    /** Returns set of accepted strings, assuming this automaton has a finite language.
	If the language is not finite, null is returned. */
    public Set getFiniteStrings()
    {
	HashSet strings = new HashSet();
	if (isSingleton())
	    strings.add(singleton);
	else if (!getFiniteStrings(initial, new HashSet(), strings, new StringBuffer(), -1))
	    return null;
	return strings;
    }

    /** Returns set of accepted strings, assuming that at most <code>limit</code> strings are accepted.
	If more than <code>limit</code> strings are accepted, null is returned. 
	If <code>limit</code>&lt;0, then this methods works like {@link #getFiniteStrings()}.
    */
    public Set getFiniteStrings(int limit)
    {
	HashSet strings = new HashSet();
	if (isSingleton()) {
	    if (limit>0)
		strings.add(singleton);
	    else
		return null;
	} else if (!getFiniteStrings(initial, new HashSet(), strings, new StringBuffer(), limit))
	    return null;
	return strings;
    }

    /** 
     * Returns a shortest accepted/rejected string.
     * If more than one string is found, the lexicographically first is returned.
     * @param accepted if true, look for accepted strings; otherwise, look for rejected strings
     * @return the string, null if none found
     */
    public String getShortestExample(boolean accepted)
    {
	if (isSingleton()) {
	    if (accepted)
		return singleton;
	    else if (singleton.length()>0)
		return "";
	    else
		return "\u0000";
		
	}
	return getShortestExample(initial, accepted, new HashMap());
    }

    String getShortestExample(State s, boolean accepted, Map map)
    {
	if (s.accept==accepted)
	    return "";
	if (map.containsKey(s))
	    return (String) map.get(s);
	map.put(s, null);
	String best = null;
	Iterator i = s.transitions.iterator();
	while (i.hasNext()) {
	    Transition t = (Transition) i.next();
	    String b = getShortestExample(t.to, accepted, map);
	    if (b!=null) {
		b = t.min+b;
		if (best==null || 
		    b.length()<best.length() || 
		    (b.length()==best.length() && b.compareTo(best)<0))
		    best = b;
	    }
	}
	map.put(s, best);
	return best;
    }

    /**
     * Returns the longest string that is a prefix of all accepted strings and 
     * visits each state at most once.
     * @return common prefix
     */
    public String getCommonPrefix()
    {
	if (isSingleton()) 
	    return singleton;
	StringBuffer b = new StringBuffer();
	HashSet visited = new HashSet();
	State s = initial;
	boolean done;
	do {
	    done = true;
	    visited.add(s);
	    if (!s.accept && s.transitions.size()==1) {
		Transition t = (Transition) s.transitions.iterator().next();
		if (t.min==t.max && !visited.contains(t.to)) {
		    b.append(t.min);
		    s = t.to;
		    done = false;
		}
	    }
	} while (!done);
	return b.toString();
    }

    /** Returns true if the language of this automaton is a subset of the
	language of the given automaton. 
	Equivalent to <code>this.intersection(a.complement()).isEmpty()</code>. */
    public boolean subsetOf(Automaton a)
    {
	if (isSingleton()) {
	    if (a.isSingleton())
		return singleton.equals(a.singleton);
	    else
		return a.run(singleton);
	}
	return intersection(a.complement()).isEmpty();
    }

    /** Returns true if the language of this automaton is equal to the
	language of the given automaton. 
	Implemented using <code>hashCode</code> and <code>subsetOf</code>. */
    public boolean equals(Object obj)
    {
	if (!(obj instanceof Automaton))
	    return false;
	Automaton a = (Automaton) obj;
	if (isSingleton() && a.isSingleton())
	    return singleton.equals(a.singleton);
	return hashCode()==a.hashCode() && subsetOf(a) && a.subsetOf(this);
    }

    /** 
     * Returns new automaton that accepts the shuffle (interleaving) of 
     * the languages of this and the given automaton.
     * As a side-effect, both this and the given automaton are determinized,
     * if not already deterministic.     
     * <p>
     * Complexity: quadratic in number of states (if already deterministic). 
     * <p>
     * <dl><dt><b>Author:</b></dt><dd>Torben Ruby 
     * &lt;<a href="mailto:ruby@daimi.au.dk">ruby@daimi.au.dk</a>&gt;</dd></dl>
     */
    public Automaton shuffle(Automaton a)
    {
	determinize();
	a.determinize();
	Transition[][] transitions1 = getSortedTransitions(getStates());
	Transition[][] transitions2 = getSortedTransitions(a.getStates());
	Automaton c = new Automaton();
	LinkedList worklist = new LinkedList();
	HashMap newstates = new HashMap();
	State s = new State();
	c.initial = s;
	StatePair p = new StatePair(s, initial, a.initial);
	worklist.add(p);
	newstates.put(p, p);
	while (worklist.size()>0) {
	    p = (StatePair) worklist.removeFirst();
	    p.s.accept = p.s1.accept && p.s2.accept;
	    Transition[] t1 = transitions1[p.s1.number];
	    for (int n1 = 0; n1<t1.length; n1++) {
		StatePair q = new StatePair(t1[n1].to, p.s2);
		StatePair r = (StatePair) newstates.get(q);
		if (r==null) {
		    q.s = new State();
		    worklist.add(q);
		    newstates.put(q, q);
		    r = q;
		}
		p.s.transitions.add(new Transition(t1[n1].min, t1[n1].max, r.s));
	    }
	    Transition[] t2 = transitions2[p.s2.number];
	    for (int n2 = 0; n2<t2.length; n2++) {
		StatePair q = new StatePair(p.s1, t2[n2].to);
		StatePair r = (StatePair) newstates.get(q);
		if (r==null) {
		    q.s = new State();
		    worklist.add(q);
		    newstates.put(q, q);
		    r = q;
		}
		p.s.transitions.add(new Transition(t2[n2].min, t2[n2].max, r.s));
	    }
	}
	c.deterministic = false;
	c.removeDeadTransitions();
	c.checkMinimizeAlways();
	return c;
    }
    
    /** Returns hash code for this automaton. 
	The hash code is based on the number of states and transitions in the minimized automaton. */
    public int hashCode()
    {
	if (hash_code==0)
	    minimize();
	return hash_code;
    }

    /** Returns a string representation of this automaton. */
    public String toString()
    {
	StringBuffer b = new StringBuffer();
	Set states = getStates();
	setStateNumbers(states);
	b.append("initial state: ").append(initial.number).append("\n");
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    b.append(s.toString());
	}
	return b.toString();
    }

    /** Returns <a href="http://www.research.att.com/sw/tools/graphviz/" 
	target="_top">Graphviz Dot</a> representation of this automaton. */
    public String toDot()
    {
	StringBuffer b = new StringBuffer("digraph Automaton {\n");
	b.append("  rankdir = LR;\n");
	Set states = getStates();
	setStateNumbers(states);
	Iterator i = states.iterator();
	while (i.hasNext()) {
	    State s = (State) i.next();
	    b.append("  ").append(s.number);
	    if (s.accept)
		b.append(" [shape=doublecircle,label=\"\"];\n");
	    else
		b.append(" [shape=circle,label=\"\"];\n");
	    if (s==initial) {
		b.append("  initial [shape=plaintext,label=\"\"];\n");
		b.append("  initial -> ").append(s.number).append("\n");
	    }
	    Iterator j = s.transitions.iterator();
	    while (j.hasNext()) {
		Transition t = (Transition) j.next();
		b.append("  ").append(s.number);
		t.appendDot(b);
	    }
	}
	return b.append("}\n").toString();
    }

    /** Returns a clone of this automaton, expands if singleton. */
    Automaton cloneExpanded()
    {
	Automaton a = (Automaton) clone();
	a.expandSingleton();
	return a;
    }

    /** Returns a clone of this automaton. */
    public Object clone()
    {
	Automaton a = new Automaton();
	if (isSingleton()) {
	    a.singleton = singleton;
	} else {
	    HashMap m = new HashMap();
	    Set states = getStates();
	    Iterator i = states.iterator();
	    while (i.hasNext())
		m.put(i.next(), new State());
	    i = states.iterator();
	    while (i.hasNext()) {
		State s = (State) i.next();
		State p = (State) m.get(s);
		p.accept = s.accept;
		if (s==initial)
		    a.initial = p;
		p.transitions = new HashSet();
		Iterator j = s.transitions.iterator();
		while (j.hasNext()) {
		    Transition t = (Transition) j.next();
		    p.transitions.add(new Transition(t.min, t.max, (State) m.get(t.to)));
		}
	    }
	}
	a.deterministic = deterministic;
	a.info = info;
	return a;
    }

    /** 
     * Retrieves a serialized <code>Automaton</code> located by a URL.
     * @param url URL of serialized automaton
     * @exception IOException if input/output related exception occurs
     * @exception OptionalDataException if the data is not a serialized object
     * @exception InvalidClassException if the class serial number does not match
     * @exception ClassCastException if the data is not a serialized <code>Automaton</code>
     * @exception ClassNotFoundException if the class of the serialized object cannot be found
     */
    public static Automaton load(URL url)
	throws IOException, OptionalDataException, ClassCastException, 
	       ClassNotFoundException, InvalidClassException
    {
	return load(url.openStream());
    }

    /** 
     * Retrieves a serialized <code>Automaton</code> from a stream.
     * @param stream input stream with serialized automaton
     * @exception IOException if input/output related exception occurs
     * @exception OptionalDataException if the data is not a serialized object
     * @exception InvalidClassException if the class serial number does not match
     * @exception ClassCastException if the data is not a serialized <code>Automaton</code>
     * @exception ClassNotFoundException if the class of the serialized object cannot be found
     */
    public static Automaton load(InputStream stream)
	throws IOException, OptionalDataException, ClassCastException, 
	       ClassNotFoundException, InvalidClassException
    {
	ObjectInputStream s = new ObjectInputStream(stream);
	return (Automaton) s.readObject();
    }

    /**
     * Writes this <code>Automaton</code> to the given stream.
     * @param stream output stream for serialized automaton
     * @exception IOException if input/output related exception occurs
     */
    public void store(OutputStream stream)
	throws IOException 
    {
	ObjectOutputStream s = new ObjectOutputStream(stream);
	s.writeObject(this);
	s.flush();
    }
}

class IntPair
{
    int n1;
    int n2;

    IntPair(int n1, int n2)
    {
	this.n1 = n1;
	this.n2 = n2;
    }
}

class StateList
{
    int size;
    StateListNode first, last;
    
    StateListNode add(State q)
    {
	return new StateListNode(q, this);
    }
}

class StateListNode
{
    State q;
    StateListNode next, prev;
    StateList sl;

    StateListNode(State q, StateList sl)
    {
	this.q = q;
	this.sl = sl;
	if (sl.size++ == 0) 
	    sl.first = sl.last = this;
	else {
	    sl.last.next = this;
	    prev = sl.last;
	    sl.last = this;
	}
    }

    void remove()
    {
	sl.size--;
	if (sl.first==this)
	    sl.first = next;
	else
	    prev.next = next;
	if (sl.last==this)
	    sl.last = prev;
	else
	    next.prev = prev;
    }
}

