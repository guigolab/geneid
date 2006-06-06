/*
 * dk.brics.automaton
 * Copyright (C) 2001-2004 Anders Moeller
 * All rights reserved.
 */

package gphase.regregex;

import java.util.*;
import java.io.*;

/** 
 * <tt>Automaton</tt> transition. 
 * <p>
 * A transition, which belongs to a source state, consists of a Unicode character interval
 * and a destination state.
 * @author Anders M&oslash;ller &lt;<a href="mailto:amoeller@brics.dk">amoeller@brics.dk</a>&gt;
 */
public class Transition implements Serializable
{ 
    static final long serialVersionUID = 40001;

    /* 
     * CLASS INVARIANT: min<=max
     */

    char min;
    char max;

    State to;

    /** 
     * Constructs new singleton interval transition. 
     * @param c transition character
     * @param to destination state
     */
    public Transition(char c, State to)
    {
	min = max = c;
	this.to = to;
    }

    /** 
     * Constructs new transition. 
     * Both end points are included in the interval.
     * @param min transition interval minimum
     * @param max transition interval maximum
     * @param to destination state
     */
    public Transition(char min, char max, State to)
    {
	if (max<min) {
	    char t = max;
	    max = min;
	    min = t;
	}
	this.min = min;
	this.max = max;
	this.to = to;
    }

    /** Returns minimum of this transition interval. */
    public char getMin()
    {
	return min;
    }

    /** Returns maximum of this transition interval. */
    public char getMax()
    {
	return max;
    }

    /** Returns destination of this transition. */
    public State getDest()
    {
	return to;
    }

    /** 
     * Checks for equality.
     * @param obj object to compare with
     * @return true if <tt>obj</tt> is a transition with same 
     *         character interval and destination state as this transition.
     */
    public boolean equals(Object obj)
    {
	if (obj instanceof Transition) {
	    Transition t = (Transition) obj;
	    return t.min==min && t.max==max && t.to==to;
	} else
	    return false;
    }
    
    /** 
     * Returns hash code.
     * The hash code is based on the character interval (not the destination state).
     * @return hash code
     */
    public int hashCode()
    {
	return min*2+max*3;
    }
    
    /** 
     * Clones this transition. 
     * @return clone with same character interval and destination state
     */
    public Object clone()
    {
	return new Transition(min, max, to);
    }

    static void appendCharString(char c, StringBuffer b)
    {
	if (c>=0x21 && c<=0x7e && c!='\\')
	    b.append(c);
	else {
	    b.append("\\u");
	    String s = Integer.toHexString((int) c);
	    if (c<0x10)
		b.append("000").append(s);
	    else if (c<0x100)
		b.append("00").append(s);
	    else if (c<0x1000)
		b.append("0").append(s);
	    else
		b.append(s);
	}
    }

    /** 
     * Returns string describing this state. Normally invoked via 
     * {@link Automaton#toString()}. 
     */
    public String toString()
    {
	StringBuffer b = new StringBuffer();
	appendCharString(min, b);
	if (min!=max){
	    b.append("-");
	    appendCharString(max, b);
	}
	b.append(" -> ").append(to.number);
	return b.toString();
    }    

    void appendDot(StringBuffer b)
    {
	b.append(" -> ").append(to.number).append(" [label=\"");
	appendCharString(min, b);
	if (min!=max){
	    b.append("-");
	    appendCharString(max, b);
	}
	b.append("\"]\n");
    }
}

class TransitionComparator implements Comparator
{
    boolean to_first;

    /** compare by (min, reverse max, to) or (to, min, reverse max) */
    public int compare(Object o1, Object o2) 
    {
	if (o1 instanceof Transition && o2 instanceof Transition) {
	    Transition t1 = (Transition) o1;
	    Transition t2 = (Transition) o2;
	    if (to_first) {
		if (t1.to.number<t2.to.number)
		    return -1;
		else if (t1.to.number>t2.to.number)
		    return 1;
	    }
	    if (t1.min<t2.min)
		return -1;
	    if (t1.min>t2.min)
		return 1;
	    if (t1.max>t2.max)
		return -1;
	    if (t1.max<t2.max)
		return 1;
	    if (!to_first) {
		if (t1.to.number<t2.to.number)
		    return -1;
		else if (t1.to.number>t2.to.number)
		    return 1;
	    }
	    return 0;
	} else
	    throw new ClassCastException();
    }
}
