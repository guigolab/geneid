/*
 * dk.brics.automaton
 * Copyright (C) 2001-2004 Anders Moeller
 * All rights reserved.
 */

package gphase.regregex;

import java.util.*;

/**
 * Regular Expression extension to <code>Automaton</code>.
 * <p>
 * Regular expressions are built from the following abstract syntax:
 * <p>
 * <table border=0>
 * <tr>
 * <td><i>regexp </i></td>
 * <td>::=</td>
 * <td><i>unionexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>unionexp </i></td>
 * <td>::=</td>
 * <td><i>interexp </i>&nbsp; <tt><b>|</b></tt> &nbsp; <i>unionexp </i>
 * </td>
 * <td>(union)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>interexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>interexp </i></td>
 * <td>::=</td>
 * <td><i>concatexp </i>&nbsp; <tt><b>&amp;</b></tt> &nbsp; <i>interexp </i>
 * </td>
 * <td>(intersection)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>concatexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>concatexp </i></td>
 * <td>::=</td>
 * <td><i>repeatexp </i>&nbsp; <i>concatexp </i></td>
 * <td>(concatenation)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>repeatexp </i></td>
 * <td>::=</td>
 * <td><i>repeatexp </i>&nbsp; <tt><b>?</b></tt></td>
 * <td>(zero or one occurrence)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i>&nbsp; <tt><b>*</b></tt></td>
 * <td>(zero or more occurrences)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i>&nbsp; <tt><b>+</b></tt></td>
 * <td>(one or more occurrences)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i>&nbsp; <tt><b>{</b><i>n</i><b>}</b></tt></td>
 * <td>(<tt><i>n</i></tt> occurrences)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i>&nbsp; <tt><b>{</b><i>n</i><b>,}</b></tt></td>
 * <td>(<tt><i>n</i></tt> or more occurrences)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>repeatexp </i>&nbsp;
 * <tt><b>{</b><i>n</i><b>,</b><i>m</i><b>}</b></tt></td>
 * <td>(<tt><i>n</i></tt> to <tt><i>m</i></tt> occurrences, including
 * both)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>complexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>complexp </i></td>
 * <td>::=</td>
 * <td><tt><b>~</b></tt> &nbsp; <i>complexp </i></td>
 * <td>(complement)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>charclassexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>charclassexp </i></td>
 * <td>::=</td>
 * <td><tt><b>[</b></tt> &nbsp; <i>charclasses </i>&nbsp; <tt><b>]</b></tt>
 * </td>
 * <td>(character class)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>[^</b></tt> &nbsp; <i>charclasses </i>&nbsp; <tt><b>]</b></tt>
 * </td>
 * <td>(negated character class)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>simpleexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>charclasses </i></td>
 * <td>::=</td>
 * <td><i>charclass </i>&nbsp; <i>charclasses </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>charclass </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>charclass </i></td>
 * <td>::=</td>
 * <td><i>charexp </i>&nbsp; <tt><b>-</b></tt> &nbsp; <i>charexp </i></td>
 * <td>(character range, including end-points)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><i>charexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * 
 * <tr>
 * <td><i>simpleexp </i></td>
 * <td>::=</td>
 * <td><i>charexp </i></td>
 * <td></td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>.</b></tt></td>
 * <td>(any single character)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>#</b></tt></td>
 * <td>(the empty language)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>@</b></tt></td>
 * <td>(any string)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>"</b></tt> &nbsp;&lt;Unicode string without
 * double-quotes&gt;&nbsp; <tt><b>"</b></tt></td>
 * <td>(a string)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>(</b></tt> &nbsp; <tt><b>)</b></tt></td>
 * <td>(the empty string)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>(</b></tt> &nbsp; <i>unionexp </i>&nbsp; <tt><b>)</b></tt>
 * </td>
 * <td>(precedence override)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>&lt;</b></tt> &nbsp;&lt;identifier&gt;&nbsp;
 * <tt><b>&gt;</b></tt></td>
 * <td>(named automaton)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>&lt;</b><i>n</i>-<i>m</i><b>&gt;</b></tt></td>
 * <td>(numerical interval)</td>
 * <td><small>[OPTIONAL] </small></td>
 * </tr>
 * 
 * <tr>
 * <td><i>charexp </i></td>
 * <td>::=</td>
 * <td>&lt;Unicode character&gt;</td>
 * <td>(a single non-reserved character)</td>
 * <td></td>
 * </tr>
 * <tr>
 * <td></td>
 * <td>|</td>
 * <td><tt><b>\</b></tt> &nbsp;&lt;Unicode character&gt;&nbsp;</td>
 * <td>(a single character)</td>
 * <td></td>
 * </tr>
 * </table>
 * <p>
 * The productions marked <small>[OPTIONAL] </small> are only allowed if
 * specified by the syntax flags passed to the <code>RegExp</code>
 * constructor. The reserved characters used in the (enabled) syntax must be
 * escaped with backslash (<tt><b>\</b></tt>) or double-quotes (
 * <tt><b>"..."</b></tt>). (In contrast to other regexp syntaxes, this is
 * required also in character classes.) Be aware that dash (<tt><b>-</b></tt>)
 * has a special meaning in <i>charclass </i> expressions. An identifier is a
 * string not containing right angle bracket (<tt><b>&gt;</b></tt>) or dash (
 * <tt><b>-</b></tt>). Numerical intervals are specified by non-negative
 * decimal integers and include both end points, and if <tt><i>n</i></tt> and
 * <tt><i>m</i></tt> have the same number of digits, then the conforming
 * strings must have that length (i.e. prefixed by 0's).
 * 
 * @author Anders M&oslash;ller &lt; <a
 *         href="mailto:amoeller@brics.dk">amoeller@brics.dk </a>&gt;
 */
public class RegExp {
	static final byte RE_UNION = 1;

	static final byte RE_CONCATENATION = 2;

	static final byte RE_INTERSECTION = 3;

	static final byte RE_OPTIONAL = 4;

	static final byte RE_REPEAT = 5;

	static final byte RE_REPEAT_MIN = 6;

	static final byte RE_REPEAT_MINMAX = 7;

	static final byte RE_COMPLEMENT = 8;

	static final byte RE_CHAR = 10;

	static final byte RE_CHAR_RANGE = 11;

	static final byte RE_ANYCHAR = 14;

	static final byte RE_EMPTY = 15;

	static final byte RE_STRING = 16;

	static final byte RE_ANYSTRING = 17;

	static final byte RE_AUTOMATON = 18;

	static final byte RE_INTERVAL = 19;

	/** Syntax flag, enables intersection (<tt>&amp;</tt>). */
	public static final int INTERSECTION = 0x0001;

	/** Syntax flag, enables complement (<tt>~</tt>). */
	public static final int COMPLEMENT = 0x0002;

	/** Syntax flag, enables empty language (<tt>#</tt>). */
	public static final int EMPTY = 0x0004;

	/** Syntax flag, enables anystring (<tt>@</tt>). */
	public static final int ANYSTRING = 0x0008;

	/**
	 * Syntax flag, enables named automata (<tt>&lt;</tt> identifier
	 * <tt>&gt;</tt>).
	 */
	public static final int AUTOMATON = 0x0010;

	/**
	 * Syntax flag, enables numerical intervals (
	 * <tt>&lt;<i>n</i>-<i>m</i>&gt;</tt>).
	 */
	public static final int INTERVAL = 0x0020;

	/** Syntax flag, enables all optional regexp syntax. */
	public static final int ALL = 0xffff;

	/** Syntax flag, enables no optional regexp syntax. */
	public static final int NONE = 0x0000;

	int kind;

	RegExp exp1, exp2;

	String s;

	char c;

	int min, max, digits;

	char from, to;

	StringBuffer b;

	int flags;

	int pos;

	RegExp() {
	}

	/**
	 * Constructs new <code>RegExp</code> from a string. Same as
	 * <code>RegExp(s, ALL)</code>.
	 * 
	 * @param s
	 *            regexp string
	 * @exception IllegalArgumentException
	 *                if an error occured while parsing the regular expression
	 */
	public RegExp(String s) throws IllegalArgumentException {
		this(s, ALL);
	}

	/**
	 * Constructs new <code>RegExp</code> from a string.
	 * 
	 * @param s
	 *            regexp string
	 * @param syntax_flags
	 *            boolean 'or' of optional syntax constructs to be enabled
	 * @exception IllegalArgumentException
	 *                if an error occured while parsing the regular expression
	 */
	public RegExp(String s, int syntax_flags) throws IllegalArgumentException {
		b = new StringBuffer(s);
		flags = syntax_flags;

		// create regexp
		RegExp e = parseUnionExp();

		if (pos < b.length()) // error
			throw new IllegalArgumentException("end-of-string expected at position " + pos);

		kind = e.kind; // copy to this
		exp1 = e.exp1;
		exp2 = e.exp2;
		this.s = e.s;
		c = e.c;
		min = e.min;
		max = e.max;
		digits = e.digits;
		from = e.from;
		to = e.to;

		b = null;
	}

	/**
	 * Constructs new <code>Automaton</code> from this <code>RegExp</code>.
	 * Same as <code>toAutomaton(null)</code> (empty automaton map).
	 */
	public Automaton toAutomaton() {
		return toAutomaton(null);
	}

	/**
	 * Constructs new <code>Automaton</code> from this <code>RegExp</code>.
	 * The constructed automaton is minimal and deterministic and has no
	 * transitions to dead states.
	 * 
	 * @param automata
	 *            a map from automaton identifiers to automata (of type
	 *            <code>Automaton</code>).
	 * @exception IllegalArgumentException
	 *                if this regular expression uses a named identifier that
	 *                does not occur in the automaton map
	 */
	public Automaton toAutomaton(Map automata) throws IllegalArgumentException {
		Automaton a = null;
		switch (kind) {
		case RE_UNION:
			a = exp1.toAutomaton(automata).union(exp2.toAutomaton(automata));
			a.minimize();
			break;
		case RE_CONCATENATION:
			a = exp1.toAutomaton(automata).concatenate(exp2.toAutomaton(automata));
			a.minimize();
			break;
		case RE_INTERSECTION:
			a = exp1.toAutomaton(automata).intersection(exp2.toAutomaton(automata));
			a.minimize();
			break;
		case RE_OPTIONAL:
			a = exp1.toAutomaton(automata).optional();
			a.minimize();
			break;
		case RE_REPEAT:
			a = exp1.toAutomaton(automata).repeat();
			a.minimize();
			break;
		case RE_REPEAT_MIN:
			a = exp1.toAutomaton(automata).repeat(min);
			a.minimize();
			break;
		case RE_REPEAT_MINMAX:
			a = exp1.toAutomaton(automata).repeat(min, max);
			a.minimize();
			break;
		case RE_COMPLEMENT:
			a = exp1.toAutomaton(automata).complement();
			a.minimize();
			break;
		case RE_CHAR:
			a = Automaton.makeChar(c);
			break;
		case RE_CHAR_RANGE:
			a = Automaton.makeCharRange(from, to);
			break;
		case RE_ANYCHAR:
			a = Automaton.makeAnyChar();
			break;
		case RE_EMPTY:
			a = Automaton.makeEmpty();
			break;
		case RE_STRING:
			a = Automaton.makeString(s);
			break;
		case RE_ANYSTRING:
			a = Automaton.makeAnyString();
			break;
		case RE_AUTOMATON:
			Automaton aa = (Automaton) automata.get(s);
			if (aa == null)
				throw new IllegalArgumentException(s + " not found");
			a = (Automaton) aa.clone();
			break;
		case RE_INTERVAL:
			a = Automaton.makeInterval(min, max, digits);
			break;
		}
		return a;
	}

	/** Constructs string from parsed regular expression. */
	public String toString() {
		return toStringBuffer(new StringBuffer()).toString();
	}

	StringBuffer toStringBuffer(StringBuffer b) {
		switch (kind) {
		case RE_UNION:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append("|");
			exp2.toStringBuffer(b);
			b.append(")");
			break;
		case RE_CONCATENATION:
			exp1.toStringBuffer(b);
			exp2.toStringBuffer(b);
			break;
		case RE_INTERSECTION:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append("&");
			exp2.toStringBuffer(b);
			b.append(")");
			break;
		case RE_OPTIONAL:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append(")?");
			break;
		case RE_REPEAT:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append(")*");
			break;
		case RE_REPEAT_MIN:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append("){").append(min).append(",}");
			break;
		case RE_REPEAT_MINMAX:
			b.append("(");
			exp1.toStringBuffer(b);
			b.append("){").append(min).append(",").append(max).append("}");
			break;
		case RE_COMPLEMENT:
			b.append("~(");
			exp1.toStringBuffer(b);
			b.append(")");
			break;
		case RE_CHAR:
			b.append("\\").append(c);
			break;
		case RE_CHAR_RANGE:
			b.append("[\\").append(from).append("-\\").append(to).append("]");
			break;
		case RE_ANYCHAR:
			b.append(".");
			break;
		case RE_EMPTY:
			b.append("#");
			break;
		case RE_STRING:
			b.append("\"").append(s).append("\"");
			break;
		case RE_ANYSTRING:
			b.append("@");
			break;
		case RE_AUTOMATON:
			b.append("<").append(s).append(">");
			break;
		case RE_INTERVAL:
			String s1 = (new Integer(min)).toString();
			String s2 = (new Integer(max)).toString();
			b.append("<");
			if (digits > 0)
				for (int i = s1.length(); i < digits; i++)
					b.append('0');
			b.append(s1).append("-");
			if (digits > 0)
				for (int i = s2.length(); i < digits; i++)
					b.append('0');
			b.append(s2).append(">");
			break;
		}
		return b;
	}

	/**
	 * Returns set of automaton identifiers that occur in this regular
	 * expression.
	 */
	public Set getIdentifiers() {
		HashSet set = new HashSet();
		getIdentifiers(set);
		return set;
	}

	void getIdentifiers(Set set) {
		switch (kind) {
		case RE_UNION:
		case RE_CONCATENATION:
		case RE_INTERSECTION:
			exp1.getIdentifiers(set);
			exp2.getIdentifiers(set);
			break;
		case RE_OPTIONAL:
		case RE_REPEAT:
		case RE_REPEAT_MIN:
		case RE_REPEAT_MINMAX:
		case RE_COMPLEMENT:
			exp1.getIdentifiers(set);
			break;
		case RE_AUTOMATON:
			set.add(s);
			break;
		}
	}

	static RegExp makeUnion(RegExp exp1, RegExp exp2) {
		RegExp r = new RegExp();
		r.kind = RE_UNION;
		r.exp1 = exp1;
		r.exp2 = exp2;
		return r;
	}

	static RegExp makeConcatenation(RegExp exp1, RegExp exp2) {
		RegExp r = new RegExp();
		r.kind = RE_CONCATENATION;
		r.exp1 = exp1;
		r.exp2 = exp2;
		return r;
	}

	static RegExp makeIntersection(RegExp exp1, RegExp exp2) {
		RegExp r = new RegExp();
		r.kind = RE_INTERSECTION;
		r.exp1 = exp1;
		r.exp2 = exp2;
		return r;
	}

	static RegExp makeOptional(RegExp exp) {
		RegExp r = new RegExp();
		r.kind = RE_OPTIONAL;
		r.exp1 = exp;
		return r;
	}

	static RegExp makeRepeat(RegExp exp) {
		RegExp r = new RegExp();
		r.kind = RE_REPEAT;
		r.exp1 = exp;
		return r;
	}

	static RegExp makeRepeat(RegExp exp, int min) {
		RegExp r = new RegExp();
		r.kind = RE_REPEAT_MIN;
		r.exp1 = exp;
		r.min = min;
		return r;
	}

	static RegExp makeRepeat(RegExp exp, int min, int max) {
		RegExp r = new RegExp();
		r.kind = RE_REPEAT_MINMAX;
		r.exp1 = exp;
		r.min = min;
		r.max = max;
		return r;
	}

	static RegExp makeComplement(RegExp exp) {
		RegExp r = new RegExp();
		r.kind = RE_COMPLEMENT;
		r.exp1 = exp;
		return r;
	}

	static RegExp makeChar(char c) {
		RegExp r = new RegExp();
		r.kind = RE_CHAR;
		r.c = c;
		return r;
	}

	static RegExp makeCharRange(char from, char to) {
		RegExp r = new RegExp();
		r.kind = RE_CHAR_RANGE;
		r.from = from;
		r.to = to;
		return r;
	}

	static RegExp makeAnyChar() {
		RegExp r = new RegExp();
		r.kind = RE_ANYCHAR;
		return r;
	}

	static RegExp makeEmpty() {
		RegExp r = new RegExp();
		r.kind = RE_EMPTY;
		return r;
	}

	static RegExp makeString(String s) {
		RegExp r = new RegExp();
		r.kind = RE_STRING;
		r.s = s;
		return r;
	}

	static RegExp makeAnyString() {
		RegExp r = new RegExp();
		r.kind = RE_ANYSTRING;
		return r;
	}

	static RegExp makeAutomaton(String s) {
		RegExp r = new RegExp();
		r.kind = RE_AUTOMATON;
		r.s = s;
		return r;
	}

	static RegExp makeInterval(int min, int max, int digits) {
		RegExp r = new RegExp();
		r.kind = RE_INTERVAL;
		r.min = min;
		r.max = max;
		r.digits = digits;
		return r;
	}

	/**
	 * 
	 * @param s
	 * @return whether s contains b[pos] and there are more 
	 * characters to process
	 */
	private boolean peek(String s) {
		return more() && s.indexOf(b.charAt(pos)) != -1;
	}

	private boolean match(char c) {
		if (pos >= b.length())
			return false;
		if (b.charAt(pos) == c) {
			pos++;
			return true;
		}
		return false;
	}

	/**
	 * 
	 * @return whether there are more chars left in the buffer
	 */
	private boolean more() {
		return pos < b.length();
	}

	private char next() throws IllegalArgumentException {
		if (!more())
			throw new IllegalArgumentException("unexpected end-of-string");
		return b.charAt(pos++);
	}

	private boolean check(int flag) {
		return (flags & flag) != 0;
	}

	RegExp parseUnionExp() throws IllegalArgumentException {
		RegExp e = parseInterExp();
		if (match('|'))
			e = makeUnion(e, parseUnionExp());
		return e;
	}

	RegExp parseInterExp() throws IllegalArgumentException {
		RegExp e = parseConcatExp();
		if (check(INTERSECTION) && match('&'))
			e = makeIntersection(e, parseInterExp());
		return e;
	}

	RegExp parseConcatExp() throws IllegalArgumentException {
		RegExp e = parseRepeatExp();
		if (more() && !peek(")&|")) {
			RegExp ee = parseConcatExp();
			if ((e.kind == RE_STRING || e.kind == RE_CHAR) && ee.kind == RE_CONCATENATION
					&& (ee.exp1.kind == RE_STRING || ee.exp1.kind == RE_CHAR)) {
				StringBuffer bb = new StringBuffer();
				if (e.kind == RE_STRING)
					bb.append(e.s);
				else
					bb.append(e.c);
				if (ee.exp1.kind == RE_STRING)
					bb.append(ee.exp1.s);
				else
					bb.append(ee.exp1.c);
				String ss = bb.toString();
				if (ss.indexOf('"') == -1) {
					e = makeString(ss);
					ee = ee.exp2;
				}
			}
			if ((e.kind == RE_STRING || e.kind == RE_CHAR)
					&& (ee.kind == RE_STRING || ee.kind == RE_CHAR)) {
				StringBuffer bb = new StringBuffer();
				if (e.kind == RE_STRING)
					bb.append(e.s);
				else
					bb.append(e.c);
				if (ee.kind == RE_STRING)
					bb.append(ee.s);
				else
					bb.append(ee.c);
				String ss = bb.toString();
				if (ss.indexOf('"') == -1)
					return makeString(ss);
			}
			e = makeConcatenation(e, ee);
		}
		return e;
	}

	RegExp parseRepeatExp() throws IllegalArgumentException {
		RegExp e = parseComplExp();
		while (peek("?*+{")) {
			if (match('?'))
				e = makeOptional(e);
			else if (match('*'))
				e = makeRepeat(e);
			else if (match('+'))
				e = makeRepeat(e, 1);
			else if (match('{')) {
				int start = pos;
				while (peek("0123456789"))
					next();
				if (start == pos)
					throw new IllegalArgumentException("integer expected at position " + pos);
				int n = Integer.parseInt(b.substring(start, pos));
				int m = -1;
				if (match(',')) {
					start = pos;
					while (peek("0123456789"))
						next();
					if (start != pos)
						m = Integer.parseInt(b.substring(start, pos));
				} else
					m = n;
				if (!match('}'))
					throw new IllegalArgumentException("expected '}' at position " + pos);
				if (m == -1)
					return makeRepeat(e, n);
				else
					return makeRepeat(e, n, m);
			}
		}
		return e;
	}

	RegExp parseComplExp() throws IllegalArgumentException {
		if (check(COMPLEMENT) && match('~'))
			return makeComplement(parseComplExp());
		else
			return parseCharClassExp();
	}

	RegExp parseCharClassExp() throws IllegalArgumentException {
		if (match('[')) {
			boolean negate = false;
			if (match('^'))
				negate = true;
			RegExp e = parseCharClasses();
			if (negate)
				e = makeIntersection(makeAnyChar(), makeComplement(e));
			if (!match(']'))
				throw new IllegalArgumentException("expected ']' at position " + pos);
			return e;
		} else
			return parseSimpleExp();
	}

	RegExp parseCharClasses() throws IllegalArgumentException {
		RegExp e = parseCharClass();
		while (more() && !peek("]"))
			e = makeUnion(e, parseCharClass());
		return e;
	}

	RegExp parseCharClass() throws IllegalArgumentException {
		char c = parseCharExp();
		if (match('-'))
			return makeCharRange(c, parseCharExp());
		else
			return makeChar(c);
	}

	RegExp parseSimpleExp() throws IllegalArgumentException {
		if (match('.'))
			return makeAnyChar();
		else if (check(EMPTY) && match('#'))
			return makeEmpty();
		else if (check(ANYSTRING) && match('@'))
			return makeAnyString();
		else if (match('"')) {
			int start = pos;
			while (more() && !peek("\""))
				next();
			if (!match('"'))
				throw new IllegalArgumentException("expected '\"' at position " + pos);
			return makeString(b.substring(start, pos - 1));
		} else if (match('(')) {
			if (match(')'))
				return makeString("");
			RegExp e = parseUnionExp();
			if (!match(')'))
				throw new IllegalArgumentException("expected ')' at position " + pos);
			return e;
		} else if ((check(AUTOMATON) || check(INTERVAL)) && match('<')) {
			int start = pos;
			while (more() && !peek(">"))
				next();
			if (!match('>'))
				throw new IllegalArgumentException("expected '>' at position " + pos);
			String s = b.substring(start, pos - 1);
			int i = s.indexOf('-');
			if (i == -1) {
				if (!check(AUTOMATON))
					throw new IllegalArgumentException("interval syntax error at position "
							+ (pos - 1));
				return makeAutomaton(s);
			} else {
				if (!check(INTERVAL))
					throw new IllegalArgumentException("illegal identifier at position "
							+ (pos - 1));
				try {
					if (i == 0 || i == s.length() - 1 || i != s.lastIndexOf('-'))
						throw new NumberFormatException();
					String smin = s.substring(0, i);
					String smax = s.substring(i + 1, s.length());
					int imin = Integer.parseInt(smin);
					int imax = Integer.parseInt(smax);
					int digits;
					if (smin.length() == smax.length())
						digits = smin.length();
					else
						digits = 0;
					if (imin > imax) {
						int t = imin;
						imin = imax;
						imax = t;
					}
					return makeInterval(imin, imax, digits);
				} catch (NumberFormatException e) {
					throw new IllegalArgumentException("interval syntax error at position "
							+ (pos - 1));
				}
			}
		} else
			return makeChar(parseCharExp());
	}

	char parseCharExp() throws IllegalArgumentException {
		match('\\');
		return next();
	}
}