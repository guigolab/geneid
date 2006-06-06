/*
 * dk.brics.automaton
 * Copyright (C) 2001-2004 Anders Moeller
 * All rights reserved.
 */

package gphase.regex;


import gphase.tools.Arrays;

import java.util.*;

/**
 * excluded: intersection (&), negation (^)
 * @author micha
 */
public class RegExp {
	
	String[] tokens;

	
	public static void main(String[] args) {
		RegExp regex= new RegExp("2(3,4)*6");
		Automaton auto= regex.toAutomaton();
		System.currentTimeMillis();
	}

	public RegExp(String s) {
	
			// init tokenizer
		StringBuffer sb= new StringBuffer(s);
		for (int i = 0; i < sb.length(); i++) {
			if (sb.charAt(i)== ',')
				continue;
			if (Character.isDigit(sb.charAt(i))) {
				if (i> 0&& (!Character.isDigit(sb.charAt(i-1)))&& sb.charAt(i-1)!= ',') 
					sb.insert((i++), ',');
				if (i+1< sb.length()&& (!Character.isDigit(sb.charAt(i+1)))&& sb.charAt(i+1)!= ',') 
					sb.insert((++i), ',');
			}
		}
		StringTokenizer toki= new StringTokenizer(sb.toString(), ",");
		tokens= new String[toki.countTokens()];
		for (int i = 0; i < tokens.length; i++) 
			tokens[i]= toki.nextToken();
		
	}

	/**
	 * Constructs new <code>Automaton</code> from this <code>RegExp</code>.
	 * Same as <code>toAutomaton(null)</code> (empty automaton map).
	 */
	public Automaton toAutomaton() {
		
		Automaton auto= new Automaton();
		Stack repeatPos= new Stack();
		repeatPos.push(auto.getInitialState());	// for skips from the beginning
		
		int pos= 0;
		int label= -1;
		State saveSourceForRepeats= null;
		boolean donor= true;	// assume most patterns to start with donors
		while (pos< tokens.length) {
			try {
				label= Integer.parseInt(tokens[pos]);
				saveSourceForRepeats= auto.getEndState();
				auto.createTransition(label, donor);
				donor= !donor;
				pos++;
			} catch (NumberFormatException e) {
				if (tokens[pos].equals("^"))
					auto.getInitialState().setStart(true);
				if (tokens[pos].equals("$"))
					auto.getEndState().setEnd(true);
				if (tokens[pos].equalsIgnoreCase("A")) 
					donor= false;
				if (tokens[pos].equalsIgnoreCase("D")) 
					donor= true;
				if (tokens[pos].equals("("))
					repeatPos.push(auto.getEndState());
				if (tokens[pos].startsWith(")")) {					
					if (tokens[pos].charAt(1)== '+'|| tokens[pos].charAt(1)== '*') {		// loop
						State loopState= (State) repeatPos.pop();
						boolean loopDonor= saveSourceForRepeats.transitions[0].isDonor();
						saveSourceForRepeats.addTransition(label, loopDonor, loopState);	
						if (tokens[pos].charAt(1)== '*')	// for reuse
							repeatPos.push(loopState);
					}
					if (tokens[pos].charAt(1)== '?'|| tokens[pos].charAt(1)== '*') {		// skip
						if (pos+ 1< tokens.length) 
							try {
								label= Integer.parseInt(tokens[pos+ 1]);
							} catch (NumberFormatException ex) {
								System.err.println("Format error, number expected: "+ tokens[pos+ 2]);
							}
						else
							label= -1;	// epsilon
									
							
						saveSourceForRepeats= auto.getEndState();
						auto.createTransition(label, donor);			// transition to next wo skip
						
						State skipState= ((State) repeatPos.pop());		// skip transition
						boolean skipDonor= skipState.getTransitions()[0].isDonor();
						skipState.addTransition(label, skipDonor, auto.getEndState());	
						pos++;
						donor= !donor;
					}
				}
				pos++;
			}
			
		}
		
		auto.sortTransitions();
		return auto;
	}
}