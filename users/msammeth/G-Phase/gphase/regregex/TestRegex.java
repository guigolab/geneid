/*
 * Created on May 5, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.regregex;

/**
 * 
 * 
 * @author msammeth
 */
public class TestRegex {

	public static void main(String[] args) {
		RegExp myRE= new RegExp("(a|c)*b");
		Automaton auto= myRE.toAutomaton();
		System.currentTimeMillis();
	}
}
