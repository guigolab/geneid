/*
 * Created on May 6, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase;

import gphase.model.Exon;
import gphase.model.Gene;

import java.awt.Point;

/**
 * 
 * 
 * @author micha
 */
public class Toolbox {

	static public String[] constraintFormater(String[] sequences, Point[][] constraints) {
		
			// clone (not loosing original)
		StringBuffer[] seqs= new StringBuffer[sequences.length];
		for (int i = 0; i < sequences.length; i++) 
			seqs[i]= new StringBuffer(sequences[i].toLowerCase());
		
			// constraints to uppercase
		for (int i = 0; i < constraints.length; i++) {
			for (int j = 0; j < constraints[i].length; j++) {
				seqs[i].replace(constraints[i][j].x, (constraints[i][j].y+ 1), 
						seqs[i].substring(constraints[i][j].x, (constraints[i][j].y+ 1)).toUpperCase());
			}
		}
		
			// convert result
		String[] result= new String[seqs.length];
		for (int i = 0; i < seqs.length; i++) 
			result[i]= seqs[i].toString();
		
		return result;
	}
	
	static public Point[] extractConstraints(Gene g, int offset) {
				
		Exon[] ex= g.getExons();
		Point[] pa= new Point[ex.length];
		for (int i = 0; i < ex.length; i++) 
			pa[i]= new Point(ex[i].getStart()- offset, ex[i].getEnd()- offset);		
		
		return pa;
	}
}
