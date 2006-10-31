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
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * 
 * 
 * @author micha
 */
public class Toolbox {

	public static String getAbsFileName(String fName) {
		return new File(fName).getAbsolutePath();
	}
	
	/**
	 * 
	 * @param fName
	 * @return absolute file name
	 */
	public static String checkFileExists(String fName) {
		File f= new File(fName);
		if (f.exists()) {
			System.out.println("File "+fName+" exists. Specify new name or press <CR> to overwrite.");
			BufferedReader r= new BufferedReader(new InputStreamReader(System.in));
			String s= "";
			try {
				s= r.readLine();
			} catch (IOException e) {
				e.printStackTrace();
			}
			s= s.trim();
			if (s.length()< 1) {
				System.out.println("Overwritten.");
				return getAbsFileName(fName);
			}
			String p= "";
			int pos= fName.lastIndexOf(File.separator);
			if (pos>= 0)
				p= fName.substring(0, pos);
			if (p.length()> 0)
				s= p+ File.separator+ s;
			System.out.println("Redirected to file "+ s);
			return getAbsFileName(s);
		}
		System.out.println("Output in file "+ fName);
		return getAbsFileName(fName);
	}
	
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
