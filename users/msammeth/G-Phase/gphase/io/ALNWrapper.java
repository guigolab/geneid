/*
 * Created on May 6, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 
 * 
 * @author micha
 */
public class ALNWrapper {
	
	File file= null;
	String[] names= null;
	String[] sequences= null;
	
	public static void main(String[] args) {
		
		ALNWrapper myaln= new ALNWrapper("D:\\projects\\g-phase\\data\\seq\\evloutiveAS3\\last_rat.aln");
		myaln.read();
		for (int i = 0; i < myaln.getNames().length; i++) 	
			System.out.println(myaln.getNames()[i]+ "\t"+ myaln.getSequences()[i]);
	}
	
	public ALNWrapper(String newFileName) {
		file= new File(newFileName);
	}
	
	public void read() {
		
		if (file== null)	// error 
			return;
		
		Pattern patty= Pattern.compile("^(\\w+)\\s+([\\-a-zA-Z]+)\\s*(\\d*)\\s*$");
		Matcher matty;
		
		Vector growingSeqs= null;
		Vector nameVec= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(file));
			matty= patty.matcher(buffy.readLine());
			String s= "";
			while (buffy.ready()) {
			
				while (!matty.matches())		// find next tupel
					matty= patty.matcher(s= buffy.readLine());
			
				if (growingSeqs== null) {		// 1st round
					growingSeqs= new Vector();
					while (matty.matches()) {
						nameVec.add(matty.group(1));
						growingSeqs.add(new StringBuffer(matty.group(2)));
						matty= patty.matcher(s= buffy.readLine());
					}
					continue;
				}
				
				for (int i = 0; i < growingSeqs.size(); i++) {	// other rounds
					((StringBuffer) growingSeqs.elementAt(i)).append(matty.group(2));
					matty= patty.matcher(s= buffy.readLine());
					matty.matches();
				}
			}			
			
		} catch (Exception e) {
			System.out.println();; // : -)
		}
		
			// convert result
		names= new String[nameVec.size()];
		sequences= new String[nameVec.size()];
		for (int i = 0; i < nameVec.size(); i++) {
			names[i]= (String) nameVec.elementAt(i);
			sequences[i]= ((StringBuffer) growingSeqs.elementAt(i)).toString();
		}		
	}
	/**
	 * @return Returns the names.
	 */
	public String[] getNames() {
		return names;
	}
	/**
	 * @param names The names to set.
	 */
	public void setNames(String[] names) {
		this.names = names;
	}
	/**
	 * @return Returns the sequences.
	 */
	public String[] getSequences() {
		return sequences;
	}
	/**
	 * @param sequences The sequences to set.
	 */
	public void setSequences(String[] sequences) {
		this.sequences = sequences;
	}
}
