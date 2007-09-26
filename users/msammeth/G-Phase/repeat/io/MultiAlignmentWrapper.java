/*
 * Created on Jul 16, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import qalign.algo.dialign.MultiFrag;
import repeat.data.Alignment;


/**
 * 
 * 
 * @author micha
 */
public class MultiAlignmentWrapper {

	public final static int TYPE_FLAT_FILE= 0;
	
	public static void main(String[] args) {

		Alignment[] ali= readAlignments(
			"D:\\Eigene Dateien\\repeats\\data\\BAliBASE\\ref6\\test2c\\sh3_2\\repeat_detection\\trust\\2c_sh3_2_detect_trust.txt",
			TYPE_FLAT_FILE
		);
		
			// tst out
		MultiFrag[] frags= ali[1].getMultiFrags();
		for (int i= 0; i < frags.length; i++) 
			System.out.println(frags[i]);	
	}
	
	public static Alignment[] readAlignments(String absFileName, int type) {
	
		try {
			switch (type) {
				case 0: return readAlignmentsFromFlatfile(absFileName); 		
			}
		} catch (Exception e) {		// FileNotFound, IO
			e.printStackTrace();
		}
		
		return null;
	}
	
	public static Alignment[] readAlignmentsFromFlatfile(String absFileName) throws Exception { 
		
		Pattern p= Pattern.compile("(\\w+)\\s+(\\d+)\\s+(\\D+)\\s+(\\d+)");
		BufferedReader buffy= new BufferedReader(new FileReader(absFileName));
		String read= "";
		Alignment ali= null;
		Vector result= new Vector();
		
		while (buffy.ready()) {
			
			read= buffy.readLine().trim();		
			if (read.length()< 1) {				// skip empty lines
				if (ali!= null) {				// save alignment
					result.add(ali);
					ali= null;
				}
				continue;
			}
			
				// read signs
			if (ali== null)
				ali= new Alignment();
			
			Matcher m= p.matcher(read);
			if (m.matches()) {
				ali.addSeqID(m.group(1).trim());
				ali.addSequence(m.group(3).trim());
				ali.addStartPos(Integer.parseInt(m.group(2).trim())- 1);
				ali.addEndPos(Integer.parseInt(m.group(4).trim())- 1);
			}
		}	// until no more can be read
		if (ali!= null)
			result.add(ali);	// add last
		
			// convert
		Alignment[] res= new Alignment[result.size()];
		for (int i= 0; i < res.length; i++) 
			res[i]= (Alignment) result.elementAt(i);	
		
		return res;
	}
	
}
