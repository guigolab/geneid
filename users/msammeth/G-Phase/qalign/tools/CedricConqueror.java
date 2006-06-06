package qalign.tools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;

import qalign.algo.dca.DCATester;
import qalign.algo.dca.extension.DCAClosureTester;
import qalign.algo.tcoffee.TCoffeeWrapper;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class CedricConqueror {
	
	public static final String DEFAULT_PATH_BASE="C:\\Programme\\balibase\\cedric\\testcases";
	public static final String DEFAULT_FASTA_EXT="pep";
	public static final String DEFAULT_TCOFFEE_EXT="cw175nw_aln";
	
	
	public static void main(String[] args) {
		
		CedricConqueror myConqueror= new CedricConqueror(DEFAULT_PATH_BASE, 2, "1cpt_ref2", 2, DEFAULT_FASTA_EXT+ ".msf");
	}

	public CedricConqueror(String pathBase, int startRef, String startAli, int endRef, String tstExt) {
		
		run(pathBase, startRef, startAli, endRef, tstExt);		
	}

	public void run(String pathBase, int startRef, String startAli, int endRef, String tstExt) {

			// retrieve sub-dirs		
		File dir= new File(pathBase);
		String[] subDir= dir.list();
		Arrays.sort(subDir);

		DCAClosureTester dcaTest= null;
//		DCATester dcaTest= null;
		AlnCompareWrapper alnCompare= null;
		float[] tmpResult= null, newResult= null;
		long start, end= 0l;
		
			// init writer for logging
		BufferedWriter writer= null;
		try {
			writer= new BufferedWriter(new FileWriter(
				DEFAULT_PATH_BASE+ "results.txt"));	
			writer.write("\n\n\nName, SP-struct, SP-total, Col-struct, Col-total, time\n\n");
		} catch (Exception e) {
			; // nothing, again..
		}
			// to sort according ref#
		Hashtable results= new Hashtable(5);
		results.put(new Integer(1), new Vector());
		results.put(new Integer(2), new Vector());
		results.put(new Integer(3), new Vector());
		results.put(new Integer(4), new Vector());
		results.put(new Integer(5), new Vector());

			// for starting point
		boolean contFlag= true;
		if (startAli.equals("*"))
			contFlag= false;
		for (int i= 0; i< subDir.length; ++i) {

				// skip if not desired
			int refNr= Integer.parseInt(
				subDir[i].substring(subDir[i].length()-1, subDir[i].length()));				
			if (contFlag&& startAli.equalsIgnoreCase(subDir[i]))
				contFlag= false;
			if ((contFlag&& refNr== startRef)|| refNr< startRef|| refNr> endRef) {
//			if (contFlag) {
				System.out.println("skipping: "+subDir[i]);
				continue;
			}
			
			System.out.println("processing: "+subDir[i]);
				// perform alignment
//			System.out.println(subDir[i]);
			dcaTest= new DCAClosureTester(
//			dcaTest= new DCATester(
				pathBase+ File.separator+
				subDir[i]+ File.separator+
				subDir[i]+ "."+ DEFAULT_FASTA_EXT);
/*			TCoffeeWrapper tcof= new TCoffeeWrapper();	
			tcof.setOutputMSF(true);
*/			
			start= System.currentTimeMillis();
			try {
				dcaTest.run();
/*				tcof.runTCoffee(
					pathBase+ File.separator+
					subDir[i]+ File.separator+
					subDir[i]+ "."+ DEFAULT_FASTA_EXT);
*/					
			} catch (OutOfMemoryError e) {
				e.printStackTrace();
				System.err.println("\n--> continuing...");
				continue;
			}
			end= System.currentTimeMillis();
			float time= (float) ((end- start)/ 1000f);
				
				// start comparison
			String tmpName= 
				pathBase+ File.separator+ 
				subDir[i]+ File.separator+ 
				subDir[i]+ "."+ tstExt;

			alnCompare= new AlnCompareWrapper(tmpName);
			
				// get results and add time
			tmpResult= alnCompare.runAll();
			newResult= new float[tmpResult.length+ 1];
			for (int j= 0; j< tmpResult.length; ++j)
				newResult[j]= tmpResult[j];
			newResult[newResult.length- 1]= time;
			
			((Vector) results.get(new Integer(refNr))).add(subDir[i]);
			((Vector) results.get(new Integer(refNr))).add(newResult);	// no need to put again

		}
		
		Vector tmpVec= null;
				// log results
		for (int i= 1; i<= 5; ++i) {
			
			tmpVec= (Vector) results.get(new Integer(i));
			for (int j= 0; j< tmpVec.size(); ++j) {
				
				try {
					writer.write(((String) tmpVec.elementAt(j))+"."+ tstExt+ ", ");
					tmpResult= (float[]) tmpVec.elementAt(j+ 1);
					for (int k= 0; k< (tmpResult.length- 1); ++k)
						writer.write(tmpResult[k]+ ", ");
					writer.write(tmpResult[tmpResult.length- 1]+ "\n");
				} catch (Exception e) {
					; // nothing, again..
				}
			}
			
			try {
				writer.write("\n\n");
			} catch (Exception e) {
				; // nothing, again..
			}
		}

			// close stream		
		try {
			writer.close();
		} catch (Exception e) {
			; // nothing, again..
		}
	}


}
