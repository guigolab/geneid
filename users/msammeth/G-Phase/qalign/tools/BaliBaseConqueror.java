package qalign.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;

import qalign.algo.dca.DCATester;
import qalign.algo.dca.extension.DCAClosureTester;
import qalign.algo.dialign.DialignWrapper;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class BaliBaseConqueror {
	
	public static String defaultPathBase= "C:\\Programme\\balibase";
	public static final String INDEX_FNAME= "00index.lst";
	public static final String[] BALIBASE_SUBDIRS= {
		"ref1"+File.separator+"test1",
		"ref1"+File.separator+"test2",
		"ref1"+File.separator+"test3",
		"ref2"+File.separator+"test",
		"ref3"+File.separator+"test",
//		"ref3"+File.separator+"test1",
		"ref4"+File.separator+"test",
		"ref5"+File.separator+"test"
	};
		

	public static void main(String[] args) {
		
		BaliBaseConqueror conqueror= new BaliBaseConqueror(defaultPathBase);
		conqueror.run(0, BALIBASE_SUBDIRS.length, "*");
	}

	/**
	 * Constructor for BaliBaseConqueror.
	 */
	public BaliBaseConqueror(String pathBase) {
		this.defaultPathBase= pathBase;
	}

	public String[] getFileNames(String directory) {
		
			// get index file
		BufferedReader idxReader= null;
		try {
			idxReader= new BufferedReader(new FileReader(
				directory+ File.separator+ INDEX_FNAME));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return null;
		}
		
			// read file names
		Vector result= new Vector();
		try {
			while (idxReader.ready())
				result.add(idxReader.readLine().trim());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
			// convert to String[]
		String[] res= new String[result.size()];
		for (int i= 0; i< result.size(); ++i)
			res[i]= (String) result.elementAt(i);
		
		return res;
	}
	
	public void run(int startDir, int stopDir, String startFName) {
		
		DCAClosureTester dcaTest= null;
		BaliScoreWrapper baliScore= null;
		String tmpFName= null;
		boolean contFlag= false;
		for (int i= startDir; i< stopDir; ++i) {
			
			String currSubDir= defaultPathBase+ File.separator+ BALIBASE_SUBDIRS[i];
			String[] fNames= getFileNames(currSubDir);
			if (i== startDir)
				contFlag= true;
			for (int j= 0; j< fNames.length; ++j) {
				
					// find name to start with
				if (contFlag) {
					if ((fNames[j].equalsIgnoreCase(startFName))
							|| (startFName.equals("*")))
						contFlag= false;
					if (contFlag)
						continue;
				}
					
				tmpFName= currSubDir+ File.separator+ fNames[j];
				dcaTest= new DCAClosureTester(tmpFName);
				try {
					dcaTest.run();
				} catch (OutOfMemoryError e) {
					e.printStackTrace();
					System.err.println("\n--> continuing...");
				}
				
				baliScore= new BaliScoreWrapper(tmpFName+ ".msf");	// written by DCA
				baliScore.runAll();
			}	// for all files in a subdir
		}	// for all subdirs of balibase
	}
	
	public void runDialign(int startDir, int stopDir, String startFName) {
		
		DialignWrapper dia= null;
		StringTokenizer st= null;
		BaliScoreWrapper baliScore= null;
		String tmpFName= null;
		boolean contFlag= false;
		for (int i= startDir; i< stopDir; ++i) {
			
			String currSubDir= defaultPathBase+ File.separator+ BALIBASE_SUBDIRS[i];
			String[] fNames= getFileNames(currSubDir);
			if (i== startDir)
				contFlag= true;
			for (int j= 0; j< fNames.length; ++j) {
				
					// find name to start with
				if (contFlag) {
					if ((fNames[j].equalsIgnoreCase(startFName))
							|| (startFName.equals("*")))
						contFlag= false;
					if (contFlag)
						continue;
				}
					
				tmpFName= currSubDir+ File.separator+ fNames[j];
				
					// start dialign
				System.out.print("running DIALIGN..");
				System.out.flush();
				dia= new DialignWrapper();
				dia.setOutputMSF(true);
				if (dia.runDialign(tmpFName)!= 0)
					System.err.println("Error in DIALIGN run: "+ tmpFName);
				st= new StringTokenizer(tmpFName, ".");
				int len= st.countTokens();
				tmpFName= "";
				for (int x= 0; x< (len- 1); ++x)
					tmpFName+= st.nextToken();
				baliScore= new BaliScoreWrapper(tmpFName+ ".ms");	// written by DCA
				baliScore.runAll();
			}	// for all files in a subdir
		}	// for all subdirs of balibase
	}
	
}
