package qalign.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Vector;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class BaliScoreParser {

	public static String defaultPathBase= "C:\\Programme\\balibase\\results\\bali_score_tst";
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
	
	/**
	 * Constructor for BaliScoreParser.
	 */
	public BaliScoreParser(String pathBase) {
		this.defaultPathBase= pathBase;
	}

	public static void main(String[] args) {
		BaliScoreParser parser= new BaliScoreParser(defaultPathBase);
		parser.run(0, BALIBASE_SUBDIRS.length, "*");
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
		
		BaliScoreReader scoreReader= null;
		TimeMSFWrapper timeReader= null;
		String tmpFName= null;
		boolean contFlag= false;
		
		BufferedWriter writer= null;
		BufferedWriter timeWriter= null;
		try {
			writer= new BufferedWriter(new FileWriter(
				defaultPathBase+ File.separator+ "table.txt"));
			writer.write("Name , SP-glob , TC-glob , SP-core, TC-core\n");
			timeWriter= new BufferedWriter(new FileWriter(
				defaultPathBase+ File.separator+ "timetable.txt"));
			writer.write("Name , Hyperspace , Consistency\n");
		} catch (Exception e) {
			e.printStackTrace();
		}
		StringBuffer tmpOut;
		for (int i= startDir; i< stopDir; ++i) {
			
			String currSubDir= defaultPathBase+ File.separator+ BALIBASE_SUBDIRS[i];
			String[] fNames= getFileNames(currSubDir);
			if (fNames== null)
				break; // end

			if (i== startDir)
				contFlag= true;
			for (int j= 0; j< fNames.length; ++j) {
				
					// find name to start with
				if (contFlag) {
					if (contFlag
						&& ((fNames[j].equalsIgnoreCase(startFName))
							|| (startFName.equals("*"))))
						contFlag= false;
					if (contFlag)
						continue;
				}
					
				tmpFName= currSubDir+ File.separator+ fNames[j];
				System.out.println(fNames[j]);
				scoreReader= new BaliScoreReader(tmpFName);
				tmpFName= tmpFName.substring(0, tmpFName.indexOf("."));
				timeReader= new TimeMSFWrapper(tmpFName+ ".fa.msf");
				
					// read & output
				try {
					scoreReader.readBoth();
					writer.write(BALIBASE_SUBDIRS[i]+ File.separator+ fNames[j]+ " , ");
					writer.write(scoreReader.spScoreGlob+ " , "+ scoreReader.tcScoreGlob+ " , ");
					writer.write(scoreReader.spScoreCore+ " , "+ scoreReader.tcScoreCore+ "\n");

					timeReader.read();
					timeWriter.write(BALIBASE_SUBDIRS[i]+ File.separator+ fNames[j]+ ".msf"+ " , ");
					timeWriter.write(Float.toString(timeReader.getHyperTime()));
					if (timeReader.isConsistentTimeFound())
						timeWriter.write(" , "+ timeReader.getConsistentTime()+ "\n");
					else
						timeWriter.write("\n");
					
				} catch (Exception e) {
					e.printStackTrace();
				}
			}	// for all files in a subdir
		}	// for all subdirs of balibase
		
		try {
			writer.close();
			timeWriter.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
}
