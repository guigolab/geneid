package qalign.tools;

import gphase.ext.DevPipeReaderThread;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.util.StringTokenizer;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AlnCompareWrapper {

	/**
	 * The executable location.
	 */
	public static final String EXEC_ABS_PATH= "C:\\Programme\\balibase\\cedric";
	/**
	 * The executable name.
	 */
	public static final String EXEC_NAME="aln_compare.exe";
	/**
	 * Option for parameter reference alignment.
	 */
	public static final String EXEC_OPT_REF_ALI= "-al1 ";
	/**
	 * Option for parameter test alignment.
	 */
	public static final String EXEC_OPT_TST_ALI= "-al2 ";
	/**
	 * Option to provide structure file.
	 */
	public static final String EXEC_OPT_STRUCT_FILE= "-st ";
	/**
	 * Option to adjust colums of output.
	 */
	public static final String EXEC_OPT_COL= "-io_cat ";
	/**
	 * Option to adjust lines of output.
	 */
	public static final String EXEC_OPT_OUTPUT_PAR= "-io_format ";
	/**
	 * Option to switch sp or column comparison.
	 */
	public static final String EXEC_OPT_COMP_MODE= "-compare_mode ";
	/**
	 * Default extension for reference alignments.
	 */
	public static final String DEF_REF_ALI_EXT= "ref_aln";
	/**
	 * Default extension for structure files.
	 */
	public static final String DEF_STRUCT_EXT= "cache_id ";

	/**
	 * test method
	 */
	public static void main(String[] args) {
		
		AlnCompareWrapper myWrapper= new AlnCompareWrapper(
			"C:\\Programme\\balibase\\cedric\\testcases\\1aab_ref1\\1aab_ref1.dca_aln");
		myWrapper.runAll();
	}



	/**
	 * Filebase directory.
	 */
	protected String filePath= null;
	/**
	 * Filebase name.
	 */
	protected String fileName= null;
	/**
	 * Filebase extension.
	 */
	protected String fileExt= null;
	/**
	 * Default column comparison extension.
	 */
	protected String colExt= "col";
	/**
	 * Default sp comparison extension.
	 */
	protected String spExt= "sp";
	/**
	 * Location of the structure file in absolute terms.
	 */
	protected String structFileAbs= null;

	/**
	 * Name of the structure comparison column.
	 */
	protected String structColName= "struc";
	/**
	 * Name of the total comparsion column.
	 */
	protected String totalColName= "Tot";
	/**
	 * Column format configuration string.
	 */
	protected String columnFormat= null;
	/**
	 * Line output configuration string.
	 */
	protected String lineFormat= null;


	/**
	 * Column score w/o structure file.
	 */
	protected float flatSPScore= 0f;
	/**
	 * SP score w structure file.
	 */
	protected float structSPScore= 0f;
	/**
	 * Column score w/o structure file.
	 */
	protected float flatColScore= 0f;
	/**
	 * SP score w structure file.
	 */
	protected float structColScore= 0f;
	
	

	/**
	 * Constructs <code>AlnCompareWrapper</code> and sets default parameters.
	 */
	public AlnCompareWrapper() {
		setDefaults();
	}
		
	/**
	 * Constructs <code>AlnCompareWrapper</code>, inits filebase and
	 * sets default parameters.
	 */
	public AlnCompareWrapper(String newTestFileAbs) {
		setFileBase(newTestFileAbs);
		setDefaults();	// afterwards to include filebase for struct file
	}
	
	/**
	 * Sets default comparison parameters:
	 * 	- activates structure comparison
	 *  - only average totals are printed
	 */
	public void setDefaults() {
		setStructure();
		setOutputLines(true, false, false);
	}
	
	/**
	 * Sets all options for output lines.
	 */
	public void setOutputLines(boolean outAvgTotal, 
								boolean outAvgSeq, 
								boolean outPairw) {
		
		lineFormat= "h";	// header
		if (outAvgTotal)
			lineFormat+= "t";
		if (outAvgSeq)
			lineFormat+= "s";
		if (outPairw)
			lineFormat+= "p";
	}
	/**
	 * Sets additional structure comparison parameters:
	 * 	- structure file included (default name)
	 *  - h and e matches are pooled in structure column
	 * 	- all residue matches are given in total column
	 */
	public void setStructure() {
		if (structFileAbs== null)
			setStructFileAbs(filePath+ File.separator+
							fileName+ "."+ DEF_STRUCT_EXT);
		this.columnFormat= threeDeeAli();
	}
	
	/**
	 * Sets the location of the structure file in absolute terms.
	 */
	public void setStructFileAbs(String absolutePath) {
		
		this.structFileAbs= absolutePath;
	}

	/**
	 * Pre-sets the default pairs to be taken into account for every column.
	 * To be specific, residue pairs of the h and e regions from structure file 
	 * are added up to form the structure comparison column, all residue pairs
	 * are regarded to sum up the total comparison column.
	 */	
	public String threeDeeAli() {

			// does not work properly in windows
			// just the struct is returned!		
/*		return "'"
			+ createColumn(align('h','h')+"+"+ align('e','e'), this.structColName)
			+ createColumn(align('*','*'), this.totalColName)
			+ "'";
*/			
		return "3d_ali";
	}
	
	/**
	 * Creates a valid alignment String for two types of structures.
	 */
	protected String align(char c1, char c2) {
		
		return "["+ c1+ "]["+ c2+ "]";
	}
	
	protected String createColumn(String alignStr, String colName) {
		
		return alignStr+ "="+ colName;
	}
	
	protected void setFileBase(String newAbsFileName) {
		
		StringTokenizer st= new StringTokenizer(newAbsFileName, File.separator);

//		System.out.println(newAbsFileName);
			// retrieve path
		this.filePath= "";// File.separator;
		int end= st.countTokens();
		String tmpToken= "";
		for (int i= 0; i< (end- 1); ++i) {
			tmpToken= st.nextToken();
			this.filePath+= tmpToken+ File.separator;
		}
			// kill last File.separator
		this.filePath= this.filePath.substring(0, this.filePath.length()- 1);
//		System.out.println("path: "+filePath);
		
			// retrieve name
		tmpToken= st.nextToken();
		st= new StringTokenizer(tmpToken, ".");
		this.fileName= st.nextToken();
//		System.out.println("name: "+fileName);
		
			// retrieve ext
		this.fileExt= tmpToken.substring(this.fileName.length()+ 1);
//		System.out.println("ext: "+fileExt);
	}	


	
	
	
	/**
	 * Runs a comparison for the given files with the parameters set.
	 */
	public int runComparison(boolean col) {
		
			// must be array command type on SUN
			// ("sh -c ..." does NOT work !)
		String[] command= {"cmd", "/c", 
			EXEC_ABS_PATH+ File.separator+ EXEC_NAME+ " "
			+ EXEC_OPT_REF_ALI
					+ filePath+ File.separator+ fileName+ "."+ DEF_REF_ALI_EXT+ " "
			+ EXEC_OPT_TST_ALI
					+ filePath+ File.separator+ fileName+ "."+ fileExt};
		
		if (structFileAbs!= null)
			command[2]+= " "+ EXEC_OPT_STRUCT_FILE+ structFileAbs;
		
		if (columnFormat!= null)
			command[2]+= " "+ EXEC_OPT_COL+ columnFormat;
		
		if (lineFormat!= null)
			command[2]+= " "+ EXEC_OPT_OUTPUT_PAR+ lineFormat;
		
		if (col)
			command[2]+= " "+ EXEC_OPT_COMP_MODE+ "column";
		
//		System.out.println(command[2]);
		OutputStream out= null;
		try {
			String outFile= filePath+ File.separator+ fileName+ ".";
			if (col)
				outFile+= colExt;
			else
				outFile+= spExt;
			out= new FileOutputStream(outFile);
		} catch (IOException e) {
			e.printStackTrace(); //nothing
		}
		Process proc= null;
		DevPipeReaderThread pipe= null;
		try {
			pipe= new DevPipeReaderThread(null, out);
			proc= Runtime.getRuntime().exec(command);
			pipe.setIn(proc.getInputStream());
			pipe.start();
		} catch (IOException e) {
			System.err.println("nope"); // nothing
		}
		
		int res= -1;
		try {
			res= proc.waitFor();
		} catch (InterruptedException e) {
			System.err.println("nope"); // nothing
		}
		try {
			while(!pipe.hasFinished()) {
				try {
					Thread.currentThread().sleep(100l);
				} catch (InterruptedException e) {
					; // nothing
				}
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return res;
	}
	
	/**
	 * Runs both, SP- and columnwise comparison of the alignments.
	 * @returns float[] SP-score(struct), SP-score(total), 
	 * Col-score(struct), Col-score(total)
	 */
	public float[] runAll() {
		
		System.out.print(filePath+ File.separator+ fileName+ "..");
		int tmpRes= -1;

		float[] result= new float[4];
		
		tmpRes= runComparison(false);
		try {
			float[] tmp= read(filePath+ File.separator+ fileName+ "."+ spExt);
			result[0]= tmp[0];
			result[1]= tmp[1];
		} catch (Exception e) {
			System.err.println("Error during parsing SP comparison file!");
			e.printStackTrace();
		}
			// Cedric returns sth like: 2009246880
/*		if (tmpRes!= 0)
			System.err.println("WARNING: "
				+ filePath+ File.separator+ fileName
				+ " failed in sp comparison. ["+ tmpRes+"]"
			);
*/
		tmpRes= runComparison(true);
		try {
			float[] tmp= read(filePath+ File.separator+ fileName+ "."+ colExt);
			result[2]= tmp[0];
			result[3]= tmp[1];
		} catch (Exception e) {
			System.err.println("Error during parsing column comparison file!");
			e.printStackTrace();
		}
/*		if (tmpRes!= 0)
			System.err.println("WARNING: "
				+ filePath+ File.separator+ fileName
				+ " failed in col comparison. ["+ tmpRes+"]"
			);
*/		
		System.out.println("done.");
		return result;
	}	
	
	
	/**
	 * Reads and parses comparison files. Just takes the first line as total scores.
	 */
	public float[] read(String absFName) throws Exception {
		
		
		BufferedReader reader= null;
		reader= new BufferedReader(
				new FileReader(absFName));
		
		String tmp= "";
		while (reader.ready()&& !tmp.trim().startsWith("seq1"))
			tmp= reader.readLine();
		tmp= reader.readLine();


		StringTokenizer st= new StringTokenizer(tmp);
		float[] res= new float[2];
				
		st.nextToken(); 	// ref aln name
		st.nextToken();		// # of seqs in tst aln name

		st.nextToken();		// avg similarity

		tmp= st.nextToken();
		res[0]= Float.parseFloat(tmp);	// structure score
		st.nextToken();		// '['
		st.nextToken();		// [total struct]

		tmp= st.nextToken();
		res[1]= Float.parseFloat(tmp);
		st.nextToken();		// '['
		st.nextToken();		// [total struct]
		
		return res;
	}
	
}
