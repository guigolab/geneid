package qalign.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class BaliScoreReader {
	
	public static final String TEST_ALI_TEXT= "Comparing test alignment in";
	public static final String REF_ALI_TEXT= "with reference alignment in";
	public static final String CORE_ANN_TEXT= "Using core blocks defined in";
	public static final String SP_SCORE_TEXT= "SP score=";
	public static final String TC_SCORE_TEXT= "TC score=";
	
	public static final String GLOB_EXT= "glob";
	public static final String CORE_EXT= "core";

	protected String fileName= null;
	protected String absPath= null;
	protected boolean readModeGlobal= true;

	protected String testAli= null;
	protected String refAli= null;
	protected String annFile= null;
	protected int seqNb= 0;
	protected StringBuffer[] sequences= null;
	protected String[] seqNames= null;
	protected int[] scoresGlob= null;
	protected int[] scoresCore= null;
	protected double spScoreGlob= 0d;
	protected double tcScoreGlob= 0d;
	protected double spScoreCore= 0d;
	protected double tcScoreCore= 0d;

	/**
	 * Constructor for BaliScoreReader.
	 */
	public BaliScoreReader() {
		super();
	}
	public BaliScoreReader(String fileBase) {
		setSeqBase(fileBase);
	}

	public static void main(String[] args) {
		
		BaliScoreReader myReader= new BaliScoreReader();
		myReader.setSeqBase("/homes/sammeth/bali/results/cac/ref1/test1/1ajsA.glob");
		myReader.readBoth();
		
		System.out.println(myReader.testAli);
		System.out.println(myReader.refAli);
		System.out.println(myReader.annFile);
		for (int i= 0; i< myReader.sequences.length; ++i)
			System.out.println("("+ myReader.sequences[i].length()+ ") "+ myReader.seqNames[i]+ ": "+ myReader.sequences[i]);
		System.out.print("("+ myReader.scoresGlob.length+ ")");
		for (int i= 0; i< myReader.scoresGlob.length; ++i)
			System.out.print(myReader.scoresGlob[i]);
		System.out.println();
		if (myReader.scoresCore!= null) {
			System.out.print("("+ myReader.scoresCore.length+ ")");
			for (int i= 0; i< myReader.scoresCore.length; ++i)
				System.out.print(myReader.scoresCore[i]);
		}
		System.out.println();
		System.out.println(myReader.spScoreGlob);
		System.out.println(myReader.tcScoreGlob);
		if (myReader.scoresCore!= null) {
			System.out.println(myReader.spScoreCore);
			System.out.println(myReader.tcScoreCore);
		}
	}

	public void setSeqBase(String newAbsFileName) {
		
		StringTokenizer st= new StringTokenizer(newAbsFileName, File.separator);
		absPath= "";
		int count= st.countTokens()- 1;
		for (int i= 0; i< count; ++i) 
			absPath+= File.separator+ st.nextToken();
//		System.out.println(absPath);
		
		st= new StringTokenizer(st.nextToken(), ".");
		fileName= st.nextToken();
//		System.out.println(fileName);
	}
	
	public void readBoth() {
	
			// read global
		readModeGlobal= true;
		try {
			read();
		} catch (Exception e) {
			e.printStackTrace();
		}
	
			// read core
		readModeGlobal= false;
		try {
			read();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
		
	public void read() throws Exception{
		
		String tmp= absPath+ File.separator+ fileName+ ".";
		if (readModeGlobal)
			tmp+= GLOB_EXT;
		else 
			tmp+= CORE_EXT;
		BufferedReader reader= new BufferedReader(
			new FileReader(tmp));
		tmp= "";
		
			// find first line
		while (tmp.equalsIgnoreCase("")) 
			tmp= reader.readLine().trim();
		
			// look for first line
			// (empty core files!)
		while (tmp.indexOf(TEST_ALI_TEXT)== (-1)) {
			if (reader.ready())
				tmp= reader.readLine().trim();
			else
				return;
		}
		testAli= tmp.substring(tmp.indexOf(TEST_ALI_TEXT)+ TEST_ALI_TEXT.length()).trim();

//		System.out.println(testAli);
		tmp= reader.readLine();
		refAli= tmp.substring(tmp.indexOf(REF_ALI_TEXT)+ REF_ALI_TEXT.length()).trim();
//		System.out.println(refAli);
		if (!readModeGlobal) {
			tmp= reader.readLine(); // blank line
			tmp= reader.readLine();
			annFile= tmp.substring(tmp.indexOf(CORE_ANN_TEXT)+ CORE_ANN_TEXT.length()).trim();
//			System.out.println(annFile);
		}
		
			// 2 blank lines
		tmp= "";
		while (tmp.equalsIgnoreCase(""))
			tmp= reader.readLine().trim();

			// read first block
		Vector vecSeq= new Vector();
		Vector vec2= new Vector();
		StringBuffer tmpSeq;
		StringTokenizer st;
		while (!tmp.equalsIgnoreCase("")) {
			
			tmpSeq= new StringBuffer();
			st= new StringTokenizer(tmp);
			
			vec2.add(st.nextToken().trim());
			
			while (st.hasMoreTokens())
				tmpSeq.append(st.nextToken().trim());
			vecSeq.add(tmpSeq);
//			System.err.println(tmpSeq);
			
			tmp= reader.readLine().trim();
		}
			// init arrays
		sequences= new StringBuffer[vecSeq.size()];
		seqNames= new String[vec2.size()];
		for (int i= 0; i< vecSeq.size(); ++i) {
			sequences[i]= (StringBuffer) vecSeq.elementAt(i);
			seqNames[i]= (String) vec2.elementAt(i);
//			System.out.println(seqNames[i]+ ": "+ sequences[i]);
		}
		
			// read scores
		vec2= new Vector();
		tmp= "";
		while (tmp.equalsIgnoreCase(""))
			tmp= reader.readLine().trim();
		st= new StringTokenizer(tmp);
		while (st.hasMoreTokens()) 
			vec2.add(new Integer(Integer.parseInt(st.nextToken().trim())));
		
			// read rest
		tmp= "";
		while (tmp.equalsIgnoreCase(""))
			tmp= reader.readLine().trim();
		while (reader.ready()) {
			
				// break condition
			if (tmp.indexOf(SP_SCORE_TEXT)!= (-1))
				break;
			
				// read in sequences
			for (int i= 0; i< sequences.length; ++i) {
				st= new StringTokenizer(tmp);
				st.nextToken(); 	// skip seq name
				while (st.hasMoreTokens())
					sequences[i].append(st.nextToken().trim());
				
				tmp= reader.readLine().trim();
			}
			
			tmp= "";
			while (tmp.equalsIgnoreCase(""))
				tmp= reader.readLine().trim();
		
				// read scores
				// bug, last block contains empty sequences (no scores!)
			if (tmp.indexOf(SP_SCORE_TEXT)== (-1)) {
				st= new StringTokenizer(tmp);
				while (st.hasMoreTokens()) 
					vec2.add(new Integer(Integer.parseInt(st.nextToken().trim())));

				tmp= "";
				while (tmp.equalsIgnoreCase(""))
					tmp= reader.readLine().trim();
			}
			
		}
		
			// copy scores, read in sp and tc
		if (readModeGlobal) {
			scoresGlob= new int[vec2.size()];
			for (int i= 0; i< scoresGlob.length; ++i)
				scoresGlob[i]= ((Integer) vec2.elementAt(i)).intValue();
			
			spScoreGlob= Double.parseDouble(tmp.substring(tmp.indexOf(SP_SCORE_TEXT)+ SP_SCORE_TEXT.length()).trim());
			
			tmp= reader.readLine().trim();
			while (tmp.indexOf(TC_SCORE_TEXT)== (-1))
				tmp= reader.readLine().trim();
			
			tcScoreGlob= Double.parseDouble(tmp.substring(tmp.indexOf(TC_SCORE_TEXT)+ TC_SCORE_TEXT.length()).trim());

		} else {
			scoresCore= new int[vec2.size()];
			for (int i= 0; i< scoresGlob.length; ++i)
				scoresCore[i]= ((Integer) vec2.elementAt(i)).intValue();
			spScoreCore= Double.parseDouble(tmp.substring(tmp.indexOf(SP_SCORE_TEXT)+ SP_SCORE_TEXT.length()).trim());

			tmp= reader.readLine().trim();
			while (tmp.indexOf(TC_SCORE_TEXT)== (-1))
				tmp= reader.readLine().trim();
			
			tcScoreCore= Double.parseDouble(tmp.substring(tmp.indexOf(TC_SCORE_TEXT)+ TC_SCORE_TEXT.length()).trim());
		}
	}
}
