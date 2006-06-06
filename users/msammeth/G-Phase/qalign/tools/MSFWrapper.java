package qalign.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class MSFWrapper extends DefaultIOWrapper {

	// wrapped information
	/**
	 * alignment type (optional)
	 * 	0= unknown (!! line not first or missing, nor Type identifier)
	 *  1= nucleic (!!NA_MULTIPLE_ALIGNMENT 1.0 or Type: N)
	 *  2= protein (!!AA_MULTIPLE_ALIGNMENT 1.0 or Type: P)
	 */
	protected byte type = 0;
	protected double version = 0d;

	/**
	 * description (optional)
	 */
	protected String description = null;

	/**
	 * dividing line identifiers
	 * (length, checksum required)
	 */
	protected String fileName = null;
	protected int length = 0;
	protected String timeStamp = null;
	protected int checksum = 0;

	/**
	 * sequence identifiers
	 * (required)
	 */
	protected String[] seqNames = null;
	protected int[] seqLengths = null;
	protected int[] seqChecks = null;
	protected float[] seqWeights = null;

	/**
	 * sequences
	 */
	protected String[] sequences = null;
	protected boolean numbersDrawn= false;
	protected int blockSize = 0;

	/**
 * Calculates the GCG checksum for a sequence of the layout.
 *
 * @param seq the sequence.
 * @return the GCG checksum.
 * @since 1.0 rc1
 */
 
public static long calcGCGchecksum(String seq) {

	int check = 0, count = 0;
	byte[] ba = seq.getBytes();
	
	for (int i = 0; i < seq.length(); i++) {
		byte b = ba[i];
		if (b >= 'a' && b <= 'z')
			b -= 32;
		count++;
		check += count * b;
		if (count == 57)
			count = 0;
	}
	check %= 10000;
	return check;
}


	/**
 * Calculates the GCG checksum for an given substring of a 
 * sequence of the layout.
 *
 * @param seq the sequence.
 * @param offset where to start coding the GCG checksum.
 * @param seqlen the length of the substring to be encoded
 * to the GCG checksum.
 * @return the GCG checksum.
 * @since 1.0 rc1
 */
 
public static long calcGCGchecksum(String seq, int offset, int seqlen) {

	int check = 0, count = 0;
	byte[] ba = seq.getBytes();
	
	for (int i = 0; i < seqlen; i++) {
		byte b = ba[i + offset];
		if (b >= 'a' && b <= 'z')
			b -= 32;
		count++;
		check += count * b;
		if (count == 57)
			count = 0;
	}
	check %= 10000;
	return check;
}


	/**
 * Encodes the multiple alignment to MSF-format.
 *
 * @param names the names of the sequences within the 
 * multiple alignment layout.
 * @param layout the multiple alignment layout.
 * @return the string representation of the MSF format
 * encoded multiple alignment layout.
 * @since 1.0 rc1
 */
 
public static String getMSFStr(String[] layout, String[] names) {

	int maxLineLen= 80;

	java.util.Random rnd= new java.util.Random();
	float[] weights= new float[layout.length];
	for (int i= 0; i< weights.length; ++i)
		weights[i]= rnd.nextFloat();
		
	long[] checks= new long[layout.length];
	int sumCRC= 0;							// sum up CRCs
	for (int i= 0; i< checks.length; ++i) {
		checks[i]= calcGCGchecksum(layout[i]);
		sumCRC+= checks[i];
	}

	int nameMax= 0;
	int nameMaxNb= -1;
	for (int i= 0; i< names.length; ++i)
		if (names[i].length()> nameMax) {
			nameMax= names[i].length();
			nameMaxNb= i;
		}

	String result= "";
	
		result+= "PileUp\tMSF: "+layout[0].length()+"\tType: N\tCheck:"+sumCRC+"\t..\n";
		result+= "\n";
		for (int i= 0; i< names.length; ++i)
			result+= "Name: "+names[i]+"\tLen: "+layout[i].length()+"\tCheck: "+checks[i]+"\tWeight: "+weights[i]+"\n";
		result+="\n//\n";
		
		for (int x= 0; x< layout[0].length(); x+= maxLineLen- nameMax) {
			for (int i= 0; i< layout.length; ++i) {
				result+= names[i];
				for (int j= names[i].length(); j< nameMax; ++j)
					result+= " ";
				if (x+(maxLineLen- nameMax)< layout[i].length())
					result+= layout[i].substring(x, x+(maxLineLen- nameMax));
				else
					result+= layout[i].substring(x, layout[i].length());
				result+= "\n";
			}
		}
			
		
		
	
	return result;
}


	/**
 * Encodes the multiple alignment to MSF-format and writes it to a file.
 *
 * @param names the names of the sequences within the 
 * multiple alignment layout.
 * @param layout the multiple alignment layout.
 * @return the string representation of the MSF format
 * encoded multiple alignment layout.
 * @since 1.0 rc1
 */
 
public boolean writeMSF(String[] layout, String[] names) throws Exception {

	int maxLineLen= 80;

//	java.util.Random rnd= new java.util.Random();
	float[] weights= new float[layout.length];
	for (int i= 0; i< weights.length; ++i)
		weights[i]= 1f; // rnd.nextFloat();
		
		// compute checksums
	long[] checks= new long[layout.length];
	int sumCRC= 0;							// sum up CRCs
	for (int i= 0; i< checks.length; ++i) {
		checks[i]= calcGCGchecksum(layout[i]);
		sumCRC+= checks[i];
	}

		// find longest name
	int nameMax= 0;
	int nameMaxNb= -1;
	for (int i= 0; i< names.length; ++i)
		if (names[i].length()> nameMax) {
			nameMax= names[i].length();
			nameMaxNb= i;
		}
	++nameMax;	// an extra space
		
//	try {
		BufferedWriter bw= new java.io.BufferedWriter(new FileWriter(fPath+File.separator+fName));
//		bw.write("!!NA_MULTIPLE_ALIGNMENT 1.0\n");
		bw.write("PileUp\tMSF: "+layout[0].length()+"\tType: N\tCheck:"+sumCRC+"\t..\n\n");
		if (description!= null)
			bw.write(description+ "\n\n");
		for (int i= 0; i< names.length; ++i) {
			bw.write("Name: "+names[i]);
			for (int j= names[i].length(); j< nameMax; ++j)
				bw.write(" ");
			bw.write("\tLen: "+layout[i].length()+"\tCheck: "+checks[i]+"\tWeight: "+weights[i]+"\n");
		}
		bw.write("\n//\n");
		
		for (int x= 0; x< layout[0].length(); x+= maxLineLen- nameMax) {
			for (int i= 0; i< layout.length; ++i) {
				bw.write(names[i]);
				for (int j= names[i].length(); j< nameMax; ++j)
					bw.write(" ");
				if (x+(maxLineLen- nameMax)< layout[i].length())
					bw.write(layout[i].substring(x, x+(maxLineLen- nameMax)));
				else
					bw.write(layout[i].substring(x, layout[i].length()));
				bw.write("\n");
			}
			bw.write("\n");
		}
		
		bw.close();
			
		
/*	} catch (IOException e) {
		System.err.println(e);
	}
*/	
	return false;
}

	public void write() throws Exception {

		writeMSF(sequences, seqNames);		
	}


	/**
 * Encodes the multiple alignment to MSF-format and writes it to a file.
 * Additionaly writes the <code>check</code> and </code>weight</code>
 * attributes (mostly unused) given.
 *
 * @param names the names of the sequences within the 
 * multiple alignment layout.
 * @param layout the multiple alignment layout.
 * @param checks check for each sequence.
 * @param weights weight for each sequence.
 * @return the string representation of the MSF format
 * encoded multiple alignment layout.
 * @since 1.0 rc1
 */
 
public boolean writeMSF(String[] layout, String[] names, long[] checks, double[] weights) throws Exception {

	int maxLineLen= 80;
	
	int sumCRC= 0;							// sum up CRCs
	for (int i= 0; i< checks.length; ++i)
		sumCRC+= checks[i];

	int nameMax= 0;
	int nameMaxNb= -1;
	for (int i= 0; i< names.length; ++i)
		if (names[i].length()> nameMax) {
			nameMax= names[i].length();
			nameMaxNb= i;
		}
		
//	try {
		BufferedWriter bw= new java.io.BufferedWriter(new FileWriter(fPath+ File.separator+ fName));
//		bw.write("!!NA_MULTIPLE_ALIGNMENT 1.0\n");
		bw.write("PileUp\tMSF: "+layout[0].length()+"\tType: N\tCheck:"+sumCRC+"\t..\n");
		bw.write("\n");
		for (int i= 0; i< names.length; ++i)
			bw.write("Name: "+names[i]+"\tLen: "+layout[i].length()+"\tCheck: "+checks[i]+"\tWeight: "+weights[i]+"\n");
		bw.write("\n//\n");
		
		for (int x= 0; x< layout[0].length(); x+= maxLineLen- nameMax) {
			for (int i= 0; i< layout.length; ++i) {
				bw.write(names[i]);
				for (int j= names[i].length(); j< nameMax; ++j)
					bw.write(" ");
				if (x+(maxLineLen- nameMax)< layout[i].length())
					bw.write(layout[i].substring(x, x+(maxLineLen- nameMax)));
				else
					bw.write(layout[i].substring(x, layout[i].length()));
				bw.write("\n");
			}
		}
			
		bw.close();
		
/*	} catch (IOException e) {
		System.err.println(e);
	}
*/	
	return false;
}

	protected int blocksPerLine = 0;


	public MSFWrapper(String newFName, String newFPath) {
		
		super(newFName, newFPath);
	}
	
	public MSFWrapper(String absFName) {
		
		super(absFName);
	}

	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void initMSF() throws Exception {

		BufferedReader fptr = null;
//		try {
			fptr = new BufferedReader(new FileReader(fPath+ File.separator+ fName));

			readHeader(fptr);
			readDescriptor(fptr);
			readSequences(fptr);

/*		} catch (Throwable e) {
			e.printStackTrace();
			// do nothing
		}
*/
	}
	
	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void read() throws Exception {
		
		initMSF();
	}

	public void readHeader(BufferedReader fptr) throws Exception {

		Exception ex;

		// find Dividing Line
		String readStr = "";
		String aStr = "";
		while ((fptr.ready())
			&& ((readStr.indexOf("MSF:") == (-1))
//				|| (readStr.indexOf("Check:") == (-1)))	not good for DIALIGN
				)) {

			// file type spotted
			if ((readStr.startsWith("!!"))
				&& (readStr.indexOf("A_MULTIPLE_ALIGNMENT") != (-1))) {
				if (readStr.charAt(2) == 'N')
					type = 1;
				if (readStr.charAt(2) == 'A')
					type = 2;
				readStr =
					readStr
						.substring(readStr.indexOf(' '), readStr.length())
						.trim();
				if (readStr.length() > 0)
					version = Double.parseDouble(readStr);

				// add to comment
			} else
				aStr += readStr.trim() + " ";

			readStr = fptr.readLine();
		}

		//  save comment (if found)
		if (aStr.trim().length() > 0)
			description = aStr;

		// check for dividing line was found
		if ((readStr.indexOf("MSF:") == (-1))
//			|| (readStr.indexOf("Check:") == (-1))
			) {
				// error("Dividing Line not found!");
			throw new Exception("Warning: Dividing Line not found.");
		} else {

			StringTokenizer tokenizer = new StringTokenizer(readStr, " ");
			aStr = "";
			while (tokenizer.hasMoreTokens()) {

				aStr = tokenizer.nextToken();

				if (aStr.trim().equals("MSF:")) {
					if (tokenizer.hasMoreTokens()) {
						try {
							length = Integer.parseInt(tokenizer.nextToken());
						} catch (NumberFormatException e) {
								// error(e.getMessage());
							ex= new Exception("Warning: Length could not be read.");
// 1.4						ex.setStackTrace(e.getStackTrace());
							throw e;
						}
					}
					continue;
				}

				// before MSF: <length>
				if (length == 0) {
					fileName = aStr;
					continue;
				}

				if (aStr.trim().equals("Type:")) {
					if (tokenizer.hasMoreTokens()) {
						aStr = tokenizer.nextToken();
						if (aStr.trim().equalsIgnoreCase("N")) {
							if ((type != 0) && (type != 1)) {
									// error("Divergent Type Descriptors");
								throw new Exception("Warning: Divergent Type Descriptors.");
							} else
								type = 1;
							continue;
						}
						if (aStr.trim().equalsIgnoreCase("P")) {
							if ((type != 0) && (type != 2)) {
									// error("Divergent Type Descriptors");
								throw new Exception("Warning: Divergent Type Descriptors.");
							} else
								type = 2;
							continue;
						}
							// error("Unknown Type");
						throw new Exception("Warning: Type is not recognised.");
					}
					continue;
				}

				if (aStr.trim().equals("Check:")) {
					if (tokenizer.hasMoreTokens()) {
						try {
							checksum = Integer.parseInt(tokenizer.nextToken());
						} catch (NumberFormatException e) {
							// error(e.getMessage());
							ex= new Exception("Warning: Checksum could not be retrieved");
// 1.4							ex.setStackTrace(e.getStackTrace());
							throw ex;
						}
					}
					continue;
				}

				if (aStr.trim().equals(".."))
					break;

				if (timeStamp == null)
					timeStamp = "";
				timeStamp += aStr;
			}
		}
	}

	public void readDescriptor(BufferedReader fptr) throws Exception {

		String readStr = "";
		while ((fptr.ready()) && (readStr.indexOf("Name:") == (-1))) {
			if (readStr.trim().length() > 0)
				//error("Undefined Description Found:\n" + readStr);
				throw new Exception(readStr);
			readStr= fptr.readLine();
		}	
		Exception ex;

		// parse sequence description
		seqNames = new String[0];
		seqLengths = new int[0];
		seqChecks = new int[0];
		seqWeights = new float[0];
		while ((fptr.ready()) && (readStr.indexOf("//") == (-1))) {
			
				// skip blank lines
			if (readStr.trim().length()== 0) {
				if (fptr.ready())
					readStr= fptr.readLine();
				continue;
			}
			
			seqNames= extentStringArray(seqNames);
			seqLengths= extentIntArray(seqLengths);
			seqChecks= extentIntArray(seqChecks);
			seqWeights= extentFloatArray(seqWeights);
			
			StringTokenizer tokenizer = new StringTokenizer(readStr, " \t");
			String aStr = "";
			while (tokenizer.hasMoreTokens()) {

				aStr = tokenizer.nextToken();
				if ((aStr.trim().equals("Name:"))
					&& (tokenizer.hasMoreTokens()))
					seqNames[seqNames.length - 1] = tokenizer.nextToken();
				if ((aStr.trim().equals("Len:"))
					&& (tokenizer.hasMoreTokens())) {
					try {
						seqLengths[seqLengths.length - 1] =
							Integer.parseInt(tokenizer.nextToken());
					} catch (NumberFormatException e) {
						// error(e.getMessage());
						ex= new Exception("Warning: Length could not be retrieved");
// 1.4						ex.setStackTrace(e.getStackTrace());
						throw ex;
					}
				}

				if ((aStr.trim().equals("Check:"))
					&& (tokenizer.hasMoreTokens())) {
					try {
						seqChecks[seqChecks.length - 1] =
							Integer.parseInt(tokenizer.nextToken());
					} catch (NumberFormatException e) {
						// error(e.getMessage());
						ex= new Exception("Warning: Checksum could not be retrieved");
// 1.4						ex.setStackTrace(e.getStackTrace());
						throw ex;
					}
				}

				if ((aStr.trim().equals("Weight:"))
					&& (tokenizer.hasMoreTokens())) {
					try {
						seqWeights[seqWeights.length - 1] =
							Float.parseFloat(tokenizer.nextToken()); 
					} catch (NumberFormatException e) {
						// error(e.getMessage());
						ex= new Exception("Warning: Weight could not be retrieved");
// 1.4						ex.setStackTrace(e.getStackTrace());
						throw ex;
					}
				}
			}

			readStr = fptr.readLine();
		}
	}
	
	public void readSequences(BufferedReader fptr) throws IOException {
		
		String readStr= "";
			// skip leading blank lines
		while ((fptr.ready())&& (readStr.trim().length()== 0))
			readStr= fptr.readLine();

		if (!fptr.ready())
			readStr= null;
		sequences= extentStringArray(sequences);
		int count= 0;
		while (readStr!= null) {
			
				
				// empty line
			if (readStr.trim().length()== 0) {
				while ((fptr.ready())&& (readStr.trim().length()== 0))
					readStr= fptr.readLine();
				count= 0;
			}			
				// reader ran empty
			if (readStr.trim().length()== 0)
				break;
				
				
				// parse line
			if (!(count< sequences.length))
				sequences= extentStringArray(sequences);
			StringTokenizer tokenizer= new StringTokenizer(readStr.trim(), " ");
				
				// read name
			String aStr= tokenizer.nextToken();
//			if (!aStr.equalsIgnoreCase(seqNames[sequences.length- 1]))
//					error("names and sequences dont match in "+ (sequences.length- 1));
			if (tokenizer.hasMoreTokens())
				aStr= tokenizer.nextToken();
			try {
				Integer.parseInt(aStr);
				numbersDrawn= true;
				continue;
			} catch (NumberFormatException e) {
			
				if (sequences[count]== null)
					sequences[count]= "";
				sequences[count]+= aStr;
				while (tokenizer.hasMoreTokens())
					sequences[count]+= tokenizer.nextToken();

			}
			
			if (fptr.ready())
				readStr= fptr.readLine();
			else
				readStr= null;
			++count;
		}
	}
	
	public String toString() {
		
		String aStr= "";
		aStr+= "File: "+fileName+"\n";
		aStr+= "Comment: "+description+"\n";
		
		if (fileName!= null)
			aStr+= fileName+" ";
		aStr+= "MSF: "+length+" ";
		if (type!= 0) {
			aStr+= "Type: ";
			if (type== 1)
				aStr+= "N";
			if (type== 2)
				aStr+= "P";
			aStr+= " ";
		}
		if (timeStamp!= null)
			aStr+= timeStamp;
		
		aStr+= "Check: "+checksum+" ..\n\n";
	
		if (seqNames!= null)	
		for (int i= 0; i< seqNames.length; ++i) {
			
			aStr+= "Name: "+seqNames[i]+" ";
			aStr+= "Len: "+seqLengths[i]+" ";
			aStr+= "Check: "+seqChecks[i]+" ";
			aStr+= "Weight: "+seqWeights[i]+" ";
			aStr+="\n";
		}
		
		aStr+="\n\n//\n\n";
		
		if((sequences!= null)&& (sequences.length> 0))
		for (int k= 0; k< sequences.length; ++k) 
				aStr+= seqNames[k]+ " "
					+ sequences[k]+"\n";
		
		return aStr;
	}
	
	
	public void error(String message) {
		System.err.println(message);
	}
	/**
	 * Returns the sequences.
	 * @return String[]
	 */
	public String[] getSequences() {
		return sequences;
	}

	/**
	 * Sets the sequences.
	 * @param sequences The sequences to set
	 */
	public void setSequences(String[] sequences) {
		this.sequences = sequences;
	}

	/**
	 * Returns the seqNames.
	 * @return String[]
	 */
	public String[] getSeqNames() {
		return seqNames;
	}

	/**
	 * Sets the seqNames.
	 * @param seqNames The seqNames to set
	 */
	public void setSeqNames(String[] seqNames) {
		this.seqNames = seqNames;
	}
	
	public SequenceWrapper[] getWrappedSequences() {
		
		SequenceWrapper[] result= new SequenceWrapper[sequences.length];
		for (int i= 0; i< result.length; ++i) {
			result[i]= new SequenceWrapper(
				seqNames[i], 
				"length: "+seqLengths[i]+", weight: "+seqWeights[i],
				fName,
				getFDate()
			);
			result[i].setSequence(sequences[i]);
			result[i].setWrapper(this);
		}
		
		return result;
	}
	
	public void setWrappedSequences(SequenceWrapper[] newData) {
		
		sequences= new String[newData.length];
		seqNames= new String[newData.length];
		
		for (int i= 0; i< newData.length; ++i) {
			
			sequences[i]= newData[i].getSequence();
			seqNames[i]= newData[i].getName();
		}
	}

	/**
	 * Returns the description.
	 * @return String
	 */
	public String getDescription() {
		return description;
	}

	/**
	 * Sets the description.
	 * @param description The description to set
	 */
	public void setDescription(String description) {
		this.description = description;
	}

}