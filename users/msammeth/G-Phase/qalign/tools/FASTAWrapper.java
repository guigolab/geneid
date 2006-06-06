package qalign.tools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class FASTAWrapper extends DefaultIOWrapper {

	protected String[] sequences= null;
	protected String[] seqNames= null;
	protected String[] seqExtNames= null;

/**
 * Line length for FASTA file output.
 *
 * @since 1.0 rc1
 */
	protected static int FASTA_FILE_SEQ_LINE = 80;
	
	public FASTAWrapper(String newFName, String newFPath) {
		
		super(newFName, newFPath);
	}
	
	public FASTAWrapper(String absFName) {
		super(absFName);
	}
	
	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void read() throws Exception {
		
		BufferedReader fptr = null;

		fptr = new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			
		String read= null;
		String seq= "";
		sequences= new String[0];
		seqNames= new String[0];
		seqExtNames= new String[0];
		while (fptr.ready()) {
			
			read= fptr.readLine().trim();

			if (read.startsWith(">")) {
				seqNames= extentStringArray(seqNames);
				seqExtNames= extentStringArray(seqExtNames);
				read= read.substring(1).trim();
				if (read.indexOf(' ')> 0) {
					seqNames[seqNames.length- 1]= read.substring(0, read.indexOf(' '));
					seqExtNames[seqExtNames.length- 1]= '>'+ read.substring(read.indexOf(' '), read.length());
/*				} else if (read.indexOf('_')> 0) {
					seqNames[seqNames.length- 1]= read.substring(1, read.indexOf('_'));
					seqExtNames[seqExtNames.length- 1]= '>'+ read.substring(read.indexOf('_'), read.length());
*/				} else {
					seqNames[seqNames.length- 1]= read.substring(0, read.length());
					seqExtNames[seqExtNames.length- 1]= read.substring(0, 1);
				}
				if (sequences.length> 0) {
					sequences[sequences.length- 1]= seq;
					seq= "";
				}
				sequences= extentStringArray(sequences);
			} else
				seq+= read;
		}
		sequences[sequences.length- 1]= seq;
			
	}
	
/**
 * Writes a FASTA-file with the given sequences.
 *
 * @param fname the name of the FASTA file.
 * @param seqs the sequences.
 * @param seqNames the names of the sequences.
 * @return <true> if the file could be written successfully,
 * <false> else.
 * @since 1.0 rc1
 */
 
public boolean writeFASTA(String[] seqNames, String[] seqs) throws Exception {
	
	BufferedWriter fptr;

//	try {
		fptr= new BufferedWriter(new FileWriter(fPath+ File.separator+ fName));
/*	} catch (IOException e) {
		System.err.println("ERROR: cant create file "+filePath+"\n"+e);
		return false;
	}
*/
	if (seqNames.length!= seqs.length) 
		throw new Exception("ERROR: name count doesnt match seq count!");

//	try {
		for (int i= 0; i< seqs.length; i++) {
			fptr.write(">"+seqNames[i]+"\n");
			for (int j= 0; j< seqs[i].length();j+= FASTA_FILE_SEQ_LINE) {
				if (FASTA_FILE_SEQ_LINE+j> seqs[i].length())
					fptr.write(seqs[i].substring(j,seqs[i].length())+"\n\n");
				else
					fptr.write(seqs[i].substring(j,j+FASTA_FILE_SEQ_LINE)+"\n");
			}
			
		}
		fptr.close();
/*	} catch (IOException e) {
		System.err.println(e);
		return false;
	}
*/
	
	return true;
	
}

public void write() throws Exception {

	writeFASTA(seqNames, sequences);
}

	public void setWrappedSequences(SequenceWrapper[] newData) {
		
		sequences= new String[newData.length];
		seqNames= new String[newData.length];
		
		for (int i= 0; i< newData.length; ++i) {
			
			sequences[i]= newData[i].getSequence();
			seqNames[i]= newData[i].getName();
		}
	}	
	
	public SequenceWrapper[] getWrappedSequences() {
		
		SequenceWrapper[] result= new SequenceWrapper[sequences.length];
		for (int i= 0; i< result.length; ++i) {
			result[i]= new SequenceWrapper(
				seqNames[i], 
				seqExtNames[i],
				fName,
				getFDate()
			);
			result[i].setSequence(sequences[i]);
			result[i].setWrapper(this);
		}
		
		return result;
	}	
	/**
	 * @return Returns the seqExtNames.
	 */
	public String[] getSeqExtNames() {
		return seqExtNames;
	}
	/**
	 * @return Returns the seqNames.
	 */
	public String[] getSeqNames() {
		return seqNames;
	}
	/**
	 * @return Returns the sequences.
	 */
	public String[] getSequences() {
		return sequences;
	}
}
