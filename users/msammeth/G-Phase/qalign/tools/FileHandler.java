package qalign.tools;

import java.util.*;
import java.io.*;

/**
 * A class which provides tools for I/O handling of the 
 * FASTA file format.
 *
 * @author Michael Sammeth
 * @version 1.0 rc1
 * @since 1.0 rc1
 */
 
public class FileHandler {
/**
 * Line length for FASTA file output.
 *
 * @since 1.0 rc1
 */
	protected static int FASTA_FILE_SEQ_LINE = 80;
/**
 * Reads the names of the FASTA-file discarding the 
 * leading <i>'>'</i> character.
 *
 * @param fname the name of the FASTA file.
 * @return the names of the sequences in the specified
 * FASTA file.
 * @since 1.0 rc1
 */
 
public static String[] readFASTAnames(String fname) {
	
	BufferedReader fptr= null;
	String s;
	Vector	nameRaw= new Vector();
	String[] names;
	
	try {
		fptr= new BufferedReader(new FileReader(fname));
	} catch (FileNotFoundException e) {
		System.err.println(e);
	}
	
	try {
		while (fptr.ready()) {
			s= fptr.readLine().trim();
			if (s.startsWith(">"))
				nameRaw.add(s.substring(1, s.length()));
		}
		fptr.close();
	} catch (IOException e) {
		System.err.println(e);
	}
			
	names= new String[nameRaw.size()];
	for (int i= 0; i< names.length; ++i)
		names[i]= (String) nameRaw.elementAt(i);
		
	return names;
}
/**
 * Reads the sequences from the FASTA-file.
 *
 * @param fname the name of the FASTA file.
 * @return the sequences within the specified
 * FASTA file.
 * @since 1.0 rc1
 */
 
public static String[] readFASTAseqs(String fname) {
	
	BufferedReader fptr= null;
	String s;
	String l= null;
	Vector	seqRaw= new Vector();
	String[] seqs;
	
	try {
		fptr= new BufferedReader(new FileReader(fname));
	} catch (FileNotFoundException e) {
		System.err.println(e);
	}
	
	try {
		while (fptr.ready()) {
			s= fptr.readLine().trim();
			if (s.startsWith(">")) {
				if (l!= null)
					seqRaw.add(l);
				l= "";
			} else
				l+= s;
		}
		seqRaw.add(l);
		fptr.close();
	} catch (IOException e) {
		System.err.println(e);
	}
			
	seqs= new String[seqRaw.size()];
	for (int i= 0; i< seqs.length; ++i)
		seqs[i]= (String) seqRaw.elementAt(i);
		
	return seqs;
}
/**
 * Writes a FASTA-file with the given sequences.
 *
 * @param fname the name of the FASTA file.
 * @param seqs array of <code>Sequence</code> instances.
 * @return <true> if the file could be written successfully,
 * <false> else.
 * @since 1.0 rc1
 */
 
public static boolean writeFASTA(String fname, qalign.object.Sequence[] seqs) {

	String[] sname= new String[seqs.length];
	String[] sseqs= new String[seqs.length];

	for (int i= 0; i< seqs.length; i++) {
		sname[i]= seqs[i].getName();
		sseqs[i]= new String(seqs[i].getSequence());
	}
	
	return writeFASTA(fname, sname, sseqs);
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
 
public static boolean writeFASTA(String filePath, String[] seqNames, String[] seqs) {
	
	BufferedWriter fptr;

	try {
		fptr= new BufferedWriter(new FileWriter(filePath));
	} catch (IOException e) {
		System.err.println("ERROR: cant create file "+filePath+"\n"+e);
		return false;
	}

	if (seqNames.length!= seqs.length) {
		System.err.println("ERROR: name count doesnt match seq count!");
		return false;
	}

	try {
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
	} catch (IOException e) {
		System.err.println(e);
		return false;
	}

	
	return true;
	
}


}
