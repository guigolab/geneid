package qalign.tools;

import java.io.*;

/**
 * A class providing useful functions for multiple alignment
 * layouts.
 *
 * @author Michael Sammeth
 * @version 1.0 rc1
 * @since 1.0 rc1
 */
 
public class LayoutModifier {
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
 * Changes all gap characters within a multiple alignment.
 *
 * @param oldGapChar the gap character to be replaced.
 * @param newGapChar the gap character to be inserted.
 * @param layout the multiple alignment layout.
 * @return the multiple alignment layout with the new gap 
 * character.
 * @since 1.0 rc1
 */
 
public static String[] changeGapChar(String[] layout, char oldGapChar, char newGapChar) {
	
	for (int i= 0; i< layout.length; ++i)
		layout[i]= layout[i].replace(oldGapChar,newGapChar);
	return layout;
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
 * Strips columns consisting of exclusively gaps from the 
 * multiple alignment layout.
 *
 * @param seqs the sequences of the layout.
 * @param quals the qualities according to the sequences of 
 * the layout.
 * @param gap the <code>char</code> representing gaps within
 * the layout.
 * @since 1.0 rc1
 */
 
public static void stripGaps(String[] seqs, int[][] quals, char gap) {
	
	for (int i= 0; i< seqs.length; ++i) {
		for (int j= 0; j< seqs[i].length();++j) {
			if (seqs[i].charAt(j)== gap) {
				seqs[i]= seqs[i].substring(0,j)+ seqs[i].substring(j+1, seqs[i].length());
				int[] tempQuals= quals[i];
				quals[i]= new int[tempQuals.length-1];
				for (int k= 0; k< quals[i].length; ++k)
					if (k< j)
						quals[i][k]= tempQuals[k];
					else if (k> j)
						quals[i][k-1]= tempQuals[k];

				--j;
			}
		}
	}
	
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
 
public static boolean writeMSF(String[] layout, String[] names) {

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
		
	try {
		BufferedWriter bw= new java.io.BufferedWriter(new FileWriter("tst.tst"));
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
			
		
		
	} catch (IOException e) {
		System.err.println(e);
	}
	
	return false;
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
 
public static boolean writeMSF(String[] layout, String[] names, long[] checks, double[] weights) {

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
		
	try {
		BufferedWriter bw= new java.io.BufferedWriter(new FileWriter("tst.tst"));
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
			
		
		
	} catch (IOException e) {
		System.err.println(e);
	}
	
	return false;
}
}
