/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


import gphase.Toolbox;
import gphase.tools.File;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Vector;
import java.util.zip.Checksum;

/**
 * Utility class that calculates a CRC64 checksum on a stream of bytes. Code was
 * based on some from BioPerl. Note that we use longs then cast them to avoid
 * the lack of an unsigned int in Java. Longs are 64-bit but we are only using
 * the bottom 32 bits. An int is 32-bit but encodes sign so we can get amusing
 * results if we don't allow for this.
 * 
 * @author Unknown. Copied from Expasy4J for convenience. See <a
 *         href="http://dev.isb-sib.ch/projects/expasy4j/">http://dev.isb-sib.ch/projects/expasy4j/</a>
 * @since 1.5
 */
public class CRC64Checksum implements Checksum {
	private static final long POLY64 = 0xD800000000000000L;

	private static final long[] crcTable = new long[256];

	private long crc;

	public static void main(String[] args) {
		File f= new File("/home/msammeth/sequences/pep_mRNA_24_fromUCSC_070731.fa");
		_070809_add_CRC64_to_Fasta(f);
	}
	
	static void _070809_add_CRC64_to_Fasta(File fastaFile) {
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fastaFile));
			File outFile= new File(fastaFile.getPathOnly()+ File.separator+ fastaFile.getFileNameWithoutExtension()+"_CRC64."+fastaFile.getExtension());
			Toolbox.checkFileExists(outFile.getAbsolutePath());
			Vector<String> seqVector= new Vector<String>();
			String label= null;		
			CRC64Checksum crc= null;
			long size= fastaFile.length();
			long bytesRead= 0l;
			int ctr= 0;
			
			while (buffy.ready()) {
				String line= buffy.readLine().trim();
				bytesRead+= line.length()+ 1;
				if (10* bytesRead/ size> ctr) {
					++ctr;
					System.out.print("*");
					System.out.flush();
				}
				if (line== null|| line.length()== 0)
					continue;
				
				if (line.charAt(0)== '>') {
					if (label!= null) {
						BufferedWriter writer= new BufferedWriter(new FileWriter(outFile, true));
						writer.write(label+" "+crc.toString()+"\n");
						for (int i = 0; i < seqVector.size(); i++) 
							writer.write(seqVector.elementAt(i)+"\n");
						writer.flush(); writer.close();
					}
					label= line;
					crc= new CRC64Checksum();
					seqVector= new Vector<String>();
					continue;
				}
				seqVector.add(line);
				char[] c= line.toUpperCase().toCharArray();
				byte[] b= new byte[c.length];
				for (int i = 0; i < b.length; i++) 
					b[i]= (byte) c[i];
				crc.update(b, 0, b.length);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static {
		for (int i = 0; i < 256; ++i) {
			long part = i;
			for (int j = 0; j < 8; ++j)
				part = ((part & 1) != 0) ? (part >>> 1) ^ POLY64 : (part >>> 1);
			crcTable[i] = part;
		}
	}

	public void update(int b) {
		long low = crc >>> 8;
		long high = crcTable[(int) ((crc ^ b) & 0xFF)];
		crc = low ^ high;
	}

	public void update(byte[] b, int offset, int length) {
		for (int i = offset; i < length; ++i)
			update(b[i]);
	}

	public void update(String s) {
		// update(s.getBytes(), 0, s.length());
		int size = s.length();
		for (int i = 0; i < size; ++i)
			update(s.charAt(i));

	}

	public long getValue() {
		return crc;
	}

	/**
	 * Returns a zero-padded 16 character wide string containing the current
	 * value of this checksum in uppercase hexadecimal format.
	 */
	public String toString() {
		StringBuffer builder = new StringBuffer();
		builder.append(Long.toHexString(crc >>> 4));
		builder.append(Long.toHexString(crc & 0xF));
		for (int i = 16 - builder.length(); i > 0; --i)
			builder.insert(0, '0');
		return builder.toString().toUpperCase();
	}

	public void reset() {
		crc = 0;
	}
}