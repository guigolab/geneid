

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class FASTAUniter {

	protected String[] sequences= null;
	protected String[] seqNames= null;
	protected String[] seqExtNames= null;
	String fName, fPath;
	
/**
 * Line length for FASTA file output.
 *
 * @since 1.0 rc1 
 */
	protected static int FASTA_FILE_SEQ_LINE = 80;
	
	public FASTAUniter(String newFName, String newFPath) {
		
		this.fPath= newFPath; 
		this.fName= newFName;
	}
	
	public FASTAUniter(String absFilePath) {
		int p= absFilePath.length();
		while ((--p>= 0)&& (absFilePath.charAt(p)!= File.separatorChar))
			; // count down
		this.fPath= absFilePath.substring(0,p);
		this.fName= absFilePath.substring((p+1), absFilePath.length());
	}
	
	/**
	 * @throws all exceptions during file handling to keep <code>this</code> independent from UI
	 */
	public void read() throws Exception {
		
		BufferedReader fptr = null;

		StringBuffer sb= null;
		try {
			fptr = new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			String read= null;
			String seq= "";
			sequences= new String[0];
			seqNames= new String[0];
			seqExtNames= new String[0];
			while (fptr.ready()) {
				
				read= fptr.readLine().trim();
		
				if (read.startsWith(">")) {
					if (sequences.length> 0) 
						sequences[sequences.length- 1]= sb.toString();
					sb= new StringBuffer();
					sequences= extentStringArray(sequences);
					seqNames= extentStringArray(seqNames);
					read= read.substring(1).trim();
					seqNames[seqNames.length- 1]= read.substring(0, read.length());
				} else
					sb.append(read);
			}
			fptr.close();
		} catch (Exception e) {
			System.err.println("PROBLEMS with reading file "+fName);
			throw e;
		}
		sequences[sequences.length- 1]= sb.toString();
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

	public void setWrapLine(int x) {
		FASTA_FILE_SEQ_LINE= x;
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

	public static String[] extentStringArray(String[] old) {
	
		if (old== null) {
			String[] extended= new String[0];
			old= extended;
			return extended;
		} else {
			String[] extended= new String[old.length+ 1];
			for (int i= 0; i< old.length; ++i)
				extended[i]= old[i];
			old= extended;
			return extended;
		}
	}

	public static void main(String[] args) {
		if (args.length< 1) {
			System.out.println("Usage: Uniter separator \ne.g. FASTAUniter .");
			System.exit(0);
		}
	
		String sep= args[0];
		String path= "./";
		
		int wrap= 50;
		if (args.length> 1)
			wrap= Integer.parseInt(args[1]);
		
		String[] list= new File(path).list();
		HashMap map= new HashMap();
		int ignoreCnt= 0;
		for (int i = 0; i < list.length; i++) {
			int p= list[i].lastIndexOf(sep);
			if (p< 0) {
				++ignoreCnt;
				continue;
			}
			String sfx= list[i].substring(p+1, list[i].length());	// ignore
			if (sfx.equalsIgnoreCase("fa")|| sfx.equalsIgnoreCase("class"))
				continue;
			String pfx= list[i].substring(0, p);
			HashMap m= (HashMap) map.get(pfx);
			if (m== null) 
				m= new HashMap();
			m.put(new Integer(Integer.parseInt(sfx)), list[i]);
			map.put(pfx, m);
		}

		Object[] o= map.keySet().toArray();		
		String[] keys= new String[o.length];
		for (int i = 0; i < keys.length; i++) 
			keys[i]= (String) o[i];
		for (int i = 0; i < keys.length; i++) 
			try {
				HashMap m= (HashMap) map.get(keys[i]);
				Object[] kk= m.keySet().toArray();
				Arrays.sort(kk);
				System.out.println("processing "+keys[i]+", "+o.length+" files.");
				list= new String[kk.length];
				for (int j = 0; j < list.length; j++) 
					list[j]= (String) m.get(kk[j]);
				StringBuffer sb= new StringBuffer();
				System.out.print("\treading..");
				for (int j = 0; j < list.length; j++) {
					FASTAUniter fasta= new FASTAUniter(new File(path+list[j]).getAbsolutePath());
					fasta.read();
					sb.append(fasta.getSequences()[0]);
				}
				System.out.print("\twriting..");
				FASTAUniter fasta= new FASTAUniter(new File(path+keys[i]+".fa").getAbsolutePath());
				fasta.seqNames= new String[] {keys[i]};
				fasta.sequences= new String[] {sb.toString()};
				FASTA_FILE_SEQ_LINE= 50;
				fasta.write();
				System.out.println("done.");
			} catch (Exception e) {
				e.printStackTrace();
				System.gc();
				System.runFinalization();
//				try {
//					System.in.read();
//				} catch (IOException e1) {;}
				Thread.yield();
				--i;
				continue;
			}
		
		
		System.out.println(keys.length+ " files created, "+ ignoreCnt+ " files unchanged.");
	}
}
