

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Vector;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class FASTAReformatter {

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
	
	public FASTAReformatter(String newFName, String newFPath) {
		
		this.fPath= newFPath; 
		this.fName= newFName;
	}
	
	public FASTAReformatter(String absFilePath) {
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

		fptr = new BufferedReader(new FileReader(fPath+ File.separator+ fName));
			
		String read= null;
		String seq= "";
		sequences= new String[0];
		seqNames= new String[0];
		seqExtNames= new String[0];
		StringBuffer sb= null;
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
			System.out.println("Usage: Renamer pattern(-a for all files in dir) \ne.g. FASTAReformatter /home/ug/root 50");
			System.exit(0);
		}
	
		String pattern= args[0];
		String path= "./";
		
		int wrap= 50;
		if (args.length> 1)
			wrap= Integer.parseInt(args[1]);
		
		String[] list= new File(path).list();
		int contain= list.length;
		int cnt= 0;
		if (!pattern.equals("-a")) {
			Vector v= new Vector(list.length);
			for (int i = 0; i < list.length; i++) {
				if (list[i].indexOf(pattern)>= 0) {
					v.add(list[i]);
					++cnt;
				}
			}
			list= new String[v.size()];
			for (int i = 0; i < list.length; i++) 
				list[i]= (String) v.elementAt(i);
		}  
		
		
		for (int i = 0; i < list.length; i++) {
			System.out.println("processing "+list.length+" files.");
			System.out.println("file "+list[i]);
			FASTAReformatter wrapper= new FASTAReformatter(new File(path+list[i]).getAbsolutePath());
			wrapper.setWrapLine(wrap);
			try {
				System.out.print("\treading..");
				wrapper.read();
				System.out.print("done.\n\twriting.."); 
				wrapper.write();
				System.out.println("done.");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		
		System.out.println(cnt+ " files changed, "+ (contain- cnt)+ " files unchanged");
	}
}
