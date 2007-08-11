/*
 * Created on Sep 2, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.ext;

import gphase.Constants;
import gphase.algo.AlignmentWrapper;

import java.io.File;
import java.io.IOException;

import qalign.OSChecker;
import qalign.tools.FASTAWrapper;
import qalign.tools.MSFWrapper;

/**
 * 
 * 
 * @author msammeth
 */
public class MLaganWrapper implements AlignmentWrapper {

	final static String ENVP= "LAGAN_DIR";
	
	String[] names= null;
	String[] seqs= null;
	String tree= null;
	String outName= null;
	
	static String execPath= System.getProperty("user.dir")+ File.separator+ "extern"+ File.separator+ "mlagan"+ File.separator+ "lagan12";
	static String execName= "mlagan";
		
	/**
	 * 
	 */
	public MLaganWrapper() {
		super();
	}
	
	public MLaganWrapper(String[] newNames, String[] newSeqs, String newOutName) {
		
		setNames(newNames);
		setSeqs(newSeqs);
		setOutName(newOutName);
	}
	
	public String[] getLayout() {
		
		FASTAWrapper fas= new FASTAWrapper(outName);
		try {
			fas.read();
		} catch (Exception e) {
			return null;
		}
		
		if (fas.getWrappedSequences()== null)
			return null;
			
		String[] sequences= new String[fas.getWrappedSequences().length];
		StringBuffer tmpSb;
		for (int i= 0; i< sequences.length; ++i) {
			tmpSb= new StringBuffer(fas.getWrappedSequences()[i].getSequence());
			for (int j= 0; j< tmpSb.length(); ++j)
				if ((Constants.IS_GAP(tmpSb.charAt(j)))&& 
					(tmpSb.charAt(j)!= '-'))
					tmpSb.setCharAt(j, '-');
						sequences[i]= 
				sequences[i]= tmpSb.toString();
		}
		
		new File(outName).delete();
		return sequences;
	}
	
	public void execute_old() {
		
		String[] fNames= writeOutTemp(names, seqs);
		String[] cmd= new String[fNames.length+ 5];
		int cctr= 0;
		cmd[cctr++]= execPath+ File.separator+ execName;
		for (int i = 0; i < fNames.length; i++) 
			cmd[cctr++]= fNames[i];
		
		cmd[cctr++]= "-tree";
		cmd[cctr++]= "\""+ getTree()+ "\"";
		cmd[cctr++]= "-out";
		cmd[cctr++]= getOutName();
		String envP= ENVP+ "="+ execPath;

		try {
			// System.out.println(cmd);
			Process proc= Runtime.getRuntime().exec(cmd, new String[]{envP});
			
//			DevPipeReaderThread pipe1= new DevPipeReaderThread(proc.getInputStream(), System.out);
//			DevPipeReaderThread pipe2= new DevPipeReaderThread(proc.getErrorStream(), System.err);
			DevNullReaderThread pipe1= new DevNullReaderThread(proc.getErrorStream());
			pipe1.start();
			
			int res= proc.waitFor(); 
			proc.destroy();
		} catch (IOException e) {
			; // nothing
		} catch (InterruptedException e) {
			; // nothing
		}
		
		for (int i = 0; i < fNames.length; i++) { 
			new File(fNames[i]).delete();
		}

}

	public void execute() {
			
			String[] fNames= writeOutTemp(names, seqs);
			String[] cmd= new String[fNames.length+ 5];
			int cctr= 0;
			cmd[cctr++]= execPath+ File.separator+ execName;
			for (int i = 0; i < fNames.length; i++) 
				cmd[cctr++]= fNames[i];
			
			cmd[cctr++]= "-tree";
			cmd[cctr++]= "\""+ getTree()+ "\"";
			cmd[cctr++]= "-out";
			cmd[cctr++]= getOutName();
			String envP= ENVP+ "="+ execPath;
	
			try {
				// System.out.println(cmd);
				Process proc= Runtime.getRuntime().exec(cmd, new String[]{envP});
				
	//			DevPipeReaderThread pipe1= new DevPipeReaderThread(proc.getInputStream(), System.out);
	//			DevPipeReaderThread pipe2= new DevPipeReaderThread(proc.getErrorStream(), System.err);
				DevNullReaderThread pipe1= new DevNullReaderThread(proc.getErrorStream());
				pipe1.start();
				
				int res= proc.waitFor(); 
				proc.destroy();
			} catch (IOException e) {
				; // nothing
			} catch (InterruptedException e) {
				; // nothing
			}
			
			for (int i = 0; i < fNames.length; i++) { 
				new File(fNames[i]).delete();
			}
	
	}

	public void startGrid() {
		
		String[] fNames= writeOutTemp(names, seqs);
		String[] cmd= new String[fNames.length+ 5];
		int cctr= 0;
		cmd[cctr++]= execPath+ File.separator+ execName;
		for (int i = 0; i < fNames.length; i++) 
			cmd[cctr++]= fNames[i];
		
		cmd[cctr++]= "-tree";
		cmd[cctr++]= "\""+ getTree()+ "\"";
		cmd[cctr++]= "-out";
		cmd[cctr++]= getOutName();
		String envP= ENVP+ "="+ execPath;


		try {
			// System.out.println(cmd);
			Process proc= Runtime.getRuntime().exec(cmd, new String[]{envP});
			
//			DevPipeReaderThread pipe1= new DevPipeReaderThread(proc.getInputStream(), System.out);
//			DevPipeReaderThread pipe2= new DevPipeReaderThread(proc.getErrorStream(), System.err);
			DevNullReaderThread pipe1= new DevNullReaderThread(proc.getErrorStream());
			pipe1.start();
			
			int res= proc.waitFor(); 
			proc.destroy();
		} catch (IOException e) {
			; // nothing
		} catch (InterruptedException e) {
			; // nothing
		}
		
		for (int i = 0; i < fNames.length; i++) { 
			new File(fNames[i]).delete();
		}

	}
	
	/**
	 * @return Returns the execName.
	 */
	public static String getExecName() {
		return execName;
	}
	/**
	 * @param execName The execName to set.
	 */
	public static void setExecName(String execName) {
		MLaganWrapper.execName = execName;
	}
	/**
	 * @return Returns the execPath.
	 */
	public static String getExecPath() {
		return execPath;
	}
	/**
	 * @param execPath The execPath to set.
	 */
	public static void setExecPath(String execPath) {
		MLaganWrapper.execPath = execPath;
	}
	/**
	 * @return Returns the names.
	 */
	public String[] getNames() {
		return names;
	}
	/**
	 * @param names The names to set.
	 */
	public void setNames(String[] names) {
		this.names = names;
	}
	/**
	 * @return Returns the seqs.
	 */
	public String[] getSeqs() {
		return seqs;
	}
	/**
	 * @param seqs The seqs to set.
	 */
	public void setSeqs(String[] seqs) {
		this.seqs = seqs;
	}
	/**
	 * @return Returns the tree.
	 */
	public String getTree() {
		return tree;
	}
	/**
	 * @param tree The tree to set.
	 */
	public void setTree(String tree) {
		this.tree = tree;
	}
	
	String[] writeOutTemp(String[] names, String[] seqs) {
		
		String[] tNames= new String[names.length];
		for (int i = 0; i < names.length; i++) {
			File tFile= null;
			try {
				tFile= File.createTempFile("GP_", null);
			} catch (Exception e) {
				return null;
			}
			tNames[i]= tFile.getAbsolutePath();
			tFile.deleteOnExit();

				// write
			FASTAWrapper fasta= new FASTAWrapper(tNames[i]);
			try {	
				fasta.writeFASTA(new String[] {names[i]}, new String[]{seqs[i]});
			} catch (Exception e) {
				return null;
			}
		}
		
		return tNames;
		
	}
	/**
	 * @return Returns the outName.
	 */
	public String getOutName() {
		return outName;
	}
	/**
	 * @param outName The outName to set.
	 */
	public void setOutName(String outName) {
		this.outName = outName;
	}
}
