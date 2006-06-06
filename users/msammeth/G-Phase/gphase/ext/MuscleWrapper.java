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
public class MuscleWrapper implements AlignmentWrapper {

	String outName= null;
	String inName= null;
	
	static String execPath= System.getProperty("user.dir")+ File.separator+ "extern"+ File.separator+ "muscle"+ File.separator+ "muscle3.52";
	static String execName= "muscle";
	
	
	/**
	 * 
	 */
	public MuscleWrapper() {
		super();
	}
	
	public MuscleWrapper(String newInName, String newOutName) {
		
		setInName(newInName);
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
	
	public void execute() {
		
		String cmd= execPath+ File.separator+ execName
			+" -in "+ getInName()
			+" -out "+ getOutName();
		
		cmd+= " -maxmb 4096";		
		cmd+= " -maxiters 1";
		cmd+= " -diags1";
		cmd+= " -sv";

		try {
			// System.out.println(cmd);
			Process proc= Runtime.getRuntime().exec(cmd);
			
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
		MuscleWrapper.execName = execName;
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
		MuscleWrapper.execPath = execPath;
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
	/**
	 * @return Returns the inName.
	 */
	public String getInName() {
		return inName;
	}
	/**
	 * @param inName The inName to set.
	 */
	public void setInName(String inName) {
		this.inName = inName;
	}
}
