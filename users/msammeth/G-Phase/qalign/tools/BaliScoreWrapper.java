package qalign.tools;

import gphase.ext.DevPipeReaderThread;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.StringTokenizer;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class BaliScoreWrapper {

	public static String bsExecAbsPath= "C:\\Programme\\balibase";
	public static String bsExecName="bali_score.exe";
	public static String refExt= "msf";
	public static String annExt= "ann";
	public static String globExt= "glob";
	public static String coreExt= "core";
	
	protected String filePath= null;
	protected String fileName= null;
	protected String fileExt= null;
	
	protected int refNb= 0;
	
	public static void main(String[] args) {
		
		BaliScoreWrapper bsw= new BaliScoreWrapper(args[0]);
		bsw.runAll();
	}

	public BaliScoreWrapper(String newAbsFileName) {
		setFileBase(newAbsFileName);
	}

	public BaliScoreWrapper(String newAbsFileName, String newAbsPath, String newExecName) {
		this(newAbsFileName);
		this.bsExecAbsPath= newAbsPath;
		this.bsExecName= newExecName;
	}
	
	protected void setFileBase(String newAbsFileName) {
		
		StringTokenizer st= new StringTokenizer(newAbsFileName, File.separator);

//		System.out.println(newAbsFileName);
			// retrieve path
		this.filePath= "";// File.separator;
		int end= st.countTokens();
		String tmpToken= "";
		for (int i= 0; i< (end- 1); ++i) {
			tmpToken= st.nextToken();
			if (tmpToken.toUpperCase().startsWith("REF")) {
				try {
					this.refNb= Integer.parseInt(tmpToken.substring(3));
				} catch (NumberFormatException e) {
					; //nothing
				}
			}
			this.filePath+= tmpToken+ File.separator;
		}
			// kill last File.separator
		this.filePath= this.filePath.substring(0, this.filePath.length()- 1);
//		System.out.println("path: "+filePath);
		
			// retrieve name
		tmpToken= st.nextToken();
		st= new StringTokenizer(tmpToken, ".");
		this.fileName= st.nextToken();
//		System.out.println("name: "+fileName);
		
			// retrieve ext
		this.fileExt= tmpToken.substring(this.fileName.length()+ 1);
//		System.out.println("ext: "+fileExt);
	}

	public int runGlobalComparison() {
		
			// must be array command type on SUN
			// ("sh -c ..." does NOT work !)
		String[] command= {"cmd", "/c", 
			bsExecAbsPath+ File.separator+ bsExecName+ " "
			+ filePath+ File.separator+ fileName+ "."+ refExt+ " "
			+ filePath+ File.separator+ fileName+ "."+ fileExt+ " "
			+ "-v"};
			
		OutputStream out= null;
		try {
			out= new FileOutputStream(filePath+ File.separator+ fileName+ "."+ globExt);
		} catch (IOException e) {
			; //nothing
		}
		Process proc= null;
		try {
			DevPipeReaderThread pipe= new DevPipeReaderThread(null, out);
			proc= Runtime.getRuntime().exec(command);
			pipe.setIn(proc.getInputStream());
			pipe.start();
		} catch (IOException e) {
			System.err.println("nope"); // nothing
		}
		
		int res= -1;
		try {
			res= proc.waitFor();
		} catch (InterruptedException e) {
			System.err.println("nope"); // nothing
		}
		
		return res;
	}
		
	public int runCoreComparison() {
		
			// must be array command type on SUN
			// ("sh -c ..." does NOT work !)
		String[] command= {"cmd", "/c", 
			bsExecAbsPath+ File.separator+ bsExecName+ " "
			+ filePath+ File.separator+ fileName+ "."+ refExt+ " "
			+ filePath+ File.separator+ fileName+ "."+ fileExt+ " "
			+ filePath+ File.separator+ fileName+ "_ref"+ refNb+ "."+ annExt
			+ " -v"};
			
		OutputStream out= null;
		try {
			out= new FileOutputStream(filePath+ File.separator+ fileName+ "."+ coreExt);
		} catch (IOException e) {
			; //nothing
		}
		Process proc= null;
		try {
			DevPipeReaderThread pipe= new DevPipeReaderThread(null, out);
			proc= Runtime.getRuntime().exec(command);
			pipe.setIn(proc.getInputStream());
			pipe.start();
		} catch (IOException e) {
			System.err.println("nope"); // nothing
		}
		
		int res= -1;
		try {
			res= proc.waitFor();
		} catch (InterruptedException e) {
			System.err.println("nope"); // nothing
		}
		
		return res;
	}
	
	public void runAll() {
		
		System.out.print(filePath+ File.separator+ fileName+ "..");
		int tmpRes= -1;
		
		tmpRes= runGlobalComparison();
		if (tmpRes!= 0)
			System.err.println("WARNING: "
				+ filePath+ File.separator+ fileName
				+ " failed in global comparison."
			);

		tmpRes= runCoreComparison();
		if (tmpRes!= 0)
			System.err.println("WARNING: "
				+ filePath+ File.separator+ fileName
				+ " failed in core comparison."
			);
		System.out.println("done.");
	}
}
