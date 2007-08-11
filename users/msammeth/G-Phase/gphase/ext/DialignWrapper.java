package gphase.ext;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;

import qalign.Constants;
import qalign.DevNullReaderThread;
import qalign.OSChecker;
import gphase.algo.AlignmentWrapper;
import qalign.gui.GUIProxy;
import qalign.tools.MSFWrapper;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class DialignWrapper implements AlignmentWrapper{

	public static String EXEC_PATH= null;
	
	protected boolean finishParse= false;
	protected boolean finished= false;
	protected String[] names= null;
	
	protected int diaNb= -1;
	protected int itNb= -1;
	protected boolean checkNb= false;
	
	protected DevParseReaderThread pipe= null, pipe2= null;
	protected String[] sequences= null;
	protected GUIProxy proxy1= null, proxy2= null;
	protected Process proc= null;
	protected String absFileName= null;

		// config
	protected static String dialignExec= OSChecker.isWindows()?"dialign2-2.exe":"dialign2-2";
	static {
		if (OSChecker.isMacOSX())
			dialignExec="dia22osx";
	}
	protected static String dialignPath= System.getProperty("user.dir")+ File.separator+ "extern"; //"C:\\Programme\\dialign\\dialign2_dir";
	
	protected String[] dialignOptions= {"ff"};
	protected static String[] dialignEnvP= {"DIALIGN2_DIR="+dialignPath+File.separator+"dialign2_dir"};
	protected boolean outputMSF= false;
	protected float threshold= 0f;
	protected int maxLen= 40;
	
	/**
	 * Creating a <code>DialignWrapper</code> with setting a new path.<br>
	 * 
	 * <b>WARNING</b>: Remember to set the shell environment variable
	 * <code>DIALIGN2_DIR</code> to the same directory!
	 */
	public DialignWrapper(String newAbsFileName, String[] newNames) {
		
		this.absFileName= newAbsFileName;
		this.names= newNames;
	}


	/**
	 * Creating a <code>DialignWrapper</code> with setting a new path
	 * and a new executable name.<br>
	 * 
	 * <b>WARNING</b>: Remember to set the shell environment variable
	 * <code>DIALIGN2_DIR</code> to the same directory!
	 */
	public DialignWrapper(String newDialignPath, String newDialignExec) {
		
//		this(newDialignPath);
		dialignExec= newDialignExec;
	}
	
	public DialignWrapper() {
	}

	/**
	 * Creating a <code>DialignWrapper</code> with setting a new path
	 * and a new executable name.<br>
	 * 
	 * <b>WARNING</b>: Remember to set the shell environment variable
	 * <code>DIALIGN2_DIR</code> to the same directory!
	 */
	public DialignWrapper(String newDialignPathBase) {
		
//		this(newDialignPath);
		File tst= new File(newDialignPathBase+ File.separator+ dialignExec);
		try {
			tst.getCanonicalPath();
			dialignPath= newDialignPathBase+ File.separator+ dialignPath;
			dialignEnvP[0]= "DIALIGN2_DIR="+ dialignPath;
		} catch (IOException e) {
			;
		}
	}	


	/**
	 * Creating a <code>DialignWrapper</code> with setting a new path,
	 * a new executable name and a new option array.<br>
	 * 
	 * Options must be given without the indication character (e.g.
	 * <code>'-'</code>.
	 * 
	 * <b>WARNING</b>: Remember to set the shell environment variable
	 * <code>DIALIGN2_DIR</code> to the same directory!
	 */
	public DialignWrapper(String newDialignPath, String newDialignExec, String[] newDialignOptions) {
		
		this(newDialignPath, newDialignExec);
		this.dialignOptions= newDialignOptions;
	}


	/**
	 * Runs QAlign in an external process and returns the program exit code.
	 */
	public int runDialign(String fName) {

			absFileName= fName;
			// init command		
			String diaCommand= 
					dialignPath+ File.separator+ 
					dialignExec;
			//	diaCommand= "\""+ diaCommand+ "\"";
			//	diaCommand+= " -ff ";	// frag file essential
			//	if (outputMSF)
			diaCommand+= " -msf ";
			//	diaCommand+= "-thr "+ threshold+ " ";
			//	diaCommand+= "-lmax "+ maxLen+ " ";
			//	diaCommand+= "-stdo ";
			//	diaCommand+= "-pst ";
			diaCommand+= "-n ";
			diaCommand+= fName;
			// System.out.println(diaCommand);
			
			// start dialign
		int res= 1;
		try {
//			System.out.println(diaCommand);
			if (OSChecker.isWindows())
				proc= Runtime.getRuntime().exec(diaCommand, dialignEnvP);
			else {			
				String[] command= {"bash", "-c", diaCommand};	// terminates always, see 1ckaA_ref4
				proc= Runtime.getRuntime().exec(command, dialignEnvP);
			}
			DevNullReaderThread null1= new DevNullReaderThread(proc.getInputStream());
			DevNullReaderThread null2= new DevNullReaderThread(proc.getErrorStream());
			null1.start();
			null2.start();
			
			res= proc.waitFor();
			proc.destroy();	// IOException: Too many open files
			
		} catch (IOException e) {
			
			System.err.println("Problems running DIALIGN: ");
			e.printStackTrace();
			return res;
			
		} catch (InterruptedException e) {
			; // nothing
		}

			// SUN reports 784416 on everythings ok...
		if (System.getProperty("os.name").toLowerCase().startsWith("sunos"))		
			res= 0;	// correct this
		
		return res;
	}
	
	public String getExecPath() {

		if (EXEC_PATH == null) {
			String execPath= Constants.getICON_PATH().substring(0,
				Constants.getICON_PATH().lastIndexOf(File.separator)
			);
			execPath= execPath.substring(0,
				execPath.lastIndexOf(File.separator)
			);		
			execPath+= File.separator+ "extern"+ File.separator+ "dialign2_dir";
			
			EXEC_PATH= execPath;
		}

		return EXEC_PATH;
	}

	public void loadCommands(DevParseReaderThread parser) {
		
		Class[] args= {String.class};
		Method par= null, inc= null;
		try {
			par= this.getClass().getMethod("parseDiagonalNumber", args);
			inc= this.getClass().getMethod("incDiagonalNumber", args);
		} catch (NoSuchMethodException e) {
			; // nothing
		}
		parser.addCommand("current diagonal =", inc);
		parser.addCommand("total number of diagonals:", par);
	}
	
	
	public void loadNames(DevParseReaderThread parser) {
		
		if (names== null)
			return;
			
		Class[] args= {String.class};
		Method as= null, fin= null;
		try {
			as= this.getClass().getMethod("addSequence", args);
			fin= this.getClass().getMethod("finishProcess", args);
		} catch (NoSuchMethodException e) {
			; // nothing
		}

		for (int i= 0; i< names.length; ++i) {
			parser.addCommand(names[i], as);
		}
		parser.addCommand("Sequence tree:", fin);
	}		
	
	public void parseIterationNumber(String update) {
		
		String up2= update.substring(
			(update.indexOf("iteration step")+ new String("iteration step").length()),
			update.lastIndexOf("in multiple alignment")
		).trim();
		
		int iter= 0;
		try  {
			iter= Integer.parseInt(up2);
		} catch (NumberFormatException e) {
			; // nothing
		}
		if (iter!= proxy1.getValue()) {
			proxy1.setValue(iter);
			System.out.println("setting new iteration");
		}
	}
	
	public void finishComputation() {
		
		proxy1.setValue(proxy1.getMaximum());
		proxy2.setValue(proxy2.getMaximum());
		proxy2.setMessage("[100 %]");
		finished= true;
	}
	
	public void finishProcess(String delMe) {
		
		if (proc!= null) {
//			System.out.println("destroying DIALIGN");
			proc.destroy();
			finishParse= true;
		}
		
		
	}

	public void addSequence(String seqToken) {

		if (finishParse)
			return;
			
		if (!finished&& pipe2.isFinished())
			finishComputation();
			
		if (sequences== null) {
			sequences= new String[names.length];
			for (int i= 0; i< sequences.length; ++i)
				sequences[i]= "";
		}
	
		StringBuffer refName;
		for (int i= 0; i< names.length; ++i) {
			if (seqToken.startsWith(names[i].toLowerCase())) {	// ignore case
				int start= 0;
				for (int j= names[i].length(); j< seqToken.length(); ++j)
					if (Character.isDigit(seqToken.charAt(j))) {
						start= j;
						break;
					}
				for (int j= (start+1); j< seqToken.length(); ++j)
					if (!Character.isDigit(seqToken.charAt(j))) {
						start= j;
						break;
					}
						
				StringBuffer sb= new StringBuffer(
					seqToken.substring(
						start,
						seqToken.length()
					).trim().toUpperCase());
				
				for (int j= 0; j< sb.length(); ++j)
					if (Character.isWhitespace(sb.charAt(j)))
						sb.deleteCharAt(j--);
						
				sequences[i]+= sb.toString();
//				System.out.println(seqToken+": "+i+"("+refName+")");
				return;
			}
		}
//		System.out.println(seqToken+": (NULL)");
	}			
	
	public void parseDiagonalNumber(String update) {

		if (!checkNb)
			return;
		
//		System.out.println("checking: "+update);
		String up2= update.substring(
			(update.indexOf("total number of diagonals:")+ new String("total number of diagonals:").length()),
			update.length()
		).trim();
		
		
		try  {
			diaNb= Integer.parseInt(up2);
		} catch (NumberFormatException e) {
			; // nothing
		}
		proxy2.setMaximum(diaNb);
		proxy2.setValue(0); 
		proxy2.setMessage(
			"["+
			((proxy2.getValue()* 100)/ diaNb)+
			" %]"
		);
		 
		proxy1.increase();
		proxy1.setMessage(
			"Iteration step "+
			proxy1.getValue()
		);
	}
	
	public void incDiagonalNumber(String update) {
		
/*		String up2= update.substring(
			(update.indexOf("current diagonal =")+ new String("current diagonal =").length()),
			update.lastIndexOf("in multiple alignment")
		).trim();
		
		
		int dia= 0;
		try  {
			dia= Integer.parseInt(up2);
		} catch (NumberFormatException e) {
			; // nothing
		}
		proxy2.setValue(dia);
*/
		proxy2.increase();
//		System.out.println(proxy2.getValue()+"/"+diaNb);
		
		proxy2.setMessage(
			"["+
			((proxy2.getValue()* 100)/ diaNb)+
			" %]"
		);
		
		if (proxy2.getValue()>= diaNb) 
			checkNb= true;
		else
			checkNb= false;

	}
					

	/**
	 * Returns the outputMSF.
	 * @return boolean
	 */
	public boolean isOutputMSF() {
		return outputMSF;
	}

	/**
	 * Sets the outputMSF.
	 * @param outputMSF The outputMSF to set
	 */
	public void setOutputMSF(boolean outputMSF) {
		this.outputMSF = outputMSF;
	}

	/**
	 * Sets the threshold.
	 * @param threshold The threshold to set
	 */
	public void setThreshold(float threshold) {
		this.threshold= threshold;
	}

	/**
	 * Returns the maxLen.
	 * @return int
	 */
	public int getMaxLen() {
		return maxLen;
	}

	/**
	 * Sets the maxLen.
	 * @param maxLen The maxLen to set
	 */
	public void setMaxLen(int maxLen) {
		this.maxLen= maxLen;
	}

	public String[] getLayout() {
			// get result
		MSFWrapper msfReader= new MSFWrapper(absFileName+ ".ms");
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		new File(msfReader.getAbsFileName()).delete();
			
		return msfReader.getSequences();
	}
	
	/**
	 * @deprecated no score for Dialign
	 */
	public int getScore() {
		return (-1);
	}
	public void setProxies(GUIProxy globalProxy, GUIProxy subProxy) {
		
		this.proxy1= globalProxy;
		this.proxy2= subProxy;
		
		proxy1.setMinimum(0);
		proxy1.setMaximum(2);
		proxy1.setValue(0);
		proxy1.setMessage(
			"Iteration step "+
			(proxy1.getValue()+1)
		);
		
		
		proxy2.setMinimum(0);
		checkNb= true;
		
	}
}
