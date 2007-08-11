package gphase.ext;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Method;

import com.sun.rmi.rmid.ExecPermission;

import qalign.OSChecker;
import gphase.Constants;
import gphase.algo.AlignmentWrapper;
import qalign.algo.CostTable;
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
public class ClustalWrapper implements AlignmentWrapper {

	int score= -1;
	public static String SCORE_MARKER= "Alignment Score";
	public static String EXEC_NAME="cl.exe";
	static {
		if (OSChecker.isWindows())
			EXEC_NAME="cl.exe";
		else if (OSChecker.isMacOSX())
			EXEC_NAME="clw183osx";
		else if (OSChecker.isLinux()|| OSChecker.isSunOS())
			EXEC_NAME="clustalw";
	}
	
	protected String[] names= null;
	protected String fileAbsName= null;
	protected Process proc= null;
	protected DevParseReaderThread pipe= null;
	protected int seqNb= -1;
	protected int costType= -1;
	
	/**
	 * Constructor for ClustalWrapper.
	 */
	public ClustalWrapper(String fName) {
		this.fileAbsName= fName;
	}
	
	public ClustalWrapper(String fName, String[] newNames, int ct) {
		
		this(fName);
		this.names= newNames;
		this.seqNb= newNames.length;
		this.costType= ct;
	}
	
	public static String getCostMatrix(int ct) {
		
		switch (ct) {
			case CostTable.BLOSUM30:
				return "BLOSUM";	//"blosum30mt";
			case CostTable.BLOSUM45:
				return "BLOSUM";	//"blosum45mt";
			case CostTable.BLOSUM62:
				return "BLOSUM";	//"blosum62mt";
			case CostTable.PAM160:
				return "PAM";		//"pam160mt";
			case CostTable.PAM250:
				return "PAM";		//"pam250mt";
			case CostTable.UNITCOST:
				return "ID";		//"idmat";
			case CostTable.DNA:
				return "CLUSTALW";	//"clustalvdnamt";
			case CostTable.RNA:
				return "CLUSTALW";	//"clustalvdnamt";
			case CostTable.DNARNA:
				return "IUB";		//"swgapdnamt";
			case CostTable.GONNET120:
				return "GONNET";	//"gon120mt";
			case CostTable.GONNET250:
				return "GONNET";	//"gon250mt";
		}
		return null;
	}
	
	/**
		 * @since 1.10T Mac OSX: deletes temp files
		 */
		public int start() {
			
			String execPath= System.getProperty("user.dir");	
			execPath+= File.separator+ "extern";
	
				// init command		
			String command= "";
			command+= execPath+ File.separator+ EXEC_NAME;
			command+= " "+ fileAbsName;
			command+= " -output=gcg";
			command+= " -outorder=input";
			if (CostTable.isProteinMatrix(costType))
				command+= " -matrix="+getCostMatrix(costType);
			else
				command+= " -dnamatrix="+getCostMatrix(costType);
	
	//		System.out.println(command);
	
				// start
			int res= 1;
			try {
				pipe= new DevParseReaderThread(null, this); 
				loadCommands(pipe);
	//			DevPipeReaderThread pipe= new DevPipeReaderThread(null, System.out);
				proc= Runtime.getRuntime().exec(command); 
				pipe.setIn(proc.getInputStream());
				pipe.setName("pipe");
				pipe.start();
	//			System.out.println("running clustal");
				res= proc.waitFor();
				if (pipe!= null) {
					try {
						pipe.join();
					} catch (InterruptedException e) {
						; // nothing
					}
				}
				if (proc!= null) {
					proc.destroy();	// IOException: Too many open files
				}
	//			System.out.println("clustal terminated, destroying...");
				
	//			for (int i= 0; ((getLayout()!= null)&& (i< getLayout().length)); ++i)
	//				System.out.println(getLayout()[i]);
				
	
			} catch (IOException e) {
				
				System.err.println("Problems running CLUSTAL: ");
				e.printStackTrace();
				return res;
				
			} catch (InterruptedException e) {
				; // nothing
			}
	
				// try to remove temp file
			try {
				File delme= new File(fileAbsName);
	//			delme.delete();
			} catch (Exception e) {
				; // nothing
			}
	
				// SUN reports 784416 on everythings ok...
			if (System.getProperty("os.name").toLowerCase().startsWith("sunos"))		
				res= 0;	// correct this
			
			return res;
		}

	/**
	 * @since 1.10T Mac OSX: deletes temp files
	 */
	public int startGrid() {
		
		String execPath= System.getProperty("user.dir");	
		execPath+= File.separator+ "extern";

			// init command		
		String command= "";
		command+= execPath+ File.separator+ EXEC_NAME;
		command+= " "+ fileAbsName;
		command+= " -output=gcg";
		command+= " -outorder=input";
		if (CostTable.isProteinMatrix(costType))
			command+= " -matrix="+getCostMatrix(costType);
		else
			command+= " -dnamatrix="+getCostMatrix(costType);

			// write command in shell script
		File script= null;
		try {
			script= File.createTempFile("123",".sh",new File(Constants.getScratchDir()));
			BufferedWriter buffy= new BufferedWriter(new FileWriter(script));
			buffy.write("qsub "+
				" -o "+ Constants.getScratchDir()+
				" -e "+ Constants.getScratchDir()+
				" "+ command);
			buffy.flush();
			buffy.close();
			Process p= Runtime.getRuntime().exec("chmod 755 "+ script.getAbsolutePath());
			p.waitFor();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

//		System.out.println(command);
		

			// start
		int res= 1;
		try {
			proc= Runtime.getRuntime().exec(script.getAbsolutePath()); 
			res= proc.waitFor();
			
		} catch (IOException e) {
			
			System.err.println("Problems running CLUSTAL: ");
			e.printStackTrace();
			return res;
			
		} catch (InterruptedException e) {
			; // nothing
		}

			// try to remove temp file
		try {
			File delme= new File(fileAbsName);
//			delme.delete();
		} catch (Exception e) {
			; // nothing
		}

			// SUN reports 784416 on everythings ok...
		if (System.getProperty("os.name").toLowerCase().startsWith("sunos"))		
			res= 0;	// correct this
		
		return res;
	}
	
	
/**
 * @since v1.10T Mac OSX
 */
public void stop() {
		
		if (proc== null)
			return;

			// quit parser threads and wait for them
/*		try {
			proc.getOutputStream().flush();	// faster than setFinished()
		} catch (IOException e) {
			; 
		}
		if (pipe!= null)
			pipe.setFinished(true);
		if (pipe2!= null)
			pipe2.setFinished(true);
		try {
			if (pipe!= null) {
				pipe.join();
				pipe= null;
			}
			if (pipe2!= null) {
				pipe2.join();
				pipe2= null;
			}
		} catch (InterruptedException e) {
			; // nothing
		}

*/
			// kill ext process
		try {
			proc.getOutputStream().close();	// faster than setFinished()
//			proc.getInputStream().close();	// do not close the other ones:
//			proc.getErrorStream().close();	// pipe2 deadlocks!
		} catch (IOException e) {
			; // nothing
		}		
		proc.destroy();
/*		try {
			Thread.currentThread().sleep(100);
		} catch (InterruptedException e) {
			;//
		}
*/		
		proc= null;	// here, in the end
	}
	
	public String[] getLayout() {
		
		if (proc== null)
			return null;
			
		String msfAbsFile= fileAbsName.substring(0, fileAbsName.lastIndexOf('.'))+ ".msf";
		MSFWrapper msf= new MSFWrapper(msfAbsFile);
		new File(msfAbsFile).deleteOnExit();
		try {
			msf.read();
		} catch (Exception e) {
			return null; // clustal has built shit
		}
		new File(msf.getAbsFileName()).delete();
		
		if (msf.getWrappedSequences()== null)
			return null;
			
		String[] sequences= new String[msf.getWrappedSequences().length];
		StringBuffer tmpSb;
		for (int i= 0; i< sequences.length; ++i) {
			tmpSb= new StringBuffer(msf.getWrappedSequences()[i].getSequence());
			for (int j= 0; j< tmpSb.length(); ++j)
				if ((Constants.IS_GAP(tmpSb.charAt(j)))&& 
					(tmpSb.charAt(j)!= '-'))
					tmpSb.setCharAt(j, '-');
						sequences[i]= 
				sequences[i]= tmpSb.toString();
		}
		
		return sequences;
	}			


	public void loadCommands(DevParseReaderThread parser) {
		
		Class[] args= {String.class};
		Method up2= null;
		try {
			up2= this.getClass().getMethod("update2", args);
		} catch (NoSuchMethodException e) {
			; // nothing
		}

		parser.addCommand(SCORE_MARKER, up2);
		
	}
	
	
	public void update2(String update) {
		
//		System.out.println(update);
		score= Integer.parseInt(
				update.substring(update.indexOf(SCORE_MARKER)+ SCORE_MARKER.length()+ 1).trim());
		
	}
	
	public int getCostType() {
		return costType;
	}
	public void setCostType(int costType) {
		this.costType = costType;
	}
	public int getScore() {
		return score;
	}
	public void setScore(int score) {
		this.score = score;
	}
}
