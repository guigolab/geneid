package gphase.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Method;
import java.util.StringTokenizer;

import qalign.OSChecker;
import qalign.algo.CostTable;
import qalign.gui.GUIProxy;
import gphase.ext.DevParseReaderThread;
import gphase.ext.DevPipeReaderThread;
import qalign.tools.MSFWrapper;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class TCoffeeWrapper  {

	public static String EXEC_NAME= "t_coffee";
	static {
		if (OSChecker.isWindows())
			EXEC_NAME="t_coffee.exe";
		else if (OSChecker.isMacOSX())
			EXEC_NAME="tco137osx";
		else if (OSChecker.isLinux())
			EXEC_NAME="t_coffee";
	}
	public static final String EXEC_PATH= "extern"; //"C:\\Programme\\t-coffee\\t_coffee4win\\t_coffee";
	
	protected int costID= 0;
	protected String fileAbsName= null;
	
	protected Process proc= null;
	protected DevPipeReaderThread pipe= null;
	protected DevPipeReaderThread pipe2= null;
	
	protected String[] names= null;
	protected String[] sequences= null;
	
	protected GUIProxy proxy1= null, proxy2= null;
	protected boolean outputMSF= false;
	protected int seqNb= 0;
	
	public static String getCostMatrix(int ct) {
		
		switch (ct) {
			case CostTable.BLOSUM30:
				return "blosum30mt";
			case CostTable.BLOSUM45:
				return "blosum45mt";
			case CostTable.BLOSUM62:
				return "blosum62mt";
			case CostTable.PAM160:
				return "pam160mt";
			case CostTable.PAM250:
				return "pam250mt";
			case CostTable.UNITCOST:
				return "idmat";
			case CostTable.DNA:
				return "idmat";
			case CostTable.RNA:
				return "idmat";
			case CostTable.DNARNA:
				return "idmat";
			case CostTable.GONNET120:
				return "gon120mt";
			case CostTable.GONNET250:
				return "gon250mt";
		}
		return null;
	}
	
	public TCoffeeWrapper(String fName) {
		this.fileAbsName= fName;
	}
	
	public TCoffeeWrapper(String fName, String[] newNames, int newCostID) {
		
		this(fName);
		this.names= newNames;
		this.seqNb= newNames.length;
		this.costID= newCostID;
	}	
	
	public int runTCoffee(String fName, String old_delme) {
		
			// init command		
		String command= 
				EXEC_PATH+ File.separator+ EXEC_NAME;
		command+= " "+ fName;
		if (outputMSF)
			command+= " -output=msf_aln";
		command+= " -outorder=input";
		command+= " -matrix=blosum62mt";

			// start t-coffee
		int res= 1;
		System.out.println(command);
		try {

			Process diaProc= Runtime.getRuntime().exec(command);
			
			res= diaProc.waitFor();
			diaProc.destroy();	// IOException: Too many open files
			
		} catch (IOException e) {
			
			System.err.println("Problems running T-COFFEE: ");
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
	
	/**
	 * @since 1.10T Mac OSX remove temp file
	 */
	public int start() {
		
		String cmd= "";
		//cmd+= "cmd -c ";	// command shell needed to help tcoffee
		String execPath= "/home/ug/msammeth/no_backup/tcoffee/445/bin";
		cmd+= execPath+ File.separator+ EXEC_NAME;
		cmd+= " "+ fileAbsName;
		//cmd+= " -in=X"+ getCostMatrix(costID);
			// otherwise OSX sets output to current (starting) directory !		
		if (OSChecker.isMacOSX()) 
			cmd+= " -outfile="+ 
				fileAbsName.substring(0, fileAbsName.lastIndexOf("."))+
				".msf";
		cmd+= " -output=gcg";
		cmd+= " -outorder=input";
		
		
//		System.out.println(cmd);

			// changed this in OSX version,
			// but token for WIN version should stay the same..
		StringTokenizer tokenizer= new StringTokenizer(cmd);
		String[] command= new String[tokenizer.countTokens()];
		for (int i= 0; i< command.length; ++i)
			command[i]= tokenizer.nextToken();
//		for (int i= 2; i< cmd.length; ++i)
//			command+= cmd[i]+ " ";			
//		command= "command /c "+ command+ " > .\\extern\\stdout.txt";

			// start t-coffee
		int res= 1;
		try {
//			loadNames(pipe);
				
			
//			DevPipeReaderThread pipe= new DevPipeReaderThread(null, System.out);
//			DevPipeReaderThread pipe2= new DevPipeReaderThread(null, buffy);
			
			proc= Runtime.getRuntime().exec(command);
			pipe= new DevPipeReaderThread(proc.getInputStream(), System.out);			
			pipe.setName("pipe");
			pipe.start();
			pipe2= new DevPipeReaderThread(proc.getErrorStream(), System.err);
			pipe2.setName("pipe2");
			pipe2.start();

//			System.out.println("running tcoffee");
			res= proc.waitFor();
//			System.out.println("tcoffee terminated, destroying...");

			if (proc!= null) {
				proc.destroy();	// IOException: Too many open files
			
//			for (int i= 0; ((getLayout()!= null)&& (i< getLayout().length)); ++i)
//				System.out.println(getLayout()[i]);
			
				
			}
			
		} catch (IOException e) {
			
//			System.err.println("Problems running T-COFFEE: ");
//			e.printStackTrace();
			return res;
			
		} catch (InterruptedException e) {
			; // nothing
		}

			// try to remove temp file
		try {
			File delme= new File(fileAbsName);
			delme.delete();
		} catch (Exception e) {
			; // nothing
		}

			// SUN reports 784416 on everythings ok...
		if (System.getProperty("os.name").toLowerCase().startsWith("sunos"))		
			res= 0;	// correct this
		
		return res;
	}
	
	public String[] getLayout() {
		
		return sequences;
	}
	
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
			

			// inform user
		proxy1.setMessage("aborted.");
		proxy1.setValue(0);
		proxy2.setMessage("");
		proxy2.setValue(0);
		
		proxy1= null;
		proxy2= null;

		sequences= null;
		proc= null;	// here, in the end
	}
	
	public void setProxies(GUIProxy globalProxy, GUIProxy subProxy) {
		
		proxy1= globalProxy;
		proxy2= subProxy;
		
		proxy1.setMinimum(0);
		proxy1.setMaximum(5);
		proxy1.setValue(0);

		proxy2.setValue(0);
	}
	
	public void loadCommands(gphase.ext.DevParseReaderThread parser) {
		
		Class[] args= {String.class};
		Method up1= null, up2= null, p2= null, inc= null;
		try {
			up1= this.getClass().getMethod("update1", args);
			up2= this.getClass().getMethod("update2", args);
			p2= this.getClass().getMethod("parsePercent2", args);
			inc= this.getClass().getMethod("incrementSeqNb", args);
		} catch (NoSuchMethodException e) {
			; // nothing
		}
		parser.addCommand("INPUT FILES", up1);
		parser.addCommand("READ/MAKE LIBRARIES", up1);
		parser.addCommand("COMPUTE PAIRWISE SIMILARITY", up1);
		parser.addCommand("MAKE NEIGHBOR JOINING DENDROGRAM", up1);	// :)
		parser.addCommand("PROGRESSIVE_ALIGNMENT", up1);
		
//		parser.addCommand(" Length= ", inc);		// count
//		parser.addCommand("][", p2);		// make libraries
		parser.addCommand("%]", p2);		// make libraries
		parser.addCommand(": score=", up2);	// pw similarity
//		parser.addCommand("\tGroup ", up2);	// prog alignment
		parser.addCommand("seq)]-->[Score=", up2);	// prog alignment
		
	}
	
	public void loadNames(DevParseReaderThread parser) {
		
		if (names== null)
			return;
			
		Class[] args= {String.class};
		Method as= null;
		try {
			as= this.getClass().getMethod("addSequence", args);
		} catch (NoSuchMethodException e) {
			e.printStackTrace();
			; // nothing
		}

		for (int i= 0; i< names.length; ++i) {
			StringBuffer tmp= new StringBuffer(names[i].toLowerCase());
/*			if (tmp.indexOf(" ")> 0)
				tmp= new StringBuffer(tmp.toString().substring(0, tmp.indexOf(" ")));
			for (int j= 0; j< tmp.length(); ++j)
				if (!Character.isLetterOrDigit(tmp.charAt(j)))
					tmp.setCharAt(j, '_');
*/					
//			System.out.println(tmp);
			parser.addCommand(tmp.toString(), as);
		}
	}
	
	public void addSequence(String seqToken) {
		
		
		if (sequences== null) {
			sequences= new String[names.length];
			for (int i= 0; i< sequences.length; ++i)
				sequences[i]= "";
		}

		StringBuffer refName;
		for (int i= 0; i< names.length; ++i) {
			refName= new StringBuffer(names[i].toLowerCase());
/*			if (refName.indexOf(" ")> 0)
				refName= new StringBuffer(refName.toString().substring(0, refName.indexOf(" ")));
			for (int j= 0; j< refName.length(); ++j)
				if (!Character.isLetterOrDigit(refName.charAt(j)))
					refName.setCharAt(j,'_');
*/					
			if (seqToken.startsWith(refName.toString())) {	// ignore case
				sequences[i]+= 
					seqToken.substring(
						seqToken.indexOf(refName.toString())
							+ refName.length(),
							seqToken.length()
					).trim().toUpperCase();
//				System.out.println(seqToken+": "+i+"("+refName+")");
				return;
			}
		}
//		System.out.println(seqToken+": (NULL)");
	}
		
	
	public void incrementSeqNb(String compabitility_only) {
		++seqNb;
	}
	
	public void update1(String update) {
		
		if (proxy1== null|| proc== null)
			return;
		
			// else
		proxy1.increase();
		proxy1.setMessage(update);
		if (proxy2!= null) {
			
			if (update.indexOf(new String("READ/MAKE LIBRARIES").toLowerCase())> -1)
				proxy2.setMaximum(100);	// percent
			if (update.indexOf(new String("COMPUTE PAIRWISE SIMILARITY").toLowerCase())> -1) 
				proxy2.setMaximum(seqNb* (seqNb- 1)/ 2);	// k* (k-1)/ 2
			if (update.indexOf(new String("PROGRESSIVE_ALIGNMENT").toLowerCase())> -1)
				proxy2.setMaximum(seqNb- 1);				// k-1 progressive joins
			proxy2.setValue(0);
			proxy2.setMessage("");
		}
	}
	
	public void update2(String update) {
		
		if (proxy2== null|| proc== null) 
			return;
		
			// else
		try {
			int i= Integer.parseInt(update);
			proxy2.setValue(i);
			proxy2.setMessage("["+i+ " %]");
		} catch (NumberFormatException e) {
			proxy2.increase();
			proxy2.setMessage("["+ (int)
				(((float) proxy2.getValue()/ (float) proxy2.getMaximum())* 100f)
				+" %]"
			);
			
		}
	}
	
	public void parsePercent2(String update) {
		
		String up2= update.substring(
			(update.lastIndexOf("][")+ 2),
			update.lastIndexOf("%]")
		).trim();
		
		update2(up2);
	}		
	
	public void setOutputMSF(boolean outputMSF) {
		this.outputMSF= outputMSF;
	}
	
	/*
	 * too unsafe with stdout-parsing
	 */
	public String[] getResult() {
		
		if (proc== null)
			return null;
			
		String msfAbsFile= fileAbsName.substring(0, fileAbsName.lastIndexOf('.'))+ ".msf";
		MSFWrapper msf= new MSFWrapper(msfAbsFile);
		try {
			msf.read();
		} catch (Exception e) {
			return null; // clustal has built shit
		}
		
		if (msf.getWrappedSequences()== null)
			return null;
			
		String[] sequences= new String[msf.getWrappedSequences().length];
		StringBuffer tmpSb;
		for (int i= 0; i< sequences.length; ++i) {
			sequences[i]= msf.getWrappedSequences()[i].getSequence();
		}
		
		return sequences;
	}

	public void waitForResult() {	

		BufferedWriter buffy= null;
		try {
			buffy= new BufferedWriter(new FileWriter(".\\extern\\t_coffee.log"));
			buffy.write("trying Result\n");
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			;//
		}
		String msfAbsFile= fileAbsName.substring(0, fileAbsName.lastIndexOf('.'))+ ".msf";
		if (!new File(msfAbsFile).exists()) {
			try {
				buffy= new BufferedWriter(new FileWriter(".\\extern\\t_coffee.log", true));
				buffy.write(msfAbsFile+" does not exist, waiting...\n");
				buffy.flush();
				buffy.close();
			} catch (Exception e) {
				;//
			}

			while (!new File(msfAbsFile).exists())
				try {
					Thread.sleep(100);
				} catch (InterruptedException ex) {
					; //
				}
		}
		try {
			buffy= new BufferedWriter(new FileWriter(".\\extern\\t_coffee.log", true));
			buffy.write("found Result\n");
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			;//
		}
	}		
}
