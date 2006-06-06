/*
 * Created on Nov 26, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.algorithm;

import gummiband.GraphPanel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.JFrame;

import org.apache.xpath.axes.PredicatedNodeTest;

import qalign.OSChecker;
import qalign.algo.CancelException;
import qalign.algo.CostTable;
import qalign.algo.dca.DCA;
import qalign.algo.dca.QDivide;
import qalign.algo.dca.extension.DCAClosure;
import qalign.algo.dialign.AliGraphClosure;
import qalign.algo.dialign.Closure;
import qalign.algo.dialign.Constants;
import qalign.algo.dialign.DialignFASTAWrapper;
import qalign.algo.dialign.DialignWrapper;
import qalign.algo.dialign.MultiFrag;
import qalign.algo.dialign.MultiFragExt;
import qalign.algo.msa.MSA;
import qalign.algo.msa.QAlign;
import qalign.algo.msa.extension.ClosureAlignment;
import qalign.algo.msa.extension.HyperClosure;
import qalign.model.MultipleAlignmentModel;
import qalign.tools.FASTAWrapper;
import qalign.tools.MSFWrapper;
import repeat.data.Alignment;
import repeat.data.Repeat;
import repeat.io.BAliParserCentre;
import repeat.io.BAliParserSchema;
import repeat.io.MultiAlignmentWrapper;
import repeat.io.RepeatMSFWrapper;
import repeat.io.RepeatPoundWrapper;
import sun.security.x509.RFC822Name;

/**
 * 
 * 
 * @author micha
 */
public strictfp class RepeatAligner {
	static int approxCut= 0;
	static int heuristicAlign= 0;
	static boolean align= false;
	static double forceFactor= -1d;
	static int FILTER_PRE_ALIGN= 60;
	int prealigned= 0;

	public static void conqerBaliBase(String pathToG6, String startdir, String startFile, String repExt) {
	
		System.out.println("Starting in "+startdir+" with "+startFile);
		System.out.println("approx "+(approxCut== 1)+", heuristic "+(heuristicAlign==1));
		// retrieve sub-dirs		
		File dir= new File(pathToG6);
		String[] subDirs= dir.list();
		Arrays.sort(subDirs);
	
			// init writer for logging
		BufferedWriter writer= null;
		try {
			String fName= "_"+approxCut+"-"+heuristicAlign+"-";
			if (mode==MODE_STRICT)
				fName+="strict";
			else
				fName+= "relaxed";
			writer= new BufferedWriter(new FileWriter(File.createTempFile("RALog-"+startdir+"-"+startFile+"-"+repExt+"_",fName,new File(pathToG6))));	
			writer.write("started "+ new Date(System.currentTimeMillis())+"\n\n");
		} catch (Exception e) {
			; // nothing, again..
		}
		
			// process subdirs	
		boolean skip= true;	
		for (int i= 0; i < subDirs.length; i++) {
			if (subDirs[i].equals(startdir))
				skip= false;
			if (skip) continue;
			System.out.println("\n== ["+subDirs[i]+"] ==");
			try {
				writer.write("\n== ["+subDirs[i]+"] ==\n");
			} catch (IOException e2) {
				e2.printStackTrace();
			}
			File sdir= new File(pathToG6+File.separator+subDirs[i]);
			String[] files= sdir.list();
			Arrays.sort(files);
			boolean skip2= true;
			if (!subDirs[i].equals(startdir)|| startFile.equals("+"))
				skip2= false;
			for (int j= 0; j < files.length; j++) {
				if (files[j].indexOf(startFile)>= 0)
					skip2= false;
				if (skip2) continue;
				int idx= files[j].indexOf("."+repExt);
				if (idx< 0)
					continue;				
				String procFile= files[j].substring(0,idx);
				System.out.println("processing "+procFile+"..");
				try {
					writer.write("processing "+procFile+"..\n");
				} catch (IOException e) {
					; // :)
				}
				try {
					RepeatAligner myRepAli= new RepeatAligner(
						sdir.getAbsolutePath()+File.separator+procFile,
						repExt);
//					System.exit(0);
				} catch (Exception e) {
					e.printStackTrace(); 
					try {
						writer.write(e.toString());
				} catch (IOException e1) {
						; //:)
					}
				}
				try {
					writer.write("done.\n");
					writer.flush();
				} catch (IOException e1) {
					; // :)
				}
				
				//System.exit(0);
			}
		}
		
	}

	protected final static int MODE_STRICT= 0;
	protected final static int MODE_SUPERSTRICT= 1;
	protected final static int MODE_RELAXED= 3;
	private static final boolean DEBUG= false;
	
	private static boolean ddebug= false;
	private static boolean numericWeightCheck= false;
	
	protected String fileBase= null;
	protected String repeatFileExtension= null;
	protected String[] baseNames= null;
	protected String[] baseSequences= null;
	protected String[] concatSeqs= null;
	protected String[] alignedRepeats= null;
	protected Repeat repeats= null;
	protected String[][] tokenizedSequences= null;
	protected boolean[][] tokenizedRepeats= null;
	
	protected MultiFrag[] fragments= null;
	protected Closure clos= null;
	protected MultiFrag[] rest= null;
	
	protected static int mode= MODE_STRICT;
	
	public RepeatAligner(String newFileBase, String newRepeatFileExtension) {
		
		this.fileBase= newFileBase; 
		this.repeatFileExtension= newRepeatFileExtension;
		init();
	}
	

	
	public String[] constrainedOutputSeqs_old(MultiFrag[] theFragments) {
			
				// only consistent fragments
			Vector consVector= new Vector();
			for (int i= 0; i < theFragments.length; i++) {
				if (theFragments[i].isConsistent())
					consVector.add(theFragments[i]);
			}
			theFragments= new MultiFrag[consVector.size()];
			System.out.println("\n--\n includedFragments:");
			for (int i= 0; i < theFragments.length; i++) {
				theFragments[i]= (MultiFrag) consVector.elementAt(i);
				System.out.println(theFragments[i]);
			}
			
				// assemble closure
			int[] lengths= new int[baseSequences.length];
			for (int i= 0; i< lengths.length; ++i)
				lengths[i]= baseSequences[i].length();
			Closure ali= AliGraphClosure.newAligGraphClosure(
				baseSequences.length, lengths, 0, null);
			for (int i= 0; i < theFragments.length; i++) 
				if (theFragments[i].isConsistent()) {
					if (!AliGraphClosure.alignableFragment(ali, theFragments[i]))
						System.err.println("inconsistency assertion");
					AliGraphClosure.addAlignedFragment(ali, theFragments[i]);
					AliGraphClosure.computeClosure(ali);
				}
		
				// construct solution 
			String[] result= new String[baseNames.length];
			int[] pos= new int[baseNames.length];
			int[] ins= new int[baseNames.length];
			for (int j= 0; j < result.length; j++) {
				result[j]= "";
				pos[j]= 0;
			}
			boolean run= true;
			while(run) {
			
					// break condition
				boolean breaks= true;
				for (int i= 0; i < ins.length; i++) 
					breaks&= (pos[i]>= baseSequences[i].length());
				if (breaks) {
					run= false;
					continue;
				}

				boolean[] skip= new boolean[pos.length];
				for (int i = 0; i < skip.length; i++) 
					skip[i]= false;
				ins= constrainedAliComputeIns(ali, pos, skip);
				
			
				int min= Integer.MAX_VALUE;					// find progress step
				for (int i= 0; i < ins.length; i++) 
					if (ins[i]< min)
						min= ins[i];
				if (min== 0)
					min= 1;
				if (min< 0) {			// insert gaps
					int[] oldIns= new int[ins.length];
					for (int i = 0; i < oldIns.length; i++) 
						oldIns[i]= ins[i];
					while (min< 0) {	// solve inconsistency
														// clone old array, init skiptable
						skip= new boolean[ins.length];
						for (int i = 0; i < oldIns.length; i++) {
							oldIns[i]= ins[i];
							if ((ins[i]< 0)&& (pos[i]< baseSequences[i].length()))
								skip[i]= false;
							else
								skip[i]= true;
						}
						
						ins= constrainedAliComputeIns(ali, pos, skip);	// compute new ins
						
						min= Integer.MAX_VALUE;	// new min
						for (int i = 0; i < ins.length; i++) {
							if (skip[i])
								continue;
							min= Math.min(min, ins[i]);
						}
					}
					
					
					for (int i= 0; i < result.length; i++) { // old
						if (oldIns[i]< 0) { 
							if (pos[i]>= baseSequences[i].length())
								result[i]+= "-";
							else {
								boolean aligned= false;
								for (int j = 0; j < oldIns.length; j++) {
									if (AliGraphClosure.predFrontier(ali,i,pos[i]+1,j)==
										AliGraphClosure.succFrontier(ali,i,pos[i]+1,j)) {
										if (i== j)
											continue;
										aligned= true;
										break;
									}
									
								}
								if (aligned)
									result[i]+= Character.toUpperCase(baseSequences[i].charAt(pos[i]++));
								else
									result[i]+= Character.toLowerCase(baseSequences[i].charAt(pos[i]++));
							}
						} else
							result[i]+= "-";
					}
				} else {
					for (int i= 0; i < result.length; i++) {
					
					
						if (min== 1) {			// one char (aligned or unaligned)
							boolean aligned= false;
							for (int j= 0; j < result.length; j++) 
								if (AliGraphClosure.predFrontier(ali,i,pos[i]+1,j)==
									AliGraphClosure.succFrontier(ali,i,pos[i]+1,j)) {
									if (i== j)
										continue;
									aligned= true;
									break;
								}
							char c= Character.toLowerCase(baseSequences[i].charAt(pos[i]));
							if (aligned)
								c= Character.toUpperCase(c);
						
							result[i]+= c;
							pos[i]++;
							continue;
						}
					
													// more chars
						result[i]+= 
							baseSequences[i].substring(pos[i], pos[i]+ min).toLowerCase();
						pos[i]+= min;
					}
				}
			
			
			}
		
			return result;
		}
		
		
	public String[] constrainedOutputSeqs(MultiFrag[] theFragments) {
			
				// only consistent fragments
			Vector consVector= new Vector();
			for (int i= 0; i < theFragments.length; i++) {
				if (theFragments[i].isConsistent())
					consVector.add(theFragments[i]);
			}
			theFragments= new MultiFrag[consVector.size()];
//			System.out.println("\n--\n includedFragments:");
			for (int i= 0; i < theFragments.length; i++) {
				theFragments[i]= (MultiFrag) consVector.elementAt(i);
//				System.out.println(theFragments[i]);
			}
			
				// assemble closure
			int[] lengths= new int[baseSequences.length];
			for (int i= 0; i< lengths.length; ++i)
				lengths[i]= baseSequences[i].length();
			Closure ali= AliGraphClosure.newAligGraphClosure(
				baseSequences.length, lengths, 0, null);
			for (int i= 0; i < theFragments.length; i++) 
				if (theFragments[i].isConsistent()) {
					if (!AliGraphClosure.alignableFragment(ali, theFragments[i]))
						System.err.println("inconsistency assertion");
					AliGraphClosure.addAlignedFragment(ali, theFragments[i]);
					AliGraphClosure.computeClosure(ali);
				}
		
				// construct solution 
			String[] result= new String[baseNames.length];
			int[] pos= new int[baseNames.length];
			int[] ins= new int[baseNames.length];
			for (int j= 0; j < result.length; j++) {
				result[j]= "";
				pos[j]= 0;
			}
			boolean run= true;
			while(run) {
			
					// break condition
				boolean breaks= true;
				for (int i= 0; i < ins.length; i++) 
					breaks&= (pos[i]>= baseSequences[i].length());
				if (breaks) {
					run= false;
					continue;
				}

				boolean[] consume= new boolean[pos.length];	// init consume
				for (int i = 0; i < consume.length; i++) 
					consume[i]= true;
				
				int state= consume(ali,pos,0,1);	
				int ref= 0;
				if (state== 2) {
					ref= 1;							// reference is always consumed
					consume[1]= true;
				} else 
					consume[0]= true;
				for (int i = 0; i < consume.length; i++) {
					if(i== ref)
						continue;
					int lb= AliGraphClosure.predFrontier(ali, ref, pos[ref]+1, i)- 1;
					int ub= AliGraphClosure.succFrontier(ali, ref, pos[ref]+1, i)- 1;
					
					state= consume(ali, pos, ref, i);
					
					if (state== 3) { 
						consume[i]= true;	// add. consuming sequence
					} else if(state== 1) 
						consume[i]= false;	// no consume
					else {
						consume[ref]= false;		// else: shit, have to use new ref
						for (int j = 1; j < i; j++) {
							if(j== ref)
								continue;
							
							if (consume(ali,pos,i,j)== 3)
								consume[j]= true;
							else if (consume(ali,pos,i,j)== 1)
								consume[j]= false;
							else
								System.err.println("error in strategy");	// assert consume(ali,pos,i,j)!= 2;
						}
						ref= i;
					}
				}
				
				for (int i = 0; i < consume.length; i++) {		// cross-check for all tupels whether alignable
					if (i== ref|| !consume[i])
						continue;
					for (int j = 0; j < consume.length; j++) {
						if (j== ref|| !consume[j]|| i== j)
							continue;
						
						state= consume(ali, pos, i, j);
						
						if (state== 1)
							consume[j]= false;
						else if (state== 2) { 
							consume[i]= false;
							break;
						}
						// 3: both stay true
					}
				}
				
				for (int i = 0; i < consume.length; i++) {
					
					boolean aligned= false;
					for (int j= 0; j < result.length; j++) 
						if (AliGraphClosure.predFrontier(ali,i,pos[i]+1,j)==
							AliGraphClosure.succFrontier(ali,i,pos[i]+1,j)) {
							if (i== j)
								continue;
							aligned= true;
							break;
						}
					
					char c= '-';
					if (consume[i]) {
						try {
						 c= Character.toLowerCase(baseSequences[i].charAt(pos[i]++));
						} catch (StringIndexOutOfBoundsException e) {
							pos[i]--;	// not nothing :)
						}
						if (aligned)
							c= Character.toUpperCase(c);
						result[i]+= c;
					} else
						result[i]+= c;
						
				}
			}
					
			return result;
		}
		
		/**
		 * 1: reference (x) moves in order to match lb]-pos[y]-[ub 
		 * 2: y has to move to match (lb,ub) on x
		 * 3: both (x and y) can move
		 * 
		 * @param ali
		 * @param pos
		 * @param x
		 * @param y
		 * @return
		 */
		int consume(Closure ali, int[] pos, int x, int y) {
			
			int lb= 0;
			int ub= 0;
			try {
				lb= AliGraphClosure.predFrontier(ali, x, pos[x]+ 1, y)- 1;
				ub= AliGraphClosure.succFrontier(ali, x, pos[x]+ 1, y)- 1;				
			} catch (ArrayIndexOutOfBoundsException e) {
				return 1;
			}
			
			if (ub== lb) {	// match
				if(ub== pos[y])
					return 3;
				else if (pos[y]< lb)
					return 2;
				else
					return 1;
			}
			
			if (pos[y]>= ub)
				return 1;
			else if (pos[y]<= lb)
				return 2;
			else
				return 3;  
		}
		
		int[] constrainedAliComputeIns(Closure ali, int[] pos, boolean[] skip) {
			
			int[] ins= new int[pos.length];
			for (int i= 0; i < ins.length; i++) 
				ins[i]= baseSequences[i].length()- pos[i];
			
			
			for (int i = 0; i < pos.length; i++) {
				if (skip[i]|| ins[i]< 0)
					continue;
				for (int j= (i+1); j < ins.length; j++) {
					if (skip[j]|| ins[j]< 0)
						continue;
					// debug
					int i2j_low= AliGraphClosure.predFrontier(ali,i,(pos[i]+1),j);
					int i2j_hi= AliGraphClosure.succFrontier(ali,i,(pos[i]+1),j);
					int j2i_low= AliGraphClosure.predFrontier(ali,j,(pos[j]+1),i);
					int j2i_hi= AliGraphClosure.succFrontier(ali,j,(pos[j]+1),i);
					if (AliGraphClosure.predFrontier(ali,i,(pos[i]+1),j)== 
							AliGraphClosure.succFrontier(ali,i,(pos[i]+1),j)) { 
						
						ins[i]= Math.min(
							AliGraphClosure.succFrontier(ali,i,(pos[i]+1),j)- 1- pos[j],	// frontier incl
							ins[i]
						);
					} else

						ins[i]= Math.min(
							AliGraphClosure.succFrontier(ali,i,(pos[i]+1),j)- 2- pos[j],	// excluded
							ins[i]
						);
				
					if (AliGraphClosure.predFrontier(ali,j,(pos[j]+1),i)== 
							AliGraphClosure.succFrontier(ali,j,(pos[j]+1),i)) 
						ins[j]= Math.min(
							AliGraphClosure.succFrontier(ali,j,(pos[j]+1),i)- 1- pos[i],	// incl
							ins[j]
						);
					else
						ins[j]= Math.min(
							AliGraphClosure.succFrontier(ali,j,(pos[j]+1),i)- 2- pos[i],	// excl
							ins[j]
						);
				}
			}
			
			return ins;
			
		}
		
		
		/**
		 * new bali init.
		 * only dependent on alignment of repeats 
		 */
		protected void init() {
	
			long myStartTime= System.currentTimeMillis();
			
				// read in unaligned sequences
			System.out.print("reading raw sequences...");			
			readSeqFasta();
			System.err.flush();
			System.out.println("done.");
	//		for (int i= 0; i < baseNames.length; i++) 
	//			System.out.println(baseNames[i]+ "\t"+ baseSequences[i]);
			
				// read in repeats
			System.out.print("reading repeats...");
			readRepeats();			
			System.err.flush();
			System.out.println("done.");
	
				// check for all repeats aligned
			System.out.print("checking for repeat alignments...");
			checkRepeatAlignment();
			System.err.flush();
			System.out.println("done.");
			
				// tokenize Sequences
			System.out.print("tokenizing sequences...");
			tokenizeSequences();
			System.err.flush();
			System.out.println("done.");
	
			long time0= System.currentTimeMillis();
			
				// retrieve all 3 kinds of fragments
			MultiFrag[] ivFragments= null;
			ivFragments= getIVFragments();
			int nextFragNb= 0;
			System.out.println("found "+ivFragments.length+" intervening frags.");
			if (ivFragments!= null)
				nextFragNb= ivFragments.length+ 1;	// 1-based
			long time1= System.currentTimeMillis();
			
			MultiFrag[] repFragments= getRepeatFragments(nextFragNb, repeats);
			System.out.println("found "+repFragments.length+" repeat frags.");
			long time2= System.currentTimeMillis();
	
			MultiFrag[] rlxFragments= new MultiFrag[0]; 
			if (mode== MODE_RELAXED) {
				nextFragNb= ivFragments.length+ repFragments.length+ 1;
				rlxFragments= getRelaxedFragments(nextFragNb);
				System.out.println("found "+rlxFragments.length+" relaxed frags.");
			}
			long time3= System.currentTimeMillis();
	
	
				// unite
			System.out.print("filtering pw optimal sets..");
// here: strategy
			if (mode== MODE_STRICT)
				repFragments= filterPW(repFragments, true, true);			// find pw consistent sets
			System.out.println("done.");

			fragments= new MultiFrag[ivFragments.length+ repFragments.length+ rlxFragments.length];
			for (int i= 0; i < ivFragments.length; i++) {
				fragments[i]= ivFragments[i]; 
				fragments[i].setConsistent(true);
			}
			for (int i= 0; i < repFragments.length; i++) {
				fragments[i+ ivFragments.length]= repFragments[i];
//				fragments[i+ ivFragments.length].setConsistent(true);
			}
			if (mode== MODE_RELAXED) {
				for (int i= 0; i < rlxFragments.length; i++) {
					fragments[i+ ivFragments.length+ repFragments.length]= rlxFragments[i];
					fragments[i+ ivFragments.length+ repFragments.length].setConsistent(true);
				}
			}
	
			if (mode== MODE_RELAXED)
				fragments= filterPW(fragments, false, true);	// strategy
			if (mode== MODE_RELAXED)
				calcOLW(fragments);								// calc OLW
			else {
				calcOLW(ivFragments);	// cannot overlap
				calcOLW(repFragments);
			}
	//		MultiFrag[] clonedFrags= new MultiFrag[fragments.length];
	//		for (int i= 0; i < clonedFrags.length; i++) 
	//			clonedFrags[i]= (MultiFrag) fragments[i].clone();
	//		if (mode== MODE_RELAXED)
	//			calcOLW_old(fragments);								// calc OLW
	//		else {
	//			calcOLW_old(ivFragments);
	//			calcOLW_old(repFragments);
	//		}
	//		for (int i= 0; i < clonedFrags.length; i++) {
	//			if (clonedFrags[i].getNumber()!= fragments[i].getNumber())
	//				System.err.println("numbers differ");
	//			if (clonedFrags[i].getOverlapWeight()!= fragments[i].getOverlapWeight())
	//				System.err.println("different olws: "+clonedFrags[i].getOverlapWeight()
	//								+"<>"+fragments[i].getOverlapWeight());
	//		}
			
							
				// sort according to criterion		
			sortFragmentsOLW(fragments);
			sortFragmentsWGT(rest);
	
				// assemble closure		
	/*		System.out.println("assembling closure..");
			int added= initClosure();
			System.out.println("--- added "+added +" fragments ---");
	*/		
	
				// split up combined fragments		
			Vector fragsExpl= new Vector(fragments.length, 10);
			int michi_tmp= 0;
			for (int i= 0; i < fragments.length; i++) {
				if ((fragments[i] instanceof MultiFragExt)
						&& (((MultiFragExt) fragments[i]).hasJumpPoints())) {
					MultiFrag[] explFrags= ((MultiFragExt) fragments[i]).getFragments();
					for (int j= 0; j < explFrags.length; j++) 
						fragsExpl.add(explFrags[j]);
				} else
					fragsExpl.add(fragments[i]);
			}
			fragments= new MultiFrag[fragsExpl.size()];		
			for (int i= 0; i < fragsExpl.size(); i++) 
				fragments[i]= (MultiFrag) fragsExpl.elementAt(i);
	
				// assemble closure -2-		
			System.out.println("assembling closure..");
			int added= initClosure();
			added+= initClosure(rest);
			System.out.println("\t--- added "+added +" fragments ---");
	
	//		for (int i= 0; i < fragments.length; i++)
	//			if (fragments[i].isConsistent())
	//				System.out.println(i+ ":\t"+ fragments[i]); 
	//		System.out.println("--- --- ---");
			System.out.println("done.");					
			
			
				
				// align rest
			long mainStartTime= System.currentTimeMillis();
			String[] alignedSeqs= null;
			if (align)
				alignedSeqs= constrainedAlignSeqs();
			else
				alignedSeqs= constrainedOutputSeqs(fragments);
//			System.out.println();
//			for (int i = 0; i < alignedSeqs.length; i++) 
//				System.out.println(alignedSeqs[i]);
			
			long myEndTime= System.currentTimeMillis();
		
				// create output file
			File outFile= null;
			try {
				String fName= fileBase+"_"+repeatFileExtension+approxCut+"-"+heuristicAlign+"-";
				if (mode==MODE_STRICT)
					fName+="strict";
				else
					fName+= "relaxed";
				fName+= ".msf";
				outFile= new File(fName);
				RepeatMSFWrapper wrapper= new RepeatMSFWrapper(outFile.getAbsolutePath());
				wrapper.setHighlight(false);
				wrapper.setSeqNames(baseNames);				
				wrapper.setRepeatFragments(fragments);
				wrapper.setSequences(alignedSeqs);
				wrapper.write();
				System.out.println("wrote to "+outFile.getAbsolutePath());
				Date date= new Date(System.currentTimeMillis());
				System.out.println(date.getDay()+" ["+date.getHours()+":"+date.getSeconds()+"]");
				BufferedWriter buffy= new BufferedWriter(new FileWriter(outFile.getAbsolutePath(), true));
				buffy.write("\n\n--\n\nrepeat alignment mode= ");
				switch (mode) {
					case MODE_STRICT: buffy.write("strict (normal) mode"); break;
					case MODE_SUPERSTRICT: buffy.write("SuperStrict mode"); break;
					case MODE_RELAXED: buffy.write("relaxed mode"); break;
				}
				buffy.write("DCA cut approx "+(approxCut== 1)+", msa opt "+(heuristicAlign== 0)+"\n");
				buffy.write("prealigned: "+prealigned+"fragments\n");
				buffy.write("\n\n--\ntotal time: "+(myEndTime- myStartTime)/1000+ " [sec]\n\n");
				buffy.write("\tInput, token, pre-align: "+(time0- mainStartTime)/1000+ " [sec]\n\n");
				buffy.write("\tget interven. fragments: "+(time1- time0)/1000+ " [sec]\n");
				buffy.write("\textract repeat fragments: "+(time2- time1)/1000+ " [sec]\n");
				buffy.write("\tget relaxed fragments: "+(time3- time2)/1000+ " [sec]\n\n");
				buffy.write("\tconstrained alignment: "+(myEndTime- mainStartTime)/1000+ " [sec]\n");
				buffy.flush();
				buffy.close();
			} catch (Exception e) {
				; //nothing
			}
			
			
				// end
			System.out.println("End.");
			
		}
		
	String[] getConstrainedDCA() {
				
		DCAClosure dca= new DCAClosure();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(this.baseNames, this.baseSequences, CostTable.BLOSUM62);
		int added= dca.setFragments(fragments);
		System.out.println("added "+added);
	
		dca.setApproximate(true);			// APPROX HERE !!!
		dca.setMSAParameters(false);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
//		dca.setOutputFormat(QDivide.MSF);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}

		return dca.getAlignmentLayout();
	}

		
	String[] constrainedAlignSeqs() {
			
		DCAClosure dca= new DCAClosure();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(this.baseNames, this.baseSequences, CostTable.BLOSUM62);
		int added= dca.setFragments(fragments);
		System.out.println("added "+added);
	
		dca.setApproximate(approxCut== 1);			// APPROX HERE !!!
		dca.setMSAParameters(heuristicAlign== 1);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
//		dca.setOutputFormat(QDivide.MSF);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
			
		return dca.getAlignmentLayout();
	}
	
	/*
	 * !!! Enhance for overlapping repeat boundaries
	 */
	 // TODO enhance for overlapping repeat boundaries (group3,4?)
	protected void tokenizeSequences() {

		tokenizedSequences= new String[baseSequences.length][];
		String[] tokSeqNames= new String[baseSequences.length];
		tokenizedRepeats= new boolean[baseSequences.length][];
		int last= 0;
		String lastSeqName= "", iS= null;
		int x= -1, y= 0;
		int nameCtr= 0;
			// assumes that repeats are ordered on the sequences
		for (Repeat tmp= repeats;tmp!= null;tmp= tmp.getNext()) {
			
				// position
			if (!tmp.getSeqName().equals(lastSeqName)) {

				tokSeqNames[nameCtr++]= tmp.getSeqName();
				
					// last seq
				if (x!= (-1)) { 	
					iS= getSequence(lastSeqName).substring(last);
					if (iS.length()> 0) {
						tokenizedRepeats[x][y]= false;
						tokenizedSequences[x][y++]= iS;
					}
				}
				
					
				last= 0;
				x++; y= 0;
				String lName= tmp.getSeqName();
				Repeat cRep= tmp;
				int count= 0;
				int lcounted= 0;
				Repeat lRepeat= null;
				while ((cRep!= null)&& (cRep.getSeqName().equals(lName))) {
					if (cRep.getStart()> lcounted)	// count intervening
						count++;
					lRepeat= cRep;
					lcounted= cRep.getStart()+ cRep.getLength();
					cRep= cRep.getNext();
					count++;
				}
				if (lcounted< getSequence(lRepeat.getSeqName()).length())
					count++;
				tokenizedSequences[x]= new String[count];
				tokenizedRepeats[x]= new boolean[count];
			}
			
				// intervening before and repeat
			if (last> tmp.getStart()) {
				System.err.println("Warning! Overlapping domain boundaries in sequence " +
					tmp.getSeqName()+ "("+ last+ ","+ tmp.getStart()+ ")");
				iS="";
			} else
				iS= getSequence(tmp.getSeqName()).substring(last, tmp.getStart());
			if (iS.length()> 0) {
				tokenizedRepeats[x][y]= false; 
				tokenizedSequences[x][y++]= iS;
			}
			tokenizedRepeats[x][y]= true;
			String hs= getSequence(tmp.getSeqName());
			int hi= tmp.getStart()+ tmp.getLength();
			if (hi> hs.length()) {
				hi= hs.length();
				tmp.setLength(hi- tmp.getStart());
				System.err.println("repeat too long..corrected.");
			}
			tokenizedSequences[x][y++]= hs.substring(tmp.getStart(), hi);
			lastSeqName= tmp.getSeqName();
			last= tmp.getStart()+ tmp.getLength();
		}
			// last seq
		iS= getSequence(lastSeqName).substring(last);
		if (iS.length()> 0) {
			tokenizedRepeats[x][y]= false;
			tokenizedSequences[x][y++]= iS;
		}
		
			// remap
		boolean[][] oldTokRep= new boolean[tokenizedRepeats.length][];
		String[][] oldTokSeq= new String[tokenizedSequences.length][];
		for (int i = 0; i < tokenizedSequences.length; i++) {
			if (tokenizedSequences[i]== null)
				continue;
			oldTokSeq[i]= new String[tokenizedSequences[i].length];
			oldTokRep[i]= new boolean[tokenizedRepeats[i].length];
			for (int j = 0; j < oldTokSeq[i].length; j++) {
				oldTokSeq[i][j]= tokenizedSequences[i][j];
				oldTokRep[i][j]= tokenizedRepeats[i][j];
			}
			tokenizedSequences[i]= null;
		}
		for (int i = 0; i < tokSeqNames.length; i++) {
			int j;
			for (j = 0; j < baseNames.length; j++) 
				if (baseNames[j].equalsIgnoreCase(tokSeqNames[i]))
					break;
			if (j>= baseNames.length) {
				if (tokSeqNames[i]!= null)
					System.err.println("remap error for tokenized sequences");
				continue;
			}
			tokenizedRepeats[j]= oldTokRep[i];	// else
			tokenizedSequences[j]= oldTokSeq[i];
		}
		
			// check for seqs without repeats
		for (int i = 0; i < tokenizedSequences.length; i++) {
			if (tokenizedSequences[i]== null) {
				tokenizedSequences[i]= new String[] {baseSequences[i]};
				tokenizedRepeats[i]= new boolean[] {false};
			}
		}

			// output tokenized sequences		
//		for (int i= 0; i< tokenizedSequences.length; ++i) {
//			for (int j= 0; j< tokenizedSequences[i].length; ++j)
//				if (tokenizedRepeats[i][j])
//					System.out.print(":"+ tokenizedSequences[i][j].toLowerCase()+ ":\t");
//				else
//					System.out.print(":"+ tokenizedSequences[i][j]+ ":\t");
//			System.out.println();
//		}
		
	}
	
	protected void checkRepeatAlignment() {
		
			// filter unaligned repeats
		boolean aligned= true;
		int len= -1;
		Repeat tmpRepeat= repeats;
		while(aligned&& tmpRepeat!= null) {
			if (tmpRepeat.getAlignedSeq()== null|| 			// one unaligned sequence OR length variation forces alignment 
				(len>= 0&& tmpRepeat.getAlignedSeq().length()!= len))
				aligned= false;
			if (tmpRepeat!= null&& len< 0)
				len= tmpRepeat.getAlignedSeq().length(); 
			tmpRepeat= tmpRepeat.getNext();
		}
		
			// all are aligned: check for position correctnes (baseSeqs match)
			// check for aligned substring equals position-substring (Balibase)
		if (aligned) {
			tmpRepeat= repeats;
			while(tmpRepeat!= null) {
				MultiFrag frag= new MultiFrag();
				frag.setSequenceStart(true, tmpRepeat.getStart());
				frag.setLength(tmpRepeat.getLength());
				int x;
				for (x= 0; x < baseNames.length; x++) 
					if (baseNames[x].equals(tmpRepeat.getSeqName()))
						break;
				if (x>= baseNames.length)
					System.err.println("Repeat Seq-Reference not found.");
				
				String unalignedSeq= frag.getSubsequence(baseSequences[x], true);
				String alignedSeq= tmpRepeat.getAlignedSeq();
				char c1= Character.toLowerCase(unalignedSeq.charAt(0));
				char c2= Character.toLowerCase(alignedSeq.charAt(0));
				x= 1;
				while (!Character.isLetter(c2))
					c2= Character.toLowerCase(alignedSeq.charAt(x++));
				char c3= Character.toLowerCase(unalignedSeq.charAt(unalignedSeq.length()- 1));
				char c4= Character.toLowerCase(alignedSeq.charAt(alignedSeq.length()- 1));
				x= alignedSeq.length()- 2;
				while (!Character.isLetter(c4))
					c4= Character.toLowerCase(alignedSeq.charAt(x--));
				
				if (c1!= c2|| c3!= c4)
					System.err.println("Repeat boundaries do not match sequence! "+tmpRepeat+"\n"+
						"\t"+unalignedSeq+"\n"+
						"\t"+alignedSeq
					);
					
				tmpRepeat= tmpRepeat.getNext();
			}
			return;											// all repeats have alignment
		}
		
			// else: there are unaligned Seqs
//		alignRepeatsDiaDCA(alignRepeatsGetAll());
		alignRepeatsConsensus(alignRepeatsCountTypes());
	}
	
	/**
	 * Counts different types of repeats according to attribute <code>type</code> 
	 * of <code>Repeat</code>. However, if there are subtypes captured in 
	 * <code>alignmentID</code>, they are further distinguished.
	 * Example: for trust detected repeats there is
	 * type= REPEAT_TYPE 1
	 * alignment_ID= O12321:REPEAT_TYPE_1
	 * @return
	 */
	Vector alignRepeatsCountTypes() {

			// sort according to subtype or leave unsorted
			// type only for unaligned seqs which have been assigned types
		int countTypes= 0;
		Repeat tmpRepeat= repeats;
		Hashtable types= new Hashtable();
		Vector values;
		while(tmpRepeat!= null) {
			Object k= tmpRepeat.getAlignmentID();
			if (k== null)
				k= tmpRepeat.getType();
			if (k== null)			// also store null
				k= k.toString();
			if (types.get(k)== null) {
				values= new Vector();
				values.add(tmpRepeat);
				types.put(k, values);
			} else  
				((Vector) types.get(k)).add(tmpRepeat);
	
			tmpRepeat= tmpRepeat.getNext();
		}
		Vector cat= new Vector(types.values());		// categories= aliID sorting
		for (int i = 0; i < cat.size(); i++) {		// prune gap-only col
			Vector tmpVec= (Vector) cat.elementAt(i);
			String[] tmpLayout= new String[tmpVec.size()];
			for (int j = 0; j <	 tmpVec.size(); j++) 
				tmpLayout[j]= ((Repeat) tmpVec.elementAt(j)).getAlignedSeq();
			for (int j = 0; j < tmpLayout[0].length(); j++) {
				int k= 0;
				for (k = 0; k < tmpLayout.length; k++) 		// gap-only col
					if (tmpLayout[k].charAt(j)!= '-')
						break;
				if (k>= tmpLayout.length) {					// eliminate col
					for (int index = 0; index < tmpLayout.length; index++) 
						tmpLayout[index]= tmpLayout[index].substring(0,j)
								+ tmpLayout[index].substring(j+1, tmpLayout[index].length());
					//System.err.println("purged gap-col");
					--j;
				}
			}
					// update
			if (tmpLayout[0].length()!= ((Repeat) tmpVec.elementAt(0)).getAlignedSeq().length())
				for (int j = 0; j < tmpLayout.length; j++) 
					((Repeat) tmpVec.elementAt(j)).setAlignedSeq(tmpLayout[j]);
		}
		
		return cat;
	}
	
	Vector alignRepeatsGetAll() {

			// get all repeat strings
		Repeat tmpRepeat= repeats;
		Vector seqVec= new Vector();
		while(tmpRepeat!= null) {
				
			seqVec.add(baseSequences[tmpRepeat.getSeqNb()].substring(
				tmpRepeat.getStart(), tmpRepeat.getStart()+ tmpRepeat.getLength()
			));
			if (tmpRepeat.getAlignedSeq()!= null) {						// check for correct string
				StringBuffer sb= new StringBuffer(tmpRepeat.getAlignedSeq());
				while(sb.indexOf("-")>= 0)
					sb.deleteCharAt(sb.indexOf("-"));
				String s1= ((String) seqVec.elementAt(seqVec.size()- 1)).trim().toLowerCase();
				String s2= sb.toString().trim().toLowerCase();
				if (!s1.equals(s2))
					System.err.println("seq retrieve mismatch ("+tmpRepeat.getAlignmentID()+"):\n"+s1+"\n"+s2);
			}
			tmpRepeat= tmpRepeat.getNext();
		}
		return seqVec;
	}
	
	
	void alignRepeatsDialign(Vector seqVec) {
		
			// align sequences Dialign
		String[] repSeqs= new String[seqVec.size()];		// convert to String[]
		String[] repNames= new String[seqVec.size()];
		for (int i= 0; i < repSeqs.length; i++) {
			repSeqs[i]= (String) seqVec.elementAt(i);
			repNames[i]= new Integer(i).toString();
		}
		
		DialignWrapper diaWrap= new DialignWrapper(System.getProperty("user.dir"));
		String fName= writeOutTemp(repNames, repSeqs);
		diaWrap.setOutputMSF(true);
		diaWrap.runDialign(fName);

		MSFWrapper msfReader= new MSFWrapper(fName+ ".ms");	// get result
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		repSeqs= msfReader.getSequences();
		repNames= msfReader.getSeqNames();
		DialignFASTAWrapper diaFasta= new DialignFASTAWrapper(fName);
		try {
			diaFasta.readFF();
		} catch (Exception e) {
			; // merkt man schon.
		}
		MultiFrag[] frags= diaFasta.getFragments();
		
		int x= Integer.MIN_VALUE;							// check for remapping
		for (int i= 0; i < repNames.length; i++) {
			if (Integer.parseInt(repNames[i])<= x)
				System.err.println("Need remapping!");
			x= Integer.parseInt(repNames[i]);
		}
		
		
			// feed in back aligned repeats
		Repeat tmpRepeat= repeats;
		x= 0;
		while(tmpRepeat!= null) {
			
			tmpRepeat.setAlignedSeq(repSeqs[x++]);
			tmpRepeat= tmpRepeat.getNext();
		}
	}
	
	
	void alignRepeatsDCA(Vector seqVec) {
		
			// align sequences Dialign
		String[] repSeqs= new String[seqVec.size()];		// convert to String[]
		String[] repNames= new String[seqVec.size()];
		for (int i= 0; i < repSeqs.length; i++) {
			repSeqs[i]= (String) seqVec.elementAt(i);
			repNames[i]= new Integer(i).toString();
		}
		
		QDivide dca= new QDivide();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(repNames, repSeqs, CostTable.BLOSUM62);
	
		dca.setApproximate(false);			// APPROX HERE !!!
		dca.setMSAParameters(false);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
		dca.setOutputFormat(QDivide.MSF);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
		repSeqs= dca.getAlignmentLayout();	// get result
		repNames= dca.getAlignmentNames();
		
		int x= Integer.MIN_VALUE;							// check for remapping
		for (int i= 0; i < repNames.length; i++) {
			if (Integer.parseInt(repNames[i].substring(1,repNames[i].length()).trim())<= x)
				System.err.println("Need remapping!");
			x= Integer.parseInt(repNames[i].substring(1,repNames[i].length()).trim());
		}
		
		
			// feed in back aligned repeats
		Repeat tmpRepeat= repeats;
		x= 0;
		while(tmpRepeat!= null) {
			
			tmpRepeat.setAlignedSeq(repSeqs[x++]);
			tmpRepeat= tmpRepeat.getNext();
		}
	}	
	
	void alignRepeatsDiaDCA(Vector seqVec) {
		
			// align sequences Dialign
		String[] repSeqs= new String[seqVec.size()];		// convert to String[]
		String[] repNames= new String[seqVec.size()];
		for (int i= 0; i < repSeqs.length; i++) {
			repSeqs[i]= (String) seqVec.elementAt(i);
			repNames[i]= new Integer(i).toString();
		}
		
		DialignWrapper diaWrap= new DialignWrapper(System.getProperty("user.dir"));
		String fName= writeOutTemp(repNames, repSeqs);
		diaWrap.setOutputMSF(true);
		diaWrap.runDialign(fName);

		MSFWrapper msfReader= new MSFWrapper(fName+ ".ms");	// get result
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		String[] repSeqsDia= msfReader.getSequences();
		String[] repNamesDia= msfReader.getSeqNames();
		DialignFASTAWrapper diaFasta= new DialignFASTAWrapper(fName);
		try {
			diaFasta.readFF();
		} catch (Exception e) {
			; // merkt man schon.
		}
		MultiFrag[] frags= diaFasta.getFragments();
		
		int x= Integer.MIN_VALUE;							// check for remapping
		for (int i= 0; i < repNames.length; i++) {
			if (Integer.parseInt(repNames[i])<= x)
				System.err.println("Need remapping!");
			x= Integer.parseInt(repNames[i]);
		}
		
			// align DCA
		Vector fragsFilter= new Vector();		// filter Frags
		int total= frags.length;
		for (int i= 0; i < frags.length; i++) 
			if (frags[i].getOverlapWeight()> FILTER_PRE_ALIGN)	// OLW>= 5
				fragsFilter.add(frags[i]);
		frags= new MultiFrag[fragsFilter.size()];
		for (int i= 0; i < frags.length; i++) 
			frags[i]= (MultiFrag) fragsFilter.elementAt(i);
		
		DCAClosure dca= new DCAClosure();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(repNames, repSeqs, CostTable.BLOSUM62);
		int added= dca.setFragments(frags);
		System.out.println("added "+added+" of "+ total+ " frags.");
	
		dca.setApproximate(false);			// APPROX HERE !!!
		dca.setMSAParameters(false);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
		dca.setOutputFormat(QDivide.MSF);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
		repSeqs= dca.getAlignmentLayout();	// get result
		repNames= dca.getAlignmentNames();
		
		x= Integer.MIN_VALUE;							// check for remapping
		for (int i= 0; i < repNames.length; i++) {
			if (Integer.parseInt(repNames[i].substring(1,repNames[i].length()).trim())<= x)
				System.err.println("Need remapping!");
			x= Integer.parseInt(repNames[i].substring(1,repNames[i].length()).trim());
		}
		
			// feed in back aligned repeats
		Repeat tmpRepeat= repeats;
		x= 0;
		while(tmpRepeat!= null) {
			
			tmpRepeat.setAlignedSeq(repSeqs[x++]);
			tmpRepeat= tmpRepeat.getNext();
		}
	}

	void alignRepeatsConsensus(Vector cat) {
		
			// collect consensus-seqs
		String[] consensusSeqs= new String[cat.size()];
		String[] consensusAliIDs= new String[cat.size()];
		String[] layout;
		Vector values;
		for (int i = 0; i < consensusSeqs.length; i++) {
			values= (Vector) cat.elementAt(i);
			consensusAliIDs[i]= ((Repeat) values.elementAt(0)).getAlignmentID();
			layout= new String[values.size()];
			int j;
			for (j = 0; j < values.size(); j++) {
				if (((Repeat) values.elementAt(j)).getAlignedSeq()== null)
					break;
				else
					layout[j]= ((Repeat) values.elementAt(j)).getAlignedSeq();
			}
			if (j< values.size()) {
				// TODO: alignType(values);
				System.err.println("FATAL: No Repeat Alignment!");
				--i;
				continue;
			}
			for (int k = 0; k < layout[0].length(); k++) {
				
				for (j = 0; j < layout.length; j++)		// a letter in the column? 
					if (Character.isLetter(layout[j].charAt(k)))
						break;
				if (j>= layout.length) {				// delete gap columns
					for (j = 0; j < layout.length; j++) 
						layout[j]= layout[j].substring(0,k)+
								layout[j].substring(k+1, layout[j].length());
					--k;
				}
			}

				// else: all seqs have alignment, determine consensus
			int[] alphabet= new int[26];
			consensusSeqs[i]= "";
			for (int k = 0; k < layout[0].length(); k++) {
				
				for (int index = 0; index < alphabet.length; index++)	// init
					alphabet[index]= 0;				
				for (int index = 0; index < layout.length; index++) { 	// count across column
					char c= layout[index].charAt(k);
					if (Character.isLetter(c))
						++alphabet[Character.toUpperCase(layout[index].charAt(k))- 65];
				}
				int max= 0;
				char c= ' ';
				for (int index = 0; index < alphabet.length; index++) {	// look for (a) majority
					if (alphabet[index]> max) {
						c= new Character((char) (index+ 65)).charValue();
						max= alphabet[index];
					}
				}
				if (c== ' ') {
					System.err.println("gap-only column!");
					continue;	// pw gaps in trust output filtered earlier
				}
				consensusSeqs[i]+= c;
			}
		}
//		System.out.println("consensus:");
//		for (int i = 0; i < consensusSeqs.length; i++) 
//			System.out.println(consensusSeqs[i]);
//		System.out.println();
		
			
				// submit to alignment, get aligned seq
		QDivide dca= new QDivide();
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(consensusAliIDs, consensusSeqs, CostTable.BLOSUM62);
	
		dca.setApproximate(approxCut== 1);			// APPROX HERE !!!
		dca.setMSAParameters(approxCut== 1);			// APPROX HERE !!!
		dca.setOutput(false);
		dca.setRecursionStop(40);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
		for (int k = 0; k < dca.getAlignmentNames().length; k++) {	// check
			StringBuffer s1= new StringBuffer(dca.getAlignmentNames()[k]);
			StringBuffer s2= new StringBuffer(consensusAliIDs[k]);
			for (int i = 0; i < s1.length(); i++)			// eliminate special chars (replaced in dca) 
				if (!Character.isLetter(s1.charAt(i)))
					s1.deleteCharAt(i--);
			for (int i = 0; i < s2.length(); i++) 
				if (!Character.isLetter(s2.charAt(i)))
					s2.deleteCharAt(i--);
			if (!s1.toString().equalsIgnoreCase(s2.toString()))
				System.err.println("("+s1+"!="+s2+"): align mismatch, need reorder!");
		}
		layout= dca.getAlignmentLayout();
//		System.out.println("alignment:");
//		for (int i = 0; i < layout.length; i++) {
//			System.out.println(layout[i]);
//		}
//		System.out.println();
				
		
				// insert gaps in all groups 
		for (int i = 0; i < layout.length; i++) 
			for (int j = 0; j < layout[i].length(); j++) 
				if (!Character.isLetter(layout[i].charAt(j))) {
					
					Vector catcat= (Vector) cat.elementAt(i);
					for (int k = 0; k < catcat.size(); k++) {
						Repeat r= (Repeat) catcat.elementAt(k);
							r.setAlignedSeq(
								r.getAlignedSeq().substring(0, j)+
								"-"+ r.getAlignedSeq().substring(j, r.getAlignedSeq().length())
							);
						// r.setType("DCA");
						r.setAlignmentID("DCA");
					}
				}
		
				// test-out
//		for (int i = 0; i < cat.size(); i++) {
//			Vector catcat= (Vector) cat.elementAt(i);
//			System.out.println(((Repeat) catcat.elementAt(0)).getAlignmentID()+ ":");
//			for (int j = 0; j < catcat.size(); j++) 
//				System.out.println(((Repeat) catcat.elementAt(j)).getAlignedSeq());
//			System.out.println();
//		}
//		System.out.println("--");	
	}
	
	/**
	 * reads repeats and syncs with sequences (i.e. call AFTER seq read)
	 *
	 */
	protected void readRepeats() {
		
		RepeatPoundWrapper reader= new RepeatPoundWrapper(fileBase+"."+repeatFileExtension);
		Repeat[] reps= reader.getRepeats();
		Arrays.sort(reps, new Repeat.SequencePositionComparator());	// sort according to seqs& positions
																	// needed lateron e.g. for tokenizing loop
			// nest repeats and get pointer on first
		repeats= reps[0];
		if (reps.length> 1)
			reps[0].setNext(reps[1]);
		for (int i= 1; i < (reps.length-1); i++) {
			reps[i].setPrev(reps[i-1]);
			reps[i].setNext(reps[i+1]);
		}
		if (reps.length> 1)
			reps[reps.length-1].setPrev(reps[reps.length- 2]);

			// assign seq#s
		Repeat tmpRepeat= repeats;						
		while(tmpRepeat!= null) {
			for (int i= 0; i < baseNames.length; i++) { 
				if (baseNames[i].equalsIgnoreCase(
						tmpRepeat.getSeqName())) {
					tmpRepeat.setSeqNb(i);
					break;				
				}	
			}
			if (tmpRepeat.getSeqNb()< 0)
				System.err.println("error in seqNb assignment for "+ tmpRepeat+ " ("+ tmpRepeat.getSeqName()+ ")");
			tmpRepeat= tmpRepeat.getNext();
		}
		
		System.currentTimeMillis();

	}
	

	
		/**
	 * original init.
	 * 
	 * see also 2 switches in BAliSchemaParser.decodeTable()
	 *
	 */
	protected void init_bali_old() {
		
			// read in unaligned sequences
		System.out.print("reading raw sequences...");			
		readSeqFasta();
		System.out.println("done.");
//		for (int i= 0; i < baseNames.length; i++) 
//			System.out.println(baseNames[i]+ "\t"+ baseSequences[i]);
		
			// read in schema
		System.out.print("reading in schema...");
		BAliParserSchema schemaParser= new BAliParserSchema(fileBase+ "_schema.html");
		repeats= schemaParser.getFirstRepeat();
						 
		Repeat tmpRepeat= repeats;						// assign seq#s
		while(tmpRepeat!= null) {
			for (int i= 0; i < baseNames.length; i++) 
				if (baseNames[i].equals(tmpRepeat.getSeqName())) {
					tmpRepeat.setSeqNb(i);
					break;				
				}	
			if (tmpRepeat.getSeqNb()< 0)
				System.err.println("error in seqNb assignment for "+ tmpRepeat+ " ("+ tmpRepeat.getSeqName()+ ")");
			tmpRepeat= tmpRepeat.getNext();
		}
		System.err.flush();
		System.out.println("done.");
/*		System.out.println("\n");
		for (Repeat tmp= repeats;tmp!= null;tmp= tmp.getNext())
			System.out.println(tmp.getID()+ ": "+ getSequence(tmp.getSeqName()).substring(tmp.getStart(), tmp.getStart()+ tmp.getLength())); 
*/		
		
			// read in aligned repeats
//		BAliAlignmentParser aliParser= new BAliAlignmentParser(fileBase+ "_centre.html");
//		Alignment ali= aliParser.getAlignment();
//		System.out.println(ali);
		
			// tokenize Sequences
		System.out.print("tokenizing sequences...");
		tokenizedSequences= new String[baseSequences.length][];
		tokenizedRepeats= new boolean[baseSequences.length][];
		int last= 0;
		String lastSeqName= "", iS= null;
		int x= -1, y= 0;
		for (Repeat tmp= repeats;tmp!= null;tmp= tmp.getNext()) {
			
				// position
			if (!tmp.getSeqName().equals(lastSeqName)) {

					// last seq
				if (x!= (-1)) { 	
					iS= getSequence(lastSeqName).substring(last);
					if (iS.length()> 0) {
						tokenizedRepeats[x][y]= false;
						tokenizedSequences[x][y++]= iS;
					}
				}
				
					
				last= 0;
				x++; y= 0;
				String lName= tmp.getSeqName();
				Repeat cRep= tmp;
				int count= 0;
				int lcounted= 0;
				Repeat lRepeat= null;
				while ((cRep!= null)&& (cRep.getSeqName().equals(lName))) {
					if (cRep.getStart()> lcounted)	// count intervening
						count++;
					lRepeat= cRep;
					lcounted= cRep.getStart()+ cRep.getLength();
					cRep= cRep.getNext();
					count++;
				}
				if (lcounted< getSequence(lRepeat.getSeqName()).length())
					count++;
				tokenizedSequences[x]= new String[count];
				tokenizedRepeats[x]= new boolean[count];
			}
			
				// intervening before and repeat
			if (last> tmp.getStart()) {
				System.err.println("Warning! Overlapping domain boundaries in sequence " +
					tmp.getSeqName()+ "("+ last+ ","+ tmp.getStart()+ ")");
				iS="";
			} else
				iS= getSequence(tmp.getSeqName()).substring(last, tmp.getStart());
			if (iS.length()> 0) {
				tokenizedRepeats[x][y]= false; 
				tokenizedSequences[x][y++]= iS;
			}
			tokenizedRepeats[x][y]= true;
			String hs= getSequence(tmp.getSeqName());
			int hi= tmp.getStart()+ tmp.getLength();
			if (hi> hs.length()) {
				hi= hs.length();
				tmp.setLength(hi- tmp.getStart());
				System.err.println("repeat too long..corrected.");
			}
			tokenizedSequences[x][y++]= hs.substring(tmp.getStart(), hi);
			lastSeqName= tmp.getSeqName();
			last= tmp.getStart()+ tmp.getLength();
		}
			// last seq
		iS= getSequence(lastSeqName).substring(last);
		if (iS.length()> 0) {
			tokenizedRepeats[x][y]= false;
			tokenizedSequences[x][y++]= iS;
		}
		System.err.flush();
		System.out.println("done.");

			// output tokenized sequences		
/*		for (int i= 0; i< tokenizedSequences.length; ++i) {
			for (int j= 0; j< tokenizedSequences[i].length; ++j)
				if (tokenizedRepeats[i][j])
					System.out.print(":"+ tokenizedSequences[i][j].toLowerCase()+ ":\t");
				else
					System.out.print(":"+ tokenizedSequences[i][j]+ ":\t");
			System.out.println();
		}
*/		
			// align intervening sequence fragments
/*		String[] layout= alignInterveningSequences();
		for (int i= 0; i < layout.length; i++) {
			System.out.println(layout[i]);
		}
		writeOut(layout);
*/		
			// retrieve both kinds of fragments
		MultiFrag[] ivFragments= getIVFragments();
		int nextFragNb= 0;
		System.out.println("found "+ivFragments.length+" intervening frags.");
		if (ivFragments!= null)
			nextFragNb= ivFragments.length+ 1;	// 1-based
		String fbase= fileBase.substring(0, fileBase.indexOf("_ref6"));
		Vector rVec= new Vector();
		for (int i= 1; i< 2; ++i) {				// read in multiple alignments
			String fibase= fbase+ /*i+*/ "_ref6_centre.html";
			MultiFrag[] rFragments= getRepeatFragments(nextFragNb, fibase);
			if (rFragments== null)
				break;
			for (int j= 0; j < rFragments.length; j++) 
				rVec.add(rFragments[j]);
			nextFragNb+= rFragments.length;
		}
		MultiFrag[] repFragments= new MultiFrag[rVec.size()];
		for (int i= 0; i < repFragments.length; i++) {
			repFragments[i]= (MultiFrag) rVec.elementAt(i);
		}
		System.out.println("found "+repFragments.length+" repeat frags.");
		MultiFrag[] rlxFragments= new MultiFrag[0]; 
		if (mode== MODE_RELAXED) {
			nextFragNb= ivFragments.length+ repFragments.length+ 1;
			rlxFragments= getRelaxedFragments(nextFragNb);
			System.out.println("found "+rlxFragments.length+" relaxed frags.");
		}

			// unite
		fragments= new MultiFrag[ivFragments.length+ repFragments.length+ rlxFragments.length];
		for (int i= 0; i < ivFragments.length; i++)
			fragments[i]= ivFragments[i]; 
		for (int i= 0; i < repFragments.length; i++)
			fragments[i+ ivFragments.length]= repFragments[i];
		if (mode== MODE_RELAXED) {
			for (int i= 0; i < rlxFragments.length; i++)
				fragments[i+ ivFragments.length+ repFragments.length]= rlxFragments[i];
		}
		fragments= filterPW(fragments, false);
		calcOLW(fragments);
			
			// sort according to criterion		
		sortFragmentsOLW(fragments);

			// assemble closure		
/*		System.out.println("assembling closure..");
		int added= initClosure();
		System.out.println("--- added "+added +" fragments ---");
*/		
			// split up combined fragments		
		Vector fragsExpl= new Vector(fragments.length, 10);
		int michi_tmp= 0;
		for (int i= 0; i < fragments.length; i++) {
			if ((fragments[i] instanceof MultiFragExt)
					&& (((MultiFragExt) fragments[i]).hasJumpPoints())) {
				MultiFrag[] explFrags= ((MultiFragExt) fragments[i]).getFragments();
				for (int j= 0; j < explFrags.length; j++) 
					fragsExpl.add(explFrags[j]);
			} else
				fragsExpl.add(fragments[i]);
		}
		fragments= new MultiFrag[fragsExpl.size()];		
		for (int i= 0; i < fragsExpl.size(); i++) 
			fragments[i]= (MultiFrag) fragsExpl.elementAt(i);

			// assemble closure -2-		
		System.out.println("assembling closure..");
		int added= initClosure();
		System.out.println("--- added "+added +" fragments ---");

		for (int i= 0; i < fragments.length; i++)
			if (fragments[i].isConsistent())
				System.out.println(i+ ":\t"+ fragments[i]); 
		System.out.println("--- --- ---");
		System.out.println("done.");					
		
			
			// align rest
		DCAClosure dca= new DCAClosure();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(this.baseNames, this.baseSequences, CostTable.BLOSUM62);
		dca.setFragments(fragments);

		dca.setApproximate(false);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
//		dca.setOutputFormat(QDivide.MSF);
//		dca.dialignTime= (float) ((System.currentTimeMillis()- myStartTime)/ 1000);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
//		myEndTime= System.currentTimeMillis();
	
			// create output file
		File outFile= null;
		try {
			outFile= File.createTempFile("delme", ".msf");
			RepeatMSFWrapper wrapper= new RepeatMSFWrapper(outFile.getAbsolutePath());
			wrapper.setSeqNames(baseNames);				
			wrapper.setSequences(dca.getAlignmentLayout());
			wrapper.setRepeatFragments(fragments);
			wrapper.write();
			System.out.println("wrote to "+outFile.getAbsolutePath());
		} catch (Exception e) {
			; //nothing
		}
		
		
			// end
		System.out.println("End.");
		
	}
	
	/**
	 * init for repeat detection
	 *
	 */
	protected void init_detect() {
		
			// read in unaligned sequences
		System.out.print("reading raw sequences...");			
		readSeqFasta();
		System.out.println("done.");
//		for (int i= 0; i < baseNames.length; i++) 
//			System.out.println(baseNames[i]+ "\t"+ baseSequences[i]);
		
			// read in schema
		System.out.print("reading in repeats...");
		Alignment[] alis= MultiAlignmentWrapper.readAlignments(
		"D:\\Eigene Dateien\\repeats\\data\\jaap\\apolipoprotein\\apo_iso_trust.txt",
//		"D:\\Eigene Dateien\\repeats\\data\\BAliBASE\\ref6\\test2c\\sh3_2\\repeat_detection\\trust\\2c_sh3_2_detect_trust.txt",
			MultiAlignmentWrapper.TYPE_FLAT_FILE
		); 
		Repeat tmpRepeat= null;						// assign seq#s
		for (int i= 0; i < alis.length; i++) {
			
			int seqNb= (-1);
			for (int j= 0; j < baseNames.length; j++)	// find seq nb 
				if (alis[i].getSeqIDs()[0].startsWith(baseNames[j])) {
					seqNb= j;
					break;
				}
			if (seqNb< 0 )
				System.err.println("seq# not found");
			
			for (int j= 0; j < alis[i].getSequences().length; j++) {
				if (tmpRepeat== null) {
					tmpRepeat= new Repeat();
					repeats= tmpRepeat;					
				} else {
					tmpRepeat.setNext(new Repeat());
					tmpRepeat.getNext().setPrev(tmpRepeat);
					tmpRepeat= tmpRepeat.getNext();
				}
				tmpRepeat.setSeqName(baseNames[seqNb]);
				tmpRepeat.setSeqNb(seqNb);
				tmpRepeat.setStart(alis[i].getStartPos()[j]- 1);	// convert to 0-based
				tmpRepeat.setLength(alis[i].getEndPos()[j]- alis[i].getStartPos()[j]+ 1);
				tmpRepeat.setType("A");
			}
		}
		while(tmpRepeat!= null) {
			for (int i= 0; i < baseNames.length; i++) 
				if (baseNames[i].equals(tmpRepeat.getSeqName())) {
					tmpRepeat.setSeqNb(i);
					break;				
				}	
			if (tmpRepeat.getSeqNb()< 0)
				System.err.println("error in seqNb assignment for "+ tmpRepeat+ " ("+ tmpRepeat.getSeqName()+ ")");
			tmpRepeat= tmpRepeat.getNext();
		}
		System.err.flush();
		System.out.println("done.");
/*		System.out.println("\n");
		for (Repeat tmp= repeats;tmp!= null;tmp= tmp.getNext())
			System.out.println(tmp.getID()+ ": "+ getSequence(tmp.getSeqName()).substring(tmp.getStart(), tmp.getStart()+ tmp.getLength())); 
*/		
		
			// read in aligned repeats
//		BAliAlignmentParser aliParser= new BAliAlignmentParser(fileBase+ "_centre.html");
//		Alignment ali= aliParser.getAlignment();
//		System.out.println(ali);
		
			// tokenize Sequences
		System.out.print("tokenizing sequences...");
		tokenizedSequences= new String[baseSequences.length][];
		tokenizedRepeats= new boolean[baseSequences.length][];
		int last= 0;
		String lastSeqName= "", iS= null;
		int x= -1, y= 0;
		Repeat tmp;
		for (tmp= repeats;tmp!= null;tmp= tmp.getNext()) {
			
				// position
			if (!tmp.getSeqName().equals(lastSeqName)) {

					// last seq
				if (tmp!= repeats) {	// first loop only 	
					iS= getSequence(lastSeqName).substring(last);
					if (iS.length()> 0) {
						tokenizedRepeats[tmp.getPrev().getSeqNb()][y]= false;
						tokenizedSequences[tmp.getPrev().getSeqNb()][y++]= iS;
					}
				}
				
				last= 0; 
				int xx= x;
				x= getSeqNo(tmp.getSeqName());	// current seq position
				for(int i= xx+ 1; i< x; ++i) {
					tokenizedRepeats[i]= new boolean[]{false};
					tokenizedSequences[i]= new String[]{baseSequences[i]};
				}
					 
				y= 0;
				String lName= tmp.getSeqName();
				Repeat cRep= tmp;
				int count= 0;
				int lcounted= 0;
				Repeat lRepeat= null;
				while ((cRep!= null)&& (cRep.getSeqName().equals(lName))) {
					if (cRep.getStart()> lcounted)	// count intervening
						count++;
					lRepeat= cRep;
					lcounted= cRep.getStart()+ cRep.getLength();
					cRep= cRep.getNext();
					count++;
				}
				if (lcounted< getSequence(lRepeat.getSeqName()).length())
					count++;
				tokenizedSequences[tmp.getSeqNb()]= new String[count];
				tokenizedRepeats[tmp.getSeqNb()]= new boolean[count];
			}
			
				// intervening before and repeat
			if (last> tmp.getStart()) {
				System.err.println("Warning! Overlapping domain boundaries in sequence " +
					tmp.getSeqName()+ "("+ last+ ","+ tmp.getStart()+ ")");
				iS="";
			} else
				iS= getSequence(tmp.getSeqName()).substring(last, tmp.getStart());
			if (iS.length()> 0) {
				tokenizedRepeats[tmp.getSeqNb()][y]= false; 
				tokenizedSequences[tmp.getSeqNb()][y++]= iS;
			}
			tokenizedRepeats[tmp.getSeqNb()][y]= true;
			String hs= getSequence(tmp.getSeqName());
			int hi= tmp.getStart()+ tmp.getLength();
			if (hi> hs.length()) {
				hi= hs.length();
				tmp.setLength(hi- tmp.getStart());
				System.err.println("repeat too long..corrected.");
			}
			tokenizedSequences[tmp.getSeqNb()][y++]= hs.substring(tmp.getStart(), hi);
			lastSeqName= tmp.getSeqName();
			last= tmp.getStart()+ tmp.getLength();
		}
			// last seq
		iS= getSequence(lastSeqName).substring(last);
		if (iS.length()> 0) {
			tokenizedRepeats[x][y]= false;
			tokenizedSequences[x][y++]= iS;
		}
		System.err.flush();
		System.out.println("done.");
		
			// fill unused
		for (int i= (x+ 1); i < baseSequences.length; i++) {
			tokenizedRepeats[i]= new boolean[]{false};
			tokenizedSequences[i]= new String[]{baseSequences[i]};
		}

			// output tokenized sequences		
		for (int i= 0; i< tokenizedSequences.length; ++i) {
			for (int j= 0; j< tokenizedSequences[i].length; ++j)
				if (tokenizedRepeats[i][j])
					System.out.print(":"+ tokenizedSequences[i][j].toUpperCase()+ ":\t");
				else
					System.out.print(":"+ tokenizedSequences[i][j].toLowerCase()+ ":\t");
			System.out.println();
		}
		
			// align intervening sequence fragments
			// retrieve intervening fragments
		MultiFrag[] ivFragments= getIVFragments();
		int nextFragNb= 0;
		System.out.println("found "+ivFragments.length+" intervening frags.");
		
			// retrieve repeat fragments
		Vector nameVec= new Vector();
		Vector seqVec= new Vector();
		Vector fragVec= new Vector();
		Vector startPosVec= new Vector(), endPosVec= new Vector();
		for (int i= 0; i < alis.length; i++) {		// collect
			MultiFrag[] tmpFrags= alis[i].getMultiFrags();	// retrieve fragments
			for (int j= 0; j < tmpFrags.length; j++) {		// ..and correct for seq#
				tmpFrags[j].setSequenceNo(true,
					tmpFrags[j].getSequenceNo(true)+ nameVec.size());
				tmpFrags[j].setSequenceNo(false,
					tmpFrags[j].getSequenceNo(false)+ nameVec.size());
				fragVec.add(tmpFrags[j]);					// .. collect frags
			}
			for (int j= 0; j < alis[i].getSeqIDs().length; j++) {
				nameVec.add(alis[i].getSeqIDs()[j]);		// collect seqIDs and seqs
				seqVec.add(alis[i].getSequences()[j]);
				startPosVec.add(new Integer(alis[i].getStartPos()[j]));
				endPosVec.add(new Integer(alis[i].getEndPos()[j]));
			}
		}
		String[] preSeqs= new String[seqVec.size()];	// prepare input
		String[] preNames= new String[nameVec.size()];
		for (int i= 0; i < preSeqs.length; i++) {
			
			StringBuffer tmpSeq= new StringBuffer((String) seqVec.elementAt(i));
			int gapIdx= tmpSeq.indexOf("-");
			while (gapIdx!= (-1)) {
				tmpSeq.deleteCharAt(gapIdx);
				gapIdx= tmpSeq.indexOf("-");
			}
			preSeqs[i]= tmpSeq.toString();
			preNames[i]= (String) nameVec.elementAt(i);
		}
		MultiFrag[] preFrags= new MultiFrag[fragVec.size()];
		for (int i= 0; i < fragVec.size(); i++) 
			preFrags[i]= (MultiFrag) fragVec.elementAt(i);
		
			// *** realign fragments new ***
			// (for different alignments are lacking global context)
			// align sequences
//		System.out.print("re-alinging repeat fragments...");
//		DialignWrapper diaWrap= new DialignWrapper();
//		String fName= writeOutTemp(preNames, preSeqs);
//		diaWrap.setOutputMSF(true);
//		diaWrap.runDialign(fName);
//		preFrags= getDialignFrags(fName);	// get dialign results
//		System.err.flush();
//		System.out.println("done.");

	
		DCAClosure dcaPre= new DCAClosure();
		//dca.setOutFileName(diaFasta.getAbsFileName());
		dcaPre.setCostTable(CostTable.BLOSUM62);
		dcaPre.setSequences(preNames, preSeqs, CostTable.BLOSUM62);
		dcaPre.setFragments(preFrags);

		dcaPre.setApproximate(false);			// APPROX HERE !!!
		dcaPre.setOutput(true);
		dcaPre.setRecursionStop(40);
		// dcaPre.setOutputFormat(QDivide.MSF);
		//dca.dialignTime= (float) ((System.currentTimeMillis()- myStartTime)/ 1000);
		try {dcaPre.run();} 
		catch (CancelException e) {
			; // nothing
		}

		Alignment repAli= new Alignment(dcaPre.getAlignmentLayout());
		String[] tmpNames= dcaPre.getAlignmentNames();	// kill ">" tags!
		for (int i= 0; i < tmpNames.length; i++) 
			tmpNames[i]= tmpNames[i].substring(1, tmpNames[i].length());
		repAli.setSeqIDs(tmpNames);
		int[] startPos= new int[startPosVec.size()];
		int[] endPos= new int[endPosVec.size()];
		for (int i= 0; i < endPos.length; i++) {
			startPos[i]= ((Integer) startPosVec.elementAt(i)).intValue();
			endPos[i]= ((Integer) endPosVec.elementAt(i)).intValue();
		}
		repAli.setStartPos(startPos);
		repAli.setEndPos(endPos);
		for (int i= 0; i < repAli.getSeqIDs().length; i++) {
			System.out.println(
				repAli.getSeqIDs()[i]+ "\t"+
				repAli.getSequences()[i]
			);
			
		}
//		System.out.print("extract fragments...");	// extract fragments
		if (ivFragments!= null)
			nextFragNb= ivFragments.length;
		MultiFragExt[] repFragments= extractFragments(repAli, nextFragNb);
		
		if (mode== MODE_SUPERSTRICT)				// cross-check with schema, assign types
			repFragments= filterMultiFragTypes(repFragments);
		repFragments= scoreFragments(repFragments);				// score fragments
		System.err.flush();
		System.out.println("done."); 
		System.out.println("found "+repFragments.length+" repeat frags.");

			// retrieve relaxed fragments
		MultiFrag[] rlxFragments= new MultiFrag[0]; 
		if (mode== MODE_RELAXED) {
			nextFragNb= ivFragments.length+ repFragments.length+ 1;
			rlxFragments= getRelaxedFragments(nextFragNb);
			System.out.println("found "+rlxFragments.length+" relaxed frags.");
		}

			// unite
		fragments= new MultiFrag[ivFragments.length+ repFragments.length+ rlxFragments.length];
		for (int i= 0; i < ivFragments.length; i++)
			fragments[i]= ivFragments[i]; 
		for (int i= 0; i < repFragments.length; i++)
			fragments[i+ ivFragments.length]= repFragments[i];
		if (mode== MODE_RELAXED) {
			for (int i= 0; i < rlxFragments.length; i++)
				fragments[i+ ivFragments.length+ repFragments.length]= rlxFragments[i];
		}
		fragments= filterPW(fragments, false);
		calcOLW(fragments);
			
			// sort according to criterion		
		sortFragmentsOLW(fragments);
//		System.out.println("* top 100 *");
//		for (int i= 0; i < fragments.length; i++) 
//			if (fragments[i].isConsistent()&& fragments[i].getTranslation()== MultiFragExt.TYPE_REPEAT)
//			System.out.println(fragments[i]);
//		System.out.println("***\n");			
	
			// assemble closure
		System.out.println("assembling closure..");
		int added= initClosure();
		System.out.println("done.");					
		System.out.println("--- added "+added +" fragments ---");
		for (int i= 0; i < fragments.length; i++) {
			if (fragments[i].isConsistent())
				System.out.println(fragments[i]+"\n"+
					fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(true)], true)+ "\n"+
					fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(false)], false)
				);
		}
			// split up combined fragments		
		Vector fragsExpl= new Vector(fragments.length, 10);
		for (int i= 0; i < fragments.length; i++) {
			if ((fragments[i] instanceof MultiFragExt)
					&& (((MultiFragExt) fragments[i]).hasJumpPoints())) {
				MultiFrag[] explFrags= ((MultiFragExt) fragments[i]).getFragments();
				for (int j= 0; j < explFrags.length; j++) 
					fragsExpl.add(explFrags[j]);
			} else 
				fragsExpl.add(fragments[i]);
		}
		
		fragments= new MultiFrag[fragsExpl.size()];	
//		System.out.println("splitted frags to be added:");	
		for (int i= 0; i < fragsExpl.size(); i++) { 
			fragments[i]= (MultiFrag) fragsExpl.elementAt(i);
//			if (fragments[i].isConsistent())
//				System.out.println(fragments[i]);
		}

			// align rest
		DCAClosure dca= new DCAClosure();
//		dca.setOutFileName(diaFasta.getAbsFileName());
		dca.setCostTable(CostTable.BLOSUM62);
		dca.setSequences(this.baseNames, this.baseSequences, CostTable.BLOSUM62);
		dca.setFragments(fragments);

		dca.setApproximate(false);			// APPROX HERE !!!
		dca.setOutput(true);
		dca.setRecursionStop(40);
//		dca.setOutputFormat(QDivide.MSF);
//		dca.dialignTime= (float) ((System.currentTimeMillis()- myStartTime)/ 1000);
		try {dca.run();} 
		catch (CancelException e) {
			; // nothing
		}
//		myEndTime= System.currentTimeMillis();
	
			// create output file
		File outFile= null;
		try {
			outFile= File.createTempFile("delme", ".msf");
			RepeatMSFWrapper wrapper= new RepeatMSFWrapper(outFile.getAbsolutePath());
			wrapper.setSeqNames(baseNames);				
			wrapper.setSequences(dca.getAlignmentLayout());
			wrapper.setRepeatFragments(fragments);
			wrapper.write();
			System.out.println("wrote to "+outFile.getAbsolutePath());
		} catch (Exception e) {
			; //nothing
		}
		
		
			// end
		System.out.println("End.");
		
	}

	
	public MultiFrag[] getRepeatFragments(int nextFragNb, Repeat firstRepeat) {
		
			// sort into groups according to type (case of superstrict also subtype)
		Repeat tmpRepeat= firstRepeat;
		Vector groupIDs= new Vector();
		Vector groups= new Vector();
		final char subtypeSep= '#';
		while (tmpRepeat!= null) {
			
			String id= tmpRepeat.getAlignmentID();
			if (mode== MODE_SUPERSTRICT) 					// check also subtypes
				id+= subtypeSep+ tmpRepeat.getType().toString();
			int i;
			for (i= 0; i < groupIDs.size(); i++) 			// pack into existing group
				if (groupIDs.elementAt(i).equals(id)) {
					Vector group= (Vector) groups.elementAt(i);
					group.add(tmpRepeat);
					break;
				}
				
			if (i>= groupIDs.size()) {						// new type / subtype
				groupIDs.add(id);
				Vector newGroup= new Vector();
				newGroup.add(tmpRepeat);
				groups.add(newGroup);
			}
			
			tmpRepeat= tmpRepeat.getNext();
		}
		
		
			// compose aligned repeats and extract fragments
		System.out.print("extract fragments...");
		Alignment ali;
		Vector tmpGroup;
		String[] seqs;
		String[] seqIDs;
		int[] starts, ends;
		MultiFragExt[] frags;
		Vector collectedFrags= new Vector();
		for (int i= 0; i < groups.size(); i++) {
			System.out.print(".");
			tmpGroup= (Vector) groups.elementAt(i);
			seqs= new String[tmpGroup.size()];
			seqIDs= new String[tmpGroup.size()];
			starts= new int[tmpGroup.size()];
			ends= new int[tmpGroup.size()];
			for (int j= 0; j < tmpGroup.size(); j++) {
				tmpRepeat= (Repeat) tmpGroup.elementAt(j);
				seqs[j]= tmpRepeat.getAlignedSeq();
				seqIDs[j]= tmpRepeat.getSeqName();
				starts[j]= tmpRepeat.getStart();
				ends[j]= tmpRepeat.getStart()+ tmpRepeat.getLength()- 1;	// last pos 
			}
			ali= new Alignment(seqs);
			ali.setSeqIDs(seqIDs);
			ali.setStartPos(starts);
			ali.setEndPos(ends);

			frags= extractFragments(ali, nextFragNb, false);		// extract frags
			nextFragNb+= frags.length;
			for (int j= 0; j < frags.length; j++) 
				collectedFrags.add(frags[j]);
		}
		for (int i = 0; i < collectedFrags.size(); i++) {	// remove self-aligned repeats
			MultiFrag tmp= (MultiFrag) collectedFrags.elementAt(i);
			if (tmp.getSequenceNo(true)== tmp.getSequenceNo(false)) 
				collectedFrags.remove(i--);
			
		}
		frags= new MultiFragExt[collectedFrags.size()];
		for (int j= 0; j < frags.length; j++) 
			frags[j]= (MultiFragExt) collectedFrags.elementAt(j);
		

				// extract fragments
		// cross-check with schema, assign types
//		if (mode== MODE_SUPERSTRICT)
//			frags= filterMultiFragTypes(frags);
					
		// score fragments
		frags= scoreFragments(frags);
		System.err.flush();
		System.out.println("done."); 
		
		return frags;
	}
	
	/**
	 * Filters the given fragments according the repeat types they contain. 
	 * All fragments with not equal repeat types are filtered off.
	 * 
	 * @param originalFrags the fragments to be filtered
	 * @return the filtered fragments
	 */
	// TODO Check for replacing MultiFragExt with MultiFrag signature
	protected MultiFragExt[] filterMultiFragTypes(MultiFragExt[] originalFrags) {
		
			// lazily create Vector for deleting single frags
//		System.out.println(originalFrags.length+ " tuples for frags");
		Vector fragVec= new Vector(originalFrags.length);
		for (int i= 0; i< originalFrags.length; ++i) 
			fragVec.add(originalFrags[i]);				

			// check for same Repeat Types
		for (int i= 0; i< fragVec.size(); ++i) {
			 
			MultiFragExt tmpFrag= (MultiFragExt) fragVec.elementAt(i);

			Repeat tmpRepeat1= repeats;
			while (tmpRepeat1!= null) {
				if (tmpRepeat1.getSeqName().equals(baseNames[tmpFrag.getSequenceNo(true)])
					&& (tmpRepeat1.getStart()== tmpFrag.getSequenceStart(true))
					&& (tmpRepeat1.getLength()== tmpFrag.getRealLength(true))
					)
					break;
				tmpRepeat1= tmpRepeat1.getNext();
			}
			if (tmpRepeat1== null) 
				System.err.println("getRepeatFragments() Repeat Type 1 not found:\n"
					+ tmpFrag
				);

			Repeat tmpRepeat2= repeats;
			while (tmpRepeat2!= null) {
				if (tmpRepeat2.getSeqName().equals(baseNames[tmpFrag.getSequenceNo(false)])
					&& (tmpRepeat2.getStart()== tmpFrag.getSequenceStart(false))
					&& (tmpRepeat2.getLength()== tmpFrag.getRealLength(false))
					)
					break;
				tmpRepeat2= tmpRepeat2.getNext();
			}
			if (tmpRepeat2== null) 
				System.err.println("getRepeatFragments() Repeat Type 2 not found:\n"
					+ tmpFrag
				);
				
			if (!tmpRepeat1.getType().equals(tmpRepeat2.getType())) {
				fragVec.remove(i--);
//				System.out.println("deleted Fragment "+ tmpFrag.getNumber()+ " joining incompatible repeats "
//									+ tmpRepeat1.getID()+ " and "+ tmpRepeat2.getID());
			}
		}
		
			// re-build MultiFrag[]
		MultiFragExt[] newFrags= new MultiFragExt[fragVec.size()];
		for (int i= 0; i < newFrags.length; i++) 
			newFrags[i]= (MultiFragExt) fragVec.elementAt(i);
		
//		System.out.println(newFrags.length+ " frags after filtering");
		return newFrags;
	}
	
	
	public MultiFragExt[] scoreFragments(MultiFragExt[] frags) {

		// calc relative weights
		for (int i= 0; i < frags.length; i++) { 
			MultiFrag[] fragComp= frags[i].getFragments();
			float partWgt= 0;
			for (int j= 0; j < fragComp.length; j++) {
				partWgt+= calcRelWgt(fragComp[j]);
			}
			
			// returns same result
			frags[i].setWeight(partWgt);
//			frags[i].setWeight(calcRelWgt_indel(frags[i]));
		}

		if (DEBUG) {
//			sortFragmentsWGT(frags);
			System.out.println("repeat frags");
			System.out.println("============");
			for (int i= 0; i < frags.length; i++) 
				System.out.println(frags[i]);
		}

		return frags;
	}
	
	public static MultiFrag[] sortFragmentsOLW(MultiFrag[] frags) {
		
		Comparator c= new Comparator() {
			/* (non-Javadoc)
			 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
			 */
			public int compare(Object o1, Object o2) {
			
				MultiFrag frag1= (MultiFrag) o1;
				MultiFrag frag2= (MultiFrag) o2;
				
				if (frag1.getOverlapWeight()< frag2.getOverlapWeight())
					return 1;	// reverse order for ascending weights
				if (frag1.getOverlapWeight()== frag2.getOverlapWeight())
					return 0;
				return (-1);	// reverse order for ascending weights
			}
			/* (non-Javadoc)
			 * @see java.util.Comparator#equals(java.lang.Object)
			 */
			public boolean equals(Object obj) {
				
				return (compare(this, obj)== 0);
			}

		};
		
		Arrays.sort(frags, c);
		
		return frags;
	}
	
	/**
	 * Sorts fragments according to ascending numbers.
	 * 
	 * @param frags
	 * @return
	 */
	public static MultiFrag[] sortFragmentsNumber(MultiFrag[] frags) {
	
		Comparator c= new Comparator() {
			/* (non-Javadoc)
			 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
			 */
			public int compare(Object o1, Object o2) {
	
				MultiFrag frag1= (MultiFrag) o1;
				MultiFrag frag2= (MultiFrag) o2;
	
				if (frag1.getNumber()> frag2.getNumber())
					return 1; // ascending numbers
				if (frag1.getNumber() == frag2.getNumber())
					return 0;
				return (-1); // ascending numbers
			}
			/* (non-Javadoc)
			 * @see java.util.Comparator#equals(java.lang.Object)
			 */
			public boolean equals(Object obj) {
	
				return (compare(this, obj) == 0);
			}
	
		};
	
		Arrays.sort(frags, c);
	
		return frags;
	}
	
	/**
	 * Sorts fragments according to sequence numbers, start positions and lengths.
	 * 
	 * @param frags
	 * @return
	 */
	public static MultiFrag[] sortFragmentsSeqNb(MultiFrag[] frags) {
		
		Comparator c= new Comparator() { 
			/* (non-Javadoc)
			 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
			 */
			public int compare(Object o1, Object o2) {
		
				MultiFrag frag1= (MultiFrag) o1;
				MultiFrag frag2= (MultiFrag) o2;
		
				if (frag1.getSequenceNo(true)> frag2.getSequenceNo(true))
					return 1; // ascending
				else if (frag1.getSequenceNo(true)< frag2.getSequenceNo(true))
					return (-1);
					// else.. first seq# is equal
				if (frag1.getSequenceNo(false)> frag2.getSequenceNo(false))
					return 1;
				else if (frag1.getSequenceNo(false)< frag2.getSequenceNo(false))
					return (-1);	
					// both seq# are equal
				if (frag1.getSequenceStart(true)> frag2.getSequenceStart(true))
					return 1;
				else if (frag1.getSequenceStart(true)< frag2.getSequenceStart(true))
					return (-1);
					// start in 1 is equal
				if (frag1.getSequenceStart(false)> frag2.getSequenceStart(false))
					return 1;
				else if (frag1.getSequenceStart(false)< frag2.getSequenceStart(false))
					return (-1);
					// start in 2 also equal
				if (frag1.getLength()> frag2.getLength())
					return 1;
				if (frag1.getLength()< frag2.getLength())
					return (-1);
					// all equal
				return 0;
				
			}
			/* (non-Javadoc)
			 * @see java.util.Comparator#equals(java.lang.Object)
			 */
			public boolean equals(Object obj) {
		
				return (compare(this, obj) == 0);
			}
		
		};
		
		Arrays.sort(frags, c);
		
		return frags;
	}

	public static void main(String[] args) {
	
	
			if (args.length!= 5)
				System.out.println("use: RepeatAligner approxCut(0/1) heuristicAlign(0/1) | noAlign(3,3) forceFactor startDir startFile repeatExt");
	
			approxCut= Integer.parseInt(args[0]);
			heuristicAlign= Integer.parseInt(args[1]);
	//		if (approxCut< 2|| heuristicAlign< 2)
	//			align= true;
			
	//		forceFactor= Double.parseDouble(args[2]);
	
			RepeatAligner myRepAli= new RepeatAligner(
				"D:\\workspace\\G-Phase\\",
				args[5]);
			System.exit(0);		
	
			String path2G6= "D:\\Eigene Dateien\\repeats\\data\\ref6";
			if (OSChecker.isLinux())
				path2G6= "/home/ug/msammeth/bali/repeats/ref6";
			else if (OSChecker.isSunOS())
				path2G6= "/homes/micha/bali/repeats/ref6";
	
			conqerBaliBase(
				path2G6,
				args[2],
				args[3],	// ktn+trk_ref6
				args[4]
			);		
		}
	
public static MultiFrag[] sortFragmentsWGT(MultiFrag[] frags) {

	Comparator c= new Comparator() {
		/* (non-Javadoc)
		 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
		 */
		public int compare(Object o1, Object o2) {

			MultiFrag frag1= (MultiFrag) o1;
			MultiFrag frag2= (MultiFrag) o2;

			if (frag1.getWeight() < frag2.getWeight())
				return 1; // reverse order for ascending weights
			if (frag1.getWeight() == frag2.getWeight())
				return 0;
			return (-1); // reverse order for ascending weights
		}
		/* (non-Javadoc)
		 * @see java.util.Comparator#equals(java.lang.Object)
		 */
		public boolean equals(Object obj) {

			return (compare(this, obj) == 0);
		}

	};

	Arrays.sort(frags, c);

	return frags;
}	
	
	public MultiFragExt[] extractFragments(Alignment ali, int nextFragNb, boolean oneBase) {
		
		String[] names= ali.getSeqIDs();
		Vector fragVec= new Vector();
		int aliLength= ali.getSequences()[0].length();
		
		for (int i= 0; i < names.length; i++) {
				// find separator
			int sep1= names[i].lastIndexOf('_');
			try {
				Integer.parseInt(names[i].substring(sep1+1, names[i].length()));
			} catch (Exception e) {
				sep1= names[i].length(); // Error while number parsing --> no number separator (eg. SPN_CHICK)
			}
				// find name
			int remap1= 0;
			String search= names[i].substring(0, sep1); 
			sep1=(sep1== (-1))?names[i].length():sep1;
			for (remap1= 0; remap1< baseNames.length; ++remap1)
				if (search.equalsIgnoreCase(baseNames[remap1]))
					break;
			if (remap1>= baseNames.length) {
//				System.err.println("Found no "+ search+ "!");	// ok, alignmnet can contain more seqs
				continue;
			}
			
			for (int j= (i+1); j < names.length; j++) {
					// find separator
				int sep2= names[j].lastIndexOf('_');
				if ((sep1!= -1)&& (sep2!= -1)) { 
					// check for number counter at end of name
					try {
						Integer.parseInt(names[j].substring(sep2+1, names[j].length()));
					} catch (Exception e) {
						sep2= names[j].length(); // Error while number parsing 
					}
					if (names[i].substring(0, sep1).equals(names[j].substring(0, sep2)))
						continue;	// identical sequences --> skip
				}
					// find name
				int remap2= 0;
				sep2=(sep2== (-1))?names[j].length():sep2;
				search= names[j].substring(0, sep2);
				for (remap2= 0; remap2< baseNames.length; ++remap2)
					if (search.equalsIgnoreCase(baseNames[remap2]))
						break;
				if (remap2>= baseNames.length) {
//					System.err.println("Found no* "+ search+ "!");	// ok, alignments can contain more seqs
					continue;
				}
			
						// create new fragment
			    MultiFragExt newFrag= new MultiFragExt();
			    int start1= ali.getStartPos()[i];
			    int start2= ali.getStartPos()[j];
			    if (oneBase) {
			    	start1--;
			    	start2--;
			    }
			    newFrag.setSequenceStarts(start1, start2);	
			    
				newFrag.setSequenceNos(remap1,remap2);

				int delta= 0, jmp1= 0, jmp2= 0, jmpTot1= 0, jmpTot2= 0;
				for (int l= 0; l < aliLength; l++) {
					if (MultipleAlignmentModel.isGapChar(ali.getSequences()[i].charAt(l)))
						if (MultipleAlignmentModel.isGapChar(ali.getSequences()[j].charAt(l)))
							++delta;	// count gap-pairs
						else {
							++jmp1;		// count gaps for jump in seq2
							if (jmp2> 0) {	// close jmp in 2
								jmpTot2+= jmp2;
								newFrag.addJumpPoint(l- jmpTot2- delta,(-1)* jmp2,false);
								jmp2= 0; 
							}
					} else
						if (MultipleAlignmentModel.isGapChar(ali.getSequences()[j].charAt(l))) {
							++jmp2;		// count gaps for jump in seq1
							if (jmp1> 0) {	// close jmp in 2
								jmpTot1+= jmp1;
								newFrag.addJumpPoint(l- jmpTot1- delta,(-1)* jmp1,true);
								jmp1= 0; 
							}
							// both are chars
						} else if (jmp1> 0) {	// close jmp in 1
							jmpTot1+= jmp1;
							newFrag.addJumpPoint(l- jmpTot1- delta,(-1)* jmp1,true);
							jmp1= 0; 
						} else if (jmp2> 0) {	// close jmp in 2
							jmpTot2+= jmp2;
							newFrag.addJumpPoint(l- jmpTot2- delta,(-1)* jmp2,false);
							jmp2= 0; 
						}
				}
						
					// close open jmp points at end
				if (jmp1> 0) {	// close jmp in 1
					jmpTot1+= jmp1;
					newFrag.addJumpPoint(aliLength- jmpTot1- delta,(-1)* jmp1,true);
					jmp1= 0; 
				} else if (jmp2> 0) {	// close jmp in 2
					jmpTot2+= jmp2;
					newFrag.addJumpPoint(aliLength- jmpTot2- delta,(-1)* jmp2,false);
					jmp2= 0; 
				}
				newFrag.setLength(aliLength- delta);	// set effective length

//				System.out.println("frag for ("+ ali.getSeqIDs()[i]+","+ ali.getSeqIDs()[j]+") "+" -["+remap1+","+remap2+"]-> ("+ baseNames[remap1]+","+baseNames[remap2]+")");
//				System.out.println(newFrag.getSubsequences(baseSequences)[0]+ "\n"+ newFrag.getSubsequences(baseSequences)[1]);

						// chk (compare 1st and last character)
				if (((!MultipleAlignmentModel.isGapChar(newFrag.getSubsequence(baseSequences[remap1], true).charAt(0)))
						&& (!MultipleAlignmentModel.isGapChar(ali.getSequences()[i].charAt(0)))
						&& (Character.toUpperCase(ali.getSequences()[i].charAt(0))!= Character.toUpperCase(newFrag.getSubsequence(baseSequences[remap1], true).charAt(0))))
					|| ((!MultipleAlignmentModel.isGapChar(newFrag.getSubsequence(baseSequences[remap1], true).charAt(newFrag.getSubsequence(baseSequences[remap1], true).length()- 1)))
						&& (!MultipleAlignmentModel.isGapChar(ali.getSequences()[i].charAt(ali.getSequences()[i].length()- 1)))
						&& (Character.toUpperCase(ali.getSequences()[i].charAt(ali.getSequences()[i].length()- 1))!= Character.toUpperCase(newFrag.getSubsequence(baseSequences[remap1], true).charAt(newFrag.getSubsequence(baseSequences[remap1], true).length()- 1))))) {
					System.out.println("MAP ERROR:");
					System.out.println(ali.getSequences()[i]+ "\t"+ newFrag.getSubsequence(baseSequences[remap1], true)); 
					System.out.println(ali.getSequences()[j]+ "\t"+ newFrag.getSubsequence(baseSequences[remap2], false)+"\n--"); 
				} if (((!MultipleAlignmentModel.isGapChar(newFrag.getSubsequence(baseSequences[remap2], false).charAt(0)))
						&& (!MultipleAlignmentModel.isGapChar(ali.getSequences()[j].charAt(0)))
						&& (Character.toUpperCase(ali.getSequences()[j].charAt(0))!= Character.toUpperCase(newFrag.getSubsequence(baseSequences[remap2], false).charAt(0))))
					|| ((!MultipleAlignmentModel.isGapChar(newFrag.getSubsequence(baseSequences[remap2], false).charAt(newFrag.getSubsequence(baseSequences[remap2], false).length()- 1)))
						&& (!MultipleAlignmentModel.isGapChar(ali.getSequences()[j].charAt(ali.getSequences()[i].length()- 1)))
						&& (Character.toUpperCase(ali.getSequences()[j].charAt(ali.getSequences()[j].length()- 1))!= Character.toUpperCase(newFrag.getSubsequence(baseSequences[remap2], false).charAt(newFrag.getSubsequence(baseSequences[remap2], false).length()- 1))))) {
					System.out.println("MAP ERROR:");
					System.out.println(ali.getSequences()[j]+ "\t"+ newFrag.getSubsequence(baseSequences[remap2], false));
					System.out.println(ali.getSequences()[i]+ "\t"+ newFrag.getSubsequence(baseSequences[remap1], true)+"\n--");
				} 
				
					// add
				fragVec.add(newFrag);
					
			}
		}
			
			// convert to array
		MultiFragExt[] frags= new MultiFragExt[fragVec.size()];
		for (int i= 0; i < frags.length; i++) {
			frags[i]= (MultiFragExt) fragVec.elementAt(i);
			frags[i].setNumber(nextFragNb++);
			frags[i].setTranslation(MultiFragExt.TYPE_REPEAT);
		}
		
		return frags;
	}
	
	public String[] alignAllSequences() {

		String[] layout= alignTCoffee(baseSequences);

		int[][] pointer= new int[tokenizedSequences.length][];
		for (int i= 0; i < pointer.length; i++) {

			int cntIV= 0; // get # of IV
			for (int j= 0; j < tokenizedSequences[i].length; ++j)
				if (!tokenizedRepeats[i][j])
					++cntIV;

			//  keep repeat-insert positions
			pointer[i]= new int[cntIV];
			int ptr= 0;
			cntIV= 0;
			for (int j= 0; j < tokenizedSequences[i].length; j++)
				if (!tokenizedRepeats[i][j]) {

					pointer[i][cntIV++]= ptr; // save starts
					ptr += tokenizedSequences[i][j].length();
				} else
					ptr += tokenizedSequences[i][j].length();
		}

		layout= reformatLayout(layout, pointer); 
		return layout;
	}

	public String[] alignInterveningSequences() {

			// concatenate Sequences
		int[][] pointer= new int[tokenizedSequences.length][];
		String[] concatSeqs= new String[tokenizedSequences.length];
		for (int i= 0; i < pointer.length; i++) {

			int cntIV= 0; // get # of IV
			for (int j= 0; j < tokenizedSequences[i].length; ++j)
				if (!tokenizedRepeats[i][j])
					++cntIV;

			// concatenate and keep repeat-insert positions
			pointer[i]= new int[cntIV];
			concatSeqs[i]= "";
			int ptr= 0;
			cntIV= 0;
			for (int j= 0; j < tokenizedSequences[i].length; j++)
				if (!tokenizedRepeats[i][j]) {

					if (((cntIV== 0)&& (j!= 0))|| (cntIV!= 0)) // only 0 if seq starts w rep
						pointer[i][cntIV++]= ptr; // save starts
					concatSeqs[i] += tokenizedSequences[i][j];
					ptr += tokenizedSequences[i][j].length();
				}
		}

/*		for (int i= 0; i < pointer.length; i++) {
			for (int j= 0; j< pointer[i].length; j++) 
				System.out.print(pointer[i][j]+ ", ");
			System.out.println();
		}
		for (int i= 0; i < concatSeqs.length; i++) {
			System.out.println(concatSeqs[i]);
		}
*/		
		String[] layout= alignDialign(concatSeqs);

		// re-insert repeats
		layout= reformatLayout(layout, pointer);

		return layout;
	}
	
public MultiFragExt[] getIVFragments() {

		// concatenate Sequences
	int[][] pointer= new int[tokenizedSequences.length][];
	concatSeqs= new String[tokenizedSequences.length];
	pointer= getIVConcatTokens(pointer);

//	System.out.println("concat IV seqs:");
//	for (int i= 0; i < concatSeqs.length; i++) 
//		System.out.println(concatSeqs[i]);
//	System.out.println();
/*	for (int i= 0; i < pointer.length; i++) {
		for (int j= 0; j < pointer[i].length; j++) {
			System.out.print(pointer[i][j]+",");
		}
		System.out.println();
	}
*/
		// align sequences
	System.out.print("alinging IV...");
	DialignWrapper diaWrap= new DialignWrapper(System.getProperty("user.dir"));
	String fName= writeOutTemp(baseNames, concatSeqs);
	diaWrap.setOutputMSF(true);
	diaWrap.runDialign(fName);
	System.err.flush();
	System.out.println("done.");
		
		// get dialign results 
	MultiFrag[] frags= getDialignFrags(fName);
	MultiFragExt[] fext= new MultiFragExt[frags.length];
	for (int i= 0; i < frags.length; i++) {			// convert
		fext[i]= new MultiFragExt(frags[i]);
		fext[i].setTranslation(MultiFragExt.TYPE_INTERVENING);
	}
	// get names (sequences necessary?)
	MSFWrapper msfReader= new MSFWrapper(fName+ ".ms");
	try {
		msfReader.read();		
	} catch (Exception e) {
		e.printStackTrace();
	}
			
	String[] dnames= msfReader.getSeqNames();		// necessary?!
	for (int i= 0; i < dnames.length; i++) {
		if (!baseNames[i].toLowerCase().startsWith(dnames[i].toLowerCase()))
			System.err.println(
				"Realign name mismatch: "+baseNames[i]+ " <> "+ dnames[i]
			);
	}

		// re-insert repeats
	fext= remapFrags(fext, pointer, concatSeqs);
	
	if (DEBUG) {
		sortFragmentsWGT(fext);
		System.out.println("iv frags");
		System.out.println("========");
		for (int i= 0; i < fext.length; i++) 
			System.out.println(fext[i]);
	}
	
	return fext;
}

protected int[][] getIVConcatTokens(int[][] pointer) {
	
	for (int i= 0; i < tokenizedSequences.length; i++) {	// for all sequences
		
		int cntRep= 0;					// count repeats to init array
		for (int j= 0; j < tokenizedSequences[i].length; j++) 
			if (tokenizedRepeats[i][j])
				++cntRep;
		pointer[i]= new int[cntRep];
		
		concatSeqs[i]= "";
		int concatPtr= 0;
		int k= 0;
		for (int j= 0; j < tokenizedSequences[i].length; j++) {
			
			if (tokenizedRepeats[i][j]) {		// repeat
				pointer[i][k++]= concatPtr;
			} else {							// IV fragment
				concatSeqs[i]+= tokenizedSequences[i][j];
				concatPtr+= tokenizedSequences[i][j].length();
			}
		}
	}	
		
	return pointer;

}

protected int[][] getIVConcatTokens_old(int[][] pointer) {
	
	for (int i= 0; i < pointer.length; i++) {

		int cntIV= 0; // get # of IV
		for (int j= 0; j < tokenizedSequences[i].length; ++j)
			if (!tokenizedRepeats[i][j])
				++cntIV;

		// concatenate and keep repeat-insert positions
		if (!tokenizedRepeats[i][0])
		if (!tokenizedRepeats[i][tokenizedRepeats[i].length- 1])
			--cntIV;	// correct if seq ends with IV
		pointer[i]= new int[cntIV];
		concatSeqs[i]= "";
		int ptr= 0;
		cntIV= 0;
		for (int j= 0; j < tokenizedSequences[i].length; j++)
			if (!tokenizedRepeats[i][j]) {

				if (((cntIV== 0)&& (j!= 0))|| (cntIV!= 0)) // only 0 if seq starts w rep
					pointer[i][cntIV++]= ptr; // save starts
				concatSeqs[i] += tokenizedSequences[i][j];
				ptr += tokenizedSequences[i][j].length();
			}
	}
	
	return pointer;
}


public MultiFrag[] getRelaxedFragments(int nextFragNb) {
	
		// align sequences
	System.out.print("alinging All...");
	DialignWrapper diaWrap= new DialignWrapper(System.getProperty("user.dir"));	// System.getProperty("user.dir"):
													// not when using wrapper 2nd time
													// change !!!
	String fName= writeOutTemp(baseNames, baseSequences);
	diaWrap.setOutputMSF(true);
	diaWrap.runDialign(fName);
	System.out.println("done.");
			
		// get dialign results 
	MultiFrag[] frags= getDialignFrags(fName);
		
		// collect all fragments x repeat
	Repeat tmpRepeat= repeats;
	Vector rlxFrags= new Vector();
	while (tmpRepeat!= null) {

			// collect/clone all fragments x repeat
		Vector repFrags= new Vector();
		for (int i= 0; i < frags.length; i++) 
			if (frags[i].intersects(tmpRepeat))
				repFrags.add(frags[i].clone());
		
			// sort acc. to 2nd seq#
		for (int i= 0; i < baseSequences.length; i++) {
			if (i== tmpRepeat.getSeqNb())
				continue;
			Vector repFragsSort= new Vector(repFrags.size());
			for (int j= 0; j < repFrags.size(); j++) {
				
				if (((MultiFrag) repFrags.elementAt(j)).getOtherSequenceNb(tmpRepeat.getSeqNb())== i)
					repFragsSort.add(repFrags.elementAt(j));
			}
				
				// submit for filtering
			MultiFrag tmpFrag= getRelaxedFragments(tmpRepeat, repFragsSort, nextFragNb);
			if (tmpFrag!= null) {
				rlxFrags.add(tmpFrag);
				++nextFragNb;
			}
		}
			
		tmpRepeat= tmpRepeat.getNext();
	}
	
		// convert
	MultiFrag[] result= new MultiFrag[rlxFrags.size()];
	for (int i= 0; i < result.length; i++) 
		result[i]= (MultiFrag) rlxFrags.elementAt(i);
	
	return result;
}

/**
 * At the moment only for normal fragments
 * @param baseRep
 * @param frags 	<code>MultiFrag</code>s which are overlapping with 
 * 					<code>baseRepeat</code>
 * @return
 */
protected MultiFragExt getRelaxedFragments(Repeat baseRep, Vector frags, int nextFragNb) {
	
		// abort
	if ((frags== null)|| (frags.size()< 1))
		return null;
		
		// cut fragments to repeat area
	for (int i= 0; i < frags.size(); i++) {		
		
			// a) cut them to overlap size
		MultiFrag tmpFrag= (MultiFrag) frags.elementAt(i);
			// start before
		if (tmpFrag.getSequenceStart(baseRep.getSeqNb())< baseRep.getStart()) {
			int delta= (baseRep.getStart()- tmpFrag.getSequenceStart(baseRep.getSeqNb()));
			tmpFrag.setLength(tmpFrag.getLength()- delta);
			tmpFrag.setSequenceStart(true, tmpFrag.getSequenceStart(true)+ delta);
			tmpFrag.setSequenceStart(false, tmpFrag.getSequenceStart(false)+ delta);
		}
			// ends after
		if (tmpFrag.getSequenceStart(baseRep.getSeqNb())+ tmpFrag.getLength()
			> baseRep.getStart()+ baseRep.getLength()) 
			
			tmpFrag.setLength(tmpFrag.getLength()- 
				((tmpFrag.getSequenceStart(baseRep.getSeqNb())+ tmpFrag.getLength())
				- (baseRep.getStart()+ baseRep.getLength()))
			);
			// check fragment integrity
		if (tmpFrag.getLength()< 1) {
			frags.remove(i--);
			continue;
		}
		
			// b) cut them if overlapping other repeats
/*		Repeat tmpRepeat= repeats;
		while (tmpRepeat!= null) {
			if (tmpRepeat.getSeqNb()== baseRep.getSeqNb()) {// repeats NON-OVERLAPPING: 
															// incl. tmpRepeat.equals(baseRep)
				tmpRepeat= tmpRepeat.getNext();
				continue;
			}
			if (tmpFrag.intersects(tmpRepeat)) {
				MultiFrag purgedFrag= getRelaxedPurgeRepeat(tmpRepeat, (MultiFrag) frags.remove(i));
				if (purgedFrag== null) {					// eliminated one
					--i;
					break;
				} else {
					frags.insertElementAt(purgedFrag, i);	// modified one (not add())s
					tmpFrag= purgedFrag;
				}
			}
			tmpRepeat= tmpRepeat.getNext();
		}
*/		
	}
	
		// assign weights
	MultiFrag[] help= new MultiFrag[frags.size()];
	for (int i= 0; i < frags.size(); i++) {
		MultiFrag tmpFrag= (MultiFrag) frags.elementAt(i);
		tmpFrag.setWeight(calcRelWgt(tmpFrag));
		help[i]= tmpFrag;
	}
	if (frags.size()== 0)
		return null;
	
		// make consistent and sort ascending to start position
	//frags= filterPWOptimal(frags);		// note: here can still be overlapping (consistent) fragments !
	MultiFrag[] help2= filterPW(help, false, false);
	frags= new Vector();
	for (int i= 0; i < help2.length; i++) {
		frags.add(help2[i]);
	}
	
	MultiFrag[] cFrags= new MultiFrag[frags.size()];
	for (int i= 0; i < cFrags.length; i++) 
		cFrags[i]= (MultiFrag) frags.elementAt(i);
	MultiFragExt.sortFragmentsStartPoints(cFrags, baseRep.getSeqNb());
		
		// prepare attributes
	String[] names= new String[2];
	names[0]= this.baseNames[baseRep.getSeqNb()];
	names[1]= this.baseNames[cFrags[0].getOtherSequenceNb(baseRep.getSeqNb())];
	
	int start1= cFrags[0].getSequenceStart(baseRep.getSeqNb());		// start and end of repeat area			
	int end1= -1;
	for (int i= 0; i < cFrags.length; i++) 
		end1= Math.max(end1, cFrags[i].getSequenceStart(baseRep.getSeqNb())+ 
								cFrags[i].getLength());
	int start2= Integer.MAX_VALUE;											// find lowest start pos
	for (int i= 0; i < cFrags.length; i++) 
		start2= Math.min(start2, cFrags[i].getSequenceStart(cFrags[i].getOtherSequenceNb(baseRep.getSeqNb())));
	int end2= -1;
	for (int i= 0; i < cFrags.length; i++) 
		end2= Math.max(end2, cFrags[i].getSequenceStart(cFrags[i].getOtherSequenceNb(baseRep.getSeqNb()))+ 
								cFrags[i].getLength());
	String[] seqs= new String[2];
	seqs[0]= baseSequences[baseRep.getSeqNb()].substring(start1, end1);
	seqs[1]= baseSequences[cFrags[0].getOtherSequenceNb(baseRep.getSeqNb())].
		substring(start2, end2);
	
	for (int i= 0; i < cFrags.length; i++) {						// adapt thy frags to cutoff repeat area
		cFrags[i].setSequenceStart(baseRep.getSeqNb(), 
			cFrags[i].getSequenceStart(baseRep.getSeqNb())- start1
		);
		cFrags[i].setSequenceStart(cFrags[i].getOtherSequenceNb(baseRep.getSeqNb()), 
			cFrags[i].getSequenceStart(cFrags[i].getOtherSequenceNb(baseRep.getSeqNb()))- start2
		);
		if (cFrags[i].getSequenceNo(true)== baseRep.getSeqNb())
			cFrags[i].setSequenceNos(0,1);
		else
			cFrags[i].setSequenceNos(1,0);
	}
	
	
		// align good
	DCAClosure dca= new DCAClosure();
//	dca.setOutFileName(diaFasta.getAbsFileName());
	dca.setCostTable(CostTable.BLOSUM62);
	dca.setSequences(names, seqs, CostTable.BLOSUM62);
	dca.setFragments(cFrags);

	dca.setApproximate(false);
	dca.setOutput(true);
	dca.setRecursionStop(40);
//	dca.setOutputFormat(QDivide.MSF);
//	dca.dialignTime= (float) ((System.currentTimeMillis()- myStartTime)/ 1000);
	dca.setOutput(false);
	try {dca.run();} 
	catch (CancelException e) {
		System.err.println(e);; // nothing
	}
	
		// extract super-frags
	Alignment ali= new Alignment();
	ali.setSeqIDs(names);
	ali.setStartPos(new int[] {start1+ 1, start2+ 1});	// 1-based
	ali.setEndPos(new int[] {end1, end2});
	ali.setSequences(dca.getAlignmentLayout());
	MultiFragExt[] extract= extractFragments(ali, nextFragNb, true);
	if (extract.length!= 1)
		System.err.println("extracted more than one super-frag!");
	scoreFragments(extract);
	extract[0].setTranslation(12);	// undefined type
	
	
	return extract[0];
}
/**
 * 
 * @param aRepeat	repeat intersecting with frag
 * @param frag
 * @return
 */
protected MultiFrag getRelaxedPurgeRepeat(Repeat aRepeat, MultiFrag frag) {
	
	int start1= aRepeat.getStart();
	int start2= frag.getSequenceStart(aRepeat.getSeqNb());
	int end1= aRepeat.getStart()+ aRepeat.getLength();
	int end2= start2+ frag.getLength();
	
	if (start2>= start1&& start2< end1)	// if starts within repeat, start after
		start2= end1;
	if (end2>= start1&& end2< end1)		// if ends within, has to end before
		end2= start1;
	
		// new characteristics
	if ((end2- start2)< 1)					// check b4 modifying reference
		return null;
	frag= (MultiFrag) frag.clone();
	int delta= frag.getSequenceStart(aRepeat.getSeqNb())+ frag.getLength()- end2; // can only fall
	if (delta> 0) 
		frag.setLength(frag.getLength()- delta);								  // just correct length
	delta= start2- frag.getSequenceStart(aRepeat.getSeqNb());		// can only rise
	if (start2!= frag.getSequenceStart(aRepeat.getSeqNb())) {
		frag.setSequenceStart(true, frag.getSequenceStart(true)+ delta);
		frag.setSequenceStart(false, frag.getSequenceStart(false)+ delta);
		frag.setLength(frag.getLength()- delta);					// correct all
	}
	
	return frag;	// else..
}

	
	public void writeOut(String[] layout) {

		MSFWrapper msf= new MSFWrapper("C:\\tcoffee.out");
		msf.setSequences(layout);
		msf.setSeqNames(baseNames);
		try {
			msf.write();
		} catch (Exception e) {
			; //
		}
	}
	
	public String[] alignDialign(String[] concatSeqs) {
		
			// align sequences
		DialignWrapper diaWrap= new DialignWrapper(System.getProperty("user.dir"));
		String fName= writeOutTemp(baseNames, concatSeqs);
		diaWrap.setOutputMSF(true);
		diaWrap.runDialign(fName);
		
			// get result
		MSFWrapper msfReader= new MSFWrapper(fName+ ".ms");
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
			
		return msfReader.getSequences();
	}
	
	public MultiFrag[] getDialignFrags(String fName) {
		
			// get result
		DialignFASTAWrapper diaFasta= new DialignFASTAWrapper(fName);
		try {
			diaFasta.readFF();
		} catch (Exception e) {
			; // merkt man schon.
		}
		MultiFrag[] frags= diaFasta.getFragments();
//		System.out.println(frags.length);
			
		return frags;
	}
	
	public String[] alignTCoffee(String[] concatSeqs) {
		
			// align sequences
		String fName= writeOutTemp(baseNames, concatSeqs);
		TCoffeeWrapper wrap= new TCoffeeWrapper(fName);
		wrap.setOutputMSF(true);
		wrap.runTCoffee(fName, "");
		
			// get result
		MSFWrapper msfReader= new MSFWrapper(fName.substring(0, fName.lastIndexOf('.'))+ ".msf_aln");
		try {
			msfReader.read();		
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		String[] layout= msfReader.getSequences();
		String[] mapNames= msfReader.getSeqNames();
		String[] mapSeqs= new String[layout.length];
		for (int i= 0; i < mapSeqs.length; i++) {
			mapSeqs[i]= layout[i];
		}
		for (int i= 0; i < mapNames.length; i++) {
			if (!mapNames[i].equalsIgnoreCase(baseNames[i]))
				for (int j= 0; j < baseNames.length; j++) {
					if (mapNames[i].equalsIgnoreCase(baseNames[j]))
						mapSeqs[i]= layout[j];
				}
		}
		
		return mapSeqs;
	}
	
	public String[] reformatLayout(String[] layout, int[][] pointer) {
		
/*		System.out.println("---");
		for (int i= 0; i< layout.length; ++i)
			System.out.println(layout[i]);
		for (int i= 0; i< pointer.length; ++i) {
			for (int j= 0; j< pointer[i].length; ++j) 
				System.out.print(pointer[i][j]+", ");
			System.out.println();
		}
		System.out.println("---");
*/		
			
		String[] newLayout= new String[layout.length];
		for (int i= 0; i < newLayout.length; i++) {
			newLayout[i]= layout[i];
		}
		String gap= ".-~0123456789|";
		for (int i= 0; i< pointer.length; ++i){
			int lastPos= 0;
			for (int j= 0; j< pointer[i].length; ++j) {
					// find position
				int pos= 0;
				for (int ctr= 0; ctr< pointer[i][j];++pos)
					if (gap.indexOf(newLayout[i].charAt(pos))== (-1))
							ctr++;
				
					// insert across all sequences
				for (int k= 0; k< newLayout.length; ++k)
					newLayout[k]= newLayout[k].subSequence(0,pos)
						+"||"+ i+"||"
						+newLayout[k].substring(pos, newLayout[k].length());
				lastPos= pos;

			}
			
		}
		return newLayout;
	}
	
	/**
	 * Re-assigns position-indices for the fragments for changes caused by inserting the repeats.
	 * Also re-calculates the relative weight <code>calcRelWgt</code> and the overlapping weight
	 * <code>calcOlWegt</code> for the new fragments (i.e., for the new length of the base sequences).
	 * 
	 * @param frags
	 * @param pointer
	 * @param concatSeqs
	 * @return
	 */
	public MultiFragExt[] remapFrags(MultiFragExt[] frags, int[][] pointer, String[] concatSeqs) {
		
/*		System.out.println("--- remap ---");
		for (int i= 0; i < baseSequences.length; i++) 
			System.out.println(baseSequences[i]);
		System.out.println("---");
		for (int i= 0; i < concatSeqs.length; i++) 
			System.out.println(concatSeqs[i]);
		System.out.println("---");
		for (int i= 0; i< frags.length; ++i)
			System.out.println(frags[i]);
		System.out.println("---");
*/		
		System.out.print("remapping fragments...");
		
			// clone frags
			// maybe cancel
		MultiFragExt[] newFrags= new MultiFragExt[frags.length];
		for (int i= 0; i < newFrags.length; i++) {
			newFrags[i]= (MultiFragExt) frags[i].clone();
			if (i> 0) {
				newFrags[i].setPred(newFrags[i-1]);
				newFrags[i-1].setNext(newFrags[i]);
			}
		}

			// for all frags
		for (int f= 0; f < newFrags.length; f++) {
			
			MultiFrag oldFrag= (MultiFrag) newFrags[f].clone();		// guard for check
			
			int a= newFrags[f].getSequenceNo(true);
			int b= newFrags[f].getSequenceNo(false);
			int x= newFrags[f].getSequenceStart(true);
			int y= newFrags[f].getSequenceStart(false);
			
				// sequence 1
				// move start
			int insRep= 0;	// count # of repeats to insert before fragment
			for (int i= 0; i < pointer[a].length; i++)	// perform changes for seqA 
				if (x>= pointer[a][i])
					++insRep;
			int insChar= 0;	// add up # of chars to insert
			for (int i= 0, reps= 0; reps< insRep; i++)
				if (tokenizedRepeats[a][i]) {
					insChar+= tokenizedSequences[a][i].length();
					++reps;
				}
			newFrags[f].setSequenceStart(true, (x+ insChar));	
			int nextRep= 0;								// jump points
			for (int i= 0; i < pointer[a].length; i++) {
				while(!tokenizedRepeats[a][nextRep])	// find next repeat (for length info)
					++nextRep;					
				if ((pointer[a][i]> x)&& (pointer[a][i]< x+ newFrags[f].getLength())) // all concat-based numbers
					newFrags[f].addJumpPoint(pointer[a][i]- x, tokenizedSequences[a][nextRep].length(), true);
				++nextRep;
			}
			
				// 2nd seq
			insRep= insChar= 0;
			for (int i= 0; i < pointer[b].length; i++)	// perform for seqB 
				if (y>= pointer[b][i])
					++insRep;
			for (int i= 0, reps= 0; reps< insRep;++i)
				if (tokenizedRepeats[b][i]) {
					insChar+= tokenizedSequences[b][i].length();
					++reps;
				}
			newFrags[f].setSequenceStart(false, (y+ insChar));
			nextRep= 0;							// jump ponts
			for (int i= 0; i < pointer[b].length; i++) {
				while(!tokenizedRepeats[b][nextRep])	// find next repeat
					++nextRep;					
				if ((pointer[b][i]> y)&& (pointer[b][i]< y+ newFrags[f].getLength())) 
					newFrags[f].addJumpPoint(pointer[b][i]- y, tokenizedSequences[b][nextRep].length(), false);
				++nextRep;
			}
			
				// recalc wgt (according to new length)
			newFrags[f].setWeight(calcRelWgt(newFrags[f]));
			
			String s1= oldFrag.getSubsequence(concatSeqs[oldFrag.getSequenceNo(true)], true);
			String s2= oldFrag.getSubsequence(concatSeqs[oldFrag.getSequenceNo(false)], false);
			String s3= newFrags[f].getSubsequence(baseSequences[newFrags[f].getSequenceNo(true)], true);
			String s4= newFrags[f].getSubsequence(baseSequences[newFrags[f].getSequenceNo(false)], false);
			if(!s1.equalsIgnoreCase(s3)|| !s2.equalsIgnoreCase(s4)) {
				System.err.println(oldFrag);
				System.err.println("oldFrag\t"+s1+"\tnewFrag\t"+s3);
				System.err.println("\t\t"+s2+"\t\t\t"+s4);
			}
			if (numericWeightCheck&& oldFrag.getWeight()> 0f&& Math.round(newFrags[f].getWeight())== 0)
				System.err.println("Warning: old w= "+oldFrag.getWeight()+"  => new w= "+newFrags[f].getWeight());

				// check calcWgt	
/*			if (calcRelWgt(newFrags[f])!= newFrags[f].getWeight())
				System.err.println(
					calcRelWgt(newFrags[f])+ "!="+
					newFrags[f].getWeight()+ "in"+
					newFrags[f]
				);
*/
/*				// check seq lengths
			String oldS= 
				concatSeqs[a].substring(
								frags[f].getSequenceStart(true), 
								frags[f].getSequenceStart(true)+ frags[f].getLength() 
								)+ "\t"+
								concatSeqs[b].substring(
								frags[f].getSequenceStart(false), 
								frags[f].getSequenceStart(false)+ frags[f].getLength()
							);
			String newS=
				newFrags[f].getSubsequence(baseSequences[a], true)
				+ "\t"+
				newFrags[f].getSubsequence(baseSequences[b], false);

			if (!oldS.equals(newS))	{
				System.out.println(frags[f]+"\n"+newFrags[f]);
				System.out.println(insRep+" -> "+insChar);							
				System.out.println(oldS+"\n"+newS);
				System.out.println("---");
			}
*/						
		}
		
			// recalc olw (according to new length)
//		calcOLW(newFrags); 
//		calcOlWegt_old(newFrags);
		
		System.err.flush();
		System.out.println("done.");
		return newFrags;
	}
	
/**
 * Adds up the weights of the respective sub-diagonals.
 * (Motivated by traditions established in mechanics and informatin theory
 * considered the neg. logarithm of a probability, the values can be added).
 * @param frag
 * @return
 */
/*
 * see rel_wgt_calc() [in wgt.c]
 * opened up to double for numerical instability
 * 
 * (with rounding in end consistent with c-code,
 * but using 2-floats, not 6-floats as in c)
 * 
 */
public strictfp float calcRelWgt_indel(MultiFrag frag) {

		// derive attributes
	int l1= baseSequences[frag.getSequenceNo(true)].length();
	int l2= baseSequences[frag.getSequenceNo(false)].length();
		
		// init basic values
	double factor= ((double) (l1* l2))/ 400d;
	String s1= frag.getSubsequence(baseSequences[frag.getSequenceNo(true)], true).toUpperCase();
	String s2= frag.getSubsequence(baseSequences[frag.getSequenceNo(false)], false).toUpperCase();
	int[] mapping= Constants.getMapBChars();
	int[][] blosum= Constants.getBlosum();

		// find single sub-diagonals and add up weight
	int pos= 0, lpos= 0;
	float sumWgt= 0f;
	while (pos< s1.length()) {
		for (;(pos< s1.length())
				&& (MultipleAlignmentModel.isGapChar(s1.charAt(pos))
				|| MultipleAlignmentModel.isGapChar(s2.charAt(pos)))
				;++pos)
				; // find new "continue"
		lpos= pos;
		for (;((pos< s1.length())
					&& (pos- lpos< 40))	// max diagonal length (lookup table)
				&& (!(MultipleAlignmentModel.isGapChar(s1.charAt(pos))
					|| MultipleAlignmentModel.isGapChar(s2.charAt(pos))));
				++pos)
				; // find new "breakpoint"
							
		int s_D= 0;
		if (pos== lpos)
			continue;				// otherwise: sim= 0, log= NaN, sumWgt+= NaN --> NaN
		for (int i= lpos; i < pos; ++i) {
			char c1= s1.charAt(i);
			char c2= s2.charAt(i);
			s_D+= blosum[mapping[c1- 65]][mapping[c2- 65]];
		}
		
			// lookup in tp
		double tpr= Constants.getTpProt400()[pos- lpos][s_D];
		
			// correct
		if (tpr> 0.0000000001)						// check rel_wgt_cut !?!
			tpr= 1d- Math.pow((1d- tpr), factor);
		else
			tpr= tpr* factor;
		
			// add check for threshold T
		// if (-log(tpr)< T) return 0f;
		
			// compute rel. Weight
		float relWgt= (-1f)* ((float) Math.log(tpr));

			// check
		if (relWgt== Float.POSITIVE_INFINITY) 
			System.err.println("relWgt infinity2: "+ frag);
		 
		sumWgt+= relWgt;	// add to result
	}
		
	return sumWgt;
}

/**
 * Adds up the weights of the respective sub-diagonals.
 * (Motivated by traditions established in mechanics and informatin theory
 * considered the neg. logarithm of a probability, the values can be added).
 * @param frag
 * @return
 */
/*
 * see rel_wgt_calc() [in wgt.c]
 * opened up to double for numerical instability
 * 
 * (with rounding in end consistent with c-code,
 * but using 2-floats, not 6-floats as in c)
 * 
 */
public strictfp float calcRelWgt_indel_old(MultiFrag frag) {

		// derive attributes
	int l1= baseSequences[frag.getSequenceNo(true)].length();
	int l2= baseSequences[frag.getSequenceNo(false)].length();
//	int l1= concatSeqs[frag.getSequenceNo(true)].length();
//	int l2= concatSeqs[frag.getSequenceNo(false)].length();
		
		// compute characteristics
	double factor= ((double) (l1* l2))/ 400d;
//	System.out.println("factor= "+ factor);
	String s1= frag.getSubsequence(baseSequences[frag.getSequenceNo(true)], true).toUpperCase();
	String s2= frag.getSubsequence(baseSequences[frag.getSequenceNo(false)], false).toUpperCase();
//	System.out.println(s1+ "\n"+ s2+ "\n--");
	int[] mapping= Constants.getMapBChars();
	int[][] blosum= Constants.getBlosum();

		// find single sub-diagonals and add up weight
	int pos= 0, lpos= 0;
	float sumWgt= 0f;
	while (pos< s1.length()) {
		for (;(MultipleAlignmentModel.isGapChar(s1.charAt(pos))
				|| MultipleAlignmentModel.isGapChar(s2.charAt(pos)))
				&& (pos< s1.length());++pos)
				; // find new "continue"
		lpos= pos;
		for (;((pos< s1.length()))	// max diagonal length (lookup table)
				&& (!(MultipleAlignmentModel.isGapChar(s1.charAt(pos))
					|| MultipleAlignmentModel.isGapChar(s2.charAt(pos))));
				++pos)
				; // find new "breakpoint"
							
		int s_D= 0;
		for (int i= lpos; i < pos; ++i) {
			char c1= s1.charAt(i);
			char c2= s2.charAt(i);
			s_D+= blosum[mapping[c1- 65]][mapping[c2- 65]];
		}
//		System.out.println("s_D= "+ s_D);
		
			// lookup in tp
		double tpr= 0d;
		try {
			tpr= Constants.getTpProt400()[pos- lpos][s_D];
			
			// if too big (repeat-frags), decompose
		} catch (ArrayIndexOutOfBoundsException e) {
			return 0f;
		}
						
//		System.out.println("tpr["+l_D+"]["+s_D+"]= "+tpr);

		
			// correct
		if (tpr> 0.0000000001)						// check rel_wgt_cut !?!
			tpr= 1d- Math.pow((1d- tpr), factor);
		else
			tpr= tpr* factor;
//		System.out.println("corrected: "+tpr);
		
			// add check for threshold T
//		if (-log(tpr)< T)
//			return 0f;
			
			// compute rel. Weight
		float relWgt= (-1f)* ((float) Math.log(tpr));
//		double relWgt= (-1f)* ((double) Math.log(tpr));
//		relWgt= Math.round(relWgt* 100f)/ 100f;		// break on 2 decimals
//		System.out.println("wgt: "+ relWgt);

		if (relWgt== Float.POSITIVE_INFINITY) 
//			System.err.println("overlap2: infinity2 "+conslen+","+match+"= "+addWgt);
			continue;	// skip if no entry in lookup-table (=0.0 -> log(0.0)= posInfinity?!!)
		 
		sumWgt+= relWgt;
	}
		
	return sumWgt;
}	/*

	 * see rel_wgt_calc() [in wgt.c]
 	 * opened up to double for numerical instability
 	 * 
	 * (with rounding in end consistent with c-code,
	 * but using 2-floats, not 6-floats as in c)
	 * 
	 */
	public strictfp float calcRelWgt_old(MultiFrag frag) {

			// derive attributes
		int l_D= frag.getLength();
		int l1= baseSequences[frag.getSequenceNo(true)].length();
		int l2= baseSequences[frag.getSequenceNo(false)].length();
//		int l1= concatSeqs[frag.getSequenceNo(true)].length();
//		int l2= concatSeqs[frag.getSequenceNo(false)].length();
		
			// compute characteristics
		double factor= ((double) (l1* l2))/ 400d;
//		System.out.println("factor= "+ factor);
		String s1= frag.getSubsequence(baseSequences[frag.getSequenceNo(true)], true).toUpperCase();
		String s2= frag.getSubsequence(baseSequences[frag.getSequenceNo(false)], false).toUpperCase();
//		System.out.println(s1+ "\n"+ s2+ "\n--");
		int[] mapping= Constants.getMapBChars();
		int[][] blosum= Constants.getBlosum();
		int s_D= 0;
		for (int i= 0; i < l_D; ++i) {
			char c1= s1.charAt(i);
			char c2= s2.charAt(i);
			if (MultipleAlignmentModel.isGapChar(c1)
				|| MultipleAlignmentModel.isGapChar(c2))
				continue;	// skip gaps (for global alignment of repeats)
			s_D+= blosum[mapping[c1- 65]][mapping[c2- 65]]; 
		}
//		System.out.println("s_D= "+ s_D);
		
			// lookup in tp
		double tpr= 0d;
		try {
			tpr= Constants.getTpProt400()[l_D][s_D];
		} catch (ArrayIndexOutOfBoundsException e) {
			return 0f;
		}
//		System.out.println("tpr["+l_D+"]["+s_D+"]= "+tpr);

			// add check for e_match		
//		if (!((max_sim_score* l)> (e_match* l))|| !(-Math.log(tpr)> threshold))
//			return 0f;	

			// correct
		if (tpr> 0.0000000001)						// check rel_wgt_cut !?!
			tpr= 1d- Math.pow((1d- tpr), factor);
		else
			tpr= tpr* factor;
//		System.out.println("corrected: "+tpr);
		
			// add check for threshold T
//		if (-log(tpr)< T)
//			return 0f;
			
			// compute rel. Weight
		float relWgt= (-1f)* ((float) Math.log(tpr));
//		double relWgt= (-1f)* ((double) Math.log(tpr));
//		relWgt= Math.round(relWgt* 100f)/ 100f;		// break on 2 decimals
//		System.out.println("wgt: "+ relWgt);

		return relWgt;
	}
	
 	/*
	 * see rel_wgt_calc() [in wgt.c]
	 * opened up to double for numerical instability
	 * 
	 * (with rounding in end consistent with c-code,
	 * but using 2-floats, not 6-floats as in c)
	 * 
	 */
	public strictfp float calcRelWgt(MultiFrag frag) {

			// derive attributes
		int l_D= frag.getLength();
		int l1= baseSequences[frag.getSequenceNo(true)].length();
		int l2= baseSequences[frag.getSequenceNo(false)].length();
//		int l1= concatSeqs[frag.getSequenceNo(true)].length();
//		int l2= concatSeqs[frag.getSequenceNo(false)].length();
		
			// init base values
		double factor= ((double) (l1* l2))/ 400d;
		String s1= frag.getSubsequence(baseSequences[frag.getSequenceNo(true)], true).toUpperCase();
		String s2= frag.getSubsequence(baseSequences[frag.getSequenceNo(false)], false).toUpperCase();
		int[] mapping= Constants.getMapBChars();
		int[][] blosum= Constants.getBlosum();
		
			// tokenize to 40mers (lookup-table!)
		int x= (l_D/ 40);	
		if ((l_D% 40)!= 0)
			++x;
		float relWgt= 0f;
		for (int l= 0; l< x; ++l) {

				// calc similarity
			int s_D= 0;
			for (int i= (l* 40); (i < l_D)&& (i< l*40+40); ++i) {
				char c1= s1.charAt(i);
				char c2= s2.charAt(i);
				if (MultipleAlignmentModel.isGapChar(c1)
					|| MultipleAlignmentModel.isGapChar(c2))
					continue;	// skip gaps (for global alignment of repeats)
				try {
					s_D+= blosum[mapping[c1- 65]][mapping[c2- 65]]; 
				} catch (ArrayIndexOutOfBoundsException e) {			// set to 0= min of matrix
					System.err.println("Non-supported char ("+ c1+ ","+ c2+ ") in sequences.");
				}
			}

				// lookup in tp
			int currentLD= 40; 						// full tokens of 40
			if ((l== (x- 1))&& (l_D% 40!= 0))		// in last (or unique) round 
				currentLD= l_D% 40; 				// is the rest
			double tpr= Constants.getTpProt400()[currentLD][s_D];	// not l_D
	
				// correct
			if (tpr> 0.0000000001)						// check rel_wgt_cut !?!
				tpr= 1d- Math.pow((1d- tpr), factor);
			else
				tpr= tpr* factor;
			
			// add check for threshold T
		// if (-log(tpr)< T) return 0f;

				// compue wgt, check and add
			float addWgt= (-1f)* ((float) Math.log(tpr));
			if (relWgt== Float.POSITIVE_INFINITY) 
				System.err.println("relWgt infinity: "+ frag);
				// compute rel. Weight
			relWgt+= addWgt;
		}

		return relWgt;
	}
	
	public void alignPWInterveningSequences() {
		
			// align intervening sequence fragments
		int maxLen= 0;
		for (int i= 0; i< tokenizedSequences.length; ++i) 
			for (int j= 0; j< tokenizedSequences[i].length; ++j) 
				if (tokenizedSequences[i][j].length()> maxLen)
					maxLen= tokenizedSequences[i][j].length();
		QAlign.setLength(maxLen);
		int[][][][] scores= new int[tokenizedSequences.length][][][];
		for (int i= 0; i< tokenizedSequences.length; ++i) {
			scores[i]= new int[tokenizedSequences[i].length][][];
			for (int j= 0; j< tokenizedSequences[i].length; ++j) { 

				scores[i][j]= new int[scores.length][];
				for (int k= 0; k< tokenizedSequences.length; ++k) {
					scores[i][j][k]= new int[tokenizedSequences[k].length];
					for (int l= 0; l< scores[i][j][k].length; ++l)
						scores[i][j][k][l]= (-1);
				}
					
				for (int k= (i+1); k< tokenizedSequences.length; ++k) 
					for (int l= 0; l< tokenizedSequences[k].length; ++l) {
					
							// dont align repeats with intervenings
						if(tokenizedRepeats[i][j]!= tokenizedRepeats[k][l])
							continue; 
						
						try {
							QAlign qalign= new QAlign();
							qalign.setCostTable(CostTable.BLOSUM62);
							qalign.setSequences(new String[] {tokenizedSequences[i][j], tokenizedSequences[k][l]});
							qalign.setSimultaneousAlignment(false);
							qalign.setOutput(false);
							qalign.run();
							scores[i][j][k][l]= qalign.getPairwiseCosts()[1][2];
						} catch (CancelException e) {
							; // :)
						}
					}
			}
		}
			// make symmetric matrix
/*		for (int i= 0; i< scores.length; ++i) 
			for (int j= 0; j< scores[i].length; ++j) 
				for (int k= (i+1); k< scores.length; ++k) 
					for (int l= 0; l< scores[k].length; ++l)
						scores[k][l][i][j]= scores[i][j][k][l];
*/
/*						 
		for (int i= 0; i< scores.length; ++i) {
			for (int j= 0; j< scores[i].length; ++j) {
				System.out.print("[");
				for (int k= 0; k< scores.length; ++k) {
					System.out.print("(");
					for (int l= 0; l< scores[k].length; ++l)
						if (scores[i][j][k][l]!= (-1))
							System.out.print(scores[i][j][k][l]+ ", ");
					System.out.print(") ; ");
				}
				System.out.print("]  :  ");
			}
			System.out.println();
		}
*/		

		// build graph
		buildGraph(scores);
	}
	
	protected void buildGraph(int[][][][] scores) {
		
//		Hashtable hash= new Hashtable();
		GraphPanel panel = new GraphPanel();
			
			// init nodes
		String[][] tokenNames= new String[tokenizedSequences.length][];
		for (int i= 0; i< tokenNames.length; ++i) {
			tokenNames[i]= new String[tokenizedSequences[i].length];
			int repCounter= 0;
			for (int j= 0; j< tokenNames[i].length; ++j) {
				tokenNames[i][j]= baseNames[i]+ "_";
				if (tokenizedRepeats[i][j]) {
					repCounter++;
					tokenNames[i][j]+= repCounter; 
				} else {
					tokenNames[i][j]+= repCounter+ ":"+ (repCounter+ 1);
				}
				if (!tokenizedRepeats[i][j])
					panel.addNode(tokenNames[i][j]);
//				if (j> 0)
//					panel.addEdge(tokenNames[i][j], tokenNames[i][j-1], 0);
			}
		}
		
		
			// init edges
		int maxScore= 0;
		for (int i= 0; i< scores.length; ++i) 
			for (int j= 0; j< scores[i].length; ++j) 
				for (int k= (i+1); k< scores.length; ++k) 
					for (int l= 0; l< scores[k].length; ++l)
						if (scores[i][j][k][l]> maxScore)
							maxScore= scores[i][j][k][l];  
		for (int i= 0; i< scores.length; ++i) 
			for (int j= 0; j< scores[i].length; ++j) 
				for (int k= (i+1); k< scores.length; ++k) 
					for (int l= 0; l< scores[k].length; ++l) {
							// skip empty
						if (scores[i][j][k][l]< 0)
							continue;
						
						if (!tokenizedRepeats[i][j])
							panel.addEdge(tokenNames[i][j], tokenNames[k][l], scores[i][j][k][l]/ 10);
					}
					
						
		JFrame frame = new JFrame("Graphtest");
		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		frame.setContentPane(panel);
		frame.setSize(500,400);
		
		frame.setVisible(true);
		panel.start();
		
	}
	
	public boolean readSeqFasta() {
		
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fileBase+ ".tfa"));
			Vector names= new Vector(), seqs= new Vector();
		
			String line;
			while (buffy.ready()) {
			
				line= buffy.readLine().trim();
				if (line.length()< 1)
					continue;
				if (line.startsWith(">")) {
					names.add(line.substring(1,line.length()));
					seqs.add("");
				} else 
					seqs.setElementAt(((String) seqs.elementAt(seqs.size()- 1))+ line, seqs.size()- 1);
			}


				// init variables, sort seqs alphabeticall according to names
			baseNames= new String[names.size()];
			baseSequences= new String[names.size()];
			String[] sortedNames= new String[names.size()];
			int[] resort= new int[names.size()];	// resort[source]= destination
			for (int i= 0; i< names.size(); ++i)  {
				baseNames[i]= (String) names.elementAt(i);
				sortedNames[i]= (String) names.elementAt(i);
				baseSequences[i]= (String) seqs.elementAt(i);
			}
			Arrays.sort(sortedNames);
			for (int i= 0; i < baseNames.length; i++) {
				int x;
				for (x= 0; x < sortedNames.length; x++) 
					if (baseNames[i]== sortedNames[x])	// here OID ok..
						break;
				resort[i]= x;
			}

			for (int i= 0; i< names.size(); ++i) {
				if (resort[i]!= i) {
					
					String t= baseNames[i];					// swap names
					baseNames[i]= baseNames[resort[i]];
					baseNames[resort[i]]= t;
					
					t= baseSequences[i];				// swap seqs
					baseSequences[i]= baseSequences[resort[i]];
					baseSequences[resort[i]]= t;
					
					int tt= resort[i];					// swap array
					resort[i]= resort[tt]; 
					resort[tt]= tt;
					
					if (resort[i]!= i)
						--i;
				}
			}
			
		} catch (Exception e) {
			return false;
		}

		return true;
	}
	
	protected String writeOutTemp(String[] names, String[] seqs) {
	
			
		File tFile= null;
		try {
//			tFile= File.createTempFile("QAL", null);
			String tmpFName= "/homes/micha/";	// for sun
			if (OSChecker.isLinux())
				tmpFName="/home/ug/msammeth/";
			File parentDir= new File(tmpFName+ "scratch_delme");
			if (OSChecker.isWindows())
				tFile= File.createTempFile("QAL", null);
			else
				tFile= File.createTempFile("QAL", null, parentDir);	// SUN does not like var/tmp
		} catch (Exception e) {
			return null;
		}
		String tName= tFile.getAbsolutePath();
	
		FASTAWrapper fasta= new FASTAWrapper(tName);
/*		SequenceWrapper[] sw= new SequenceWrapper[names.length];
		for (int i= 0; i< sw.length; ++i) {
			sw[i].setName(names[i]);
			sw[i].setSequence(seqs[i]);
		}
		fasta.setWrappedSequences(sw);
*/	
		try {	
			fasta.writeFASTA(names, seqs);
		} catch (Exception e) {
			return null;
		}
		
		return tName;
	}
	
	/**
	 * Gets fragments of all pairs of sequences and delegates to 
	 * <code>filterPWOptimal(Vector)</code> for getting an optimal set
	 * of pw fragments.
	 *  
	 * @param frags
	 * @return
	 */
	protected MultiFrag[] filterPW(MultiFrag[] frags, boolean optimal, boolean output) {
		
			// create dynamic structure
		Vector allFrags= new Vector(frags.length);
		for (int i= 0; i < frags.length; i++) 
			allFrags.add(frags[i]);
				
			// filter for pw fragments
		Vector restVec= new Vector();	// pruned fragments
		int total= (baseSequences.length* (baseSequences.length- 1))/2;
		if (output)
			System.out.print("["+total+"]: ");
		int cnt= 0;
		for (int i= 0; i< baseSequences.length; ++i) {
			if (output)
				System.out.print(cnt+",");
			for (int j= (i+1); j< baseSequences.length; ++j) {
				
				cnt++;
				Vector pwFrags= new Vector();
				for (int k= 0; k< allFrags.size(); ++k) {
					MultiFrag frag= (MultiFrag) allFrags.elementAt(k); 
					if (((frag.getSequenceNo(true)== i)&& (frag.getSequenceNo(false)== j))
						|| ((frag.getSequenceNo(false)== j)&& (frag.getSequenceNo(true)== i))) {
						
						pwFrags.add(allFrags.remove(k--));	// temporarily remove from set
					}
				}
				
				if (optimal) {
					pwFrags= filterOptimalRepeatsPW(pwFrags, restVec);		// exhaustive search ..
				} else
					pwFrags= filterPWGreedy(pwFrags);		// .. or greedy selection
					
				allFrags.addAll(pwFrags);			// re-unite
			}
		}
			
			// re-convert to array
		rest= new MultiFrag[restVec.size()];
		for (int i= 0; i < rest.length; i++) 
			rest[i]= (MultiFrag) restVec.elementAt(i);
		MultiFrag[] result= new MultiFrag[allFrags.size()];
		for (int i= 0; i < result.length; i++) 
			result[i]= (MultiFrag) allFrags.elementAt(i);
		
		return result;

	}
	
	
	/**
	 * Finds for the fragments given an optimal consistent subset
	 * according to the rel. weights assigned.
	 * WARNING: changes <code>consistent</code> flag of fragments
	 * 
	 * @param someID
	 * @return
	 */
	protected Vector filterPWOptimal(Vector fragVec) {
		
			// submit only conflicting fragments (keep rest for solution)
		Vector conflictingFrags= new Vector();
		for (int i= 0; i < fragVec.size(); ++i) {
			
			boolean conflict= false;
			for (int j= (i+1); j< fragVec.size(); ++j) {
				if (((MultiFrag) fragVec.elementAt(i)).contradicts((MultiFrag) fragVec.elementAt(j))) {
					((MultiFrag) fragVec.elementAt(j)).setConsistent(false);
					conflict= true;		// there happened a conflict
				} else {
					((MultiFrag) fragVec.elementAt(j)).setConsistent(
						((MultiFrag) fragVec.elementAt(j)).isConsistent()&& true
					);					// mark fragment
				}
			}
				
			
			if (conflict) {  			// at least one conflict
				((MultiFrag) fragVec.elementAt(i)).setConsistent(false);
			} else {
				((MultiFrag) fragVec.elementAt(i)).setConsistent(
					((MultiFrag) fragVec.elementAt(i)).isConsistent()&& true
				);						// mark fragment
			}
		}
		for (int i= 0; i < fragVec.size(); i++) {
//			System.out.println(fragVec.elementAt(i));
			conflictingFrags.add(fragVec.remove(i--));
		}
		
			// recursionally generate all orders
		float maxScore= -1;
		Vector maxVec= null;
		Hashtable bHash= new Hashtable(conflictingFrags.size());
		for (int i= 0; i< conflictingFrags.size(); i++) {
			
			Vector tmpVec= (Vector) conflictingFrags.clone();	// new copy
			
				// prepare closure
			int[] lengths= new int[baseSequences.length];
			for (int j= 0; j< lengths.length; ++j)
				lengths[j]= baseSequences[j].length();
			Closure clos= AliGraphClosure.newAligGraphClosure(
				baseSequences.length, lengths, 0, null);
			
				// first one is always consistent
			Vector tmpIn= new Vector();		// for fragments aligned
			tmpIn.add(tmpVec.remove(i));
			MultiFrag tmpFrag= (MultiFrag) tmpIn.lastElement();
			float score= tmpFrag.getWeight();
			if (tmpIn.lastElement() instanceof MultiFragExt) {
		// 		deprecated
		//		AliGraphClosure.addAlignedFragment(clos, (MultiFragExt) tmpIn.lastElement());
				MultiFrag[] tmpFrags= ((MultiFragExt) tmpIn.lastElement()).getFragments();
				for (int j= 0; j < tmpFrags.length; j++) 
					AliGraphClosure.addAlignedFragment(clos, tmpFrags[j]);
			} else
				AliGraphClosure.addAlignedFragment(clos, (MultiFrag) tmpIn.lastElement());
			AliGraphClosure.computeClosure(clos);
			
				// recursion
			score+= filterPWRecursion(clos, tmpVec, tmpIn, score, bHash, 0);
			tmpVec.insertElementAt(tmpFrag, 0);
			if (score> maxScore) {			// TODO: for == count repeats used
				maxScore= score;
				maxVec= tmpVec;
			}
		}
		
			// prepare results
		fragVec.addAll(maxVec);								// assuming a non-null result
		for (int i= 0; i < fragVec.size(); i++) 			// mark as solution
			((MultiFrag) fragVec.elementAt(i)).setConsistent(true);
		for (int i= 0; i < conflictingFrags.size(); i++) {	// find rest
			int j;
			for (j= 0; j< maxVec.size(); ++j)				// TODO: make mor efficient, binary search..
				if (conflictingFrags.elementAt(i).equals(maxVec.elementAt(j)))
					break;
			if (j>= maxVec.size()) {
				((MultiFrag) conflictingFrags.elementAt(i)).setConsistent(false);
				fragVec.add(conflictingFrags.elementAt(i));	// mark as not in opt. solution
			}
		}
		
		return fragVec;
	}
	
	protected Vector filterPWGreedy(Vector fragVec) {

			// sort frags
		MultiFrag[] sortedFrags= new MultiFrag[fragVec.size()];
		for (int i= 0; i < sortedFrags.length; i++) 
			sortedFrags[i]= (MultiFrag) fragVec.elementAt(i);
		sortFragmentsWGT(sortedFrags);
		fragVec= new Vector();
		for (int i= 0; i < sortedFrags.length; i++) 
			fragVec.add(sortedFrags[i]);
					
			// prepare closure
		int[] lengths= new int[baseSequences.length];
		for (int j= 0; j< lengths.length; ++j)
			lengths[j]= baseSequences[j].length();
		Closure clos= AliGraphClosure.newAligGraphClosure(
			baseSequences.length, lengths, 0, null);

			// greedily search for next-best consistent fragment and add
		for (int i= 0; i< fragVec.size(); i++) {
			
			MultiFrag tmpFrag= (MultiFrag) fragVec.elementAt(i);
			
				// check consistency
			Vector allFrags= new Vector();
			if (tmpFrag instanceof MultiFragExt) {
				MultiFrag[] tmpFrags= ((MultiFragExt) tmpFrag).getFragments();
				boolean alignable= true;
				for (int j= 0; alignable&& j < tmpFrags.length; j++) 
					alignable&= AliGraphClosure.alignableFragment(clos, tmpFrags[j]);
				if (!alignable) {
//				if (!AliGraphClosure.alignableFragment(clos, (MultiFragExt) tmpFrag)) {
					tmpFrag.setConsistent(false);
					continue;
				}
			} else if (!AliGraphClosure.alignableFragment(clos, tmpFrag)) {
				tmpFrag.setConsistent(false);
				continue;
			}

				// if alignment is possible				
			tmpFrag.setConsistent(true);
			if (tmpFrag instanceof MultiFragExt) {
// deprecated
//				AliGraphClosure.addAlignedFragment(clos, (MultiFragExt) tmpFrag);
				MultiFrag[] tmpFrags= ((MultiFragExt) tmpFrag).getFragments();
				for (int j= 0; j < tmpFrags.length; j++) 
					AliGraphClosure.addAlignedFragment(clos, tmpFrags[j]);
			} else
				AliGraphClosure.addAlignedFragment(clos, (MultiFrag) tmpFrag);
			AliGraphClosure.computeClosure(clos);
		}
		
		return fragVec;		// result in consistent flags
	}
	
	
	/**
	 * Recursionally tries all combinations and returns best score and vector 
	 * from its sub-tree.
	 * 
	 * @param clos
	 * @param tmpVec
	 * @param inVec
	 * @return
	 */
	private float filterPWRecursion(Closure origClos, Vector origVec, Vector inVec, float score, Hashtable bHash, int depth) {
		
			// abort 
		if (origVec== null)
			return 0f;
				
			// filter off no longer alignable fragments
		for (int i= 0; i< origVec.size(); i++) {
			
				// check revelance
				// remove not alignable fragments from all subtrees
			MultiFrag tmpFrag= (MultiFrag) origVec.elementAt(i);
			if (tmpFrag instanceof MultiFragExt) {
// deprecated
//				if (!AliGraphClosure.alignableFragment(origClos, (MultiFragExt) tmpFrag))
				boolean consistent= true;
				MultiFrag[] tmpFrags= ((MultiFragExt) tmpFrag).getFragments();
				for (int j= 0; j < tmpFrags.length; j++) 
					consistent&= AliGraphClosure.alignableFragment(origClos, tmpFrags[j]);
				if (!consistent)
					origVec.remove(i--);
			} else if (!AliGraphClosure.alignableFragment(origClos, tmpFrag))
					origVec.remove(i--);
		}		
			// abort condition
		if (origVec.size()< 1)
			return 0f;
			
			// recursionally generate all orders 
		float maxScore= -1;
		Vector maxVec= null;
		for (int i= 0; i< origVec.size(); i++) {
			
			MultiFrag tmpFrag= (MultiFrag) origVec.elementAt(i);
			float tmpScore= 0; 
			
				// init
			Vector tmpVec= (Vector) origVec.clone();		// new copy
			Vector tmpIn= new Vector(inVec);					// new vector for results 
			tmpIn.add(tmpVec.remove(i));
			ResultHash result= filterGetHash(bHash, tmpIn);

				// compute recursionally solution
			if (result== null || result.getScore()< 0) {
				Closure clos= (Closure) origClos.clone();	
				
				if (tmpFrag instanceof MultiFragExt) {
//					deprecated
//					AliGraphClosure.addAlignedFragment(clos, (MultiFragExt) tmpFrag);
					MultiFrag[] tmpFrags= ((MultiFragExt) tmpFrag).getFragments();
					for (int j= 0; j < tmpFrags.length; j++) 
						AliGraphClosure.addAlignedFragment(clos, tmpFrags[j]);
				} else
					AliGraphClosure.addAlignedFragment(clos, tmpFrag);
				AliGraphClosure.computeClosure(clos);
			
					// recursion
				tmpScore+= filterPWRecursion(clos, tmpVec, tmpIn, score, bHash, depth+ 1);
				filterPutHash(bHash, tmpIn, tmpVec, tmpScore);
				tmpScore+= tmpFrag.getWeight();

				// or lookup in hash of hashes
			} else {
				tmpScore= result.getScore();
				tmpVec= result.getResult(); 
			}
			
				// new maximum
			if (tmpScore> maxScore) {				// TODO: for == count repeats used				
				maxScore= tmpScore;
				maxVec= (Vector) tmpVec.clone();	// there the result comes back, clone for not change hash
				maxVec.insertElementAt(tmpFrag, 0);	// now write to result
			} 
		}
		
			// build score and result vector from depth to top
		origVec.removeAllElements();
		if (maxVec!= null)							// if all inconsistent
			origVec.addAll(maxVec);
		return maxScore;
	}
	
	/**
	 * Writes in the hash of hashes to the ascending order of the given path
	 * the given solution.
	 * 
	 * @param bHash	current hash of the recursional step
	 * @param path current path to be permutated for writing
	 * @param solution the solution to be written
	 */
	private void filterPutHash(Hashtable bHash, Vector path, Vector solution, float score) {
		
		if (path== null)
			return;
			
			// sort for putting at defined position (instead of mirroring) 
		int[] ids= new int[path.size()];
		for (int i= 0; i< ids.length; ++i)
			ids[i]= ((MultiFrag) path.elementAt(i)).getNumber();
		Arrays.sort(ids);

			// iterate hash of hashes to get last hashtable for entry
		for (int i= 0; i < ids.length; ++i) {
			
			Integer key= new Integer(ids[i]);
			if (bHash.get(key)== null) {
				ResultHash tmpHash= new ResultHash();
				bHash.put(key, tmpHash);
				bHash= tmpHash;
			} else
				bHash= (Hashtable) bHash.get(key);
		}
		
			// put in the result (and additional info)
		ResultHash rHash= (ResultHash) bHash;
		rHash.setResult(solution);
		rHash.setScore(score);
		rHash.setPath(path);
	}
	
	/**
	 * Retrieves a result from the hash of hashes, if there has one be saved
	 * 
	 * @param bHash	current hash of the recursional step
	 * @param path current path to searched
	 * @return the <code>Vector</code> with the solution or <code>null</code>
	 */	 
	private ResultHash filterGetHash(Hashtable bHash, Vector path) {
			
			// abort
		if ((path== null)|| (bHash== null))
			return null;
		
			// sort for getting from defined position
		int[] ids= new int[path.size()];
		for (int i= 0; i< ids.length; ++i)
			ids[i]= ((MultiFrag) path.elementAt(i)).getNumber();
		Arrays.sort(ids);
			
			// iterate hash of hashes for correct hashtable
		for (int i= 0; i < ids.length; i++) {
			
			if (bHash== null)
				return null;
			Integer key= new Integer(ids[i]);
			bHash= (Hashtable) bHash.get(key);
		}

			// return hashtable or null if not found
		if (bHash== null)
			return null;
		else
			return (ResultHash) bHash;	
	}

	
	public int getSeqNo(String someID) {
			
		for (int i= 0; i < baseNames.length; i++) 
			if (baseNames[i].equals(someID))
				return i;
				
		return -1;
	}
	
	public String getSequence(String someID) {
		
		for (int i= 0; i < baseNames.length; i++) 
			if (baseNames[i].equals(someID))
				return baseSequences[i];
		
		return null;
	}
	
	
	/**
	 * works now not on baseSeq but on subseq from diagonals
	 * (for jumping frags)
	 * 
	 * see: 
	 * ow_add in functions.c
	 * Morgenstern, Dress, Werner: Multiple DNA and protein sequence alignment 
	 * based on segment-to-segment comparison.
	 * 
	 * WARNING: not consistent with c-code!
	 * 2-decimal floats (c) instead of 6-decimal (c)
	 *  
	 * @param Dl
	 * @param Dm
	 * 
	 * @deprecated Not capable of comparing normal fragments to ones with jpoints
	 * @see calcOlWegt_old()
	 */
	public strictfp void addOW_old(MultiFragExt Dl, MultiFragExt Dm) {
		
		if ((Dl.getNumber()== 124&& Dm.getNumber()== 37)||
			(Dl.getNumber()== 37&& Dm.getNumber()== 124))
			System.currentTimeMillis();
			
		int[] mapping= Constants.getMapBChars();
		for (int i= 0; i< 2; ++i) 
			for (int j= 0; j< 2; ++j)
				if ((Dl.getSequenceNos()[i]== Dm.getSequenceNos()[j])	// common sequence (sj)
					&& (Dl.getSequenceNos()[j]!= Dm.getSequenceNos()[i])// different sequences (si,sk)
						// start before end of other diagonal on common sequence
						// (real length matters for comparing normal frag to one w jpoints, e.g. inconsistent IVs) 
					&& ((Dl.getSequenceStarts()[i]< Dm.getSequenceStarts()[j]+ Dm.getRealLength(j))
						&& (Dm.getSequenceStarts()[j]< Dl.getSequenceStarts()[i]+ Dl.getRealLength(i)))
					) {
						
						
							// length of overlap (real length does matter !)
						int conslen= Math.min(Dl.getSequenceStarts()[i]+ Dl.getRealLength(i), 
												Dm.getSequenceStarts()[j]+ Dm.getRealLength(j))
									 - Math.max(Dl.getSequenceStarts()[i], Dm.getSequenceStarts()[j]);
						
							// calc start points of Dn
						int s1= Dl.getSequenceNos()[(i+1)%2];	// = si
						int s2= Dm.getSequenceNos()[(j+1)%2];	// = sk
/*						boolean micha_dbg= (Dl.getWeight()>60&& Dl.getWeight()< 61);
						if (micha_dbg)
							System.out.println(
								s1+","+Dl.getSequenceNos()[i]+","+s2+":"+
								Dl.getWeight()+","+Dm.getWeight()
							);
*/						
//						int b1= Dl.getSequenceStarts()[(i+1)%2];
						int b1= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
						int diff= Dm.getSequenceStarts()[j]- Dl.getSequenceStarts()[i];
						if (diff> 0)
							b1+= diff;
						
//						int b2= Dm.getSequenceStarts()[(j+1)% 2];
						int b2= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
						diff= Dl.getSequenceStarts()[i]- Dm.getSequenceStarts()[j];
						if (diff> 0)
							b2+= diff;

						// ?!! change here to sum(all little diagonals) !!?
							
							// calc s_Dn
						int match= 0;
//						String s1Up= baseSequences[s1].toUpperCase();
//						String s2Up= baseSequences[s2].toUpperCase();
							// now we work on the subsequences ...
						String s1Up= Dl.getSubsequence(baseSequences[s1],((i+1)%2)).toUpperCase();
						String s2Up= Dm.getSubsequence(baseSequences[s2],((j+1)%2)).toUpperCase();
						for (int k= 0; k < conslen; ++k) 
							match+= Constants.getBlosum()
								[mapping[s1Up.charAt(b1+ k)- 65]]
								[mapping[s2Up.charAt(b2+ k)- 65]];
//						if (micha_dbg)
//							System.out.println("\tsim= "+ match+", len= "+conslen+", b1= "+b1+", b2= "+b2);
						
							// calc w(Dn)
						double addWgt= 					// lookup
							Constants.getTpProt400()[conslen][match];
//						System.out.print(" -> "+ addWgt);

						double factor= ((double)
						(baseSequences[s1].length()
						* baseSequences[s2].length()))
//						((double) concatSeqs[s1].length()
//						* (double) concatSeqs[s2].length()))
							/ 400d;
						if (addWgt> 0.0000000001)		// calc
							addWgt= 1d- Math.pow((1d- addWgt), factor);
						else
							addWgt= addWgt* factor;
//						if (micha_dbg) 
//							System.out.println("\tfactor= "+factor+"("+concatSeqs[s1].length()+","+concatSeqs[s2].length()+")\n\tcorrected: "+addWgt);
		
						float relWgt= 					// compute rel. Weight
							(-1f)* ((float) Math.log(addWgt)); 
//						relWgt= Math.round(relWgt* 100f)/ 100f;
//						if (micha_dbg)
//							System.out.println("\trelWgt= "+relWgt);	
						
							// add to olw
								// !!! not very nice :(
						if (relWgt== Float.POSITIVE_INFINITY) { 
							System.err.println("overlap: infinity!");
							continue;
						}
//						if (micha_dbg)
//							System.out.print("\t"+ Dl.getOverlapWeight()+ "+"+relWgt+"= ");
						Dl.setOverlapWeight(Dl.getOverlapWeight()+ (float) relWgt);
//						if (micha_dbg)
//							System.out.println(Dl.getOverlapWeight());
//						if (micha_dbg)
//							System.out.print("\t"+ Dm.getOverlapWeight()+"+ "+relWgt+"= ");
						Dm.setOverlapWeight(Dm.getOverlapWeight()+ (float) relWgt);
//						if (micha_dbg)
//							System.out.println(Dm.getOverlapWeight());
					}
	}
	
    /**
	 * works now on realLengths / assumes full overlap of indel frags
	 * (for indel frags)
	 * 
	 * see: 
	 * ow_add in functions.c
	 * Morgenstern, Dress, Werner: Multiple DNA and protein sequence alignment 
	 * based on segment-to-segment comparison.
	 * 
	 * REMARK: not consistent with c-code (2-decimal floats instead of 6-decimal (c))
	 * REM2: computes full overlap if gaps in common sequence
	 * 		e.g.: ~~~~~~~~~~~~~~~~~~
	 * 			  ~~~~~~-~~~~---~~~~
	 * 			  ~~~~~~~~~~~~~~~~~~
	 * -->
	 * 			  ~~~~~~~~~~~~~~~~~~
	 * 			  ~~~~~~~~~~~~~~~~~~
	 * 
	 * in contrast to
	 * 
	 * 			  ~~~~~   ~~~~   ~~~~
	 *            ~~~~~ + ~~~~ + ~~~~
	 *  
	 * @param Dl
	 * @param Dm
	 */
	public strictfp void addOW_indel_old(MultiFragExt Dl, MultiFragExt Dm) {
		
		int[] mapping= Constants.getMapBChars();
		float oldOLW= Dl.getOverlapWeight();
		for (int i= 0; i< 2; ++i) 
			for (int j= 0; j< 2; ++j)
				if ((Dl.getSequenceNos()[i]== Dm.getSequenceNos()[j])	// common sequence (sj)
					&& (Dl.getSequenceNos()[j]!= Dm.getSequenceNos()[i])// different sequences (si,sk)
						// start before end of other diagonal on common sequence
						// (!! real length !!) 
					&& ((Dl.getSequenceStarts()[i]< Dm.getSequenceStarts()[j]+ Dm.getRealLength(j))
						&& (Dm.getSequenceStarts()[j]< Dl.getSequenceStarts()[i]+ Dl.getRealLength(i)))
					) {
						
							// subsequences (maybe stringbuffers ?!..)
						int si= Dl.getSequenceNos()[(i+1)%2];	// = si
						int sj= Dl.getSequenceNos()[i]; 		// = Dm.getSequenceNos()[j] 
						int sk= Dm.getSequenceNos()[(j+1)%2];	// = sk
						String useq1= Dl.getSubsequence(baseSequences[si], j).toUpperCase();	
						String cseq1= Dl.getSubsequence(baseSequences[sj], i).toUpperCase();
						String useq2= Dm.getSubsequence(baseSequences[sk], i).toUpperCase();
						String cseq2= Dm.getSubsequence(baseSequences[sj], j).toUpperCase();
						
							// make common sequences "equal" (ie. introduce gaps where not corresponding)
							// [bit stupid, optimize]
						String lcommon= "", scommon= "", lunique="", sunique= "";
						if (cseq1.length()< cseq2.length()) {	// start with shorter one
							scommon= cseq1; lcommon= cseq2; sunique= useq1; lunique= useq2; 
						} else {
							scommon= cseq2; lcommon= cseq1; sunique= useq2; lunique= useq1;
						}
						
						for (int k= 0; k < scommon.length(); k++) 
							if (scommon.charAt(k)!= lcommon.charAt(k))
								if (MultipleAlignmentModel.isGapChar(lcommon.charAt(k))) {
									scommon= scommon.substring(0,k)+ "-"+ scommon.substring(k,scommon.length());
									sunique= sunique.substring(0,k)+ "-"+ sunique.substring(k,sunique.length());
								} else { // gap char in cseq1 (another if ?!!)
									lcommon= lcommon.substring(0,k)+ "-"+ lcommon.substring(k,lcommon.length());
									lunique= lunique.substring(0,k)+ "-"+ lunique.substring(k,lunique.length());
								}
								// evtl. trailing gaps
						for (int k= scommon.length(); k < lcommon.length(); k++) 
							if (MultipleAlignmentModel.isGapChar(lcommon.charAt(k))) {
								scommon+= "-";
								sunique+= "-";
							} // cannot be else, other sequence shorter

							// chk
						if ((scommon.length()!= lcommon.length())
							|| (scommon.length()!= sunique.length())
							|| (lcommon.length()!= lunique.length()))	// all 4 must be equal
							System.err.println("shit. consensus not equal.");

							// common attributes of all sub-diagonals
						double factor= ((double)
										(baseSequences[si].length()
										* baseSequences[sk].length()))
										/ 400d;

							// find single sub-diagonals and add up weight
						int pos= 0, lpos= 0; 
						while (pos< scommon.length()) {
							for (;(MultipleAlignmentModel.isGapChar(sunique.charAt(pos))
									|| MultipleAlignmentModel.isGapChar(lunique.charAt(pos)))
									&& (pos< sunique.length());++pos)
									; // find new "continue"
							lpos= pos;
							for (;((pos< sunique.length())
										&& (pos- lpos< 40))	// max diagonal length (lookup table)
									&& (!(MultipleAlignmentModel.isGapChar(sunique.charAt(pos))
										|| MultipleAlignmentModel.isGapChar(lunique.charAt(pos))));
									++pos)
									; // find new "breakpoint"
							
							String subfrag1= sunique.substring(lpos, pos);
							String subfrag2= lunique.substring(lpos, pos);
							 
							int conslen= subfrag1.length();	// consensus length
							
								// calc s_Dn
							int match= 0;					// no more gaps in diagonal
							for (int k= 0; k< conslen; k++) {
								match+= Constants.getBlosum()
									[mapping[subfrag1.charAt(k)- 65]]
									[mapping[subfrag2.charAt(k)- 65]];
							}
//							if (micha_dbg)
//								System.out.println("\tsim= "+ match+", len= "+conslen+", b1= "+b1+", b2= "+b2);
							
							double addWgt= 												// lookup
								Constants.getTpProt400()[conslen][match];		// substract gaps !!
//							System.out.print(" -> "+ addWgt);
							if (addWgt> 0.0000000001)		// calc
								addWgt= 1d- Math.pow((1d- addWgt), factor);
							else
								addWgt= addWgt* factor;
//							if (micha_dbg) 
//								System.out.println("\tfactor= "+factor+"("+concatSeqs[s1].length()+","+concatSeqs[s2].length()+")\n\tcorrected: "+addWgt);
							

							float relWgt= 					// compute rel. Weight
								(-1f)* ((float) Math.log(addWgt)); 
//							relWgt= Math.round(relWgt* 100f)/ 100f;
//							if (micha_dbg)
//								System.out.println("\trelWgt= "+relWgt);	
						
								// add to olw
							// !!! not very nice :(
							if (relWgt== Float.POSITIVE_INFINITY) {
								System.err.println("overlap2: infinity2 "+conslen+" for "+ Dl.getNumber()+"x"+Dm.getNumber());
								continue;	// skip if no entry in lookup-table (=0.0 -> log(0.0)= posInfinity?!!) 
							}
//							if (micha_dbg)
//								System.out.print("\t"+ Dl.getOverlapWeight()+ "+"+relWgt+"= ");
							Dl.setOverlapWeight(Dl.getOverlapWeight()+ (float) relWgt);
//							if (micha_dbg)
//								System.out.println(Dl.getOverlapWeight());
//							if (micha_dbg)
//								System.out.print("\t"+ Dm.getOverlapWeight()+"+ "+relWgt+"= ");
							Dm.setOverlapWeight(Dm.getOverlapWeight()+ (float) relWgt);
//							if (micha_dbg)
//								System.out.println(Dm.getOverlapWeight());

						}
						
							// calc w_Dn
							// ?!! change here to sum(all little diagonals) !!?
					}
		float newOLW= Dl.getOverlapWeight()- oldOLW;
		System.currentTimeMillis();
	}	
	
	/**
	 * ATTENTION: Weights have to be already been computed!
	 * @param frags
	 */
	public void calcOLW(MultiFrag[] frags) {
	
			// init olw with wgt
		for (int i= 0; i < frags.length; i++) 
			frags[i].setOverlapWeight(frags[i].getWeight());
			
			// add up olw
		System.out.println("calculating OLW ["+frags.length+"]:");
		System.out.print("\t");
		for (int i= 0; i < frags.length; i++) {
			if ((i< 10)
				|| (i< 100 && i%10== 0) 
				|| (i< 1000 && i%100== 0) 
				|| (i< 10000 && i%1000== 0) 
				|| (i< 100000 && i%10000== 0)) 
				System.out.print(i+",");
			for (int j= (i+1); j < frags.length; j++) 
				calcOLW(frags[i], frags[j]);
		}
		System.out.println("done.");
	}
		/**
	 * ATTENTION: Weights have to be already been computed!
	 * @param frags
	 */
	public void calcOLW_old(MultiFrag[] frags) {
	
			// init olw with wgt
		for (int i= 0; i < frags.length; i++) 
			frags[i].setOverlapWeight(frags[i].getWeight());
			
			// add up olw
		System.out.println("calculating OLW ["+frags.length+"]:");
		System.out.print("\t");
		for (int i= 0; i < frags.length; i++) {
			if ((i< 10)
				|| (i< 100 && i%10== 0) 
				|| (i< 1000 && i%100== 0) 
				|| (i< 10000 && i%1000== 0) 
				|| (i< 100000 && i%10000== 0)) 
				System.out.print(i+",");
			for (int j= (i+1); j < frags.length; j++) 
				calcOLW_old(frags[i], frags[j]);
		}
		System.out.println("done.");
	}

	
	
	/**
	 * returns frags of best solution
	 */
	Vector filterOptimalRepeatsPW(Vector fragVec, Vector pruned) {
		
		
			// sort according to weight
		MultiFrag[] pwFrags= new MultiFrag[fragVec.size()];
		for (int i= 0; i < pwFrags.length; i++) 
			pwFrags[i]= (MultiFrag) fragVec.elementAt(i);
		sortFragmentsWGT(pwFrags);
			
			// construct greedy solution such that all have been tried
		Vector tryVec= new Vector();
		Vector solution, rest, nbAdded, bestSolution= null, bestRest= null;
		double value= 0d, bestValue= -1d;
		Closure clos;
		int[] lengths= new int[baseSequences.length];
		for (int j= 0; j< lengths.length; ++j)
			lengths[j]= baseSequences[j].length();
		
		boolean run= true;
		while(run|| tryVec.size()!= 0) {

			solution= new Vector();		
			rest= new Vector();
			nbAdded= new Vector();
			clos= AliGraphClosure.newAligGraphClosure(
				baseSequences.length, lengths, 0, null);

			if (!run) {
				solution.add(tryVec.remove(0));		// try a new fragment
				value= ((MultiFrag) tryVec.elementAt(0)).getWeight();
			}
		
			for (int i = 0; i < pwFrags.length; i++) {	// greedily add
				if (AliGraphClosure.alignableFragment(clos, pwFrags[i])) {
					pwFrags[i].setAccepted(true);
					solution.add(pwFrags[i]);
					int j;
					for (j = 0; j < nbAdded.size(); j++) 
						if (pwFrags[i].getNumber()==
							((Integer) nbAdded.elementAt(j)).intValue())
							break;						// see if superfrag already added
					if (j== nbAdded.size())
						value+= pwFrags[i].getWeight();		// not yet added
					
				} else {
					pwFrags[i].setAccepted(false);
					rest.add(pwFrags[i]);					
				}
			}
			
			if (value> bestValue) {		// compare to optimum
				bestSolution= solution;
				bestRest= rest;
				bestValue= value;
			}
			
			if (!run) {				// compare with still untried frags
				for (int i = 0; i < tryVec.size(); i++) 
					for (int j = 0; j < solution.size(); j++) 
						if (tryVec.elementAt(i).equals(solution.elementAt(i)))
							tryVec.remove(i);	// remove frags already included in an greedy solution
			} else
				tryVec= (Vector) rest.clone();
			
			run= false; 					// force 1st round
		}
		
			// result
		pruned= bestRest;	// return via parameter
		for (int i = 0; i < bestSolution.size(); i++) 			// remove double added
			for (int j = (i+1); j < bestSolution.size(); j++) 
				if (bestSolution.elementAt(i).equals(bestSolution.elementAt(j)))
					bestSolution.remove(i--);
		return bestSolution;
	}

	/**
	 * Calculates overlap weights for two given fragments. The olw
	 * no longer is symetrical and is only added if the respective
	 * fragment which is compared has been marked <code>consistent
	 * </code>. Therefore it is calculated for:<br>
	 * 
	 * frag1 frag2 	(frag1+olw) (frag2+olw)
	 * cons. cons.	yes			yes
	 * cons. incon. no			yes
	 * incon.cons.	yes			no
	 * incon.incon. no			no
	 * 
	 * REMARK: adds up olw of subfrags if same gaps in common seq (s. addOW_old)
	 * TODO check if addOW_indel_old is better
	 * 
	 * @param frag1 
	 * @param frag2  
	 */
	public void calcOLW(MultiFrag frag1, MultiFrag frag2) {
		
			// check
		if ((frag1== null)|| (frag2== null)
			|| ((!frag1.isConsistent())&& (!frag2.isConsistent()))	// if one inconsistent -> upweight donor
		)
			return;

		boolean doit= false;
		for (int i= 0; i< 2; ++i) { 
			for (int j= 0; j< 2; ++j) {
				if ((frag1.getSequenceNos()[i]== frag2.getSequenceNos()[j])	// common sequence (si,sj) f.a. combis (i,j)
						// start before end of other diagonal on common sequence
						// (real length doesn't matter for overlap in same sequence) 
					&& ((frag1.getSequenceStarts()[i]< frag2.getSequenceStarts()[j]+ frag2.getLength())
						&& (frag2.getSequenceStarts()[j]< frag1.getSequenceStarts()[i]+ frag1.getLength()))
					) { 
						doit= true;
						break;
					}
			}
			if (doit)
				break;
		}
		if (!doit)
			return;
		
			// decompose 
		MultiFrag[] decompFrags1= null;
		MultiFrag[] decompFrags2= null;
		if (frag1 instanceof MultiFragExt)
			decompFrags1= ((MultiFragExt) frag1).getFragments();
		if (frag2 instanceof MultiFragExt)
			decompFrags2= ((MultiFragExt) frag2).getFragments();
		
			// iterate sub-frags and add olw
//		if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) {
//			System.out.println(frag1.getSubsequence(baseSequences[frag1.getSequenceNo(true)], true));
//			System.out.println(frag1.getSubsequence(baseSequences[frag1.getSequenceNo(false)], false));
//			System.out.println(frag2.getSubsequence(baseSequences[frag2.getSequenceNo(true)], true));
//			System.out.println(frag2.getSubsequence(baseSequences[frag2.getSequenceNo(false)], false));
//		}
		float olw= 0f;
		for (int i= 0; i < decompFrags1.length; i++) 
			for (int j= 0; j < decompFrags2.length; j++) { 
				float olw_up= getOLW(decompFrags1[i], decompFrags2[j]); 
//				if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) {
//					System.out.println(decompFrags1[i].getSubsequence(baseSequences[decompFrags1[i].getSequenceNo(true)], true));
//					System.out.println(decompFrags1[i].getSubsequence(baseSequences[decompFrags1[i].getSequenceNo(false)], false));
//					System.out.println(decompFrags2[j].getSubsequence(baseSequences[decompFrags2[j].getSequenceNo(true)], true));
//					System.out.println(decompFrags2[j].getSubsequence(baseSequences[decompFrags2[j].getSequenceNo(false)], false));
//					System.out.println("==> "+ olw_up);
//				}
				olw+= olw_up;
//				if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) 
//					System.out.println();
			}
		
			// add olw to frags		
		if (frag1.isConsistent())										// if the other is consistent...
			frag2.setOverlapWeight(frag2.getOverlapWeight()+ olw);
		if (frag2.isConsistent())										// ... it can give an upweight
			frag1.setOverlapWeight(frag1.getOverlapWeight()+ olw);
	}
	
public void calcOLW_old(MultiFrag frag1, MultiFrag frag2) {
		
		// check
	if ((frag1== null)|| (frag2== null)
		|| ((!frag1.isConsistent())&& (!frag2.isConsistent()))	// if one inconsistent -> upweight donor
	)
		return;

		// decompose 
	MultiFrag[] decompFrags1= null;
	MultiFrag[] decompFrags2= null;
	if (frag1 instanceof MultiFragExt)
		decompFrags1= ((MultiFragExt) frag1).getFragments();
	if (frag2 instanceof MultiFragExt)
		decompFrags2= ((MultiFragExt) frag2).getFragments();
		
		// iterate sub-frags and add olw
//	if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) {
//		System.out.println(frag1.getSubsequence(baseSequences[frag1.getSequenceNo(true)], true));
//		System.out.println(frag1.getSubsequence(baseSequences[frag1.getSequenceNo(false)], false));
//		System.out.println(frag2.getSubsequence(baseSequences[frag2.getSequenceNo(true)], true));
//		System.out.println(frag2.getSubsequence(baseSequences[frag2.getSequenceNo(false)], false));
//	}
	float olw= 0f;
	for (int i= 0; i < decompFrags1.length; i++) 
		for (int j= 0; j < decompFrags2.length; j++) { 
			float olw_up= getOLW_old(decompFrags1[i], decompFrags2[j]); 
//			if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) {
//				System.out.println(decompFrags1[i].getSubsequence(baseSequences[decompFrags1[i].getSequenceNo(true)], true));
//				System.out.println(decompFrags1[i].getSubsequence(baseSequences[decompFrags1[i].getSequenceNo(false)], false));
//				System.out.println(decompFrags2[j].getSubsequence(baseSequences[decompFrags2[j].getSequenceNo(true)], true));
//				System.out.println(decompFrags2[j].getSubsequence(baseSequences[decompFrags2[j].getSequenceNo(false)], false));
//				System.out.println("==> "+ olw_up);
//			}
			olw+= olw_up;
//			if ((frag1.getNumber()== 94/*|| frag1.getNumber()== 85*/)&& frag2.getNumber()== 149) 
//				System.out.println();
		}
		
		// add olw to frags		
	if (frag1.isConsistent())										// if the other is consistent...
		frag2.setOverlapWeight(frag2.getOverlapWeight()+ olw);
	if (frag2.isConsistent())										// ... it can give an upweight
		frag1.setOverlapWeight(frag1.getOverlapWeight()+ olw);
}	

	/**
	 * Calculates and returns the weight of the overlapping area of both diagonals (s_Dn). 
	 * Note that only elementary fragments are input.
	 * 
	 * @param Dl
	 * @param Dm
	 * @return
	 */
	public float getOLW(MultiFrag Dl, MultiFrag Dm) {
		
		int[] mapping= Constants.getMapBChars();
		boolean tryCombi= true;
		float relWgt= 0f;
		for (int i= 0; tryCombi&& (i< 2); ++i) 
			for (int j= 0; tryCombi&& (j< 2); ++j)
				if (
//					(Dl.getSequenceNos()[i]== Dm.getSequenceNos()[j])	// common sequence (sj)
//					&& (Dl.getSequenceNos()[j]!= Dm.getSequenceNos()[i])// different sequences (si,sk)
//					&&
						// start before end of other diagonal on common sequence
						// (real length doesn't matter for overlap in same sequence)
					((Dl.getSequenceStarts()[i]< Dm.getSequenceStarts()[j]+ Dm.getLength())
						&& (Dm.getSequenceStarts()[j]< Dl.getSequenceStarts()[i]+ Dl.getLength()))
					) {
						
						tryCombi= false;
						
							// length of overlap (real length doesnt matter, same jpoints)
						int conslen= Math.min(Dl.getSequenceStarts()[i]+ Dl.getLength(), 
												Dm.getSequenceStarts()[j]+ Dm.getLength())
									 - Math.max(Dl.getSequenceStarts()[i], Dm.getSequenceStarts()[j]);

							// calc start points of Dn
						int s1= Dl.getSequenceNos()[(i+1)%2];	// = si
						int s2= Dm.getSequenceNos()[(j+1)%2];	// = sk

						int b1= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
						int diff= Dm.getSequenceStarts()[j]- Dl.getSequenceStarts()[i];
						if (diff> 0)
							b1+= diff;
						
						int b2= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
						diff= Dl.getSequenceStarts()[i]- Dm.getSequenceStarts()[j];
						if (diff> 0)
							b2+= diff;

							// calc s_Dn
						String s1Up=			// the sub-sequences 
							Dl.getSubsequence(baseSequences[s1],((i+1)%2)).toUpperCase();
						String s2Up= 
							Dm.getSubsequence(baseSequences[s2],((j+1)%2)).toUpperCase();
						int x= (conslen/ 40);	// tokenize to 40mers (lookup-table!)
						if ((conslen% 40)!= 0)
							++x;
						for (int l= 0; l< x; ++l) {
				
							int match= 0;
							int k;
							for (k= (l* 40); (k < conslen)&& (k< l*40+40); ++k) 
								try {
									match+= Constants.getBlosum()
										[mapping[s1Up.charAt(b1+ k)- 65]]
										[mapping[s2Up.charAt(b2+ k)- 65]];
								} catch (ArrayIndexOutOfBoundsException e) {
									System.err.println("Not supported char in {"+ 
										s1Up.charAt(b1+ k)+ ","+ s2Up.charAt(b2+ k)+ "}, skipped.");
								}
							
								// calc w(Dn)
							double addWgt= 					// lookup
								Constants.getTpProt400()[k- (l* 40)][match];
	
							double factor= ((double)
							(baseSequences[s1].length()
							* baseSequences[s2].length()))
								/ 400d;
							if (addWgt> 0.0000000001)		// calc
								addWgt= 1d- Math.pow((1d- addWgt), factor);
							else
								addWgt= addWgt* factor;
			
							addWgt= (-1f)* ((float) Math.log(addWgt)); 
							relWgt+= addWgt;				// compute rel. Weight

						}
					}
			
			return relWgt;
	}
	
public float getOLW_old(MultiFrag Dl, MultiFrag Dm) {
		
	int[] mapping= Constants.getMapBChars();
	boolean tryCombi= true;
	float relWgt= 0f;
	for (int i= 0; tryCombi&& (i< 2); ++i) 
		for (int j= 0; tryCombi&& (j< 2); ++j)
			if (
				(Dl.getSequenceNos()[i]== Dm.getSequenceNos()[j])	// common sequence (sj)
				&& (Dl.getSequenceNos()[j]!= Dm.getSequenceNos()[i])// different sequences (si,sk)
				&&
					// start before end of other diagonal on common sequence
					// (real length doesn't matter for overlap in same sequence)
				((Dl.getSequenceStarts()[i]< Dm.getSequenceStarts()[j]+ Dm.getLength())
					&& (Dm.getSequenceStarts()[j]< Dl.getSequenceStarts()[i]+ Dl.getLength()))
				) {
						
					tryCombi= false;
						
						// length of overlap (real length doesnt matter, same jpoints)
					int conslen= Math.min(Dl.getSequenceStarts()[i]+ Dl.getLength(), 
											Dm.getSequenceStarts()[j]+ Dm.getLength())
								 - Math.max(Dl.getSequenceStarts()[i], Dm.getSequenceStarts()[j]);

						// calc start points of Dn
					int s1= Dl.getSequenceNos()[(i+1)%2];	// = si
					int s2= Dm.getSequenceNos()[(j+1)%2];	// = sk

					int b1= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
					int diff= Dm.getSequenceStarts()[j]- Dl.getSequenceStarts()[i];
					if (diff> 0)
						b1+= diff;
						
					int b2= 0;	// offset-based on the diagonal now! (no longer on base-seq!)
					diff= Dl.getSequenceStarts()[i]- Dm.getSequenceStarts()[j];
					if (diff> 0)
						b2+= diff;

						// calc s_Dn
					String s1Up=			// the sub-sequences 
						Dl.getSubsequence(baseSequences[s1],((i+1)%2)).toUpperCase();
					String s2Up= 
						Dm.getSubsequence(baseSequences[s2],((j+1)%2)).toUpperCase();
					int x= (conslen/ 40);	// tokenize to 40mers (lookup-table!)
					if ((conslen% 40)!= 0)
						++x;
					for (int l= 0; l< x; ++l) {
				
						int match= 0;
						int k;
						for (k= (l* 40); (k < conslen)&& (k< l*40+40); ++k) 
							try {
								match+= Constants.getBlosum()
									[mapping[s1Up.charAt(b1+ k)- 65]]
									[mapping[s2Up.charAt(b2+ k)- 65]];
							} catch (ArrayIndexOutOfBoundsException e) {
								System.err.println("Not supported char in {"+ 
									s1Up.charAt(b1+ k)+ ","+ s2Up.charAt(b2+ k)+ "}, skipped.");
							}
							
							// calc w(Dn)
						double addWgt= 					// lookup
							Constants.getTpProt400()[k- (l* 40)][match];
	
						double factor= ((double)
						(baseSequences[s1].length()
						* baseSequences[s2].length()))
							/ 400d;
						if (addWgt> 0.0000000001)		// calc
							addWgt= 1d- Math.pow((1d- addWgt), factor);
						else
							addWgt= addWgt* factor;
			
						addWgt= (-1f)* ((float) Math.log(addWgt)); 
						relWgt+= addWgt;				// compute rel. Weight

					}
				}
			
		return relWgt;
}	
	
	/**
	  * WARNING:
	  * numerically not stable: should sort results ascending before adding
	  * !!! not consistent with c-code !!!
 	  * 2-decimal floats instead of 6-decimal (c)
	  *  
	  * @deprecated Not capable of comparing normal fragments to ones with jpoints
	  * @see addOW_old
	  */
	public void calcOlWegt_old(MultiFragExt[] frags) {

		if(frags== null)
			return;   

			// reset olw to weight, save for chk
		float[] oldOLW= new float[frags.length];
		for (int i= 0; i < frags.length; i++) {
			oldOLW[i]= frags[i].getOverlapWeight();
			frags[i].setOverlapWeight(frags[i].getWeight());
		}

		for (int d1= 0; d1< frags.length; ++d1) {
          
			 	// pw comparison and add to olw
			MultiFragExt diagonal1= frags[d1];
			for (int d2= (d1+1); d2< frags.length; ++d2) {
				MultiFragExt diagonal2= frags[d2];
				addOW_old(diagonal1 , diagonal2); 
			}
		}
		
			// chk
/*		for (int i= 0; i < frags.length; i++) 
			if (frags[i].getOverlapWeight()!= oldOLW[i])
				System.out.println(
					frags[i].getOverlapWeight()+ "!= "
					+oldOLW[i]
				);
			else
				System.out.println(
					frags[i].getOverlapWeight()+ "== "
					+oldOLW[i]
					+"--- correct!!! ---"
				);
*/				
	}
	
protected int initClosure() {
		
	boolean output= false;
	if (DEBUG)
		output= true;
		
		// get sequence lengths and create 'closure'
	int[] lengths= new int[baseSequences.length];
	for (int i= 0; i< lengths.length; ++i)
		lengths[i]= baseSequences[i].length();
	clos= AliGraphClosure.newAligGraphClosure(
		baseSequences.length, lengths, 0, null);
	int counter= 0;

		// if there is a foce factor, first align the consistent iv until threshold
	if (forceFactor>= 0) {
		for (int i = 0; i < fragments.length; i++) {
			if (fragments[i].isIntervening()&&
				(fragments[i].getOverlapWeight()> forceFactor)&&
				AliGraphClosure.alignableFragment(clos, fragments[i])) {
					fragments[i].setConsistent(true);
					AliGraphClosure.addAlignedFragment(clos, fragments[i]);
					++counter;
					AliGraphClosure.computeClosure(clos);			// this one is necessary...
						
			} else
				fragments[i].setConsistent(false);
	
		}
	}
	if (counter> 0) {
		prealigned= counter;
		System.out.println("prealigned iv fragments: "+prealigned);
	}
			
		// add fragments
	for (int i= 0; i< fragments.length; ++i) {
		
		
		if (AliGraphClosure.alignableFragment(clos, fragments[i])) { 

			fragments[i].setConsistent(true);
			AliGraphClosure.addAlignedFragment(clos, fragments[i]);
			++counter;
			AliGraphClosure.computeClosure(clos);			// this one is necessary...
//			System.out.println("ACCEPT: "+ fragments[i]+ "\n"+
//				fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(true)], true)+ "\n"+
//				fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(false)], false)
//			);

		} else {
//			if (fragments[i].getTranslation()== MultiFragExt.TYPE_REPEAT) {
//				System.out.println("DISCARD "+ fragments[i]+ "\n"+
//					fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(true)], true)+ "\n"+
//					fragments[i].getSubsequence(baseSequences[fragments[i].getSequenceNo(false)], false)
//				);
//			}
			fragments[i].setConsistent(false);
//			if(fragments[i].isRepeat())
//				System.currentTimeMillis();
		}
		
		
		if (DEBUG)
			System.out.println("added: "+ fragments[i]);
	}
		

		
	if (output) {
		System.out.println("Closure Test:");
		int pos1= 1;
		int pos2= 1;
		int a= AliGraphClosure.predFrontier(clos, 0, pos1, 1);
		int b= AliGraphClosure.succFrontier(clos, 0, pos1, 1);
		int c= AliGraphClosure.predFrontier(clos, 1, pos2, 0);
		int d= AliGraphClosure.succFrontier(clos, 1, pos2, 0); 
		int[] b1= AliGraphClosure.getTransitivityBounds(clos, 0, pos1-1, 1);
		int[] b2= AliGraphClosure.getTransitivityBounds(clos, 1, pos2-1, 0);
		System.out.println("["+a+","+b+"], ["+c+","+d+"]");
		System.out.println("["+b1[0]+","+b1[1]+"], ["+b2[0]+","+b2[1]+"]");
		System.out.println("-------------");
	}
					
	return counter;
}

int initClosure(MultiFrag[] restFrags) {
	
	int counter= 0;
	
		// add fragments
	for (int i= 0; i< restFrags.length; ++i) {
		
		if (AliGraphClosure.alignableFragment(clos, restFrags[i])) { 

			fragments[i].setConsistent(true);
			AliGraphClosure.addAlignedFragment(clos, restFrags[i]);
			++counter;
			AliGraphClosure.computeClosure(clos);			// this one is necessary...

		} else {
			fragments[i].setConsistent(false);
		}
		
	}
		
	return counter;	
}
}	

