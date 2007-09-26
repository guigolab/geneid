/*
 * Created on Oct 14, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JSeparator;

import repeat.data.AlignedTupel;
import repeat.data.Alignment;
import repeat.data.Block;
import repeat.data.Repeat;
import sun.rmi.runtime.GetThreadPoolAction;

/**
 * Amalgamates the information of the <i>schema-</i> (repeat types, bounds) 
 * and the <i>centre-</i> files (repeat alignments, features, bounds). The 
 * intermediately created <code>Repeat</code> objects may be written out to
 * suitable file formats.
 * 
 * @author micha
 */
public class BAliManager {
	
	public static final String GAP_CHARS= "-.~ ";
	protected String fileBase= null;
	protected String fileTestcase= null;
	protected String filePath= null;
	
	protected Repeat[] repeats= null;
	protected Alignment[] alis= null;
	protected String[] seqNames= null;
	
	protected AlignedTupel[][][][][] tupels= null;
	protected int[][][][] coreTupelCount= null;	
	public static void converterGroup6(String pathToG6) {

			// retrieve sub-dirs		
		File dir= new File(pathToG6);
		String[] subDirs= dir.list();
		Arrays.sort(subDirs);

			// init writer for logging
		BufferedWriter writer= null;
		try {
			writer= new BufferedWriter(new FileWriter(File.createTempFile("amalgamator","log")));	
			writer.write("started "+ new Date(System.currentTimeMillis())+"\n\n");
		} catch (Exception e) {
			; // nothing, again..
		}
		
			// process subdirs		
		for (int i= 0; i < subDirs.length; i++) {
//			if (!subDirs[i].trim().equals("test"))
//				continue;
			System.out.println("\n== ["+subDirs[i]+"] ==");
			try {
				writer.write("\n== ["+subDirs[i]+"] ==\n");
			} catch (IOException e2) {
				e2.printStackTrace();
			}
			File sdir= new File(pathToG6+File.separator+subDirs[i]);
			String[] files= sdir.list();
			for (int j= 0; j < files.length; j++) {
				int idx= files[j].indexOf("_schema");
				if (idx< 0)
					continue;				
//				if (files[j].indexOf("lrr_ref6")< 0)
//					continue;
				String procFile= files[j].substring(0,idx);
				System.out.print("processing "+procFile+"..");
				try {
					writer.write("processing "+procFile+"..");
				} catch (IOException e) {
					; // :)
				}
				try {
					BAliManager converter= new BAliManager(
						sdir.getAbsolutePath()+File.separator+procFile);
					RepeatPoundWrapper repper=
						new RepeatPoundWrapper(converter.fileBase+".rep", converter.getRepeats());
					repper.write();
				} catch (Exception e) {
					e.printStackTrace(); 
					try {
						writer.write(e.toString());
						writer.write(e.getMessage());
					} catch (IOException e1) {
						; //:)
					}
				}
				System.out.println("done.");
				try {
					writer.write("done.\n");
					writer.flush();
				} catch (IOException e1) {
					; // :)
				}
			}
		}
		
	}
	
	public static void main(String[] args) {

			converterGroup6("D:\\Eigene Dateien\\repeats\\data\\ref6");
	// single test
//			BAliManager balu= 
//				new BAliManager("D:\\Eigene Dateien\\repeats\\data\\ref6\\test\\ank_ref6");
//	
//			RepeatPoundWrapper repper=
//				new RepeatPoundWrapper(balu.fileBase+".rep", balu.getRepeats());
//			repper.write();
//			System.out.println("Wrote "+balu.getRepeats().length+" repeats:");
//			System.out.println(balu.getRepeats()[0]);			
//			System.out.println(balu.getRepeats()[balu.getRepeats().length-1]);			
	// read test
//			RepeatPoundWrapper repper2=
//				new RepeatPoundWrapper(balu.fileBase+".rep");
//			repper2.getRepeats();
//			System.out.println("Read "+repper2.getRepeats().length+" repeats:");
//			System.out.println(repper2.getRepeats()[0]);			
//			System.out.println(repper2.getRepeats()[repper2.getRepeats().length-1]);			
			
			

	}
	
	public BAliManager(String newFileBase) {
	
		this.fileBase= newFileBase;		// "<path>/Testcase_ref6"
		int pos1= fileBase.lastIndexOf(File.separator);
		int pos2= fileBase.lastIndexOf('_');
		this.fileTestcase= newFileBase.substring(pos1+1, pos2);
	}
	
	public Repeat[] getRepeats() {
		if (repeats == null) {
			repeats= readRepeats();
		}
		return repeats;
	}
	
	protected Repeat[] readRepeats() {
		
			// preprocessed file available
		File preprocFile= new File(fileBase+".rep");
		if (preprocFile.exists()) {
			RepeatPoundWrapper pWrapper= new RepeatPoundWrapper(preprocFile.getAbsolutePath());
			repeats= pWrapper.getRepeats();
			return repeats;
		}
			
			// amalgamate information from html-files
		repeats= readRepeatSchema();
		alis= readRepeatAlignments();
		for (int i= 0; i < alis.length; i++) 
			alis[i]= formatAli(alis [i]);
		
		repeats= mergeRepeats(repeats, alis);
		return repeats;
	}
	
	public int[][][][] getCoreTupelCount() {

		if (coreTupelCount == null) {
			
			tupels= getAlignedTupels();
			coreTupelCount= new int[tupels.length][][][];
			for (int i= 0; (tupels!= null)&& i < tupels.length; i++) {
				coreTupelCount[i]= new int[tupels[i].length][][];
				for (int j= 0; (tupels[i]!= null)&& j < tupels[i].length; j++) {
					if (i== j)		// same sequences -> null
						continue;
					coreTupelCount[i][j]= new int[tupels[i][j].length][];
					for (int k= 0; (tupels[i][j]!= null)&& k < tupels[i][j].length; k++) {
						coreTupelCount[i][j][k]= new int[tupels[i][j][k].length];
						for (int l= 0; (tupels[i][j][k]!=null)&& l < tupels[i][j][k].length; l++) {
							int coreCounter= 0;
							for (int m= 0; (tupels[i][j][k][l]!= null)&& m < tupels[i][j][k][l].length; m++) 
								if (tupels[i][j][k][l][m].isCore())
									++coreCounter;
							coreTupelCount[i][j][k][l]= coreCounter;
						}
					}
				}
			}
			
		}

		return coreTupelCount;
	}
	
	/**
	 * get aligned character pairs
	 * sorted according to seq names
	 * @return 3D tupel array [seqNameA][seqNameB][repNbA][repNbB][]
	 */		
	public AlignedTupel[][][][][] getAlignedTupels() {
		
		if (tupels == null) {
			repeats= getRepeats();
			if (repeats== null)
				return null;
		
				// get all seq names
			Vector nameVec= new Vector();
			for (int i= 0; i < repeats.length; i++) {
				int j;
				for (j= 0; j < nameVec.size(); j++) {
					if (repeats[i].getSeqName().equalsIgnoreCase((String) nameVec.elementAt(j)))
						break;
				}
				if (j>= nameVec.size())
					nameVec.add(repeats[i].getSeqName());
			}
			seqNames= new String[nameVec.size()];
			for (int i= 0; i < seqNames.length; i++) 
				seqNames[i]= (String) nameVec.elementAt(i);
			Arrays.sort(seqNames);
		
				// get tupels for all seq pairs
			tupels= new AlignedTupel[seqNames.length][seqNames.length][][][];
			for (int i= 0; i < tupels.length; i++) 
				for (int j= 0; j < tupels.length; j++) { 
					if (i== j)
						continue;
					tupels[i][j]= getAlignedTupels(seqNames[i], seqNames[j]);	// mirror
//					AlignedTupel[][][] mirror= 
//						new AlignedTupel[tupels[i][j][0].length][tupels[i][j].length][];
//					for (int k= 0; k < mirror.length; k++) {
//						for (int l= 0; l < mirror[k].length; l++) {
//							mirror[k][l]= new AlignedTupel[tupels[i][j][l][k].length];
//							for (int m= 0; m < mirror[k][l].length; m++) {
//								mirror[k][l][m]= new AlignedTupel();
//								mirror[k][l][m].setCharA(tupels[i][j][l][k][m].getCharB());
//								mirror[k][l][m].setCharB(tupels[i][j][l][k][m].getCharA());
//								mirror[k][l][m].setCore(tupels[i][j][l][k][m].isCore());
//								mirror[k][l][m].setPositionA(tupels[i][j][l][k][m].getPositionB());
//								mirror[k][l][m].setPositionB(tupels[i][j][l][k][m].getPositionA());
//								mirror[k][l][m].setRepNbA(tupels[i][j][l][k][m].getRepNbB());
//								mirror[k][l][m].setRepNbB(tupels[i][j][l][k][m].getRepNbA());
//								mirror[k][l][m].setSeqNameA(tupels[i][j][l][k][m].getSeqNameB());
//								mirror[k][l][m].setSeqNameB(tupels[i][j][l][k][m].getSeqNameA());
//							}
//						}
//					}
//					tupels[j][i]= mirror;
				}
		}
		
		return tupels;
	}
	
	/**
	 * returns all alignment tupels for a seq pair
	 * 
	 * @param seqA
	 * @param seqB
	 * @return AlignedTupel[repeatA][repeatB][]
	 */
	protected AlignedTupel[][][] getAlignedTupels(String seqA, String seqB) {
		
			// retrieve repeats for the sequences
		Vector seqARepVec= new Vector();
		Vector seqBRepVec= new Vector();
		for (int i= 0; i < repeats.length; i++) {
			if (repeats[i].getSeqName().equals(seqA))
				seqARepVec.add(repeats[i]);
			if (repeats[i].getSeqName().equals(seqB))
				seqBRepVec.add(repeats[i]);
		}
		Repeat[] seqARep= new Repeat[seqARepVec.size()];
		for (int i= 0; i < seqARepVec.size(); i++) 
			seqARep[i]= (Repeat) seqARepVec.elementAt(i);
		Repeat[] seqBRep= new Repeat[seqBRepVec.size()];
		for (int i= 0; i < seqBRepVec.size(); i++) 
			seqBRep[i]= (Repeat) seqBRepVec.elementAt(i);
		
			// get tupels for each repeat pair (of the seq pair)
		AlignedTupel[][][] tupels= new AlignedTupel[seqARep.length][seqBRep.length][];
		for (int i= 0; i < tupels.length; i++) {
			for (int j= 0; j < tupels[i].length; j++) {
				seqARep[i].setNb(i);
				seqBRep[j].setNb(j);
				tupels[i][j]= getAlignedTupels(seqARep[i], seqBRep[j]);
			}
		}
		
		return tupels;
	}
	
	/**
	 * tupels for a certain repeat pair
	 * @param seqARep
	 * @param seqBRep
	 * @param repA
	 * @param repB
	 * @return 	AlignedTupel[] with the aligned positions
	 * 			<code>null</code> if both repeat types are not aligned
	 */
	protected AlignedTupel[] getAlignedTupels(Repeat repA, Repeat repB) {
		
		if (!repA.getAlignmentID().equals(repB.getAlignmentID()))
			return null;			// not alignable / not aligned in reference
		
			// get sequences and tokenize in tupels
		String aliA= repA.getAlignedSeq();
		String aliB= repB.getAlignedSeq();
		if (aliA.length()!= aliB.length())
			System.err.println("alignements of unequal length: "+repA.getSeqName()+"_"+repA.getNb()+
				" and "+ repB.getSeqName()+"_"+repB.getNb());
		
		AlignedTupel tmpTupel;
		Vector tupelVec= new Vector();
		int gapCharsA= 0;
		int gapCharsB= 0;
		for (int i= 0; i < aliA.length(); i++) {		// iterate
			
			if (GAP_CHARS.indexOf(aliA.charAt(i))> 0		// skip gaps
					|| GAP_CHARS.indexOf(aliB.charAt(i))> 0) {
				if (GAP_CHARS.indexOf(aliA.charAt(i))> 0)
					++gapCharsA;
				if (GAP_CHARS.indexOf(aliB.charAt(i))> 0)
					++gapCharsB;
				continue;
			}
			
			tmpTupel= new AlignedTupel();
			tmpTupel.setCore(Character.isUpperCase(aliA.charAt(i))&& 
								Character.isUpperCase(aliB.charAt(i)));
			tmpTupel.setCharA(aliA.charAt(i));
			tmpTupel.setCharB(aliB.charAt(i));
			tmpTupel.setPositionA(repA.getStart()+ i- gapCharsA);
			tmpTupel.setPositionB(repB.getStart()+ i- gapCharsB);
			tmpTupel.setRepNbA(repA.getNb());
			tmpTupel.setRepNbB(repB.getNb());
			tmpTupel.setSeqNameA(repA.getSeqName());
			tmpTupel.setSeqNameB(repB.getSeqName());
			tupelVec.add(tmpTupel);
		}
		
			// convert
		AlignedTupel[] tupels= new AlignedTupel[tupelVec.size()];
		for (int i= 0; i < tupels.length; i++) {
			tupels[i]= (AlignedTupel) tupelVec.elementAt(i);
		}
		Arrays.sort(tupels,new AlignedTupel.NaturalOrderComparator());
		
		return tupels;
	}
	
	protected Alignment formatAli(Alignment ali) {
		
			// init: all lowercase
		for (int i= 0; i < ali.getSequences().length; i++) 
			ali.setSequence(i, ali.getSequences()[i].toLowerCase());
		
			// only core blocks to uppercase..
		LinkedList blocks= ali.getBlocks();
		Iterator iter= blocks.iterator();
		while (iter.hasNext()) {
			
			Block tmpBlock= (Block) iter.next();
			if (tmpBlock.isCBlock()) {
				
				String seq= ali.getSequence(tmpBlock.getSeqID()[0]);		// tokenize
				String pfx= seq.substring(0, tmpBlock.getStartPos()[0]);
				String core= seq.substring(tmpBlock.getStartPos()[0], 
										tmpBlock.getStartPos()[0]+ tmpBlock.getLength());
				String sfx= seq.substring(tmpBlock.getStartPos()[0]+ tmpBlock.getLength(),
										seq.length());
				
				core= core.toUpperCase();		// uppercase for core regions
				seq= pfx+ core+ sfx;
				
				ali.setSequence(tmpBlock.getSeqID()[0], seq);	// set new Sequence
			}
		}
		
//		System.out.println(ali);
		return ali;
	}
	
	protected Repeat[] readRepeatSchema() {

//		System.out.print("reading in schema...");
		BAliParserSchema schemaParser= new BAliParserSchema(fileBase+ "_schema.html");
		System.err.flush();
//		System.out.println("done.");
		
		return schemaParser.getRepeats();
	}
	
	protected Alignment[] readRepeatAlignments() {
		
		String fileName= fileBase+ "_centre.html";
		if (new File(fileName).exists()) 				// easy, just one alignment
			return (new Alignment[] {readRepeatAlignment(fileName)});
		else {
			int separate= fileName.lastIndexOf(File.separator);
			String path= fileName.substring(0, separate+1);
			String fName= fileName.substring(separate+1, fileName.length());
			separate= fName.indexOf('_');		// try to read several alignments
			String prefix= fName.substring(0, (separate+ 1));			// pfx+"_"
			String suffix= fName.substring(separate, fName.length());	// _sfx
			
			Vector collectedAlis= new Vector();
			for (int i= 1;;++i) {
				fName= path+ prefix+ i+ suffix;
				if (new File(fName).exists()) {
					String optSchema= fName.substring(0, fName.lastIndexOf('_'))+ "_schema.html";
					if (new File(optSchema).exists()) 
						readAdditionalSubschema(optSchema, repeats);
					collectedAlis.add(readRepeatAlignment(fName));
				} else
					break;	// no more alignments
			}
			Alignment[] result= new Alignment[collectedAlis.size()];
			for (int i= 0; i < result.length; i++) 
				result[i]= (Alignment) collectedAlis.elementAt(i);
			return (result.length> 0) ? result: null;
		}
	}
	
	
	/**
	 * Reads an secondary schema file (see test4/sh3).
	 * 
	 * @param fName
	 */
	protected Repeat[] readAdditionalSubschema(String fName, Repeat[] repeats) {
		
		BAliParserSchema parser= new BAliParserSchema(fName);
		Repeat[] rep= parser.getRepeats();
		
			// merge with existing repeats
		Vector addReps= new Vector();
		for (int i= 0; i < rep.length; i++) {
			for (int j= 0; (repeats!= null)&& j < repeats.length; j++) {
				if (rep[i].getSeqName().equalsIgnoreCase(repeats[j].getSeqName())
					&& rep[i].getStart()== repeats[j].getStart()
					&& rep[i].getLength()== repeats[j].getLength()
					) {
						if (repeats[j].getType()!= null)			// merge by concat types
							repeats[j].setType(repeats[j].getType().toString()+ 
										"_"+rep[i].getType().toString());
						
					} else
						addReps.add(rep[i]);						// add new repeat
			}
		}
		
		if (addReps.size()> 0) {									// extent array if necessary
			Repeat[] newRepeats= new Repeat[repeats.length+ addReps.size()];
			for (int i= 0; i < repeats.length; i++) 
				newRepeats[i]= repeats[i];
			for (int i= 0; i < addReps.size(); i++) 
				newRepeats[i+ repeats.length]= (Repeat) addReps.elementAt(i);
			repeats= newRepeats;
		}
		
		return repeats;
	}
	
	protected Alignment readRepeatAlignment(String fileName) {
		
		BAliParserCentre aliParser= new BAliParserCentre(fileName);
		Alignment ali= aliParser.getAlignment();
//		System.out.println(ali);

		return ali;
	}
	
	
	protected Repeat[] mergeRepeats(Repeat[] repeats, Alignment[] alis) {
		
		for (int i= 0; i < alis.length; i++) {
			for (int j= 0; j < alis[i].getSequences().length; j++) {
				int k;
				for (k= 0; k < repeats.length; k++) {
					if (repeats[k].getStart()== alis[i].getStartPos()[j]
// endpos is not a marker -> bug
//						&& (repeats[k].getStart()+ repeats[k].getLength()- 1)== alis[i].getEndPos()[j]
						&& alis[i].getSeqIDs()[j].startsWith(repeats[k].getSeqName())
						) {
						
							if (repeats[k].getAlignedSeq()!= null)
								System.err.println("Repeat "+ repeats[k]+ " double choosen!");
							repeats[k].setAlignedSeq(alis[i].getSequences()[j]);
							repeats[k].setAlignmentID(fileTestcase+ (i+1));
							if ((repeats[k].getStart()+ repeats[k].getLength()- 1)!= alis[i].getEndPos()[j]) {
								System.err.println("Length Bug for sequence with gap end corrected! "+repeats[k].getSeqName());
								repeats[k].setLength(alis[i].getEndPos()[j]- 
														alis[i].getStartPos()[j]+ 1);
							}
							break;
					}
				}
				if (k>= repeats.length)
					System.err.println(alis[i].getSeqIDs()[j]+
							"("+alis[i].getStartPos()[j]+","+alis[i].getEndPos()[j]+")"+ 
							" not found in schema, purged.");
			}
		}
		
			// eliminate sequences w/o alignment
			// e.g. all additional domains of group 3
		Vector result= new Vector();
		for (int i= 0; i < repeats.length; i++) {
			if (repeats[i].getAlignedSeq()== null)
				System.err.println("Repeat: "+repeats[i]+" has no aligned sequence, purged.");
			else
				result.add(repeats[i]);
		}
		Repeat[] repeats2= new Repeat[result.size()];
		for (int i= 0; i < result.size(); i++) 
			repeats2[i]= (Repeat) result.elementAt(i);
		
		return repeats2;
	}
	
	
	/**
	 * @return
	 */
	public String[] getSeqNames() {
		return seqNames;
	}
	
	public Alignment[] getAlis() {
		alis= readRepeatAlignments();
		return alis;
	}

}
