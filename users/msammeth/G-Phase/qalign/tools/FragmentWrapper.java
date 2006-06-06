package qalign.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

import qalign.algo.dialign.MultiFrag;
import qalign.model.Block;

/**
 * @author sammeth
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class FragmentWrapper implements DecorationWrapper {

	protected String absFName= null;
	
		// information of ff-File
	protected String comment= null;
	protected int[] lengths= null;
	protected String[] names= null;
	static final String FF_LENGTH_RECCOGNITION="seq_len:";
	static final String FF_NAME_RECOGNITION="sequences:";
	protected MultiFrag[] fragments= null;
	protected Block[] blocks= null;
	
	/**
	 * Constructor for DialignFASTAWrapper.
	 * @param newFName
	 * @param newFPath
	 */
	public FragmentWrapper(String newFName, String newFPath) {
		
		this.absFName= newFPath+ File.separator+ newFName;
	}

	/**
	 * Constructor for DialignFASTAWrapper.
	 * @param absFName
	 */
	public FragmentWrapper(String absFName) {
		
		this.absFName= absFName;
	}

	public static void main(String[] args) {
	}
	
	public void readFF() {
		
//		System.out.println(fPath+":::"+fName);

		String tmpStr;
		
		try {
			BufferedReader fptr= new BufferedReader(
						new FileReader(absFName));
						
				// read comment		
			tmpStr= fptr.readLine().trim();
			comment= "";
			while (!tmpStr.startsWith(FF_LENGTH_RECCOGNITION)) {
				comment+= tmpStr+ "\n";
				tmpStr= fptr.readLine().trim();
			}
			
				// read header
			readFFSeqLen(tmpStr);
			tmpStr= fptr.readLine().trim();
			while (!tmpStr.startsWith(FF_NAME_RECOGNITION))
				tmpStr= fptr.readLine().trim();
			readFFSeqNames(tmpStr);
			
				// read fragments
			tmpStr= fptr.readLine().trim();
			while (tmpStr.length()< 2)
				tmpStr= fptr.readLine().trim();		// blank lines
			Vector frgs= new Vector();
			frgs.add(readFFFragment(tmpStr));
			while (fptr.ready()) {
				tmpStr= fptr.readLine().trim();		// blank lines
				frgs.add(readFFFragment(tmpStr));
			}			
				
				// convert & save
			fragments= new MultiFrag[frgs.size()];
			for (int i= 0; i< frgs.size(); ++i) {
				fragments[i]= (MultiFrag) frgs.elementAt(i);
//				System.out.println(fragments[i]);
			}
			
				// end
			fptr.close();
		
		} catch (Exception e) {	// catch all IO's
			e.printStackTrace();
		}
		
	}
	
	public void readFFSeqLen(String read) {
		
		StringTokenizer st= new StringTokenizer(read, " ");
		lengths= new int[st.countTokens()- 1];
	
		st.nextToken(); 	// "seq_len: "	
		try {

			for (int i= 0; i< lengths.length; ++i) 
				lengths[i]= Integer.parseInt(st.nextToken());

		} catch (NumberFormatException e) {
			e.printStackTrace();
		}
	}
	
	public MultiFrag readFFFragment(String read) {
		
		MultiFrag result= new MultiFrag();
		
		StringTokenizer st= new StringTokenizer(read, " ");
	
			// read seq no's
		String nb= st.nextToken().trim(); 	// "23) "	
		result.setNumber(Integer.parseInt(nb.substring(0,nb.length()-1)));
		st.nextToken(); 	// "seq: "	
		try {
			result.setSequenceNos(
				Integer.parseInt(st.nextToken())- 1,	// 1-based
				Integer.parseInt(st.nextToken())- 1);
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}
		
			// read seq starts
		st.nextToken(); 	// "beg: "	
		try {
			result.setSequenceStarts(
				Integer.parseInt(st.nextToken())- 1,	// 1-based
				Integer.parseInt(st.nextToken())- 1);
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

			// read seq length
		st.nextToken(); 	// "len: "	
		try {
			result.setLength(
				Integer.parseInt(st.nextToken()));
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

			// read seq weight
		st.nextToken(); 	// "wgt: "	
		try {
			result.setWeight(
				Float.parseFloat(st.nextToken()));
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

			// read seq overlapping weight
		st.nextToken(); 	// "olw: "	
		try {
			result.setOverlapWeight(
				Float.parseFloat(st.nextToken()));
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

			// read iteration no
		st.nextToken(); 	// "it: "	
		try {
			result.setIteration(
				Integer.parseInt(st.nextToken()));	// 1-based ??!
		} catch (NumberFormatException e) {
			e.printStackTrace();
		}

			// read seq consistency
		if (st.nextToken().startsWith("cons"))
			result.setConsistent(true);
		else
			result.setConsistent(false);

			// read seq consistency
			// ??! check with Burkhard
		if (st.nextToken().startsWith("N-frag"))
			result.setTranslation(0);
		else
			result.setTranslation(1);
			
		return result;
	}	
	
	public void readFFSeqNames(String read) {
		
		StringTokenizer st= new StringTokenizer(read, " ");
		names= new String[st.countTokens()- 1];
	
		st.nextToken(); 	// "sequences: "	
		for (int i= 0; i< names.length; ++i) 
			names[i]= st.nextToken().replace('_', ' ');
	}		
	/**
	 * Returns the fragments.
	 * @return MultiFrag[]
	 */
	public MultiFrag[] getFragments() {
		
		return fragments;
	}
	
	/**
	 * Returns the lengths.
	 * @return int[]
	 */
	public int[] getLengths() {
		return lengths;
	}

	/**
	 * Returns the names.
	 * @return String[]
	 */
	public String[] getNames() {
		return names;
	}


	/**
	 * returns the blocks induced by the (consistent) fragments
	 * and creates them if they are not already initialised (converted).
	 */
	public Block[] getBlocks() {
		
		if (blocks== null) {
		
			blocks= convertToBlocks();
		}
		
		return blocks;
	}


	/**
	 * converts the MultiFrags read in to Blocks.
	 */
	public Block[] convertToBlocks() {
		
			// error, abort
		if ((fragments== null)|| (names== null))	// names needed for block size init !
			return null;
			
			// convert MultiFrags to Blocks
		Vector blockVector= new Vector();	// to store temporary results
		for (int i= 0; i< fragments.length; ++i) {
			blockVector= addMultiFrag(fragments[i], blockVector);
		}
		
			// convert to Block[]
		Block[] result= new Block[blockVector.size()];
		for (int i= 0; i< result.length; ++i)
			result[i]= (Block) blockVector.elementAt(i);
			
		return result;
	}
	
	/**
	 * adds a MultiFrag to an existing Block vector.
	 */
	public Vector addMultiFrag(MultiFrag newFrag, Vector blockVector) {
		
			// only consistent fragments are added
		if (!newFrag.isConsistent())
			return blockVector;
		
			// create Block from MultiFrag
		int[] start= new int[names.length];
		int[] end= new int[names.length];
		for (int i= 0; i< start.length; ++i) 
			start[i]= end[i]= (-1);
		start[newFrag.getSequenceNo(true)]= newFrag.getSequenceStart(true);
		start[newFrag.getSequenceNo(false)]= newFrag.getSequenceStart(false);
		end[newFrag.getSequenceNo(true)]= newFrag.getSequenceStart(true)+ newFrag.getLength();
		end[newFrag.getSequenceNo(false)]= newFrag.getSequenceStart(false)+ newFrag.getLength();
		Block newBlock= new Block(start, end);
		newBlock.setAnnotation(newFrag.getNumber()+ ") wgt: "+ newFrag.getWeight()+ ", olw: "+ newFrag.getOverlapWeight());
	
			// check for overlaps
		if (checkOverlaps(newBlock, blockVector))
			return blockVector;
			
			// no Block to append to found, create new
		blockVector.add(newBlock);
		return blockVector;
	}
	
	
	public boolean checkOverlaps(Block newBlock, Vector blockVector) {
		
		for (int i= 0; i< blockVector.size(); ++i) {
			
			if ( ((Block) blockVector.elementAt(i)).isIdenticalArea(newBlock) ) {

				((Block) blockVector.elementAt(i)).extend(newBlock);
				return true;
			}
		}
		
		return false;
	}
	
}
