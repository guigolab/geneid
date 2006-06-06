/*
 * Created on Nov 27, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.data;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.Vector;

import qalign.algo.dialign.MultiFrag;
import qalign.model.MultipleAlignmentModel;


/**
 * 
 * 
 * @author micha
 */
public class Alignment {

	public static void main(String[] args) {
	}
	
	protected LinkedList blocks= null;
	protected int[] startPos= null;
	protected int[] endPos= null;
	protected String[] sequences= null;
	protected String[] seqIDs= null;
	
	public Alignment() {
	}
	
	public Alignment(String[] newSequences) {
		setSequences(newSequences);
	}
	
	public String getSequence(String someID) {
		
		for (int i= 0; i < seqIDs.length; i++) 
			if (seqIDs[i].equals(someID))
				return sequences[i];
		
		return null;
	}
	
	public void setSequence(String someID, String newSeq) {
		
		for (int i= 0; i < seqIDs.length; i++) 
			if (seqIDs[i].equals(someID))
				sequences[i]= newSeq;
	}
	
	public void setSequence(int index, String newSeq) {
		
		sequences[index]= newSeq;
	}	
	
	public boolean addBlock(Block newBlock) {
		if (blocks== null)
			blocks= new LinkedList();
		
			// look if already in list
		Iterator biter= blocks.iterator();
		while (biter.hasNext())
			if (newBlock.equals((Block) biter.next()))
				return false;
				
		blocks.add(newBlock);
		return true;
	}
	
	/**
	 * @return
	 */
	public int[] getEndPos() {
		return endPos;
	}

	/**
	 * @return
	 */
	public int[] getStartPos() {
		return startPos;
	}

	/**
	 * @param is
	 */
	public void setEndPos(int[] is) {
		endPos= is;
	}

	/**
	 * @param is
	 */
	public void setStartPos(int[] is) {
		startPos= is;
	}

	/**
	 * @return
	 */
	public String[] getSeqIDs() {
		return seqIDs;
	}

	/**
	 * @return
	 */
	public String[] getSequences() {
		return sequences;
	}

	/**
	 * @param strings
	 */
	public void setSeqIDs(String[] strings) {
		seqIDs= strings;
	}

	/**
	 * @param strings
	 */
	public void setSequences(String[] strings) {
		
			// init
		this.sequences= strings;
		if (sequences== null)
			return;
		
			// derive attributes
		if (startPos== null&& endPos== null) {	
			this.startPos= new int[sequences.length];
			this.endPos= new int[sequences.length];
			int counter;
			for (int i= 0; i < endPos.length; i++) {
				counter= 0;
				for (int j= 0; j < sequences[i].length(); j++) 
					if (!MultipleAlignmentModel.isGapChar(sequences[i].charAt(j)))
						++counter;
				if (counter> 0) {
					startPos[i]= 1;
					endPos[i]= counter;
				} else
					startPos[i]= endPos[i]= 0;
			}
		}
	}

	public String toString() {
		
		String result= "--- Alignment ---\n";
		for (int i= 0; i< seqIDs.length; ++i) 
			result+= seqIDs[i]+ " ("+ startPos[i]+ ","+ endPos[i]+ ")\t"+ sequences[i]+"\n";
/*		if(blocks!= null) {
			result+= "\nassociated Blocks: "+ blocks.size()+ "\n";
			Iterator biter= blocks.iterator();
			while (biter.hasNext()) {
//				result+= biter.next()+ "\n";
				Block tmpB= (Block) biter.next();
				if (!tmpB.isRArea())
					continue;				
				result+= tmpB+ " : "+ getSequence(tmpB.getSeqID()[0]).substring(
				tmpB.getStartPos()[0], tmpB.getStartPos()[0]+ tmpB.getLength())+ "\n"; 
			}
/*			for (int i= 0; i< 1; ++i) {
			
				Block tmpB= (Block) biter.next();
				if (!tmpB.isRArea()) {
					i--;
					continue;
				}
				result+= " "+ tmpB;
				try {
					result+= " : "+ getSequence(tmpB.getSeqID()[0]).substring(
						tmpB.getStartPos()[0], tmpB.getStartPos()[0]+ tmpB.getLength()); 
				} catch (Exception e) {
					; // :)
				}
				result+= "\n";
			}
		}
*/			
		
		result+= "--- ---";
		
		return result;
	}
	
	/**
	 * Decomposes a alignment layout into pairwise, gap-free diagonals.
	 * Uppercase for aligned positions is mandatory!
	 * @return
	 */
	public MultiFrag[] getMultiFrags() {
		
		if (sequences== null)
			return null;
			
		int[] startIdx= null, actIdx= null, gapCnt= new int[sequences.length], oldGapCnt= new int[gapCnt.length];
		for (int i= 0; i < gapCnt.length; i++) 
			gapCnt[i]= oldGapCnt[i]= 0;
		Vector results= new Vector();
		for (int j= 0; j< sequences[0].length(); ++j) {	// iterate over all positions

			actIdx= new int[sequences.length];	// the index of active sequences 
			for (int i= 0; i < actIdx.length; i++) 
				actIdx[i]= j;
			boolean skipCol= false;
			
				// find active positions
			for (int i= 0; i < sequences.length; ++i) {		// iterate over all sequences
				
				char c= sequences[i].charAt(j);
				
				if (Character.isLetter(c)) {
					if (Character.isUpperCase(c)) 
						actIdx[i]= j;	// mark sequences participating
					else {	// lowercase letter -> skip col
						skipCol= true;		
					}
				} else {
					if (MultipleAlignmentModel.isGapChar(c)) {
						actIdx[i]= (-1);
						++gapCnt[i];
					} else {	// no letter and no gap char -> skip column
						skipCol= true;
						break;
					}
				}
			}
			
				// compare against start of current batch
			boolean newBatch= false;
			if (startIdx== null) 
				startIdx= actIdx;
			else {
				for (int i= 0; i < actIdx.length; i++) 
					if ((actIdx[i]< 0 && startIdx[i]>= 0)
						|| (actIdx[i]>= 0 && startIdx[i]< 0)) 
					
						newBatch= true;
			}
			
			int gapCount= 0;						// don't start batch for only one sequence
			for (int i= 0; i < actIdx.length; i++) 
				if (actIdx[i]< 0)
					++gapCount; 
			if (gapCount>= actIdx.length- 1)
				skipCol= true;
			
				// create new batch
			if (skipCol || newBatch) {
				
				int x= 0;
				for (;x< startIdx.length; ++x)		// find first seq
					if (startIdx[x]>= 0)
						break;
				for (int y= (x+1);y< startIdx.length; ++y) {	// find all other sequences and chain
					
					if (startIdx[y]< 0)
						continue;
					
					MultiFrag newFrag= new MultiFrag();
					newFrag.setSequenceNos(x, y);
					newFrag.setSequenceStarts(startIdx[x]- oldGapCnt[x], startIdx[y]- oldGapCnt[y]);
					newFrag.setLength(j- startIdx[x]);		// +1 -1
					newFrag.setConsistent(true);
					results.add(newFrag);					
				}
				
				if (skipCol)				// init new batch
					startIdx= null;
				else 
					for (int i= 0; i < startIdx.length; i++) 
						startIdx[i]= actIdx[i];
				
				for (int i= 0; i< oldGapCnt.length; i++) 
					oldGapCnt[i]= gapCnt[i];
			}
			
		}	// end iterate over all positions
		
		if (startIdx!= null) {		// close last batch

			int x= 0;
			for (;x< startIdx.length; ++x)		// find first seq
				if (startIdx[x]>= 0)
					break;
			for (int y= (x+1);y< startIdx.length; ++y) {	// find all other sequences and chain
					
				if (startIdx[y]< 0)
					continue;
					
				MultiFrag newFrag= new MultiFrag();
				newFrag.setSequenceNos(x, y);
				newFrag.setSequenceStarts(startIdx[x]- oldGapCnt[x], startIdx[y]- oldGapCnt[y]);
				newFrag.setLength(actIdx[x]- startIdx[x]+ 1);		// +1
				newFrag.setConsistent(true);
				results.add(newFrag);					
			}
				
		}
		
		MultiFrag[] res= new MultiFrag[results.size()];
		for (int i= 0; i < res.length; i++) 
			res[i]= (MultiFrag) results.elementAt(i);
		return res;
	}
	
	public void addSequence(String newSequence) {
	
		if (sequences== null)							// extend array 
			sequences= new String[1];
		else {
			String[] tmp= sequences;
			sequences= new String[tmp.length+ 1];
			for (int i= 0; i < tmp.length; i++) 
				sequences[i]= tmp[i];
		}
		
		sequences[sequences.length- 1]= newSequence;	// perform extension
	}		

	public void addSeqID(String newSeqID) {
		
		if (seqIDs== null)							// extend array 
			seqIDs= new String[1];
		else {
			String[] tmp= seqIDs;
			seqIDs= new String[tmp.length+ 1];
			for (int i= 0; i < tmp.length; i++) 
				seqIDs[i]= tmp[i];
		}
			
		seqIDs[seqIDs.length- 1]= newSeqID;	// perform extension
	}		
	
	public void addStartPos(int newStartPos) {
		
		if (startPos== null)							// extend array 
			startPos= new int[1];
		else {
			int[] tmp= startPos;
			startPos= new int[tmp.length+ 1];
			for (int i= 0; i < tmp.length; i++) 
				startPos[i]= tmp[i];
		}
			
		startPos[startPos.length- 1]= newStartPos;		// perform extension
	}		

	public void addEndPos(int newEndPos) {
			
		if (endPos== null)							// extend array 
			endPos= new int[1];
		else {
			int[] tmp= endPos;
			endPos= new int[tmp.length+ 1];
			for (int i= 0; i < tmp.length; i++) 
				endPos[i]= tmp[i];
		}
				
		endPos[endPos.length- 1]= newEndPos;		// perform extension
	}		

	/**
	 * @return
	 */
	public LinkedList getBlocks() {
		return blocks;
	}
	
	public AlignedTupel[][][] getAlignedTupels() {
		
			// assuming the alignment has been read
		AlignedTupel[][][] result= new AlignedTupel[seqIDs.length][seqIDs.length][];
		Vector tupelVec= null;
		final String GAP_CHARS= "-.~ ";
		for (int i= 0; i < result.length; i++) {
			for (int j= 0; j < result[i].length; j++) {
				if (j== i)
					continue;
				int ctrA= 0;
				int ctrB= 0;
				tupelVec= new Vector();
				for (int k= 0; k < sequences[i].length(); k++) {
					
					char charA= sequences[i].charAt(k);		// 0-based pos counters
					char charB= sequences[j].charAt(k);
					
					if (GAP_CHARS.indexOf(charA)>= 0) {		// gaps
						if(GAP_CHARS.indexOf(charB)< 0)
							++ctrB;
						continue;
					}
					if (GAP_CHARS.indexOf(charB)>= 0) {
						if(GAP_CHARS.indexOf(charA)< 0)
							++ctrA;
						continue;
					}
					
					AlignedTupel newTupel= new AlignedTupel();	// new tupel
					newTupel.setCharA(charA);
					newTupel.setCharB(charB);
					newTupel.setSeqNameA(seqIDs[i]);
					newTupel.setSeqNameB(seqIDs[j]);
					newTupel.setPositionA(ctrA);
					newTupel.setPositionB(ctrB);
					tupelVec.add(newTupel);
					
					++ctrA;++ctrB;
				}
				result[i][j]= new AlignedTupel[tupelVec.size()];	// convert
				for (int k= 0; k < result[i][j].length; k++) 
					result[i][j][k]= (AlignedTupel) tupelVec.elementAt(k);
				
			}
		}
		
		
		return result;
	}

}