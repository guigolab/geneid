package gphase.regex;

import java.util.EmptyStackException;
import java.util.Vector;

import gphase.model.SpliceSite;
import gphase.tools.Arrays;
import gphase.tools.IntStack;
import gphase.tools.IntVector;
import gphase.tools.Stack;

public class ChessMaster {

	public static void main(String[] args) {
		
		SpliceSite[][] sc= new SpliceSite[][] {
				{new SpliceSite(null, 1006, false)},
				{new SpliceSite(null, 1001, false), new SpliceSite(null, 1002, true), new SpliceSite(null, 1003, false), new SpliceSite(null, 1004, true), new SpliceSite(null, 1005, false)}
		};
		
		RegExp reg1= new RegExp("*,A,i");
		RegExp reg2= new RegExp("*,A,(i+1)|");
		Automaton[] autos= new Automaton[] {reg1.toAutomaton(), reg2.toAutomaton()};
		
		ChessMaster kasparov= new ChessMaster(autos, sc);
		System.out.println(kasparov.match());
	}
	
	Automaton[] autos;
	SpliceSite[][] spliceChains;
	
	private Stack[] branchAutos;
	private Stack[] branchSC;
	private IntStack branchVal;
	public ChessMaster(Automaton[] newAutos, SpliceSite[][] newSpliceChains) {
		this.autos= newAutos;
		this.spliceChains= newSpliceChains;
	}
	
	
	void pushBranchPoint(Automaton[] autos, int min, Stack[] sc) {
		
		for (int i = 0; i < autos.length; i++) { 
			branchAutos[i].push(autos[i].getCurrentState());
			branchSC[i].push(sc[i].clone());
		}
		branchVal.push(min);
	}
	
	int popBranchPoint(Automaton[] autos, Stack[] sc) {
		if (branchVal.isEmpty())	// dont eat yellow snow!
			return (-1);
		
		for (int i = 0; i < autos.length; i++) 
			autos[i].setCurrentState((State) branchAutos[i].pop());
		for (int i = 0; i < sc.length; i++) 
			sc[i]= (Stack) branchSC[i].pop();
		
		return branchVal.pop();
	}
	
	boolean match(Automaton[] autos, int min, Stack[] sc) {
		
		int[] minAutos= findAutos(autos, min);	// determine current minimum
		if (min< 0|| minAutos== null|| minAutos.length< 1) {
			int x;
			for (x = 0; x < autos.length; x++) 
				if ((!autos[x].isEndState()) 
						|| (autos[x].currentState.isEnd()	// $
							&& !sc[x].isEmpty())) 
					break;
			if (x== autos.length) 
				return true; 	// success: a hit
			
			try {
				min= popBranchPoint(autos, sc);
			} catch (EmptyStackException e) {
				return false;	// failed: no more branch points
			}
			if (min< 0 || autos== null)
				return false; 	// failed: no more branch points
			return match(autos, findMin(autos, min), sc);	// try last branch point
		}
		
		int x;	// if all minima provide another way, create a branch for backtracking
		for (x = 0; x < minAutos.length; x++) 
			if (!autos[minAutos[x]].getCurrentState().hasHigherTransition(min))
				break;
		if (x== minAutos.length)
			pushBranchPoint(autos, min, sc);
		
			// get chunks of equally labelled transitions
		Transition[][] chunk= new Transition[minAutos.length][];
		int[] chunkPos= new int[chunk.length];
		for (int i = 0; i < minAutos.length; i++) { 
			chunk[i]= autos[minAutos[i]].getChunk(min);
			if (chunk[i]== null)
				chunkPos[i]= -1;
			else
				chunkPos[i]= 0;
		}
		
			// try positional move
		int pos= -1;	// check for consistency amongst min autos
		for (x = 0; x < minAutos.length; x++) {
			SpliceSite s;
			try {
				s= ((SpliceSite) sc[minAutos[x]].pop());
			} catch (Exception e) {
				break;
			}
			if (pos< 0) { 
				pos= s.getPos();
				if (s.isDonor()== chunk[x][0].isDonor())
					continue;
				else
					break;
			}
			if ((s.isDonor()!= chunk[x][0].isDonor())
					|| (pos!= s.getPos())
					|| (autos[minAutos[x]].getCurrentState().isStart()&&	// ^
							s!= spliceChains[minAutos[x]][0]))
				break;
		}
		if (x< minAutos.length) 	// failed, not equal splice sites
			return false;			
				
		int[] succAutos= findAutos(autos, findMin(autos, min));	// check successors (not in same sequence -> always true)
		for (int i = 0; succAutos!= null&& i < succAutos.length; i++) {
			SpliceSite s;
			try {
				s= ((SpliceSite) sc[succAutos[i]].top());
			} catch (Exception e) {
				return false;
			}
			if (s.getPos()>= 0&& s.getPos()<= pos) 	// failed
				return false;
		}


			// positional move is possible, lets try all structural different possibilities
		State[] saveMinAutos= new State[autos.length];
		Stack[] saveMinSC= new Stack[autos.length];
		for (int i = 0; i < saveMinAutos.length; i++) {  
			saveMinAutos[i]= autos[i].getCurrentState();
			saveMinSC[i]= (Stack) sc[i].clone();
		}
		boolean first= true;
		while (first|| incrementChunks(chunk, chunkPos)) {
			first= false;
			for (x = 0; x < minAutos.length; x++)	// make move 
				autos[minAutos[x]].move(chunk[x][chunkPos[x]]);
			boolean result= match(autos, findMin(autos, -1), sc);	// recursion, continue with next min
			if (result== true)
				return true;	// return first hit found			
			for (x = 0; x< saveMinAutos.length; x++) {	// recover position 
				autos[x].setCurrentState(saveMinAutos[x]);
				sc[x]= saveMinSC[x];
			}
		}
		
		return false;	// chunk possibilities exhausted wo finding hit
	}
	
	private boolean incrementChunks(Transition[][] chunks, int[] chunkPos) {
		
		int i;
		for (i = 0; i < chunkPos.length; i++) 
			if (chunkPos[i]+ 1== chunks[i].length) 
				chunkPos[i]= 0;
			else {
				++chunkPos[i];
				break;
			}
		
		if (i== chunkPos.length)
			return false;
		return true;
	}


	/**
	 * @deprecated
	 * @param autos
	 * @param min
	 * @param sc
	 * @param blackList
	 * @return
	 */
	boolean match_old(Automaton[] autos, int min, Stack[] sc) {
		
		int[] minAutos= findAutos(autos, min);	// determine current minimum
		if (minAutos== null|| minAutos.length< 1) {
			for (int x = 0; x < autos.length; x++) 
				if (autos[x].isEndState())
					return false; 	// failed
			return true;	// success
		}
		
		int x;	// if all minima provide another way, create a branch for backtracking
		for (x = 0; x < minAutos.length; x++) 
			if (!autos[minAutos[x]].getCurrentState().hasHigherTransition(min))
				break;
		if (x== minAutos.length)
			pushBranchPoint(autos, min, sc);
		
			// try the move
		int pos= -1;	// check for consistency amongst min autos
		for (x = 0; x < minAutos.length; x++) {
			autos[minAutos[x]].move(min);
			if (pos< 0)
				pos= ((SpliceSite) sc[minAutos[x]].pop()).getPos();
			else if (pos!= ((SpliceSite) sc[minAutos[x]].pop()).getPos())
				break;
		}
		if (x< minAutos.length) {	// failed, not equal splice sites
			min= popBranchPoint(autos, sc);
			if (autos== null)
				return false;			
			return match(autos, findMin(autos, min), sc);
		}
		
		int[] succAutos= findAutos(autos, findMin(autos, min));	// check successors
		for (int i = 0; succAutos!= null&& i < succAutos.length; i++) {
			if (((SpliceSite) sc[succAutos[i]].top()).getPos()<= pos) {	// failed
				min= popBranchPoint(autos, sc);
				if (autos== null)
					return false;
				min= findMin(autos, min);					// recursion, continue with next min
				return match(autos, findMin(autos, min), sc);
			}
		}
		return match(autos, findMin(autos, min), sc);	// recursion, continue with next min		
	}
	
	/**
	 * 
	 * @param autos
	 * @param threshold
	 * @return the next minimum value in the current outgoing transitions
	 * in the automatons <code>autos</code>, or <code>-1</code> if all 
	 * no such minimum is found (all automatons in end state or dead end).
	 */
	 int findMin(Automaton[] autos, int threshold) {
		int min = Integer.MAX_VALUE;
		for (int i = 0; i < autos.length; i++) {
			Transition[] trans= autos[i].getCurrentState().getTransitions();
			for (int j = 0; trans!= null&& j < trans.length; j++) 
				if (trans[j].getLabel()< min&& trans[j].getLabel()> threshold) 
					min= trans[j].getLabel();
		}
		
		if (min== Integer.MAX_VALUE)
			return (-1);
		return min;
	}
	
	 /**
	  * 
	  * @param autos
	  * @param val
	  * @return automatons found in <code>autos</code> with the 
	  * <code>val</code> value in one of the current outgoing edges
	  * or <code>null</code>.
	  */
	 int[] findAutos(Automaton[] autos, int val) {
		IntVector result= new IntVector(autos.length);
		for (int i = 0; i < autos.length; i++) {
			Transition[] trans= autos[i].getCurrentState().getTransitions();
			for (int j = 0; trans!= null&& j < trans.length; j++) {
				if (trans[j].getLabel()== val) {
					result.add(i);
					break;
				}
			}
		}
		
		if (result.size()< 1)
			return null;		
		return result.toIntArray();
	}
		
	boolean recursivePermute(Integer[] chosen) {
		
		if (chosen.length== autos.length) {	// abort
			Automaton[] autos= map(chosen);
			int min= findMin(autos, -1);
			
			branchVal= new IntStack();	// re-init branch points
			for (int i = 0; i < branchAutos.length; i++) {
				branchAutos[i]= new Stack();
				branchSC[i]= new Stack();
			}
			
			Stack[] sc= new Stack[spliceChains.length];		// init ss stacks
			for (int i = 0; i < spliceChains.length; i++) {
				sc[i]= new Stack();
				for (int j = spliceChains[i].length- 1; j >= 0; j--) 
					sc[i].push(spliceChains[i][j]);
			}
			
			return match(autos, min, sc);
		}
	
			// genrate new combination
		for (int i = 0; i < autos.length; i++) { 	// check for free positions
			int j;
			for (j = 0; j < chosen.length; j++) 
				if (i== chosen[j].intValue())
					break;
			if (j< chosen.length)
				continue;
			
			boolean found= recursivePermute((Integer[]) Arrays.extendField(chosen, new Integer(i)));	// recursion
			if (found)
				return found;
			else
				for (int k = 0; k < autos.length; k++) 
					autos[k].init();	// re-init automatons				
		}
		
		return false;
	}
	
	Automaton[] map(Integer[] mapTable) {
		Automaton[] result= new Automaton[autos.length];
		for (int i = 0; i < result.length; i++) 
			result[i]= autos[mapTable[i].intValue()];
		return result;
	}
	
	public boolean match() {
		
		branchAutos= new Stack[autos.length];
		branchSC= new Stack[spliceChains.length];
		for (int i = 0; i < branchAutos.length; i++) { 
			branchAutos[i]= new Stack();
			branchSC[i]= new Stack();
		}
		branchVal= new IntStack();
		
		return recursivePermute(new Integer[]{});
	}
}
