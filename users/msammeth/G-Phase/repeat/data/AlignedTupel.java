/*
 * Created on Oct 15, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.data;

import java.util.Comparator;

/**
 * 
 * 
 * @author micha
 */
public class AlignedTupel {

	protected boolean core= false;
	protected int positionA= -1;
	protected int positionB= -1;
	protected char charA= '*';
	protected char charB= '*';
	protected int repNbA= -1;
	protected int repNbB= -1;
	protected String seqNameA= "";
	protected String seqNameB= "";
	
	/**
	 * Sorts aligned tupels according to natural order with decision hierarchy:
	 * 		positionA
	 * 		positionB
	 */
	public static class NaturalOrderComparator implements Comparator {
		public int compare(Object o1, Object o2) {
				
			AlignedTupel t1= (AlignedTupel) o1;
			AlignedTupel t2= (AlignedTupel) o2;

				// check for equal sequences, swap if necessary
			if (t1.getSeqNameA().equals(t2.getSeqNameA())) {		
				if (!t1.getSeqNameB().equals(t2.getSeqNameB()))	{	// unequal
					// throw new ClassCastException();
					return (t1.getSeqNameB().compareTo(t2.getSeqNameB()));
				}
				// else check postions
			} else if (t1.getSeqNameA().equals(t2.getSeqNameB())) { // over-cross
				AlignedTupel t= new AlignedTupel();
				t.setPositionA(t2.getPositionB());
				t.setPositionB(t2.getPositionA());
				t2= t;
				if (!t1.getSeqNameB().equals(t2.getSeqNameA())){// swap condition
					// throw new ClassCastException();
					return t1.getSeqNameB().compareTo(t2.getSeqNameA());
				} 
				// else check position
			} else 
				return t1.getSeqNameA().compareTo(t2.getSeqNameA());
			
				// delta in 1st position
			if (t1.getPositionA()< t2.getPositionA())
				return (-1);
			else if (t1.getPositionA()> t2.getPositionA())
				return 1;
				
				// else: 1st pos equal					
			if (t1.getPositionB()< t2.getPositionB())
				return (-1);
			else if (t1.getPositionB()> t2.getPositionB())
				return 1;
					
			return 0;
		}
		
		public boolean equals(Object obj) {
			return compare(this, obj)== 0;
		}

	}
	
	/**
	 * @return
	 */
	public char getCharA() {
		return charA;
	}

	/**
	 * @return
	 */
	public char getCharB() {
		return charB;
	}

	/**
	 * @return
	 */
	public boolean isCore() {
		return core;
	}

	/**
	 * @return
	 */
	public int getPositionA() {
		return positionA;
	}

	/**
	 * @return
	 */
	public int getPositionB() {
		return positionB;
	}

	/**
	 * @return
	 */
	public int getRepNbA() {
		return repNbA;
	}

	/**
	 * @return
	 */
	public int getRepNbB() {
		return repNbB;
	}

	/**
	 * @return
	 */
	public String getSeqNameA() {
		return seqNameA;
	}

	/**
	 * @return
	 */
	public String getSeqNameB() {
		return seqNameB;
	}

	/**
	 * @param c
	 */
	public void setCharA(char c) {
		charA= c;
	}

	/**
	 * @param c
	 */
	public void setCharB(char c) {
		charB= c;
	}

	/**
	 * @param b
	 */
	public void setCore(boolean b) {
		core= b;
	}

	/**
	 * @param i
	 */
	public void setPositionA(int i) {
		positionA= i;
	}

	/**
	 * @param i
	 */
	public void setPositionB(int i) {
		positionB= i;
	}

	/**
	 * @param i
	 */
	public void setRepNbA(int i) {
		repNbA= i;
	}

	/**
	 * @param i
	 */
	public void setRepNbB(int i) {
		repNbB= i;
	}

	/**
	 * @param string
	 */
	public void setSeqNameA(String string) {
		seqNameA= string;
	}

	/**
	 * @param string
	 */
	public void setSeqNameB(String string) {
		seqNameB= string;
	}
	
	public String toString() {
		String res="["+seqNameA+","+seqNameB+"] "+positionA+","+positionB+" ("+charA+","+charB+")";
		if (core)
			res+="*";
		if (repNbA>= 0)
			res+= "-"+repNbA+","+repNbB+"-";
			
		return res;
	}

}
