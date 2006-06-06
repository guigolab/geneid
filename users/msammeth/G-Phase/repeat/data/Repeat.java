/*
 * Created on Nov 26, 2003
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
public class Repeat {

	public static class SequencePositionComparator implements Comparator {

		public int compare(Object o1, Object o2) {
			Repeat r1= (Repeat) o1;
			Repeat r2= (Repeat) o2;
			
			if (r1.getSeqName().equalsIgnoreCase(r2.getSeqName())) {
				if (r1.getStart()< r2.getStart())
					return -1;
				else  
					if (r1.getStart()> r2.getStart())
						return 1;
					else return 0;
			} else	
				return r1.getSeqName().compareTo(r2.getSeqName());
				
		}

	}
	
	protected static int maxFreeID= 0;
	
	protected int ID= -1;
	protected int start= -1, length= -1;
	protected Object type= null;
	protected Repeat next= null, prev= null;
	protected String SeqName= null;
	protected int SeqNb= -1;
	
	protected String alignedSeq= null;
	protected String alignmentID= null;
	
	protected int nb= -1;	// number of repeat within seq 
	
	public static void main(String[] args) {
	}
	
	public Repeat(int newStart, int newLength, Object newType) {
		
		this();
		this.start= newStart;
		this.length= newLength;
		this.type= newType;
	}
	
	public Repeat(int newSeqNb, int newStart, int newLength) {
		
		this();	// get ID
		this.SeqNb= newSeqNb;
		this.start= newStart;
		this.length= newLength;
	}
	
	public Repeat() {
		ID= maxFreeID++;
	}
	
	/**
	 * @return
	 */
	public Object getType() {
		return type;
	}

	/**
	 * @return
	 */
	public int getLength() {
		return length;
	}

	/**
	 * @return
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @return
	 */
	public Repeat getNext() {
		return next;
	}

	/**
	 * @return
	 */
	public Repeat getPrev() {
		return prev;
	}

	/**
	 * @param string
	 */
	public void setType(Object string) {
		type= string;
	}

	/**
	 * @param i
	 */
	public void setLength(int i) {
		length= i;
	}

	/**
	 * @param repeat
	 */
	public void setNext(Repeat repeat) {
		next= repeat;
	}

	/**
	 * @param repeat
	 */
	public void setPrev(Repeat repeat) {
		prev= repeat;
	}

	/**
	 * @param i
	 */
	public void setStart(int i) {
		start= i;
	}

	/**
	 * @return
	 */
	public String getSeqName() {
		return SeqName;
	}

	/**
	 * @param string
	 */
	public void setSeqName(String string) {
		SeqName= string;
	}

	public String toString() {
		
		return getType()+ " ("+ getSeqName()+ ") : "+ getStart()+ " ["+ getLength()+"]";
		
	}
	
	public boolean contains(int pos) {
		
		return (pos>= start&& pos< start+length);
	}
	/**
	 * @return
	 */
	public int getSeqNb() {
		return SeqNb;
	}

	/**
	 * @param i
	 */
	public void setSeqNb(int i) {
		SeqNb= i;
	}

	/**
	 * @return
	 */
	public int getID() {
		return ID;
	}

	/**
	 * @param i
	 */
	public void setID(int i) {
		ID= i;
	}

	/**
	 * @return
	 */
	public String getAlignedSeq() {
		return alignedSeq;
	}

	/**
	 * @param string
	 */
	public void setAlignedSeq(String string) {
		alignedSeq= string;
	}

	/**
	 * @return
	 */
	public String getAlignmentID() {
		return alignmentID;
	}

	/**
	 * @param string
	 */
	public void setAlignmentID(String string) {
		alignmentID= string;
	}

	/**
	 * @return
	 */
	public int getNb() {
		return nb;
	}

	/**
	 * @param i
	 */
	public void setNb(int i) {
		nb= i;
	}

}
