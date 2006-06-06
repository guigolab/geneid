/*
 * Created on Nov 27, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.Serializable;
import java.util.Comparator;

/**
 * 
 * 
 * @author msammeth
 */
public abstract class AbstractRegion implements Region {

	public static final int REGION_5UTR= 1;
	public static final int REGION_CDS= 2;
	public static final int REGION_3UTR= 3;

	static final long serialVersionUID=  5443375142823871946L;
	public abstract Species getSpecies();

	public abstract String getChromosome();

	/* assumption: exons share same strand
	 * 
	 * @author micha
	 */
	public static class PositionComparator implements Comparator {
		
		public boolean equals(Object o1, Object o2) {
			
			if (compare(o1, o2)== 0)
				return true;
			return false;
		}
			
		public int compare(Object o1, Object o2) {
				
			if (o1.equals(o2))	// object identity
				return 0;
			
				// non-overlapping, one before the other
			if (((Region) o1).getEnd()< ((Region) o2).getStart())
				return -1;		// one stops before the other
			if (((Region) o2).getEnd()< ((Region) o1).getStart())
				return 1;
			
				// overlapping: none stops before the other
			if (((Region) o1).getStart()< ((Region) o2).getStart())
				return -1;
			if (((Region) o2).getStart()< ((Region) o1).getStart())
				return 1;
			
				// overlapping and same start position
			if (((Region) o1).getEnd()< ((Region) o2).getEnd())
				return -1;
			if (((Region) o2).getEnd()< ((Region) o1).getEnd())
				return 1;
			
			return 0;	// identical positions
			
		}
	}

	public boolean overlaps(Region anotherRegion) {
		
		if (getChromosome()!= null&& anotherRegion.getChromosome()!= null
			&& (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome())))
			return false;
		
		if ((getStart()>= anotherRegion.getStart()&& getStart()< anotherRegion.getEnd())
				|| (anotherRegion.getStart()>= getStart()&& anotherRegion.getStart()< getEnd()))
			return true;
		return false;
	}
	
	/**
	 * Returns <code>true</code> if <code>this</code> region contains <code>anotherRegion</code>.
	 */
	public boolean contains(Region anotherRegion) {

		if (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome()))
			return false;
		
		if ((getStart()<= anotherRegion.getStart())&& (getEnd()>= anotherRegion.getEnd()))
			return true;
		return false;
	}

	int end = Integer.MIN_VALUE;
	int start = Integer.MAX_VALUE;
	/**
	 * @return Returns the end.
	 */
	public int getEnd() {
		return end;
	}

	/**
	 * @return Returns the start.
	 */
	public int getStart() {
		return start;
	}

	/**
	 * @param end The end to set.
	 */
	public void setEnd(int end) {
		this.end = end;
	}

	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		this.start = start;
	}

	public String toString() {
		return "["+getStart()+";"+getEnd()+"]";
	}
}
