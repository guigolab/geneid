package gphase.model;

import java.util.Comparator;

public abstract class DirectedRegion extends AbstractRegion {

	/**
	 * 
	 * @author micha (written in the "Café de Indias" coffee shop 
	 * in Seville)
	 *
	 */
	public static class PositionComparator extends AbstractRegion.PositionComparator {
		// allows exons on different strands
		public int compare(Object o1, Object o2) {
			
			DirectedRegion r1, r2;
			try {
				r1= (DirectedRegion) o1;
				r2= (DirectedRegion) o2;
			} catch (ClassCastException e) {
				return super.compare(o1, o2);
			}
			
			int start1= r1.isForward()? r1.getStart():-r1.getStart();	// reverse strand exons are still given backward:
			int end1= r1.isForward()? r1.getEnd():-r1.getEnd();			// (start is the 5'position= end of rev exons)
			int start2= r2.isForward()? r2.getStart():-r2.getStart();
			int end2= r2.isForward()? r2.getEnd():-r2.getEnd();
			
			if (start1> start2)
				return 1;
			if (start1< start2)
				return -1;
			if (end1> end2)
				return 1;
			if (end1< end2)
				return -1;
			return 0;
		}
	}
	public int strand = 0;

	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== (strand> 0));
	}
	
	public boolean isForward() {
		return strand== 1;
	}

	public boolean isLaggingStrand() {
		return isReverse();
	}

	public boolean isLeadingStrand() {
		return isForward();
	}

	public boolean isReverse() {
		return strand== -1;
	}

	/**
	 * @param end The end to set.
	 */
	public void setEnd(int end) {
		if (strand== 0)	// not inited 
			throw new RuntimeException("Error: set strand before start/end!");
		
		if (strand> 0&& end> 0)
			this.end = Math.abs(end);	// force pos
		else if (strand< 0)
			this.end= -Math.abs(end);	// force neg
	}

	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		if (strand== 0)	// not inited 
			throw new RuntimeException("Error: set strand before start/end!");
		
		if (strand> 0&& start> 0)
			this.start = Math.abs(start);	// force pos
		else if (strand< 0)
			this.start= -Math.abs(start);	// force neg
	}

	public void setStrand(int strand) {
		this.strand = strand;
	}

	public int getStrand() {
		return strand;
	}

	public boolean overlaps(DirectedRegion anotherRegion) {
		
		if (getChromosome()!= null&& anotherRegion.getChromosome()!= null
			&& (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome())))
			return false;
		
		int start1, end1, start2, end2;
		if (isForward()) {
			start1= getStart();
			end1= getEnd();
		} else if (isReverse()) {
			start1= -getEnd();
			end1= -getStart();
		} else
			throw new RuntimeException("Error: strand not set!");
		if (anotherRegion.isForward()) {
			start2= anotherRegion.getStart();
			end2= anotherRegion.getEnd();
		} else if (anotherRegion.isReverse()) {
			start2= -anotherRegion.getEnd();
			end2= -anotherRegion.getStart();
		} else
			throw new RuntimeException("Error: strand not set!");
		
		if ((start1>= start2&& start1< end2)
				|| (start2>= start1&& start2< end1))
			return true;
		return false;
	}

}
