package gphase.model.intransparent.init;

import gphase.io.gtf.GTFObject;
import gphase.tools.Arrays;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

public class DirectedRegion extends DefaultRegion {

	static final long serialVersionUID = 4346170163999111167l;
	public static class StartComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			DirectedRegion reg0= (DirectedRegion) o1;
			DirectedRegion reg1= (DirectedRegion) o2;
			if (reg0.get5PrimeEdge()< reg1.get5PrimeEdge())
				return -1;
			if (reg0.get5PrimeEdge()> reg1.get5PrimeEdge())
				return 1;
			return 0;
		}
	}

	public static class EndComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			DirectedRegion reg0= (DirectedRegion) o1;
			DirectedRegion reg1= (DirectedRegion) o2;
			if (reg0.get3PrimeEdge()< reg1.get3PrimeEdge())
				return -1;
			if (reg0.get3PrimeEdge()> reg1.get3PrimeEdge())
				return 1;
			return 0;
		}
	}

	public static class EndsBeforeComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			DirectedRegion reg0= (DirectedRegion) o1;
			DirectedRegion reg1= (DirectedRegion) o2;
			if (reg0.get3PrimeEdge()< reg1.get5PrimeEdge())
				return -1;
			if (reg1.get3PrimeEdge()< reg0.get5PrimeEdge())
				return 1;
			return 0;
		}
	}

	public static class LengthComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			DirectedRegion reg0= (DirectedRegion) arg0;
			DirectedRegion reg1= (DirectedRegion) arg1;
			int len0= reg0.getLength();
			int len1= reg1.getLength();
			if (len0< 0|| len1< 0) {
				System.err.println("Invalid region length "+len0+", "+len1);
			}
			if (len0< len1)
				return -1;
			if (len0> len1)
				return 1;
			return 0;
		}
	}
	public static DirectedRegion[] intersect(DirectedRegion[] dir1 , DirectedRegion[] dir2) {
		Vector interV= new Vector();
		for (int i = 0; i < dir1.length; i++) {
			for (int j = 0; j < dir2.length; j++) {
				if (dir1[i].overlaps(dir2[j]))
					interV.add(dir1[i].intersect(dir2[j]));
			}
		}
		return (DirectedRegion[]) Arrays.toField(interV);
	}
	
	
	/**
	 * returns sites contained in ANY of the regions. inefficient
	 * 
	 * @param regs
	 * @param s
	 * @return
	 */
	public static AbstractSite[] contained(DirectedRegion[] regs, AbstractSite[] s) {
		if (regs== null|| s== null)
			return null;
		
		Vector v= new Vector(s.length);
		Comparator compi= new AbstractSite.PositionComparator();
		for (int i = 0; regs!= null&& i < regs.length; i++) {
			AbstractSite[] as= contained(regs[i], s);
			v= Arrays.addUnique(v, as, compi);
		}
		
		return (AbstractSite[]) Arrays.toField(v);
	}
	
	public static DirectedRegion getUnion(DirectedRegion[] regs) {
		if (regs== null)
			return null;
		int min= Integer.MAX_VALUE;
		int max= Integer.MIN_VALUE;
		for (int i = 0; i < regs.length; i++) {
			if (regs[i].get5PrimeEdge()< min)
				min= regs[i].get5PrimeEdge();
			if (regs[i].get3PrimeEdge()> max)
				max= regs[i].get3PrimeEdge();
		}
		
		DirectedRegion result= new DirectedRegion(min, max, regs[0].getStrand());
		result.setChromosome(regs[0].getChromosome());
		result.setSpecies(regs[0].getSpecies());
		return result;
	}
	
	
	/**
	 * deletes contained regions, concatenates overlapping regions, st every nt is only covered once
	 * @param regs
	 * @return
	 */
	public static DirectedRegion[] getUniqueRegions(DirectedRegion[] regs) {
		if (regs== null)
			return null;
		regs= filterContainedRegions(regs);
		Comparator compi= new PositionComparator();
		java.util.Arrays.sort(regs, compi);
		Vector remV= new Vector(regs.length);
		Vector addV= new Vector();
		for (int i = 0; i < regs.length; i++) {
			int j;
			for (j = 0; j < remV.size(); j++) 
				if (remV.elementAt(j)== regs[i])
					break;
			if (j< remV.size())	// already concatenated
				continue;
			
			DirectedRegion concatReg= (DirectedRegion) regs[i].clone();
			for (j = (i+1); j < regs.length; j++) {
				if (concatReg.overlaps(regs[j])) {
					concatReg= concatenate(concatReg, regs[j]);
					remV.add(regs[j]);
				}
			}
			if (compi.compare(regs[i], concatReg)== 0) {
				remV.add(regs[i]);
				addV.add(concatReg);
			}
		}
		
		Vector keepV= new Vector(regs.length- remV.size());
		for (int i = 0; i < regs.length; i++) {
			int j;
			for (j = 0; j < remV.size(); j++) 
				if (remV.elementAt(j)== regs[i])
					break;
			if (j== remV.size())
				keepV.add(regs[i]);
		}
		
		for (int j = 0; j < addV.size(); j++) 
			keepV.add(addV.elementAt(j));
		
		return (DirectedRegion[]) Arrays.toField(keepV);
	}

	public static DirectedRegion[] filterIdenticalRegions(DirectedRegion[] regs) {
		if (regs== null)
			return null;
		Comparator compi= new PositionComparator();
		Vector remV= new Vector(regs.length);
		for (int i = 0; i < regs.length; i++) 
			for (int j = i+1; j < regs.length; j++) 
				if (compi.compare(regs[i], regs[j])== 0)
					remV.add(regs[j]);
		
		Vector keepV= new Vector(regs.length- remV.size());
		for (int i = 0; i < regs.length; i++) {
			int j;
			for (j = 0; j < remV.size(); j++) 
				if (regs[i]== remV.elementAt(j))
					break;
			if (j== remV.size())
				keepV.add(regs[i]);
		}
		return (DirectedRegion[]) Arrays.toField(keepV);
	}
	
	public static DirectedRegion[] filterContainedRegions(DirectedRegion[] regs) {
		if (regs== null)
			return null;
		regs= filterIdenticalRegions(regs);
		java.util.Arrays.sort(regs, new LengthComparator());
		
		Comparator compi= new PositionComparator();
		Vector remV= new Vector(regs.length);
		for (int i = 0; i < regs.length; i++) 
			for (int j = i+1; j < regs.length; j++) 
				if (regs[i].contains(regs[j]))
					remV.add(regs[j]);
		
		Vector keepV= new Vector(regs.length- remV.size());
		for (int i = 0; i < regs.length; i++) {
			int j;
			for (j = 0; j < remV.size(); j++) 
				if (regs[i]== remV.elementAt(j))
					break;
			if (j== remV.size())
				keepV.add(regs[i]);
		}
		return (DirectedRegion[]) Arrays.toField(keepV);
	}
	
	public static DirectedRegion concatenate(DirectedRegion reg0, DirectedRegion reg1) {
		DirectedRegion reg= new DirectedRegion(
				Math.min(reg0.getStart(), reg1.getStart()),
				Math.max(reg0.getEnd(), reg1.getEnd()),
				reg0.getStrand());
		reg.setChromosome(reg0.getChromosome());
		reg.setSpecies(reg0.getSpecies());
		return reg;
	}
	
	public static AbstractSite[] contained(DirectedRegion dir, AbstractSite[] s) {
		if (dir== null|| s== null)
			return null;
		Vector v= new Vector();
		for (int i = 0; i < s.length; i++) 
			if (dir.contains(s[i]))
				v.add(s[i]);
		return (AbstractSite[]) Arrays.toField(v);
	}
	
	public boolean contains(SpliceSite ss) {
		if (!ss.getGene().getChromosome().equals(chromosome))
			return false;
		
		int apos= Math.abs(ss.getPos());
		if (apos>= Math.abs(getStart())&& apos<= Math.abs(getEnd()))
			return true;
		return false;
	}


	public static DirectedRegion[] unite_old(DirectedRegion[] dir1, DirectedRegion[] dir2) {
		
		Comparator compi= new PositionComparator();
		java.util.Arrays.sort(dir1, compi);
		java.util.Arrays.sort(dir2, compi);
		
		
		Vector interV= new Vector();
		int i= 0, j= 0;
		DirectedRegion lastCluster= null;
		while (i < dir1.length&& j < dir2.length) {
			
			if (lastCluster== null) {
				if (dir1[i].overlaps(dir2[j])) {
					lastCluster= dir1[i].merge(dir2[j]);
					++i; ++j;
				} else {
					if (Math.abs(dir1[i].getStart())< Math.abs(dir2[j].getStart())) {
						interV.add(dir1[i++]);	// give up to merge this
					} else {
						interV.add(dir2[j++]);	// give up to merge this
					}
				}
			} else {	// try to extend last cluster
				if (lastCluster.overlaps(dir1[i])) {
					lastCluster= lastCluster.merge(dir1[i]);
					++i;
				} else if (lastCluster.overlaps(dir2[j])) {
					lastCluster= lastCluster.merge(dir2[j]);
					++j;
				} else {
					interV.add(lastCluster);
					lastCluster= null;
				}
			}
			
		}
		if (lastCluster!= null)
			interV.add(lastCluster);
		for (; i < dir1.length; i++) 
			interV.add(dir1[i]);
		for (; j < dir2.length; j++) 
			interV.add(dir2[j]);
		
		return (DirectedRegion[]) Arrays.toField(interV);
	}
	
	public Object clone() {
		DirectedRegion reg= new DirectedRegion();
		reg.setStrand(getStrand());
		reg.setStart(getStart());
		reg.setEnd(getEnd());
		reg.setChromosome(getChromosome());
		reg.setSpecies(getSpecies());
		reg.setID(getID());
		return reg;
	}
	
	@Override
	public boolean equals(Object obj) {
		
		if (obj== null|| !(obj instanceof DirectedRegion))
			return false;
		
		DirectedRegion otherReg= (DirectedRegion) obj;
		if (otherReg.getStart()== getStart()&& otherReg.getEnd()== getEnd()
				&& otherReg.getChromosome().equals(getChromosome()))
			return true;
		return false;
	}
	
	public static class OrderComparator implements Comparator {
		/**
		 * sorts regions according to their order on pos/neg strand 
		 */
		public int compare(Object o1, Object o2) {
			DirectedRegion r1, r2;
			try {
				r1= (DirectedRegion) o1;
				r2= (DirectedRegion) o2;
			} catch (ClassCastException e) {
				return -1;
			}
			if (r1.get5PrimeEdge()< r2.get5PrimeEdge())
				return -1;
			if (r1.get5PrimeEdge()> r2.get5PrimeEdge())
				return 1;
			if (r1.get3PrimeEdge()< r2.get3PrimeEdge())
				return -1;
			if (r1.get3PrimeEdge()> r2.get3PrimeEdge())
				return 1;
			return 0;
		}
		
		
	}
	/**
	 * 
	 * @author micha (written in the "Caf? de Indias" coffee shop 
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
	/**
	 * for ordering regions ascending according to their order in -1/+1
	 * genes
	 * @author micha 
	 *
	 */
	public static class DirectedPositionComparator extends AbstractRegion.PositionComparator {
		// allows exons on different strands
		public int compare(Object o1, Object o2) {
			
			DirectedRegion r1, r2;
			try {
				r1= (DirectedRegion) o1;
				r2= (DirectedRegion) o2;
			} catch (ClassCastException e) {
				return super.compare(o1, o2);
			}
			
			int start1= r1.getStart();	
			int end1= r1.getEnd();			// (start is the 5'position= end of rev exons)
			int start2= r2.getStart();
			int end2= r2.getEnd();
			
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
	public DirectedRegion(int newStart, int newEnd, int strand) {
		setStrand(strand);
		// now checked in setStart(), setEnd()
//		if(Math.abs(newStart)> Math.abs(newEnd)) {
//			int h= newStart;
//			newStart= newEnd;
//			newEnd= h;
//		}
		setStart(newStart);
		setEnd(newEnd);
	}
	public DirectedRegion(int newStart, int newEnd, int strand, String newChr) {
		this(newStart, newEnd, strand);
		setChromosome(newChr);
	}
	
	public DirectedRegion(DirectedRegion source) {
		setStrand(source.getStrand());
		setStart(source.getStart());
		setEnd(source.getEnd());
	}
	
	public DirectedRegion(GTFObject obj) {
		setStrand(obj.getStrand());
		setStart(obj.getStart());
		setEnd(obj.getEnd());
		setChromosome(obj.getChromosome());		
	}
	
	public DirectedRegion() {
	}
	
	/**
	 * Returns <code>true</code> if <code>this</code> region contains <code>anotherRegion</code>.
	 */
	public boolean contains(Region anotherRegion) {

		if (!anotherRegion.getChromosome().equalsIgnoreCase(getChromosome()))
			return false;
		
		if ((Math.abs(getStart())<= Math.abs(anotherRegion.getStart()))&& 
				(Math.abs(getEnd())>= Math.abs(anotherRegion.getEnd())))
			return true;
		return false;
	}
	

	
	public boolean isUpstream(int pos) {
		if ((isForward()&& pos< start)||
			(!isForward()&& pos< end))
			return true;
		return false;
			
	}
	
	public boolean isDownstream(int pos) {
		if ((isForward()&& pos> end)||
			(!isForward()&& pos> start))
			return true;
		return false;
			
	}
	
	public int getLength() {
		return (get3PrimeEdge()- get5PrimeEdge()+ 1);
	}
	
	public boolean contains(int pos) {
		if ((isForward()&& pos>= start&& pos<= end)||
			(!isForward()&& pos<= start&& pos>= end))
			return true;
		return false;
	}
	
	public boolean contains(AbstractSite ss) {
		if (ss== null)
			return false;
		if (ss.getGene().getStrand()!= getStrand()||
				ss.getGene().getChromosome()!= getChromosome())
			return false;
		return contains(ss.getPos());
	}
	
	
	
	public DirectedRegion intersect(DirectedRegion anotherRegion) {
		
		if (!overlaps(anotherRegion))
			return null;
		
		DirectedRegion dir= new DirectedRegion(Math.max(Math.abs(getStart()), Math.abs(anotherRegion.getStart())),
				Math.min(Math.abs(getEnd()), Math.abs(anotherRegion.getEnd())), getStrand());
		dir.setChromosome(getChromosome());
		return dir;
	}
	
	public DirectedRegion merge(DirectedRegion anotherRegion) {
		
		if (!overlaps(anotherRegion))
			return null;
		
		DirectedRegion dir= new DirectedRegion(Math.min(Math.abs(getStart()), Math.abs(anotherRegion.getStart())),
				Math.max(Math.abs(getEnd()), Math.abs(anotherRegion.getEnd())), getStrand());
		dir.setChromosome(getChromosome());
		return dir;
	}
	
	public int strand = 0;

	public boolean overlaps(DirectedRegion anotherRegion) {
		if (anotherRegion== null)
			return false;
		if (getStrand()== 0|| anotherRegion.getStrand()== 0||
				getChromosome()== null|| anotherRegion.getChromosome()== null||
				getStart()== 0|| anotherRegion.getStart()== 0||
				getEnd()== 0 || anotherRegion.getEnd()== 0) {
			System.err.println("Not inited for overlap "+this+", "+anotherRegion);
			return false;
		}
			

		if (!(anotherRegion.getChromosome().equalsIgnoreCase(getChromosome()))||
				getStrand()!= anotherRegion.getStrand())
			return false;
		
		// includes contains
		if (this.contains(anotherRegion)|| anotherRegion.contains(this))
			return true;
		
		// overlaps
		if ((Math.abs(getStart())>= Math.abs(anotherRegion.getStart())&& 
				Math.abs(getStart())<= Math.abs(anotherRegion.getEnd()))
				|| (Math.abs(anotherRegion.getStart())>= Math.abs(getStart())
						&& Math.abs(anotherRegion.getStart())<= Math.abs(getEnd())))
			return true;
		return false;
	}
	
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
		
		if (this.start!= 0&& Math.abs(end)< Math.abs(start)) {
			int h= this.start;
			this.start= this.end;
			this.end= h;
		}
	}

	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		if (strand== 0)	// not inited 
			throw new RuntimeException("Error: set strand before start/end!");
		
		if (strand> 0&& start>= 0)
			this.start = Math.abs(start);	// force pos
		else if (strand< 0)
			this.start= -Math.abs(start);	// force neg
		
		if (this.end!= 0&& Math.abs(end)< Math.abs(start)) {
			int h= this.start;
			this.start= this.end;
			this.end= h;
		}
	}
	
	public Region getAbsoluteRegion() {
		int p1= Math.abs(getStart());
		int p2= Math.abs(getEnd());
		if (p1> p2) {
			int h= p1;
			p1= p2;
			p2= h;
		}
		DefaultRegion reg= new DefaultRegion(p1, p2);
		reg.setChromosome(getChromosome());
		return reg;
	}
	
	public int get5PrimeEdge() {
		if (isForward())
			return getStart();
		return getEnd();
	}
	
	public void set5PrimeEdge(int val) {
		if (isForward())
			setStart(val);
		else
			setEnd(val);
	}
	
	public void set3PrimeEdge(int val) {
		if (isForward())
			setEnd(val);
		else
			setStart(val);
	}
	
	public String toUCSCString() {
		return getChromosome()+":"+Math.abs(getStart())+"-"+Math.abs(getEnd());
	}
	
	public int get3PrimeEdge() {
		if (isForward())
			return getEnd();
		return getStart();
	}

	public void setStrand(int strand) {
		this.strand = strand;
	}

	public int getStrand() {
		return strand;
	}

	// copy back to DirectedRegion some day
	public static DirectedRegion[] unite(DirectedRegion[] dir) {
		Comparator compi= new AbstractRegion.PositionComparator(); 
		java.util.Arrays.sort(dir, compi);
		
		
		Vector interV= new Vector();
		int i= 0;
		DirectedRegion lastCluster= null;
		while (i < dir.length) {
			
			if (lastCluster== null) {
				if (i+1< dir.length&& dir[i].overlaps(dir[i+1])) {
					lastCluster= dir[i].merge(dir[i+1]);
					i+= 2;
				} else {
					interV.add(new DirectedRegion(dir[i++]));	// give up to merge this
				}
			} else {	// try to extend last cluster
				if (lastCluster.overlaps(dir[i])) {
					lastCluster= lastCluster.merge(dir[i]);
					++i;
				} else {
					interV.add(lastCluster);	// end cluster
					lastCluster= null;
				}
			}
			
		}
		if (lastCluster!= null)
			interV.add(lastCluster);
		
		return (DirectedRegion[]) Arrays.toField(interV);
	}

	// copy back to DirectedRegion some day
	public static DirectedRegion[] unite(DirectedRegion[][] dir2) {
		if (dir2== null)
			return null;
		int len= 0;
		for (int i = 0; i < dir2.length; i++) 			
			len+= dir2[i].length;
		DirectedRegion[] dir= new DirectedRegion[len];
		int pos= 0;
		for (int i = 0; i < dir2.length; i++) {
			for (int j = 0; j < dir2[i].length; j++) 
				dir[pos+ j]= dir2[i][j];
			pos+= dir2[i].length;
		}
		
		return unite(dir);
	}

}
