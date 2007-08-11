package gphase.graph;

import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Stack;
import java.util.Vector;

public class SpliceBubble {
	
	boolean iBubble= false;
	
	public static class PositionComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			SpliceBubble bub0= (SpliceBubble) arg0;
			SpliceBubble bub1= (SpliceBubble) arg1;
			if (bub0.getSource().getSite().getPos()< bub1.getSource().getSite().getPos())
				return -1;
			if (bub0.getSource().getSite().getPos()> bub1.getSource().getSite().getPos())
				return 1;
			if (bub0.getSink().getSite().getPos()< bub1.getSink().getSite().getPos())	
				return -1;
			if (bub0.getSink().getSite().getPos()> bub1.getSink().getSite().getPos())
				return 1;
			return 0;
		}
	}
	
	public static class SizePartitionSizeComparator extends SpliceBubble.SizeComparator {
		public int compare(Object o1, Object o2) {
			int val= super.compare(o1, o2);
			if (val!= 0)
				return val;
			SpliceBubble bub1= (SpliceBubble) o1;
			SpliceBubble bub2= (SpliceBubble) o2;
			Transcript[][] t1= bub1.getTranscriptPartitions();
			Transcript[][] t2= bub2.getTranscriptPartitions();
			
			//if (t1.length> t2.length)
			if (containedPreservesPartitions(t2, t1))
				return 1;
			//if (t2.length> t1.length)
			if (containedPreservesPartitions(t1, t2))
				return -1;
			return 0;
		}
	}
	
	public static class SizeComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			SpliceBubble bub1= (SpliceBubble) o1;
			SpliceBubble bub2= (SpliceBubble) o2;
			
			int size1= bub1.getSize();
			int size2= bub2.getSize();
			
			if (size1< size2)
				return -1;
			if (size1> size2)
				return 1;
			return 0;
		}
	}

	public static class SizeASComparator extends SpliceBubble.SizeComparator {
		public int compare(Object o1, Object o2) {
			int val= super.compare(o1, o2);
			if (val!= 0)
				return val;
			SpliceBubble bub1= (SpliceBubble) o1;
			SpliceBubble bub2= (SpliceBubble) o2;
			boolean src1SS= (bub1.getSource().getSite() instanceof SpliceSite);
			boolean snk1SS= (bub1.getSink().getSite() instanceof SpliceSite);
			boolean src2SS= (bub2.getSource().getSite() instanceof SpliceSite);
			boolean snk2SS= (bub2.getSink().getSite() instanceof SpliceSite);
			Transcript[][] t1= bub1.getTranscriptPartitions();
			Transcript[][] t2= bub2.getTranscriptPartitions();
			
//			if (SpliceBubble.contained(t2, t1))
//				return 1;
//			if (SpliceBubble.contained(t1, t2))
//				return -1;
			
				// AS bubble is bigger
			if ((!src1SS)&& (!snk1SS))
				return 1;
			if ((!src2SS)&& (!snk2SS))
				return -1;
			
			if (((!src1SS)&& (src2SS)&& (!snk2SS))||
					((!snk1SS)&& (snk2SS)&& (!src2SS)))
				return 1;
			if (((!src2SS)&& (src1SS)&& (!snk1SS))||
					((!snk2SS)&& (snk1SS)&& (!src1SS)))
				return -1;
			
			assert(false);
			return 0;
		}
	}

	SpliceNode source= null;
	SpliceNode sink= null;
	SpliceNode[][] pathes= null;
	SplicePath[] savPathes= null;
	HashMap transHash= null; // maps SpliceNode[]#Transcript[]
	SpliceBubble[] parents= null;
	SpliceBubble[] children= null;
	/**
	 * @deprecated now children, parents
	 */
	SpliceBubble containerBubble= null;
	/**
	 * @deprecated now children, parents
	 */
	SpliceBubble[] containedBubbles= null;
	
	public SpliceBubble(SpliceNode src, SpliceNode snk, SplicePath[] newPathes) {
		this.source= src;
		this.sink= snk;
		setPathes(newPathes);
	}
	
	public boolean hasParents() {
		if (parents== null|| parents.length< 1)
			return false;
		return true;
	}
	
	public boolean hasChildren() {
		if (children== null|| children.length< 1)
			return false;
		return true;
	}
	
	public Transcript[][] getTranscriptPartitions() {
		Transcript[][] parts= new Transcript[pathes.length][];
		Object[] o= transHash.values().toArray();
		for (int i = 0; i < o.length; i++) 
			parts[i]= (Transcript[]) o[i];
		return parts;
	}
	
	protected Object clone() throws CloneNotSupportedException {
		
		SplicePath[] cPathes= null;
		if (savPathes!= null) {
			cPathes= new SplicePath[savPathes.length];
			for (int i = 0; i < cPathes.length; i++) {
				cPathes[i]= (SplicePath) savPathes[i].clone();
			}
		}
		
		SpliceBubble clone= new SpliceBubble(getSource(), getSink(), cPathes);
		return clone;
	}
	
	public void setPathes(SplicePath[] newPathes) {
		pathes= new SpliceNode[newPathes.length][];
		savPathes= newPathes;
		transHash= new HashMap();
		for (int i = 0; i < newPathes.length; i++) {
			Vector v= newPathes[i].getNodeV();
			if (v.size()< 2)
				pathes[i]= new SpliceNode[0];
			else {
				pathes[i]= new SpliceNode[v.size()- 2];	// omit src and snk
				for (int j = 1; j < v.size()- 1; j++) 
					pathes[i][j-1]= (SpliceNode) v.elementAt(j);
			}
			transHash.put(pathes[i], newPathes[i].getTranscripts());
		}
	}
	/**
	 * @deprecated now children, parents
	 */
	public void addContainedBubble(SpliceBubble blob) {
		
		blob.setContainerBubble(this);
		if (containedBubbles== null) 
			containedBubbles= new SpliceBubble[] {blob};
		else {
			containedBubbles= (SpliceBubble[]) Arrays.extendField(containedBubbles, blob);
		}
	}
	
	public SpliceNode getSource() {
		return source;
	}
	public SpliceNode getSink() {
		return sink;
	}
	
	public int getPathSetSize() {
		return ((Integer) sink.getFromList().get(source.getSite())).intValue();
	}
	
	public SpliceNode[][] getPathes() {
		return pathes;
	}
	
	public SpliceNode[] removePath(Transcript[] part) {
		Object[] o= transHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) {
			Transcript[] partmp= (Transcript[]) transHash.get(o[i]);
			if (identical(partmp, part)) {
				transHash.remove(o[i]);

				//pathes= (SpliceNode[][]) Arrays.remove(pathes, o[i]);
				SpliceNode[][] nuPathes= new SpliceNode[pathes.length- 1][];
				int pos= 0;
				for (int j = 0; pathes!= null&& j < pathes.length; j++) {
					if (identical(pathes[j], (SpliceNode[]) o[i]))
						continue;
					nuPathes[pos++]= pathes[j];
				}
				pathes= nuPathes;
				
				for (int j = 0; j < savPathes.length; j++) {	// TODO second search ineficient, make data structure better
					if (identical(part, savPathes[j].getTranscripts())) {
						savPathes= (SplicePath[]) Arrays.remove(savPathes, savPathes[j]);
						break;
					}
				}
				return (SpliceNode[]) o[i];
			} else if (contained(partmp, part)) {
				transHash.remove(o[i]);
				Vector v= new Vector(partmp.length);
				for (int j = 0; j < partmp.length; j++) {
					int k;
					for (k = 0; k < part.length; k++) {
						if (partmp[j]== part[k])
							break;
					}
					if (k== part.length)
						v.add(partmp[j]);	// not found in remove array
				}
				Transcript[] nuPart= (Transcript[]) Arrays.toField(v);
				transHash.put(o[i], nuPart);
				
				// pathes hasnt to be changed, no path is completely lost
				
					// correct savPathes
				for (int j = 0; j < savPathes.length; j++) {
					if (identical(partmp, savPathes[j].getTranscripts())) {
						savPathes[j].setTranscripts(nuPart);
						break;	// can only happen once
					}
				}
			}
		}
		return null;
	}

	public SpliceNode[][] getPathes_old() {
		if (pathes == null) {
			Vector v= new Vector();	// for result
			Stack backtrackStack= new Stack();
			backtrackStack.push(source);

			getPath(backtrackStack, new Vector(), v);
			
			pathes= new SpliceNode[v.size()][];
			for (int i = 0; i < pathes.length; i++) {
				pathes[i]= (SpliceNode[]) Arrays.toField(v.elementAt(i));
				if (pathes[i]== null)
					pathes[i]= new SpliceNode[0];
				transHash.put(pathes[i], transHash.remove(v.elementAt(i)));
			}
		}

		return pathes;
	}
	
	int getPartitionNr(Transcript t) {
		Object[] o= transHash.values().toArray();
		for (int i = 0; i < o.length; i++) { 
			Transcript[] tr= (Transcript[]) o[i];
			for (int j = 0; j < tr.length; j++) {
				if (tr[j]== t)
					return i;
			}
		}
		return -1;
	}
	
	void getPath(Stack backtrackStack, Vector path, Vector v) {
		
		if (backtrackStack.isEmpty())
			return;
		
		SpliceNode parent= (SpliceNode) backtrackStack.pop();
		if (parent== source)
			transHash= new HashMap();
		SpliceEdge[] outEdges= parent.getOutEdges();
		for (int i = 0; i < outEdges.length; i++) {
			SpliceNode tmpNode= outEdges[i].getHead();
			if (tmpNode== sink) {
				v.add(path);	// do not add sink to path
				transHash.put(path, outEdges[i].getTranscripts());	// store transcripts for path
				continue;
			}
			Vector newPath= (Vector) path.clone();
			newPath.add(tmpNode);
			backtrackStack.push(tmpNode);	// else
			getPath(backtrackStack, newPath, v);
		}
		
	}
	
	public Transcript[] getTranscripts_old() {
			// intersect outgoing from src with incoming from sink
		SpliceEdge[] edges= source.getOutEdges();
		Vector outV= new Vector();
		for (int i = 0; i < edges.length; i++) 
			for (int j = 0; j < edges[i].getTranscripts().length; j++) 
				outV.add(edges[i].getTranscripts()[j]);	// doubles impossible
		edges= sink.getInEdges();
		Vector inV= new Vector();
		for (int i = 0; i < edges.length; i++) 
			for (int j = 0; j < edges[i].getTranscripts().length; j++) 
				inV.add(edges[i].getTranscripts()[j]);	// doubles impossible
		
		Vector resultV= new Vector();
		for (int i = 0; i < outV.size(); i++) 	// sort for more efficient ?!
			for (int j = 0; j < inV.size(); j++) {
				if (i>= outV.size())
					break;
				if (outV.elementAt(i)== inV.elementAt(j)) {
					resultV.add(outV.remove(i));
					inV.remove(j);
				}
			}
			
		return (Transcript[]) Arrays.toField(resultV);
	}
	
	public String toString() {
		
		return "["+source+" ==> "+sink+ "]";
	}
	
	public boolean comprises(SpliceBubble anotherBubble) {
		if (anotherBubble.getSink().getSite().getPos()== getSource().getSite().getPos()
				&& anotherBubble.getSource().getSite().getPos()== getSink().getSite().getPos())
			return true;	// check for transcript set of anotherBubble is smaller or equal ?!
		return false;
	}

	public boolean containsGeometrically(SpliceBubble anotherBubble) {
		// src/snk changed on neg strand, dont need abs()
		if (anotherBubble.getSource().getSite().getPos()>= getSource().getSite().getPos()&&
				anotherBubble.getSink().getSite().getPos()<= getSink().getSite().getPos())
			return true;
		return false;
	}
	
	public boolean overlapsGeometrically(SpliceBubble anotherBubble) {
		// src/snk changed on neg strand, dont need abs()
		if (anotherBubble.getSink().getSite().getPos()< getSource().getSite().getPos()
				|| anotherBubble.getSource().getSite().getPos()> getSink().getSite().getPos())
			return false;
		return true;
	}
	
	public boolean contains(SpliceBubble anotherBubble) {
		
		if (!containsGeometrically(anotherBubble))
			return false;
		
			// transcript set not included
		Transcript[][] t1= getTranscriptPartitions();
		Transcript[][] t2= anotherBubble.getTranscriptPartitions();
		if (t2.length> t1.length)
			return false;
		if (SpliceBubble.containedByTranscript(t2, t1))
			return true;
		return false;
	}

	public boolean containsPreservesPartitions(SpliceBubble anotherBubble) {
		
		if (!containsGeometrically(anotherBubble))
			return false;
		
			// transcript set not included
		Transcript[][] t1= getTranscriptPartitions();
		Transcript[][] t2= anotherBubble.getTranscriptPartitions();
		if (SpliceBubble.preservesPartitions(t1, t2))
			return true;
		return false;
	}
	
	public void addParent(SpliceBubble newParent) {
		parents= (SpliceBubble[]) Arrays.addUnique(parents, newParent);
	}
	
	public void removeParent(SpliceBubble toBremoved) {
		parents= (SpliceBubble[]) Arrays.remove(parents, toBremoved);
	}
	
	public void removeChild(SpliceBubble toBremoved) {
		children= (SpliceBubble[]) Arrays.remove(children, toBremoved);
	}
	
	public void addChild(SpliceBubble newChild) {
		children= (SpliceBubble[]) Arrays.addUnique(children, newChild);
	}
	// -1 s1 transcript set included in s2, 0 equal, +1 s2 tset included in s1, 2 none included
	public static int compareTranscriptSets_old(SpliceBubble s1, SpliceBubble s2) {
		// transcript set not included
		Transcript[] t1= s1.getTranscripts();
		Transcript[] t2= s2.getTranscripts();
		boolean ovlp= false;
		for (int i = 0; i < t2.length; i++) {	// sort ?!
			int j;
			for (j = 0; j < t1.length; j++) { 
				if (t1[j]== t2[i])
					break;
				else
					ovlp= true;
			}
			if (j== t1.length) {
				if (ovlp)
					return 2;
				else
					return -1;
			}
		}
		ovlp= false;
		for (int i = 0; i < t1.length; i++) {	// sort ?!
			int j;
			for (j = 0; j < t2.length; j++) { 
				if (t2[j]== t1[i])
					break;
				else
					ovlp= true;
			}
			if (j== t2.length) {
				if (ovlp)
					return 2;
				else
					return 1;
			}
		}
		return 0;
	}
	
	public HashMap getTransHash() {
		return transHash;
	}

	public int getSize() {
		return sink.getSite().getPos()- 
			source.getSite().getPos();
	}

	/**
	 * @deprecated now children, parents
	 */
	public void setContainerBubble(SpliceBubble containerBubble) {
		this.containerBubble = containerBubble;
	}

	public boolean equals(Object obj) {
		SpliceBubble b= null;
		try {
			b= (SpliceBubble) obj;
		} catch (ClassCastException e) {
			return super.equals(obj);
		}

		if (b.getSource()== getSource()&& b.getSink()== getSink()&&
				identical(getTranscriptPartitions(), b.getTranscriptPartitions())) {
			return true;
		} else
			return false;
	}
	
	
	public static boolean intersect_old(SpliceBubble bub0, SpliceBubble bub1, Vector chkBubV) {
		
			// find wider bubble
		SpliceBubble wideBub= null, tinyBub= null;
		if (bub0.getSource().getSite().getPos()<= bub1.getSource().getSite().getPos()
				&& bub0.getSink().getSite().getPos()>= bub1.getSink().getSite().getPos()) {
			wideBub= bub0;
			tinyBub= bub1;
		} else if (bub1.getSource().getSite().getPos()<= bub0.getSource().getSite().getPos()
				&& bub1.getSink().getSite().getPos()>= bub0.getSink().getSite().getPos()) {
			wideBub= bub1;
			tinyBub= bub0;
		} else
			return false;	// geometrically not overlapping
		
			// create intersection bubble
		Transcript[][] tinyPart= tinyBub.getTranscriptPartitions();
		Transcript[][] widePart= wideBub.getTranscriptPartitions();
		Vector v= new Vector();	// collects transcripts found in the intersection
		for (int i = 0; i < tinyPart.length; i++) {
			int j;
			for (j = 0; j < widePart.length; j++) {
				Transcript[][] iPart= intersect(tinyPart[i], widePart[j]);	// [0] sobre tiny, [1] sobre wide, [2] intersect
				if (iPart!= null&& iPart[2]!= null) 
						v= (Vector) Arrays.addUnique(v, iPart[2]);	// partitions a sobre in tiny
			}
		}
		Vector vv= new Vector();	// collect transcript NOT found in any intersection
		for (int i = 0; i < tinyPart.length; i++) 
			for (int j = 0; j < tinyPart[i].length; j++) {
				int k;
				for (k = 0; k < v.size(); k++) 
					if (v.elementAt(k)== tinyPart[i][j])
						break;
				if (k== v.size())
					vv.add(tinyPart[i][j]);
			}
		
		SpliceBubble interBub= null;
		try {
			interBub = (SpliceBubble) tinyBub.clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		for (int i = 0; i < vv.size(); i++) 
			interBub.removePath(new Transcript[] {(Transcript) vv.elementAt(i)});	
					// works on [], not very efficient
		int x;
		for (x = 0; x < chkBubV.size(); x++) 
			if (interBub.equals(chkBubV.elementAt(x))) {
				interBub= (SpliceBubble) chkBubV.elementAt(x);
				break;
			}
		if (x== chkBubV.size())
			chkBubV.add(interBub);
		
			// insert into hierarchy
		if (vv.size()> 0&& interBub.getPathes()!= null&& interBub.getPathes().length> 0) {
			SpliceBubble[] c= wideBub.getChildren();
			for (int i = 0; c!= null&& i < c.length; i++) 
				if (interBub.contains(c[i])) {
					try {
						c[i].removeParent(wideBub);
					} catch (Exception e) {;}
					interBub.addChild(c[i]);
					c[i].addParent(interBub);
					wideBub.removeChild(c[i]);
				}
			wideBub.addChild(interBub);
			interBub.addParent(wideBub);
			
			c= tinyBub.getChildren();
			for (int i = 0; c!= null&& i < c.length; i++) 
				if (interBub.contains(c[i])) {
					try {
						c[i].removeParent(tinyBub);
					} catch (Exception e) {;}
					interBub.addChild(c[i]);
					c[i].addParent(interBub);
					tinyBub.removeChild(c[i]);
				}
			tinyBub.addChild(interBub);
			interBub.addParent(tinyBub);
		} else 
			return false;
		
		return true;
	}
	
	public static boolean identical(SpliceNode[] p0, SpliceNode[] p1) {
		if ((p0== null&& p1!= null)|| (p0!= null&& p1== null))
			return false;
		if (p0== null&& p1== null) 
			return true;
		
		if (p0.length!= p1.length)
			return false;
		for (int i = 0; i < p1.length; i++) {
			if (p0[i]!= p1[i])	// SpliceNode identity maintained
				return false;
		}
		return true;
	}
		
	public static boolean identical(Transcript[][] part0, Transcript[][] part1) {
		if (part0.length!= part1.length)
			return false;
		
		for (int i = 0; i < part0.length; i++) {
			int j;
			for (j = 0; j < part1.length; j++) {
				if (identical(part0[i], part1[j]))
					break;
			}
			if (j== part1.length)	// not found
				return false;
		}
		return true;
	}
	
	public static boolean identical(Transcript[] part0, Transcript[] part1) {
		
		if (part0== null|| part1== null)
			return false;
		
		if (part0.length!= part1.length)
			return false;
		
		for (int i = 0; i < part0.length; i++) {
			int j;
			for (j = 0; j < part1.length; j++) {
				if (part0[i]== part1[j])
					break;
			}
			if (j== part1.length)
				return false;	// 1 trans not found
		}
		
		return true;
	}
	
	/**
	 * 
	 * @param part0
	 * @param part1
	 * @return <code>true</code> iff there is at least one transcript in both
	 * partitions
	 */
	public static boolean intersects(Transcript[][] part0, Transcript[][] part1) {
		for (int i = 0; i < part0.length; i++) {
			for (int j = 0; j < part1.length; j++) {
				boolean b= intersects(part0[i], part1[j]);
				if (b)
					return true;
			}
		}
		return false;
	}
	
	/**
	 * 
	 * @param part0
	 * @param part1
	 * @return <code>true</code> iff there is at least one transcript in both
	 * partitions
	 */
	public static boolean intersects(Transcript[] part0, Transcript[] part1) {
		for (int i = 0; i < part0.length; i++) {
			for (int j = 0; j < part1.length; j++) {
				if (part0[i]== part1[j])
					return true;
			}
		}
		return false;
	}
	/**
	 * 
	 * @param part0
	 * @param part1
	 * @return <code>Transcript[3][]</code> (sobre in part0, sobre in part1, intersection)
	 * 		   iff there was an intersection possible,
	 * 		   <code>null</code> else 
	 */
	public static Transcript[][] intersect(Transcript[] part0, Transcript[] part1) {
		if (part0== null|| part1== null) {
			System.err.println("assertion failed: empty transcript set partition");
			return null;
		}
		Vector v0= new Vector();
		Vector v2= new Vector();
		for (int i = 0; i < part0.length; i++) {
			int j;
			for (j = 0; j < part1.length; j++) {
				if (part0[i]== part1[j]) {
					v2.add(part0[i]);
					break;
				}
			}
			if (j== part1.length)	// part0[i] not found
				v0.add(part0[i]);
		}
		
		if (v2.size()< 1)
			return null;	// no intersection
		
		Vector v1= new Vector();
		for (int i = 0; i < part1.length; i++) {
			int j;
			for (j = 0; j < v2.size(); j++) {
				if (v2.elementAt(j)== part1[i])
					break;
			}
			if (j== v2.size())
				v1.add(part1[i]);	// sobre in part1
		}
		
		Transcript[][] result= new Transcript[][] 
		         {(Transcript[]) Arrays.toField(v0), (Transcript[]) Arrays.toField(v1), (Transcript[]) Arrays.toField(v2)};
		return result;
	}
	
	public static boolean contained(Transcript[][] containedTrans, Transcript[][] containingTrans) {
		int i;
		for (i = 0; i < containedTrans.length; i++) {
			int j;
			for (j = 0; j < containingTrans.length; j++) {
				if (contained(containedTrans[i], containingTrans[j]))
					break;
			}
			if (j == containingTrans.length)	// one partition of contained not found
				break;
		}
		
		return (i== containedTrans.length);
	}
	
	public static boolean containedByTranscript(Transcript[][] containedTrans, Transcript[][] containingTrans) {
		int i;
		for (i = 0; i < containedTrans.length; i++) {
			int j;
			for (j = 0; j < containedTrans[i].length; j++) {
				int k;
				for (k = 0; k < containingTrans.length; k++) {
					int m;
					for (m = 0; m < containingTrans[k].length; m++) 
						if (containedTrans[i][j]== containingTrans[k][m])
							break;
					if (m< containingTrans[k].length)	// found
						break;
				}
				if (k== containingTrans.length)
					return false;
			}
		}
		
		return true;
	}
	
	/**
	 * all partitions matched by at least 1 transcript support.
	 * no partition matches 2 of the other
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static boolean preservesPartitions(Transcript[][] t1, Transcript[][] t2) {
	
		for (int i = 0; i < t1.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t1[i].length; j++) {
				for (int k = 0; k < t2.length; k++) {
					for (int m = 0;m < t2[k].length; m++) {
						if (t1[i][j]== t2[k][m]) {
							if (found== null)
								found= t2[k];
							else if (found!= t2[k])
								return false;
						}
					}
				}
			}
		}
		
		return true;
	}

	/**
	 * t1 shows in each partition at least one tID of each partition of t2 (exclusively?).
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static boolean atLeastOneContained(Transcript[][] t1, Transcript[][] t2) {

		for (int i = 0; i < t1.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t1[i].length; j++) {
				for (int k = 0; k < t2.length; k++) {
					for (int m = 0;m < t2[k].length; m++) {
						if (t1[i][j]== t2[k][m]) {
							if (found== null)
								found= t2[k];
							else if (found!= t2[k])
								return false;
						}
					}
				}
			}
			if (found== null)
				return false;	// all tIDs of a partition not found
		}
		
		return true;
	}

	/**
	 * all partitions matched by at least 1 transcript support.
	 * no partition matches 2 of the other
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static boolean preservesAllPartitions(Transcript[][] t1, Transcript[][] t2) {
		if (t1.length!= t2.length)
			return false;
		
		for (int i = 0; i < t1.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t1[i].length; j++) {
				for (int k = 0; k < t2.length; k++) {
					for (int m = 0;m < t2[j].length; m++) {
						if (t1[i][j]== t2[k][m]) {
							if (found== null)
								found= t2[k];
							else if (found!= t2[k])
								return false;
						}
					}
				}
			}
			if (found== null)
				return false;
		}
		
			// necessary
		for (int i = 0; i < t2.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t2[i].length; j++) {
				for (int k = 0; k < t1.length; k++) {
					for (int m = 0;m < t1[j].length; m++) {
						if (t1[i][j]== t1[k][m]) {
							if (found== null)
								found= t1[k];
							else if (found!= t1[k])
								return false;
						}
					}
				}
			}
			if (found== null)
				return false;
		}
		
		return true;
	}

	/**
	 * all transcripts of t1 contained in t2 AND partitions only narrowed t2->t1, eg. [5],[3] -> [5,3].
	 * no partition matches 2 of the other
	 * @param t1
	 * @param t2
	 * @return
	 */
	public static boolean containedPreservesPartitions(Transcript[][] t1, Transcript[][] t2) {
		for (int i = 0; i < t1.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t1[i].length; j++) {
				for (int k = 0; k < t2.length; k++) {
					for (int m = 0;m < t2[j].length; m++) {
						if (t1[i][j]== t2[k][m]) {
							if (found== null)
								found= t2[k];
							else if (found!= t2[k])
								return false;
						}
					}
				}
				if (found== null)
					return false;
			}
		}
		
			// necessary
		for (int i = 0; i < t2.length; i++) {
			Transcript[] found= null;
			for (int j = 0; j < t2[i].length; j++) {
				for (int k = 0; k < t1.length; k++) {
					for (int m = 0;m < t1[j].length; m++) {
						if (t1[i][j]== t1[k][m]) {
							if (found== null)
								found= t1[k];
							else if (found!= t1[k])
								return false;
						}
					}
				}
			}
			if (found== null)
				return false;
		}
		
		return true;
	}

	public static boolean contained(Transcript[] containedTrans, Transcript[] containingTrans) {
		int k;
		for (k = 0; k < containedTrans.length; k++) {
			int m;
			for (m = 0; m < containingTrans.length; m++) 
				if (containedTrans[k].getTranscriptID().equals(containingTrans[m].getTranscriptID()))
					break;
			if (m== containingTrans.length)	// tID of contained not found
				break;
		}
		if (k== containedTrans.length)	// all tID of contained found 
			return true;
		return false;
	}

	public SpliceBubble[] getChildren() {
		return children;
	}

	public SpliceBubble[] getParents() {
		return parents;
	}
	
	Vector getAncestors(Vector v, Vector resV) {

		v.add(this);
		
		if (!hasParents()) {
			//resV.add(v);
			return resV;	
		}
		
		for (int i = 0; i < getParents().length; i++) {
			Vector vv= (Vector) v.clone();
			resV.add(vv);
			getParents()[i].getAncestors(vv, resV);
		}
		return resV;
	}
	
	public SpliceBubble[][] getAncestors() {
		Vector resV= new Vector();
		Vector v= new Vector();
		resV.add(v);
		resV= getAncestors(v, resV);
		return (SpliceBubble[][]) Arrays.toField(resV); 
	}
	
	public void setParents(SpliceBubble[] parents) {
		this.parents = parents;
	}

	public boolean isIBubble() {
		return iBubble;
	}

	public void setIBubble(boolean bubble) {
		iBubble = bubble;
	}

	public SplicePath[] getSavPathes() {
		return savPathes;
	}
}
