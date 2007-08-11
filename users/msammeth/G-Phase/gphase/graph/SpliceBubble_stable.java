package gphase.graph;

import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

import java.util.Comparator;
import java.util.HashMap;
import java.util.Stack;
import java.util.Vector;

public class SpliceBubble_stable {
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
	
	SpliceNode source= null;
	SpliceNode sink= null;
	SpliceNode[][] pathes= null;
	HashMap transHash= null;
	SpliceBubble containerBubble= null;
	SpliceBubble[] containedBubbles= null;
	
	public SpliceBubble_stable(SpliceNode src, SpliceNode snk, SplicePath[] newPathes) {
		this.source= src;
		this.sink= snk;
		setPathes(newPathes);
	}
	
	public Transcript[][] getTranscriptPartitions() {
		Transcript[][] parts= new Transcript[pathes.length][];
		Object[] o= transHash.values().toArray();
		for (int i = 0; i < o.length; i++) 
			parts[i]= (Transcript[]) o[i];
		return parts;
	}
	
	public void setPathes(SplicePath[] newPathes) {
		pathes= new SpliceNode[newPathes.length][];
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
	public void addContainedBubble(SpliceBubble blob) {
		
		blob.setContainerBubble(this);
		if (containedBubbles== null) 
			containedBubbles= new SpliceBubble[] {blob};
		else {
			containedBubbles= (SpliceBubble[]) Arrays.extendField(containedBubbles, blob);
		}
	}
	
	public void removeContainedBubble(SpliceBubble blob) {
		blob.setContainerBubble(null);
		containedBubbles= (SpliceBubble[]) Arrays.remove(containedBubbles, blob);
		if (containedBubbles.length< 1)
			containedBubbles= null;
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
	
	public boolean contains(SpliceBubble anotherBubble) {
			// positions not included
		if (anotherBubble.getSink().getSite().getPos()< getSource().getSite().getPos()
				|| anotherBubble.getSource().getSite().getPos()> getSink().getSite().getPos())
			return false;
		
			// transcript set not included
		Transcript[] t1= getTranscripts();
		Transcript[] t2= anotherBubble.getTranscripts();
		for (int i = 0; i < t2.length; i++) {	// sort ?!
			int j;
			for (j = 0; j < t1.length; j++) 
				if (t1[j]== t2[i])
					break;
			if (j== t1.length)
				return false;
		}
		return true;
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

	public SpliceBubble[] getContainedBubbles() {
		return containedBubbles;
	}

	public void setContainedBubbles(SpliceBubble[] containedBubbles) {
		this.containedBubbles = containedBubbles;
	}

	public SpliceBubble getContainerBubble() {
		return containerBubble;
	}

	public void setContainerBubble(SpliceBubble containerBubble) {
		this.containerBubble = containerBubble;
	}

	public static boolean contained(Transcript[][] containedTrans, Transcript[][] containingTrans) {
		int i;
		for (i = 0; i < containedTrans.length; i++) {
			int j;
			for (j = 0; j < containingTrans.length; j++) {
				int k;
				for (k = 0; k < containedTrans[i].length; k++) {
					int m;
					for (m = 0; m < containingTrans[j].length; m++) 
						if (containedTrans[i][k].getTranscriptID().equals(containingTrans[j][m].getTranscriptID()))
							break;
					if (m== containingTrans[j].length)	// tID of contained not found
						break;
				}
				if (k== containedTrans[i].length)	// all tID of contained found 
					break;
			}
			if (j == containingTrans.length)	// one partition of contained not found
				break;
		}
		
		return (i== containedTrans.length);
	}
}
