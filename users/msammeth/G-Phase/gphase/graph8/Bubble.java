package gphase.graph8;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class Bubble {

	HashMap<long[], Vector<long[]>> partitionMap;
	Node src, snk;
	Vector<Bubble> children;
	Vector<Bubble> parents;
	Bubble visitedBy;
	
	public Bubble(Node src, Node snk) {
		this.src= src;
		this.snk= snk;
		initPartitions();
	}
	
	Bubble() {		
	}
	
	private void initPartitions() {
		Object[] out= src.getOutEdges().toArray();
		Object[] in= snk.getInEdges().toArray();
		Edge e, f;
		for (int i = 0; i < out.length; i++) {
			e= (Edge) out[i];
			for (int j = 0; j < in.length; j++) {
				f= (Edge) in[j]; 
				long[] inter= Graph.intersect(e.getTranscripts(), f.getTranscripts());
				if (Graph.isNull(inter))
					continue;
				partitionMap.put(inter, null);
			}
		}
	}

	public Node getSrc() {
		return src;
	}

	public Node getSnk() {
		return snk;
	}

	public HashMap<long[], Vector<long[]>> getPartitionMap() {
		return partitionMap;
	}

	public boolean satisfiesSink() {
		
		Object[] inPart= snk.getInEdges().toArray();
		Object[] part= partitionMap.
		for (int i = 0; i < inPart.length; i++) {
			
		}
	}

	public Vector<Bubble> getChildren() {
		return children;
	}
	
	public void addChild(Bubble child) {
		if (children== null)
			children= new Vector<Bubble>(2,2);
		children.add(child);
	}
	
	public void addParent(Bubble parent) {
		if (parents== null)
			parents= new Vector<Bubble>(2,2);
		parents.add(parent);
	}

	public Vector<Bubble> getParents() {
		return parents;
	}
	
	public boolean contains(Bubble b) {
		if (src.getSite().getPos()<= b.getSrc().getSite().getPos()&& 
				snk.getSite().getPos()>= b.getSnk().getSite().getPos())
			return true;
		return false;
	}
	
	public boolean intersects(Bubble b) {
		if ((src.getSite().getPos()>= b.getSrc().getSite().getPos()&&
				src.getSite().getPos()<= b.getSnk().getSite().getPos())||
				(b.getSrc().getSite().getPos()>= src.getSite().getPos()&&
						b.getSrc().getSite().getPos()<= snk.getSite().getPos()))
				return true;
		return false;
	}

	public boolean isVisitedBy(Bubble b) {
		return (b== visitedBy);
	}

	public void setVisitedBy(Bubble visitedBy) {
		this.visitedBy = visitedBy; 
	}
	
	void getCombinations(int n) {
		HashMap<long[], long[]> map= new HashMap<long[], long[]>();
		Iterator<long[]> iterOuterPart= partitionMap.keySet().iterator();
		while (iterOuterPart.hasNext()) {
			long[] outerPart= iterOuterPart.next();
			
			for (int i = 0; i < getChildren().size(); i++) {
				HashMap<long[], Vector<long[]>> innerMap= getChildren().elementAt(i).getPartitionMap();
				Iterator<long[]> iterInnerPart= innerMap.keySet().iterator();
				while (iterInnerPart.hasNext()) {
					long[] innerPart= iterInnerPart.next();
					long[] inter= Graph.intersect(outerPart, innerPart);
				}
			}
		}
	}
}
