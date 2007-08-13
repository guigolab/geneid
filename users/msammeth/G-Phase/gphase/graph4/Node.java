package gphase.graph4;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import gphase.model.SpliceSite;
import gphase.model.Transcript;

public class Node {
	
	public static class PositionTypeComparator extends SpliceSite.PositionTypeComparator {
		public int compare(Object o1, Object o2) {
			return super.compare(((Node) o1).getSite(),(((Node) o2).getSite()));			
		}
	}
	
	SpliceSite site;
	long[] transcripts= null;
	Vector<Edge> outEdges= new Vector<Edge>(2);	// 2,1f
	Vector<Edge> inEdges= new Vector<Edge>(2);	// 2,1f
	boolean processed= false;	// coloring for contracting graph
	HashMap<Node, HashMap<Edge, Vector<Path>>> fromNodeMap= new HashMap<Node, HashMap<Edge,Vector<Path>>>();
	
	public Node(SpliceSite newSite, long[] newTranscripts) {
		this.site= newSite;
		this.transcripts= newTranscripts;
	}
	
	@Override
	public int hashCode() {
		return getSite().getPos();
	}
	
	@Override
	public boolean equals(Object obj) {		
		return getSite().equals(((Node) obj).getSite());
	}
	
	public void addOutEdge(Edge e) {
		outEdges.add(e);
	}
	
	public void addInEdge(Edge e) {
		inEdges.add(e);
	}

	public Vector<Edge> getInEdges() {
		return inEdges;
	}

	public Vector<Edge> getOutEdges() {
		return outEdges;
	}
	
	public void removeInEdge(Edge e) {
		for (int i = 0; i < inEdges.size(); i++) {
			if (inEdges.elementAt(i)== e) {
				inEdges.remove(i);	// do not use equals()
				return;
			}
		}
	}
	
	public void removeOutEdge(Edge e) {
		for (int i = 0; i < outEdges.size(); i++) {
			if (outEdges.elementAt(i)== e) {
				outEdges.remove(i);	// do not use equals()
				return;
			}
		}
	}
	
	public long[] getTranscripts() {
		return transcripts;
	}

	public SpliceSite getSite() {
		return site;
	}
	
	public String toString() {
		return getSite().toString();
	}
	
	public void addPath(Path p, Edge e) {
		HashMap<Edge, Vector<Path>> map= fromNodeMap.get(p.getSourceNode());
		if (map== null) {
			map= new HashMap<Edge,Vector<Path>>(getInEdges().size());
			fromNodeMap.put(p.getSourceNode(),map);
		}
		Vector<Path> v= map.get(e);
		if (v== null) {
			v= new Vector<Path>(1);
			map.put(e, v);
		}
		v.add(p);
	}

	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

	public Set<Node> getFromNodes() {
		return fromNodeMap.keySet();
	}
	
	public Vector<Path> getPathesFrom(Node src) {
		
		Iterator<Vector<Path>> iter= fromNodeMap.get(src).values().iterator();
		Vector<Path> v= new Vector<Path>(2);
		while (iter.hasNext()) {
			Vector<Path> vv= iter.next();	// dont trust addAll()
			for (int i = 0; i < vv.size(); i++) {
				v.add(vv.elementAt(i));
			}
		}
		
		return v;
	}

	public HashMap<Node, HashMap<Edge, Vector<Path>>> getFromNodeMap() {
		return fromNodeMap;
	}
}
