package gphase.graph3;

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
	HashMap<Edge,Edge> outEdges= new HashMap<Edge,Edge>(2,1f);	// 2,1f
	HashMap<Edge,Edge> inEdges= new HashMap<Edge,Edge>(2,1f);	// 2,1f
	HashMap<Node, Vector<Path>> bucketHash= new HashMap<Node, Vector<Path>>(2,1f); // 2,1f
	boolean processed= false;	// coloring for contracting graph
	PartitionSet outPartitionSet, inPartitionSet;
	
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
		outEdges.put(e,e);
	}
	
	public void addInEdge(Edge e) {
		inEdges.put(e,e);
	}

	public Set<Edge> getInEdges() {
		return inEdges.keySet();
	}

	public Set<Edge> getOutEdges() {
		return outEdges.keySet();
	}
	
	public void removeInEdge(Edge e) {
		inEdges.remove(e);
	}
	
	public void removeOutEdge(Edge e) {
		outEdges.remove(e);
	}
	
	public long[] getTranscripts() {
		return transcripts;
	}

	public SpliceSite getSite() {
		return site;
	}
	
	public Collection<Vector<Path>> getBuckets() {
		return bucketHash.values(); 
	}
	
	public boolean removeBucket(Vector<Path> bucket) {
		Object o= bucketHash.remove(bucket.elementAt(0).getSource());
		if (o== null)
			return false;
		return true;
	}
	
	public String toString() {
		return getSite().toString();
	}
	
	public void addToBucket(Node v, Path p) {
		Vector<Path> vec= bucketHash.get(v);
		if (vec== null) {
			vec= new Vector<Path>();
			bucketHash.put(v, vec);
		}
		vec.add(p);
	}
	
	public void addPath(Path p) {
		Vector<Path> v= bucketHash.get(p.getSource());
		if (v== null)
			v= new Vector<Path>();
		v.add(p);
		bucketHash.put(p.getSource(),v);
	}

	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

	public PartitionSet getInPartitionSet() {
		
		if (inPartitionSet == null) {
			inPartitionSet = getPartitionSet(getInEdges());
		}

		return inPartitionSet;
	}

	private PartitionSet getPartitionSet(Set<Edge> edges) {
		Vector<long[]> psVec= new Vector<long[]>(edges.size());
		Iterator<Edge> iter= edges.iterator();
		while (iter.hasNext()) {
			Edge e= iter.next();
			psVec.add(e.getTranscripts());	// NOT: Graph.intersect(e.getTail().getTranscripts(), e.getHead().getTranscripts()), see skipped exons
		}
		return new PartitionSet(psVec);
	}
	
	public PartitionSet getOutPartitionSet() {
		if (outPartitionSet == null) {
			outPartitionSet= getPartitionSet(getOutEdges());
		}

		return outPartitionSet;
	}
	
	/*
	 * will be very inefficient
	 */
	public Vector<Vector<Path>> collectPathes(Node target, Vector<long[]> interPartVec) {
		Set<Edge> outEdges= getOutEdges();
		Vector<Vector<Path>> pathes= new Vector<Vector<Path>>(interPartVec.size());
		Iterator<Edge> iter= outEdges.iterator();
		while (iter.hasNext()) {
			Edge e= iter.next();
			Vector<long[]> inter= new Vector<long[]>(interPartVec);
			Graph.intersect(e.getTranscripts(), interPartVec, inter, false);
			Path p= new Path(e.getTail());
			p.addNode(e.getHead());
			Vector<Path> v= collectPathes(target, Graph.unite(inter), p);
		}
		
		Graph.intersect(ps1, ps2, interSet, onlyFullIntersect)
	}
}
