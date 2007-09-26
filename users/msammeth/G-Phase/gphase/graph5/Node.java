package gphase.graph5;

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
	HashMap<Node, Vector<Path>> fromNodeMap= new HashMap<Node, Vector<Path>>();
	HashMap<Edge, Vector<Edge>> outPartitionMap= null;
	int outPartitionSize= 0;
	
	public Node(SpliceSite newSite, long[] newTranscripts) {
		this.site= newSite;
		this.transcripts= newTranscripts;
	}
	
	public Edge getOutPartition(Edge rootEdge) {
		return getOutPartitionMap().get(rootEdge);
	}
	
	public Vector<Vector<Path>> splitPathes(Vector<Path> pathes) {
		HashMap<Vector<Edge>, Vector<Path>> xRefMap= new HashMap<Vector<Edge>, Vector<Path>>(getOutPartitionSize());
		Vector<Vector<Path>> v= new Vector<Vector<Path>>(getOutPartitionSize());
		for (int i = 0; i < pathes.size(); i++) {
			Edge srcEdge= pathes.elementAt(i).getSourceEdge();
			Vector<Edge> chk= getOutPartitionMap().get(srcEdge);
			Vector<Path> pv= xRefMap.get(chk);
			if (pv== null) { 
				pv= new Vector<Path>(4);
				xRefMap.put(outPartitionMap.get(srcEdge), pv);
				v.add(pv);
			}
			pv.add(pathes.elementAt(i));
		}
		return v;
	}
	
	public HashMap<Edge, Vector<Edge>> getOutPartitionMap() {
		if (outPartitionMap == null) {
			outPartitionMap = new HashMap<Edge, Vector<Edge>>(outEdges.size());
			Iterator<Edge> iter= getOutEdges().iterator();
			while (iter.hasNext()) {
				Edge e= iter.next();
				Vector<Edge> v= new Vector<Edge>(1);
				v.add(e);
				outPartitionMap.put(e,v);
			}
			outPartitionSize= getOutEdges().size();
		}

		return outPartitionMap;
	}
	
	public void mergePartitions(Vector<Vector<Path>> splitPathes) {
				
		HashMap<Edge, Edge> map= new HashMap<Edge, Edge>(splitPathes.size(),1f);
		for (int i = 0; i < splitPathes.size(); i++) {
			map.put(splitPathes.elementAt(i).elementAt(0).getSourceEdge(),
					splitPathes.elementAt(i).elementAt(0).getSourceEdge());
		}
		
		Iterator<Edge> iter= map.keySet().iterator();
		Vector<Edge> v= new Vector<Edge>(map.values());
		while (iter.hasNext()) {
			Edge e= iter.next();
			outPartitionMap.remove(e);
			outPartitionMap.put(e, v);
		}
		
		HashMap<Vector<Edge>, Vector<Edge>> m= new HashMap<Vector<Edge>, Vector<Edge>>(outPartitionMap.size(),1f);
		Iterator<Vector<Edge>> iter2= outPartitionMap.values().iterator();
		while (iter.hasNext())
			m.put(iter2.next(), null);
		outPartitionSize= m.size();
	}
	
	public int getOutPartitionSize() {
		if (outPartitionSize== 0)
			getOutPartitionMap();
		return outPartitionSize;
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

	public HashMap<Node, Vector<Path>> getFromNodeMap() {
		return fromNodeMap;
	}
}
