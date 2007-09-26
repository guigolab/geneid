package gphase.graph6;
/**
 * invested a lot in new graph contraction (now 2 edges between two vertexes can exist)
 * mergePartitions not correct -- cannot merge partitions has to save n-combinations realized
 * time benchmark ok -> 15 sec RefSeq, <5min mRNA
 */

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
	private static HashMap<Edge, Edge> mapEdgeEdge;
	private static Iterator<Edge> iterEdge;
	private static Iterator<Vector<Edge>> iterVecEdge;
	private static HashMap<Vector<Edge>, Vector<Edge>> mapVecEdgeVecEdge;
	private static HashMap<Vector<Edge>, Vector<Path>> mapVecEdgeVecPath;
	private static Vector<Vector<Path>> vecVecPath;	
	
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
	
	public Vector<Vector<Path>> splitPathes(Vector<Path> pathes) {
		mapVecEdgeVecPath= new HashMap<Vector<Edge>, Vector<Path>>(getOutPartitionSize());
		vecVecPath= new Vector<Vector<Path>>(getOutPartitionSize());
		for (int i = 0; i < pathes.size(); i++) {
			Edge srcEdge= pathes.elementAt(i).getSourceEdge();
			Vector<Edge> chk= getOutPartitionMap().get(srcEdge);
			Vector<Path> pv= mapVecEdgeVecPath.get(chk);
			if (pv== null) { 
				pv= new Vector<Path>(4);
				mapVecEdgeVecPath.put(outPartitionMap.get(srcEdge), pv);
				vecVecPath.add(pv);
			}
			pv.add(pathes.elementAt(i));
		}
		return vecVecPath;
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
				
		mapEdgeEdge= new HashMap<Edge, Edge>(splitPathes.size(),1f);
		for (int i = 0; i < splitPathes.size(); i++) {
				mapEdgeEdge.put(splitPathes.elementAt(i).elementAt(0).getSourceEdge(),
						splitPathes.elementAt(i).elementAt(0).getSourceEdge());
		}
		
		iterEdge= mapEdgeEdge.keySet().iterator();
		Vector<Edge> v= new Vector<Edge>(mapEdgeEdge.values());
		while (iterEdge.hasNext()) {
			Edge e= iterEdge.next();
			outPartitionMap.remove(e);
			outPartitionMap.put(e, v);
		}
		
		mapVecEdgeVecEdge= new HashMap<Vector<Edge>, Vector<Edge>>(outPartitionMap.size(),1f);
		iterVecEdge= outPartitionMap.values().iterator();
		while (iterVecEdge.hasNext())
			mapVecEdgeVecEdge.put(iterVecEdge.next(), null);
		outPartitionSize= mapVecEdgeVecEdge.size();
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
	
	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

	public Set<Node> getFromNodes() {
		return fromNodeMap.keySet();
	}
	
	public HashMap<Node, Vector<Path>> getFromNodeMap() {
		return fromNodeMap;
	}

	public void setFromNodeMap(HashMap<Node, Vector<Path>> fromNodeMap) {
		this.fromNodeMap = fromNodeMap;
	}
}
