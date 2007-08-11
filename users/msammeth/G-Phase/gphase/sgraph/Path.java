package gphase.sgraph;

import gphase.model.SpliceSite;
import gphase.tools.IntVector;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class Path {
	long[] transcripts= null;
	Vector nodesAndEdges= new Vector(2,1);	// 2,2
	HashMap<Path,Path> superPathes= new HashMap<Path, Path>(2,1f);	//2,1f
	IntVector bubbleIDs= new IntVector();
	int hc= 0;
	boolean hasSpliceSite= false, lastSpliceSite= false;
	
	public Path(Node newSrc) {
		setTranscripts(newSrc.getTranscripts());
		addNode(newSrc);
	}
	
	private Path() {
		
	}

	public Vector<Node> getAllNodes() {
		Vector<Node> v= new Vector<Node>();
		for (int i = 0; i < nodesAndEdges.size(); i++) {
			if (nodesAndEdges.elementAt(i) instanceof Edge) {
				Path p= ((Edge) nodesAndEdges.elementAt(i)).getPath();
				v.addAll(p.getAllNodes());
			} else {
				v.add((Node) nodesAndEdges.elementAt(i));
			}
		}
		return v;
	}

	public void addNode(Node newNode) {
		//nodes.add(newNode);		
		nodesAndEdges.add(newNode);
		hasSpliceSite|= lastSpliceSite;
		if (!(nodesAndEdges.size()== 1))
			lastSpliceSite= newNode.getSite().isSpliceSite();
		//hc|= newNode.getSite().getPos();
	}
	
	public void addEdge(Edge newEdge) {
		nodesAndEdges.add(newEdge);
		hasSpliceSite|= lastSpliceSite;
		lastSpliceSite= newEdge.getPath().hasSpliceSite;
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
	
	public Node getSource() {
		if (nodesAndEdges.size()== 0)
			return null;
		return (Node) nodesAndEdges.elementAt(0);
	}
	
	public Node getSecondNode() {
		if (nodesAndEdges.elementAt(1) instanceof Node)
			return (Node) nodesAndEdges.elementAt(1);
		return ((Edge) nodesAndEdges.elementAt(1)).getPath().getSource();
	}
	
	public Node getSecondLastNode() {
		if (nodesAndEdges.elementAt(nodesAndEdges.size()-2) instanceof Node)
			return (Node) nodesAndEdges.elementAt(nodesAndEdges.size()-2);
		return ((Edge) nodesAndEdges.elementAt(nodesAndEdges.size()-2)).getPath().getSink();
	}
	
	public Node getSink() {
		if (nodesAndEdges.size()== 0)
			return null;
		return (Node) nodesAndEdges.elementAt(nodesAndEdges.size()-1);
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		Path newPath= new Path();
		//newPath.nodes= (Vector<Node>) nodes.clone();
		newPath.nodesAndEdges= (Vector) nodesAndEdges.clone();
		newPath.transcripts= new long[transcripts.length];
		for (int i = 0; i < transcripts.length; i++) 
			newPath.transcripts[i]= transcripts[i];
		newPath.bubbleIDs= this.bubbleIDs.cloneIntVector();
		newPath.hasSpliceSite= hasSpliceSite;
		newPath.lastSpliceSite= lastSpliceSite;
		return newPath;
	}
	
	public Path clonePath() {
		try {
			return (Path) clone();
		} catch (Exception e) {
			return null;
		}
	}
	
	@Override
	public int hashCode() {		
	return super.hashCode(); //hc;
	}
	
	@Override
	public boolean equals(Object obj) {
		return this== obj;
//		Path otherPath= (Path) obj;
//		if (length()!= otherPath.length())
//			return false;		
//		if (!getSource().equals(otherPath.getSource()))
//			return false;
//		if (!getSink().equals(otherPath.getSink()))
//			return false;
//		for (int i = 0; i < getTranscripts().length; i++) 
//			if (getTranscripts()[i]!= otherPath.getTranscripts()[i])
//				return false;
//		for (int i = 1; i < nodes.size()-1; i++) 
//			if (!nodes.elementAt(i).equals(otherPath.getNodes().elementAt(i)))
//				return false;
//		return true;
	}
	
	public Set<Path> getSuperPathes() {
		return superPathes.keySet();
	}
	
	public void addSuperPath(Path p) {
		superPathes.put(p,p);
	}
	
	public void removeSuperPath(Path p) {
		superPathes.remove(p);
	}
	
	public void removeSource() {
		assert(nodesAndEdges.elementAt(0) instanceof Node);
		nodesAndEdges.remove(0);
	}
	
	public void removeSink() {
		assert(nodesAndEdges.elementAt(nodesAndEdges.size()-1) instanceof Node);
		nodesAndEdges.remove(nodesAndEdges.size()-1);
	}
	
	public void replaceSuperPath(Path p1, Path p2) {
		superPathes.remove(p1);
		superPathes.put(p2, p2);
	}

	public int[] getBubbleIDs() {
		return bubbleIDs.toIntArray();
	}

	public void addBubbleID(int newBid) {
		this.bubbleIDs.add(newBid);
		Iterator<Path> iter= getSuperPathes().iterator();
		while (iter.hasNext())
			iter.next().addBubbleID(newBid);
	}

	public Vector getNodesAndEdges() {
		return nodesAndEdges;
	}
	
	public boolean hasSpliceSite() {
		return hasSpliceSite;
//		for (int i = 1; i < (nodesAndEdges.size()-1); i++) {
//			if (nodesAndEdges.elementAt(i) instanceof Node) {
//				if (i== 0|| i== nodesAndEdges.size()- 1)
//					continue;
//				if (((Node) nodesAndEdges.elementAt(i)).getSite().isSpliceSite())
//					return true;
//			} else {
//				Path p= ((Edge) nodesAndEdges.elementAt(i)).getPath();
//				for (int j = 0; j < p.getNodesAndEdges().size(); j++) {
//					// redundant, src and snk are to be nodes!
////					if ((i== 0&& j== 0)|| (i== nodesAndEdges.size()- 1&& j== p.getNodesAndEdges().size()))
////						continue;
//					if (((Node) p.getNodesAndEdges().elementAt(j)).getSite().isSpliceSite())
//						return true;
//				}
//			}
//		}
//		return false;
	}

}
