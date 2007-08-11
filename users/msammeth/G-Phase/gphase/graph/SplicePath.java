package gphase.graph;

import gphase.model.AbstractSite;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.IntVector;

import java.util.Vector;

import javax.swing.tree.TreePath;

public class SplicePath {
	Vector nodeV= null;
	Vector edgeV= null;
	Transcript[] minFlowTrans= null;
	SpliceNode src= null;
	
	
	public SplicePath() {
		nodeV= new Vector();
		edgeV= new Vector();
	}
	
	
	public int getExonicLength() {
		int len= 0;
		for (int i = 0; i < edgeV.size(); i++) {
			SpliceEdge e= (SpliceEdge) edgeV.elementAt(i);
			if (e.isExonic())
				len+= e.getLength();
		}
		
		return len;
	}
	
	public int[] getExonicLengthes() {
		IntVector v= new IntVector();
		for (int i = 0; i < edgeV.size(); i++) {
			SpliceEdge e= (SpliceEdge) edgeV.elementAt(i);
			if (e.isExonic())
				v.add(e.getLength());
		}
		
		return v.toIntArray();
	}


	public int getExonicPieces() {
		int sum= 0;
		for (int i = 0; i < edgeV.size(); i++) {
			SpliceEdge e= (SpliceEdge) edgeV.elementAt(i);
			if (e.isExonic())
				++sum;
		}
		
		return sum;
	}
	
	/**
	 * creates new path from current one
	 * @param newEdge
	 * @return
	 */
	public SplicePath createPath(SpliceEdge newEdge) {
		
			// TODO: ??? chk
		if (!(newEdge.getTail().getSite() instanceof SpliceSite))	// forbid ...->AS->... transitions
			return null;
		
			// update minFlow, intersect edges
		Transcript[] newTrans= newEdge.getTranscripts();
		Vector intersection= new Vector();
		for (int i = 0; i < minFlowTrans.length; i++) 
			for (int j = 0; j < newTrans.length; j++) 
				if (minFlowTrans[i]== newTrans[j])
					intersection.add(minFlowTrans[i]);
		
		if (intersection.size()== 0)
			return null;
		
		SplicePath path= null;
		try {
			path= (SplicePath) this.clone();
		} catch (CloneNotSupportedException e) {
			System.err.println(e);
		}
		path.nodeV.add(newEdge.getHead());	// add node
		path.edgeV.add(newEdge);
		path.minFlowTrans= (Transcript[]) Arrays.toField(intersection);
		if (path.src== null)
			path.src= newEdge.getTail(); 
		return path;
	}
	
	
	public String toString() {
		String res= "[";
		for (int i = 0; i < minFlowTrans.length; i++) 
			res+= minFlowTrans[i].toString()+ ", ";
		return res.substring(0, res.length()- 2);
	}
	
	public SplicePath(SpliceEdge newEdge) {
		this();
		nodeV.add(newEdge.getTail());	// root
		nodeV.add(newEdge.getHead());	
		edgeV.add(newEdge);
		minFlowTrans= newEdge.getTranscripts();
	}
	public Vector getNodeV() {
		return nodeV;
	}
	
	public SpliceEdge[] getEdges() {
		if (edgeV== null)
			return null;
		return (SpliceEdge[]) Arrays.toField(edgeV);
	}
	/**
	 * creates new path from current one, and extends it by the new edge
	 * @param newEdge
	 * @return
	 */
	public SplicePath exendPath(SpliceEdge newEdge) {
		
			// TODO: ??? chk
		// now excluded in graph
//		if (!(newEdge.getTail().getSite() instanceof SpliceSite))	// forbid ...->AS->... transitions
//			return null;
		
			// update minFlow, intersect edges
		Transcript[] newTrans= newEdge.getTranscripts();
		Vector intersection= new Vector();
		for (int i = 0; i < minFlowTrans.length; i++) 
			for (int j = 0; j < newTrans.length; j++) 
				if (minFlowTrans[i]== newTrans[j])
					intersection.add(minFlowTrans[i]);
		
		if (intersection.size()== 0)
			return null;
		
		SplicePath path= null;
		try {
			path= (SplicePath) this.clone();
		} catch (CloneNotSupportedException e) {
			System.err.println(e);
		}
		
			//extend
		path.nodeV.add(newEdge.getHead());	// add node
		path.edgeV.add(newEdge);
		path.minFlowTrans= (Transcript[]) Arrays.toField(intersection);
		if (path.src== null)
			path.src= newEdge.getTail(); 
		return path;
	}
	
	public int edgeLength() {
		if (edgeV== null)
			return 0; 
		return edgeV.size();
	}
	
protected Object clone() throws CloneNotSupportedException {
		SplicePath path= new SplicePath();
		path.nodeV= (Vector) nodeV.clone();
		path.edgeV= (Vector) edgeV.clone();
		path.minFlowTrans= minFlowTrans;
		return path; 
	}

public Transcript[] getTranscripts() {
	return minFlowTrans;
}

public void setTranscripts(Transcript[] nuTrans) {
	minFlowTrans= nuTrans;
}


public SpliceNode getSrc() {
	return src;
}
}
