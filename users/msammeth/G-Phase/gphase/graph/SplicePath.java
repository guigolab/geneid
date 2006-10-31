package gphase.graph;

import gphase.model.AbstractSite;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

import java.util.Vector;

import javax.swing.tree.TreePath;

public class SplicePath {
	Vector nodeV= null;
	Transcript[] minFlowTrans= null;
	
	private SplicePath() {
	}
	
	public SplicePath(SpliceEdge newEdge) {
		nodeV= new Vector();
		nodeV.add(newEdge.getTail());	// root
		nodeV.add(newEdge.getHead());	
		minFlowTrans= newEdge.getTranscripts();
	}
	
	public Vector getNodeV() {
		return nodeV;
	}
	
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
		path.minFlowTrans= (Transcript[]) Arrays.toField(intersection);
		return path;
	}
	
protected Object clone() throws CloneNotSupportedException {
		SplicePath path= new SplicePath();
		path.nodeV= (Vector) nodeV.clone();
		path.minFlowTrans= minFlowTrans;
		return path; 
	}

public Transcript[] getTranscripts() {
	return minFlowTrans;
}
}
