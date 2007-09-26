package gphase.graph7;

import gphase.model.SpliceSite;
import gphase.tools.IntVector;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class Path {
	long[] transcripts= null;
	Edge sourceEdge, sinkEdge;
	
	
	
	public Path() {
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
	
	public String toString() {
		return getSourceNode().toString()+"->"+getSinkNode().toString();
	}
	
	public Node getSourceNode() {
		return sourceEdge.getTail();
	}
	
	public Edge getSourceEdge() {
		return sourceEdge;
	}
	
	public Edge getSinkEdge() {
		return sinkEdge;
	}
	
	public Node getSinkNode() {
		return sinkEdge.getHead();
	}

	public void setSourceEdge(Edge sourceEdge) {
		this.sourceEdge = sourceEdge;
	}

	public void setSinkEdge(Edge sinkEdge) {
		this.sinkEdge = sinkEdge;
	}

}
