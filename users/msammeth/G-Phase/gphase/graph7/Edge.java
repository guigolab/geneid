/**
 * test
 */
package gphase.graph7;

public class Edge {
	
	boolean contracted= false, processed= false;	
	
	public static String getStringRep(Node v, Node w) {
		return v.getSite().toString()+w.getSite().toString();
	}
	
	Node tail, head;
	String stringRep;
	long[] transcripts;
	
	
	public Edge(Node newTail, Node newHead) {
		this.tail= newTail;
		this.head= newHead;
		tail.addOutEdge(this);
		head.addInEdge(this);
	}

	public Node getHead() {
		return head;
	}

	public Node getTail() {
		return tail;
	}
	
	@Override
	public int hashCode() {		
		return toString().hashCode();
	}
	
	public String toString() {
		if (stringRep == null) 
			stringRep = getStringRep(getTail(), getHead());

		return stringRep;
	}
	
	@Override
	public boolean equals(Object obj) {
		Edge e= (Edge) obj;
		if (getTail().equals(e.getTail())&& getHead().equals(e.getHead())&& Graph.equalSet(getTranscripts(), e.getTranscripts()))
			return true;
		return false;
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}

	public boolean isContracted() {
		return contracted;
	}

	public void setContracted(boolean contracted) {
		this.contracted = contracted;
	}

	public boolean isProcessed() {
		return processed;
	}

	public void setProcessed(boolean processed) {
		this.processed = processed;
	}
}
