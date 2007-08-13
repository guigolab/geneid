package gphase.graph3;

public class Edge {
	public static String getStringRep(Node v, Node w) {
		return v.getSite().toString()+w.getSite().toString();
	}
	
	Path path= null;
	Node tail, head;
	String stringRep;
	int hc= 0;
	long[] transcripts;
	
	
	public Edge(Node newTail, Node newHead) {
		this.tail= newTail;
		this.head= newHead;
		tail.addOutEdge(this);
		head.addInEdge(this);
	}

	public Path getPath() {
		return path;
	}

	public void setPath(Path path) {
		this.path = path;
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
	
	public int getHC() {
		if (hc == 0) {
			hc = getTail().getSite().getPos()+ getHead().getSite().getPos();
		}

		return hc;
	}

	@Override
	public boolean equals(Object obj) {
		Edge e= (Edge) obj;
		if (getTail().equals(e.getTail())&& getHead().equals(e.getHead()))
			return true;
		return false;
	}

	public long[] getTranscripts() {
		return transcripts;
	}

	public void setTranscripts(long[] transcripts) {
		this.transcripts = transcripts;
	}
}
