/*
 * Created on Nov 26, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.graph;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

/**
 * 
 * 
 * @author micha
 */
public class SequenceNode {

	public static void main(String[] args) {
	}

	protected int sequenceID= -1;
	protected int startPosition= -1;
	protected int length= -1;
	protected Hashtable inEdges= null, outEdges= null, biEdges= null; 
	
	public SequenceNode(int newSeqID, int newStartPos, int newLength) {
		
		this.sequenceID= newSeqID;
		this.startPosition= newStartPos;
		this.length= newLength;
	}
	/**
	 * @return
	 */
	public int getLength() {
		return length;
	}

	/**
	 * @return
	 */
	public int getSequenceID() {
		return sequenceID;
	}

	/**
	 * @return
	 */
	public int getStartPosition() {
		return startPosition;
	}

	public boolean addEdge(SequenceEdge newEdge) {
	
		if (newEdge.isOutgoing(this)) {
			addOutEdge(newEdge);
			return true;
		} else if (newEdge.isIncoming(this)) {
			addInEdge(newEdge);
			return true;
		} else if (newEdge.isBidirectional()) {
			addBiEdge(newEdge);
			return true;
		}
		return false;	
	}
	
	public boolean addOutEdge(SequenceEdge newOutEdge) {
	
		if (outEdges== null) 						// create new
			outEdges= new Hashtable();
			
		return addEdge(newOutEdge, outEdges);
	}

	public boolean addInEdge(SequenceEdge newInEdge) {
	
		if (inEdges== null) 						// create new
			inEdges= new Hashtable();
			
		return addEdge(newInEdge, inEdges);
	}

	public boolean addBiEdge(SequenceEdge newBiEdge) {
	
		if (biEdges== null) 						// create new
			biEdges= new Hashtable();
			
		return addEdge(newBiEdge, biEdges);
	}

	public boolean addEdge(SequenceEdge newEdge, Hashtable newHash) {
		
		if (newHash.get(newEdge)!= null)	// do not add if already in
			return false;

		newHash.put(newEdge, newEdge);		// add
		return true;
	}
	
	public Enumeration getOutEdges() {
		
		return outEdges.elements();
	}
	
}
