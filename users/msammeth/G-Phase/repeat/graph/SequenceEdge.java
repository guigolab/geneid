/*
 * Created on Nov 26, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.graph;

/**
 * 
 * 
 * @author micha
 */
public class SequenceEdge {

	public static void main(String[] args) {
	}
	
	protected boolean bidirectional= false;
	protected SequenceNode tailNode= null, headNode= null;
	
	public SequenceEdge(SequenceNode newTailNode, SequenceNode newHeadNode) {
		
		this.tailNode= newTailNode;
		this.headNode= newHeadNode;
	}
	
	public SequenceEdge(SequenceNode newTailNode, SequenceNode newHeadNode, boolean newBidirectional) {
		
		this(newTailNode, newHeadNode);
		this.bidirectional= newBidirectional;
	}
	
	/**
	 * @return
	 */
	public boolean isBidirectional() {
		return bidirectional;
	}
	
	public boolean isHeadNode(SequenceNode testNode) {
		
		return (testNode== headNode); 
	}

	public boolean isTailNode(SequenceNode testNode) {
		
		return (testNode== tailNode); 
	}
	
	public SequenceNode getOtherEnd(SequenceNode aNode) {
		
		if (isTailNode(aNode))
			return headNode;
		if (isHeadNode(aNode))
			return tailNode;
		return null;
	}
	
	public boolean isOutgoing(SequenceNode aNode) {
		
		if (bidirectional)
			return false;
		return isTailNode(aNode);
	}
	
	public boolean isIncoming(SequenceNode aNode) {
		
		if (bidirectional)
			return false;
		return isHeadNode(aNode);
	}
}
