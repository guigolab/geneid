package gphase.graph7;

public class BubbleGraph {
	
	Bubble root;
	
	public BubbleGraph() {
		root= new Bubble();
	}
	
	
	/*
	 * parent of targetB is always the root, by iteration order
	 */
	public void insert(Bubble targetB, Bubble currentB) {
		if (currentB.isVisitedBy(targetB))
			return;
		currentB.setVisitedBy(targetB); 
		
		if (targetB.contains(currentB)) {
			targetB.addChild(currentB);
			currentB.addParent(targetB);
			return;
		}			
		if (targetB.intersects(currentB)) {
			for (int i = 0; i < currentB.getChildren().size(); i++) {
				insert(targetB, currentB.getChildren().elementAt(i));
			}
		}
	}
	
}
