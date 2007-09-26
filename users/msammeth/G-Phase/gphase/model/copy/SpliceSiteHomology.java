/*
 * Created on Apr 18, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.copy;

/**
 * 
 * 
 * @author msammeth
 */
public class SpliceSiteHomology {

	SpliceSite splice1= null;
	SpliceSite splice2= null;
	
	float percID= -1f;
	
	String[] ali= null;
	int cost= -1;

	public SpliceSiteHomology(SpliceSite s1, SpliceSite s2) {
		splice1= s1;
		splice2= s2;
	}
	public float getPercID() {
		if ((percID< 0f)&& (ali!= null)) {
			int ctr= 0;
			for (int i = 0; i < ali[0].length(); i++) 
				if (ali[0].charAt(i)== ali[1].charAt(i))
					++ctr;
			percID= (float) ctr/ (float) ali[0].length();
		}
		return percID;
	}
	public SpliceSite getSplice1() {
		return splice1;
	}
	public SpliceSite getSplice2() {
		return splice2;
	}
	public String[] getAlignment() {
		return ali;
	}
	public void setAlignment(String[] ali) {
		this.ali = ali;
	}
	public int getCost() {
		return cost;
	}
	public void setCost(int cost) {
		this.cost = cost;
	}
}
