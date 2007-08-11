/*
 * Created on Feb 21, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/**
 * 
 * 
 * @author msammeth
 */
public class SuperExon implements Serializable {
	
	HashMap exons= new HashMap();
	transient HashMap hitparade = null;
	
	public PWHit[] getHits(Gene g) {
		
		if (hitparade== null|| hitparade.get(g)== null)
			return null;
		
		Vector v= (Vector) hitparade.get(g);
		PWHit[] result= new PWHit[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (PWHit) v.elementAt(i);
		
		return result;
	}

	/**
	 * merges exons of another <code>SuperExon</code> to <code>this</code>.
	 * @param anotherSE
	 * @return
	 */
	public boolean merge(SuperExon anotherSE) {
	
		Iterator iter= anotherSE.getExons().iterator();
		while (iter.hasNext())
			add((Exon) iter.next());
		return true;
	}

	/**
	 * 
	 * @param anExon
	 * @return <code>false</code> if exon is already contained
	 */
	public boolean add(Exon anExon) {
	
		// add exon to exonlist
		if (exons.get(anExon.getStableID())== null)
			return false;
		exons.put(anExon.getStableID(), anExon);
		
		// add superexon to exon
		anExon.setSuperExon(this);
		return true;		
	}
	
	public Collection getExons() {
		return exons.values();
	}

	public boolean addHit(Gene g, PWHit h) {
			
			Vector v= null;
			if (hitparade== null)	// create new
				hitparade= new HashMap();
			else 
				v= (Vector) hitparade.get(g);	// lookup existing
			
			if (v== null)
				v= new Vector();
			else 
				for (int i = 0; i < v.size(); i++) {	// check for already added hit
					PWHit hit= (PWHit) v.elementAt(i);
					if (hit.getOtherObject(this)== h.getOtherObject(this))	// same hit between exons
						if (hit.getScore()< h.getScore()) 					// can be for superexons
							v.remove(hit);
						else
							return false;	// not added
				}
			
				// keep sorted with ascending costs
	//		int i= 0;
	//		for (i = 0; i < v.size(); i++) 
	//			if (((PWHit) v.elementAt(i)).getCost()> h.getCost())
	//				break;
	//		v.insertElementAt(h, i);	// add
			v.add(h);
			hitparade.put(g, v);
			return true;
		}
	
}
