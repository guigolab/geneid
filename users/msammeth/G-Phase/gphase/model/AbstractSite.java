/*
 * Created on Mar 12, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.Serializable;
import java.util.Comparator;

/**
 * 
 * 
 * @author msammeth
 */
public class AbstractSite implements Serializable {
	int pos= -1;
	Gene gene= null;
	
	static final long serialVersionUID = 3169139368723074072L;
	
	public static class PositionToSpliceSiteComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			AbstractSite s= null;
			Integer i= null;
			try {
				s= (AbstractSite) arg0;
				i= (Integer) arg1;
			} catch (ClassCastException e) {
				try {
					i= (Integer) arg0;
					s= (AbstractSite) arg1;
				} catch (ClassCastException ex) {
					ex.printStackTrace();
					return -1;
				}
			}
			
			if (s.getPos()< i.intValue())
				return -1;
			if (s.getPos()> i.intValue())
				return 1;
			return 0;
		}
		
	}
	
	public static class PositionComparator implements Comparator {
	
		public int compare(Object arg0, Object arg1) {
			
			AbstractSite s1= null, s2= null;
			try {
				s1= (AbstractSite) arg0;
				s2= (AbstractSite) arg1;
			} catch (ClassCastException e) {
				e.printStackTrace();
				return -1;
			}
			
			if (s1.getPos()< s2.getPos())
				return -1;
			if (s1.getPos()> s2.getPos())
				return 1;
			return 0;
		}
		
	}

	public String toString() {
		return Integer.toString(getPos());
	}
	
	public int getPos() {
		return pos;
	}
	public void setPos(int pos) {
		this.pos = pos;
	}
	public Gene getGene() {
		return gene;
	}
	public void setGene(Gene gene) {
		this.gene = gene;
	}
}
