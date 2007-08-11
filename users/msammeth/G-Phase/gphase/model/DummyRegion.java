/*
 * Created on Mar 6, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

/**
 * 
 * 
 * @author msammeth
 */
public class DummyRegion extends AbstractRegion {
	/* (non-Javadoc)
	 * @see gphase.model.AbstractRegion#getSpecies()
	 */
	public Species getSpecies() {
		// TODO Auto-generated method stub
		return null;
	}
	
	/* (non-Javadoc)
	 * @see gphase.model.AbstractRegion#getChromosome()
	 */
	public String getChromosome() {
		// TODO Auto-generated method stub
		return null;
	}
	
	public DummyRegion(int newStart, int newEnd) {
		this.start= newStart;
		this.end= newEnd;
	}
}
