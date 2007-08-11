/*
 * Created on Nov 27, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.Serializable;

/**
 * 
 * 
 * @author msammeth
 */
public interface Region extends Serializable{
	
	public Species getSpecies();
	public int getStart();
	public int getEnd();
	public String getChromosome();
	public boolean contains(Region anotherRegion);
}
