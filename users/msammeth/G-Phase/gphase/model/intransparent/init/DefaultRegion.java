/*
 * Created on Nov 27, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.intransparent.init;

import java.io.Serializable;

/**
 * 
 * 
 * @author msammeth
 */
public class DefaultRegion extends AbstractRegion implements Region {

	static long serialVersionUID= 6558964245030911377l;
	String chromosome = null;
	Species species = null;
	double score= 0d;
	
	public boolean contains(SpliceSite ss) {
		if (!ss.getGene().getChromosome().equals(chromosome))
			return false;
		
		int apos= Math.abs(ss.getPos());
		if (apos>= getStart()&& apos<= getEnd())
			return true;
		return false;
	}
	
	public Object clone() {
		DefaultRegion clone= new DefaultRegion(getStart(), getEnd());
		clone.setID(getID());
		clone.setChromosome(getChromosome());
		clone.setSpecies(getSpecies());
		return clone;
	}
	
	public DefaultRegion(int newStart, int newEnd) {
		setStart(newStart);
		setEnd(newEnd);
	}
	
	public DefaultRegion() {
	}
	
	/**
	 * @return Returns the chromosome.
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @param chromosome The chromosome to set.
	 */
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	/**
	 * @return Returns the species.
	 */
	public Species getSpecies() {
		return species;
	}
	/**
	 * @param species The species to set.
	 */
	public void setSpecies(Species species) {
		this.species = species;
	}
	
	public String toString() {
		String str= super.toString();
		if (getID()!= null)
			str= getID()+ " "+ str;
		if (getChromosome()!= null)
			str= getChromosome()+ ":"+ str;
		if (getSpecies()!= null)
			str= getSpecies().getCommonName()+ "-"+ str;
		if (!Double.toString(getScore()).equals(Double.toString(Double.NaN)))
			str+= ": "+getScore();
			
		return str;
	}

	public double getScore() {
		return score;
	}

	public void setScore(double score) {
		this.score = score;
	}
}
