package gphase.model;

public class DefaultDirectedRegion extends DirectedRegion {
	String chromosome = null;
	Species species = null;
	
	public DefaultDirectedRegion(int newStart, int newEnd, int strand) {
		setStrand(strand);
		setStart(newStart);
		setEnd(newEnd);
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
		String borders= super.toString();
		return getSpecies().getCommonName()+ "-"+ getChromosome()+ ":"+ borders;
	}
}
