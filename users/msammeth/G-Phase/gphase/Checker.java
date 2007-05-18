package gphase;

import gphase.model.DirectedRegion;
import gphase.model.Graph;
import gphase.model.Species;

public class Checker {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		checkSequenceExtraction();

	}

	static void checkSequenceExtraction() {
		Species spec= new Species("human");
		spec.setGenomeVersion("hg18");
		DirectedRegion reg= new DirectedRegion(
				3071, 3083, 1
				);
		reg.setChromosome("chrX");
		reg.setSpecies(spec);
		String seq= Graph.readSequence(reg);
		
		spec= new Species("fruitfly");
		spec.setGenomeVersion("bdgp43");
		reg= new DirectedRegion(
				1229794,1229797, 1
				);
		reg.setChromosome("chr4");
		reg.setSpecies(spec);
		seq= Graph.readSequence(reg);
		
		System.currentTimeMillis();
	}
}
