/*
 * Created on Nov 15, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.algo;

import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.Species;

import java.util.Vector;

/**
 * k-partite graph to capture multiple gene homologs 
 * 
 * @author msammeth
 */
public class GeneSet {

	Vector[] set;	// same species order like in graph
	Graph graph;
	
	/**
	 * Inits a set with the given number of species.
	 */
	public GeneSet(Graph g) {
		graph= g;
		set= new Vector[g.getSpecies().length];
		for (int i = 0; i < set.length; i++) 
			set[i]= new Vector();
	}


	/**
	 * @param g
	 * @return all homologs of <code>g</code> within the same
	 * species
	 */
	public Gene[] getGeneParalogs(Gene g) {

		int specNb= getSpeciesNumber(g.getSpecies());
		Vector resV= (Vector) set[specNb].clone();
		resV.remove(g);
		
		return Gene.toGeneArray(resV);
	}
	
	/**
	 * 
	 * @param spec
	 * @return
	 */ 
	int getSpeciesNumber(Species spec) {
		int specNb;
		for (specNb = 0; specNb < graph.getSpecies().length; specNb++) {
			if (graph.getSpecies()[specNb]== spec)	
				break;
		}
		return specNb;
	}

	
	/**
	 * Adds new gene and recursionally adds new homologs.
	 * @param g
	 * @return
	 */
	public boolean add(Gene g) {

		int specNb= getSpeciesNumber(g.getSpecies());
		
		for (int i = 0; i < set[specNb].size(); i++) 
			if (set[specNb].elementAt(i)== g)
				return false;
		
		set[specNb].add(g);	// add
		for (int i = 0; i < graph.getSpecies().length; i++) {	
			// also in same species!
			GeneHomology[] homols= g.getHomologies(graph.getSpecies()[i]);
			if (homols== null) {
				if (graph.getSpecies()[i]!= g.getSpecies())
					System.err.println("Gene "+g+" has no homolog in "+graph.getSpecies()[i]);
				continue;
			}
			for (int j = 0; j < homols.length; j++) {
				add(homols[j].getOtherGene(g));
			}
		}
		
		return true;
	}
	
	public Gene[][] getSet() {
		if (set== null)
			return null;
		
		Gene[][] result= new Gene[set.length][];
		for (int i = 0; i < result.length; i++) {
			result[i]= new Gene[set[i].size()];
			for (int j = 0; j < result.length; j++) 
				result[i][j]= (Gene) set[i].elementAt(j);
		}
		
		return result;
	}
}
