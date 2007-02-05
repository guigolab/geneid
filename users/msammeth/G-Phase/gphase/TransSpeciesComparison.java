/*
 * Created on Feb 4, 2007
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase;

import java.io.File;
import java.io.PrintStream;

import gphase.algo.ASAnalyzer;
import gphase.db.EnsemblDBAdaptor;
import gphase.model.Graph;
import gphase.model.Species;

public class TransSpeciesComparison {

	public static void main(String[] args) {
		for (int i = 14; i < Species.SP_NAMES_COMMON.length; i++) {
			System.out.println(Species.SP_NAMES_COMMON[i]);
			Graph g= null;
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			try {
				g= adaptor.getGraphAllGenes(Species.SP_NAMES_COMMON[i]);
				String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
				String sfx= ".landscape";
				ASAnalyzer.test04_determineVariations(g, iname, sfx);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
}
