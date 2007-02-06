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
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;

public class TransSpeciesComparison {

	public static void _00_mainLoop() {
		for (int i = 0; i < Species.SP_NAMES_COMMON.length; i++) {
			System.out.println(Species.SP_NAMES_COMMON[i]);
			Graph g= null;
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			try {
				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(new Species(Species.SP_NAMES_COMMON[i]))+"_filtNonsense");
				String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
				String sfx= ".landscape";
				PrintStream p= null;
				p= new PrintStream("graph"+File.separator+Species.SP_NAMES_COMMON[i]+"_length_distr.txt");
				//ASAnalyzer.test04_determineVariations(g, iname, sfx);
				ASAnalyzer.test03c_statisticAll(g, p);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public static void main(String[] args) {
		//_00_mainLoop();
		_01_testDroso(); 
	}

	public static void _01_testDroso() {
		Graph g= null;
		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
		try {
			g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(new Species("fruitfly"))+"_download");
			ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);

			for (int i = 0; i < vars.length; i++) {
				if (!vars[i][0].isIntronRetention())
					continue;
				for (int j = 0; j < vars[i].length; j++) {
					if (vars[i][j].is_affecting_5UTR()&& !(vars[i][j].is_contained_in_5UTR())&&
							!(vars[i][j].is_affecting_CDS()))
						System.out.println(vars[i][j].getSsRegionIDs()+"\t"+vars[i][j].toStringUCSC());
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
