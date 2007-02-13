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
import gphase.io.gtf.ArabidopsisGFFReader;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;

public class TransSpeciesComparison {

	public static Graph getGraph(String commonName) {
		for (int i = 0; i < Species.SP_NAMES_COMMON.length; i++) {
			if (!commonName.equalsIgnoreCase(Species.SP_NAMES_COMMON[i]))		
				continue;
			return GraphHandler.readIn(GraphHandler.getGraphAbsPath(new Species(Species.SP_NAMES_COMMON[i]))+"_download");
		}
		System.err.println("Species "+commonName+" not known.");
		return null;
	}
	
	public static void _00_mainLoop() {
		for (int i = 0; i < Species.SP_NAMES_COMMON.length; i++) {
			if ((!Species.SP_NAMES_COMMON[i].equals("seasquirt"))&&
					(!Species.SP_NAMES_COMMON[i].equals("seqsquirt2")))		// !!! REMOVE
				continue;
			System.out.println(Species.SP_NAMES_COMMON[i]);
			Graph g= null;
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			try {
				String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
				Species dummySpec= new Species(Species.SP_NAMES_COMMON[i]);
				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
//				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_download");
//				g.getSpecies()[0].filter();
//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
//				g.filterNMDTranscripts();
//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
				String sfx= ".landscape.GTAGonly.NMD";
				ASAnalyzer.test04_determineVariations(g, iname, sfx, false);
				
//				String sfx= ".landscape";
//				ASAnalyzer.test04_determineVariations(g, iname, sfx, false);
//				sfx= ".landscape.GTAGonly";
//				ASAnalyzer.test04_determineVariations(g, iname, sfx, true);
//				PrintStream p= null;
				//p= new PrintStream("graph"+File.separator+Species.SP_NAMES_COMMON[i]+"_length_distr.txt");
				//ASAnalyzer.test03c_statisticAll(g, p);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public static void main(String[] args) {
		
//		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);
//		g.filterNonCodingTranscripts();
//		g.filterNMDTranscripts();
//		
//		g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);
//		g.filterNMDTranscripts();
//		
//		g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_REFSEQ_CODING_FROM_UCSC);
//		g.filterNMDTranscripts();
//		
//		g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENSEMBL_CODING_FROM_UCSC);
//		g.filterNMDTranscripts();
		
//		Graph g= getGraph("seasquirt");
//		g.getSpecies()[0].filter();
//		g.filterNMDTranscripts();
		
		_01_testCiona();
		//_00_mainLoop();
		//_01_testDroso();
		//_01_testAthaliana();
	}

	public static void _01_testCiona() {
		Graph g= null;
		try {
			//"C_intestinalis_RefSeqGenes_from_UCSC.gtf"
			//"C_intestinalis_RefSeqGenes_from_UCSC.gtf"
			// "C_intestinalis_splicedESTs_from_UCSC.gtf"
			String fName= "/home/msammeth/annotations"+File.separator+
			"C_intestinalis_splicedESTs_from_UCSC_200503.gtf";
			System.out.println("Reading annotation: "+ fName);
			g= ASAnalyzer.getGraph(fName);
			//g.getSpecies()[0].setBuildVersion(200503);
			//g.getSpecies()[0].spNumber= 17;	
			//g.filter();	// reading errors, EOF
			g.filterUnInformative();	// ESTs
//			try {
//				GTFWrapper wrapper= new GTFWrapper("Ciona_filtered.gtf");
//				wrapper.addGTFObject(g.getGenes());
//				wrapper.writeGTF();
//			} catch (Exception e) {
//				e.printStackTrace();
//			}
			ASAnalyzer.test04_determineVariations(g, System.out);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		String path= "/home/msammeth/annotations";
		String[] list= new File(path).list();
		for (int i = 0; i < list.length; i++) {
			try {
				String fName= path+File.separator+list[i];
				System.out.println("Reading annotation: "+fName);
				g= ASAnalyzer.getGraph(fName);
				g.filterUnInformative();	// ESTs
				ASAnalyzer.test04_determineVariations(g, System.out);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
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
	
	public static void _01_testAthaliana() {
		ArabidopsisGFFReader reader= new ArabidopsisGFFReader("A_thaliana.gff");
		Graph g= reader.getGraph();
		ASAnalyzer.test04_determineVariations(g, System.out);
	}
}
