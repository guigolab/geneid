/*
 * Created on Feb 4, 2007
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase;

import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Method;

import javax.swing.JFrame;

import gphase.algo.ASAnalyzer;
import gphase.db.EnsemblDBAdaptor;
import gphase.gui.SpliceOSigner;
import gphase.io.gtf.ArabidopsisGFFReader;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFWrapper;
import gphase.io.gtf.JGI_GFFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.SpliceSite;

public class TransSpeciesComparisonExil {

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
				if ((!Species.SP_NAMES_COMMON[i].equals("seasquirt")))		// !!! REMOVE
					continue;
				System.out.println(Species.SP_NAMES_COMMON[i]);
				Graph g= null;
				EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
				try {
					String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
					Species dummySpec= new Species(Species.SP_NAMES_COMMON[i]);
					g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
					//g.filterScaffold();	// REMOVE !!!
	//				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_download");
	//				g.getSpecies()[0].filter();
	//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
	//				g.filterNMDTranscripts();
	//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
					String sfx= ".landscape.GTAGonly";
					ASAnalyzer.test04_determineVariations(g, iname, sfx, true);
					
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



	/**
		 * target method must comply with signature (Graph, Printstream)
		 * @param target
		 */
		public static void _00_mainLoop(Method[] target) {
			String[] spe= Species.SP_NAMES_METAZOA;
			for (int i = 0; i < spe.length; i++) {
				if (!spe[i].equals("worm"))
					continue;
				System.out.println(Species.SP_NAMES_METAZOA[i]);
				Graph g= null;
				try {
					String iname= "annotation"+ File.separator+ "ensembl"+ File.separator+
						Species.SP_NAMES_METAZOA[i]+"_42_filtDNA.gtf";
					EncodeWrapper wrapper= new EncodeWrapper(iname);
					g= wrapper.getGraph(false);
	//				String iname= "graph"+ File.separator+ Species.SP_NAMES_METAZOA[i];
	//				Species dummySpec= new Species(Species.SP_NAMES_METAZOA[i]);
	//				File f= new File(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
	//				if (f.exists())
	//					g= GraphHandler.readIn(f.getAbsolutePath());
	//				else {
	//					g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_download");
	//					g.getSpecies()[0].filter(); 
	//					GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
	//				}
					
					for (int j = 0; j < target.length; j++) {
						try {
							PrintStream p= new PrintStream(iname+"."+target[j].getName());
							target[j].invoke(null, new Object[] {g, p});
							p.flush();
							p.close();
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		}



	/**
	 * target method must comply with signature (Graph, Printstream)
	 * @param target
	 */
	public static void _00_mainLoopEnsemblChrom(Method[] target) {
		String[] spe= Species.SP_NAMES_METAZOA;
		for (int i = 0; i < spe.length; i++) {
			if (!spe[i].equals("worm"))
				continue;
			System.out.println(Species.SP_NAMES_METAZOA[i]);
			Graph g= null;
			try {
				String iname= "annotation"+ File.separator+ "ensembl"+ File.separator+
					Species.SP_NAMES_METAZOA[i]+"_42_filtDNA.gtf";
				EncodeWrapper wrapper= new EncodeWrapper(iname);
				g= wrapper.getGraph(false);
//				String iname= "graph"+ File.separator+ Species.SP_NAMES_METAZOA[i];
//				Species dummySpec= new Species(Species.SP_NAMES_METAZOA[i]);
//				File f= new File(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
//				if (f.exists())
//					g= GraphHandler.readIn(f.getAbsolutePath());
//				else {
//					g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_download");
//					g.getSpecies()[0].filter(); 
//					GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
//				}
				
				for (int j = 0; j < target.length; j++) {
					try {
						PrintStream p= new PrintStream(iname+"."+target[j].getName());
						target[j].invoke(null, new Object[] {g, p});
						p.flush();
						p.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}



	public static void _00_RefSeq_Gencode_GTAGonly() {
			String id= ASAnalyzer.INPUT_ENCODE;	// INPUT_REFSEQ_CODING_FROM_UCSC;
			Graph g= ASAnalyzer.getGraph(id);
			g.filterCodingClustersOnly();
			try {
				String iname= "graph"+ File.separator+ id+"-codingCluster";
				//g.filterScaffold();	// REMOVE !!!
//				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_download");
//				g.getSpecies()[0].filter();
//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA");
//				g.filterNMDTranscripts();
//				GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
				String sfx= ".landscape.GTAGonly";
				ASAnalyzer.test04_determineVariations(g, iname, sfx, true);
				
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
	
	public static void _01_analyzeCiona() {
		for (int i = 0; i < Species.SP_NAMES_COMMON.length; i++) {
			if ((!Species.SP_NAMES_COMMON[i].equals("seasquirt")))		// !!! REMOVE
				continue;
			System.out.println(Species.SP_NAMES_COMMON[i]);
			Graph g= null;
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			try {
				String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
				Species dummySpec= new Species(Species.SP_NAMES_COMMON[i]);
				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
				ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);			
				
				
				int eq= 0, neq= 0, slip= 0, slp2= 0;
				for (int j = 0; j < vars.length; j++) {
					if (!vars[j][0].toString().equals("1^3- , 2^4-"))
						continue;
					System.out.println(vars[j][0].toString());
					for (int k = 0; k < vars[j].length; k++) {
						SpliceSite[] sc1= vars[j][k].getSpliceChain1();
						SpliceSite[] sc2= vars[j][k].getSpliceChain2();
						int dif1= Math.abs(sc1[0].getPos()- sc2[0].getPos());
						int dif2= Math.abs(sc1[1].getPos()- sc2[1].getPos());
						int itr1= Math.abs(sc1[1].getPos()- sc1[0].getPos());
						int itr2= Math.abs(sc2[1].getPos()- sc2[0].getPos());
						if (dif1== dif2) {
							++eq;
							if (dif1< 10)
								++slp2;
						} else {							
							++neq;
							if (dif1< 10|| dif2< 10|| vars[j][k].getGene().getChromosome().contains("scaffold"))
								++slip;
							else {
								if (vars[j][k].isBorderEvent())
									System.out.print("*** ");
								System.out.println(dif1+","+dif2+"\t"+itr1+","+itr2);
								System.out.println(vars[j][k].toStringUCSC()+"\t"+vars[j][k].getTranscript1().getTranscriptID()+
										"\t"+ vars[j][k].getTranscript2().getTranscriptID()+"\t"+vars[j][k].getTranscript1().getGene().getGeneID());
							}
						}
					}
					System.out.println("equal= "+eq+"("+slp2+"),\tnot eq= "+neq+" with dubious= "+slip);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static void _01_analyzeCiona2() {
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
				ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
				int eq= 0, neq= 0, slip= 0, slp2= 0;
				
				for (int j = 0; j < vars.length; j++) {
					if (!vars[j][0].toString().equals("1^ , 2^"))
						continue;
					System.out.println(vars[j][0].toString());
					for (int k = 0; k < vars[j].length; k++) {
						SpliceSite[] sc1= vars[j][k].getSpliceChain1();
						SpliceSite[] sc2= vars[j][k].getSpliceChain2();
						int dif1= Math.abs(sc1[0].getPos()- sc2[0].getPos());
						if (dif1< 10) 
							++slp2;
					}
					System.out.println("little= "+slp2+"\tof "+vars[j].length);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	public static void _01_ciona() {
		for (int i = 0; i < Species.SP_NAMES_COMMON.length; i++) {
			if ((!Species.SP_NAMES_COMMON[i].equals("seasquirt")))		// !!! REMOVE
				continue;
			System.out.println(Species.SP_NAMES_COMMON[i]);
			Graph g= null;
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			try {
				String iname= "graph"+ File.separator+ Species.SP_NAMES_COMMON[i];
				Species dummySpec= new Species(Species.SP_NAMES_COMMON[i]);
				g= GraphHandler.readIn(GraphHandler.getGraphAbsPath(dummySpec)+"_filtDNA_NMD");
				ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
				for (int j = 0; j < vars.length; j++) {
					if (!vars[j][0].toString().equals("1-2^ , 0"))
						continue;
					for (int k = 0; k < vars[j].length; k++) 
						System.out.println(vars[j][k].getGene().getGeneID());
				}
				String[] ids= new String[] {
						"ENSCING00000003762",
						"ENSCING00000009296",
						"ENSCING00000000150",
						"ENSCING00000008202",
						"ENSCING00000015904",
						"ENSCING00000008665",
						"ENSCING00000004415",
						"ENSCING00000004739",
						"ENSCING00000005479",
						
						"ENSCING00000006076",
						"ENSCING00000002676",
						"ENSCING00000011628",
						"ENSCING00000008014",
						"ENSCING00000014604",
						"ENSCING00000001957",
						"ENSCING00000001957",
						"ENSCING00000005356",
						"ENSCING00000005356",
						"ENSCING00000003848",
						"ENSCING00000004147",
						"ENSCING00000004147",
						
						"ENSCING00000002430",
						"ENSCING00000001017",
						"ENSCING00000005722",
						"ENSCING00000002514",
						"ENSCING00000005377",
						"ENSCING00000004304",
						"ENSCING00000009011",
						"ENSCING00000003244"
						
				};
				for (int j = 0; j < ids.length; j++) {
					Gene ge= g.getSpecies()[0].getGene(ids[j]);
					SpliceOSigner designer= new SpliceOSigner(ge);
					JFrame myFrame= new JFrame();
					myFrame.getContentPane().add(designer);
					myFrame.pack();
					myFrame.setVisible(true);
				}
				
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
		
		Method[] m= new Method[2];
		try {
			m[0]= ASAnalyzer.class.getMethod("test03c_statisticExonsWithRegulatorSequences", new Class[] {Graph.class, PrintStream.class});
			m[1]= ASAnalyzer.class.getMethod("test03c_statisticIntronsWithGCcontent", new Class[] {Graph.class, PrintStream.class});
		} catch (Exception e) {
			e.printStackTrace();
		}
		_00_mainLoop(m);
		
		
		//_00_RefSeq_Gencode_GTAGonly();
		//_01_analyzeCiona2();
		//_01_analyzeCiona();
		//_01_testCiona();
		//_00_mainLoop();
		//_01_testDroso();
		//_01_testAthaliana();
	}

	public static void _01_testCiona() {
		Graph g= null;
		try {
			// "/home/msammeth/annotations"+File.separator+
			//"C_intestinalis_RefSeqGenes_from_UCSC.gtf"
			//"C_intestinalis_RefSeqGenes_from_UCSC.gtf"
			// "C_intestinalis_splicedESTs_from_UCSC.gtf"
			// 			"C_intestinalis_splicedESTs_from_UCSC_200503.gtf";
			String fName= "encode"+ File.separator+ "ciona_genes_small.gff";
			System.out.println("Reading annotation: "+ fName);
			//g= ASAnalyzer.getGraph(fName);
			JGI_GFFWrapper wrapper= new JGI_GFFWrapper(new File(fName).getAbsolutePath());
			g= wrapper.getGraph();
			//g.getSpecies()[0].setBuildVersion(200503);
			//g.getSpecies()[0].spNumber= 17;	
			g.filter();	// reading errors, EOF
			g.filterNonGTAGIntrons();
			//g.filterUnInformative();	// ESTs
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
					if (vars[i][j].isAffecting5UTR()&& !(vars[i][j].isContained5UTR())&&
							!(vars[i][j].isAffectingCDS()))
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
