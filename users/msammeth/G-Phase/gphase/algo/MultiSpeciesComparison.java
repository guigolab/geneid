/*
 * Created on Jul 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.algo;

import java.io.File;
import java.io.PrintStream;
import java.util.Vector;

import javax.swing.JFrame;

import qalign.algo.CostTable;
import qalign.tools.FASTAWrapper;
import qalign.tools.MSFWrapper;
import qalign.tools.SequenceWrapper;
import qalign2.algo.CancelException;
import qalign2.algo.sequence.dca.QDivide;
import qalign2.algo.sequence.dca.QDivideWrapper;

import gphase.db.EnsemblDBAdaptor;
import gphase.graph.SpliceGraph;
import gphase.graph.alignment.GraphAligner;
import gphase.graph.alignment.Mapping;
import gphase.graph.gui.GraphView;
import gphase.graph.gui.PFGraph;
import gphase.io.TCoffeeWrapper;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class MultiSpeciesComparison {

	static void testGenomicSequences(Graph g) {
	
		try { 
			PrintStream p= System.out;
			Gene tst= g.getSpecies()[9].getGenes()[0];
			//Species[] spec= g.getSpecies();
			p.println(tst.getSpecies().getCommonName()+"\t"+tst.getGeneID());
			for (int i = 0; i < tst.getExons().length; i++) {
				p.println(">"+ tst.getExons()[i].getExonID());
				p.println(Graph.readSequence(tst.getExons()[i]));
			}
			p.flush(); p.close();
		} catch (Exception e) {
			e.printStackTrace();  
		}
	
	}

	public static void main(String[] args) {
					
					EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
					Graph g= adaptor.getGraphAllHomologs(SPECIES_ISMB);
					g.filterNonsense();
					g.filterNonCodingTranscripts();
//					Species[] spe= g.getSpecies();
//					for (int i = 0; i < spe.length; i++) {
//						Gene[] ge= spe[i].getGenes();
//						for (int j = 0; j < ge.length; j++) {
//							Transcript[] trans= ge[j].getTranscripts();
//							for (int k = 0; k < trans.length; k++) {
//								SpliceSite[] ss= trans[k].getSpliceChain();
//								for (int m = 0; m < ss.length; m++) {
//									ss[m].addTranscripts(new Transcript[] {trans[k]});
//								}
//							}
//						}
//					}
//					GraphHandler.writeOut(g);
		
					try { 
						PrintStream p= System.out;
						Gene tst= g.getSpecies()[0].getGenes()[0];
						Gene[] hGenes= (Gene[]) Arrays.add(tst.getHomologGenes(), tst);
						SpliceGraph[] graphs= new SpliceGraph[hGenes.length];
						for (int i = 0; i < 4; i++) {	// graphs.length
							graphs[i]= new SpliceGraph(hGenes[i].getTranscripts());
							graphs[i].init();
						}
						
							// align
						for (int i = 0; i < 1/*graphs.length*/; i++) {
							for (int j = (i+2); j < 3/*graphs.length*/; j++) {
								p.println(graphs[i]+ "\n"+ graphs[j]);
								Mapping[] maps= GraphAligner.align(graphs[i], graphs[j]);
						        //JFrame frame = GraphView.demo(new PFGraph(graphs[i]), "i");
						        JFrame frame = GraphView.demo("i", new PFGraph(graphs[i]), "j", new PFGraph(graphs[j]));
						        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
						        //JFrame frame2 = GraphView.demo(new PFGraph(graphs[j]), "j");
						        //frame2.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
								p.println(hGenes[i]+" x "+hGenes[j]);
								for (int k = 0; k < maps.length; k++) 
									p.println(maps[k]);
								p.println();
							}
						}
						
						p.flush(); p.close();
					} catch (Exception e) {
						e.printStackTrace();
					}
		}

	public static void main_ismb(String[] args) {
				
				EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
				Graph g= adaptor.getGraphAllHomologs(SPECIES_ISMB);
				g.filterNonsense();
				g.filterNonCodingTranscripts();
				GraphHandler.writeOut(g);
	
				try { 
					PrintStream p= System.out;
					Gene tst= g.getSpecies()[0].getGenes()[0];
					Gene[] hGenes= (Gene[]) Arrays.add(tst.getHomologGenes(), tst);
					//String[] ali= align(hGenes);

					MSFWrapper msf= new MSFWrapper("/home/ug/msammeth/workspace/G-Phase/gphase64881.msf");	//52161
					msf.read();
					String[] ali= msf.getSequences();
					String[] names= msf.getSeqNames();
						// reorder
					for (int i = 0; i < ali.length; i++) {
						int j= 0;
						for (j = 0; j < names.length; j++) 
							if (names[j].equalsIgnoreCase(hGenes[i].getSpecies().getCommonName()))
								break;
						if (i!= j) {
							String s= names[i];
							names[i]= names[j];
							names[j]= s;
							s= ali[i];
							ali[i]= ali[j];
							ali[j]= s;
						}
					}
					int[][] aliCoords= getAliCoordinates(ali, hGenes);
//					for (int i = 0; i < ali.length; i++) 
//						p.println(hGenes[i].getSpecies().getCommonName().substring(0, 3)+ "\t"+ ali[i]);
					
					p.flush(); p.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
	}
	
	
	private static int[][] getAliCoordinates(String[] ali, Gene[] hGenes) {
		
		int[][] ssCoords= new int[hGenes.length][];
		for (int i = 0; i < ssCoords.length; i++) 
			ssCoords[i]= hGenes[i].getExonicRegionsSSCoords();
		
		int[][] result= new int[ssCoords.length][]; 	
		for (int i = 0; i < ssCoords.length; i++) {
			java.util.Arrays.sort(ssCoords[i]);
			result[i]= new int[ssCoords[i].length];
			int lettNb= -1;
			int pos= -1;
			for (int j = 0; j < ssCoords[i].length; j++) {
				while (lettNb< ssCoords[i][j]) 
					if (Character.isLetter(ali[i].charAt(++pos)))
						++lettNb;
				result[i][j]= pos;
			}
		}
		
			// test out
		Vector[] outV= new Vector[result.length];
		for (int i = 0; i < outV.length; i++) {
			outV[i]= new Vector(result[i].length);
			outV[i].add(hGenes[i].getSpecies().getCommonName());
		}
		int maxNam= 0;	// make names equal length
		for (int i = 0; i < outV.length; i++) 
			if (((String) outV[i].elementAt(0)).length()> maxNam)
				maxNam= ((String) outV[i].elementAt(0)).length();
		maxNam+= 3; 	// spacers
		for (int i = 0; i < outV.length; i++) {
			while (((String) outV[i].elementAt(0)).length()< maxNam)
				outV[i].add(((String) outV[i].remove(0))+ ".");
		}
		
		int[] pos= new int[result.length];
		for (int i = 0; i < pos.length; i++) 
			pos[i]= 0;
		
		boolean cont= false;
		for (int i = 0; !cont&& i < pos.length; i++) 	// check for pos's left
			cont|= (pos[i]< result[i].length);
		
		while (cont) {
			int min= Integer.MAX_VALUE;	// find next min
			for (int i = 0; i < pos.length; i++) {
				if (pos[i]< result[i].length&& result[i][pos[i]]< min)
					min= result[i][pos[i]]; 
			}
			for (int i = 0; i < pos.length; i++) {	
				if (pos[i]< result[i].length&& result[i][pos[i]]== min) 
					outV[i].add(Integer.toString(result[i][pos[i]++])+ "\t");	// aligned bunch (column)
				else 
					outV[i].add("\t");						// empty
			}
			
			cont= false;	// check for pos's left
			for (int i = 0; !cont&& i < pos.length; i++) 
				cont|= (pos[i]< result[i].length);
		}
		
		for (int i = 0; i < outV.length; i++) {
			for (int j = 0; j < outV[i].size(); j++) 
				System.out.print(outV[i].elementAt(j));
			System.out.println();
		}
		return result;
	}

	static String[] align(Gene[] hGenes) {
		
		String[] names= new String[hGenes.length];
		String[] seqs= new String[hGenes.length];
		for (int i = 0; i < hGenes.length; i++) {
			names[i]= hGenes[i].getSpecies().getCommonName();
			seqs[i]= hGenes[i].getExonicRegionsSplicedSequence();
			System.out.println(names[i].substring(0,3)+ "\t"+ seqs[i]);
		}
		
//		QDivide dca= new QDivide();
//		dca.setAll(names, seqs, CostTable.DNARNA,
//				false, true, 40, false, 5);
//		dca.setApproximate(true);
//		try {
//			dca.run();
//		} catch (CancelException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		String fName= null;
		try {
			fName= File.createTempFile("gphase", ".tmp").getAbsolutePath();
			FASTAWrapper wrapper= new FASTAWrapper(fName);
			wrapper.writeFASTA(names, seqs);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		TCoffeeWrapper tco= new TCoffeeWrapper(fName);
		tco.start();
		return tco.getResult();
	}

	public static final String[] SPECIES_ISMB = { "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish", "fruitfly", "mosquito" };

}
