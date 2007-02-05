package gphase;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import javax.swing.table.DefaultTableModel;
import javax.swing.text.DefaultEditorKit.CutAction;

import gphase.algo.ASAnalyzer;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.model.TranslationLocus;
import gphase.tools.Arrays;
import gphase.tools.IntVector;

public class StopCodons {

	public static void _01_numericalAnalysis(Graph g) {
		Gene[] ge= g.getGenes();
		int tlnStarts= 0, ncTrpts= 0, cnt= 0;
		for (int i = 0; i < ge.length; i++) {
			int x= ge[i].getTranslationInitSites().length;
			if (x> 1)
				++cnt;
		}
		System.out.println(cnt+","+ge.length);
	}
	public static void _02_ASDifferences(Graph g) {
		Gene[] ge= g.getGenes();
		int tlnStarts= 0, ncTrpts= 0, cnt= 0;
		for (int i = 0; i < ge.length; i++) {
			TranslationLocus[] tls= ge[i].getTranslationLoci();
			if (tls.length< 2)
				continue;
			
			System.out.println(ge[i].getTranscripts()[0].getTranscriptID());
			ASVariation[][] vars= new ASVariation[tls.length][];
			for (int j = 0; j < tls.length; j++) 
				vars[j]= tls[j].getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
			
			
				// tls[] x tls[], common, different
			ASVariation[][][][] cvars= new ASVariation[tls.length][][][];
			for (int j = 0; j < cvars.length; j++) {
				cvars[j]= new ASVariation[tls.length][][];
				for (int k = 0; k < cvars[j].length; k++) {	// tls[] x tls[]
					cvars[j][k]= ASVariation.commonVariations(vars[j], vars[k]);
					for (int m = 0; m < cvars[j][k].length; m++) {
						System.out.print((cvars[j][k][m]==null?"0":cvars[j][k][m].length)+"/");
					}
					System.out.print("\t");
				}
				System.out.println();
			}
			System.out.println();
		}
		System.out.println(cnt+","+ge.length);
	}
	public static void _02_ASPatternDifferences(Graph g) {
		Gene[] ge= g.getGenes();
		int tlnStarts= 0, ncTrpts= 0, cnt= 0;
		for (int i = 0; i < ge.length; i++) {
			TranslationLocus[] tls= ge[i].getTranslationLoci();
			if (tls.length< 2)
				continue;
			TranslationLocus tlAll= new TranslationLocus(ge[i], ge[i].getTranscripts());
			tls= (TranslationLocus[]) Arrays.add(tls, tlAll);
			
			String id= ge[i].getTranscripts()[0].getTranscriptID();
			System.out.println(id.substring(0, id.indexOf('-')));
			System.out.print("       \t");
			for (int j = 0; j < tls.length; j++) {
				if (tls[j].getTranscripts().length== ge[i].getTranscripts().length) {
					System.out.print("all\t");
					continue;
				}
				Translation[] tln= tls[j].getTranscripts()[0].getTranslations();
				if (tln!= null&& tln[0]!= null) {
					String pos= new Integer(Math.abs(tln[0].get5PrimeEdge())).toString();
					pos= ".."+ pos.substring(Math.max(0, pos.length()- 5), pos.length());
					System.out.print(pos+ "\t");
				} else
					System.out.print("nc\t");
			}
			System.out.println();
			
			
			ASVariation[][] vars= new ASVariation[tls.length][];
			for (int j = 0; j < tls.length; j++) 
				vars[j]= tls[j].getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
			
			
				// tls[] x tls[], common, different
			ASVariation[][][][] cvars= new ASVariation[tls.length][][][];
			for (int j = 0; j < cvars.length; j++) {
				cvars[j]= new ASVariation[tls.length][][];
				
				Translation[] tln= tls[j].getTranscripts()[0].getTranslations();
				if (tls[j].getTranscripts().length== ge[i].getTranscripts().length) 
					System.out.print("all\t");
				else if (tln!= null&& tln[0]!= null) {
					String pos= new Integer(Math.abs(tln[0].get5PrimeEdge())).toString();
					pos= ".."+ pos.substring(Math.max(0, pos.length()- 5), pos.length());
					System.out.print(pos+ "\t");
				} else
					System.out.print("nc\t");
				
				for (int k = 0; k < cvars[j].length; k++) {	// tls[] x tls[]
					cvars[j][k]= ASVariation.commonVariations(vars[j], vars[k]);
					for (int m = 0; m < cvars[j][k].length; m++) {
						System.out.print((cvars[j][k][m]==null?"0":cvars[j][k][m].length));
						if (m< cvars[j][k].length- 1)
							System.out.print("/");
					}
					System.out.print("\t");
				}
				System.out.println();
			}
			System.out.println();
		}
		System.out.println(cnt+","+ge.length);
	}
	public static void _02_DifferencesSStructureCDS(Graph g) {
		g.filterNonCodingTranscripts();	// no nc
		Gene[] ge= g.getGenes();
		int tlnStarts= 0, ncTrpts= 0, cnt= 0;
		for (int i = 0; i < ge.length; i++) {
			TranslationLocus[] tls= ge[i].getTranslationLoci();
			if (tls.length< 2)
				continue;
			TranslationLocus tlAll= new TranslationLocus(ge[i], ge[i].getTranscripts());
			tls= (TranslationLocus[]) Arrays.add(tls, tlAll);
			
			String id= ge[i].getTranscripts()[0].getTranscriptID();
			System.out.println(id.substring(0, id.indexOf('-')));
			System.out.print("       \t");
			for (int j = 0; j < tls.length; j++) {
				if (tls[j].getTranscripts().length== ge[i].getTranscripts().length) {
					System.out.print("all\t");
					continue;
				}
				Translation[] tln= tls[j].getTranscripts()[0].getTranslations();
				if (tln!= null&& tln[0]!= null) {
					String pos= new Integer(Math.abs(tln[0].get5PrimeEdge())).toString();
					pos= ".."+ pos.substring(Math.max(0, pos.length()- 5), pos.length());
					System.out.print(pos+ "\t");
				} else
					System.out.print("nc\t");
			}
			System.out.println();
			
			
			ASVariation[][] vars= new ASVariation[tls.length][];
			for (int j = 0; j < tls.length; j++) 
				vars[j]= tls[j].getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
			
			
				// tls[] x tls[], common, different
			ASVariation[][][][] cvars= new ASVariation[tls.length][][][];
			Vector remV= new Vector();
			for (int j = 0; j < cvars.length; j++) {
				cvars[j]= new ASVariation[tls.length][][];
				
				Translation[] tln= tls[j].getTranscripts()[0].getTranslations();
				if (tls[j].getTranscripts().length== ge[i].getTranscripts().length) 
					System.out.print("all\t");
				else if (tln!= null&& tln[0]!= null) {
					String pos= new Integer(Math.abs(tln[0].get5PrimeEdge())).toString();
					pos= ".."+ pos.substring(Math.max(0, pos.length()- 5), pos.length());
					System.out.print(pos+ "\t");
				} else
					System.out.print("nc\t");
				
				for (int k = 0; k < cvars[j].length; k++) {	// tls[] x tls[]
					cvars[j][k]= ASVariation.commonVariations(vars[j], vars[k]);
					for (int m = 0; m < cvars[j][k].length; m++) {
							// get variations common to something
						if (m== 0&& j!= k&&
								tls[j].getTranscripts().length== ge[i].getTranscripts().length) {
							for (int x = 0; x < cvars[j][k][m].length; x++) 
								remV.add(cvars[j][k][m][x]);
						}
						System.out.print((cvars[j][k][m]==null?"0":cvars[j][k][m].length));
						if (m< cvars[j][k].length- 1)
							System.out.print("/");
					}
					System.out.print("\t");
				}
				System.out.println();
			}
			
			Vector keepV= new Vector();
			for (int j = 0; j < tls.length; j++) {
				if (tls[j].getTranscripts().length!= ge[i].getTranscripts().length)
					continue;
				for (int k = 0; vars[j]!= null&& k < vars[j].length; k++) {
					int m;
					for (m = 0; m < remV.size(); m++) 
						if (vars[j][k]== remV.elementAt(m))
							break;
					if (m== remV.size())
						keepV.add(vars[j][k]);
				}
			}
			System.out.println("==> "+keepV.size()+" AS events due to differences in translational clusters.");
			for (int j = 0; j < keepV.size(); j++) {
				ASVariation var= (ASVariation) keepV.elementAt(j);
				System.out.println(var.toStringUCSC()+"\t"+var);
			}
			System.out.println();
		}
		System.out.println(cnt+","+ge.length);
	}
	public static void _03_utrAS(Graph g) {
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as= ge[i].getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
			int cnt= 0;
			for (int j = 0; as!= null&& j < as.length; j++) 
				if (as[j].is_contained_in_5UTR())
					++cnt;
			if (cnt> 0)
				System.out.println(cnt+"\t"+ge[i].getTranscripts()[0]);
		}
	}
	
	public static void _03_utrLongExons(Graph g) {
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] trpts= ge[i].getTranscripts();
			int longestExon= -1;
			for (int j = 0; j < trpts.length; j++) {
				if (trpts[j].isNonCoding())
					continue;
				Exon[] exons= trpts[j].getExons();
				for (int k = 0; k < exons.length; k++) {	// not 1st exon
					if (exons[k].isCoding())
						break;
					if (k> 0&& exons[k].getLength()> longestExon)
						longestExon= exons[k].getLength();
				}
			}
			if (longestExon> 150)
				System.out.println(longestExon+"\t"+ge[i].getTranscripts()[0]);
		}
	}

	static void _03_ncExonsORF_outputStops(Vector v, PrintStream p) {
		for (int i = 0; i < v.size(); i++) {
			DirectedRegion dir= (DirectedRegion) v.elementAt(i);
			String seq= Graph.readSequence(dir);
			int[] stopCnt= Translation.getStopCount(seq);
			java.util.Arrays.sort(stopCnt);
			for (int j = 0; j < stopCnt.length; j++) 
				p.print(stopCnt[j]+ "\t");
			p.print(dir.getLength()+"\t("+(dir.getLength()/64f)+")");	//3/64 /3 for codonns..
			if (dir instanceof Exon&& stopCnt[0]> 0)
				p.print("\t"+ dir.toUCSCString());
			p.println();
		}
	}
	
	public static void _03_ncGenes(Graph g) {
		Gene[] ge= g.getGenes();
		for (int i = 0; i < ge.length; i++) {
			if (!ge[i].isProteinCoding())
				System.out.println(ge[i].getTranscripts()[0]);
		}
	}

	public static void _03_ncExonsORF(Graph g) {
		Gene[] ge= g.getGenes();
		Vector ncxV= new Vector();
		Vector nciV= new Vector();
		Vector utrxV= new Vector();
		Vector utriV= new Vector();
		Vector cdsxV= new Vector();
		Vector cdsiV= new Vector();
		Comparator compi= new DirectedRegion.PositionComparator(); 
		for (int i = 0; i < ge.length; i++) { 
			Transcript[] trpts= ge[i].getTranscripts();
			for (int j = 0; j < trpts.length; j++) {
				Exon[] exons= trpts[j].getExons();
				for (int k = 1; k < exons.length- 1; k++) {	// never 1st or last exon
					if (!ge[i].isProteinCoding()) {
						ncxV= Arrays.addUnique(ncxV, new Object[] {exons[k]}, compi);
						DirectedRegion dir= new DirectedRegion(exons[k-1].get3PrimeEdge(), exons[k].get5PrimeEdge(), ge[i].getStrand());
						dir.setSpecies(ge[i].getSpecies());
						dir.setChromosome(ge[i].getChromosome());
						nciV= Arrays.addUnique(nciV, new Object[] {dir}, compi);
						if (k== exons.length- 2) { 
							dir= new DirectedRegion(exons[k].get3PrimeEdge(), exons[k+1].get5PrimeEdge(), ge[i].getStrand());
							dir.setSpecies(ge[i].getSpecies());
							dir.setChromosome(ge[i].getChromosome());
							if (ge[i].getChromosome()== null)
								System.out.println();
							nciV= Arrays.addUnique(nciV, new Object[] {dir}, compi);
						}
					} else {	// protein coding gene
						if (exons[k].get3PrimeEdge()< ge[i].getMinCDSStart()) {
							utrxV= Arrays.addUnique(utrxV, new Object[] {exons[k]}, compi);
							DirectedRegion dir= new DirectedRegion(exons[k-1].get3PrimeEdge(), exons[k].get5PrimeEdge(), ge[i].getStrand());
							dir.setSpecies(ge[i].getSpecies());
							dir.setChromosome(ge[i].getChromosome());
							utriV= Arrays.addUnique(utriV, new Object[] {dir}, compi);
						}
						if (exons[k].get5PrimeEdge()> ge[i].getMaxCDSEnd()) {
							utrxV= Arrays.addUnique(utrxV, new Object[] {exons[k]}, compi);
							DirectedRegion dir= new DirectedRegion(exons[k].get3PrimeEdge(), exons[k+1].get5PrimeEdge(), ge[i].getStrand());
							dir.setSpecies(ge[i].getSpecies());
							dir.setChromosome(ge[i].getChromosome());
							utriV= Arrays.addUnique(utriV, new Object[] {dir}, compi);
						}
						if (trpts[j].isCoding()&& exons[k].get5PrimeEdge()> trpts[j].getTranslations()[0].get5PrimeEdge()
								&& exons[k].get3PrimeEdge()< trpts[j].getTranslations()[0].get3PrimeEdge()) {
							cdsxV= Arrays.addUnique(cdsxV, new Object[] {exons[k]}, compi);
							DirectedRegion dir= new DirectedRegion(exons[k].get3PrimeEdge(), exons[k+1].get5PrimeEdge(), ge[i].getStrand());
							dir.setSpecies(ge[i].getSpecies());
							dir.setChromosome(ge[i].getChromosome());
							cdsiV= Arrays.addUnique(cdsiV, new Object[] {dir}, compi);
							dir= new DirectedRegion(exons[k-1].get3PrimeEdge(), exons[k].get5PrimeEdge(), ge[i].getStrand());
							dir.setSpecies(ge[i].getSpecies());
							dir.setChromosome(ge[i].getChromosome());
							cdsiV= Arrays.addUnique(cdsiV, new Object[] {dir}, compi);
						}

							
					}
						
				}
			}
		}
		
		PrintStream p= null;
		try {
			p= new PrintStream("exon_intron_3frames.txt");
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		p.println("nc exons:");
		_03_ncExonsORF_outputStops(ncxV, p);
		
		p.println("nc introns:");
		_03_ncExonsORF_outputStops(nciV, p);
		
		p.println("utr exons:");
		_03_ncExonsORF_outputStops(utrxV, p);
		
		p.println("utr introns:");
		_03_ncExonsORF_outputStops(utriV, p);

		p.println("cds exons:");
		_03_ncExonsORF_outputStops(cdsxV, p);
		
		p.println("cds introns:");
		_03_ncExonsORF_outputStops(cdsiV, p);

	}
	
	public static void _04_atg_density(Graph g) {
		Gene[] ge= g.getGenes();
		int c5UTRi= 0, c5UTRe= 0, cCDSi= 0, cCDSe= 0;
		int c5UTRiLen= 0, c5UTReLen= 0, cCDSiLen= 0, cCDSeLen= 0;
		for (int i = 0; i < ge.length; i++) {
			DirectedRegion[] utr5i= ge[i].getIntrons(Gene.REGION_REAL_5UTR);
			utr5i= DirectedRegion.getUniqueRegions(utr5i);
			for (int j = 0; utr5i!= null&& j < utr5i.length; j++) {
				int[] all= Translation.getStartCount(Graph.readSequence(utr5i[j]));
				for (int k = 0; k < all.length; k++) 
					c5UTRi+= all[k];
				c5UTRiLen+= utr5i[j].getLength();
			}
			
			DirectedRegion[] utr5e= ge[i].getExons(Gene.REGION_REAL_5UTR);
			utr5e= DirectedRegion.getUniqueRegions(utr5e);
			for (int j = 0; utr5e!= null&& j < utr5e.length; j++) {
				int[] all= Translation.getStartCount(Graph.readSequence(utr5e[j]));
				for (int k = 0; k < all.length; k++) 
					c5UTRe+= all[k];
				c5UTReLen+= utr5e[j].getLength();
			}
			
			DirectedRegion[] cdsi= ge[i].getIntrons(Gene.REGION_REAL_CDS);
			cdsi= DirectedRegion.getUniqueRegions(cdsi);
			for (int j = 0; cdsi!= null&& j < cdsi.length; j++) {
				int[] all= Translation.getStartCount(Graph.readSequence(cdsi[j]));
				for (int k = 0; k < all.length; k++) 
					cCDSi+= all[k];
				cCDSiLen+= cdsi[j].getLength();
			}
			
			DirectedRegion[] cdse= ge[i].getExons(Gene.REGION_REAL_CDS);
			cdse= DirectedRegion.getUniqueRegions(cdse);
			for (int j = 0; cdse!= null&& j < cdse.length; j++) {
				int[] all= Translation.getStartCount(Graph.readSequence(cdse[j]));
				for (int k = 0; k < all.length; k++) 
					cCDSe+= all[k];
				cCDSeLen+= cdse[j].getLength();
			}
		}
		System.out.println("5UTR intron:\t"+c5UTRi+"\t"+((float)c5UTRi/ c5UTRiLen));
		System.out.println("5UTR exon:\t"+c5UTRe+"\t"+((float)c5UTRe/ c5UTReLen));
		System.out.println("CDS intron:\t"+cCDSi+"\t"+((float)cCDSi/ cCDSiLen));
		System.out.println("CDS exon:\t"+cCDSe+"\t"+((float)cCDSe/ cCDSeLen));
	}
	
	public static void _04_atg_distribution(Graph g) {
		Gene[] ge= g.getGenes();
		HashMap c5UTRiHash= new HashMap(), c5UTReHash= new HashMap(), cCDSiHash= new HashMap(), cCDSeHash= new HashMap();
		int c5UTRiLen= 0, c5UTReLen= 0, cCDSiLen= 0, cCDSeLen= 0;
		for (int i = 0; i < ge.length; i++) {
			IntVector vec= new IntVector();
			for (int j = 0; j < ge[i].getTranscripts().length; j++) {
				if (!ge[i].getTranscripts()[j].isCoding())
					continue;
				String seq= ge[i].getTranscripts()[j].getSplicedSequence();
				seq= seq.substring(0,	// only 5'UTR 
						ge[i].getTranscripts()[j].getExonicPosition(ge[i].getTranscripts()[j].getTranslations()[0].get5PrimeEdge()));
				int x;
				int[] v= vec.toIntArray();
				for (x = 0; x < v.length; x++) 
					if (seq.length()== v[x])
						break;
				if (x< v.length)	// only unique UTRs
					continue;
				int[] all= Translation.getCodonPositions(
						new String[] {Translation.START_CODON}, seq);
				for (int k = 0; k < all.length; k++) {
					Integer key= new Integer(seq.length()- all[k]); 
					Object val= c5UTReHash.get(key);
					if (val== null)
						val= new Integer(1);
					else
						val= new Integer(((Integer) val).intValue()+ 1);
					c5UTReHash.put(key, val);
				}
			}
		}
		
		Object[] keys= c5UTReHash.keySet().toArray();
		java.util.Arrays.sort(keys);
		System.out.println("5UTR exon");
		int cumul= 0;
		for (int i = 0; i < keys.length; i++) {
			cumul+= ((Integer) c5UTReHash.get(keys[i])).intValue();
			System.out.println(keys[i]+"\t"+c5UTReHash.get(keys[i]) +"\t"+ cumul);
		}
	}
	
	public static void main(String[] args) {
		Graph g= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);
		//_01_numericalAnalysis(g);
		_02_DifferencesSStructureCDS(g);
		//_03_utrLongExons(g);
		//_03_ncExonsORF(g);
		//_03_ncGenes(g);
		//_04_atg_distribution(g);
	}
}
