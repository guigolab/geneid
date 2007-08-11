package gphase;

import java.util.Vector;

import gphase.io.gtf.GTFChrReader;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Intron;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.DoubleVector;
import gphase.tools.Formatter;
import gphase.tools.IntVector;

public class ExonIntronDefinition {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		exonIntronLengthDistribution();

	}

	public static void exonIntronLengthDistribution() {
		String speName= "human";
		String annoName= "Gencode";
		String[] keywords= null;
		int eventFilter= ASMultiVariation.FILTER_STRUCTURALLY;
		
		// get basic transcripts
		String absFName= Species.getAnnotation(speName, null, annoName, keywords); 
			//Constants.getLatestUCSCAnnotation(speName, annoName, keywords);
		System.out.println("Searching for reference transcripts in "+absFName+".");
		GTFChrReader reader= new GTFChrReader(absFName);
		reader.setChromosomeWise(true);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		
		Vector<DirectedRegion> rIntrons= new Vector<DirectedRegion>(), riExon5prime= new Vector<DirectedRegion>(),
		riExon3prime= new Vector<DirectedRegion>(), sExons= new Vector<DirectedRegion>(), seIntron5prime= new Vector<DirectedRegion>(),
		seIntron3prime= new Vector<DirectedRegion>(), restExons= new Vector<DirectedRegion>(), restIntrons= new Vector<DirectedRegion>(),
		reIntron5prime= new Vector<DirectedRegion>(), reIntron3prime= new Vector<DirectedRegion>(), reiExon5prime= new Vector<DirectedRegion>(), 
		reiExon3prime= new Vector<DirectedRegion>();
		IntVector seV= new IntVector(), seIn5V= new IntVector(), seIn3V= new IntVector(), 
		reV= new IntVector(), reIn5V= new IntVector(), reIn3V= new IntVector(), 
		riV= new IntVector(), riEx5V= new IntVector(), riEx3V= new IntVector(), 
		ciV= new IntVector(), ciEx5V= new IntVector(), ciEx3V= new IntVector();
		DoubleVector seSc= new DoubleVector(), seIn5Sc= new DoubleVector(), seIn3Sc= new DoubleVector(), 
		ceSc= new DoubleVector(), ceIn5Sc= new DoubleVector(), ceIn3Sc= new DoubleVector(), 
		riSc= new DoubleVector(), riEx5Sc= new DoubleVector(), riEx3Sc= new DoubleVector(), 
		ciSc= new DoubleVector(), ciEx5Sc= new DoubleVector(), ciEx3Sc= new DoubleVector();
		Noboru nob= null;
		while (genes!= null) {
			if (nob== null)
				nob= new Noboru(speName, annoName);
			for (int i = 0; i < genes.length; i++) {
				Intron[] introns= genes[i].getIntrons();
				Exon[] exons= genes[i].getExons();
				for (int j = 0; introns!= null&& j < introns.length; j++) {
					if (!introns[j].isInternal())
						continue;
					if (genes[i].hasVariation(introns[j], ASVariation.ID_INTRON_RETENTION, eventFilter)) {
						rIntrons.add(introns[j]);
						riV.add(introns[j].getLength());
						double[] sc = nob.scoreSpliceSite(new SpliceSite[] {introns[j].getAcceptor(), introns[j].getDonor()});
						riSc.add(Math.max(sc[0], sc[1]));	//sc[0]+sc[1]);
						
						Vector<DirectedRegion> res= genes[i].getAdjacentStructures(introns[j], true, true);
						for (int k = 0; k < res.size(); k++) {
							Exon e= (Exon) res.elementAt(k);
							if (e.isInternal()) {
								riEx5V.add(e.getLength());
								riExon5prime.add(e);
								sc = nob.scoreSpliceSite(new SpliceSite[] {e.getAcceptor(), e.getDonor()});
								riEx5Sc.add(sc[0]);	//sc[0]+sc[1]);
							}
						}
						res= genes[i].getAdjacentStructures(introns[j], false, true);
						for (int k = 0; k < res.size(); k++) {
							Exon e= (Exon) res.elementAt(k);
							if (e.isInternal()) {
								riEx3V.add(e.getLength());
								riExon3prime.add(e);
								sc = nob.scoreSpliceSite(new SpliceSite[] {e.getAcceptor(), e.getDonor()});
								riEx3Sc.add(sc[1]);		//sc[0]+sc[1]);
							}
						}
					} else if (genes[i].isConstitutive(introns[j], eventFilter)) {
						restIntrons.add(introns[j]);
						ciV.add(introns[j].getLength());
						double[] sc = nob.scoreSpliceSite(new SpliceSite[] {introns[j].getAcceptor(), introns[j].getDonor()});
						ciSc.add(Math.max(sc[0],sc[1]));	//sc[0]+sc[1]);
						
						Vector<DirectedRegion> res= genes[i].getAdjacentStructures(introns[j], true, true); 
						for (int k = 0; k < res.size(); k++) {
							Exon e= (Exon) res.elementAt(k);
							if (((Exon) res.elementAt(k)).isInternal()) {
								ciEx5V.add(res.elementAt(k).getLength());
								reiExon5prime.add(res.elementAt(k));
								sc = nob.scoreSpliceSite(new SpliceSite[] {e.getAcceptor(), e.getDonor()});
								ciEx5Sc.add(sc[0]);	//sc[0]+sc[1]);
							}
						}
						
						res= genes[i].getAdjacentStructures(introns[j], false, true);
						for (int k = 0; k < res.size(); k++) { 
							Exon e= (Exon) res.elementAt(k);
							if (e.isInternal()) {
								ciEx3V.add(e.getLength());
								reiExon3prime.add(e);
								sc = nob.scoreSpliceSite(new SpliceSite[] {e.getAcceptor(), e.getDonor()});
								ciEx3Sc.add(sc[1]);	//sc[0]+sc[1]);
							}
						}
					}
				}
				for (int j = 0; j < exons.length; j++) {
					if (!exons[j].isInternal())
						continue;
					if (exons[j].hasVariation(ASVariation.ID_SKIPPED, eventFilter)) {
						sExons.add(exons[j]);
						seV.add(exons[j].getLength());
						double[] sc = nob.scoreSpliceSite(new SpliceSite[] {exons[j].getDonor(), exons[j].getAcceptor()});
						seSc.add(Math.max(sc[0],sc[1]));	//sc[0]+sc[1]);

						Vector<DirectedRegion> res= genes[i].getAdjacentStructures(exons[j], true, false); 
						for (int k = 0; k < res.size(); k++) {
							Intron intr= (Intron) res.elementAt(k);
							if (intr.isInternal()) {
								seIn5V.add(intr.getLength());
								seIntron5prime.add(intr);
								sc = nob.scoreSpliceSite(new SpliceSite[] {intr.getAcceptor(), intr.getDonor()});
								seIn5Sc.add(sc[0]);	//sc[0]+sc[1]);
							}
						}
						res= genes[i].getAdjacentStructures(exons[j], false, false);
						for (int k = 0; k < res.size(); k++) { 
							Intron intr= (Intron) res.elementAt(k);
							if (intr.isInternal()) {
								seIn3V.add(intr.getLength());
								seIntron3prime.add(intr);
								sc = nob.scoreSpliceSite(new SpliceSite[] {intr.getAcceptor(), intr.getDonor()});
								seIn3Sc.add(sc[1]);	//sc[0]+sc[1]);
							}
						}
					} else if (exons[j].isConstitutive(eventFilter)){
						restExons.add(exons[j]);
						reV.add(exons[j].getLength());
						double[] sc = nob.scoreSpliceSite(new SpliceSite[] {exons[j].getDonor(), exons[j].getAcceptor()});
						ceSc.add(Math.max(sc[0],sc[1]));	//sc[0]+sc[1]);
						
						Vector<DirectedRegion> res= genes[i].getAdjacentStructures(exons[j], true, false); 
						for (int k = 0; k < res.size(); k++) {
							Intron intr= (Intron) res.elementAt(k);
							if (((Intron) res.elementAt(k)).isInternal()) {
								reIn5V.add(res.elementAt(k).getLength());
								reIntron5prime.add(res.elementAt(k));
								sc = nob.scoreSpliceSite(new SpliceSite[] {intr.getAcceptor(), intr.getDonor()});
								ceIn5Sc.add(sc[0]);	//sc[0]+sc[1]);
							}
						}
						res= genes[i].getAdjacentStructures(exons[j], false, false);
						for (int k = 0; k < res.size(); k++) { 
							Intron intr= (Intron) res.elementAt(k);
							if (((Intron) res.elementAt(k)).isInternal()) {
								reIn3V.add(res.elementAt(k).getLength());
								reIntron3prime.add(res.elementAt(k));
								sc = nob.scoreSpliceSite(new SpliceSite[] {intr.getAcceptor(), intr.getDonor()});
								ceIn3Sc.add(sc[1]);	//sc[0]+sc[1]);
							}
						}
					}
				}
			}
			
			try {reader.read();}
			catch (Exception e) {e.printStackTrace();}
			genes= reader.getGenes();
		}
		
		System.out.println(sExons.size()+" skipped exons, "+seIntron5prime.size()+" 5' introns and "+seIntron3prime.size()+" 3' introns."+
				"\n\tmed lens: "+new Distribution(seV.toIntArray()).getMedian()+", "+new Distribution(seIn5V.toIntArray()).getMedian()+", "+ new Distribution(seIn3V.toIntArray()).getMedian()+
				"\n\tmed scores: "+new Distribution(seSc.toDoubleArray()).getMedian()+", "+new Distribution(seIn5Sc.toDoubleArray()).getMedian()+", "+ new Distribution(seIn3Sc.toDoubleArray()).getMedian());
		System.out.println(restExons.size()+" constit exons, "+reIntron5prime.size()+" 5' introns and "+reIntron3prime.size()+" 3' introns."+
				"\n\tmed lens: "+new Distribution(reV.toIntArray()).getMedian()+", "+new Distribution(reIn5V.toIntArray()).getMedian()+", "+ new Distribution(reIn3V.toIntArray()).getMedian()+
				"\n\tmed scores: "+new Distribution(ceSc.toDoubleArray()).getMedian()+", "+new Distribution(ceIn5Sc.toDoubleArray()).getMedian()+", "+ new Distribution(ceIn3Sc.toDoubleArray()).getMedian());
		System.out.println(rIntrons.size()+" retained introns, "+riExon5prime.size()+" 5' exons and "+riExon3prime.size()+" 3' exons."+
				"\n\tmed lens: "+new Distribution(riV.toIntArray()).getMedian()+", "+new Distribution(riEx5V.toIntArray()).getMedian()+", "+ new Distribution(riEx3V.toIntArray()).getMedian()+
				"\n\tmed scores: "+new Distribution(riSc.toDoubleArray()).getMedian()+", "+new Distribution(riEx5Sc.toDoubleArray()).getMedian()+", "+ new Distribution(riEx3Sc.toDoubleArray()).getMedian());
		System.out.println(restIntrons.size()+" constit introns, "+reiExon5prime.size()+" 5' exons and "+reiExon3prime.size()+" 3' exons."+
				"\n\tmed lens: "+new Distribution(ciV.toIntArray()).getMedian()+", "+new Distribution(ciEx5V.toIntArray()).getMedian()+", "+ new Distribution(ciEx3V.toIntArray()).getMedian()+
				"\n\tmed scores: "+new Distribution(ciSc.toDoubleArray()).getMedian()+", "+new Distribution(ciEx5Sc.toDoubleArray()).getMedian()+", "+ new Distribution(ciEx3Sc.toDoubleArray()).getMedian());

		double[] riScores= riSc.toDoubleArray(), riEx5Scores= riEx5Sc.toDoubleArray(), riEx3Scores= riEx3Sc.toDoubleArray();
		int[] riLengths= riV.toIntArray(), riEx5Lengths= riEx5V.toIntArray(), riEx3Lengths= riEx3V.toIntArray();
		Vector v= new Vector();
		v.add(riEx5Scores); v.add(riEx3Scores); v.add(riLengths); v.add(riEx5Lengths); v.add(riEx3Lengths);
		Arrays.synchroneousSort(riScores, v);
		for (int i = 0; i < riScores.length; i++) {
			System.out.println(Formatter.fprint(riScores[i],2)+"("+riLengths[i]+") 3' "+
					Formatter.fprint(riEx5Scores[i], 2)+"("+riEx5Lengths[i]+") 5' "+
					Formatter.fprint(riEx3Scores[i], 2)+"("+riEx3Lengths[i]+")\t"
			);
		}

	}

}
