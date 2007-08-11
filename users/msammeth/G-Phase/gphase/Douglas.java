package gphase;

import gphase.Analyzer.VectorSizeComparator;
import gphase.algo.ASAnalyzer;
import gphase.io.DomainWrapper;
import gphase.io.DouglasDomainWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASEventold;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.DoubleVector;
import gphase.tools.Formatter;
import gphase.tools.IntVector;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Vector;

import org.apache.commons.collections.BidiMap;
import org.apache.commons.collections.bidimap.DualHashBidiMap;


public class Douglas {
	
	
	public static void main(String[] args) {
		//_05_01_entireDistribution();
		
		//_05_00_pfamDomains_all();
		//_05_00_pfamDomains_all_contained();
		
		//_05_exonsStatistic();
		
		_070716_getAllAAsequences();
	}
	
	static void _070716_getAllAAsequences() {
		try {
			GTFChrReader reader= new GTFChrReader("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_CDS_chr8.gtf");
			reader.setReadGene(true);
			reader.setChromosomeWise(true);
			reader.read();
			Gene[] g= reader.getGenes();
			BufferedWriter writer= new BufferedWriter(new FileWriter("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_CDS_aaSequences.txt"));
			Species spe= new Species("human");
			spe.setGenomeVersion("hg18");
			while (g!= null) {
				
				for (int i = 0; i < g.length; i++) {
					g[i].setSpecies(spe);
					Transcript[] t= g[i].getTranscripts();
					for (int j = 0; j < t.length; j++) {
						t[j].setGene(g[i]);		// TODO repair this !!!
						if (!t[j].isCoding())
							continue;
						String s= t[j].getTranslations()[0].translate();
						if (s== null)
							continue;
						writer.write(">"+t[j].getTranscriptID()+"\n");
						for (int x = 0; x < s.length(); x+=50) 
							writer.write(s.substring(x, Math.min(x+50, s.length()))+"\n");						
					}
				}
				
				reader.read();
				g= reader.getGenes();
			}
			
			writer.flush();
			writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static int[] countExonDomOvl(Vector<Exon> eVector, HashMap regMap, IntVector lenCDom, IntVector lenSDom, IntVector lenTot, 
			Vector collectHits) {
		int counterHit= 0, counterNoHit= 0;		
		for (int j = 0; j < eVector.size(); j++) {
			lenTot.add(eVector.elementAt(j).getLength());
			Transcript[] trpts= eVector.elementAt(j).getTranscripts();
			boolean found= false;
			for (int k = 0; k < trpts.length; k++) {	
				Vector<DirectedRegion> v= (Vector) regMap.get(trpts[k].getTranscriptID());
				for (int m = 0; v!= null&& m < v.size(); m++) {
					if (eVector.elementAt(j).contains(v.elementAt(m))) {
						++counterHit;
						lenCDom.add(eVector.elementAt(j).getLength());
						collectHits.add(eVector.elementAt(j));
						found= true;
						break;
					}
				}
				if (found)
					break;
			}
			if (found)
				;	// break; why?
			else {
				++counterNoHit;
				lenSDom.add(eVector.elementAt(j).getLength());				
			}
		}
		
		return new int[] {counterHit, counterNoHit};

	}
	
	public static strictfp void _05_00_pfamDomains() {
			String speName= "human";
			String annoName= "RefSeq";
			String[] keywords= new String[] {"_xref"}; 
			
			System.out.print("Loading domains..");
			System.out.flush();
			DouglasDomainWrapper wrapper= new DouglasDomainWrapper(
					"douglas"+File.separator+"H.sapiens.All.mRNAs.with.CDS.pfam.parsed");
			try {
				wrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}		
			HashMap map= wrapper.getMap(); 
			Object[] o= map.keySet().toArray();
			String[] someIDs= new String[o.length];
			int cntDomains= 0;
			for (int i = 0; i < o.length; i++) { 
	//			System.out.print(o[i]+"\t");
				Vector v= (Vector) map.get(o[i]); 
				cntDomains+= v.size(); 
	//			for (int j = 0; j < v.size(); j++) 
	//				System.out.print(v.elementAt(j)+",");
	//			System.out.println();
				someIDs[i]= (String) o[i];
			}
			System.out.println("found "+someIDs.length+" genes, "+cntDomains+" domains.\n");
			
			
				// get basic transcripts
			String absFName= Species.getAnnotation(speName, null, annoName, keywords); 
				//Constants.getLatestUCSCAnnotation(speName, annoName, keywords);
			System.out.println("Searching for reference transcripts in "+absFName+".");
			GTFChrReader reader= new GTFChrReader(absFName);
			reader.setChromosomeWise(false);
			reader.setFiltSomeIDs(someIDs);
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			Gene[] genes= reader.getGenes();
			int cntMultiTrpt= 0;
			int totNTcovByORFs= 0;
			for (int i = 0; i < genes.length; i++) {
				totNTcovByORFs+= genes[i].getTranscripts()[0].getTranslations()[0].getLength();
				if (genes[i].getTranscriptCount()> 1)
					++cntMultiTrpt;
			}
			System.out.println("Found "+genes.length+" genes ("+cntMultiTrpt+"multi-transcript).\n");		
	//		System.out.print("Not found: ");
	//		String[] notFound= reader.getFiltSomeIDsNotFound();
	//		for (int i = 0; i < notFound.length; i++) 
	//			System.out.print(notFound[i]+" ");
	//		System.out.println();
			
			
			
				// map AA -> NT
			System.out.print("Mapping domain info to genome..");
			System.out.flush();
			Vector vRegs= new Vector();
			HashMap mapProtTID= new HashMap(someIDs.length);
			int totNTcovByDom= 0;
			for (int i = 0; i < someIDs.length; i++) {
				int x;
				for (x = 0; x < genes.length; x++) {
					//assert(genes[x].getTranscripts().length== 1);
					String[] protIDs= genes[x].getTranscripts()[0].getTranslations()[0].getProteinIDsAll();
					int k;
					for (k = 0; k < protIDs.length; k++) 
						if (protIDs[k].equals(someIDs[i]))
							break;
					if (k< protIDs.length)
						break;
				}
				if (x== genes.length) 
					continue;
				mapProtTID.put(someIDs[i], genes[x].getTranscripts()[0].getTranscriptID());
				Vector v= (Vector) map.get(someIDs[i]);
				for (int j = 0; j < v.size(); j++) {
					DefaultRegion reg= (DefaultRegion) v.elementAt(j);
					Transcript trpt= genes[x].getTranscripts()[0];				
					int start= trpt.getTranslations()[0].getGenomicPosition((reg.getStart()- 1)*3);		// TODO check for 0-based
					int end= trpt.getTranslations()[0].getGenomicPosition((reg.getEnd())*3);		// -1 +1
					DirectedRegion genReg= new DirectedRegion(start, end, trpt.getStrand());
					genReg.setChromosome(trpt.getChromosome());
					genReg.setID(reg.getID());
					genReg.setScore(reg.getScore());
					vRegs.add(genReg);
					totNTcovByDom+= genReg.getLength();
				}
			}
			System.out.println("mapped "+vRegs.size()+" domains.");
			
			
			
				// get all transcripts in the regions
			System.out.print("Getting all transcripts in locus of domain..");
			System.out.flush();
			reader= new GTFChrReader(Constants.getLatestUCSCAnnotation(speName, annoName, null));
			DirectedRegion[] diregs= (DirectedRegion[]) Arrays.toField(vRegs);
			reader.setFiltRegs(diregs);
			reader.setSilent(false);
			reader.setChromosomeWise(false);
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
			int cntDomSingle= 0, cntSingle= 0;
			for (int i = 0; i < someIDs.length; i++) {	// another time, map to protIDs
				String tid= (String) mapProtTID.get(someIDs[i]);
				int x;
				for (x = 0; x < genes.length; x++) {
					int j;
					for (j = 0; j < genes[x].getTranscripts().length; j++) {
						if (genes[x].getTranscripts()[j].getTranscriptID().equals(tid))
							break;
					}
					if (j< genes[x].getTranscripts().length)
						break;
				}
				if (x== genes.length) 
					continue;
				
				int a= 0;
				for (int j = 0; j < genes[x].getTranscripts().length; j++) 
					if (genes[x].getTranscripts()[j].isCoding())
						++a;
				int b= 0;
				for (int j = 0; j < genes[x].getExons().length; j++) 
					if (genes[x].getExons()[j].isCoding())
						++b;
				
				if (a< 2|| b< 2) {
					++cntSingle;
					cntDomSingle+= ((Vector) map.get(someIDs[i])).size();
				}
			}
			System.out.println("found "+(genes.length- cntSingle)+" genes ("+cntDomSingle+
					" domains) with multiple transcripts/exons.\n");
			
				// retrieve events
			System.out.print("Retrieving events..");
			System.out.flush();
			HashMap localMap= new HashMap();
			Vector vASNDom= new Vector();
			int cntASgenes= 0, cntASEvents= 0, cntASinDom= 0;
			HashMap mapLandscape= new HashMap(), mapLandscapeNonDom= new HashMap(), mapLandscapeAll= new HashMap();
			for (int i = 0; i < genes.length; i++) { 
				ASVariation[][] vars= Analyzer.getASVariations(genes[i], ASMultiVariation.FILTER_HIERARCHICALLY, ASVariation.TYPE_ALL);
				int oldCntASEvents= cntASEvents;
				for (int x = 0; vars!= null&& x < vars.length; x++) {
					for (int xx = 0; xx < vars[x].length; xx++) {
						if (!vars[x][xx].isContainedCDS())
							continue;
						++cntASEvents;
						DirectedRegion asReg= vars[x][xx].getRegion();
						int y;
						for (y = 0; y < diregs.length; y++) 
							//if (asReg.overlaps(diregs[y])) {
							if (diregs[y].contains(asReg)) {
								++cntASinDom;
								Vector v= (Vector) localMap.get(diregs[y]);
								if (v== null)
									v= new Vector();
								v.add(vars[x][xx]);
								localMap.put(diregs[y], v);
								
								v= (Vector) mapLandscape.get(vars[x][xx].toString());
								if (v== null)
									v= new Vector();
								v.add(vars[x][xx]);
								mapLandscape.put(vars[x][xx].toString(), v);
								break;
							}
						if (y== diregs.length) {
							vASNDom.add(vars[x]);
							
							Vector v= (Vector) mapLandscapeNonDom.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscapeNonDom.put(vars[x][xx].toString(), v);
						}
						Vector v= (Vector) mapLandscapeAll.get(vars[x][xx].toString());
						if (v== null)
							v= new Vector();
						v.add(vars[x][xx]);
						mapLandscapeAll.put(vars[x][xx].toString(), v);
					}
				}
				if (cntASEvents!= oldCntASEvents)
					++cntASgenes;
			}
			System.out.println("found "+cntASgenes+" genes with AS.\n\n");
			
			
			System.out.println("Domain\tevents");
			o= localMap.keySet().toArray();
			int cntDomAS= 0;
			for (int i = 0; i < o.length; i++) {
				System.out.print(o[i]+"\t");
				Vector v= (Vector) localMap.get(o[i]);
				if (v== null|| v.size()== 0) 
					continue;
				++cntDomAS;
				for (int j = 0; j < v.size(); j++) 
					System.out.print(v.elementAt(j)+",");
				System.out.println();
			}
			System.out.println(cntDomAS+" domains containing at least 1 AS event.\n");
			
			System.out.println("Landscape Domains");
			Vector vv= new Vector();
			Object[] keys= mapLandscape.keySet().toArray();
			HashMap mapDiffLenDom= new HashMap(mapLandscape.size());
			IntVector diffAllDom= new IntVector();
			int cntAsymDon= 0, cntAsymAcc= 0, cntDonExtExon= 0, cntDonRedExon= 0,
				cntAccExtExon= 0, cntAccRedExon= 0, cntSmallDonStronger= 0, cntBigDonStronger= 0,
				cntSmallAccStronger= 0, cntBigAccStronger= 0;
			IntVector donExt= new IntVector(), donRed= new IntVector(), accExt= new IntVector(),
				accRed= new IntVector();
			Noboru nobi= new Noboru(speName, annoName);
			HashMap mapLittleDonCodonMap= new HashMap(), mapGrandeDonCodonMap= new HashMap(), 
				mapLittleAccCodonMap= new HashMap(), mapGrandeAccCodonMap= new HashMap();
			Vector transAAdonExt= new Vector(), transAAdonRed= new Vector(), transAAaccExt= new Vector(), transAAaccRed= new Vector();
			DoubleVector percStopDE= new DoubleVector(), percStopDR= new DoubleVector(), percStopAE= new DoubleVector(), percStopAR= new DoubleVector(),   
			percStericDE= new DoubleVector(), percStericDR= new DoubleVector(), percStericAE= new DoubleVector(), percStericAR= new DoubleVector(), 
			percPolarDE= new DoubleVector(), percPolarDR= new DoubleVector(), percPolarAE= new DoubleVector(), percPolarAR= new DoubleVector(), 
			percPosDE= new DoubleVector(), percPosDR= new DoubleVector(), percPosAE= new DoubleVector(), percPosAR= new DoubleVector(), 
			percNegDE= new DoubleVector(), percNegDR= new DoubleVector(), percNegAE= new DoubleVector(), percNegAR= new DoubleVector(); 
			for (int i = 0; i < keys.length; i++) { 
				vv.add(mapLandscape.get(keys[i]));
				Vector v= (Vector) mapLandscape.get(keys[i]);
				for (int j = 0; j < v.size(); j++) {
					ASVariation var= (ASVariation) v.elementAt(j);
					int val= var.getLengthDiff(true);
					IntVector ivec= (IntVector) mapDiffLenDom.get(keys[i]);
					if (ivec== null)
						ivec= new IntVector();
					ivec.add(val);
					diffAllDom.add(val);
					mapDiffLenDom.put(keys[i], ivec);
					
					if (keys[i].equals("1^ , 2^")) {
						if (val% 3!= 0)
							++cntAsymDon;
						int posLittle= var.getSpliceUniverse()[0].getPos();
						int posGrande= var.getSpliceUniverse()[1].getPos();
						Gene g= var.getTranscript1().getGene();
						int cntLittle= 0, cntGrande= 0;
						for (int k = 0; k < g.getTranscripts().length; k++) 
							for (int m = 0; m < g.getTranscripts()[k].getSpliceChain().length; m++) 
								if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posLittle) {
									++cntLittle;
									break;
								} else if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posGrande) {
									++cntGrande;
									break;
								}
						if (cntLittle> cntGrande) {
							++cntDonExtExon;
							donExt.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								Integer cnt= (Integer) mapLittleDonCodonMap.get(codons[k]);
								aa+= Translation.getTranslatedAA(codons[k]);
								if (cnt== null)
									mapLittleDonCodonMap.put(codons[k], new Integer(1));
								else
									mapLittleDonCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAdonExt.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopDE.add((double) cntStop/ aa.length());
							percStericDE.add((double) cntSteric/ aa.length());
							percPolarDE.add((double) cntPolar/ aa.length());
							percPosDE.add((double) cntPos/ aa.length());
							percNegDE.add((double) cntNeg/ aa.length());
						} else if (cntGrande> cntLittle) {
							++cntDonRedExon;
							donRed.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								Integer cnt= (Integer) mapGrandeDonCodonMap.get(codons[k]);
								aa+= Translation.getTranslatedAA(codons[k]);
								if (cnt== null)
									mapGrandeDonCodonMap.put(codons[k], new Integer(1));
								else
									mapGrandeDonCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAdonRed.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopDR.add((double) cntStop/ aa.length());
							percStericDR.add((double) cntSteric/ aa.length());
							percPolarDR.add((double) cntPolar/ aa.length());
							percPosDR.add((double) cntPos/ aa.length());
							percNegDR.add((double) cntNeg/ aa.length());
						}
						
	//					double[] strengths= nobi.scoreSpliceSite(var.getSpliceUniverse());
	//					if (strengths[0]> strengths[1])
	//						++cntSmallDonStronger;
	//					else if (strengths[0]< strengths[1])
	//						++cntBigDonStronger;
					}
					if (keys[i].equals("1- , 2-")) {
						if (val% 3!= 0)
							++cntAsymAcc;
						int posLittle= var.getSpliceUniverse()[0].getPos();
						int posGrande= var.getSpliceUniverse()[1].getPos();
						Gene g= var.getTranscript1().getGene();
						int cntLittle= 0, cntGrande= 0;
						for (int k = 0; k < g.getTranscripts().length; k++) 
							for (int m = 0; m < g.getTranscripts()[k].getSpliceChain().length; m++) 
								if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posLittle) {
									++cntLittle;
									break;
								} else if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posGrande) {
									++cntGrande;
									break;
								}
						if (cntLittle> cntGrande) {
							++cntAccExtExon;
							accExt.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								aa+= Translation.getTranslatedAA(codons[k]);
								Integer cnt= (Integer) mapLittleAccCodonMap.get(codons[k]);
								if (cnt== null)
									mapLittleAccCodonMap.put(codons[k], new Integer(1));
								else
									mapLittleAccCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAaccExt.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopAE.add((double) cntStop/ aa.length());
							percStericAE.add((double) cntSteric/ aa.length());
							percPolarAE.add((double) cntPolar/ aa.length());
							percPosAE.add((double) cntPos/ aa.length());
							percNegAE.add((double) cntNeg/ aa.length());
						} else if (cntGrande> cntLittle) {
							++cntAccRedExon;
							accRed.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								Integer cnt= (Integer) mapGrandeAccCodonMap.get(codons[k]);
								aa+= Translation.getTranslatedAA(codons[k]);
								if (cnt== null)
									mapGrandeAccCodonMap.put(codons[k], new Integer(1));
								else
									mapGrandeAccCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAaccRed.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopAR.add((double) cntStop/ aa.length());
							percStericAR.add((double) cntSteric/ aa.length());
							percPolarAR.add((double) cntPolar/ aa.length());
							percPosAR.add((double) cntPos/ aa.length());
							percNegAR.add((double) cntNeg/ aa.length());
						}
	
	//					double[] strengths= nobi.scoreSpliceSite(var.getSpliceUniverse());
	//					if (strengths[0]> strengths[1])
	//						++cntSmallAccStronger;
	//					else if (strengths[0]< strengths[1])
	//						++cntBigAccStronger;
					}
					
				}
				int[] distArray= ((IntVector) mapDiffLenDom.get(keys[i])).toIntArray();
				java.util.Arrays.sort(distArray);
				Distribution distr= new Distribution(distArray);
				System.out.print(keys[i]+"\t"+distr.getMedian()+"\t");
				for (int j = 0; j < distArray.length; j++) 
					System.out.print(distArray[j]+",");
				System.out.println();
			}
			Distribution distr= new Distribution(diffAllDom.toIntArray());
			System.out.println("ALL\t"+distr.getMedian());
			System.out.println("asymDon: "+cntAsymDon+", asymAcc: "+cntAsymAcc+ "\t"+
					"don ext "+cntDonExtExon+", don red "+cntDonRedExon+", acc ext "+
					cntAccExtExon+", acc red "+cntAccRedExon+ ";\t+smDon "+cntSmallDonStronger+
					", +bigDon "+cntBigDonStronger+", +smAcc "+cntSmallAccStronger+", +bigAcc "+
					cntBigAccStronger);
			Distribution dist; 
			int[] h= donExt.toIntArray();
			dist= new Distribution(h);
			System.out.print(h.length+" donExt "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			Object[] kys= mapLittleDonCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapLittleDonCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAdonExt.size(); i++) 
				System.out.print(transAAdonExt.elementAt(i)+", ");
			System.out.println();
			double[] hull= percStopDE.toDoubleArray();
			dist= new Distribution(hull);
			int cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			
			h= donRed.toIntArray();
			dist= new Distribution(h);
			System.out.println(h.length+" donRed "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapGrandeDonCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapGrandeDonCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAdonRed.size(); i++) 
				System.out.print(transAAdonRed.elementAt(i)+", ");
			System.out.println();
			hull= percStopDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			
			h= accExt.toIntArray();
			dist= new Distribution(h);
			System.out.println(h.length+" accExt "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapLittleAccCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapLittleAccCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAaccExt.size(); i++) 
				System.out.print(transAAaccExt.elementAt(i)+", ");
			System.out.println();
			hull= percStopAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
			h= accRed.toIntArray();
			System.out.println(h.length+" accRed "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapGrandeAccCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapGrandeAccCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAaccRed.size(); i++) 
				System.out.print(transAAaccRed.elementAt(i)+", ");
			System.out.println();
			hull= percStopAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
					
			ASVariation[][] vars= (ASVariation[][]) Arrays.toField(vv);
			Arrays.sort2DFieldRev(vars);
			ASAnalyzer.outputVariations(vars, false, false, System.out);
			System.out.println("\n");
			
			System.out.println("Landscape NonDomains");
			vv= new Vector();
			keys= mapLandscapeNonDom.keySet().toArray();
			HashMap mapDiffLenNDom= new HashMap(mapLandscapeNonDom.size());
			IntVector diffAllNDom= new IntVector();
			cntAsymDon= 0; cntAsymAcc= 0; cntDonExtExon= 0; cntDonRedExon= 0;
			cntAccExtExon= 0; cntAccRedExon= 0;  cntSmallDonStronger= 0; cntBigDonStronger= 0;
			cntSmallAccStronger= 0; cntBigAccStronger= 0;
			donExt= new IntVector(); donRed= new IntVector(); accExt= new IntVector();
			accRed= new IntVector();
			mapLittleDonCodonMap= new HashMap(); mapGrandeDonCodonMap= new HashMap(); 
			mapLittleAccCodonMap= new HashMap(); mapGrandeAccCodonMap= new HashMap();
			transAAdonExt= new Vector(); transAAdonRed= new Vector(); transAAaccExt= new Vector(); transAAaccRed= new Vector();
			percStopDE= new DoubleVector(); percStopDR= new DoubleVector(); percStopAE= new DoubleVector(); percStopAR= new DoubleVector();   
			percStericDE= new DoubleVector(); percStericDR= new DoubleVector(); percStericAE= new DoubleVector(); percStericAR= new DoubleVector(); 
			percPolarDE= new DoubleVector(); percPolarDR= new DoubleVector(); percPolarAE= new DoubleVector(); percPolarAR= new DoubleVector(); 
			percPosDE= new DoubleVector(); percPosDR= new DoubleVector(); percPosAE= new DoubleVector(); percPosAR= new DoubleVector(); 
			percNegDE= new DoubleVector(); percNegDR= new DoubleVector(); percNegAE= new DoubleVector(); percNegAR= new DoubleVector(); 
			for (int i = 0; i < keys.length; i++) { 
				vv.add(mapLandscapeNonDom.get(keys[i]));
				Vector v= (Vector) mapLandscapeNonDom.get(keys[i]);
				for (int j = 0; j < v.size(); j++) {
					ASVariation var= (ASVariation) v.elementAt(j);
					int val= var.getLengthDiff(true);
					IntVector ivec= (IntVector) mapDiffLenNDom.get(keys[i]);
					if (ivec== null)
						ivec= new IntVector();
					ivec.add(val);
					diffAllNDom.add(val);
					mapDiffLenNDom.put(keys[i], ivec);
	
					if (keys[i].equals("1^ , 2^")) {
						if (val% 3!= 0)
							++cntAsymDon;
						int posLittle= var.getSpliceUniverse()[0].getPos();
						int posGrande= var.getSpliceUniverse()[1].getPos();
						Gene g= var.getTranscript1().getGene();
						int cntLittle= 0, cntGrande= 0;
						for (int k = 0; k < g.getTranscripts().length; k++) 
							for (int m = 0; m < g.getTranscripts()[k].getSpliceChain().length; m++) 
								if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posLittle) {
									++cntLittle;
									break;
								} else if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posGrande) {
									++cntGrande;
									break;
								}
						if (cntLittle> cntGrande) {
							++cntDonExtExon;
							donExt.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								Integer cnt= (Integer) mapLittleDonCodonMap.get(codons[k]);
								aa+= Translation.getTranslatedAA(codons[k]);
								if (cnt== null)
									mapLittleDonCodonMap.put(codons[k], new Integer(1));
								else
									mapLittleDonCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAdonExt.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopDE.add((double) cntStop/ aa.length());
							percStericDE.add((double) cntSteric/ aa.length());
							percPolarDE.add((double) cntPolar/ aa.length());
							percPosDE.add((double) cntPos/ aa.length());
							percNegDE.add((double) cntNeg/ aa.length());
						} else if (cntGrande> cntLittle) {
							++cntDonRedExon;
							donRed.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								aa+= Translation.getTranslatedAA(codons[k]);
								Integer cnt= (Integer) mapGrandeDonCodonMap.get(codons[k]);
								if (cnt== null)
									mapGrandeDonCodonMap.put(codons[k], new Integer(1));
								else
									mapGrandeDonCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAdonRed.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopDR.add((double) cntStop/ aa.length());
							percStericDR.add((double) cntSteric/ aa.length());
							percPolarDR.add((double) cntPolar/ aa.length());
							percPosDR.add((double) cntPos/ aa.length());
							percNegDR.add((double) cntNeg/ aa.length());
						}
	
	//					double[] strengths= nobi.scoreSpliceSite(var.getSpliceUniverse());
	//					if (strengths[0]> strengths[1])
	//						++cntSmallDonStronger;
	//					else if (strengths[0]< strengths[1])
	//						++cntBigDonStronger;
					}
					if (keys[i].equals("1- , 2-")) {
						if (val% 3!= 0)
							++cntAsymAcc;
						int posLittle= var.getSpliceUniverse()[0].getPos();
						int posGrande= var.getSpliceUniverse()[1].getPos();
						Gene g= var.getTranscript1().getGene();
						int cntLittle= 0, cntGrande= 0;
						for (int k = 0; k < g.getTranscripts().length; k++) 
							for (int m = 0; m < g.getTranscripts()[k].getSpliceChain().length; m++) 
								if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posLittle) {
									++cntLittle;
									break;
								} else if (g.getTranscripts()[k].getSpliceChain()[m].getPos()== posGrande) {
									++cntGrande;
									break;
								}
						if (cntLittle> cntGrande) {
							++cntAccExtExon;
							accExt.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								aa+= Translation.getTranslatedAA(codons[k]);
								Integer cnt= (Integer) mapLittleAccCodonMap.get(codons[k]);
								if (cnt== null)
									mapLittleAccCodonMap.put(codons[k], new Integer(1));
								else
									mapLittleAccCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAaccExt.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopAE.add((double) cntStop/ aa.length());
							percStericAE.add((double) cntSteric/ aa.length());
							percPolarAE.add((double) cntPolar/ aa.length());
							percPosAE.add((double) cntPos/ aa.length());
							percNegAE.add((double) cntNeg/ aa.length());
						} else if (cntGrande> cntLittle) { 
							++cntAccRedExon;
							accRed.add(val);
							String[] codons= var.getVariableCodons();
							String aa= "";
							for (int k = 0; k < codons.length; k++) {
								codons[k]= codons[k].toUpperCase();
								aa+= Translation.getTranslatedAA(codons[k]);
								Integer cnt= (Integer) mapGrandeAccCodonMap.get(codons[k]);
								if (cnt== null)
									mapGrandeAccCodonMap.put(codons[k], new Integer(1));
								else
									mapGrandeAccCodonMap.put(codons[k], new Integer(cnt.intValue()+ 1));
							}
							transAAaccRed.add(aa);
							int cntStop= 0, cntSteric= 0, cntPolar= 0, cntPos= 0, cntNeg= 0; 
							for (int k = 0; k < aa.length(); k++) {
								if (Translation.isStop(aa.charAt(k)))
									++cntStop;
								if (Translation.isAAsteric(aa.charAt(k)))
									++cntSteric;
								if (Translation.isAApolar(aa.charAt(k)))
									++cntPolar;
								if (Translation.isAAposCharge(aa.charAt(k)))
									++cntPos;
								if (Translation.isAAnegCharge(aa.charAt(k)))
									++cntNeg;
							}
							percStopAR.add((double) cntStop/ aa.length());
							percStericAR.add((double) cntSteric/ aa.length());
							percPolarAR.add((double) cntPolar/ aa.length());
							percPosAR.add((double) cntPos/ aa.length());
							percNegAR.add((double) cntNeg/ aa.length());
						}
	
	//					double[] strengths= nobi.scoreSpliceSite(var.getSpliceUniverse());
	//					if (strengths[0]> strengths[1])
	//						++cntSmallAccStronger;
	//					else if (strengths[0]< strengths[1])
	//						++cntBigAccStronger;
					}
					
				}
				int[] distArray= ((IntVector) mapDiffLenNDom.get(keys[i])).toIntArray();
				java.util.Arrays.sort(distArray);
				distr= new Distribution(distArray);
				System.out.print(keys[i]+"\t"+distr.getMedian()+"\t");
				for (int j = 0; j < distArray.length; j++) 
					System.out.print(distArray[j]+",");
				System.out.println();
			}
			distr= new Distribution(diffAllNDom.toIntArray());
			System.out.println("ALL\t"+distr.getMedian());
			System.out.println("asymDon: "+cntAsymDon+", asymAcc: "+cntAsymAcc+ "\t"+
					"don ext "+cntDonExtExon+", don red "+cntDonRedExon+", acc ext "+
					cntAccExtExon+", acc red "+cntAccRedExon+ ";\t+smDon "+cntSmallDonStronger+
					", +bigDon "+cntBigDonStronger+", +smAcc "+cntSmallAccStronger+", +bigAcc "+
					cntBigAccStronger);
			h= donExt.toIntArray();
			dist= new Distribution(h);
			System.out.print(h.length+" donExt "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapLittleDonCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapLittleDonCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAdonExt.size(); i++) 
				System.out.print(transAAdonExt.elementAt(i)+", ");
			System.out.println();
			hull= percStopDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosDE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegDE.toDoubleArray();
			dist= new Distribution(hull);
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
			h= donRed.toIntArray();
			dist= new Distribution(h);
			System.out.println(h.length+" donRed "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapGrandeDonCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapGrandeDonCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAdonRed.size(); i++) 
				System.out.print(transAAdonRed.elementAt(i)+", ");
			System.out.println();
			hull= percStopDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegDR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
			h= accExt.toIntArray();
			dist= new Distribution(h);
			System.out.println(h.length+" accExt "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapLittleAccCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapLittleAccCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAaccExt.size(); i++) 
				System.out.print(transAAaccExt.elementAt(i)+", ");
			System.out.println();
			hull= percStopAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegAE.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
			h= accRed.toIntArray();
			System.out.println(h.length+" accRed "+dist.getMedian()+" : ");
			java.util.Arrays.sort(h);
			for (int i = 0; i < h.length; i++) 
				System.out.print(h[i]+",");
			System.out.println();
			kys= mapGrandeAccCodonMap.keySet().toArray();
			java.util.Arrays.sort(kys);
			for (int i = 0; i < kys.length; i++) 
				System.out.print(kys[i]+" "+mapGrandeAccCodonMap.get(kys[i])+",");
			System.out.println();
			for (int i = 0; i < transAAaccRed.size(); i++) 
				System.out.print(transAAaccRed.elementAt(i)+", ");
			System.out.println();
			hull= percStopAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("stops ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percStericAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("steric ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPolarAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("polar ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percPosAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("pos ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
			hull= percNegAR.toDoubleArray();
			dist= new Distribution(hull);
			cc= 0;
			for (int i = 0; i < hull.length; i++) 
				if (hull[i]> 0d)
					++cc;
			System.out.print("neg ("+Formatter.fprint(100d* cc/ h.length, 2)+"% med "+Formatter.fprint(100d* dist.getMedian(), 2)+"%): ");
			for (int i = 0; i < hull.length; i++) 
				System.out.print(Formatter.fprint(hull[i]* 100d, 2)+"%, ");
			System.out.println();
	
			vars= (ASVariation[][]) Arrays.toField(vv);
			Arrays.sort2DFieldRev(vars);
			ASAnalyzer.outputVariations(vars, false, false, System.out);
			System.out.println("\n");
			
			System.out.println("Landscape All");
			vv= new Vector();
			keys= mapLandscapeAll.keySet().toArray();
			for (int i = 0; i < keys.length; i++) 
				vv.add(mapLandscapeAll.get(keys[i]));
			vars= (ASVariation[][]) Arrays.toField(vv);
			Arrays.sort2DFieldRev(vars);
			ASAnalyzer.outputVariations(vars, false, false, System.out);
			System.out.println("\n");
	
			System.out.println();
			float percGenSingle= ((genes.length- cntSingle)*100f)/ genes.length;
			float percGenAS= (cntASgenes*100f)/ (genes.length- cntSingle);
			System.out.println(Formatter.fprint(percGenSingle, 2)+"% multi transcript genes, of them "+
					Formatter.fprint(percGenAS, 2)+"% with AS.");
			
			float percDomSingle= ((cntDomains- cntDomSingle)*100f)/ cntDomains;  
			float percDomAS= (cntDomAS*100f)/ (cntDomains- cntDomSingle);  
			System.out.println(Formatter.fprint(percDomSingle, 2)+"% domains in loci with multiple transcripts, of them "+
					Formatter.fprint(percDomAS, 2)+"% with AS.");
			
			float percASinDom= (cntASinDom*100f)/ cntASEvents;  
			float percNTcov= (totNTcovByDom*100f)/ totNTcovByORFs;  
			System.out.println(Formatter.fprint(percASinDom, 2)+"% of AS events are in domains, which occupy a fraction of NT of "+
					Formatter.fprint(percNTcov, 2)+"%");
			
	
		}

		public static strictfp void _05_00_pfamDomains_all() {
			long t0= System.currentTimeMillis();
			String input= "H.sapiens.All.mRNAs.with.CDS.pfam.parsed";
			String speName= "human";
			String annoName= "RefSeq";
			String[] keywords= new String[] {"_xref"}; 
			
			System.out.print("Loading domains..");
			System.out.flush();
			DomainWrapper wrapper= new DomainWrapper(
					"douglas"+File.separator+input);
			try {
				wrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}		
			HashMap map= wrapper.getMap(); 
			Object[] o= map.keySet().toArray();
			String[] someIDs= new String[o.length];
			int cntDomains= 0;
			DoubleVector scores= new DoubleVector();
			for (int i = 0; i < o.length; i++) { 
	//			out.print(o[i]+"\t");
				Vector v= (Vector) map.get(o[i]); 
				cntDomains+= v.size(); 
				for (int j = 0; j < v.size(); j++) { 
	//				out.print(v.elementAt(j)+",");
					DefaultRegion reg= (DefaultRegion) v.elementAt(j);
					scores.add(reg.getScore());
				}
	//			out.println();
				someIDs[i]= (String) o[i];
			}		
			Distribution dist= new Distribution(scores.toDoubleArray());
			double median= dist.getMedian();
			double threshold= -1d;
			PrintStream out= null;
			try {
				out= new PrintStream("douglas"+File.separator+"pfam"+File.separator+
						speName+"_"+annoName+"_"+threshold+"."+input);
			} catch (FileNotFoundException e1) {			
				e1.printStackTrace();
			}
			out.println("found "+someIDs.length+" genes, "+cntDomains+" domains.\n");
			out.println("score med "+median+", min "+dist.getMin()+", max "+dist.getMax());
			
			
				// get basic transcripts
			String absFName= Species.getAnnotation(speName, null, annoName, keywords); 
				//Constants.getLatestUCSCAnnotation(speName, annoName, keywords);
			out.println("Searching for reference transcripts in "+absFName+".");
			GTFChrReader reader= new GTFChrReader(absFName);
			reader.setChromosomeWise(true);
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			Gene[] genes= reader.getGenes();
			
			int cntFoundTrpts= 0;
			int cntDomSingle= 0, cntSingle= 0;
			IntVector asGenesV= new IntVector(), asInDomV= new IntVector(), asInNonDomV= new IntVector(),
			asDomGenes= new IntVector(), asNonDomGenes= new IntVector();
			HashMap mapLandscapeDom= new HashMap(), mapLandscapeNonDom= new HashMap(), 
				mapLandscapeAllDomProts= new HashMap(), mapLandscapeAll= new HashMap(),
				mapLandscapeNonDomProts= new HashMap();
			int cntAllGenes= 0, cntNonCoding= 0, cntNonASAllGenes= 0, cntNonASDomGenes= 0, 
					cntNonASDom= 0, cntNonASNonDom= 0, cntNonASNonDomGenes= 0;
			HashMap<DirectedRegion,Vector> domainSpecMap= new HashMap<DirectedRegion,Vector>();
			HashMap<String,Vector<DirectedRegion>> domGroupMap= new HashMap<String,Vector<DirectedRegion>>();
			Vector<ASVariation> mutExDom= new Vector<ASVariation>();
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					++cntAllGenes;
					if (!genes[i].isProteinCoding())
						++cntNonCoding;
					Transcript[] trpts= genes[i].getTranscripts();
					boolean found= false;
					Vector diregV= new Vector();
					for (int j = 0; j < trpts.length; j++) {
						Vector v= (Vector) map.get(trpts[j].getTranscriptID());
						for (int k = 0; v!= null&& k < v.size(); k++) {
							DefaultRegion regD= (DefaultRegion) v.elementAt(k);
							if (threshold> 0d&& regD.getScore()> threshold)
								continue;
							++cntFoundTrpts;
							found= true;
							Translation tln= trpts[j].getTranslations()[0];
							DirectedRegion reg= new DirectedRegion();
							reg.setStrand(trpts[j].getStrand());
							reg.setChromosome(trpts[j].getChromosome());
							reg.setStart(tln.getGenomicPosition((regD.getStart()- 1)*3));		// TODO check for 0-based
							reg.setEnd(tln.getGenomicPosition((regD.getEnd())*3));		// -1 +1
							reg.setScore(regD.getScore());
							reg.setID(regD.getID());
							diregV.add(reg);
							
							Vector<DirectedRegion> vv= domGroupMap.get(reg.getID());
							if (vv== null)
								vv= new Vector<DirectedRegion>();
							vv.add(reg);
							domGroupMap.put(reg.getID(), vv);
						}
					}
		
						// events
					DirectedRegion[] diregs= null;
					if (found) 
						diregs= (DirectedRegion[]) Arrays.toField(diregV);
					ASVariation[][] vars= Analyzer.getASVariations(genes[i], ASMultiVariation.FILTER_HIERARCHICALLY, ASVariation.TYPE_ALL);
					int cntASEvents= 0, cntASinDom= 0, cntASinNonDom= 0, cntASinDomGenes= 0, cntASinNonDomGenes= 0;
					if (vars== null|| vars.length== 0) {
						++cntNonASAllGenes;
						if (found)
							++cntNonASDomGenes;
						else
							++cntNonASNonDomGenes;
					} 
					for (int x = 0; vars!= null&& x < vars.length; x++) {
						for (int xx = 0; xx < vars[x].length; xx++) {
							if (!vars[x][xx].isContainedCDS())
								continue;
							++cntASEvents;
							Vector v= (Vector) mapLandscapeAll.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscapeAll.put(vars[x][xx].toString(), v);
							if (!found) {
								v= (Vector) mapLandscapeNonDomProts.get(vars[x][xx].toString());
								if (v== null)
									v= new Vector();
								v.add(vars[x][xx]);
								mapLandscapeNonDomProts.put(vars[x][xx].toString(), v);
								++cntASinNonDomGenes;
								continue;
							}
							++cntASinDomGenes;
							DirectedRegion asReg= vars[x][xx].getRegion();
							int y;
							for (y = 0; diregs!= null&& y < diregs.length; y++) 
								if (asReg.overlaps(diregs[y])) {	// case of AS can be bigger than domain
								//if (diregs[y].contains(asReg)) {
								//if (diregs[y].contains(asReg)|| asReg.contains(diregs[y])) {
									++cntASinDom;
									
									String asID= vars[x][xx].toString();
									v= (Vector) mapLandscapeDom.get(asID);
									if (v== null)
										v= new Vector();
									v.add(vars[x][xx]);
									mapLandscapeDom.put(asID, v);
									
									v= (Vector) domainSpecMap.get(diregs[y]);
									if (v== null)
										v= new Vector();
									v.add(vars[x][xx]);
									domainSpecMap.put(diregs[y], v);
									
									if(asID.equals(ASVariation.ID_MUTEX))
										mutExDom.add(vars[x][xx]);
									break;	// fucks up the domainSpecMap
								}
							if (y== diregs.length) {	// hit in non-domain
								++cntASinNonDom;
								v= (Vector) mapLandscapeNonDom.get(vars[x][xx].toString());
								if (v== null)
									v= new Vector();
								v.add(vars[x][xx]);
								mapLandscapeNonDom.put(vars[x][xx].toString(), v);
							}
							v= (Vector) mapLandscapeAllDomProts.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscapeAllDomProts.put(vars[x][xx].toString(), v);
						}
					}
					
					if (vars!= null&& vars.length> 0) {
						asGenesV.add(cntASEvents);
						if (found) {
							asDomGenes.add(cntASinDomGenes);
							for (int j = 0; j < diregs.length; j++) {
								Vector v= domainSpecMap.get(diregs[j]);
								if (v== null)
									++cntNonASDom;
								else
									asInDomV.add(v.size());
							}
							if (cntASinNonDom== 0)
								++cntNonASNonDom;
							else
								asInNonDomV.add(cntASinNonDom); 	// can be 0
						} else
							asNonDomGenes.add(cntASinNonDomGenes);
					} else if (found) {
						++cntNonASDom;
						++cntNonASNonDom;
					}
	
				}
				
				try {
					reader.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				genes= reader.getGenes();
			}
			
				// output
			out.println(cntAllGenes+" genes, "+cntNonCoding+" non-coding genes.");
			out.println(
				cntNonASAllGenes+ " non AS genes ("+Formatter.fprint(100d* cntNonASAllGenes/cntAllGenes, 2)+"%)\n"+
				(cntNonASDomGenes+asDomGenes.length)+ " genes with annotated domains, "+
					asDomGenes.length+" alt. spliced ("+ Formatter.fprint(100d* asDomGenes.length/ (cntNonASDomGenes+ asDomGenes.length), 2)+"%)\n"+
				(cntNonASNonDomGenes+asNonDomGenes.length)+ " genes w/no annotated domains, "+
					asNonDomGenes.length+" alt. spliced ("+ Formatter.fprint(100d* asNonDomGenes.length/ (cntNonASNonDomGenes+ asNonDomGenes.length), 2)+"%)\n"+
				(cntNonASDom+ asInDomV.length)+" domains, "+
					asInDomV.length+" alt. spliced ("+ Formatter.fprint(100d* asInDomV.length/(cntNonASDom+ asInDomV.length), 2)+ "%)\n"+
				(cntNonASNonDom+ asInNonDomV.length)+" genes w domains, "+
					asInNonDomV.length+" alt. spliced outside of domains ("+ Formatter.fprint(100d* asInNonDomV .length/(cntNonASNonDom+ asInNonDomV.length), 2)+ "%)\n"
			);
			
			
			out.println("\n==================");
			out.println("All Genes:");
			dist= new Distribution(asGenesV.toIntArray());
			out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
			Analyzer.outputLandscape(mapLandscapeAll, out);
			
			out.println("\n==================");
			out.println("All Domain Genes:");
			dist= new Distribution(asDomGenes.toIntArray());
			out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
			Analyzer.outputLandscape(mapLandscapeAllDomProts, out);
			
			out.println("\n==================");
			out.println("All Non-Domain Genes:");
			dist= new Distribution(asNonDomGenes.toIntArray());
			out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
			Analyzer.outputLandscape(mapLandscapeNonDomProts, out);
	
			out.println("\n==================");
			out.println("in Domains:");
			dist= new Distribution(asInDomV.toIntArray());
			out.println(dist.getSum()+" events, med "+dist.getMedian()+" per domain.");
			Analyzer.outputLandscape(mapLandscapeDom, out);
			for (int i = 0; i < mutExDom.size(); i++) 
				System.out.println(ASVariation.ID_MUTEX+"\t"+mutExDom.elementAt(i).toCoordinates());
	
			out.println("\n==================");
			out.println("in Non-Domains:");
			dist= new Distribution(asInNonDomV.toIntArray());
			out.println(dist.getSum()+" events, med "+dist.getMedian()+" per non-domain.");
			Analyzer.outputLandscape(mapLandscapeNonDom, out);
	
	
			HashMap<String,Vector> domGroupEventMap= new HashMap<String,Vector>(domGroupMap.size());
			Iterator<String> iter= domGroupMap.keySet().iterator();
			while(iter.hasNext()) {
				String oID= (String) iter.next();
				Vector<ASVariation> vv= domGroupEventMap.get(oID);
				if (vv== null)
					vv= new Vector<ASVariation>();
				Vector<DirectedRegion> v= (Vector<DirectedRegion>) domGroupMap.get(oID);
				for (int i = 0; i < v.size(); i++) { 
					Vector<ASVariation> inV= domainSpecMap.get(v.elementAt(i));	// directedRegx V(ASevents)
					if (inV== null)
						continue;
					for (int j = 0; j < inV.size(); j++) 
						vv.add(inV.elementAt(j));
					if (vv.size()> 0)
						domGroupEventMap.put(oID, vv);
				}
				
			}
			out.println("\n\nDomain Group Events");
			DualHashBidiMap mapD= new DualHashBidiMap(domGroupEventMap);
			Object[] values= domGroupEventMap.values().toArray();
			Comparator compi= new VectorSizeComparator();
			java.util.Arrays.sort(values, compi);
			BidiMap revMap= mapD.inverseBidiMap();
			for (int i = values.length- 1; i >= 0; --i) {
				Vector<ASVariation> varV= (Vector<ASVariation>) values[i];
				out.print(varV.size()+ "\t"+ revMap.get(values[i])+"\t|");
				for (int j = 0; j < varV.size(); j++) 
					out.print(varV.elementAt(j)+"|");
				out.println();
			}
			
				
			out.println("\n\nDomain Groups");
			mapD= new DualHashBidiMap(domGroupMap);
			values= domGroupMap.values().toArray();
			java.util.Arrays.sort(values, compi);
			revMap= mapD.inverseBidiMap();
			for (int i = values.length- 1; i >= 0; --i) {
				Vector<ASVariation> varV= (Vector<ASVariation>) values[i];
				out.print(varV.size()+ "\t"+ revMap.get(values[i])+"\t|");
	//			for (int j = 0; j < varV.size(); j++) 
	//				out.print(varV.elementAt(j)+"|");
				out.println();
			}
			
			out.println("\n\nSingle Domains");
			IntVector cntEvPerDomain= new IntVector();
			Iterator<Vector> iterVec= domainSpecMap.values().iterator();
			while (iter.hasNext()) 
				cntEvPerDomain.add(iterVec.next().size());		
			dist= new Distribution(cntEvPerDomain.toIntArray());
			out.println("events per domain: med "+dist.getMedian()+", min "+dist.getMin()+", max "+dist.getMax());
			Analyzer.outputLandscape(domainSpecMap, out);
			
			out.println("\n\n[took "+(System.currentTimeMillis()- t0)/1000+" sec.]");
			out.flush();
			out.close();
		}

		public static strictfp void _05_exonsStatistic() {
			long t0= System.currentTimeMillis();
			String input= "H.sapiens.All.mRNAs.with.CDS.pfam.parsed";
			String speName= "human";
			String annoName= "RefSeq";
			String[] keywords= new String[] {"_xref"}; 
			
			System.out.print("Loading domains..");
			System.out.flush();
			DomainWrapper wrapper= new DomainWrapper(
					"douglas"+File.separator+input);
			try {
				wrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}		
			HashMap map= wrapper.getMap(); 
			Object[] o= map.keySet().toArray();
			String[] someIDs= new String[o.length];
			int cntDomains= 0;
			DoubleVector scores= new DoubleVector();
			for (int i = 0; i < o.length; i++) { 
	//			out.print(o[i]+"\t");
				Vector v= (Vector) map.get(o[i]); 
				cntDomains+= v.size(); 
				for (int j = 0; j < v.size(); j++) { 
	//				out.print(v.elementAt(j)+",");
					DefaultRegion reg= (DefaultRegion) v.elementAt(j);
					scores.add(reg.getScore());
				}
	//			out.println();
				someIDs[i]= (String) o[i];
			}		
			Distribution dist= new Distribution(scores.toDoubleArray());
			double median= dist.getMedian();
			double threshold= -1d;
			
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
			
			int cntMutExCDom= 0, cntMutExSDom= 0, cntSkipCDom= 0, cntSkipSDom= 0, cntConstCDom= 0, cntConstSDom= 0;
			IntVector lenMutEx= new IntVector(), lenMutExCDom= new IntVector(), lenMutExSDom= new IntVector(),
				lenSkipEx= new IntVector(), lenSkipExCDom= new IntVector(), lenSkipExSDom= new IntVector(),
				lenConstEx= new IntVector(), lenConstExCDom= new IntVector(), lenConstExSDom= new IntVector();
			Vector<Exon> mutexExon= new Vector<Exon>();
			int cntMutexEv= 0;
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					if (!genes[i].isProteinCoding())
						continue;
						
						// cnt events
					ASVariation[] vars= genes[i].getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
					int localCntMutEx= 0;
					for (int j = 0; vars!= null&& j < vars.length; j++) {
						if (vars[j].toString().equals(ASVariation.ID_MUTEX)&& vars[j].isContainedCDS())
							++localCntMutEx;
					}
					cntMutexEv+= localCntMutEx;
					
						// get genomic domain regions
					Transcript[] trpts= genes[i].getTranscripts();
					for (int j = 0; j < trpts.length; j++) {
						Vector v= (Vector) map.remove(trpts[j].getTranscriptID());
						Vector diregV= new Vector();
						for (int k = 0; v!= null&& k < v.size(); k++) {
							DefaultRegion regD= (DefaultRegion) v.elementAt(k);
							if (threshold> 0d&& regD.getScore()> threshold)
								continue;
							Translation tln= trpts[j].getTranslations()[0];
							DirectedRegion reg= new DirectedRegion();
							reg.setStrand(trpts[j].getStrand());
							reg.setChromosome(trpts[j].getChromosome());
							reg.setStart(tln.getGenomicPosition((regD.getStart()- 1)*3));		// TODO check for 0-based
							reg.setEnd(tln.getGenomicPosition((regD.getEnd())*3));		// -1 +1
							reg.setScore(regD.getScore());
							reg.setID(regD.getID());
							diregV.add(reg);
						}
						map.put(trpts[j].getTranscriptID(),diregV);
					}

						// get exons groups
					Exon[] exons= genes[i].getCodingExons(true);	// fully coding
					Vector<Exon> exMutex= new Vector<Exon>(), exSkipped= new Vector<Exon>(), exConstit= new Vector<Exon>();
					int localCntMEexons= 0;
					for (int j = 0; exons!= null&& j < exons.length; j++) {
						if (exons[j].isConstitutive(ASMultiVariation.FILTER_CONTAINED_IN_CDS))
							exConstit.add(exons[j]);
						else if (exons[j].hasVariation(ASVariation.ID_MUTEX, ASMultiVariation.FILTER_CONTAINED_IN_CDS)) {
							exMutex.add(exons[j]);
							++localCntMEexons;
						} else if (exons[j].hasVariation(ASVariation.ID_SKIPPED, ASMultiVariation.FILTER_CONTAINED_IN_CDS))	// else for the not really mutex?
							exSkipped.add(exons[j]);
					}
					if (localCntMutEx* 2!= localCntMEexons) 
						System.currentTimeMillis();
					
						// overlap
					int[] res= countExonDomOvl(exMutex, map, lenMutExCDom, lenMutExSDom, lenMutEx,
							mutexExon);
					cntMutExCDom+= res[0];
					cntMutExSDom+= res[1];
					if (localCntMutEx* 2!= (res[0]+res[1]))
						System.currentTimeMillis();
					res= countExonDomOvl(exSkipped, map, lenSkipExCDom, lenSkipExSDom, lenSkipEx, new Vector<Exon>());
					cntSkipCDom+= res[0];
					cntSkipSDom+= res[1];
					res= countExonDomOvl(exConstit, map, lenConstExCDom, lenConstExSDom, lenConstEx, new Vector<Exon>());
					cntConstCDom+= res[0];
					cntConstSDom+= res[1];
				}
				
				try {
					reader.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				genes= reader.getGenes();
				System.gc();
			}
			
				// output
			System.out.println("Exons constit. "+(cntConstCDom+cntConstSDom)+", contain dom "+cntConstCDom+" ("
					+Formatter.fprint(100d* cntConstCDom/(cntConstCDom+cntConstSDom), 2)+"%) med len ");
			System.out.println("Exons skip "+(cntSkipCDom+cntSkipSDom)+", contain dom "+cntSkipCDom+" ("
					+Formatter.fprint(100d* cntSkipCDom/(cntSkipCDom+cntSkipSDom), 2)+"%)");
			System.out.println("MutEx events "+cntMutexEv);
			System.out.println("Exons mutex "+(cntMutExCDom+cntMutExSDom)+", contain dom "+cntMutExCDom+" ("
					+Formatter.fprint(100d* cntMutExCDom/(cntMutExCDom+cntMutExSDom), 2)+"%)");
			System.out.println("Med lengths:\t(all\tcDom\tsDom)");
			System.out.println("constit.\t"+ 
					new Distribution(lenConstEx.toIntArray()).getMedian()+"\t"+
					new Distribution(lenConstExCDom.toIntArray()).getMedian()+"\t"+
					new Distribution(lenConstExSDom.toIntArray()).getMedian()+"\t");
			System.out.println("skip\t"+ 
					new Distribution(lenSkipEx.toIntArray()).getMedian()+"\t"+
					new Distribution(lenSkipExCDom.toIntArray()).getMedian()+"\t"+
					new Distribution(lenSkipExSDom.toIntArray()).getMedian()+"\t");
			System.out.println("mutex\t"+ 
					new Distribution(lenMutEx.toIntArray()).getMedian()+"\t"+
					new Distribution(lenMutExCDom.toIntArray()).getMedian()+"\t"+
					new Distribution(lenMutExSDom.toIntArray()).getMedian()+"\t");
			System.out.println("\nMutEx exons:");
			for (int i = 0; i < mutexExon.size(); i++) {
				System.out.print(mutexExon.elementAt(i)+"\t"+mutexExon.elementAt(i).toUCSCString()+"   :   \t");
				for (int j = 0; j < mutexExon.elementAt(i).getTranscripts().length; j++) 
					System.out.print(mutexExon.elementAt(i).getTranscripts()[j]+",");
				System.out.println();
			}
		}

	public static strictfp void _05_00_pfamDomains_all_contained() {
		long t0= System.currentTimeMillis();
		String input= "H.sapiens.All.mRNAs.with.CDS.pfam.parsed";
		String speName= "human";
		String annoName= "RefSeq";
		String[] keywords= new String[] {"_xref"}; 
		
		System.out.print("Loading domains..");
		System.out.flush();
		DomainWrapper wrapper= new DomainWrapper(
				"douglas"+File.separator+input);
		try {
			wrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}		
		HashMap map= wrapper.getMap(); 
		Object[] o= map.keySet().toArray();
		String[] someIDs= new String[o.length];
		int cntDomains= 0;
		DoubleVector scores= new DoubleVector();
		for (int i = 0; i < o.length; i++) { 
//			out.print(o[i]+"\t");
			Vector v= (Vector) map.get(o[i]); 
			cntDomains+= v.size(); 
			for (int j = 0; j < v.size(); j++) { 
//				out.print(v.elementAt(j)+",");
				DefaultRegion reg= (DefaultRegion) v.elementAt(j);
				scores.add(reg.getScore());
			}
//			out.println();
			someIDs[i]= (String) o[i];
		}		
		Distribution dist= new Distribution(scores.toDoubleArray());
		double median= dist.getMedian();
		double threshold= -1d;
		PrintStream out= null;
		try {
			out= new PrintStream("douglas"+File.separator+"pfam"+File.separator+
					speName+"_"+annoName+"_"+threshold+"_2contained"+
					"."+input);
		} catch (FileNotFoundException e1) {			
			e1.printStackTrace();
		}
		out.println("found "+someIDs.length+" genes, "+cntDomains+" domains.\n");
		out.println("score med "+median+", min "+dist.getMin()+", max "+dist.getMax());
		
		
			// get basic transcripts
		String absFName= Species.getAnnotation(speName, null, annoName, keywords); 
			//Constants.getLatestUCSCAnnotation(speName, annoName, keywords);
		out.println("Searching for reference transcripts in "+absFName+".");
		GTFChrReader reader= new GTFChrReader(absFName);
		reader.setChromosomeWise(true);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		
		int cntFoundTrpts= 0;
		int cntDomSingle= 0, cntSingle= 0;
		IntVector asGenesV= new IntVector(), asInDomV= new IntVector(), asInNonDomV= new IntVector(),
		asDomGenes= new IntVector(), asNonDomGenes= new IntVector();
		HashMap mapLandscapeDom= new HashMap(), mapLandscapeNonDom= new HashMap(), 
			mapLandscapeAllDomProts= new HashMap(), mapLandscapeAll= new HashMap(),
			mapLandscapeNonDomProts= new HashMap();
		int cntAllGenes= 0, cntNonCoding= 0, cntNonASAllGenes= 0, cntNonASDomGenes= 0, 
				cntNonASDom= 0, cntNonASNonDom= 0, cntNonASNonDomGenes= 0;
		HashMap<DirectedRegion,Vector> domainSpecMap= new HashMap<DirectedRegion,Vector>();
		HashMap<String,Vector<DirectedRegion>> domGroupMap= new HashMap<String,Vector<DirectedRegion>>();
		Vector<ASVariation> mutExDom= new Vector<ASVariation>();
		while (genes!= null) {
			for (int i = 0; i < genes.length; i++) {
				++cntAllGenes;
				if (!genes[i].isProteinCoding())
					++cntNonCoding;
				Transcript[] trpts= genes[i].getTranscripts();
				boolean found= false;
				Vector diregV= new Vector();
				for (int j = 0; j < trpts.length; j++) {
					Vector v= (Vector) map.get(trpts[j].getTranscriptID());
					for (int k = 0; v!= null&& k < v.size(); k++) {
						DefaultRegion regD= (DefaultRegion) v.elementAt(k);
						if (threshold> 0d&& regD.getScore()> threshold)
							continue;
						++cntFoundTrpts;
						found= true;
						Translation tln= trpts[j].getTranslations()[0];
						DirectedRegion reg= new DirectedRegion();
						reg.setStrand(trpts[j].getStrand());
						reg.setChromosome(trpts[j].getChromosome());
						reg.setStart(tln.getGenomicPosition((regD.getStart()- 1)*3));		// TODO check for 0-based
						reg.setEnd(tln.getGenomicPosition((regD.getEnd())*3));		// -1 +1
						reg.setScore(regD.getScore());
						reg.setID(regD.getID());
						diregV.add(reg);
						
						Vector<DirectedRegion> vv= domGroupMap.get(reg.getID());
						if (vv== null)
							vv= new Vector<DirectedRegion>();
						vv.add(reg);
						domGroupMap.put(reg.getID(), vv);
					}
				}
	
					// events
				DirectedRegion[] diregs= null;
				if (found) 
					diregs= (DirectedRegion[]) Arrays.toField(diregV);
				ASVariation[][] vars= Analyzer.getASVariations(genes[i], ASMultiVariation.FILTER_HIERARCHICALLY, ASVariation.TYPE_ALL);
				int cntASEvents= 0, cntASinDom= 0, cntASinNonDom= 0, cntASinDomGenes= 0, cntASinNonDomGenes= 0;
				if (vars== null|| vars.length== 0) {
					++cntNonASAllGenes;
					if (found)
						++cntNonASDomGenes;
					else
						++cntNonASNonDomGenes;
				} 
				boolean overlapASdom= false;
				for (int x = 0; vars!= null&& x < vars.length; x++) {
					for (int xx = 0; xx < vars[x].length; xx++) {
						if (!vars[x][xx].isContainedCDS())
							continue;
						++cntASEvents;
						Vector v= (Vector) mapLandscapeAll.get(vars[x][xx].toString());
						if (v== null)
							v= new Vector();
						v.add(vars[x][xx]);
						mapLandscapeAll.put(vars[x][xx].toString(), v);
						if (!found) {
							v= (Vector) mapLandscapeNonDomProts.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscapeNonDomProts.put(vars[x][xx].toString(), v);
							++cntASinNonDomGenes;
							continue;
						}
						++cntASinDomGenes;
						DirectedRegion asReg= vars[x][xx].getRegion();
						int y;
						boolean counted= false, overlap= false;
						for (y = 0; diregs!= null&& y < diregs.length; y++) 
							//if (asReg.overlaps(diregs[y])) {	// case of AS can be bigger than domain
							if (diregs[y].contains(asReg)) {
							//if (diregs[y].contains(asReg)|| asReg.contains(diregs[y])) {
								++cntASinDom;
								String asID= vars[x][xx].toString();
								if (!counted) {
									v= (Vector) mapLandscapeDom.get(asID);
									if (v== null)
										v= new Vector();
									v.add(vars[x][xx]);
									mapLandscapeDom.put(asID, v);
									if (asID.equals(ASVariation.ID_MUTEX))
										mutExDom.add(vars[x][xx]);
								}
								counted= true;
								
								if ((!overlap)&& asReg.overlaps(diregs[y])) {
									overlap= true;
									overlapASdom= true;
								}
								
								v= (Vector) domainSpecMap.get(diregs[y]);
								if (v== null)
									v= new Vector();
								v.add(vars[x][xx]);
								domainSpecMap.put(diregs[y], v);
								
								//	break;	// fucks up the domainSpecMap
							}
						if (!overlap) {	// hit in non-domain
							++cntASinNonDom;
							v= (Vector) mapLandscapeNonDom.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscapeNonDom.put(vars[x][xx].toString(), v);
						}
						v= (Vector) mapLandscapeAllDomProts.get(vars[x][xx].toString());
						if (v== null)
							v= new Vector();
						v.add(vars[x][xx]);
						mapLandscapeAllDomProts.put(vars[x][xx].toString(), v);
					}
				}
				
				if (vars!= null&& vars.length> 0) {	// events
					asGenesV.add(cntASEvents);
					if (found) {	// domains
						asDomGenes.add(cntASinDomGenes);
						for (int j = 0; j < diregs.length; j++) {
							Vector v= domainSpecMap.get(diregs[j]);
							if (v== null)
								++cntNonASDom;
							else
								asInDomV.add(v.size());
						}
						if (cntASinNonDom== 0)
							++cntNonASNonDom;
						else
							asInNonDomV.add(cntASinNonDom); 	// can be 0
					} else
						asNonDomGenes.add(cntASinNonDomGenes);
				} else if (found) {	// no events
					cntNonASDom+= diregs.length;
					++cntNonASNonDom;	// ??
				}

			}
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
		}
		
			// output
		out.println(cntAllGenes+" genes, "+cntNonCoding+" non-coding genes.");
		out.println(
			cntNonASAllGenes+ " non AS genes ("+Formatter.fprint(100d* cntNonASAllGenes/cntAllGenes, 2)+"%)\n"+
			(cntNonASDomGenes+asDomGenes.length)+ " genes with annotated domains, "+
				asDomGenes.length+" alt. spliced ("+ Formatter.fprint(100d* asDomGenes.length/ (cntNonASDomGenes+ asDomGenes.length), 2)+"%)\n"+
			(cntNonASNonDomGenes+asNonDomGenes.length)+ " genes w/no annotated domains, "+
				asNonDomGenes.length+" alt. spliced ("+ Formatter.fprint(100d* asNonDomGenes.length/ (cntNonASNonDomGenes+ asNonDomGenes.length), 2)+"%)\n"+
			(cntNonASDom+ asInDomV.length)+" domains, "+
				asInDomV.length+" alt. spliced ("+ Formatter.fprint(100d* asInDomV.length/(cntNonASDom+ asInDomV.length), 2)+ "%)\n"+
			(cntNonASNonDom+ asInNonDomV.length)+" genes w domains, "+
				asInNonDomV.length+" alt. spliced outside of domains ("+ Formatter.fprint(100d* asInNonDomV .length/(cntNonASNonDom+ asInNonDomV.length), 2)+ "%)\n"
		);
		
		
		out.println("\n==================");
		out.println("All Genes:");
		dist= new Distribution(asGenesV.toIntArray());
		out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
		Analyzer.outputLandscape(mapLandscapeAll, out);
		
		out.println("\n==================");
		out.println("All Domain Genes:");
		dist= new Distribution(asDomGenes.toIntArray());
		out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
		Analyzer.outputLandscape(mapLandscapeAllDomProts, out);
		
		out.println("\n==================");
		out.println("All Non-Domain Genes:");
		dist= new Distribution(asNonDomGenes.toIntArray());
		out.println(dist.getSum()+" events, med "+dist.getMedian()+" per gene.");
		Analyzer.outputLandscape(mapLandscapeNonDomProts, out);

		out.println("\n==================");
		out.println("in Domains:");
		dist= new Distribution(asInDomV.toIntArray());
		out.println(dist.getSum()+" events, med "+dist.getMedian()+" per domain.");
		Analyzer.outputLandscape(mapLandscapeDom, out);
		for (int i = 0; i < mutExDom.size(); i++) 
			System.out.println(ASVariation.ID_MUTEX+"\t"+mutExDom.elementAt(i).toCoordinates());
		
		
		out.println("\n==================");
		out.println("in Non-Domains:");
		dist= new Distribution(asInNonDomV.toIntArray());
		out.println(dist.getSum()+" events, med "+dist.getMedian()+" per non-domain.");
		Analyzer.outputLandscape(mapLandscapeNonDom, out);


		HashMap<String,Vector> domGroupEventMap= new HashMap<String,Vector>(domGroupMap.size());
		Iterator<String> iter= domGroupMap.keySet().iterator();
		while(iter.hasNext()) {
			String oID= (String) iter.next();
			Vector<ASVariation> vv= domGroupEventMap.get(oID);
			if (vv== null)
				vv= new Vector<ASVariation>();
			Vector<DirectedRegion> v= (Vector<DirectedRegion>) domGroupMap.get(oID);
			for (int i = 0; i < v.size(); i++) { 
				Vector<ASVariation> inV= domainSpecMap.get(v.elementAt(i));	// directedRegx V(ASevents)
				if (inV== null)
					continue;
				for (int j = 0; j < inV.size(); j++) 
					vv.add(inV.elementAt(j));
				if (vv.size()> 0)
					domGroupEventMap.put(oID, vv);
			}
			
		}
		out.println("\n\nDomain Group Events");
		DualHashBidiMap mapD= new DualHashBidiMap(domGroupEventMap);
		Object[] values= domGroupEventMap.values().toArray();
		Comparator compi= new VectorSizeComparator();
		java.util.Arrays.sort(values, compi);
		BidiMap revMap= mapD.inverseBidiMap();
		for (int i = values.length- 1; i >= 0; --i) {
			Vector<ASVariation> varV= (Vector<ASVariation>) values[i];
			out.print(varV.size()+ "\t"+ revMap.get(values[i])+"\t|");
			for (int j = 0; j < varV.size(); j++) 
				out.print(varV.elementAt(j)+"|");
			out.println();
		}
		
			
		out.println("\n\nDomain Groups");
		mapD= new DualHashBidiMap(domGroupMap);
		values= domGroupMap.values().toArray();
		java.util.Arrays.sort(values, compi);
		revMap= mapD.inverseBidiMap();
		for (int i = values.length- 1; i >= 0; --i) {
			Vector<ASVariation> varV= (Vector<ASVariation>) values[i];
			out.print(varV.size()+ "\t"+ revMap.get(values[i])+"\t|");
//			for (int j = 0; j < varV.size(); j++) 
//				out.print(varV.elementAt(j)+"|");
			out.println();
		}
		
		out.println("\n\nSingle Domains");
		IntVector cntEvPerDomain= new IntVector();
		Iterator<Vector> iterVec= domainSpecMap.values().iterator();
		while (iter.hasNext()) 
			cntEvPerDomain.add(iterVec.next().size());		
		dist= new Distribution(cntEvPerDomain.toIntArray());
		out.println("events per domain: med "+dist.getMedian()+", min "+dist.getMin()+", max "+dist.getMax());
		Analyzer.outputLandscape(domainSpecMap, out);
		
		out.println("\n\n[took "+(System.currentTimeMillis()- t0)/1000+" sec.]");
		out.flush();
		out.close();
	}
	
	public static void _05_01_entireDistribution() {
		String fName= Constants.getLatestUCSCAnnotation("mouse", "Ensembl", null);
		//fName= "D:\\workspace\\G-Phase\\annotation\\ensembl\\fruitfly_EnsEMBL43_BDGP43_dlMart.gff";
		
		System.out.println("Retrieving events for "+fName);
		GTFChrReader reader= new GTFChrReader(fName);
		reader.setChromosomeWise(true);
		reader.setSilent(false);
		Method m = null, m2= null;
		try {
			m= ASVariation.class.getMethod("is_contained_in_CDS", null);
			m2= ASVariation.class.getMethod("has_onlyGTAG_Intron", null);
		} catch (Exception e) {
			e.printStackTrace();
		}
		Method[] filters= new Method[] {m};	// m2
		
		ASVariation[][] vars= Analyzer.getASVariations(reader, ASMultiVariation.FILTER_HIERARCHICALLY, filters);	//		
		ASAnalyzer.outputVariations(vars, false, false, System.out);
		
//		for (int i = 0; i < vars.length; i++) {
//			if (!vars[i][0].toString().equals("1^ , 2^"))
//				continue;
//			for (int j = 0; j < vars[i].length; j++) {
//				String[] codons= vars[i][j].getVariableCodons();
//				System.currentTimeMillis();
//			}
//		}
	}
	
}
