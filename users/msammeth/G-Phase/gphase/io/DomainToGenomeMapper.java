package gphase.io;

import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Arrays;
import gphase.tools.File;

import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

public class DomainToGenomeMapper {

	public final static String DOMAIN_ANNOTATION_ATTRIBUTE_ID= "domain_annotation";
	public final static String DOMAIN_ANNOTATION_DOM_ANNOTATED= "annotated";	// everything fine
	public final static String DOMAIN_ANNOTATION_DOM_OMITTED= "removed";	// discarded due to no CDS
	public final static String DOMAIN_ANNOTATION_DOM_NA= "NA";	// no hits above threshold
	public final static String DOMAIN_ANNOTATION_DOM_NONE= "none";	// not run

	public static final String GTF_DOMAIN_FEATURE_TAG= "domain";
	public static final String GTF_DOMAIN_ID_TAG= GTF_DOMAIN_FEATURE_TAG+ "_id";
	static boolean outConflicts= false;
	static boolean domainsOnly= false;

	
	static final String errorMsg= "usage: DomainToGenomeMapper [options] <domainFile> <genomeAnnotationFile>\n\n"+
		"where\n"+
		"<domainFile>\ta file containing the domains (Douglas format)\n"+
		"<genomeAnnotationFile>\ta GTF file with the reference annotation\n\n"+
		"and [options] may be:\n"+
		"-writeAll\twrite merge the domain coordinates with the input GTF file (default is writing a GTF with the domain coordinates only)"+
		"-conflicts\toutput a GTF file with the coordinates of overlapping domains"+
		"\n\nmicha, may 07";
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length< 2) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		File domFile= new File(args[args.length- 2]);
		if (!domFile.exists())
			System.err.println("Domain file not valid.");
		File annFile= new File(args[args.length- 1]);
		if (!annFile.exists())
			System.err.println("Annotation file not valid.");
		for (int i = 0; i < args.length- 2; i++) {
			if (args[i].equalsIgnoreCase("-writeAll")) {
				domainsOnly= false;
				continue;
			}
			if (args[i].equalsIgnoreCase("-conflicts")) {
				outConflicts= true;
				continue;
			}
			
		}
		String outFileName= annFile.getFileNameWithoutExtension();
		if (outConflicts) {
			outFileName+="_conflicts";
		} else {
			if (domainsOnly)
				outFileName+="_domains";
			else
				outFileName+="_addDomains";
		}
		outFileName= annFile.getPathOnly()+ File.separator+ outFileName+"."+ annFile.getExtension();
		File outFile= new File(outFileName);
		File.checkForOverwrite(System.out, outFile);
		
	
		System.out.println("Reading domains..");
		DomainWrapper wrapper= new DomainWrapper(domFile.getAbsolutePath());
		try {
			wrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		if (outConflicts) {
			System.out.println("Mapping conflicting domains..");
			mapToGenomicCoordinates(wrapper.getMapConflict(), annFile, outFile, true);
		} else {
			System.out.println("Mapping (non-conflicting) domains..");
			// tID x domainRegions (AA coordinates)
			mapToGenomicCoordinates(wrapper.getMap(), annFile, outFile, domainsOnly);
			// tID x domainsRegions (genomic coordinates)
		}
	}
	
	public static HashMap mapToGenomicCoordinates(HashMap mapAA, File annFile, File outFile, boolean domainsOnly) {
		GTFChrReader reader= new GTFChrReader(annFile.getAbsolutePath());
		GTFChrReader writer= new GTFChrReader(outFile.getAbsolutePath());
		reader.setChromosomeWise(true);
		reader.setReadGene(true);
		reader.setReadGTF(true);
		reader.setReadAllLines();
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		Vector outOfRangeV= new Vector();
		int cntDomains= 0;
		Vector nonCodTrptV= new Vector();
		GTFObject[] obs= reader.getGtfObj();
		HashMap mapGTF= new HashMap();
		for (int i = 0; i < obs.length; i++) {
			Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
			if (v== null)
				v= new Vector();
			v.add(obs[i]);
			mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
		}
		while (genes!= null) {
			for (int i = 0; i < genes.length; i++) {
					// get genomic domain regions
				Transcript[] trpts= genes[i].getTranscripts();
				for (int j = 0; j < trpts.length; j++) {
					GTFObject[] gtfRest= new GTFObject[0];
					if (!domainsOnly)
						gtfRest= (GTFObject[]) Arrays.toField(mapGTF.get(trpts[j].getTranscriptID()));
					
					Vector aaRegV= (Vector) mapAA.remove(trpts[j].getTranscriptID());
					if (!trpts[j].isCoding()&& aaRegV!= null&& aaRegV.size()> 0) {
						nonCodTrptV.add(trpts[j].getTranscriptID());
						for (int k = 0; k < gtfRest.length; k++) {
							if (DOMAIN_ANNOTATION_DOM_ANNOTATED.equals(((DefaultRegion) aaRegV.elementAt(0)).getAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID)))
								gtfRest[k].addAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID, DOMAIN_ANNOTATION_DOM_OMITTED);
							//else
								//gtfRest[k].addAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID, DOMAIN_ANNOTATION_DOM_NONE);
						}
						//continue;
					}

					GTFObject[] gtfDoms= new GTFObject[0];
					Vector diregV= new Vector();
					for (int k = 0; aaRegV!= null&& k < aaRegV.size(); k++) {
						DefaultRegion regD= (DefaultRegion) aaRegV.elementAt(k);
						if (regD.getAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID).equals(DOMAIN_ANNOTATION_DOM_NONE)) {
							for (int m = 0; m < gtfRest.length; m++) 
								gtfRest[m].addAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID, DOMAIN_ANNOTATION_DOM_NONE);
							//continue;
							gtfDoms= gtfRest;
						} else if (trpts[j].getTranslations()!= null&& trpts[j].getTranslations().length> 0) {
							Translation tln= trpts[j].getTranslations()[0];
							DirectedRegion reg= new DirectedRegion();
							reg.setStrand(trpts[j].getStrand());
							reg.setChromosome(trpts[j].getChromosome());
							reg.setStart(tln.getGenomicPosition((regD.getStart()- 1)*3));		// TODO check for 0-based
							reg.setEnd(tln.getGenomicPosition((regD.getEnd())*3));		// -1 +1
							reg.setScore(regD.getScore());
							reg.setID(regD.getID());
							if (reg.get5PrimeEdge()< tln.get5PrimeEdge()|| reg.get3PrimeEdge()> tln.get3PrimeEdge()) {
								regD.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, trpts[j].getTranscriptID());
								outOfRangeV.add(regD);						
							} else
								diregV.add(reg);
						}
					}
					DirectedRegion[] genDomains= null;
					if (aaRegV== null) {
						for (int m = 0; m < gtfRest.length; m++) 
							gtfRest[m].addAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID, DOMAIN_ANNOTATION_DOM_NA);
						//continue;
						gtfDoms= gtfRest;
					} else {
						genDomains= (DirectedRegion[]) Arrays.toField(diregV);
						if (genDomains!= null)
							cntDomains+= genDomains.length;
						mapAA.put(trpts[j], genDomains);
					}
					
					if (outFile!= null) {	//always add now && genDomains!= null&& genDomains.length> 0) {
						// get domain objects
						if (genDomains!= null) {
							Exon[] exons= trpts[j].getExons();
							DirectedRegion[] exDomRegs= DirectedRegion.intersect(genDomains, exons);	// ensure is sorted!
							gtfDoms= new GTFObject[exDomRegs.length];
							for (int k = 0; k < gtfDoms.length; k++) {
								gtfDoms[k]= GTFObject.createGTFObjects(exDomRegs[k], trpts[j])[0];
								//gtfDoms[k].addAttribute(GTFObject.TRANSCRIPT_ID_TAG, trpts[j].getTranscriptID());
								gtfDoms[k].setSeqname(trpts[j].getChromosome());
								gtfDoms[k].setFeature(GTF_DOMAIN_FEATURE_TAG);
								//gtfDoms[k].setSource(trpts[j].getSpecies().getAnnotationVersion());
								HashMap doubleNames= new HashMap();
								for (int m = 0; m < genDomains.length; m++) {
									Integer cnt= (Integer) doubleNames.get(genDomains[m].getID());
									if (cnt== null)
										cnt= new Integer(1);
									else
										cnt= new Integer(cnt.intValue()+ 1);
									doubleNames.put(genDomains[m].getID(), cnt);
									if (genDomains[m].overlaps(exDomRegs[k])) {
										gtfDoms[k].setScore(Double.toString(genDomains[m].getScore()));
										gtfDoms[k].addAttribute(GTF_DOMAIN_ID_TAG, genDomains[m].getID()+"-"+cnt);
										break;
									}
								}
							}
							if (gtfDoms.length> 0)
								for (int m = 0; m < gtfRest.length; m++) 
									gtfRest[m].addAttribute(DOMAIN_ANNOTATION_ATTRIBUTE_ID, DOMAIN_ANNOTATION_DOM_ANNOTATED);
						
							if (!domainsOnly) {
								GTFObject[] gtfAll= new GTFObject[gtfDoms.length+ gtfRest.length];
								int posRest= 0, posDom= 0, pos= 0;	// merge
								while (true) {
									if (posRest== gtfRest.length|| posDom== gtfDoms.length)
										break;
									if ((trpts[j].isForward()&& gtfDoms[posDom].getStart()< gtfRest[posRest].getStart())||
											((!trpts[j].isForward())&& gtfDoms[posDom].getEnd()> gtfRest[posRest].getEnd()))
										gtfAll[pos++]= gtfDoms[posDom++];
									else
										gtfAll[pos++]= gtfRest[posRest++];
								}
								for (int k = posRest; k < gtfRest.length; k++) 
									gtfAll[pos++]= gtfRest[posRest++];
								for (int k = posDom; k < gtfDoms.length; k++) 
									gtfAll[pos++]= gtfDoms[posDom++];
								
								gtfDoms= gtfAll;
							}
						}
						
						writer.setGtfObj(gtfDoms);
						try {
							writer.write(true);
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
				}

					
			}
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
			obs= reader.getGtfObj();
			mapGTF= new HashMap();
			for (int i = 0; obs!= null&& i < obs.length; i++) {
				Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
				if (v== null)
					v= new Vector();
				v.add(obs[i]);
				mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
			}
			System.gc();
			Thread.yield();
		}
		
		// eliminate AA regions without tID found
		Object[] keys= mapAA.keySet().toArray();
		Vector notFound= new Vector();
		int cntTrpts= 0;
		for (int i = 0; i < keys.length; i++) {
			try {
				Object o= (DirectedRegion[]) mapAA.get(keys[i]);
				++cntTrpts;
			} catch (ClassCastException e) {
				notFound.add(keys[i]);
				mapAA.remove(keys[i]);
			}
		}
		System.out.println("Found "+cntTrpts+" transcripts with annotated domains (in total "+cntDomains+" domains).");
		if (notFound.size()> 0) {
			System.out.println("WARNING: didnt find "+notFound.size()+" transcript IDs:");
			for (int i = 0; i < notFound.size(); i++) 
				System.out.print(notFound.elementAt(i)+" ");
			System.out.println();
		}
		if (nonCodTrptV.size()> 0) {
			System.out.println("WARNING: "+nonCodTrptV.size()+" transcripts with domains but without CDS.");
			for (int i = 0; i < nonCodTrptV.size(); i++) 
				System.out.print(nonCodTrptV.elementAt(i)+" ");
			System.out.println();
		}
		if (outOfRangeV.size()> 0) {
			String fName= "domains_out_of_range.txt";
			System.out.println("WARNING: "+outOfRangeV.size()+" domains out of the translated area, wrote details to "+fName);
			try {
				PrintStream p= new PrintStream(fName);
				for (int i = 0; i < outOfRangeV.size(); i++) {
					DefaultRegion reg= (DefaultRegion) outOfRangeV.elementAt(i);
				p.println(reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG)+"\t"+
							reg.getID()+"\t"+
							reg.getStart()+"\t"+
							reg.getEnd()+"\t"+
							reg.getScore()
					);
				}
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	
		return mapAA;
	} 
		
}
