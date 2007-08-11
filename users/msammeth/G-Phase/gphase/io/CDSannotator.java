package gphase.io;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.File;

// H:annotations\human_hg17_mRNA_cds.txt
public class CDSannotator {

	public final static String CDS_ANNOTATION_ATTRIBUTE_ID= "cds_annotation";
	public final static String CDS_ANNOTATION_CDS_ANNOTATED= "annotated";
	public final static String CDS_ANNOTATION_CDS_OMITTED= "removed";
	public final static String CDS_ANNOTATION_CDS_NA= "NA";
	public final static String CDS_ANNOTATION_CDS_NONE= "none";
	
	final static String UCSC_DB_MRNA_QSTART= "all_mrna.qStart";
	final static String UCSC_DB_MRNA_QBLOCKSIZES= "all_mrna.blockSizes"; 
	final static String UCSC_DB_MRNA_QBLOCKSTARTS= "all_mrna.qStarts";   
	final static String UCSC_DB_ACCESSION= "gbCdnaInfo.acc";	//"all_mrna.qName";//     
	final static String UCSC_DB_CDS= "cds.name";
	final static String UCSC_DB_STRAND= "all_mrna.strand";
	final static String UCSC_DB_QSIZE= "all_mrna.qSize";
	
	
	static int cntBadAlignedCDS= 0;
	
	static final String errorMsg = "usage: CDSannotator [-out <outFileName> -CDSonly] [CDS file] <GTF file> \n\n" + 
	"where\n" + "<GTF file>\ta GTF file containing annotated transcripts\n" + 
	"<CDS file>\ta tab-delimited file with transcript_id \t CDS coordinates on mRNA (Genbank format)\n\n"+
	"and optional\n"+
	"-out\t specifies the output file name (default is <input file name>_CDS.gtf)\n"+
	"-CDSonly\t writes only the gtf lines of the inserted cds objects.\n"+
	"-chkIntrons\t checks and replaces introns of size 0"
	+"\n\nmicha, may 07";

	static boolean writeAll= true, chkIntrons= false;
	static int cntSuccess= 0, cntMisaligned= 0, cntUnspliced= 0, cntParseError= 0;
	static int corrIntrons= 0;
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		if (args.length< 2) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		File gtfFile= new File(args[args.length- 1]);
		if (!gtfFile.exists()) {
			System.err.println("GTF file not valid.");
			System.exit(0);
		}
		File annFile= new File(args[args.length- 2]);
		if (!annFile.exists()) {
			System.err.println("CDS file not valid.");
			annFile= null;
		}
		
		String outFileName= gtfFile.getPathOnly()+ File.separator+
			gtfFile.getFileNameWithoutExtension()+ "_CDS.gtf";
		
			// options
		for (int i = 0; i < args.length- 2; i++) {
			if (args[i].equalsIgnoreCase("-out")|| args[i].equalsIgnoreCase("-outfile")) {
				outFileName= args[i+1];
				++i;
				continue;
			}
			if (args[i].equalsIgnoreCase("-CDSonly")) {
				writeAll= false;
				continue;
			}
			if (args[i].equalsIgnoreCase("-chkIntrons")) {
				chkIntrons= true;
				continue;
			}

		}
		File outFile= new File(outFileName);
		if (!File.checkForOverwrite(System.out, outFile))
			System.exit(0);

		insert(gtfFile, annFile, outFile);
		System.out.println("[took "+(System.currentTimeMillis()- t0)/1000+" sec]");
	}
	
	public static void insert(File gtfFile, File annFile, File outFile) {
	
		try {
			HashMap cdsMap= null;
			if (annFile!= null) {
				System.out.println("Reading CDS annotation "+annFile.getFileNameOnly()+"..");
				cdsMap= getCDSMap(annFile);
			}
			System.out.println("successfully mapped: "+cntSuccess);
			System.out.println("misaligned mRNAs: "+cntMisaligned);
			System.out.println("unspliced CDSs: "+cntUnspliced);
			System.out.println("parsing errors: "+cntParseError);
			System.out.println("badly aligned CDSs: "+cntBadAlignedCDS);
			
			System.out.println("Mapping to GTF file "+gtfFile.getFileNameOnly()+"..");
			GTFChrReader reader= new GTFChrReader(gtfFile.getAbsolutePath());
			reader.setChromosomeWise(true);
			if (writeAll)
				reader.setReadGTF(true);
			reader.setReadGene(true);
			reader.setReadAllLines();

			reader.read();
			GTFObject[] obs= reader.getGtfObj();
			HashMap mapGTF= new HashMap();
			for (int i = 0; i < obs.length; i++) {
				Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
				if (v== null)
					v= new Vector();
				v.add(obs[i]);
				mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
			}
			Gene[] genes= reader.getGenes();
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					for (int j = 0; j < genes[i].getTranscriptCount(); j++) {
						
						Transcript trpt= genes[i].getTranscripts()[j];
						GTFObject[] gtfRest= null;
						if (writeAll)
							gtfRest= (GTFObject[]) Arrays.toField(mapGTF.get(trpt.getTranscriptID()));
						
						if (chkIntrons) {							
							Exon[] ex= trpt.getExons();
							java.util.Arrays.sort(ex, new AbstractRegion.PositionComparator());
							Vector exonV= new Vector();
							for (int k = 1; k < ex.length; k++) {
								if (ex[k-1].get3PrimeEdge()+1>= ex[k].get5PrimeEdge()) {	// edge-to-edge
									Exon lastExon= ex[k-1];
									if (exonV.size()> 0&& lastExon.get3PrimeEdge()== ((Exon) exonV.elementAt(exonV.size()- 1)).get3PrimeEdge())
										lastExon= (Exon) exonV.elementAt(exonV.size()-1);	// multi-chains
									Exon e= new Exon(trpt, null, lastExon.get5PrimeEdge(), ex[k].get3PrimeEdge());
									e.setChromosome(trpt.getChromosome());
									exonV.remove(lastExon);
									exonV.add(e);
									//System.out.println("killed Intron "+trpt.getTranscriptID()+" "+ex[k-1].get3PrimeEdge()+" -> "+ex[k].get5PrimeEdge());
									++corrIntrons;
								} else if (exonV.size()== 0|| (((Exon) exonV.elementAt(exonV.size()-1)).get3PrimeEdge()!= ex[k-1].get3PrimeEdge()))	// already merged
									exonV.add(ex[k-1]);
							}
							if (exonV.size()> 0&& ((Exon) exonV.elementAt(exonV.size()- 1)).get3PrimeEdge()!= ex[ex.length-1].get3PrimeEdge())
								exonV.add(ex[ex.length- 1]);
							if (exonV.size()>0&& exonV.size()!= ex.length) {
								trpt.setExons((Exon[]) Arrays.toField(exonV));
								if (writeAll) {
									Vector restGTFv= new Vector();
									for (int k = 0; k < gtfRest.length; k++) 
										if (!gtfRest[k].getFeature().equals(GTFObject.EXON_FEATURE_TAG))
											restGTFv.add(gtfRest[k]);
									for (int k = 0; k < exonV.size(); k++) 
										restGTFv.add(GTFObject.createGTFObject((Exon) exonV.elementAt(k), trpt));
									gtfRest= (GTFObject[]) Arrays.toField(restGTFv);
								}

							}
						}
						
						GTFObject[] gtfs= new GTFObject[0];
						if (annFile!= null) {
							DefaultRegion reg= (DefaultRegion) cdsMap.get(trpt.getTranscriptID());
							if (reg!= null&& CDS_ANNOTATION_CDS_ANNOTATED.equals(reg.getAttribute(CDS_ANNOTATION_ATTRIBUTE_ID))) {	// get gtfs for cds | exons
								DirectedRegion regD= new DirectedRegion(
										trpt.getGenomicPosition(reg.getStart()- 1), 
										trpt.getGenomicPosition(reg.getEnd()- 1), 
										trpt.getStrand());
								regD.setChromosome(trpt.getChromosome());
								DirectedRegion[] regs= DirectedRegion.intersect(
										new DirectedRegion[] {regD}, 
										trpt.getExons());
								Vector obsV= new Vector();
								for (int k = 0; regs!= null&& k < regs.length; k++) {
									GTFObject o= GTFObject.createGFFObject(regs[k]);
									o.setFeature(GTFObject.CDS_FEATURE_TAG);
									o.addAttribute(GTFObject.GENE_ID_TAG, trpt.getGene().getGeneID());
									o.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, trpt.getTranscriptID());	
									o.setSource(trpt.getSource());
									obsV.add(o);
								}
								for (int m = 0; m < gtfRest.length; m++)
									if (gtfRest[m].getFeature().equals(GTFObject.EXON_FEATURE_TAG))
										gtfRest[m].addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_ANNOTATED);
								if (obsV.size()> 0)
									gtfs= (GTFObject[]) Arrays.toField(obsV); 
							} else if (gtfRest!= null) {
								for (int m = 0; m < gtfRest.length; m++) {
									if (!gtfRest[m].getFeature().equals(GTFObject.EXON_FEATURE_TAG))
										continue;
									if (reg== null)
										gtfRest[m].addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_NONE);
									else
										gtfRest[m].addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, (String) reg.getAttribute(CDS_ANNOTATION_ATTRIBUTE_ID));
								}
							}
						}
						
						
						if (writeAll) {	// merge
							GTFObject[] gtfAll= new GTFObject[gtfs.length+ gtfRest.length];
							//DefaultRegion reg= (DefaultRegion) cdsMap.get(trpt.getTranscriptID());
							
							Comparator compi= new GTFObject.PositionComparator();
							java.util.Arrays.sort(gtfRest, compi);
							java.util.Arrays.sort(gtfs, compi);
							int posRest= 0, posCDS= 0, pos= 0;	// merge
							while (true) {
								if (posRest== gtfRest.length|| posCDS== gtfs.length)
									break;
								if (compi.compare(gtfs[posCDS], gtfRest[posRest])< 0) {	// cds after exon
									gtfs[posCDS].setSource(gtfRest[0].getSource());
									gtfAll[pos++]= gtfs[posCDS++];
								} else
									gtfAll[pos++]= gtfRest[posRest++];
							}
							for (int k = posRest; k < gtfRest.length; k++) 
								gtfAll[pos++]= gtfRest[posRest++];
							for (int k = posCDS; k < gtfs.length; k++) {
								gtfs[posCDS].setSource(gtfRest[0].getSource());
								gtfAll[pos++]= gtfs[posCDS++];
							}
							
							gtfs= gtfAll;
						}
						
						
						if (gtfs!= null&& gtfs.length> 0) {	// write
							GTFChrReader writer= new GTFChrReader(outFile.getAbsolutePath());
							writer.setGtfObj(gtfs);
							writer.write(true);
						}
					}
				}
				
					// read
				reader.read();
				obs= reader.getGtfObj();
				mapGTF= new HashMap();
				for (int i = 0; obs!= null&& i < obs.length; i++) {
					Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
					if (v== null)
						v= new Vector();
					v.add(obs[i]);
					mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
				}
				genes= reader.getGenes();

			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println("Corrected 0-length introns: "+corrIntrons);

	}

	private static int getFieldNumber(String[] fields, String id) {
		if (fields== null)
			return -1;
		for (int i = 0; i < fields.length; i++) {
			if (fields[i].contains(id))
				return i;
		}
		return -1;
	}
	
	
	
	public static HashMap getCDSMap(File annFile) {
		
		TabDelimitedFormatWrapper tabReader= new TabDelimitedFormatWrapper(annFile.getAbsolutePath());
		try {
			tabReader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		String[][] tab= tabReader.getTable();
		HashMap map= new HashMap();
		// #hg17.all_mrna.qStart   hg17.all_mrna.blockSizes        hg17.all_mrna.qStarts   hg17.gbCdnaInfo.acc     hg17.cds.name
		int field_blockSizes= getFieldNumber(tab[0], UCSC_DB_MRNA_QBLOCKSIZES);
		int field_blockStarts= getFieldNumber(tab[0], UCSC_DB_MRNA_QBLOCKSTARTS);
		int field_cds= getFieldNumber(tab[0], UCSC_DB_CDS);
		int field_acc= getFieldNumber(tab[0], UCSC_DB_ACCESSION);
		int field_strand= getFieldNumber(tab[0], UCSC_DB_STRAND);
		int field_qsize= getFieldNumber(tab[0], UCSC_DB_QSIZE);
		for (int i = 1; tab!= null&& i < tab.length; i++) {
			DefaultRegion reg= parseGenbankCDSpos(tab[i][field_cds], tab[i][field_blockStarts], tab[i][field_blockSizes],
					tab[i][field_strand], tab[i][field_qsize]);
			map.put(tab[i][field_acc],reg);
		}
		
		return map;
	}
	
	
	private static int getGapsUntil(int pos, int[] starts, int[] sizes) {
		
		int gaps= 0;
		int lastEnd= 0;
		for (int i = 0; i < sizes.length; i++) {
			if (starts[i]> pos) {
				if (i== 0)
					gaps= Math.min(starts[i],pos);
//				if (pos> lastEnd)
//					System.out.println("WARNING: pos in not-aligned area.");
				break;	// assuming that pos is in an aligned area
			}
			if (lastEnd!= starts[i])
				gaps+= starts[i]- lastEnd;
			
			lastEnd= starts[i]+ sizes[i];
		}
		
		return gaps;
	}
	
	
	public static DefaultRegion parseGenbankCDSpos(String label, String blockStarts, String blockSizes, String strand, String qSize) {
		
		DefaultRegion dummy= new DefaultRegion();		
		
		if (label.equals("n/a")) {
			dummy.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_NA);
			return dummy;
		}
		
//		if (label.contains("(")) {	// join(), complement(), unclosed reading frames..
//			System.out.println("WARNING: not understood "+label);
//			return null;
//		}

			// remove open ends
		StringBuffer sb= new StringBuffer(label);
		for (int i = 0; i < sb.length(); i++) {
			if (sb.charAt(i)== '<'|| sb.charAt(i)== '>')
				sb.deleteCharAt(i--);
		}
		label= sb.toString();
		
		String compl= "complement(";
		if (label.startsWith(compl)) {
			label= label.substring(compl.length(), label.length()- 1);
		}
		

		String[] s= blockStarts.split(",");
		int[] starts= new int[s.length];
		for (int i = 0; i < starts.length; i++) 
			starts[i]= Integer.parseInt(s[i]);
		s= blockSizes.split(",");
		int[] sizes= new int[s.length];
		for (int i = 0; i < sizes.length; i++) 
			sizes[i]= Integer.parseInt(s[i]);
		
			// reverse for the reverse strand
		int st= GTFObject.parseStrand(strand);
		if (st< 0) {
			int qsiz= Integer.parseInt(qSize);
			
			int[] newStarts= new int[starts.length];
			int[] newSizes= new int[sizes.length];
			for (int i = 0; i < newStarts.length; i++) {
				newStarts[i]= qsiz- (starts[starts.length- 1- i]+ sizes[sizes.length- 1- i]);
				newSizes[i]= sizes[sizes.length- 1- i];
			}
			sizes= newSizes;
			starts= newStarts;
		}
		
		if (starts.length!= sizes.length)
			System.out.println("WARNING: block starts/sizes do not match "+ blockStarts+"  "+blockSizes);

		
		
		Matcher matty= Pattern.compile("^(\\d+)[.]{2}(\\d+)$").matcher(label);
		
		String jID= "join(";
		if (matty.matches()) {
			DefaultRegion reg= new DefaultRegion();
			int first= Integer.parseInt(matty.group(1));
			int gaps= getGapsUntil(first, starts, sizes);
			reg.setStart(first- gaps);
			int last= Integer.parseInt(matty.group(2));
			int gap2= getGapsUntil(last, starts, sizes);
			if (gaps!= gap2) {
				dummy.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_OMITTED);
				++cntBadAlignedCDS;
				return dummy;
			}
			reg.setEnd(last- gap2);
			reg.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_ANNOTATED);
			++cntSuccess;
			return reg;
		} else if (label.startsWith(jID)){ 	// join statement
			label= label.substring(jID.length(), label.length()- 1);
			matty= Pattern.compile("^(\\d+)[.]{2}(\\d+).*").matcher(label);
			int prevEnd= -1, begin= -1, last= -1;
			int spliceDiff= 0;
			while (matty.matches()) {
				int start= Integer.parseInt(matty.group(1));
				int end= Integer.parseInt(matty.group(2));
				if (begin< 0) 
					begin= start;
				last= end;
				if (prevEnd> 0&& start> prevEnd+ 1) {	// unspliced, not fully spliced
					spliceDiff= start- prevEnd- 1;
					break;
				}
				prevEnd= end;
				if (matty.end(2)+1> label.length())
					break;
				matty.region(matty.end(2)+1, matty.regionEnd());
			}
			if (begin>= 0&& last>= begin) {
				if (spliceDiff!= 0) {
					if (spliceDiff>= 33) {
						System.out.println("unspliced CDS\t"+label);
						++cntUnspliced;
					} else { 
						System.out.println("misaligned mRNA\t"+label);
						++cntMisaligned;
					}
					dummy.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_OMITTED);
					return null;
				} else {
					DefaultRegion reg= new DefaultRegion();
					int gaps= getGapsUntil(begin, starts, sizes);
					reg.setStart(begin- gaps);
					int gap2= getGapsUntil(last, starts, sizes);
					if (gaps!= gap2) {
						dummy.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_OMITTED);
						++cntBadAlignedCDS;
						return dummy;
					}
					reg.setEnd(last- gap2);
					++cntSuccess;
					return reg;
				}
			}
		} 
		
		System.out.println("pbs with "+label);
		++cntParseError;
		dummy.addAttribute(CDS_ANNOTATION_ATTRIBUTE_ID, CDS_ANNOTATION_CDS_OMITTED);
		return dummy;
	}

}
