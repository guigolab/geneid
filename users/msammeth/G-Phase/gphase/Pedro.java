package gphase;

import gphase.io.TabPedroReader;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASEvent;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.model.VariantGroup_SpliceChain;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.Formatter;
import gphase.tools.IntVector;
import gphase.tools.Sequence;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Date;
import java.util.HashMap;
import java.util.Observable;
import java.util.Vector;

import org.apache.commons.collections.BidiMap;
import org.apache.commons.collections.bidimap.DualHashBidiMap;

import com.sun.media.sound.RmfFileReader;

public class Pedro {
	
	static final String ID_CNT_ALL= "cntAll";
	static final String ID_CNT_3INCOMPLETE= "cnt3Pincl";
	static final String ID_CNT_REAL= "cntReal";
	static final String ID_CNT_MULT= "cntMulti";
	static final String ID_CNT_SINGLE= "cntSingle";
	static final String ID_CNT_MULT_TRPT= "cntMultTrpt";
	static final String ID_CNT_MULT_TAG= "cntMultTag";
	static final String ID_CNT_MULT_TAG_SAGE= "cntMultTagSAGE";
	static final String ID_CNT_REAL_TRPT= "cntRealTrpt";
	static final String ID_CNT_REAL_SAGE= "cntRealSAGE";
	static final String ID_CNT_MULT_SAGE= "cntMultiSAGE";
	static final String ID_CNT_SINGLE_SAGE= "cntSingleSAGE";
	static final String ID_CNT_MULT_TRPT_SAGE= "cntMultTrptSAGE";
	static final String ID_CNT_REAL_TRPT_SAGE= "cntRealTrptSAGE";
	static final String ID_CNT_SGLE_TRUNC= "cntSglTrunc";
	static final String ID_CNT_MULTI_TRUNC= "cntMultiTrunc";
	static final String ID_CNT_SGLE_REPEAT= "cntSglRepeat";
	static final String ID_CNT_MULTI_REPEAT= "cntMultiRepeat";
	static final String ID_CNT_SGLE_TRUNC_SAGE= "cntSglTruncSAGE";
	static final String ID_CNT_MULTI_TRUNC_SAGE= "cntMultiTruncSAGE";
	static final String ID_CNT_SGLE_REPEAT_SAGE= "cntSglRepeatSAGE";
	static final String ID_CNT_MULTI_REPEAT_SAGE= "cntMultiRepeatSAGE";
	static final String ID_DIST_FOUND= "distFound";
	static final String ID_DIST_NOT_FOUND= "distNotFound";
	static final String ID_CNT_SGLE_TRPT= "cntSingleTrpt";
	static final String ID_CNT_SGLE_TRPT_SAGE= "cntSingleTrptSAGE";
	static final String ID_GENE_WO_TAG_SAGE= "cntGeneWOtagSAGE";
	static final String ID_GENE_WO_TAG_MPSS= "cntGeneWOtagMPSS";
	static final String ID_GENE_INCOMPLETE_MPSS= "cntGeneIncompMPSS";
	static final String ID_GENE_INCOMPLETE_SAGE= "cntGeneIncompSAGE";
	static final String ID_GENE_COMPLETE_MPSS= "cntGeneComplMPSS";
	static final String ID_GENE_COMPLETE_SAGE= "cntGeneComplSAGE";
	
	
	public final static String MPSS_DIGEST_SITE_DPNII= "GATC";
	public final static String MPSS_TAG_QUALIFIER= "mpss_tag";
	public final static int MPSS_TAG_LENGTH= 13;
	public final static String SAGE_DIGEST_SITE_NLAIII= "CATG";
	public final static int SAGE_TAG_LENGTH_LONG= 17;
	public final static int SAGE_TAG_LENGTH_SHORT= 10;
	public final static String SAGE_TAG_QUALIFIER= "sage_tag";
	public final static String ID_TAG_TYPE= "tag_type";
	
	
	public final static String ID_DIGEST= "tag_id"; 
	public final static String ID_TAGLEN= "tag_length"; 

	static HashMap realTagHash= null;
	static DualHashBidiMap virtTagHash= null;

	public static void main(String[] args) {
		String speName= "human";
		String annoName= "ASD";
		String fName= Species.getAnnotation(speName, null, annoName, null);
		//_00_convertToGTF();
		//_01_grepVirtualTags();
		//_02_findTagsInTranscripts(fName);
		
		String[] marker= new String[] {"_tagged"};
		fName= Species.getAnnotation(speName, null, annoName, marker);
		_03_getStatistics(true, fName);
		//_04_diversityPerGene(fName);
		//_05_overlapASEvents(fName);
		_05_assignASevents(fName);
	}
	
	/**	  
	 *	Also checks the input/genome ver compliance
	 */
	public static void _00_convertToGTF() {
		String speName= "human";
		String annoName= "_UCSCGenes";	//"Known";
		String genVer= "200603";	// hg17
		boolean checkSeq= true;
		
		String fName= //Constants.getUCSCAnnotation(speName, annoName, genVer);
			Constants.getLatestUCSCAnnotation(speName, annoName, null);
		System.out.println("Attaching to annotation file "+fName);
		
		TabPedroReader reader= 
			new TabPedroReader("pedro"+File.separator+"mpss_lib_distinct_tags_chr_position.txt");
		//reader.sweepToChromosome("chr1");
		GTFChrReader gtfReader= new GTFChrReader(fName);
		gtfReader.setChromosomeWise(true);
		gtfReader.setReadGene(false);
		gtfReader.setReadGTF(true);
		gtfReader.setSilent(true);
		
		Species spec= new Species(speName);
		spec.setGenomeVersion(genVer);
		
		try {reader.read();} 
		catch (Exception e) {e.printStackTrace();}
		String[][] table= reader.getTable();
		while (table!= null) {
			Date date= new Date(System.currentTimeMillis());
			System.out.println(date.getHours()+":"+date.getMinutes()+" processing "+table[0][0]);
			HashMap tagHash= new HashMap(table.length);
			int cntFalse= 0;
			for (int i = 0; i < table.length; i++) {
				int start= Integer.parseInt(reader.getTable()[i][2]);
				int end= Integer.parseInt(reader.getTable()[i][3]);
				int strand= 1;
				if (end< start)
					strand= -1;
				DirectedRegion reg= new DirectedRegion(start,end,strand);	
				reg.setChromosome(reader.getTable()[i][0]);
				reg.setSpecies(spec);
				
					// check genome version
				String seq= null;
				if (checkSeq) {
					seq= Graph.readSequence(reg);
					if (seq!= null)
						seq= seq.toUpperCase();
				}
				if (checkSeq&& (seq== null|| ((!seq.equals(table[i][1]))&& (!seq.equals("NNNNNNNNNNNNN"))))) {
					++cntFalse;
					System.err.println("Line "+i+" ("+table[i][0]+"): sequence does not match to genome:"
							+ seq+" <> "+table[i][1]+"\t"+start+","+end);
				} else
					tagHash.put(reg, table[i][1]);	// seq
			}
			System.out.println(cntFalse+" of "+table.length+" tag seqs not matching ("+
					Formatter.fprint((cntFalse* 100d)/ table.length, 2)+"%).");
			
				// append info to gtf
			gtfReader.sweepToChromosome(table[0][0]);
			try {gtfReader.read();} 
			catch (Exception e) {e.printStackTrace();}
			GTFObject[] gtfs= gtfReader.getGtfObj();
			boolean[] mapped= new boolean[tagHash.size()];
			for (int i = 0; i < mapped.length; i++) 
				mapped[i]= false;
			Object[] keys= tagHash.keySet().toArray();
			HashMap doubleHash= new HashMap();
			for (int j = 0; gtfs!= null&& j < gtfs.length; j++) {
				if (!gtfs[j].getFeature().equalsIgnoreCase("exon"))
					continue;
				DirectedRegion reg= new DirectedRegion(gtfs[j].getStart(), gtfs[j].getEnd(), gtfs[j].getStrand());
				reg.setChromosome(gtfs[j].getSeqname());
				reg.setSpecies(spec);
				for (int i = 0; i < keys.length; i++) {
					if (reg.overlaps((DirectedRegion) keys[i])) {
						gtfs[j].addAttribute("mpss_tag", (String) tagHash.get(keys[i]));
						mapped[i]= true;
						Integer val= (Integer) doubleHash.get(tagHash.get(keys[i]));
						if (val== null)
							val= new Integer(1);
						else
							val= new Integer(val.intValue()+ 1);
						doubleHash.put(tagHash.get(keys[i]), val);
					}
				}
			}
			int cntMapped= 0;
			for (int i = 0; i < mapped.length; i++) 
				if (mapped[i])
					++cntMapped;
			double val= -1d;
			if (tagHash.size()> 0)
				val= (cntMapped* 100d)/ tagHash.size();
			System.out.println("Mapped "+cntMapped+" of "+tagHash.size()+" tags to genes ("+
					Formatter.fprint(val, 2)+"%), "+doubleHash.size()+ " of them possibly have multiple hits ("+
					Formatter.fprint((doubleHash.size()* 100d)/ cntMapped, 2)+"%).");
			GTFWrapper gtfWriter= new GTFWrapper(fName+"_mpss");	// append
			gtfWriter.setGtfObj(gtfs);
			try {gtfWriter.write(true);}
			catch (Exception e) {e.printStackTrace();}
			
				// next
			try {reader.read();} 
			catch (Exception e) {e.printStackTrace();}
			table= reader.getTable();
		}
	}
	
	public static void _01_grepVirtualTags() {
		String subdir= "pedro";
		String outfile= subdir+ File.separator+
							"mpss_lib_distinct_tags_chr_position.txt";
		File outFile= new File(outfile);
		if (outFile.exists())
			outFile.delete();
		System.out.print("Reading distinct tags..");
		System.out.flush();
		HashMap tagHash= getRealTagHash();
		System.out.println("done. Found "+tagHash.size()+" tags.");
		try {
			String[] dir= new File(subdir).list();
			System.out.println("Scanning..");
			for (int i = 0; i < dir.length; i++) {
				BufferedWriter writer= new BufferedWriter(new FileWriter(outfile, true));
				if (!dir[i].startsWith("result_"))
					continue;
				System.out.println("\t"+dir[i]);
				BufferedReader buffy= new BufferedReader(new FileReader(subdir+ File.separator+ dir[i]));
				while (buffy.ready())  {
					String line= buffy.readLine();
					String[] tokens= line.split("\t");
					if (tagHash.get(tokens[1])!= null)
						writer.write(line+"\n");
				}
				writer.flush(); writer.close();
				buffy.close();
			}
			System.out.println("done.");

		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static HashMap getRealTagHash() {
		String fName= "pedro"+ File.separator+ "mpss_lib_distinct_tags.sql";
		if (realTagHash == null) {
			try {
				realTagHash= new HashMap(341000);
				BufferedReader buffy= new BufferedReader(new FileReader(fName));
				buffy.readLine();	// skip first line				
				while (buffy.ready()) {
					String line= buffy.readLine().trim();
					realTagHash.put(line, line);
				}
				buffy.close();
				int totalTags= realTagHash.size();
				
				Object[] keys= realTagHash.keySet().toArray();
				for (int i = 0; i < keys.length; i++) {
					String seqR= Arrays.reverseComplement((String) keys[i]);
					if (realTagHash.get(seqR)!= null) {
						realTagHash.remove(seqR);
						realTagHash.remove(keys[i]);
					}
						
				}
				System.out.println("Read "+totalTags+" real MPSS tags, "+
						realTagHash.size()+ " of them unique.");
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return realTagHash;
	}

	public static GTFObject[] _02_findTagInTranscripts(Gene g, String digSite, Integer tagLen, int nbTags, HashMap localHash, boolean countMPSS) {
		Integer cntAll= (Integer) localHash.get(ID_CNT_ALL);
		if (cntAll== null)
			cntAll= new Integer(0);
		Integer cnt3Pincl= (Integer) localHash.get(ID_CNT_3INCOMPLETE);
		if (cnt3Pincl== null)
			cnt3Pincl= new Integer(0);
		IntVector distNotFound= (IntVector) localHash.get(ID_DIST_NOT_FOUND);
		if (distNotFound== null)
			distNotFound= new IntVector();
		IntVector distFound= (IntVector) localHash.get(ID_DIST_FOUND);
		if (distFound== null)
			distFound= new IntVector();
		
		
		Integer cntReal;
		if (countMPSS)
			cntReal= (Integer) localHash.get(ID_CNT_REAL);
		else
			cntReal= (Integer) localHash.get(ID_CNT_REAL_SAGE);
		if (cntReal== null)
			cntReal= new Integer(0);
		Integer cntMulti;
		if (countMPSS)
			cntMulti= (Integer) localHash.get(ID_CNT_MULT);
		else
			cntMulti= (Integer) localHash.get(ID_CNT_MULT_SAGE);
		if (cntMulti== null)
			cntMulti= new Integer(0);
		Integer cntSingle;
		if (countMPSS)
			cntSingle= (Integer) localHash.get(ID_CNT_SINGLE);
		else
			cntSingle= (Integer) localHash.get(ID_CNT_SINGLE_SAGE);
		if (cntSingle== null)
			cntSingle= new Integer(0);
		Integer cntSingleTrpt;
		if (countMPSS)
			cntSingleTrpt= (Integer) localHash.get(ID_CNT_SGLE_TRPT);
		else
			cntSingleTrpt= (Integer) localHash.get(ID_CNT_SGLE_TRPT_SAGE);
		if (cntSingleTrpt== null)
			cntSingleTrpt= new Integer(0);
		Integer cntMultiTrpt;
		if (countMPSS)
			cntMultiTrpt= (Integer) localHash.get(ID_CNT_MULT_TRPT);
		else
			cntMultiTrpt= (Integer) localHash.get(ID_CNT_MULT_TRPT_SAGE);
		if (cntMultiTrpt== null)
			cntMultiTrpt= new Integer(0);
		Integer cntRealTrpt;
		if (countMPSS)
			cntRealTrpt= (Integer) localHash.get(ID_CNT_REAL_TRPT);
		else
			cntRealTrpt= (Integer) localHash.get(ID_CNT_REAL_TRPT_SAGE);
		if (cntRealTrpt== null)
			cntRealTrpt= new Integer(0);
		Integer cntSingleTrunc;
		if (countMPSS)
			cntSingleTrunc= (Integer) localHash.get(ID_CNT_SGLE_TRUNC);
		else
			cntSingleTrunc= (Integer) localHash.get(ID_CNT_SGLE_TRUNC_SAGE);
		if (cntSingleTrunc== null)
			cntSingleTrunc= new Integer(0);
		Integer cntMultiTrunc;
		if (countMPSS)
			cntMultiTrunc= (Integer) localHash.get(ID_CNT_MULTI_TRUNC);
		else
			cntMultiTrunc= (Integer) localHash.get(ID_CNT_MULTI_TRUNC_SAGE);
		if (cntMultiTrunc== null)
			cntMultiTrunc= new Integer(0);
		Integer cntSingleRepeat;
		if (countMPSS)
			cntSingleRepeat= (Integer) localHash.get(ID_CNT_SGLE_REPEAT);
		else
			cntSingleRepeat= (Integer) localHash.get(ID_CNT_SGLE_REPEAT_SAGE);
		if (cntSingleRepeat== null)
			cntSingleRepeat= new Integer(0);
		Integer cntMultiRepeat;
		if (countMPSS)
			cntMultiRepeat= (Integer) localHash.get(ID_CNT_MULTI_REPEAT);
		else
			cntMultiRepeat= (Integer) localHash.get(ID_CNT_MULTI_REPEAT_SAGE);
		if (cntMultiRepeat== null)
			cntMultiRepeat= new Integer(0);
		Integer cntMultiTags;
		if (countMPSS)
			cntMultiTags= (Integer) localHash.get(ID_CNT_MULT_TAG);
		else
			cntMultiTags= (Integer) localHash.get(ID_CNT_MULT_TAG_SAGE);
		if (cntMultiTags== null)
			cntMultiTags= new Integer(0);
		Integer cntGeneWOtag;
		if (countMPSS)
			cntGeneWOtag= (Integer) localHash.get(ID_GENE_WO_TAG_MPSS);
		else
			cntGeneWOtag= (Integer) localHash.get(ID_GENE_WO_TAG_SAGE);
		if (cntGeneWOtag== null)
			cntGeneWOtag= new Integer(0);
		Integer cntGenesIncompl;
		if (countMPSS)
			cntGenesIncompl= (Integer) localHash.get(ID_GENE_INCOMPLETE_MPSS);
		else
			cntGenesIncompl= (Integer) localHash.get(ID_GENE_INCOMPLETE_SAGE);
		if (cntGenesIncompl== null)
			cntGenesIncompl= new Integer(0);
		Integer cntGenesComplete;
		if (countMPSS)
			cntGenesComplete= (Integer) localHash.get(ID_GENE_COMPLETE_MPSS);
		else
			cntGenesComplete= (Integer) localHash.get(ID_GENE_COMPLETE_SAGE);
		if (cntGenesComplete== null)
			cntGenesComplete= new Integer(0);

		
		
		boolean found= false;
		Vector obs= new Vector();
		if (g!= null) {
//			String digSite= (String) map.get(ID_DIGEST);
//			Integer tagLen= (Integer) map.get(ID_TAGLEN);
			if (digSite== null|| tagLen== null)
				return null;
			Vector v= new Vector();
			int trptNotFound= 0;
			int trptIncomplete= 0;
			boolean completeTranscript= false;
			for (int i = 0; i < g.getTranscripts().length; i++) {
				cntAll= new Integer(cntAll+ 1);
				Transcript trpt= g.getTranscripts()[i];
				String seq= trpt.getSplicedSequence();
					
					// check for 3'completeness
				// check complete3PSandro()
//				int pos;
//				for (pos = seq.length()- 1; pos >0; --pos) 
//					if (seq.charAt(pos)!= 'A')
//						break;		// trim polyA
//				
//				if (pos+3> seq.length())	// pA signal ends with "AAA"
//					pos= seq.length();
//				else 
//					pos= pos+3;
//				
//				String window= seq.substring(pos-30, pos);
//				if (!window.contains("AATAAA")&& !window.contains("ATTAAA")) {
//					cnt3Pincl= new Integer(cnt3Pincl+ 1);
//					continue;	// not 3' complete
//				}
				if (!trpt.is3Pcomplete()) {
					cnt3Pincl= new Integer(cnt3Pincl+ 1);
					++trptIncomplete;
					continue;
				}
				
				completeTranscript= true;
				int pos= seq.length();
				int cntTags= 0;
				for (int j = pos-4; j >0; --j) {	// KMP?
					if (seq.substring(j, j+4).equalsIgnoreCase(digSite)) {
						int end= Math.min(j+4+tagLen.intValue(), seq.length());
						String tag= seq.substring(j+4, end);
						while (tag.length()< tagLen.intValue()) 
							tag+= "X";	// fill with unknown
						
						GTFObject obj= new GTFObject();
						obj.setStart(trpt.getGenomicPosition(j+4));
						obj.setEnd(trpt.getGenomicPosition(end- 1));
						obj.setStrand(trpt.getStrand());
						obj.setFeature(getDigSiteQualifier(digSite));
						obj.setSeqname(trpt.getChromosome());
						obj.setSource(trpt.getSource());
						obj.addAttribute(GTFObject.GENE_ID_TAG, trpt.getGene().getGeneID());
						obj.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, trpt.getTranscriptID());
						obj.addAttribute(GTFObject.ID_ATTRIBUTE_SEQUENCE, tag);
						obs.add(obj);
						
						trpt.addAttribute(digSite, tag);
						v= (Vector) Arrays.addUnique(v, tag);
						found= true;
						++cntTags;
						if (cntTags== nbTags)
							break;
					}
				}
				if (found)
					distFound.add(seq.length());
				else {
					//System.out.println("No tag found for "+trpt.getTranscriptID());
					distNotFound.add(seq.length());
					++trptNotFound;
				}
			}
			
			if (completeTranscript)
				cntGenesComplete= new Integer(cntGenesComplete+ 1);
			
				// count stats
			if (v.size()> 1) {
//				System.out.print("Gene with transcript "+g.getTranscripts()[0].getTranscriptID()+
//						" has "+v.size()+" different "+digSite+" tags.");
				cntMulti= new Integer(cntMulti.intValue()+ 1);
				cntMultiTrpt= new Integer(cntMultiTrpt.intValue()+ g.getTranscriptCount()- trptIncomplete);
				cntMultiTags= new Integer(cntMultiTags.intValue()+ v.size());
				
					// filter real ones
				int cntRealOnes= 0;				
				for (int i = 0; i < v.size(); i++) {
					if (getRealTagHash().get(v.elementAt(i))!= null)
						++cntRealOnes;
					String tag= (String) v.elementAt(i);
					if (tag.contains("X"))
						cntMultiTrunc= new Integer(cntMultiTrunc+ 1);
					for (int j = 0; j < tag.length(); j++) 
						if (Character.isLowerCase(tag.charAt(j))) {
							cntMultiRepeat= new Integer(cntMultiRepeat+ 1);
							break;
						}
				}
				if (cntRealOnes> 1) {
					cntReal= new Integer(cntReal+ 1);
					cntRealTrpt= new Integer(cntRealTrpt.intValue()+ cntRealOnes);
				} 
			} else if (v.size()== 1) {
				cntSingle= new Integer(cntSingle+ 1);
				cntSingleTrpt= new Integer(cntSingleTrpt+ g.getTranscriptCount()- trptIncomplete);
				String tag= (String) v.elementAt(0);
				if (tag.contains("X"))
					cntSingleTrunc= new Integer(cntSingleTrunc+ 1);
				for (int i = 0; i < tag.length(); i++) 
					if (Character.isLowerCase(tag.charAt(i))) {
						cntSingleRepeat= new Integer(cntSingleRepeat+ 1);
						break;
					}
			} else if (completeTranscript) {
				if (v.size()== 0)
					cntGeneWOtag= new Integer(cntGeneWOtag+ 1);
			} else 
				++cntGenesIncompl;
				
		}
		
		if (countMPSS) {
			localHash.put(ID_CNT_ALL, cntAll);
			localHash.put(ID_CNT_3INCOMPLETE, cnt3Pincl);
			localHash.put(ID_DIST_FOUND, distFound);
			localHash.put(ID_DIST_NOT_FOUND, distNotFound);
			localHash.put(ID_CNT_MULT, cntMulti);
			localHash.put(ID_CNT_MULT_TAG, cntMultiTags);
			localHash.put(ID_CNT_REAL, cntReal);
			localHash.put(ID_CNT_SINGLE, cntSingle);
			localHash.put(ID_CNT_SGLE_TRPT, cntSingleTrpt);
			localHash.put(ID_CNT_MULT_TRPT, cntMultiTrpt);
			localHash.put(ID_CNT_REAL_TRPT, cntRealTrpt);
			localHash.put(ID_CNT_SGLE_TRUNC, cntSingleTrunc);
			localHash.put(ID_CNT_MULTI_TRUNC, cntMultiTrunc);
			localHash.put(ID_CNT_SGLE_REPEAT, cntSingleRepeat);
			localHash.put(ID_CNT_MULTI_REPEAT, cntMultiRepeat);
			localHash.put(ID_GENE_WO_TAG_MPSS, cntGeneWOtag);
			localHash.put(ID_GENE_INCOMPLETE_MPSS, cntGenesIncompl);
			localHash.put(ID_GENE_COMPLETE_MPSS, cntGenesComplete);
		} else {
			localHash.put(ID_DIST_FOUND, distFound);
			localHash.put(ID_DIST_NOT_FOUND, distNotFound);
			localHash.put(ID_CNT_MULT_SAGE, cntMulti);
			localHash.put(ID_CNT_MULT_TAG_SAGE, cntMultiTags);
			localHash.put(ID_CNT_REAL_SAGE, cntReal);
			localHash.put(ID_CNT_SINGLE_SAGE, cntSingle);
			localHash.put(ID_CNT_SGLE_TRPT_SAGE, cntSingleTrpt);
			localHash.put(ID_CNT_MULT_TRPT_SAGE, cntMultiTrpt);
			localHash.put(ID_CNT_REAL_TRPT_SAGE, cntRealTrpt);
			localHash.put(ID_CNT_SGLE_TRUNC_SAGE, cntSingleTrunc);
			localHash.put(ID_CNT_MULTI_TRUNC_SAGE, cntMultiTrunc);
			localHash.put(ID_CNT_SGLE_REPEAT_SAGE, cntSingleRepeat);
			localHash.put(ID_CNT_MULTI_REPEAT_SAGE, cntMultiRepeat);
			localHash.put(ID_GENE_WO_TAG_SAGE, cntGeneWOtag);
			localHash.put(ID_GENE_INCOMPLETE_SAGE, cntGenesIncompl);
			localHash.put(ID_GENE_COMPLETE_SAGE, cntGenesComplete);
		}
		return (GTFObject[]) Arrays.toField(obs);
	}
	
	/**
	 * extracts the tags according to a given annotation and writes 
	 * a gtf file with the tags
	 *
	 */
	public static void _02_findTagsInTranscripts(String fName) {
		try {
			for (int x = 1; x < 11; x++) {
				System.out.println("Tagging "+fName);
				String outName= fName+"_tagged"+x;
				File outFile= new File(outName);
				if (outFile.exists()) {
					System.out.println("Confirm removing file "+outName);
					int c= System.in.read();
					if (c== 'y'|| c== 'Y')
						outFile.delete();
					else
						System.exit(0);
				}
				GTFChrReader reader= new GTFChrReader(fName);
				reader.setChromosomeWise(true);
				reader.setReadGene(true);
				reader.read();
				Gene[] genes= reader.getGenes();
				int cnt= 0;
				HashMap localHash= new HashMap(4);
				while (genes!= null) {
					cnt+= genes.length;
					System.out.print("+");
					System.out.flush();
					Vector v= new Vector();
					for (int i = 0; i < genes.length; i++) {
						GTFObject[] obs= _02_findTagInTranscripts(genes[i], MPSS_DIGEST_SITE_DPNII, new Integer(MPSS_TAG_LENGTH), x, localHash, true);
						HashMap mpssHash= new HashMap();	
						for (int j = 0; obs!= null&& j < obs.length; j++) 
							mpssHash.put(obs[j].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), obs[j]);
						
						obs= _02_findTagInTranscripts(genes[i], SAGE_DIGEST_SITE_NLAIII, new Integer(SAGE_TAG_LENGTH_LONG), x, localHash, false);
						HashMap sageHash= new HashMap();	
						for (int j = 0; obs!= null&& j < obs.length; j++) 
							sageHash.put(obs[j].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), obs[j]);
						
						for (int j = 0; j < genes[i].getTranscripts().length; j++) {
							GTFObject[] trObs= GTFWrapper.getGTFObjects(genes[i].getTranscripts()[j]);
							for (int k = 0; k < trObs.length; k++) 
								v.add(trObs[k]);
							if (mpssHash.get(genes[i].getTranscripts()[j].getTranscriptID())!= null)
								v.add(mpssHash.get(genes[i].getTranscripts()[j].getTranscriptID()));
							if (sageHash.get(genes[i].getTranscripts()[j].getTranscriptID())!= null)
								v.add(sageHash.get(genes[i].getTranscripts()[j].getTranscriptID()));
						}
					}
					
					if (v.size()> 0) {
						GTFWrapper writer= new GTFWrapper(outName);
						writer.setGtfObj((GTFObject[]) Arrays.toField(v));
						writer.setSortAttributes(new String[] {GTFObject.GENE_ID_TAG, GTFObject.TRANSCRIPT_ID_TAG});
						writer.write(true);					
					}
					
					reader.read();
					genes= reader.getGenes();
				}			
				
				
				Integer cntAll= (Integer) localHash.get(ID_CNT_ALL);
				if (cntAll== null)
					cntAll= new Integer(0);
				Integer cnt3Pincl= (Integer) localHash.get(ID_CNT_3INCOMPLETE);
				if (cnt3Pincl== null)
					cnt3Pincl= new Integer(0);
				Integer cntReal= (Integer) localHash.get(ID_CNT_REAL);
				Integer cntMulti= (Integer) localHash.get(ID_CNT_MULT);
				Integer cntSingle= (Integer) localHash.get(ID_CNT_SINGLE);
				Integer cntMultiTrpt= (Integer) localHash.get(ID_CNT_MULT_TRPT);
				Integer cntMultiTags= (Integer) localHash.get(ID_CNT_MULT_TAG);
				Integer cntMultiTagsSAGE= (Integer) localHash.get(ID_CNT_MULT_TAG_SAGE);
				Integer cntRealTrpt= (Integer) localHash.get(ID_CNT_REAL_TRPT);
				Integer cntRealSAGE= (Integer) localHash.get(ID_CNT_REAL_SAGE);
				Integer cntMultiSAGE= (Integer) localHash.get(ID_CNT_MULT_SAGE);
				Integer cntSingleSAGE= (Integer) localHash.get(ID_CNT_SINGLE_SAGE);
				Integer cntSingleTrpt= (Integer) localHash.get(ID_CNT_SGLE_TRPT);
				Integer cntSingleTrptSAGE= (Integer) localHash.get(ID_CNT_SGLE_TRPT_SAGE);
				Integer cntMultiTrptSAGE= (Integer) localHash.get(ID_CNT_MULT_TRPT_SAGE);
				Integer cntRealTrptSAGE= (Integer) localHash.get(ID_CNT_REAL_TRPT_SAGE);
				Integer cntSingleTrunc= (Integer) localHash.get(ID_CNT_SGLE_TRUNC);
				Integer cntSingleTruncSAGE= (Integer) localHash.get(ID_CNT_SGLE_TRUNC_SAGE);
				Integer cntSingleRepeatSAGE= (Integer) localHash.get(ID_CNT_SGLE_REPEAT_SAGE);
				Integer cntSingleRepeat= (Integer) localHash.get(ID_CNT_SGLE_REPEAT);
				Integer cntMultiRepeat= (Integer) localHash.get(ID_CNT_MULTI_REPEAT);
				Integer cntMultiRepeatSAGE= (Integer) localHash.get(ID_CNT_MULTI_REPEAT_SAGE);
				Integer cntMultiTrunc= (Integer) localHash.get(ID_CNT_MULTI_TRUNC);
				Integer cntMultiTruncSAGE= (Integer) localHash.get(ID_CNT_MULTI_TRUNC_SAGE);
				IntVector distFound= (IntVector) localHash.get(ID_DIST_FOUND);
				IntVector distNotFound= (IntVector) localHash.get(ID_DIST_NOT_FOUND);
				Integer geneWOtagMPSS= (Integer) localHash.get(ID_GENE_WO_TAG_MPSS);
				Integer geneWOtagSAGE= (Integer) localHash.get(ID_GENE_WO_TAG_SAGE);
				Integer cntGenesIncomplMPSS= (Integer) localHash.get(ID_GENE_INCOMPLETE_MPSS);
				Integer cntGenesIncomplSAGE= (Integer) localHash.get(ID_GENE_INCOMPLETE_SAGE);
				Integer cntGenesCompleteMPSS= (Integer) localHash.get(ID_GENE_COMPLETE_MPSS);
				Integer cntGenesCompleteSAGE= (Integer) localHash.get(ID_GENE_COMPLETE_SAGE);
				System.out.println("\n\nChecked "+cnt+" genes ("+cntAll.intValue()+" transcripts)\n\t"
						+cnt3Pincl.intValue()+" trpts 3'incomplete\nMPSS\t"
						+cntMultiTrpt.intValue()+" transcripts with multiple tags, " 
						+cntMultiTags+" tags, of them "
						+cntRealTrpt.intValue()+" found in real tags.\n\t"
						+geneWOtagMPSS+ " genes without any tag, "+cntGenesIncomplMPSS+" incomplete ones, "+cntGenesCompleteMPSS+" complete ones\n\t"
						+cntSingle.intValue()+" genes with a single tag ("
						+cntSingleTrpt+ " trpts), of them "
						+cntSingleTrunc+" truncated and "
						+cntSingleRepeat+" in repeat regions\n\t"
						+cntMulti.intValue()+" genes with multiple tags, of them "
						+cntMultiTrunc.intValue()+" truncated and "
						+cntMultiRepeat+ " in repeats and "
						+cntReal.intValue()+" found in real tags.\nSAGE\t"
				+cntMultiTrptSAGE.intValue()+" transcripts with multiple tags, "
				+cntMultiTags+" tags, of them "
				+cntRealTrptSAGE.intValue()+" found in real tags.\n\t"
				+geneWOtagSAGE+ " genes without any tag "+cntGenesIncomplSAGE+" incomplete ones, "+cntGenesCompleteSAGE+" complete ones\n\t"
				+cntSingleSAGE.intValue()+" genes with a single tag (" 
				+cntSingleTrptSAGE+ " trpts), of them "
				+cntSingleTruncSAGE+" truncated and "
				+cntSingleRepeatSAGE+" in repeat regions\n\t"
				+cntMultiSAGE.intValue()+" genes with multiple tags, of them "
				+cntMultiTruncSAGE.intValue()+" truncated and "
				+cntMultiRepeatSAGE+ " in repeats and "
				+cntRealSAGE.intValue()+" found in real tags.\nSAGE\t");
				
				System.out.println("In total "+ distNotFound.size()+ 
						" transcripts had no MPSS or SAGE tag, med len "+ new Distribution(distNotFound.toIntArray()).getMedian()+
						", of the ones with tags med len "+ new Distribution(distFound.toIntArray()).getMedian());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void _03_getStatistics(boolean filter, String fName) {
			
		DualHashBidiMap vMap= getVirtualTagHash(fName);
		HashMap rMap= getRealTagHash();
		
			// check virtual found in real
		Object[] keys= vMap.keySet().toArray();
		int cntMPSS= 0, cntSAGE= 0, cntMPSStot= 0, cntSAGEtot= 0;
		for (int i = 0; i < keys.length; i++) {
			if (isMPSStag(((String) keys[i]).length()))
				++cntMPSStot;
			else {
				++cntSAGEtot;
				continue;
			}
			String seqR= Arrays.reverseComplement((String) keys[i]);
			if (rMap.get(keys[i])== null&& rMap.get(seqR)== null) {
				if (filter) {
					vMap.remove(keys[i]);
					vMap.remove(seqR);
				}
				continue;
			}
			if (isMPSStag(((String) keys[i]).length()))
				++cntMPSS;
			else
				++cntSAGE;
		}
		System.out.println("Found "+cntMPSS+" virtual MPSS tags in real. ("+Formatter.fprint(100d* cntMPSS/cntMPSStot, 2)+"%)");
		System.out.println("Found "+cntSAGE+" virtual SAGE tags in real. ("+Formatter.fprint(100d* cntSAGE/cntSAGEtot, 2)+"%)");

			// check real found in virtual
		keys= rMap.keySet().toArray();
		cntMPSS= 0; cntSAGE= 0; cntMPSStot= 0; cntSAGEtot= 0;
		for (int i = 0; i < keys.length; i++) {
			if (isMPSStag(((String) keys[i]).length()))
				++cntMPSStot;
			else {
				++cntSAGEtot;
				continue;
			}
			String seqR= Arrays.reverseComplement((String) keys[i]);
			if (vMap.get(keys[i])== null&& vMap.get(seqR)== null) {
				if (filter) {
					rMap.remove(keys[i]);
					rMap.remove(seqR);
				}
				continue;
			}
			if (isMPSStag(((String) keys[i]).length()))
				++cntMPSS;
			else
				++cntSAGE;
		}
		System.out.println("Found "+cntMPSS+" real MPSS tags in virtuals. ("+Formatter.fprint(100d* cntMPSS/cntMPSStot, 2)+"%)");
		System.out.println("Found "+cntSAGE+" real SAGE tags in virtuals. ("+Formatter.fprint(100d* cntSAGE/cntSAGEtot, 2)+"%)");
	}
		
	public static void _04_diversityPerGene(String fName) {
			
		try {
			DualHashBidiMap vTags= getVirtualTagHash(fName);
			HashMap tIDhash= new HashMap(vTags.size());
			Object[] keys= vTags.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				if (!isMPSStag(((String) keys[i]).length()))
					continue;
				DirectedRegion reg= (DirectedRegion) vTags.get(keys[i]);
				Vector v= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
				for (int j = 0; j < v.size(); j++) 
					tIDhash.put(v.elementAt(j), keys[i]);
			}
	
			System.out.println("Reading "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			int cntGenes= 0;
			IntVector cntMultiplesV= new IntVector();
			
			IntVector distTrptGrps= new IntVector();
			while (genes!= null) {
				cntGenes+= genes.length;
				for (int i = 0; i < genes.length; i++) {
					Transcript[] trpts= genes[i].getTranscripts();
					int cntTrptUniqueTag= 0;
					HashMap tagsInGene= new HashMap();
					for (int j = 0; j < trpts.length; j++) 
						if (tIDhash.get(trpts[j].getTranscriptID())!= null&&
								tagsInGene.get(tIDhash.get(trpts[j].getTranscriptID()))== null)
							tagsInGene.put(tIDhash.get(trpts[j].getTranscriptID()), vTags.get(tIDhash.get(trpts[j].getTranscriptID())));
					
					cntMultiplesV.putValue(tagsInGene.size(), cntMultiplesV.getValue(tagsInGene.size())+1);
					if (tagsInGene.size()> 1) {
						Object[] trptPerTag= tagsInGene.values().toArray();
						for (int j = 0; j < trptPerTag.length; j++) {
							DirectedRegion reg= (DirectedRegion) trptPerTag[j];
							Vector idV= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
							distTrptGrps.putValue(idV.size(), distTrptGrps.getValue(idV.size())+ 1);
						}
					}
				}
				reader.read();
				genes= reader.getGenes();
			}
			
			System.out.println("Unique tags per gene:");
			int[] cntMultiples= cntMultiplesV.toIntArray();
			for (int i = 0; i < cntMultiples.length; i++) {
				System.out.println(i+"\t"+cntMultiples[i]);
			}
			System.out.println("with transcripts per unique tag in a gene (only for genes with >1 tag):");
			IntVector allGroupV= new IntVector();
			int sum= 0;
			for (int i = 1; i < distTrptGrps.length; i++) {
				System.out.print(i+" "+distTrptGrps.get(i)+",");
				for (int j = 0; j < distTrptGrps.get(i); j++) 
					allGroupV.add(i);
				sum+= distTrptGrps.get(i);
			}
			System.out.println("Median "+ new Distribution(allGroupV.toIntArray()).getMedian()+", sum "+sum);	
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void _05_assignASevents(String fName) {
			
		try {
			DualHashBidiMap vTags= getVirtualTagHash(fName);	// tag x reg(tIDs..)
			HashMap tIDhash= new HashMap(vTags.size());	// tID x tag
			Object[] keys= vTags.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				if (!isMPSStag(((String) keys[i]).length()))
					continue;
				DirectedRegion reg= (DirectedRegion) vTags.get(keys[i]);
				Vector v= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
				for (int j = 0; j < v.size(); j++) 
					tIDhash.put(v.elementAt(j), keys[i]);
			}
	
			System.out.println("Reading "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			HashMap eventToGroupMap= new HashMap(); // asEv x V(grps)
			int cntIncInGroup= 0, cntIncOutGroup= 0, cntTotGrp= 0, cntHalfUnique= 0, cntUnique= 0, cntNoneUnique= 0;
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					HashMap grpHash= genes[i].getVariantGroups();	// var x Vector[]{varTIDs}
					keys= grpHash.keySet().toArray();
					for (int j = 0; keys!= null&& j < keys.length; j++) {
						Vector[] grps= (Vector[]) grpHash.get(keys[j]);
						for (int x = 0; x < grps.length; x++) {
							String[] tIDs= (String[]) Arrays.toField(grps[x]);
							String currTag= null;
							HashMap littleHash= new HashMap(tIDs.length);
							int k;
							for (k = 0; k < tIDs.length; k++) {	// all of the group must have same unique tag
								littleHash.put(tIDs[k], tIDs[k]);
								if (tIDhash.get(tIDs[k])== null)
									continue;
								if (currTag== null) { 
									currTag= (String) tIDhash.get(tIDs[k]);
									continue;
								}
								if (!currTag.equals(tIDhash.get(tIDs[k]))) 
									break;
							}
							if (currTag== null|| k< tIDs.length) {	// none or diverging tags found
								++cntIncInGroup;
								grps[x]= null;
								continue;
							}
							for (k = 0; k < genes[i].getTranscripts().length; k++) {	// all other trpts must have a distinct tag
								if (littleHash.get(genes[i].getTranscripts()[k].getTranscriptID())!= null||
										tIDhash.get(genes[i].getTranscripts()[k].getTranscriptID())== null)	// trpt in group or has a null tag
									continue;
								if (currTag.equals(tIDhash.get(genes[i].getTranscripts()[k].getTranscriptID()))) {
									++cntIncOutGroup;
									break;
								}
							}
							if (k< genes[i].getTranscripts().length)
								grps[x]= null;
							
						}
						if (grps[0]== null^grps[1]== null)
							++cntHalfUnique;
						else if (grps[0]== null&& grps[1]== null)
							++cntNoneUnique;
						else
							++cntUnique;
						eventToGroupMap.put(keys[j], grps);
					}
	
				}
				reader.read();
				genes= reader.getGenes();
			}
			System.out.println(cntIncInGroup+" incons in group, "+cntIncOutGroup+" incons out group.");
			System.out.println(cntNoneUnique+" none unique ev, "+cntHalfUnique+" half unique ev, "+cntUnique+" unique events.");
			
			keys= eventToGroupMap.keySet().toArray();
			HashMap evClassMap= new HashMap();
			int cntDouble= 0;			
			for (int i = 0; i < keys.length; i++) {
				Vector[] v= (Vector[]) eventToGroupMap.get(keys[i]);
				if (v[0]!= null&& v[1]!= null) {
					Vector vv= (Vector) evClassMap.get(keys[i].toString());
					if (vv== null)
						evClassMap.put(keys[i].toString(), vv= new Vector());
					vv.add(keys[i]);
					System.out.println("Event "+keys[i]);
					System.out.println(((ASVariation) keys[i]).toStringCoordinates());
					for (int j = 0; j < v.length; j++) {
						System.out.print("\t");
						for (int k = 0; k < v[j].size(); k++) 
							System.out.print(v[j].elementAt(k)+" ");						
						System.out.println(" =>  "+tIDhash.get(v[j].elementAt(0)));
					} 
					System.out.println();
				}
			}
			
			keys= evClassMap.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				System.out.println(((Vector) evClassMap.get(keys[i])).size()+"\t"+keys[i]);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	
	/**
	 * @deprecated
	 */
	public static void _05_assignASevents_spliceChains(String fName) {
			
		try {
			DualHashBidiMap vTags= getVirtualTagHash(fName);	// tag x reg(tIDs..)
			HashMap tIDhash= new HashMap(vTags.size());	// tID x tag
			Object[] keys= vTags.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				if (!isMPSStag(((String) keys[i]).length()))
					continue;
				DirectedRegion reg= (DirectedRegion) vTags.get(keys[i]);
				Vector v= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
				for (int j = 0; j < v.size(); j++) 
					tIDhash.put(v.elementAt(j), keys[i]);
			}
	
			System.out.println("Reading "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			HashMap eventToGroupMap= new HashMap(); // asEv x V(grps)
			int cntSameTag= 0, cntRedFilt= 0, cntIncInGroup= 0, cntIncOutGroup= 0, cntTotGrp= 0;
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					VariantGroup_SpliceChain[] varGrps= genes[i].getVariantGroups_spliceChains();
					Vector candidateGrps= new Vector();
					for (int j = 0; varGrps!= null&& j < varGrps.length; j++) {
						++cntTotGrp;
						String[] tIDs= varGrps[j].getTranscriptIDs();
						String currTag= null;
						HashMap littleHash= new HashMap(tIDs.length);
						int k;
						for (k = 0; k < tIDs.length; k++) {	// all of the group must have same unique tag
							littleHash.put(tIDs[k], tIDs[k]);
							if (k== 0) {
								currTag= (String) tIDhash.get(tIDs[k]);
								if (currTag== null) 
									break;
								continue;
							}
							if (!currTag.equals(tIDhash.get(tIDs[k]))) 
								break;
						}
						if (k< tIDs.length) {
							++cntIncInGroup;
							continue;
						}
						for (k = 0; k < genes[i].getTranscripts().length; k++) {	// all other trpts must have a distinct tag
							if (littleHash.get(genes[i].getTranscripts()[k].getTranscriptID())!= null)
								continue;
							if (currTag.equals(genes[i].getTranscripts()[k].getTranscriptID())) {
								++cntIncOutGroup;
								break;
							}
						}
						if (k== genes[i].getTranscripts().length)
							candidateGrps.add(varGrps[j]);
					}
	
						// intersect events of candidate groups
					HashMap evToGrpMap= new HashMap(); // asEv x V2(grps)
					for (int j = 0; j < candidateGrps.size(); j++) {
						VariantGroup_SpliceChain grp= (VariantGroup_SpliceChain) candidateGrps.elementAt(j);
						String grpTag= (String) tIDhash.get(grp.getTranscriptIDs()[0]);
						if (grpTag== null || !isMPSStag(grpTag.length()))
							continue;	// throw out in-specific ones
						ASVariation[] vars= grp.getASevents();
						for (int k = 0; k < vars.length; k++) {
							Vector grpV= (Vector) evToGrpMap.get(vars[k]);
							if (grpV== null) {
								grpV= new Vector();
								grpV.add(grp);
								evToGrpMap.put(vars[k], grpV);
							} else {
									// check for non-identical IDs between the groups
								int m;
								for (m = 0; m < grpV.size(); m++) {
									String tag= (String) tIDhash.get(((VariantGroup_SpliceChain) grpV.elementAt(m)).getTranscriptIDs()[0]);
									if (grpTag.equals(tag))
										break;
								}
								if (m< grpV.size()) {
									evToGrpMap.remove(vars[k]);
									++cntSameTag;
								} else
									grpV.add(grp);
							}
						}
					}
					
						// remove redundancy and add to result
					ASVariation[] vars= (ASVariation[]) Arrays.toField(evToGrpMap.keySet().toArray());
					int beforeRedFilt= 0;
					if (vars!= null)
						beforeRedFilt= vars.length;
					vars= ASMultiVariation.removeRedundancy(
							vars, 
							new ASVariation.StructureComparator());
					if (vars== null)
						cntRedFilt+= (beforeRedFilt- 0);
					else
						cntRedFilt+= (beforeRedFilt- vars.length);
					
					for (int j = 0; vars!= null&& j < vars.length; j++) 
						eventToGroupMap.put(vars[j], evToGrpMap.get(vars[j]));
				}
				reader.read();
				genes= reader.getGenes();
			}
			System.out.println(cntTotGrp+" variant groups, "+cntIncInGroup+" incons in group, "+cntIncOutGroup+" incons out group.");
			System.out.println("Removed "+cntSameTag+" (redundant) events due to non-distinguishable tags.");
			System.out.println("Removed "+cntRedFilt+" redundant events.");
			
			keys= eventToGroupMap.keySet().toArray();
			System.out.println("Unique events: "+keys.length);
			HashMap evClassMap= new HashMap();
			int cntDouble= 0;			
			for (int i = 0; i < keys.length; i++) {
				Vector trptGrpV= (Vector) eventToGroupMap.get(keys[i]);
				boolean print= false;
				if (trptGrpV.size()> 1) {
					print= true;
					++cntDouble;					
					Vector v= (Vector) evClassMap.get(keys[i].toString());
					if (v== null)
						evClassMap.put(keys[i].toString(), v= new Vector());
					v.add(keys[i]);
					System.out.println("Event "+keys[i]);
					System.out.println(((ASVariation) keys[i]).toStringCoordinates());
					for (int j = 0; j < trptGrpV.size(); j++) {
						System.out.print("\t");
						String[] tIDs= ((VariantGroup_SpliceChain) trptGrpV.elementAt(j)).getTranscriptIDs();
						for (int k = 0; k < tIDs.length; k++) 
							System.out.print(tIDs[k]+" ");						
						System.out.println(" =>  "+tIDhash.get(tIDs[0]));
					} 
				}
				if (print)
					System.out.println();
			}
			
			System.out.println("Double events: "+cntDouble);
			keys= evClassMap.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				System.out.println(((Vector) evClassMap.get(keys[i])).size()+"\t"+keys[i]);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void _05_overlapASEvents(String fName) {
			
		try {
			DualHashBidiMap vTags= getVirtualTagHash(fName);	// tag x reg(tIDs..)
			HashMap tIDhash= new HashMap(vTags.size());	// tID x tag
			Object[] keys= vTags.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				if (!isMPSStag(((String) keys[i]).length()))
					continue;
				DirectedRegion reg= (DirectedRegion) vTags.get(keys[i]);
				Vector v= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
				for (int j = 0; j < v.size(); j++)  {
					tIDhash.put(v.elementAt(j), keys[i]);
				}
			}

			System.out.println("Reading "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			HashMap eventToGroupMap= new HashMap(); // asEv x V(grps)
			int cntOvlEventsSAGE= 0, cntOvlEJSAGE= 0, 
				cntOvlEventsMPSS= 0, cntOvlEJMPSS= 0, 
				cntTotSAGE= 0, cntTotMPSS= 0;
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					ASVariation[] vars= genes[i].getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
					DirectedRegion[] regs= null;
					if (vars!= null) {
						regs= new DirectedRegion[vars.length];
						for (int j = 0; j < regs.length; j++) 
							regs[j]= vars[j].getRegion();
					}
					Transcript[] trpts= genes[i].getTranscripts();
					for (int j = 0; j < trpts.length; j++) {
							// check ovl events
						String tag= (String) tIDhash.get(trpts[j].getTranscriptID());
						if (tag== null)
							continue;
						if (isMPSStag(tag.length()))
							++cntTotMPSS;
						else
							++cntTotSAGE;
						
						DirectedRegion tagReg= (DirectedRegion) vTags.get(tag);
						for (int k = 0; regs!= null&& k < regs.length; k++) 
							if (tagReg.overlaps(regs[k])) {
								if (isMPSStag(tag.length()))
									++cntOvlEventsMPSS;
								else
									++cntOvlEventsSAGE;
								break;
							}
						
							// chk overlap exon junction
						Exon[] exons= trpts[j].getExons();
						for (int k = 1; k < exons.length; k++) {
							if ((k> 0&& tagReg.contains(exons[k].get5PrimeEdge()))||
									(k< exons.length- 2&& tagReg.contains(exons[k].get3PrimeEdge()))) {
								if (isMPSStag(tag.length()))
									++cntOvlEJMPSS;
								else
									++cntOvlEJSAGE;
								break;
							}
						}
					}
				}
				reader.read();
				genes= reader.getGenes();
			}
			
			System.out.println("MPSS total: "+cntTotMPSS+", "+cntOvlEventsMPSS+" overlapped by ASevents and "+cntOvlEJMPSS+" by EJ.");
			System.out.println("SAGE total: "+cntTotSAGE+", "+cntOvlEventsSAGE+" overlapped by ASevents and "+cntOvlEJSAGE+" by EJ.");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static DualHashBidiMap getVirtualTagHash(String fName) {
		if (virtTagHash == null) {
			virtTagHash = new DualHashBidiMap();
			try {
				System.out.println("Reading "+fName);
				GTFChrReader reader= new GTFChrReader(fName);
				reader.setChromosomeWise(false);
				reader.setReadFeatures(new String[] {MPSS_TAG_QUALIFIER, SAGE_TAG_QUALIFIER});
				reader.setReadGene(false);
				reader.setReadGTF(true);
				reader.read();
				GTFObject[] obs= reader.getGtfObj();
				
				int cntMPSS= 0, cntSAGE= 0, cntMPSStot= 0, cntSAGEtot= 0;
				HashMap remMap= new HashMap();				
				for (int i = 0; i < obs.length; i++) {
					if (obs[i].getFeature().equals(MPSS_TAG_QUALIFIER))
						++cntMPSStot;
					else if (obs[i].getFeature().equals(SAGE_TAG_QUALIFIER))
						++cntSAGEtot;					
					String seq= obs[i].getAttribute(GTFObject.ID_ATTRIBUTE_SEQUENCE);
					String seqR= Arrays.reverseComplement(seq);
					DirectedRegion reg= new DirectedRegion(obs[i]);	
					if (virtTagHash.get(seq)!= null|| virtTagHash.get(seqR)!= null) {
						DirectedRegion regTag= (DirectedRegion) virtTagHash.get(seq);
						if (regTag== null)
							regTag= (DirectedRegion) virtTagHash.get(seqR);
						if (regTag.equals(reg)) {
							Vector v= (Vector) regTag.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
							v.add(obs[i].getTranscriptID());
						} else { // not unique
							Integer cnflictCnt= (Integer) remMap.remove(seq);
							if (cnflictCnt== null) 
								cnflictCnt= (Integer) remMap.remove(seqR);
							if (cnflictCnt== null) 
								cnflictCnt= new Integer(2);
							else
								cnflictCnt= new Integer(cnflictCnt.intValue()+ 1);
							remMap.put(seq, cnflictCnt);
						}
						continue;
					}
					Vector v= new Vector();
					v.add(obs[i].getTranscriptID());
					reg.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, v);
					virtTagHash.put(obs[i].getAttribute(GTFObject.ID_ATTRIBUTE_SEQUENCE), reg);
					if (obs[i].getFeature().equals(MPSS_TAG_QUALIFIER))
						++cntMPSS;
					else if (obs[i].getFeature().equals(SAGE_TAG_QUALIFIER))
						++cntSAGE;					
				}
				Object[] keys= remMap.keySet().toArray();		// remove firsts		
				for (int i = 0; i < keys.length; i++) {
					virtTagHash.remove(keys[i]);
					String seqR= Arrays.reverseComplement((String) keys[i]);
					virtTagHash.remove(seqR);
					if (isMPSStag(((String) keys[i]).length()))
						--cntMPSS;
					else
						--cntSAGE;
				}
				
				System.out.println("Read "+ cntMPSStot+ " MPSS tags, "+cntMPSS+" of them unique.");
				System.out.println("Read "+ cntSAGEtot+ " SAGE tags, "+cntSAGE+" of them unique.");

				IntVector distSAGEv= new IntVector(), distMPSSv= new IntVector();
				keys= virtTagHash.keySet().toArray();
				for (int i = 0; i < keys.length; i++) {
					DirectedRegion reg= (DirectedRegion) virtTagHash.get(keys[i]);
					Vector variants= (Vector) reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
					if (isMPSStag(((String) keys[i]).length()))
						distMPSSv.putValue(variants.size(), distMPSSv.getValue(variants.size())+ 1);
					else 
						distSAGEv.putValue(variants.size(), distSAGEv.getValue(variants.size())+ 1);
				}
				System.out.println("Distribution for variants causing same MPSS tag:");
				int[] dist= distMPSSv.toIntArray();
				int sum= 0;
				for (int j = 1; j < dist.length; j++) {
					if (j== 1) {
						System.out.println("("+j+ "\t"+ dist[j]+")");
						continue;
					}
					System.out.println(j+ "\t"+ dist[j]);
					sum+= dist[j];
				}
				System.out.println("total\t"+sum);
				System.out.println("Distribution for variants causing same SAGE tag:");
				dist= distSAGEv.toIntArray();
				sum= 0;
				for (int j = 1; j < dist.length; j++) {
					if (j== 1) {
						System.out.println("("+j+ "\t"+ dist[j]+")");
						continue;
					}
					System.out.println(j+ "\t"+ dist[j]);
					sum+= dist[j];
				}
				System.out.println("total\t"+sum);
				
				System.out.println("Distribution for not unique tags:");
				Object[] vals= remMap.values().toArray();
				java.util.Arrays.sort(vals);
				sum= 0;
				for (int i = 0; i < vals.length; i++) {
					sum+= ((Integer) vals[i]).intValue();
					System.out.print(vals[i]+ " ");
				}
				System.out.println("= "+sum);
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		return virtTagHash;
	}
	
	
	
	public static boolean isMPSStag(int len) {
		if (len== MPSS_TAG_LENGTH)
			return true;
		return false;
	}

	void checkComplete3PSandro() {
		
		// check polyA:
			// at least 8nt A at the end
		
		// check internal priming: 
		// poly-T for mRNA generation hybridizes with A-rich region
		// in last exon
			// sliding window, size 10, 8nt A
		
	}
	
	public static String getDigSiteQualifier(String digSite) {
		if (digSite.equalsIgnoreCase(MPSS_DIGEST_SITE_DPNII))
			return MPSS_TAG_QUALIFIER;
		if (digSite.equalsIgnoreCase(SAGE_DIGEST_SITE_NLAIII))
			return SAGE_TAG_QUALIFIER;
		return null;
	}
}


