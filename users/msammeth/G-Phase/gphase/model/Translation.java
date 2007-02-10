/*
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.NMDSimulator;
import gphase.StopCodons;
import gphase.tools.Arrays;
import gphase.tools.IntVector;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Vector;

/**
 * 
 * 
 * @author micha
 */
public class Translation extends DirectedRegion {

	static final long serialVersionUID= 8996021902187779155L;
	
	public static final String START_CODON= "ATG";
	public static final String[] STOP_CODONS= new String[] {"TAA", "TGA", "TAG"};
	
	String translationID= null;
	Transcript transcript= null;
	int splicedLength= -1;
	int frame= -1;
	
	public static int[] getCodonCount(String[] codons, String seq) {
		int[] res= new int[3];
		for (int i = 0; i < res.length; i++) 
			res[i]= getCodonCount(codons, seq, i);
		return res;
	}
	
	public static int getCodonCount(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		int cnt= 0;
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				++cnt;
		}
		return cnt;
	}

	/**
	 * condenses all 3 frames to one position array
	 * @param codons
	 * @param seq
	 * @return
	 */
	public static int[] getCodonPositions(String[] codons, String seq) {
		IntVector res= new IntVector();
		for (int i = 0; i < 3; i++) { 
			int[] pos= getCodonPositions(codons, seq, i);
			for (int j = 0; j < pos.length; j++) {
				int ins= java.util.Arrays.binarySearch(res.toIntArray(), pos[j]);
				if (ins< 0) 
					res.insert(pos[j], ins);
			}
		}
		return res.toIntArray();
	}
	
	public static int[] getCodonPositions(String[] codons, String seq, int frame) {
		seq= seq.toUpperCase();
		IntVector pos= new IntVector();
		for (int i = frame; i < seq.length()- 2; i+= 3) {
			String cod= seq.substring(i, i+3);
			int j;
			for (j = 0; j < codons.length; j++) 
				if (codons[j].equals(cod))
					break;
			if (j< codons.length)
				pos.add(i);
		}
		return pos.toIntArray();
	}
	
	public static int[] getStartCount(String seq) {
		return getCodonCount(new String[] {START_CODON}, seq);
	}
	
	public static int[] getStopCount(String seq) {
		return getCodonCount(STOP_CODONS, seq);
	}
	
	public Translation(Transcript newTranscript, int newStart, int newEnd, int newStrand) {
		super(newStart, newEnd, newStrand);
		this.transcript= newTranscript;
	}
	
	public Translation(Transcript newTranscript) {
		this.transcript= newTranscript;
		this.strand= getTranscript().getStrand();
		setID("translation");
		setStrand(getTranscript().getGene().getStrand());
	}
	
	public int getSplicedLength() {
		if (splicedLength< 0) {
			DirectedRegion[] regs= transcript.getCDSRegions();
			splicedLength= 0;
			for (int i = 0; i < regs.length; i++) 
				splicedLength+= regs[i].getLength();
		}
		return splicedLength;
	}
	
	public DirectedRegion[] getExonicRegions() {
		Exon[] ex= getTranscript().getExons();
		Vector regV= new Vector(ex.length);
		for (int i = 0; i < ex.length; i++) {
			if (!ex[i].isCoding())
				continue;
			if (ex[i].isCoding5Prime()&& ex[i].isCoding3Prime()) {
				regV.add(ex[i]);
				continue;
			}
				
			DirectedRegion reg= new DirectedRegion(
					ex[i].get5PrimeCDS(), ex[i].get3PrimeCDS(), ex[i].getStrand());
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			regV.add(reg);
		}
		return (DirectedRegion[]) Arrays.toField(regV);
	}
	

	public Translation(Transcript newTranscript, String stableTranslationID) {

		this(newTranscript);
		this.translationID= stableTranslationID;
	}

	public String getChromosome() {
		return getTranscript().getChromosome();
	}
	/**
	 * @return Returns the transcript.
	 */
	public Transcript getTranscript() {
		return transcript;
	}
	
	public Species getSpecies() {
		return getTranscript().getSpecies();
	}
	/**
	 * @return Returns the translationID.
	 */
	public String getTranslationID() {
		return translationID;
	}
	/**
	 * @param translationID The translationID to set.
	 */
	public void setTranslationID(String newTranslationID) {
		this.translationID = newTranslationID;
	}

	public void setSplicedLength(int splicedLength) {
		this.splicedLength = splicedLength;
	}

	/**
	 * Never annotate an ATG starting internal of another CDS > 35 aa upstream
	 * of the ATG as is subject to NMD. [HAVANA]
	 * 
	 * @param trans
	 * @param maxDistAA
	 * @return
	 */
	public Translation[] getUsORF() {
		
		Translation[] trns= getTranscript().getAllORFs();
		Vector uOrfV= new Vector();
		for (int i = 0; trns!= null&& i < trns.length; i++) {
			if (trns[i].get3PrimeEdge()< get5PrimeEdge())
				uOrfV.add(trns[i]);
		}
		
		Translation[] uOrf= (Translation[]) Arrays.toField(uOrfV);
		return uOrf;
	}

	/**
	 * @return genomic pos of premature stops, sorted
	 */
	public int[] getPrematureStops() {
		Exon startEx= getFirstCodingExon();
		String seq= getSplicedSequence();
		int[] pos= getCodonPositions(STOP_CODONS, seq, getFrame(startEx));
		
			// correct to genomic coordinates
		for (int i = 0; i < pos.length; i++) 
			pos[i]= getGenomicPosition(pos[i]+ 2);	// genomic pos of end of stop-codon 
		return pos;
	}

	public int getFrame(Exon ex) {
		int cdsStart= getTranscript().getExonicPosition(get5PrimeEdge());
		int exStart= getTranscript().getExonicPosition(ex.get5PrimeEdge());
		int frame= getFrame();
		if (exStart> cdsStart) {	// cds starts before exon 5'boundary
			frame+= (exStart- cdsStart)% 3;	// frameshift
			if (frame> 2)
				frame%= 3;
		} else
			assert(frame>= 0);	// frame detection must have worked!
		return frame;
	}
	
	/**
	 * @return starting frame of translation
	 */
	public int getFrame() {
		
		if (frame< 0) {
			frame= getFirstCodingExon().getFrame();	// not reliable, also not in Havana
			if (frame< 0) 	// if initialized..
				frame= 0;	// guess 0, good for predicted reading frames
			
			int tlnExStart= getTranscript().getExonicPosition(get5PrimeEdge())+ frame;
			String exSeq= getTranscript().getSplicedSequence();		
				String startCodon= exSeq.substring(tlnExStart, tlnExStart+3);
				if (startCodon.equalsIgnoreCase(START_CODON))
					return frame;	// if there is a start codon here, trust the annotated frame
			
				// otherwise have to guess, look for frames without stop codons
			tlnExStart-= frame;
			int tlnExEnd= getTranscript().getExonicPosition(get3PrimeEdge());
			exSeq= exSeq.substring(tlnExStart, tlnExEnd+ 1);	// max tln seq, evtlly includes term stop
			int[][] stops= new int[3][];
			for (int i = 0; i < stops.length; i++) 
				stops[i]= getCodonPositions(STOP_CODONS, exSeq, i);	// stops in all 3 frames
			
				// and now? take the one with least stops? the one with longest ORF? ???
			int min= Integer.MAX_VALUE;
			IntVector v= new IntVector();
			for (int i = 0; i < stops.length; i++) {
				if (stops[i]== null) {
					if (min> 0)
						v= new IntVector();
					min= 0;
					v.add(i);
					continue;
				}
				if (stops[i].length< min) {
					min= stops[i].length;
					v= new IntVector();
					v.add(i);
				} else if (stops[i].length== min) 
					v.add(i);
			}
			
				// unique one
			if (v.size()== 1) {
				frame= v.get(0);
				return frame;
			}
			
			for (int i = 0; i < v.size(); i++) {
				if (v.get(i)== frame)
					return frame;	// trust frame annotation, if no inframe stop
			}
			
			assert(true);
				// what if (the 2 not as "frame" annotated) frames have equal count of stop codons?		
//			int maxORFFrame= -1;
//			int maxORFEnd= Integer.MIN_VALUE;
//			for (int i = 0; i < v.size(); i++) 
//				if (stops[v.get(i)][0]> maxORFEnd) {	// actually, this should be the longest reading frame.. 
//					maxORFEnd= stops[v.get(i)][0];			// maybe not, we keep for the moment this tln, so assume start here/us
//					maxORFFrame= v.get(i);
//				}
//			return maxORFFrame;
		}
		
		return frame;
	}
	
	public String getSplicedSequence() {
		DirectedRegion[] regs= getExonicRegions();	// not sorted?
		if (regs== null)
			return "";
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) 
			sb.append(Graph.readSequence(regs[i]));
		return sb.toString();
	}
	
	public int getGenomicPosition(int exonPos) {

		Exon[] exons= getTranscript().getExons();
		
			// find containing exon
		int x;
		int dist= 0;
		for (x = 0; dist<= exonPos&& x < exons.length; x++) 
			dist+= getCDSLength(exons[x]);
		if (x> 0) {
			--x;
			dist-= getCDSLength(exons[x]);
		}
		
		int genPos= Math.max(exons[x].get5PrimeEdge(), get5PrimeEdge())+ (exonPos- dist);		
		return genPos;
	}
	
	/**
	 * to also deal with predicted reading frames
	 */
	public int getCDSLength(Exon ex) {
		if (!this.overlaps(ex))
			return 0;
		DirectedRegion reg= this.intersect(ex);
		return reg.getLength();
	}
	

	public Exon getFirstCodingExon() {
		Exon[] ex= getTranscript().getExons();
		int i;
		for (i = 0; i < ex.length; i++) 
			if (ex[i].get3PrimeEdge()> get5PrimeEdge())		// better than exon.isCodin()
				break;
		if (i< ex.length)
			return ex[i];
		return null;
	}
	
	public boolean isOpenEnded5() {
		return (get5PrimeEdge()== transcript.get5PrimeEdge());
	}

	public boolean isOpenEnded3() {
		return (get3PrimeEdge()== transcript.get3PrimeEdge());
	}
	
	public boolean isOpenEnded() {
		// TODO: check seq for ATG/stop
		if (isOpenEnded5()|| isOpenEnded3())
			return true;
		return false;
	}
	
	

	
}
