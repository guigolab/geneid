/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.Constants;
import gphase.NMDSimulator;
import gphase.SpliceSiteConservationComparator;
import gphase.tools.ENCODE;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import com.mysql.jdbc.UpdatableResultSet;
import com.sun.corba.se.spi.legacy.connection.GetEndPointInfoAgainException;
import com.sun.org.apache.xerces.internal.impl.xs.opti.DefaultDocument;

/**
 * 
 * 
 * @author micha
 */
public class Transcript extends DirectedRegion {

	static final long serialVersionUID = 2863324463934791891L;
	int type = Constants.NOINIT;
	int confidence = Constants.NOINIT;
	TU[] tu= null;
	byte nmd= 0;
	Translation predORF= null;
	String HUGO= null;
	
	public static int REGION_5UTR= 1;
	public static int REGION_3UTR= 2;
	public static int REGION_CDS= 3;
	public static int REGION_COMPLETE_TRANSCRIPT= 4;
	
	
	public boolean containsSS(SpliceSite ss) {
		for (int i = 0; i < spliceChain.length; i++) {
			if (ss.getPos()== spliceChain[i].getPos()&&
					ss.isDonor()== spliceChain[i].isDonor())
				return true;
		}
		return false;
	}
	/**
	 * returns genomic position of 
	 * @param pos genomic position, positive or negative
	 * @param offset pos or neg (5', 3')
	 * @return exonic pos, positive
	 */
	public int getExonicPosition(int pos) {
		
		if (pos> 0&& !isForward())
			pos= -pos;
		
			// find containing exon
		int x;
		int dist= 0;
		for (x = 0; x < exons.length; x++) {
			if (exons[x].contains(pos))
				break;
			else
				dist+= exons[x].getLength();
		}
	
		if (x== exons.length)
			return -1;
		
		return dist+ (pos- exons[x].get5PrimeEdge());	// dunno understand why not +1
	}

	public void addTU(TU newTU) {
			// break multiple TUs by #SSs?
		tu= (TU[]) gphase.tools.Arrays.add(tu, newTU);
		//System.currentTimeMillis();
	}
	
	public void removeTU(TU newTU) {
		tu= (TU[]) gphase.tools.Arrays.remove(tu, newTU);
	}
	
	public int getGenomicPosition(int exonPos) {
		
			// find containing exon
		int x;
		int dist= 0;
		for (x = 0; dist<= exonPos&& x < exons.length; x++) 
			dist+= exons[x].getLength();
		if (x> 0) {
			--x;
			dist-= exons[x].getLength();
		}
		
		int genPos= exons[x].get5PrimeEdge()+ (exonPos- dist);		
		return genPos;
	}
	
	/**
	 * 
	 * @param newType type descriptor
	 * @return <code>false</code> if the type is already set or the type
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setType(String newType) {
		
		if (newType== null|| type!= Constants.NOINIT) 
			return false;
		
		int i= Constants.findIgnoreCase(newType, Constants.TYPES);
		if (i< 0) {
			System.err.println("Unknown type "+ newType);
			return false;
		}
		type= i;
		return true;
	}
	
	public String toString() {
		return getTranscriptID();
	}
	
	public String toStringNMD() {
		String s= transcriptID+ "\t";
		if ((translations!= null)&& (!translations[0].isOpenEnded()))
			s+= translations[0].getStart()+ "\t"+ translations[0].getEnd()+ "\t";
		else
			s+= ".\t.\t";
		if (predORF!= null)
			s+= predORF.getStart()+ "\t"+ predORF.getEnd()+ "\t";
		else
			s+= ".\t.\t";
		s+= "NMD"+ nmd;
		return s;
	}
	
	public static String toStringSChain(SpliceSite[] spliceChain) {
		String s= "{";
		for (int i = 0; spliceChain!= null&& i < spliceChain.length; i++) 
			s+= spliceChain[i].getPos()+ ",";
		if (spliceChain!= null&& spliceChain.length> 0)
			s= s.substring(0, s.length()- 1);
		s+= "}";
		return s;
	}

	/**
	 * 
	 * @param newConfidence confidence descriptor
	 * @return <code>false</code> if confidence is already set or the confidence
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setConfidence(String newConfd) {
		
		if (newConfd== null|| type!= Constants.NOINIT)
			return false;
		
		int i= Constants.findIgnoreCase(newConfd, Constants.CONFIDENCES);
	
		if (i< 0) {
			System.err.println("Unknown confidence "+ newConfd);
			return false;
		}
	
		confidence= i;
		return true;
	}
	Translation[] translations= null;
	SpliceSite[] spliceChain= null;	// sorted!!
	public final static Transcript[] toTranscriptArray(Vector v) {
		if (v== null)
			return null;
		Transcript[] result= new Transcript[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Transcript) v.elementAt(i);
		return result;
	}	
	
	/**
	 * 
	 * @return the first translation registered for <code>this</code> transcript.
	 */
	public Translation getDefaultTranslation() {
		if (translations== null) {
			translations= new Translation[]{new Translation(this)};
		}
		return translations[0];
	}
	
	public SpliceSite[] getSpliceChain(){
		
		spliceChain= null;
		if (spliceChain== null) {
			if (exons== null)
				return null;
			
			spliceChain= new SpliceSite[exons.length* 2- 2];
			int pos= 0;
			for (int i = 0; i < exons.length; i++) {
				if(i> 0)
					spliceChain[pos++]= exons[i].getAcceptor();
				if(i< exons.length- 1)
					spliceChain[pos++]= exons[i].getDonor();
			}
			
			Arrays.sort(spliceChain, new AbstractSite.PositionComparator());
		}
		return spliceChain;
	}
	
	public DirectedRegion[] get5UTRRegion(boolean beforeSplicing) {
		Translation[] trans= getTranslations();
		if (trans== null|| trans.length< 1)
			return null;
		if (beforeSplicing) {
			return new DirectedRegion[] {
					new DirectedRegion(
						Math.min(Math.abs(trans[0].get5PrimeEdge()), Math.abs(get5PrimeEdge())),
						Math.max(Math.abs(trans[0].get5PrimeEdge()), Math.abs(get5PrimeEdge())),
						getStrand()
					)};
		} else {
			int atg= trans[0].get5PrimeEdge();
			Vector v= new Vector();
			for (int i = 0; i < exons.length; i++) {
				if (exons[i].get3PrimeEdge()>= atg) {	// atg containing exon
					int min= exons[i].get5PrimeEdge();
					int max= atg- 1;
					if (!isForward()) {
						int h= min;
						min= max;
						max= h;
					}
					DirectedRegion dir= null;
					if (Math.abs(max)- Math.abs(min)>= 0) {	// non-empty region
						dir= new DirectedRegion(min, max, getStrand());
						dir.setChromosome(getChromosome());
						v.add(dir);
					}
					break;
				}
				DirectedRegion dir= new DirectedRegion(exons[i].getStart(), exons[i].getEnd(), exons[i].getStrand());
				dir.setChromosome(getChromosome());
				v.add(dir);
			}
			return (DirectedRegion[]) gphase.tools.Arrays.toField(v);
		}
	}

	public int get5UTR(boolean beforeSplicing) {
		Translation[] trans= getTranslations();
		if (trans== null|| trans.length< 1)
			return -1;
		if (beforeSplicing) {
			return (trans[0].get5PrimeEdge()- get5PrimeEdge()+ 1);
		} else {
			int atg= trans[0].get5PrimeEdge();
			int sum= 0;
			for (int i = 0; i < exons.length; i++) {
				if (exons[i].get3PrimeEdge()> atg) {	// atg containing exon
					sum+= atg- exons[i].get5PrimeEdge()+ 1;	// add rest
					break;
				}
				sum+= exons[i].getLength();
			}
			return sum;
		}
	}
	
	public DirectedRegion[] getExonicRegions() {
		DirectedRegion[] regs= new DirectedRegion[exons.length];
		for (int i = 0; i < regs.length; i++) {
			regs[i]= new DirectedRegion(exons[i].getStart(), exons[i].getEnd(), exons[i].getStrand());
			regs[i].setChromosome(getChromosome());
			regs[i].setSpecies(getSpecies());
		}
		return regs;
	}
	
	public String getSplicedSequence() {
		DirectedRegion[] regs= getExonicRegions();	// not sorted
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) 
			sb.append(Graph.readSequence(regs[i]));
		return sb.toString();
	}
		
	public int getCDSLength(boolean spliced) {
		if (translations== null|| translations.length< 1)
			return -1;
		if (spliced) {
			DirectedRegion[] reg= getCDSRegions();
			int sum= 0;
			for (int i = 0; i < reg.length; i++) 
				sum+= reg[i].getLength();
			return sum;
		} else {
			return translations[0].get3PrimeEdge()- translations[0].get5PrimeEdge()+ 1;
		}
	}
	
	public int get5UTRLength(boolean spliced) {
		if (translations== null|| translations.length< 1)
			return -1;
		if (spliced) {
			DirectedRegion[] reg= get5UTRRegion(false);
			int sum= 0;
			for (int i = 0; reg!= null&& i < reg.length; i++) 
				sum+= reg[i].getLength();
			return sum;
		} else {
			return translations[0].get5PrimeEdge()- getStart();
		}
	}
	
	/**
	 * 	working with translation
	 * @return
	 */
	public DirectedRegion[] getCDSRegions_old() {
		
		int tlStart= getTranslations()[0].get5PrimeEdge();
		int tlEnd= getTranslations()[0].get3PrimeEdge();
		Vector v= new Vector();
		int i;
		
			// first region;
		for (i = 0; i < exons.length; i++) 
			if (exons[i].contains(tlStart)) {
				int min= tlStart;
				int max= exons[i].get3PrimeEdge();
				if (Math.abs(min)> Math.abs(max)) {
					int h= min;
					min= max;
					max= h;
				}
				DirectedRegion reg= new DirectedRegion(
						min, 
						max, 
						exons[i].getStrand());
				reg.setChromosome(getChromosome());
				v.add(reg);
				break;
			}
		
		for (; i < exons.length; i++) {
			
				// last exon
			if (exons[i].contains(tlEnd)) {
				int min= exons[i].get5PrimeEdge();
				int max= tlEnd;
				if (Math.abs(min)> Math.abs(max)) {
					int h= min;
					min= max;
					max= h;
				}
				DirectedRegion reg= new DirectedRegion(
						min, 
						max, 
						exons[i].getStrand());
				reg.setChromosome(getChromosome());
				v.add(reg);
				break;
			}
			
				// intervening exons
			DirectedRegion reg= new DirectedRegion(
					exons[i].getStart(), 
					exons[i].getEnd(), 
					exons[i].getStrand());
			reg.setChromosome(getChromosome());
			v.add(reg);
		}
			
		return (DirectedRegion[]) gphase.tools.Arrays.toField(v);
	}
	
	public String get5UTRSequence() {
		DirectedRegion[] reg= get5UTRRegion(false);
		
		String region= "";
		int min= 0, max= 0;
		for (int j = 0; reg!= null&& j < reg.length; j++) {
			if (j== 0)
				min= Math.abs(reg[j].get5PrimeEdge());
			if (j== reg.length- 1)
				max= Math.abs(reg[j].get3PrimeEdge());
			DirectedRegion regSav= (DirectedRegion) reg[j].clone();
			reg[j]= ENCODE.convertToEncodeCoord(reg[j]);
			if (reg[j]== null) {
				System.err.println("outside encode, skipping transcript "+getTranscriptID());
				return null;
			}
			String s= SpliceSiteConservationComparator.getSubstring( 
					reg[j].getID(),
					"human", 
					Math.abs(reg[j].getStart()), 
					Math.abs(reg[j].getEnd())+ 1);	// end excl
			if (!reg[j].isForward())
				s= gphase.tools.Arrays.reverseComplement(s);
			region+= s;
		}
		
		return region;
		
	}

	public String getSequence(int regCode) {
		DirectedRegion[] reg= null;
		if (regCode== REGION_5UTR)
			reg= get5UTRRegion(false);
		else if (regCode== REGION_COMPLETE_TRANSCRIPT)
			reg= getExonicRegions();
		else if (regCode== REGION_CDS)
			reg= getCDSRegions();
		
		String region= "";
		int min= 0, max= 0;
		for (int j = 0; reg!= null&& j < reg.length; j++) {
			if (j== 0)
				min= Math.abs(reg[j].get5PrimeEdge());
			if (j== reg.length- 1)
				max= Math.abs(reg[j].get3PrimeEdge());
			DirectedRegion regSav= (DirectedRegion) reg[j].clone();
			reg[j]= ENCODE.convertToEncodeCoord(reg[j]);
			if (reg[j]== null) {
				System.err.println("outside encode, skipping transcript "+getTranscriptID());	// skip fragments
				return null;	// continue?
			}
			String s= SpliceSiteConservationComparator.getSubstring( 
					reg[j].getID(),
					"human", 
					Math.abs(reg[j].getStart()), 
					Math.abs(reg[j].getEnd())+ 1);	// end excl
			if (!reg[j].isForward())
				s= gphase.tools.Arrays.reverseComplement(s);
			region+= s;
		}
		
		return region;
		
	}
	

	
	public SpliceSite getPredSpliceSite(SpliceSite s) {
		
		SpliceSite[] ss= getSpliceChain();
		int p= Arrays.binarySearch(ss, s, new SpliceSite.PositionComparator());
		if ((p< 0)|| (p== 0)) 
			return null;
		return ss[p-1];
	}
	
	public SpliceSite getSuccSpliceSite(SpliceSite s) {
		
		SpliceSite[] ss= getSpliceChain();
		int p= Arrays.binarySearch(ss, s, new SpliceSite.PositionComparator());
		if ((p< 0)|| (p> (ss.length-2))) 
			return null;
		return ss[p+1];
	}
	
	public SpliceSite getSuccSpliceSite(int pos) {
		SpliceSite[] ss= getSpliceChain();
		return getSuccSpliceSite(ss, pos);
	}
	public static SpliceSite getSuccSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos+1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p].getPos()> pos)
			return ss[p];
		else
			return null;
	}
	
	public static int getSuccPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos+1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1);	// insertionpoint= point of successor
		if (p< ss.length&& ss[p]> pos)
			return ss[p];
		else
			return 0;
	}
	
	public SpliceSite getPredSpliceSite(int pos) {
		SpliceSite[] ss= getSpliceChain();
		return getPredSpliceSite(ss, pos);
	}
	
	public static SpliceSite getPredSpliceSite(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos-1), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) {	// found ?!
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p].getPos()< pos)
			return ss[p];
		else
			return null;
	}
	
	public SpliceSite[] getLastUTRIntron() {
		if (translations== null|| exons.length< 3|| exons[1].get3PrimeCDS()>= translations[0].get5PrimeEdge())
			return null;
		int i;
		for (i = 0; i < exons.length; i++) 
			if (exons[i].get3PrimeEdge()>= translations[0].get5PrimeEdge())
				break;
		SpliceSite[] ss= new SpliceSite[2];
		ss[0]= exons[i-1].getDonor();
		ss[1]= exons[i].getAcceptor();
		return ss;
	}
	
	public static int getPredPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos-1);	// abstract site for only comparing pos
		if(p>= 0) {
			return ss[p];
		}
		p= -(p+ 1)- 1;	// insertion point- 1= point of predecessor
		if (p>= 0&& p< ss.length&& ss[p]< pos)
			return ss[p];
		else
			return 0;
	}
	
	public static SpliceSite getSpliceSiteByPos(SpliceSite[] ss, int pos) {
		int p= Arrays.binarySearch(ss, new Integer(pos), new AbstractSite.PositionToSpliceSiteComparator());	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return null;
	}
	
	public static int getPos(int[] ss, int pos) {
		int p= Arrays.binarySearch(ss, pos);	// abstract site for only comparing pos
		if(p>= 0) 
			return ss[p];
		return -1;
	}
	
	public boolean addTranslation(Translation newTrans) {
		
			// new transcipt array
		if (translations== null) {
			translations= new Translation[] {newTrans};
			return true;
		}
			
			// search transcript for same translation, not necessary
			
			// add translation
		Translation[] nTranslations= new Translation[translations.length+ 1];
		for (int i= 0; i < translations.length; i++) 
			nTranslations[i]= translations[i];
		nTranslations[nTranslations.length- 1]= newTrans;
		translations= nTranslations;
		return true;
	}
	
	public Translation[] getTranslation(int start, int end) {
		if (translations== null)
			return null;
		Vector result= new Vector();
		for (int i = 0; i < translations.length; i++) {
			if ((translations[i].get5PrimeEdge()<= start)&&
				(translations[i].get3PrimeEdge()>= end))
				result.add(translations[i]);
		}
		
		Object o= gphase.tools.Arrays.toField(result);
		if (o== null|| result.size()< 1)
			return null;
		return (Translation[]) o;
	}
	
	public Translation[] getTranslations() {
		return translations;
	}
	
	public boolean isNonCoding() {
		return (translations== null|| translations.length== 0);
	}
	
	public boolean is5UTR(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) {
			if (translations[i].get5PrimeEdge()> pos)
				break;
		}
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}
	
	public boolean is3UTR(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) {
			if (translations[i].get3PrimeEdge()< pos)
				break;
		}
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}

	/**
	 * when 
	 * @param pos
	 * @return
	 */
	public boolean isCDS(int pos) {
		
		if (translations== null)
			return false;
		
		int i;
		for (i = 0; translations!= null&& i < translations.length; i++) 
			if (pos>= translations[i].get5PrimeEdge()&& pos<= translations[i].get3PrimeEdge())
				break;
		
		if (translations!= null&& i< translations.length)
			return true;
		return false;
	}
	
	public boolean isExonic(int pos) {
		for (int i = 0; exons!= null&& i < exons.length; i++) 
			if (pos>= exons[i].get5PrimeEdge()&& pos<= exons[i].get3PrimeEdge())
				return true;
		return false;
	}
	
	public boolean isCoding() {
		return (translations!= null&& translations.length> 0);
	}
	
	public boolean isForward() {
		if (gene!= null)
			return gene.isForward();
		return super.isForward();
	}
	
	Gene gene= null;

	String transcriptID= null;
	Exon[] exons= null;	// sorted !!
	public Transcript(Gene newGene, String stableTranscriptID) {

		this.strand= newGene.getStrand();
		this.gene= newGene;
		this.transcriptID= stableTranscriptID;
		setID("transcript");
	}
	
	public Transcript(String newID) {
		this.transcriptID= newID;
		setID("transcript");
	}
	
	/**
	 * @return
	 */
	public Exon[] getExons() {
		return exons;
	}
	
	public DirectedRegion[] getIntrons() {
		if (exons== null|| exons.length< 2)
			return null;
		DirectedRegion[] introns= new DirectedRegion[exons.length- 1];
		for (int i = 0; i < introns.length; i++) {
			int pos1= exons[i].get3PrimeEdge()+ 1;
			int pos2= exons[i+1].get5PrimeEdge()- 1;
			if (isForward())
				introns[i]= new DirectedRegion(pos1, pos2, getStrand());
			else
				introns[i]= new DirectedRegion(pos2, pos1, getStrand());
		}
		return introns;
	}
	
	/**
	 * @return
	 */
	public Gene getGene() {
		return gene;
	}

	/**
	 * @return
	 */
	public String getTranscriptID() {
		return transcriptID;
	}
	
	/**
	 * @return
	 */
	public String getStableID() {
		
		return transcriptID;
	}	

	/**
	 * @param exons
	 */
	public void setExons(Exon[] exons) {
		this.exons= exons;
	}

	/**
	 * @deprecated mandatory in constructor now
	 * @param gene
	 */
	public void setGene(Gene gene) {
		this.gene= gene;
		setStrand(getGene().getStrand());
	}

	/**
	 * @deprecated not valid for HOX-genes, AY... genes, ..
	 * @param i
	 */
	public void setTranscriptID(String i) {
		transcriptID= i;
	}

	public String getChromosome() {
		return getGene().getChromosome();
	}
	
	public Species getSpecies() {
		return getGene().getSpecies();
	}
	
	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== getGene().isForward());
	}
	
	public boolean checkStrand(String newStrand) {
		
		String nStrand= newStrand.trim();
		if (nStrand.equals("1")) 	// || nStrand.equals("forward")
			return checkStrand(true); 
		else if (nStrand.equals("-1"))	// || nStrand.equals("reverse")
			return checkStrand(false);
		
		return false; // error			
	}


	public void addCDS(int start, int end) {
		if (!isForward()) {
			start= -start;
			end= -end;
		}
			
		if (translations== null) {
			translations= new Translation[] {new Translation(this)};
			translations[0].setStart(start);
			translations[0].setEnd(end);
			translations[0].setChromosome(getChromosome());
			translations[0].setSpecies(getSpecies());
			return;
		}
			// else
		if (Math.abs(start)< Math.abs(translations[0].getStart()))
			translations[0].setStart(start);
		if (Math.abs(end)> Math.abs(translations[0].getEnd()))
			translations[0].setEnd(end);
	}
	
	/**
		 * Inserts the exons in an array sorted according to ascending order
		 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
		 * transcripts to ensure the correct init of AS types.
		 * 
		 * @param newExon
		 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
		 */
		public boolean addExon_new(Exon newExon) {
	
				// new exon array
			if (exons== null) 
				exons= new Exon[] {newExon};
			else {
				
					// search for identical exon (HERE necessary)
				int p= Arrays.binarySearch(
						exons,
						newExon,
						new AbstractRegion.PositionComparator()	
					);
		
				if (p>= 0) 
					return false;	// already contained, not added
				
					// new Exon: search for overlapping exons 
					// and accordingly construct SuperExon
	//			Exon[] exs= getGene().getExons();		// in ALL transcripts
	//			for (int i = 0; i < exons.length; i++) {
	//				if ((exs[i].getStart()>= newExon.getStart()&& exs[i].getStart()< newExon.getEnd())||
	//						newExon.getStart()>= exs[i].getStart()&& newExon.getStart()< exs[i].getEnd()) { // intersecting
	//					if (exs[i].getSuperExon()!= null)
	//						if (newExon.getSuperExon()!= null)
	//							newExon.getSuperExon().merge(exs[i].getSuperExon()); 	// merge two super-exs
	//						else
	//							exs[i].getSuperExon().add(newExon);	// add to super-exon of the other
	//					else {
	//						SuperExon superEx= newExon.getSuperExon();
	//						if (superEx== null) {
	//							superEx= new SuperExon();
	//							superEx.add(newExon);
	//						}
	//						superEx.add(exs[i]);
	//					}
	//				}
	//			}
				
					// add exon
				exons= (Exon[]) gphase.tools.Arrays.insert(this.exons, newExon, p);
			}
	
				// init splice sites 
			int p= Arrays.binarySearch(
					this.exons,
					newExon,
					new AbstractRegion.PositionComparator()	
			);
			
			if ((p== 0)&& exons.length>1) {	// ex-first exon now has an acceptor
				Exon ex= this.exons[1];
				int posAcceptor= isForward()?ex.getStart()-2:ex.getEnd()+1;
				SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false);
				SpliceSite ss= getGene().checkSpliceSite(acceptor);
				if (ss== null) 
					getGene().addSpliceSite(acceptor);
				else
					acceptor= ss;
				ex.setAcceptor(acceptor);
			}
			
			if (p> 0) {										// has acceptor
				int posAcceptor= isForward()?newExon.getStart()-2:newExon.getEnd()+1;
				SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false);
				SpliceSite ss= getGene().checkSpliceSite(acceptor);
				if (ss== null) 
					getGene().addSpliceSite(acceptor);
				else
					acceptor= ss;
				exons[p].setAcceptor(acceptor);
			}
			
			if (p< this.exons.length- 1) {					// has donor
				int posDonor= isForward()?newExon.getEnd()+1:newExon.getStart()-2;
				SpliceSite donor= new SpliceSite(getGene(), posDonor, true);
				SpliceSite ss= getGene().checkSpliceSite(donor);
				if (ss== null) 
					getGene().addSpliceSite(donor);
				else
					donor= ss;
				exons[p].setDonor(donor);
			}
	
			if ((p== this.exons.length- 1)&& exons.length>1) {	// ex-last exon now has an donor
				Exon ex= this.exons[exons.length- 2];
				int posDonor= isForward()?ex.getEnd()+1:ex.getStart()-2;
				SpliceSite donor= new SpliceSite(getGene(), posDonor, true);
				SpliceSite ss= getGene().checkSpliceSite(donor);
				if (ss== null) 
					getGene().addSpliceSite(donor);
				else
					donor= ss;
				ex.setDonor(donor);
			}
			
			return true;
		}

			/**
			 * Inserts the exons in an array sorted according to ascending order
			 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
			 * transcripts to ensure the correct init of AS types.
			 * 
			 * @param newExon
			 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
			 */
			public boolean addExon(Exon newExon) {
		
					// new exon array
				if (exons== null) 
					exons= new Exon[] {newExon};
				else {
					
						// search for identical exon (HERE necessary)
					int p= Arrays.binarySearch(
							exons,
							newExon,
							new AbstractRegion.PositionComparator()	// has to be directed, end< start for two ident. regions	
						);
			
					if (p>= 0) 
						return false;	// already contained, not added
					
						// new Exon: search for overlapping exons 
						// and accordingly construct SuperExon
		//			Exon[] exs= getGene().getExons();		// in ALL transcripts
		//			for (int i = 0; i < exons.length; i++) {
		//				if ((exs[i].getStart()>= newExon.getStart()&& exs[i].getStart()< newExon.getEnd())||
		//						newExon.getStart()>= exs[i].getStart()&& newExon.getStart()< exs[i].getEnd()) { // intersecting
		//					if (exs[i].getSuperExon()!= null)
		//						if (newExon.getSuperExon()!= null)
		//							newExon.getSuperExon().merge(exs[i].getSuperExon()); 	// merge two super-exs
		//						else
		//							exs[i].getSuperExon().add(newExon);	// add to super-exon of the other
		//					else {
		//						SuperExon superEx= newExon.getSuperExon();
		//						if (superEx== null) {
		//							superEx= new SuperExon();
		//							superEx.add(newExon);
		//						}
		//						superEx.add(exs[i]);
		//					}
		//				}
		//			}
					
						// add exon
					exons= (Exon[]) gphase.tools.Arrays.insert(this.exons, newExon, p);
				}
		
					// init splice sites 
				int p= Arrays.binarySearch(
						this.exons,
						newExon,
						new AbstractRegion.PositionComparator()		// here also abstractregion possible	
				);
				
				if (p== 0) {
					//int tss= isForward()?newExon.getStart():newExon.getEnd();
					if (isForward())
						setStart(newExon.getStart());
					else
						setEnd(newExon.getEnd());
					if (exons.length>1) {	// ex-first exon now has an acceptor
						Exon ex= this.exons[1];
						int posAcceptor= isForward()?ex.getStart():ex.getEnd();
						SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false, ex);
						SpliceSite ss= getGene().checkSpliceSite(acceptor);
						if (ss== null) 
							getGene().addSpliceSite(acceptor);
						else {
							acceptor= ss;
							ss.addExon(ex);
						}
						ex.setAcceptor(acceptor);
					}
				}
				
				if (p> 0) {										// has acceptor
					int posAcceptor= isForward()?newExon.getStart():newExon.getEnd();
					SpliceSite acceptor= new SpliceSite(getGene(), posAcceptor, false, newExon);
					SpliceSite ss= getGene().checkSpliceSite(acceptor);
					if (ss== null) 
						getGene().addSpliceSite(acceptor);
					else {
						acceptor= ss;
						ss.addExon(newExon);
					}
					exons[p].setAcceptor(acceptor);
				}
				
				if (p< this.exons.length- 1) {					// has donor
					int posDonor= isForward()?newExon.getEnd():newExon.getStart();
					SpliceSite donor= new SpliceSite(getGene(), posDonor, true, newExon);
					SpliceSite ss= getGene().checkSpliceSite(donor);
					if (ss== null) 
						getGene().addSpliceSite(donor);
					else {
						donor= ss;
						ss.addExon(newExon);
					}
					exons[p].setDonor(donor);
				}
		
				if (p== this.exons.length- 1) {
					//int tes= isForward()?newExon.getEnd():newExon.getStart();
					if (isForward())
						setEnd(newExon.getEnd());
					else
						setStart(newExon.getStart());
					if (exons.length>1) {	// ex-last exon now has an donor
						Exon ex= this.exons[exons.length- 2];
						int posDonor= isForward()?ex.getEnd():ex.getStart();
						SpliceSite donor= new SpliceSite(getGene(), posDonor, true, ex);
						SpliceSite ss= getGene().checkSpliceSite(donor);
						if (ss== null) 
							getGene().addSpliceSite(donor);
						else {
							donor= ss;
							ss.addExon(ex);
						}
						ex.setDonor(donor);
					}
				}
				
				int sstart= start;
				int send= end;
				setBoundaries(newExon);
				if (start!= sstart|| send!= end)
					System.currentTimeMillis();
				return true;
			}

		/**
		 * Inserts the exons in an array sorted according to ascending order
		 * of their start/stop position. <b>IMPORTANT</b>: add exons AFTER adding 
		 * transcripts to ensure the correct init of AS types.
		 * 
		 * @param newExon
		 * @return the exon already contained or <code>newExon</code> case of the exon was added successfully
		 */
		public void setBoundaries(Exon newExon) {
	
			if (Math.abs(newExon.getStart())< Math.abs(getStart()))
				setStart(newExon.getStart());
			if (Math.abs(newExon.getEnd())> Math.abs(getEnd()))
				setEnd(newExon.getEnd());
//			if (newExon.getStart()< getStart())
//				setStart(newExon.getStart());
//			if (newExon.getEnd()> getEnd())
//				setEnd(newExon.getEnd());
		}
		
		public void setStart(int v) {
			super.setStart(v);
		}
		
		public void setEnd(int v) {
			super.setEnd(v);
		}

	/**
	 * Finds the first exon containing the corresponding position
	 * 
	 * @deprecated inconsistent for overlapping exons		 
	 * @param absPos
	 * @return
	 */
	public Exon getExon(int absPos) {
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].contains(absPos))
				return exons[i];
		
		return null;
	}

	public int getOtherSideOfExon(int pos) {
		for (int i = 0; i < exons.length; i++) {
			if (exons[i].getStart()== pos)
				return exons[i].getEnd();
			else if (exons[i].getEnd()== pos)
				return exons[i].getStart();
		}
		
		return -1;
	}
	
	public Translation[][] findORFs(boolean openEnded) {
		String seq= getSplicedSequence();
		Translation[][] orfs= new Translation[3][];
		for (int i = 0; i < 3; i++) {
			orfs[i]= findORFs(seq, i, openEnded);
		}
		// skipped for the moment
//		int cntORFs= 0;
//		for (int i = 0; i < orfs.length; i++) 
//			if (orfs[i]!= null)
//				cntORFs+= orfs[i].length;	// TODO min threshold
//		if (cntORFs== 0) {
//			for (int i = 0; i < orfs.length; i++) {
//				orfs[i]= forceORFs(seq, i);
//			}
//		}
		return orfs;
	}
	
	/**
	 * Use the most upstream ATG as start codon. If there is no stop codon 
	 * upstream of the first ATG, check genomic sequence upstream beyond the 
	 * gene for in-frame stop. If one is found and the start of the gene is
	 * in a CpG island, use the ATG, otherwise start the open-ended CDS at the
	 * start of the gene.
	 * @return
	 */
	public Translation findHavanaORF() {
		
		Translation[][] tlns= findORFs(false);
		if (tlns== null)
			return null;
		
			// find annotated start atgs in the gene
		Transcript[] trpts= getGene().getTranscripts();
		Vector atgPosV= new Vector();
		for (int i = 0; i < trpts.length; i++) {
			if (trpts[i]== this)
				continue;
			if (trpts[i].isCoding()) {
				Translation tmpTln= trpts[i].getTranslations()[0];
				// now, predict open ended, if already reported
				if (tmpTln.isOpenEnded())	
					continue;
				Integer pos= new Integer(tmpTln.get5PrimeEdge());
				int j; 
				for (j = 0; j < atgPosV.size(); j++) 
					if (atgPosV.elementAt(j).equals(pos))
						break;
				if (j== atgPosV.size())
					atgPosV.add(pos);
			}
		}
		
			// find ORFs at the corresponding starts
		Vector tlnV= new Vector();
		for (int i = 0; i < tlns.length; i++) {
			for (int j = 0; tlns[i]!= null&& j < tlns[i].length; j++) {
				int k;
				for (k = 0; k < atgPosV.size(); k++) 
					if (tlns[i][j].get5PrimeEdge()== ((Integer) atgPosV.elementAt(k)).intValue())
						break;
				if (k< atgPosV.size())
					tlnV.add(tlns[i][j]);
			}
		}

			// if coincidence with annotated ORF, get longest of this
		if (tlnV.size()< 1)
			return findLongestORF(tlns);
		else
			return findLongestORF(new Translation[][] {(Translation[]) gphase.tools.Arrays.toField(tlnV)});
	}
	
	public Translation findLongestORF() {
		return findLongestORF(findORFs(false));
	}
	
	Translation findLongestORF(Translation[][] predORFs) {
		Translation longestORF= null;
		if (predORFs== null)
			return null;
		for (int i = 0; i < predORFs.length; i++) 
			for (int j = 0; predORFs[i]!= null&& j < predORFs[i].length; j++) 
				if (longestORF== null|| predORFs[i][j].getSplicedLength()> longestORF.getSplicedLength())
					longestORF= predORFs[i][j];
		
		if (longestORF!= null&& longestORF.getSplicedLength()< (NMDSimulator.MIN_ORF_LENGTH_AA* 3))
			return null;
		return longestORF;
	}
	
	public Translation[] getAllORFs() {
		Translation[][] predORFs= findORFs(false);
		if (predORFs== null)
			return null;
		Vector allOrfV= new Vector();
		for (int i = 0; i < predORFs.length; i++) 
			for (int j = 0; predORFs[i]!= null&& j < predORFs[i].length; j++) 
				allOrfV.add(predORFs[i][j]);
		
		Translation[] allORFs= (Translation[]) gphase.tools.Arrays.toField(allOrfV);
		return allORFs;
	}
	
	public Translation[] findORFs(String seq, int frame, boolean openEnded) {
			final int minORFlen= NMDSimulator.MIN_ORF_LENGTH_AA* 3;
			seq= seq.substring(frame).toUpperCase();
			Vector startV= new Vector(), stopV= new Vector();
			int pos= 0;
			for (int i = 0; (i+3) < seq.length(); i+= 3) {
				String codon= seq.substring(i, i+3);
				if (codon.equals(Translation.START_CODON))
					startV.add(new Integer(i));
				else if (codon.equals(Translation.STOP_CODONS[0])||
						codon.equals(Translation.STOP_CODONS[1])||
						codon.equals(Translation.STOP_CODONS[2]))
					stopV.add(new Integer(i));
			}
			int[] startPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(startV));
			int[] stopPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(stopV));
	
				// determine ORFs
			Vector v= new Vector();
			Translation reg;
			for (int i = 0; startPos!= null&& stopPos!= null&& i < startPos.length; i++) {
				int j= Arrays.binarySearch(stopPos, startPos[i]);
				if (j>= 0)
					System.err.println("Assertion failed: start/stop at same pos!");
				j= -j- 1;
				if (j== stopPos.length) {
					reg= new Translation(this, getGenomicPosition(frame+ startPos[i]), get3PrimeEdge(), getStrand());
					reg.setSplicedLength(seq.length()- startPos[i]+ 1);
				} else {
					reg= new Translation(this, getGenomicPosition(frame+ startPos[i]), getGenomicPosition(frame+ stopPos[j]+ 2), getStrand());
					reg.setSplicedLength(stopPos[j]+ 2- startPos[i]+ 1);
				}
				reg.setChromosome(getChromosome());
				reg.setSpecies(getSpecies());
				//if (reg.getLength()> minORFlen)
					v.add(reg);
			}
			
			if (!openEnded)
				return (Translation[]) gphase.tools.Arrays.toField(v);

				// 3'OpenEnded
			if (stopPos!= null&& stopPos.length> 0) {	// ORF starting outside transcript
				reg= new Translation(this, get5PrimeEdge(), getGenomicPosition(frame+ stopPos[0]+ 2), getStrand());	// +2 to include stop_codon
				reg.setSplicedLength(stopPos[0]+ 2+ frame+ 1);	// here, q si frame ... goes outside transcript..
				reg.setChromosome(getChromosome());
				reg.setSpecies(getSpecies());
				v.add(reg);
			} else {	// bi-OpenEnded
				reg= new Translation(this, get5PrimeEdge(), get3PrimeEdge(), getStrand());	// +2 to include stop_codon
				reg.setSplicedLength(seq.length());	// here, q si frame ... goes outside transcript..
				reg.setChromosome(getChromosome());
				reg.setSpecies(getSpecies());
				v.add(reg);
			}
				// 5'OpenEnded
			if (startPos!= null&& startPos.length> 0){	
				reg= new Translation(this, getGenomicPosition(frame+ startPos[0]), get3PrimeEdge(), getStrand());
				reg.setSplicedLength(seq.length()- startPos[0]+ 1);
				reg.setChromosome(getChromosome());
				reg.setSpecies(getSpecies());
				v.add(reg);
			} 
			
			return (Translation[]) gphase.tools.Arrays.toField(v);
		}
	
	/**
	 * 
	 * @param seq
	 * @param frame
	 * @deprecated not in use, check inframe stop
	 * @return
	 */
	public Translation[] forceORFs(String seq, int frame) {
		seq= seq.substring(frame).toUpperCase();
		Vector startV= new Vector(), stopV= new Vector();
		for (int i = 0; (i+3) < seq.length(); i+= 3) {
			String codon= seq.substring(i, i+3);
			if (codon.equals(Translation.START_CODON))
				startV.add(new Integer(i));
			else if (codon.equals(Translation.STOP_CODONS[0])||
					codon.equals(Translation.STOP_CODONS[1])||
					codon.equals(Translation.STOP_CODONS[2]))
				stopV.add(new Integer(i));
		}
		int[] startPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(startV));
		int[] stopPos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(stopV));

			// determine ORFs
		Vector v= new Vector();
		Translation reg;
		if (stopPos!= null&& stopPos.length> 0) {	// ORF starting outside transcript
			reg= new Translation(this, get5PrimeEdge(), getGenomicPosition(frame+ stopPos[0]+ 2), getStrand());	// +2 to include stop_codon
			reg.setSplicedLength(stopPos[0]+ 2+ frame+ 1);	// here, q si frame ... goes outside transcript..
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			
			v.add(reg);
		} 
		
		if (startPos!= null&& startPos.length> 0){	
			reg= new Translation(this, getGenomicPosition(frame+ startPos[0]), get3PrimeEdge(), getStrand());
				// TODO check for no inframe stop
			reg.setSplicedLength(seq.length()- startPos[0]+ 1);
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			v.add(reg);
		}
		
		return (Translation[]) gphase.tools.Arrays.toField(v);
	}
	

	
	public int getTSSPos() {
		return exons[0].get5PrimeEdge();
	}

	public AbstractSite getTSS() {
		return new TSSite(this, getTSSPos(), true);
	}
	
	public AbstractSite getTES() {
		return new TSSite(this, getTESPos(), false);
	}
	
	public int getTESPos() {
		return exons[exons.length- 1].get3PrimeEdge();
	}
	
	public Exon getExon(String stableID) {
		
		if (exons== null)
			return null;
		
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].getExonID().equals(stableID))
				return exons[i];
	
		return null;
	}

	public Exon getLastExon() {
		
		if (exons== null|| exons.length< 1)
			return null;
		return exons[exons.length- 1];
	}

	public SpliceSite getSpliceSite(int pos) {
		
		Comparator compi= new SpliceSite.PositionComparator();
		SpliceSite ss= new SpliceSite(null, pos, true);
		int p= Arrays.binarySearch(spliceChain, ss, compi);
		if (p>= 0)
			return spliceChain[p];
		return null;
	}

	public int getDistFromATG(int gPos) {

		if (getTranslations()== null|| getTranslations().length!= 1)
			return -1;
		
		int atg= getTranslations()[0].get5PrimeEdge();
		int min= Math.min(gPos, atg);
		int max= Math.max(gPos, atg);
		
		int i;
		for (i = 0; i < exons.length; i++) 
			if(exons[i].get3PrimeEdge()>= min)
				break;
		if (i== exons.length)	// outa range
			return -1;

		int j;
		for (j = i; j < exons.length; j++) 
			if(exons[j].get3PrimeEdge()>= max)
				break;
		if (j== exons.length)	// outa range
			return -1;
		
		int dist= 0;
		for (int k = i+1; k < j; k++) 
			dist+= exons[k].getLength();
		
		if (i== j)
			dist+= max- min;
		else {
			dist+= exons[i].get3PrimeEdge()- min;	// 1st exon, not +-1
			dist+= max- exons[j].get5PrimeEdge();	// last exon, not +-1
			++dist;	// and in both cases +1
		}

		return dist;
	}
	public TU[] getTu() {
		return tu;
	}

	/**
	 * working with cds in exons
	 * @return
	 */
	public DirectedRegion[] getCDSRegions() {
		Exon[] exons= getExons();
		Vector v= new Vector(exons.length);
		DirectedRegion reg;
		for (int i = 0; i < exons.length; i++) {
			if (!exons[i].isCoding())
				continue;
			reg= new DirectedRegion(exons[i].get5PrimeCDS(), exons[i].get3PrimeCDS(), exons[i].getStrand());
			reg.setChromosome(getChromosome());
			reg.setSpecies(getSpecies());
			v.add(reg);
		}
		
		DirectedRegion[] regs= (DirectedRegion[]) gphase.tools.Arrays.toField(v);
		// no longer necessary, stop in VEGA st included, st not
		//regs[regs.length- 1].set3PrimeEdge(regs[regs.length- 1].get3PrimeEdge()+ 3);	// to include stop codon
		return regs;
	}
	
	public String getCDSSequenceNt() {
		DirectedRegion[] regs= getCDSRegions();	// not sorted
		java.util.Arrays.sort(regs, new DirectedRegion.DirectedPositionComparator());
		StringBuffer sb= new StringBuffer();
		for (int i = 0; i < regs.length; i++) 
			sb.append(Graph.readSequence(regs[i]));
		return sb.toString();
	}
	public byte getNmd() {
		return nmd;
	}
	public void setNmd(byte nmd) {
		this.nmd = nmd;
	}
	public Translation getPredORF() {
		return predORF;
	}
	public void setPredORF(Translation predORF) {
		this.predORF = predORF;
	}
	public String getHUGO() {
		return HUGO;
	}
	public void setHUGO(String hugo) {
		HUGO = hugo;
	}



}
