/*
 * Created on Feb 22, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.io.gtf.GTFObject;
import gphase.tools.Arrays;

import java.io.Serializable;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import org.freehep.graphicsio.swf.SWFAction.RemoveSprite;

import com.sun.net.ssl.SSLContext;

/**
 * 
 * 
 * @author msammeth
 */
public class SpliceSite extends AbstractSite {

	public static final String ATTRIBUTE_ID_SS_TYPE= "ss_type";
	
	public final static byte TYPE_NOT_INITED= -1;
	public final static byte TYPE_TSS= 0;
	public final static byte TYPE_ACCEPTOR= 1;
	public final static byte TYPE_DONOR= 2;
	public final static byte TYPE_TES= 3;
	
	
	public final static byte NOTYPE_SS= 0;
	public final static byte CONSTITUTIVE= 1;
	public final static byte ALTERNATIVE= 2;
	
	public final static int DEFAULT_DONOR_5FLANK= 2;
	public final static int DEFAULT_DONOR_3FLANK= 6;
	public final static int DEFAULT_ACCEPTOR_5FLANK= 15;
	public final static int DEFAULT_ACCEPTOR_3FLANK= 2;
	
	static final long serialVersionUID = 422169949942461293L;
	static final int SPLICE_SITE_FLANK_SIZE= 20; // Alternative splicing in the human, mouse and rat genomes is associated with an increased frequency of exon creation and/or loss
												 // Barmak Modrek & Christopher J Lee
	public static final int DELTA_RANGE= 20;
	
	public static final int NOBORU_DON5_EXTENT= 20;
	public static final int NOBORU_DON3_EXTENT= 18;
	public static final int NOBORU_ACC5_EXTENT= 21;
	public static final int NOBORU_ACC3_EXTENT= 17;

	public static final int SHAPIRO_DON5_FIRST= -3;
	public static final int SHAPIRO_DON3_LAST= +6;
	public static final int SHAPIRO_ACC5_FIRST= -13;
	public static final int SHAPIRO_ACC3_LAST= +1;

	public static final String CANONICAL_DONOR = "GT";

	public static final String CANONICAL_ACCEPTOR = "AG";


	static PositionTypeComparator defaultPositionTypeComparator= null;

	public static char CHAR_ACCEPTOR = '-';

	public static char CHAR_DONOR = '^';

	public static char CHAR_TES = ';';

	public static char CHAR_TSS = '*';
	public static char CHAR_UNDEFINED = '?';
	
	public static PositionTypeComparator getDefaultPositionTypeComparator() {
		if (defaultPositionTypeComparator == null) {
			defaultPositionTypeComparator = new PositionTypeComparator();			
		}

		return defaultPositionTypeComparator;
	}
	
	
	float scoreU12= Float.NaN;
	
	HashMap hitparade= null;
	HashMap homologs= null;	// maps genes to exon homologs
	
	Exon[] exons= null;
	byte type= TYPE_NOT_INITED;
	byte modality= NOTYPE_SS;
	boolean constitutive= true;
	String dinucleotide= null, stringRep= null;

	
	ASVariation[] asVars = null;
	
	public static class PositionComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			SpliceSite ss1= (SpliceSite) arg0;
			SpliceSite ss2= (SpliceSite) arg1;
			if (ss1.getPos()< ss2.getPos())
				return -1;
			if (ss2.getPos()< ss1.getPos())
				return 1;
			return 0;
		}
	}
	
	public static class PositionTypeComparator extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}

		public int compare_old(Object arg0, Object arg1) {
				// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}
	}

	public static class PositionEqualSSTypeComparator extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			// positions equal: TSS before donor, acc, TES
			SpliceSite s0= (SpliceSite) arg0;
			SpliceSite s1= (SpliceSite) arg1;
			if (s0.isTSS()&& !s1.isTSS())
				return -1;
			if ((!s0.isTSS())&& s1.isTSS())
				return 1;
			if (s0.isDonor()&& !s1.isDonor())
				return -1;
			if ((!s0.isDonor())&& s1.isDonor())
				return 1;
			if (s0.isAcceptor()&& !s1.isAcceptor())
				return -1;
			if ((!s0.isAcceptor())&& s1.isAcceptor())
				return 1;
			if (s0.isTES()&& !s1.isTES())
				return -1;
			if ((!s0.isTES())&& s1.isTES())
				return 1;
			
			return 0; 	// everything equal
		}
	
		public int compare_old(Object arg0, Object arg1) {
				// positions equal
			if (arg0 instanceof SpliceSite) {
				if (arg1 instanceof SpliceSite) {
					SpliceSite s0= (SpliceSite) arg0;
					SpliceSite s1= (SpliceSite) arg1;
					if (s0.isDonor()&& !s1.isDonor())
						return 1;	// for 1-lentgth exons, sort acceptor before donor (same pos)
					if (s1.isDonor()&& !s0.isDonor())
						return -1;
					return 0;	// same type and position					
				} else 
					return 1;	// arg1 is no SS					
			} else {
				if (arg1 instanceof SpliceSite) 
					return -1;	// arg0 is no SS
				else {	// 2 no SS
					boolean tss0= ((AbstractSite) arg0).isTSS();
					boolean tss1= ((AbstractSite) arg1).isTSS();
					if (tss0== tss1)
						return 0;
					if (tss0)
						return -1;
					return 1;
				}
			}
		}
	}

	public static boolean isU2ss(String id) {
		if (id.contains("U2"))
			return true;
		return false;
	}
	
	public static int[] getPositions(SpliceSite[] ss) {
		if (ss== null)
			return null;
		int[] pos= new int[ss.length];
		for (int i = 0; i < pos.length; i++) 
			pos[i]= ss[i].getPos();
		return pos;
	}
	
	public SpliceSite(int pos, byte newType) {

		super(pos);
		setType(newType);		
	}
	
	public SpliceSite(int pos, byte newType, Exon exon) {
		this(pos,newType);
		addExon(exon);
	}
	
	public void removeExon(Exon ex) {
		if (ex== null)
			return;
		
			// remove exon
		Exon[] newExons= new Exon[exons.length- 1];
		int p= 0;
		for (int i = 0; i < exons.length; i++) {
			if (exons[i]!= ex)
				newExons[p++]= exons[i];
		}
		exons= newExons;
	}
	
	public boolean addExon(Exon newExon) {
		addTranscripts(newExon.getTranscripts());
		if (exons== null) {
			exons= new Exon[] {newExon};
			return true;
		}
		
		Comparator compi= new AbstractRegion.PositionComparator();
		int i;
		for (i = 0; i < exons.length; i++) 
			//if (compi.compare(newExon, exons[i])== 0)
			if (newExon== exons[i])	// have to hold redundantly ALL exons, for SpliceSite.remove(Exon)
				break;
		if (i== exons.length) {
			exons= (Exon[]) Arrays.add(exons, newExon);
			return true;
		}
		return false;
	}
	
	/**
	 * checks for identity, not for homology
	 */
	public boolean equals(Object obj) {
		
		SpliceSite otherSS= (SpliceSite) obj;
		if (getPos()== otherSS.getPos()&& getType()== otherSS.getType())
			return true;
		return false;
		
//		if (!super.equals(anotherSS))
//			return false;
//		
//		if (!(anotherSS instanceof SpliceSite))
//			return false;
//		
//		SpliceSite s2= (SpliceSite) anotherSS;
//		if (s2.isDonor()!= isDonor())
//			return false;
////		SpliceSite aSS= (SpliceSite) anotherSS;
////		if (gene!= aSS.getGene()|| pos!= aSS.getPos())
////			return false;
//		return true;
	}
	
//	---------------------------------------------- 
//	hashCode 
//	public int hashCode() 
//	Returns a hash code value for the object. This method is supported for the 
//	benefit of hashtables such as those provided by java.util.Hashtable. 
//	The general contract of hashCode is: 
//	* Whenever it is invoked on the same object more than once during an 
//	execution of a Java application, the hashCode method must consistently 
//	return the same integer, provided no information used in equals comparisons 
//	on the object is modified. This integer need not remain consistent from one 
//	execution of an application to another execution of the same application. 
//	* If two objects are equal according to the equals(Object) method, then 
//	calling the hashCode method on each of the two objects must produce the same 
//	integer result. 
//	* It is not required that if two objects are unequal according to the 
//	equals(java.lang.Object) method, then calling the hashCode method on each of 
//	the two objects must produce distinct integer results. However, the 
//	programmer should be aware that producing distinct integer results for 
//	unequal objects may improve the performance of hashtables. 
//
//
//	As much as is reasonably practical, the hashCode method defined by class 
//	Object does return distinct integers for distinct objects. (This is 
//	typically implemented by converting the internal address of the object into 
//	an integer, but this implementation technique is not required by the JavaTM 
//	programming language.) 
//	----------------------------------------------
	public int hashCode() {
		// not required, read text
		// bullshit, gotta change otherwise ss overwrite as in hash
		// and the other way
		return getPos(); //super.hashCode();
		
	}
	
	public boolean isDonor() {
		return (type== TYPE_DONOR);
	}

	public boolean isTSS() {
		return (type== TYPE_TSS);
	}
	
	public boolean isTES() {
		return (type== TYPE_TES);
	}
	
	public boolean isSpliceSite() {
		return (isDonor()|| isAcceptor());
	}

	public boolean isAcceptor() {
		return (type== TYPE_ACCEPTOR);
	}
	
	public boolean isCanonical() {
		if ((!this.isDonor())&& (!this.isAcceptor()))
			System.err.println("Not canonical() called for non-splicesite");
		
		String seq= getDinucleotide();		
		if ((isDonor()&& seq.equalsIgnoreCase(CANONICAL_DONOR))||
				(isAcceptor()&& seq.equalsIgnoreCase(CANONICAL_ACCEPTOR)))
			return true;
		return false;				
	}
	
	public String getDinucleotide() {
		if (dinucleotide == null) {
			dinucleotide = Graph.readSequence(this, 0, 0);
		}

		return dinucleotide;
	}
	
	public void setType(byte newType) {
		this.type = newType;
		if (isDonor())
			setID(GTFObject.FEATURE_ID_DONOR);
		else if (isAcceptor())
			setID(GTFObject.FEATURE_ID_ACCEPTOR);
	}
	public Exon[] getExons() {
		return exons;
	}
	public void setExons(Exon[] exons) {
		this.exons = exons;
	}
	public String toString() {
		if (stringRep == null) {	// || stringRep.charAt(stringRep.length()- 1)== '?'
			if (isDonor())
				stringRep= getPos()+ "^";
			else if (isAcceptor())
				stringRep= getPos()+ "-";
			else if (isTSS())
				stringRep= getPos()+ "*";
			else if (isTES())
				stringRep= getPos()+ ";";
			else 
				stringRep= getPos()+ "?";
		}

		return stringRep;
	}
	
	public String toOutString() {
		String result= getGene().getChromosome()+ " "+ Math.abs(getPos())+ " "+ getGene().getStrand();
		return result;
	}

	public boolean isCDSRealTranscript() {
		
				// net gut, twilight-frontier events !
	//		for (int i = 0; i < asVars.length; i++) 
	//			if (asVars[i].isProteinCoding())
	//				return true;
			
			int ctr= 0;
			for (int i = 0; i < transcripts.length; i++) {
				if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
					continue;
				if (transcripts[i].getTranslations()[0].contains(this))
					ctr++;	// at least in the CDS of one transcript
				else 
					return false;
			}
			return (ctr> 0);	// if all are non-annotated, ok, its never in CDS
		}

	public boolean isCDSMaxTranscript() {
		
		// net gut, twilight-frontier events !
//		for (int i = 0; i < asVars.length; i++) 
//			if (asVars[i].isProteinCoding())
//				return true;
	
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
				continue;
			if (transcripts[i].getTranslations()[0].contains(getPos()))
				return true;
		}
		return false;	// if all are non-annotated, ok, its never in CDS
	}
	
	public boolean is5UTRRealTranscript() {
		int ctr= 0;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1)
				continue;
			if (this.getPos()>= transcripts[i].getTranslations()[0].get5PrimeEdge())
				return false;	// at least outside of 5UTR
			++ctr;
		}
		return (ctr> 0);	// at least one transcript w annotated CDS found
	}
	
	public String getGenicLocation() {
		Transcript[] t= getTranscripts();
		boolean cds= false, utr5= false, utr3= false;
		int pos= getRealSSpos();
		for (int i = 0; i < t.length; i++) {
			if (!t[i].isCoding())
				continue;
			if (t[i].getTranslations()[0].contains(pos))
				cds= true;
			else if (pos < t[i].getTranslations()[0].get5PrimeEdge())
				utr5= true;
			else if (pos > t[i].getTranslations()[0].get3PrimeEdge())
				utr3= true;
		}

		if (cds)
			return "CDS";
		if (utr5&& !utr3)
			return "5UTR";
		if (utr3&& !utr5)
			return "3UTR";
		return "UNDEFINED";
	}
	
	public boolean is5UTRMaxTranscript() {
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1
					|| transcripts[i].getTranslations()[0].isOpenEnded5())
				continue;
			if (this.getPos()< transcripts[i].getTranslations()[0].get5PrimeEdge())
				return true;	// at least outside of 5UTR
		}
		return false;
	}

	public boolean is3UTRMaxTranscript() {
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranslations()== null|| transcripts[i].getTranslations().length< 1||
					transcripts[i].getTranslations()[0].isOpenEnded3())
				continue;
			if (this.getPos()> transcripts[i].getTranslations()[0].get3PrimeEdge())
				return true;	// at least outside of 5UTR
		}
		return false;
	}


	/**
	 * 
	 * @param constitutive
	 * @deprecated as to link SSs with ASVars
	 */
	public void setConstitutive(boolean constitutive) {
		this.constitutive = constitutive;
	}
	
	public boolean addHit(Gene g, PWHit h) {
			
			Vector v= null;
			if (hitparade== null)	// create new
				hitparade= new HashMap();
			else 
				v= (Vector) hitparade.get(g);	// lookup existing
			
			if (v== null)
				v= new Vector();
			else 
				for (int i = 0; i < v.size(); i++) {	// check for already added hit
					PWHit hit= (PWHit) v.elementAt(i);
					if (hit.getOtherObject(this)== h.getOtherObject(this))
						return false;
				}
			
				// keep sorted with ascending costs
	//		int i= 0;
	//		for (i = 0; i < v.size(); i++) 
	//			if (((PWHit) v.elementAt(i)).getCost()> h.getCost())
	//				break;
	//		v.insertElementAt(h, i);	// add
			v.add(h);
			hitparade.put(g, v);
			return true;
		}

	/**
	 * returns hits to exons of all homolog genes
	 * @return
	 */
	public PWHit[] getHits() {
		
			// count values
		Object[] oa= hitparade.values().toArray();
		int size= 0;
		for (int i = 0; i < oa.length; i++) 
			size+= ((Vector) oa[i]).size();
		
		PWHit[] result= new PWHit[size];
		int x= 0;
		for (int i = 0; i < oa.length; i++) {
			Vector v= (Vector) oa[i];
			for (int j = 0; j < v.size(); j++) 
				result[x++]= (PWHit) v.elementAt(j);
		}
		
		return result;
	}
	
	/**
	 * from (-x to + y) number of exonic/intronic positions to
	 * number of positions before and after ss
	 * before and after ss
	 * @param beforeSS
	 * @param afterSS
	 * @return
	 */
	public int[] convertToExonIntronPositions(int left, int right) {
		int[] result= new int[2];
		if (isDonor()) {
			result[0]= Math.abs(left);
			result[1]= right- 2;
		} else {
			result[0]= Math.abs(left)- 2;
			result[1]= right;
		}
		return result;
	}
	
	public DirectedRegion getShapiroRegion() {
		
		if (isDonor()) {
			int[] extent= convertToExonIntronPositions(SHAPIRO_DON5_FIRST, SHAPIRO_DON3_LAST);
			return getRegion(extent[0], extent[1]);
		} else {
			int[] extent= convertToExonIntronPositions(SHAPIRO_ACC5_FIRST, SHAPIRO_ACC3_LAST);
			return getRegion(extent[0], extent[1]);
		}
			
	}
	
	public DirectedRegion getRegion(int downstream, int upstream) {
		int start= getPos();
		int end= getPos();
		
		if (isDonor()) {
			start-= downstream- 1;	// -2
			end+= upstream+ 2;	// +4
		} else {
			start-= downstream+ 2;	// +1
			end+= upstream- 1;	// 
		}
			
		DirectedRegion reg= new DirectedRegion(start, end, getTranscripts()[0].getStrand());
		reg.setChromosome(getTranscripts()[0].getChromosome());
		reg.setSpecies(getTranscripts()[0].getSpecies());
		return reg;
	}

	public DirectedRegion getNoboruRegion() {
		int start= getPos();
		int end= getPos();
		if (isDonor()) {
			start-= NOBORU_DON5_EXTENT- 1;	// -2
			end+= NOBORU_DON3_EXTENT+ 2;	// +4
		} else {
			start-= NOBORU_ACC5_EXTENT+ 2;	// +1
			end+= NOBORU_ACC3_EXTENT- 1;	// 
		}
			
		DirectedRegion reg= new DirectedRegion(start, end, getTranscripts()[0].getStrand());
		reg.setChromosome(getTranscripts()[0].getChromosome());
		reg.setSpecies(getTranscripts()[0].getSpecies());
		return reg;
	}

	public PWHit[] getHits(Gene g) {
		
		if (hitparade== null|| hitparade.get(g)== null)
			return null;
		
		Vector v= (Vector) hitparade.get(g);
		PWHit[] result= new PWHit[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (PWHit) v.elementAt(i);
		
		return result;
	}

	// assuming just one homolog
	public boolean addHomolog(Gene g, SpliceSite e) {
	
		if (homologs== null)	// create new
			homologs= new HashMap();
	
		if (homologs.get(g)!= null) {
			if (homologs.get(g)!= e)
				System.err.println("Multiple exon homology: "+ e+ ", "+homologs.get(g));	
			return false;
		}
	
		homologs.put(g, e);
		return true;
	}

	public float getScoreU12() {
		return scoreU12;
	}

	public void setScoreU12(float scoreU12) {
		this.scoreU12 = scoreU12;
	}

	public byte getType() {
		return type;
	}

	public char getSiteSymbol() {
		if (isDonor())
			return CHAR_DONOR;
		if (isAcceptor())
			return CHAR_ACCEPTOR;
		if (isTSS())
			return CHAR_TSS;
		if (isTES())
			return CHAR_TES;
		
		return CHAR_UNDEFINED;
	}
	
	/**
	 * +1 for donors, -1 for acceptors
	 * @return
	 */
	public int getRealSSpos() {
		if (isDonor())
			return getPos()+ 1;
		if (isAcceptor())
			return getPos()- 1;
		return getPos();
	}

	public void addASVar(ASVariation newASVar) {
		if (asVars== null) 
			asVars= new ASVariation[] {newASVar};
		else 
			asVars= (ASVariation[]) Arrays.add(asVars, newASVar);
	}

	public boolean isAlternative() {
		//return (!isConstitutive());
		return (modality== ALTERNATIVE);
	}

	public boolean isConstitutive() {
		//return constitutive;
		//return (asVars== null|| asVars.length< 1);
		return (modality== CONSTITUTIVE);
	}

	public ASVariation[] getAsVars() {
		return asVars;
	}

	public void removeASVar(ASVariation var) {
		asVars= (ASVariation[]) Arrays.remove(asVars, var);
	}

	public byte getModality() {
		return modality;
	}

	public void setModality(byte modality) {
		this.modality = modality;
	}
}
