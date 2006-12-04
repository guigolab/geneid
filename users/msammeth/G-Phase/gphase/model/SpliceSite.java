/*
 * Created on Feb 22, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.tools.Arrays;

import java.io.Serializable;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import org.freehep.graphicsio.swf.SWFAction.RemoveSprite;

/**
 * 
 * 
 * @author msammeth
 */
public class SpliceSite extends AbstractSite {

	public final static int NOTYPE_SS= 0;
	public final static int CONSTITUTIVE_SS= 1;
	public final static int ALTERNATE_SS= 2;
	public final static int DEFAULT_DONOR_5FLANK= 2;
	public final static int DEFAULT_DONOR_3FLANK= 6;
	public final static int DEFAULT_ACCEPTOR_5FLANK= 15;
	public final static int DEFAULT_ACCEPTOR_3FLANK= 2;
	
	static final long serialVersionUID = 422169949942461293L;
	static final int SPLICE_SITE_FLANK_SIZE= 20; // Alternative splicing in the human, mouse and rat genomes is associated with an increased frequency of exon creation and/or loss
												 // Barmak Modrek & Christopher J Lee
	public static final int DELTA_RANGE= 20;
	
	
	HashMap hitparade= null;
	HashMap homologs= null;	// maps genes to exon homologs
	
	Exon[] exons= null;
	boolean donor= false;
	boolean constitutive= true;
	ASVariation[] asVars= null;
	public static class PositionComparator extends AbstractSite.PositionComparator {
		public int compare(Object arg0, Object arg1) {
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			SpliceSite s0= (SpliceSite) arg0;
			SpliceSite s1= (SpliceSite) arg1;
			if (s0.isDonor()!= s1.isDonor())
				return -1;
			return 0;
		}
	}
	
	public SpliceSite(Gene gene, int pos, boolean donor) {

		super(pos);
		this.donor= donor;
		this.gene= gene;
	}
	
	public SpliceSite(Gene gene, int pos, boolean donor, Exon exon) {
		this(gene,pos,donor);
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
	
	public void addASVar(ASVariation newASVar) {
		if (asVars== null) 
			asVars= new ASVariation[] {newASVar};
		
		asVars= (ASVariation[]) Arrays.add(asVars, newASVar);
	}
	
	/**
	 * checks for identity, not for homology
	 */
	public boolean equals(Object anotherSS) {
		
		if (!super.equals(anotherSS))
			return false;
		
		if (!(anotherSS instanceof SpliceSite))
			return false;
		
		SpliceSite s2= (SpliceSite) anotherSS;
		if (s2.isDonor()!= isDonor())
			return false;
//		SpliceSite aSS= (SpliceSite) anotherSS;
//		if (gene!= aSS.getGene()|| pos!= aSS.getPos())
//			return false;
		return true;
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
		return super.hashCode();
		
	}
	
	public boolean isDonor() {
		return donor;
	}
	public boolean isAcceptor() {
		return !donor;
	}
	public void setDonor(boolean donor) {
		this.donor = donor;
	}
	public Exon[] getExons() {
		return exons;
	}
	public void setExons(Exon[] exons) {
		this.exons = exons;
	}
	public String toString() {
		String result= "";
		if (!isDonor())
			result+= ">";				
		result+= getPos();
		if (isDonor())
			result+= ">";
		return result;
	}
	
	public String toOutString() {
		String result= getGene().getChromosome()+ " "+ Math.abs(getPos())+ " "+ getGene().getStrand();
		return result;
	}

	public boolean isConstitutive() {
		//return constitutive;
		return (asVars== null|| asVars.length< 1);
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
	
	public void removeASVar(ASVariation var) {
		asVars= (ASVariation[]) Arrays.remove(asVars, var);
	}

	public ASVariation[] getAsVars() {
		return asVars;
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
}
