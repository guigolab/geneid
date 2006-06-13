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

/**
 * 
 * 
 * @author msammeth
 */
public class SpliceSite extends AbstractSite {

	public final static int NOTYPE_SS= 0;
	public final static int CONSTITUTIVE_SS= 1;
	public final static int ALTERNATE_SS= 2;
	static final long serialVersionUID = 422169949942461293L;
	static final int SPLICE_SITE_FLANK_SIZE= 20; // Alternative splicing in the human, mouse and rat genomes is associated with an increased frequency of exon creation and/or loss
												 // Barmak Modrek & Christopher J Lee

	
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

		this.pos= pos;
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
		Exon[] newExons= new Exon[exons.length- 1];
		int p= 0;
		for (int i = 0; i < exons.length; i++) {
			if (exons[i]!= ex)
				newExons[p++]= exons[i];
		}
		exons= newExons;
	}
	
	public boolean addExon(Exon newExon) {
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
		
		Arrays.add(asVars, newASVar);
	}
	
	/**
	 * checks for identity, not for homology
	 */
	public boolean equals(Object anotherSS) {
		
		if (!(anotherSS instanceof SpliceSite))
			return false;
		
		SpliceSite aSS= (SpliceSite) anotherSS;
		if (gene!= aSS.getGene()|| pos!= aSS.getPos())
			return false;
		return true;
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
		return constitutive;
	}

	public void setConstitutive(boolean constitutive) {
		this.constitutive = constitutive;
	}

	public ASVariation[] getAsVars() {
		return asVars;
	}
}
