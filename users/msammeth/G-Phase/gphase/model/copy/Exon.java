 /*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.copy;

import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.ASMultiVariation.SpliceChainComparator;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import prefuse.util.UpdateListener;

/**
 * See <code>http://www.genomicglossaries.com/content/gene_def.asp</code>:
 * "What are the rules for deciding whether two partially overlapping mRNAs should be 
 * declared to be  alternative transcripts of the same gene or products of different genes? 
 * We have none."
 * <br><br>
 * Ensembl says in <code>http://whatislife.com/gene-definition.html</code>: 
 * A gene is a set of connected transcripts. A transcript is a set of exons via transcription 
 * followed (optionally) by pre-mRNA splicing. Two transcripts are connected if they share at 
 * least part of one exon in the genomic coordinates. At least one transcript must be expressed 
 * outside of the nucleus and one transcript must encode a protein (see footnotes).  
 * <br><br>
 * Hence, we (Sylvain, Tyler and me) define an exon to belong to strictly one gene.
 * <br><br>
 * 
 * 
 * @author micha
 */
public class Exon extends DirectedRegion {

	SpliceSite donor= null;
	SpliceSite acceptor= null;
	int frame= -1;
	int cds5Prime= 0;
	int cds3Prime= 0;
	EqualComparator defaultEqualComparator= new EqualComparator();
	
	public static class EqualComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			Exon ex1= (Exon) o1;
			Exon ex2= (Exon) o2;

			int start1= ex1.getStart();
			int end1= ex1.getEnd();
			String chr1= ex1.getChromosome();
			SpliceSite don1= ex1.getDonor();
			SpliceSite acc1= ex1.getAcceptor();
			int start2= ex2.getStart();
			int end2= ex2.getEnd();
			String chr2= ex2.getChromosome();
			SpliceSite don2= ex2.getDonor();
			SpliceSite acc2= ex1.getAcceptor();
			
			if (start1< start2)
				return -1;
			if (start2< start1)
				return 1;
			if (end1< end2)
				return -1;
			if (end2< end1)
				return 1;
			if (acc1.getType()< acc2.getType())
				return -1;
			if (acc2.getType()< acc1.getType())
				return 1;
			if (don1.getType()< don2.getType())
				return -1;
			if (don2.getType()< don1.getType())
				return 1;
			return 0;
		}
	}
	public static class IdentityComparator implements Comparator {

		/**
		 * checks chr, pos (start, end), types (don, acc)
		 */
		public int compare(Object o1, Object o2) {
		
				Exon ex1= (Exon) o1;
				Exon ex2= (Exon) o2;

				int start1= ex1.getStart();
				int end1= ex1.getEnd();
				String chr1= ex1.getChromosome();
				SpliceSite don1= ex1.getDonor();
				SpliceSite acc1= ex1.getAcceptor();
				int start2= ex2.getStart();
				int end2= ex2.getEnd();
				String chr2= ex2.getChromosome();
				SpliceSite don2= ex2.getDonor();
				SpliceSite acc2= ex1.getAcceptor();
				
				if (chr1.equals(chr2)&& start1== start2&& end1== end2&&
						don1== don2&& acc1== acc2)	// no object identity
					return 0;
				
					// non-overlapping, one before the other
				// cancelled, not working for neg. strand (clustering, sort array asc with start, end pos)
	//			if (end1< start2)
	//				return -1;		// one stops before the other
	//			if (end2< start1)
	//				return 1;
				
				if (!chr1.equals(chr2))
					return chr1.compareTo(chr2);
				
					// overlapping: none stops before the other
				if (start1< start2)
					return -1;
				if (start2< start1)
					return 1;
				
					// overlapping and same start position
				if (start1< end2)
					return -1;
				if (end2< start1)
					return 1;
				
				if (don1== null) {
					if (don2!= null)
						return -1; 	// abstract sites before real SSs
				} else {
					if (don2== null)
						return 1;
					else {
						if (don1.getPos()< don2.getPos())
							return -1;
						else if (don2.getPos()< don1.getPos())
							return 1;
					}
				}
				if (acc1== null) {
					if (acc2!= null)
						return 1; 	// abstract sites after real SSs
				} else {
					if (acc2== null)
						return -1;
					else {
						if (acc1.getPos()< acc2.getPos())
							return -1;
						else if (acc2.getPos()< acc1.getPos())
							return 1;
					}
				}
				
				//System.err.println("assertion in abstractregion.positioncomparator failed");
				return 0;	// identical positions --> never reached
					
				}
		
	}
	
	
	public static class PositionSSComparator extends AbstractRegion.PositionComparator {
	
		/**
		 * checks chr, pos (start, end), types (don, acc)
		 */
		public int compare(Object o1, Object o2) {
				
				int res= super.compare(o1, o2);
				if (res!= 0)
					return res;
			
				Exon ex1= (Exon) o1;
				Exon ex2= (Exon) o2;
	
				SpliceSite don1= ex1.getDonor();
				SpliceSite acc1= ex1.getAcceptor();
				SpliceSite don2= ex2.getDonor();
				SpliceSite acc2= ex1.getAcceptor();
				
				if (don1== null&& don2!= null)	// no object identity
					return -1;
				if (don1!= null&& don2== null)
					return 1;
				if (acc1== null&& acc2!= null)
					return -1;
				if (acc1!= null&& acc2== null)
					return 1;
				
				if (don1!= null&& don2!= null) {
					if (don1.getPos()< don2.getPos())
						return -1;
					if (don1.getPos()> don2.getPos())
						return 1;
				}
				if (acc1!= null&& acc2!= null) {
					if (acc1.getPos()< acc2.getPos())
						return -1;
					if (acc1.getPos()> acc2.getPos())
						return 1;
				}
				return 0;
					
			}
		
	}


	public boolean isCodingCompletely() {
		if (get5PrimeCDS()!= 0&& get3PrimeCDS()!= 0)
			return true;
		return false;
	}

	public boolean isCoding() {
		if (get5PrimeCDS()!= 0|| get3PrimeCDS()!= 0)
			return true;
		return false;
	}
	
	public boolean setPhase(int phase) {
		
		int f= 3- phase;
		if (frame>= 0&& frame!= f) {
			System.out.println("Exon "+this+" already used in another frame "+frame+"!");
			return false;
		}
		
		frame= f;
		return true;
	}
	
	public int get5PrimeCDS() {
		if (cds5Prime == 0) {
			for (int i = 0; i < transcripts.length; i++) {
				if (!transcripts[i].isCoding())
					continue;
				Translation trans= transcripts[i].getTranslations()[0];
				if (this.overlaps(trans)) {
					int p5= Math.max(this.get5PrimeEdge(), trans.get5PrimeEdge());
					if (cds5Prime== 0)
						cds5Prime= p5;
					else
						cds5Prime= Math.min(p5, cds5Prime);
				} 
			}
			
		}

		return cds5Prime;
	}
	
	public int get3PrimeCDS() {
		if (cds3Prime == 0) {
			for (int i = 0; i < transcripts.length; i++) {
				if (!transcripts[i].isCoding())
					continue;
				Translation trans= transcripts[i].getTranslations()[0];
				if (this.overlaps(trans)) {
					int p3= Math.min(this.get3PrimeEdge(), trans.get3PrimeEdge());
					if (cds3Prime== 0)
						cds3Prime= p3;
					else
						cds3Prime= Math.max(p3, cds5Prime);
				} 
			}
			
		}

		return cds3Prime;
	}
	
	public boolean isCoding5Prime() {
		int x= get5PrimeCDS();
		if (x!= 0&& get5PrimeCDS()== get5PrimeEdge())
			return true;
		return false;
	}
	
	/**
	 * 
	 * @return <true> if some transcript shows this exon as internal exon
	 */
	public boolean isInternal() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].isInternalExon(this))
				return true;
		}
		return false;
	}
	
	
	public boolean isCodingSomewhere5Prime() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.get5PrimeEdge()<= get5PrimeEdge()&& tln.get3PrimeEdge()>= get5PrimeEdge())
				return true;
		}
		return false;	 
	}
	public boolean isUpstreamRegion(DirectedRegion reg) {
		if (reg.get3PrimeEdge()+ 1== get5PrimeEdge())
			return true;
		return false;
	}
	
	public boolean isDownstreamRegion(DirectedRegion reg) {
		if (get3PrimeEdge()+1== reg.get5PrimeEdge())
			return true;
		return false;
	}

	public boolean overlapsCDS() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.overlaps(this))
				return true;
		}
		return false;	 
	}
	
	public boolean isCodingSomewhere3Prime() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (!getTranscripts()[i].isCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.get5PrimeEdge()<= get3PrimeEdge()&& tln.get3PrimeEdge()>= get3PrimeEdge())
				return true;
		}
		return false;	 
	}
	
	
	public boolean isCoding3Prime() {
		int x= get3PrimeCDS();
		if (x!= 0&& get3PrimeCDS()== get3PrimeEdge())
			return true;
		return false;
	}
	
	public boolean setFrame(int newFrame) {
		
		if (frame!= 0&& frame!= newFrame) {
			System.out.println("Exon "+this+" already used in another frame "+frame+"!");
			return false;
		}
				
		frame= newFrame;
		return true;
	}
	
	public boolean removeTranscript(Transcript trans) {
		
		if (getGene().isConstruct()) {
			Object o= transHash.remove(trans);
			if (o== null)
				return false;
			return true;
		} else {
			Transcript[] newTranscripts= new Transcript[transcripts.length- 1];
			int pos= 0;
			boolean flag= false;
			for (int i = 0; i < transcripts.length; i++) 
				if (transcripts[i]!= trans)
					newTranscripts[pos++]= transcripts[i];
				else
					flag= true;
			if (flag)
				transcripts= newTranscripts;
			return flag;
		}
	}
	
	public boolean addTranscript(Transcript trans) {
		
		if (trans.getGene().isConstruct()) {
			if (transHash== null)
				transHash= new HashMap<Transcript, Transcript>(1,1f);
			if (transHash.get(trans)!= null)
				return false;
			transHash.put(trans, trans);
			return true;
		} else {
			if (transcripts== null) {
				transcripts= new Transcript[] {trans};
				return true;
			}
			
				// search transcript
			for (int i = 0; i < transcripts.length; i++) 
				if (transcripts[i].getStableID().equalsIgnoreCase(trans.getStableID()))
					return false;
			
				// add transcript
			Transcript[] nTranscripts= new Transcript[transcripts.length+ 1];
			for (int i= 0; i < transcripts.length; i++) 
				nTranscripts[i]= transcripts[i];
			nTranscripts[nTranscripts.length- 1]= trans;
			transcripts= nTranscripts;
			return true;
		}
	}


	/**
	 * @param b
	 */
	public boolean checkStrand(boolean b) {
		
		return (b== getGene().isForward());
	}
	public boolean checkStrand(String newStrand) {
		
		String nStrand= newStrand.trim();
		if (nStrand.equals("1"))	// || nStrand.equals("forward")
			return checkStrand(true);
		else if (nStrand.equals("-1"))	// || nStrand.equals("reverse")
			return checkStrand(false);
		
		return false; // error			
	}
	
	Transcript[] transcripts= null;
	HashMap<Transcript, Transcript> transHash= null;
	
	String exonID= null;	
	protected Exon() {
		// for subclasses
	}
	public Exon(Transcript newTranscript, String stableExonID, int start, int end) {

		this.strand= newTranscript.getStrand();
		this.chromosome= newTranscript.getChromosome();
		setStart(start);
		setEnd(end);
		
			// decompose ID
		//this.exonID= stableExonID;
		setID("exon");
		if (stableExonID!= null)
			exonID= stableExonID;
		
		// check for duplicate
		addTranscript(newTranscript);
		newTranscript.addExon(this);
					
		for (int i = 0; i < getTranscripts().length; i++) 
			getTranscripts()[i].updateBoundaries(this);
	}	
	/**
	 * @return
	 */
	public String getExonID() {
		return exonID;
	}

	/**
	 * @return
	 */
	public Gene getGene() {
		if (transcripts== null)
			return transHash.keySet().iterator().next().getGene();
		return transcripts[0].getGene();
	}
	
	public Transcript[] getTranscripts() {
		if (transcripts== null) {
			transcripts= new Transcript[transHash.size()];
			Object[] o= transHash.keySet().toArray();
			for (int i = 0; i < transcripts.length; i++) {
				transcripts[i]= (Transcript) o[i];
			}
			if (getGene().isConstruct()) {
				Transcript[] t= transcripts;
				transcripts= null;
				return t;
			} else
				transHash= null;
		}
		return transcripts;
	}
	

	public int getStrand() {
		if (strand!= 0)
			return strand;
		if (getTranscripts()!= null&& transcripts[0]!= null) {
			if (transcripts[0].getStrand()!= 0)
				return transcripts[0].getStrand();
			if (transcripts[0].getGene()!= null)
				return transcripts[0].getGene().getStrand();
		}
		return 0;
	}
	
	/**
	 * @param transcripts
	 */
	public void setTranscripts(Transcript[] transcripts) {
		this.transcripts = transcripts;
	}
	
	public String toString() {
		// return getExonID();
		String res= "";
		if (getAcceptor()!= null)
			res+= getAcceptor();
		else
			res+=0;
		res+="==";
		if (getDonor()!= null)
			res+= getDonor();
		else
			res+=0;
		
		return res;
	}

	public String toPosString() {
		return getGene().getChromosome()+" "+getStart()+" "+getEnd();
	}

	public String getChromosome() {
		if (getGene()== null) {
			if (chromosome!= null)
				return chromosome.toUpperCase();
			return chromosome;
		}
		return getGene().getChromosome();
	}
	
	public Species getSpecies() {
		if (species!= null) {
			return species;
		} else if (transcripts!= null&& transcripts[0].getSpecies()!= null) {
			return transcripts[0].getSpecies();
		} else if (getGene()!= null&& getGene().getSpecies()!= null)
			return getGene().getSpecies();
		
		return null;
	}

	static final long serialVersionUID = 8914674126313232057L;
	public SpliceSite getAcceptor() {
		return acceptor;
	}
	public void setAcceptor(SpliceSite acceptor) {
		acceptor.addExon(this);
		this.acceptor = acceptor;
	}
	public SpliceSite getDonor() {
		//donor.addExon(this);
		return donor;
	}
	public void setDonor(SpliceSite donor) {
		donor.addExon(this);
		this.donor = donor;
	}
	
	public boolean replaceSite(SpliceSite oldSite, SpliceSite newSite) {
		if (oldSite== this.acceptor) {
			this.acceptor= newSite;
			return true;
		}
		if (oldSite== this.donor) {
			this.donor= newSite;
			return true;
		}
		
		return false;
	}
	
	public int get3PrimeFrame() {
		int fr= getFrame();
		int len= get3PrimeCDS()- get5PrimeCDS()+ 1;
		if (fr>= len- 1)	 
			fr-= (len- 1);
		else 
			fr= (fr+ ((len- 1)% 3))% 3;		// simpler?, inverse op of %??
		
		return fr;		
	}
	
	
	public int getFrame() {
		if (frame== -1) {
			
			for (int i = 0; i < transcripts.length; i++) {
				if ((!transcripts[i].isCoding())|| (!transcripts[i].getTranslations()[0].contains(get5PrimeCDS())))
					continue;
				int pos= transcripts[i].getTranslations()[0].getTranslatedPosition(get5PrimeCDS());
				int f= pos% 3;
				if (frame== -1)
					frame= f;
				else
					assert(frame== f);
			}
			
		}
		return frame;
	}

	public int getEndCDS() {
		if (isForward())
			return get3PrimeCDS();
		return get5PrimeCDS();
	}

	public void setEndCDS(int endCDS) {
		if (getGene().getStrand()< 0)
			this.cds5Prime= -Math.abs(endCDS);
		else
			this.cds3Prime = endCDS;
	}

	public int getStartCDS() {
		if (isForward())
			return get5PrimeCDS();	// here not getter method, avoid init for buildup
		return get3PrimeCDS();
	}

	public void setStartCDS(int startCDS) {
		if (getGene().getStrand()< 0)
			cds3Prime= -Math.abs(startCDS);
		else
			cds5Prime = startCDS;
	}
	
	public void extendStartCDS(int nuStartCDS) {
		if ((this.getStartCDS()== 0)|| (Math.abs(nuStartCDS)< Math.abs(this.getStartCDS())))
			setStartCDS(nuStartCDS);
	}

	public void extendEndCDS(int nuEndCDS) {
		if ((this.getEndCDS()== 0)|| (Math.abs(nuEndCDS)> Math.abs(this.getEndCDS()))) 
			setEndCDS(nuEndCDS);
	}

	public HashMap<Transcript, Transcript> getTransHash() {
		return transHash; 
	}

	public EqualComparator getDefaultEqualComparator() {
		return defaultEqualComparator;
	}

}
