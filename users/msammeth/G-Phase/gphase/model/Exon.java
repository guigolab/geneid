 /*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

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

	transient HashMap hitparade= null;	// do not serialize
	HashMap homologs= null;	// maps genes to exon homologs
							// also captures the alternatively spliced exon
							// in same gene for: a3, a5, me
	
	SuperExon superExon= null;
	SpliceSite donor= null;
	SpliceSite acceptor= null;
	int frame= -1;
	int cds5Prime= 0;
	int cds3Prime= 0;
	
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
	
	public boolean isConstitutive(int codingCode) {
		getGene().getASVariations(codingCode);
		if (getDonor()!= null&& (!getDonor().isConstitutive()))			
			return false;
		if (getAcceptor()!= null&& (!getAcceptor().isConstitutive()))			
			return false;
		return true;
		
//		for (int i = 0; i < getTranscripts().length; i++) {
//			Exon[] ex= getTranscripts()[i].getExons();
//			int j;
//			for (j = 0; j < ex.length; j++) {
//				if (this.equals(ex[j]))  	// o ident?!
//					break;
//			}
//			if (j== ex.length)
//				return false;
//		}
//		return true;	 
	}
	
	public boolean hasVariation(String varCode, int codingCode) {
		if (getDonor()== null|| getAcceptor()== null)
			return false;
		ASVariation[] vars= getGene().getASVariations(codingCode);
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			if (!vars[i].toString().equals(varCode))
				continue;
			Comparator compi= new SpliceSite.PositionComparator();
			SpliceSite[] sc1= vars[i].getSpliceChain1();
			if (sc1!= null) {
				int p1= Arrays.binarySearch(sc1, getAcceptor(), compi);
				int p2= Arrays.binarySearch(sc1, getDonor(), compi);
				if (p1>= 0&& p2>= 0&& p2==p1+1)
					return true;
			}
			SpliceSite[] sc2= vars[i].getSpliceChain2();
			if (sc2!= null) {
				int p1= Arrays.binarySearch(sc2, getAcceptor(), compi);
				int p2= Arrays.binarySearch(sc2, getDonor(), compi);
				if (p1>= 0&& p2>= 0&& p2==p1+1)
					return true;
			}
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
		
		if (superExon!= null)
			return superExon.getHits(g);
		
		if (hitparade== null|| hitparade.get(g)== null)
			return null;
		
		Vector v= (Vector) hitparade.get(g);
		PWHit[] result= new PWHit[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (PWHit) v.elementAt(i);
		
		return result;
	}
	
	public boolean removeHit(Gene g, PWHit h) {
	
		Vector v;
		if (hitparade== null|| (v= (Vector) hitparade.get(g))== null)
			return false;

		return v.remove(h);
	}
	
	public boolean removeTranscript(Transcript trans) {
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
	
	// assuming just one homolog
	public boolean addHomolog(Gene g, Exon e) {

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
	
	// assuming just one homolog
	public Exon getHomolog(Gene g) {
		
		if (homologs== null)
			return null;
		return (Exon) homologs.get(g);
	}
	
	public boolean addHit(Gene g, PWHit h) {
		
		if (superExon!= null)
			superExon.addHit(g, h);	// delegate to super-exon if there is some
		
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
	
	public PWHit[] getBestHits(Gene g) {
		
		if (hitparade== null|| hitparade.get(g)== null)
			return null;
		
		Vector v= (Vector) hitparade.get(g);
		Vector bestHits= new Vector();
		int bestScore= (-1);
		for (Iterator iter = v.iterator(); iter.hasNext();) {
			PWHit tmpHit = (PWHit) iter.next();
			if (tmpHit.getScore()== bestScore) 			// more hits
				bestHits.add(tmpHit);
			else if (tmpHit.getScore()> bestScore) {	// new best score
				bestHits.removeAllElements();
				bestHits.add(tmpHit);
				bestScore= tmpHit.getScore();
			}
		}
		
		PWHit[] result= new PWHit[bestHits.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (PWHit) bestHits.elementAt(i);
		return result;
	}
	
	public boolean addTranscript(Transcript trans) {
		
			// new transcipt array
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
	
	String exonID= null;	
	public Exon(Transcript newTranscript, String stableExonID, int start, int end) {

		this.strand= newTranscript.getStrand();
		
			// check for duplicate
		addTranscript(newTranscript);		
					
		setStart(start);
		setEnd(end);
		for (int i = 0; i < getTranscripts().length; i++) 
			getTranscripts()[i].updateBoundaries(this);
		
		
			// decompose ID
		this.exonID= stableExonID;
		setID("exon");
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
		return transcripts[0].getGene();
	}
	
	/**
	 * @return
	 */
	public Transcript[] getTranscripts() {
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
		return getExonID();
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
	/**
	 * @return Returns the homologs.
	 */
	public HashMap getHomologs() {
		return homologs;
	}
	
	public SuperExon getSuperExon() {
		return superExon;
	}
	public void setSuperExon(SuperExon superExon) {
		this.superExon = superExon;
	}
	public SpliceSite getAcceptor() {
		return acceptor;
	}
	public void setAcceptor(SpliceSite acceptor) {
		this.acceptor = acceptor;
	}
	public SpliceSite getDonor() {
		return donor;
	}
	public void setDonor(SpliceSite donor) {
		this.donor = donor;
	}
	
	public AbstractSite getEndSite() {
		return getGene().getSite(getEnd());
	}

	public AbstractSite getStartSite() {
		return getGene().getSite(getStart());
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

}
