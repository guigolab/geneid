package gphase.model;

import gphase.Constants;
import gphase.algo.AlignmentGenerator;

import java.io.File;
import java.io.Serializable;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.plaf.SplitPaneUI;

import com.p6spy.engine.logging.appender.StdoutLogger;
import com.sun.corba.se.impl.ior.OldPOAObjectKeyTemplate;
import com.sun.corba.se.spi.ior.Identifiable;

import qalign.algo.CancelException;
import qalign.algo.CostTable;
import qalign.algo.msa.QAlign;

/*
 * Created on Mar 3, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * 
 * 
 * @author micha
 */
public class Gene extends DirectedRegion {

	String geneID_triv= null;	// trivial name for HOX genes in GenoScan tetraodon
	SpliceSite[] spliceSites= null; // ordered according to position
	ASMultiVariation[] asComplexes= null;
	AbstractSite[] sites= null;
	
	public boolean isProteinCoding() {
		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) 
			if (!transcripts[i].isNonCoding())
				// transcripts[i].getTranslations()!= null&& transcripts[i].getTranslations().length> 0)
				return true;
		return false;
	}
	
	public int getMinCDSStart() {
		int min= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get5PrimeEdge())< Math.abs(min))
					min= tl[j].get5PrimeEdge();
			}
		}
		if (min== Integer.MAX_VALUE)
			return 0;
		return min;
	}
	
	public int getMinCDSEnd() {
		int min= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get3PrimeEdge())< Math.abs(min))
					min= tl[j].get3PrimeEdge();
			}
		}
		if (min== Integer.MAX_VALUE)
			return 0;
		return min;
	}
	

	
	public Transcript[][] recluster() {
		Arrays.sort(transcripts, new AbstractRegion.PositionComparator());
		Vector v= null;
		int max= Integer.MIN_VALUE;
		Vector clusters= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			if (Math.abs(transcripts[i].getStart())> Math.abs(max)) {
				if (v!= null)
					clusters.add(v);
				v= new Vector();
			} 
			v.add(transcripts[i]);
			if (Math.abs(transcripts[i].getEnd())> Math.abs(max))
				max= transcripts[i].getEnd();
		}
		if (v!= null)
			clusters.add(v);
		
		return (Transcript[][]) gphase.tools.Arrays.toField(clusters);
	}
	
	public SpliceSite[] getSpliceSites(int regionType, int type) {
		DirectedRegion[] regions= new DirectedRegion[1];
		if (regionType== AbstractRegion.REGION_5UTR)
			regions[0]= getReal5UTR();
		else if (regionType== AbstractRegion.REGION_CDS)
			regions[0]= getRealCDS();
		else if (regionType== AbstractRegion.REGION_3UTR)
			regions[0]= getReal3UTR();
		else if (regionType== AbstractRegion.REGION_COMPLETE_CLUSTER)
			regions[0]= new DefaultDirectedRegion(getStart(), getEnd(), getStrand());
		
		if (regions[0]== null) {
			return null;
		}
		return getSpliceSites(regions, type);
	}
	
	public SpliceSite[] getSpliceSites(DirectedRegion[] target, int type) {
		Vector v= new Vector();
		for (int i = 0; spliceSites!= null&& i < spliceSites.length; i++) {
			if ((type== SpliceSite.ALTERNATE_SS&& spliceSites[i].isConstitutive())
				|| (type== SpliceSite.CONSTITUTIVE_SS&& !spliceSites[i].isConstitutive()))
				continue;
			for (int j = 0; j < target.length; j++) {
				
				if (spliceSites[i].getPos()< target[j].get5PrimeEdge()|| 
						spliceSites[i].getPos()> target[j].get3PrimeEdge())
					continue;
				v.add(spliceSites[i]);
				break;
			}
		}
		return (SpliceSite[]) gphase.tools.Arrays.toField(v);
	} 
	
	public DirectedRegion getReal5UTR() {
		if (isForward()) {
			int x= getMinCDSStart();
			if (x== 0)
				return null;
			return new DefaultDirectedRegion(getStart(), x- 1, getStrand());	// utr starts outside CDS			
		} else {
			int x= getMaxCDSStart();
			if (x== 0)
				return null;
			return new DefaultDirectedRegion(x- 1, getEnd(), getStrand());
		}
		
	}
	
	public DirectedRegion getRealCDS() {
		
		if (isForward())
			return new DefaultDirectedRegion(getMaxCDSStart(), getMinCDSEnd(), getStrand());
		else
			return new DefaultDirectedRegion(getMaxCDSEnd(), getMinCDSStart(), getStrand());
	}
	
	public DirectedRegion getReal3UTR() {
		
		if (isForward()) {
			int x= getMaxCDSEnd();
			if (x== 0)
				return null;
			return new DefaultDirectedRegion(x+1, getEnd(), getStrand());	// utr starts outside CDS			
		} else {
			int x= getMinCDSEnd();
			if (x== 0)
				return null;
			return new DefaultDirectedRegion(getStart(), x+ 1, getStrand());	// neg strand -(-1)
		}
		
		
	}
	
	public int getMaxCDSEnd() {
		int max= 0;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get3PrimeEdge())> Math.abs(max))
					max= tl[j].get3PrimeEdge();
			}
		}
		if (max== 0)
			return 0;
		return max;
	}

	public int getMaxCDSStart() {
		int max= 0;
		for (int i = 0; i < transcripts.length; i++) {
			Translation[] tl= transcripts[i].getTranslations();
			if (tl== null)
				continue;
			for (int j = 0; j < tl.length; j++) {
				if (Math.abs(tl[j].get5PrimeEdge())> Math.abs(max))
					max= tl[j].get5PrimeEdge();
			}
		}
		if (max== 0)
			return 0;
		return max;
	}
	
	public AbstractRegion getCDSRegion() {
		return new DefaultRegion(getMinCDSStart(), getMaxCDSEnd());
	}
	
	/**
	 * compares arbitrary two instances of String or Gene 
	 * 
	 * 
	 * @author msammeth
	 */
	public static class StableIDComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			String sID1= null;
			try {
				sID1= ((Gene) arg0).getStableID();
			} catch (ClassCastException e) {
				sID1= (String) arg0;
			}
			
			String sID2= null;
			try {
				sID2= ((Gene) arg1).getStableID();
			} catch (ClassCastException e) {
				sID2= (String) arg1;
			}
			
			return sID1.compareToIgnoreCase(sID2);
		}
	}
	
	final static int getArrayPosition(String query, String[] array) {
 
		if (query== null)
			return -1;
		
		int i;
		for (i = 0; i< array.length; i++) 
			if (query.equalsIgnoreCase(array[i]))
				break;
		
		if (i>= array.length) 
			return -1;
		
		return i;
	}
	
	int type= Constants.NOINIT; 
	String confidence= null;


	/**
	 * Retrieves <b>all</b> homolog genes for <code>this</code>'s homologies
	 * and eliminates doubles.
	 * @param array
	 * @return
	 */
	public Gene[] getHomologGenes() {
		return getHomologGenes(getHomologies());
	}
	
	public Exon getExon(String exonID) {
		
		if (transcripts== null)
			return null;
		
		Exon e= null;
		for (int i = 0; i < transcripts.length; i++) {
			e= transcripts[i].getExon(exonID);
			if (e!= null)
				return e;
		}
		return e;
	}
	
	public Exon getExon(int start, int end) {
		Exon[] exons= getExons();
		for (int i = 0; i < exons.length; i++) 
			if (exons[i].getStart()== start&& exons[i].getEnd()== end)
				return exons[i];
		return null;
	}
	
	/**
	 * Retrieves homolog genes for the given species
	 * and eliminates doubles.
	 * @param array
	 * @return
	 */
	public Gene[] getHomologGenes(Species refSpec) {
		return getHomologGenes(getHomologies(refSpec));
	}
	
	/**
	 * Retrieves the homolog genes for the homologies given
	 * and eliminates doubles.
	 * @param array
	 * @return
	 */
	public Gene[] getHomologGenes(GeneHomology[] homols) {

		if (homols== null)
			return null; 
		
			//	get genes
		Vector geneV= new Vector(homols.length);		
		for (int i = 0; i < homols.length; i++) 
			geneV.add(homols[i].getOtherGene(this));
		
			// check for doubles (gene conversion/duplication)
		for (int i = 0; i < geneV.size(); i++) 
			for (int j = (i+1); j < geneV.size(); j++) 
				if (geneV.elementAt(i)== geneV.elementAt(j))
					geneV.remove(j--);
		
		return toGeneArray(geneV);
	}

	
	public SpliceSite checkSpliceSite(SpliceSite ss) {
		if (spliceSites== null)
			return null;
		
		int p= Arrays.binarySearch(spliceSites, ss, new SpliceSite.PositionComparator());	// no abstract site, check don/acc !!
		if (p>= 0) {
			if (spliceSites[p].isDonor()!= ss.isDonor()) 
				System.err.println("Mismatching splice site type (donor/acceptor): "+ getStableID()+" pos "+ss.getPos());
			return spliceSites[p]; 	// bullshit, add first and remove afterwards when filtering!
		}
		return null;
	}
	
	public boolean addSpliceSite(SpliceSite ss) {
		
		if (checkSpliceSite(ss)!= null)
			return false;
		
			// insert new splice site
		if(spliceSites== null) 
			spliceSites= new SpliceSite[] {ss};
		else {
			int p= Arrays.binarySearch(spliceSites, ss, new SpliceSite.PositionComparator());
			if (p< 0)
				spliceSites= (SpliceSite[]) gphase.tools.Arrays.insert(spliceSites, ss, p);
		}
		return true;
	}

	/**
	 * @deprecated no longer in use
	 *
	 */
	public void initSpliceSites() {
		
		Exon[] exons= getExons();
		if (exons== null)
			return;
		
		for (int i = 0; i < exons.length; i++) {
			SpliceSite acceptor= new SpliceSite(this, exons[i].getStart(), false, exons[i]);
			SpliceSite donor= new SpliceSite(this, exons[i].getEnd(), true, exons[i]);
		
			if (spliceSites== null)
				spliceSites= new SpliceSite[] {acceptor, donor};
			else {
				int x= Arrays.binarySearch(spliceSites, acceptor, new SpliceSite.PositionComparator());
				if (x< 0) {	// insert
					x=(x+ 1)* -1;
					SpliceSite[] ss= new SpliceSite[spliceSites.length+ 1];
					for (int j = 0; j < x; j++) 
						ss[j]= spliceSites[j];
					ss[x]= acceptor;
					for (int j = x+1; j < ss.length; j++) 
						ss[j]= spliceSites[j-1];
				} else
					acceptor= spliceSites[x];
				
				x= Arrays.binarySearch(spliceSites, donor, new SpliceSite.PositionComparator());
				if (x< 0) {	// insert
					x=(x+ 1)* -1;
					SpliceSite[] ss= new SpliceSite[spliceSites.length+ 1];
					for (int j = 0; j < x; j++) 
						ss[j]= spliceSites[j];
					ss[x]= donor;
					for (int j = x+1; j < ss.length; j++) 
						ss[j]= spliceSites[j-1];
				} else
					donor= spliceSites[x];
			}
			
			exons[i].setAcceptor(acceptor);
			exons[i].setDonor(donor);			
		}
		
	}
	
	/**
	 * 
	 * @param newConfidence confidence descriptor
	 * @return <code>false</code> if confidence is already set or the confidence
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public void setConfidence(String newConfd) {
		this.confidence= newConfd;
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
	
	/**
	 * 
	 * @param newType type descriptor
	 * @return <code>false</code> if the type is not set or the type
	 * does not equal to the descriptor, <code>true</code> otherwise.  
	 */
	public boolean isType(String matchType) {
		
		if (type== Constants.NOINIT)
			return false;
		
		return (Constants.TYPES[type].equalsIgnoreCase(matchType));
	}
	
	/**
	 * 
	 * @param newConfd confidence descriptor
	 * @return <code>false</code> if confidence is not set or the confidence
	 * does not equal to the descriptor, <code>true</code> otherwise.  
	 */
	public boolean isConfidence(String matchType) {
		
		if (type== Constants.NOINIT)
			return false;
		
		return (Constants.CONFIDENCES[type].equalsIgnoreCase(matchType));
	}	
	
	
	public final static Gene[] toGeneArray(Vector v) {

		if (v== null)
			return null;
		return toGeneArray(v.toArray());
	}
	
	public final static Gene[] toGeneArray(Object[] o) {

		Gene[] result= new Gene[o.length];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) o[i];
		return result;
	}
		public boolean isExonicPosition(int absPos) {
			
			Exon[] exons= getExons();
			for (int i = 0; i < exons.length; i++) 
				if (exons[i].contains(absPos))
					return true;
			
			return false;
		}
		
		// all exons identical, one/some missing		
		Exon[][] filterAlternativeExons(Transcript t1, Transcript t2) {
			
				// get& check
			Exon[] e1= t1.getExons();
			Exon[] e2= t2.getExons();
			if (e1== null|| e2== null|| e1.length< 1|| e2.length< 1)
				return null;
			
				// compare
			Vector e1Vec= new Vector();
			for (int i = 0; i < e1.length; i++) 
				e1Vec.add(e1[i]);
			Vector e2Vec= new Vector();
			for (int i = 0; i < e2.length; i++) 
				e2Vec.add(e2[i]);
			for (int i = 0; i < e1Vec.size(); i++)	// remove identical exon pairs 
				for (int j = 0; i < e2Vec.size(); j++) 
					if (e1Vec.elementAt(i).equals(e2Vec.elementAt(j))) {
						e1Vec.remove(i);
						e2Vec.remove(j);
						--i; --j;
					}
			
			Exon[][] concat= new Exon[2][];
			concat[0]= new Exon[e1Vec.size()];
			for (int i = 0; i < e1Vec.size(); i++) 
				concat[0][i]= (Exon) e1Vec.elementAt(i);
			for (int i = 0; i < e2Vec.size(); i++) 
				concat[1][i]= (Exon) e2Vec.elementAt(i);
            return concat;
		}
		
		String assembly= null;
	String geneID= null;
	HashMap homologies= null; // species to vector of genes
	
	Transcript[] transcripts= null;

	public Gene(Species spec, String stableGeneID) {

		setSpecies(spec);
		
		geneID= stableGeneID;
	}
	
	/**
	 * Checks for duplicates and adds gene to the list of homologs.
	 * 
	 * @param homol the homolog <code>Gene</code> to add
	 * @return <code>false</code> if a gene with same stable ID already has been added,
	 * <code>true</code> otherwise
	 */
	public boolean addHomology(GeneHomology homol) {
		
		if (homologies== null) {
			homologies= new HashMap();
		}
		
			// else
		Vector homs= (Vector) homologies.get(homol.getOtherGene(this).getSpecies());
		if (homs== null)
			homs= new Vector();
		else 				
			for (int i = 0; i < homs.size(); i++) {		// check whether exists
				if (((GeneHomology) homs.elementAt(i)).getOtherGene(this).getStableID().equals(
						homol.getOtherGene(this).getStableID()))
					return false;
			}
		
		homs.add(homol);	// add if not already in there
		
		homologies.put(homol.getOtherGene(this).getSpecies(), homs);
		return true;
	}
	
	public String toString() {
		return getStableID();
	}

	public void removeSpliceSite(SpliceSite ss) {
		
		if (spliceSites== null|| spliceSites.length< 1)
			return;
		
		SpliceSite[] newSpliceSites= new SpliceSite[spliceSites.length- 1];
		int pos= 0;
		for (int i = 0; i < spliceSites.length; i++) {
			if (spliceSites[i]!= ss)
				newSpliceSites[pos++]= spliceSites[i];
		}
		spliceSites= newSpliceSites;
	}
	
	public boolean removeTranscript(Transcript trans) {
		
		//System.out.print("Removing transcript "+trans.getStableID()+":");
			// remove from shared exons
		for (int i = 0; i < trans.getExons().length; i++) {
			Exon ex= trans.getExons()[i]; 
			ex.removeTranscript(trans);
			if (ex.getTranscripts().length== 0) {	// remove exon/ss
				if (ex.getDonor()!= null) {
					ex.getDonor().removeExon(ex);
					if (ex.getDonor().getExons()== null|| ex.getDonor().getExons().length< 1)
						removeSpliceSite(ex.getDonor());
				} 
				if (ex.getAcceptor()!= null) {
					ex.getAcceptor().removeExon(ex);
					if (ex.getAcceptor().getExons()== null|| ex.getAcceptor().getExons().length< 1)
						removeSpliceSite(ex.getAcceptor());
				}
			}
		}
		
			// remove from gene
		Transcript[] newTranscripts= new Transcript[transcripts.length- 1];
		int pos= 0;
		boolean flag= false;
		for (int i = 0; i < transcripts.length; i++) 
			if (transcripts[i]!= trans)
				newTranscripts[pos++]= transcripts[i];
			else
				flag= true;
		
		if (flag) {
			transcripts= newTranscripts;
			//System.out.println(" ok.");
		} else 
			;//System.out.println(" failed!");
		
			// update region
		int minStart= Integer.MAX_VALUE, maxEnd= 0;
		for (int i = 0; i < transcripts.length; i++) {
			if (Math.abs(transcripts[i].getStart())< Math.abs(minStart))
				minStart= transcripts[i].getStart();
			if (Math.abs(transcripts[i].getEnd())> Math.abs(maxEnd))
				maxEnd= transcripts[i].getEnd();
		}
		setStart(minStart);
		setEnd(maxEnd);
		
		return flag;
	}
	
	public boolean removeHomology(Gene hgene) {
		
		if (homologies== null|| homologies.get(hgene.getSpecies())== null)
			return false;
		
		Vector v= (Vector) homologies.get(hgene.getSpecies());
		for (int i = 0; i < homologies.size(); i++) 
			if (((GeneHomology) v.elementAt(i)).getOtherGene(this).equals(hgene)) {
				v.remove(i);
				return true;
			}				
			
		return false;
	}
	
	/**
	 * @return
	 */
	public String getAssembly() {
		return assembly;
	}

	/**
	 * @return
	 */
	public Exon[] getExons() {

		Vector v= new Vector();
		for (int i = 0; i < transcripts.length; i++) 
			for (int j = 0; transcripts[i].getExons()!= null&&
					j < transcripts[i].getExons().length; j++) 
				v.add(transcripts[i].getExons()[j]);
			
		Exon[] exons= new Exon[v.size()];
		for (int i = 0; i < v.size(); i++) 
			exons[i]= (Exon) v.elementAt(i);
		
		return exons;
	}
	
	/**
	 * 
	 * @param seqs	return by parameter !!!
	 * @return
	 * @deprecated too mem-intensive
	 */
	int align_qalign(String[] seqs) {
		
		QAlign qalign= new QAlign();
		try {
			qalign.setAll(
			    0,	// weighting tree
			    1,	// output console
			    1,	// simultaneous alignment
			    seqs,
			    CostTable.DNARNA,
			    5,	// minimal epsilon
			    null,
			    null,
			    null);
			qalign.run();
			seqs[0]= qalign.getSimultaneousLayout()[0];
			seqs[1]= qalign.getSimultaneousLayout()[1];
		} catch (CancelException e) {
			e.printStackTrace();
		}
		
		return qalign.getCost();
	}
	
	/**
	 * @return
	 */
	public String getGeneID() {
		return geneID;
	}
	
	public AbstractSite[] getSites() {

		if (sites == null) {
				// assemble signals (SS, TSS, TES) 
			SpliceSite[] sUniverse= getSpliceSites();
			Vector v= new Vector(sUniverse.length);
			Comparator compi= new AbstractSite.PositionComparator();
			for (int i = 0; i < sUniverse.length; i++) 
				v.add(sUniverse[i]);
			
			for (int i = 0; i < getTranscripts().length; i++) {	// get tss
				AbstractSite tss= getTranscripts()[i].getTSS();
				int j;
				for (j = 0; j < v.size(); j++) 
					if (compi.compare(v.elementAt(j), tss)== 0)
						break;
				if (j>= v.size()) 
					v.add(tss);
			}
			for (int i = 0; i < getTranscripts().length; i++) {	// get tes
				AbstractSite tes= getTranscripts()[i].getTES();
				int j;
				for (j = 0; j < v.size(); j++) 
					if (compi.compare(v.elementAt(j), tes)== 0)
						break;
				if (j>= v.size()) 
					v.add(tes);
			}
			
			sites= new AbstractSite[v.size()];
			for (int i = 0; i < sites.length; i++) 
				sites[i]= (AbstractSite) v.elementAt(i);
			Arrays.sort(sites, compi);
		}

		return sites;
	}
	
	public AbstractSite getSite(int pos) {
		AbstractSite s= null;
		AbstractSite o= new AbstractSite();
		o.setPos(pos);
		int b= Arrays.binarySearch(getSites(), o, new AbstractSite.PositionComparator());
		if (b>= 0)
			s= getSites()[b];
		
		return s;
	}

	/**
	 * Gets <b>all</b> homologies, acroos all other species.
	 * @return
	 */
	public GeneHomology[] getHomologies() {
		
		if (homologies== null)
			return null;
		Collection v= homologies.values();
		if (v== null)
			return null;
		
		Vector result= new Vector();
		Iterator iter= v.iterator();
		while (iter.hasNext())
			result.addAll((Vector) iter.next());
		
		return GeneHomology.toGeneHomologyArray(result);
	}

	/**
	 * @return
	 */
	public Transcript[] getTranscripts() {
		return transcripts;
	}
	
	public int getTranscriptCount() {
		if (transcripts== null)
			return 0;
		return transcripts.length;
	}
	
	public ASVariation[] getASVariations(int codingCode) {
		ASMultiVariation[] asm= getASMultiVariations();
		if (asm== null)
			return null;
		Vector resVec= new Vector(asm.length);
		for (int i = 0; i < asm.length; i++) {
			ASVariation[] as= null;
			if (codingCode== ASMultiVariation.FILTER_NONE)
				as= asm[i].getASVariationsAll();
			else if (codingCode== ASMultiVariation.FILTER_CODING_REDUNDANT)
				as= asm[i].getASVariationsHierarchicallyFiltered();
			else if (codingCode== ASMultiVariation.FILTER_HIERARCHICALLY)
				as= asm[i].getASVariationsClusteredCoding();
			for (int j = 0; j < as.length; j++) 
				resVec.add(as[j]);
		}
			
		return (ASVariation[]) gphase.tools.Arrays.toField(resVec);
	}
	
	/**
	 * @deprecated does not work that easy for constitutive exons, 
	 * bette take non-involvement in ASVariation as criterion
	 * @param constitutive
	 * @return
	 */
	public Exon[] getExons(boolean constitutive) {
		
		Vector v= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			for (int j = 0; j < transcripts.length; j++) {
				Exon[] ex1= transcripts[i].getExons();
				for (int k = 1; k < ex1.length- 1; k++) {	// no terminal exons for the comparison
					if (!transcripts[j].contains(ex1[k]))
						continue;
					Exon[] ex2= transcripts[j].getExons();
					int m;
					for (m = 1; m < ex2.length- 2; m++) {
						if (ex1[k].overlaps(ex2[m])) {
							if (ex1[k].get)
						}
					}
					if (m== ex2.length&& !constitutive)
						v.add(ex1);
				}
			}
		}
	}
	
	public ASMultiVariation[] getASMultiVariations() {
		
		asComplexes= null;
		if (asComplexes == null) {
			if (transcripts== null|| transcripts.length< 2) 
				return null;			
			
			SpliceSite[][] spliceChains= new SpliceSite[transcripts.length][];
			for (int i = 0; i < transcripts.length; i++) 
				spliceChains[i]= transcripts[i].getSpliceChain();
			
			Vector spliceClusters= tokenizeASClusters(spliceChains, false, false);	// get splice clusters across all sequences
																			//TODO warning! TSS/TES flanked events included !!
			
				// determine pw splice vars
			Vector asComp= new Vector(spliceClusters.size());
			int[] rightEdge= new int[transcripts.length];
			for (int x = 0; x < spliceClusters.size(); x++) {	// iterate complexes
				
				SpliceSite[][] cluster= (SpliceSite[][]) spliceClusters.elementAt(x);
				
				Vector as2Events= new Vector();
				for (int i = 0; i < transcripts.length; i++) {		// iterate pw combinations for a complex 
					for (int j = (i+1); j < transcripts.length; j++) {
						SpliceSite[][] ss2= new SpliceSite[2][];
						ss2[0]= cluster[i];
						ss2[1]= cluster[j];

						Transcript[] tt= new Transcript[2];
						tt[0]= transcripts[i];
						tt[1]= transcripts[j];
						
						ss2= ASVariation.trim(ss2, tt);	// trim left in first, last at right
						
						if (ss2[0].length< 1&& ss2[1].length< 1)	// skip when both are empty
							continue;
						
						Vector cc2= tokenizeASClusters(ss2, false, false);
						for (int k = 0; k < cc2.size(); k++) {			// maybe more than one AS variation for a pair
							SpliceSite[][] pwEvent= (SpliceSite[][]) cc2.elementAt(k);
							ASVariation as2= new ASVariation(transcripts[i], transcripts[j],
									pwEvent[0], pwEvent[1]);
							as2Events.add(as2);
						}
					}
				}
				
				ASVariation[] as2Evs= new ASVariation[as2Events.size()];
				for (int i = 0; i < as2Events.size(); i++) 
					as2Evs[i]= (ASVariation) as2Events.elementAt(i);
				if (as2Evs.length!= 0)		// no pw var without a TSS/TES in multi-cluster
					asComp.add(new ASMultiVariation(as2Evs));
			}
			
			asComplexes = new ASMultiVariation[asComp.size()];
			for (int i = 0; i < asComplexes.length; i++) 
				asComplexes[i]= (ASMultiVariation) asComp.elementAt(i);
		}

		return asComplexes;
	}
	
	
	Vector tokenizeASClusters(SpliceSite[][] spliceChains, boolean omitStart, boolean omitEnd) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionComparator();	// here! SpliceSite.PositionComparator
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			SpliceSite s= spliceChains[0][i];
			int[] pos= new int[spliceChains.length];
			pos[0]= i;
			int j;
			for (j = 1; j < spliceChains.length; j++) {
				int p= Arrays.binarySearch(spliceChains[j], s, c);
				if (p< 0)
					break;	// as soon as not found in one schain, give up..
				else
					pos[j]= p;
			}
			if (j>= spliceChains.length)	// found
				tokens.add(pos);
		}
		
			// tokenize by these conserved splice sites
		SpliceSite[][] ss= new SpliceSite[spliceChains.length][];
		int[] posOld= new int[spliceChains.length];
		for (int i = 0; i < posOld.length; i++) 
			posOld[i]= (-1);	// anchor before sequence start	
		int[] pos;
		Vector result= new Vector();
		for (int i = 0; i < tokens.size(); i++) {
			pos= (int[]) tokens.elementAt(i);
			if (omitStart&& i== 0) {	// omit splice chains containing start
				posOld= pos;
				continue;
			}
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) {
					result.add(ass);
					break;
				}
			posOld= pos;
		}
		
		if (!omitEnd) {
			pos= new int[spliceChains.length];		// anchor after sequence end
			for (int i = 0; i < pos.length; i++) 
				pos[i]= spliceChains[i].length;
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				ass[j]= new SpliceSite[pos[j]- posOld[j]- 1];
				for (int k = 0; k < (pos[j]- posOld[j]- 1); k++) 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) {
					result.add(ass);
					break;
				}
		}
		
		return result;
	}

	/**
	 * @param string
	 */
	public void setAssembly(String string) {
		assembly= string;
	}

	/**
	 * @param i
	 */
	public void setGeneID(String i) {
		geneID= i; 
	}

	/**
	 * @param transcripts
	 */
	public boolean addTranscript(Transcript newTranscript) {
		
			// update boundaries
		if (Math.abs(newTranscript.getStart())< Math.abs(getStart()))
			setStart(newTranscript.getStart());
		if (Math.abs(newTranscript.getEnd())> Math.abs(getEnd()))
			setEnd(newTranscript.getEnd());
		
		if(transcripts== null) {
			transcripts= new Transcript[] {newTranscript};
			return true;
		}
		
			// search transcript
		for (int i = 0; i < transcripts.length; i++) 
			if (transcripts[i].getStableID().equalsIgnoreCase(newTranscript.getStableID()))
				return false;

			// add transcript		
		Transcript[] nTranscripts= new Transcript[transcripts.length+ 1];
		for (int i= 0; i < transcripts.length; i++) 
			nTranscripts[i]= transcripts[i];
		nTranscripts[nTranscripts.length- 1]= newTranscript;
		transcripts= nTranscripts;
		
		sites= null;	// are to be re-inited
		
		return true;
	}

	/**
	 * @param transcripts
	 */
	public void setTranscripts(Transcript[] transcripts) {
		this.transcripts= transcripts;
	}
	
	public static final String getSpeciesPfx(String ensemblString) {
		
		Pattern patty= Pattern.compile("^(\\D+)\\d+$");
		Matcher matty= patty.matcher(ensemblString);
		
		if (matty.matches())
			return matty.group(1).substring(0, matty.group(1).length()- 1);	// chop marker ("G", "T", "E",..)
		return null;	// error
	}
	
	public static final String getStableID(String pfx, char marker, int nb) {
		
		StringBuffer result= new StringBuffer(""+ nb); // the int nb
		
			// add leading '0's
		char[] c= new char[11- result.length()];
		Arrays.fill(c, '0');
		result.insert(0,c);		
		
		result.insert(0, marker);	// add marker (gene, transcript or exon)
		result.insert(0, pfx);	// add species pfx
		
		return result.toString();
	}
	
	/**
	 * @return
	 */
	public String getStableID() {
		
		return geneID;
	}	

	static final long serialVersionUID = 8933737248601221991L;


	/**
	 * @param species2
	 * @return
	 */
	public GeneHomology[] getHomologies(Species refSpec) {
		
		if ((homologies== null)|| (homologies.get(refSpec)== null))
			return null;
		
		return GeneHomology.toGeneHomologyArray((Vector) homologies.get(refSpec));
	}
	
	public String printGeneStructure() {
		
		String result= getStart()+ "----";
		for (int i = 0; i < spliceSites.length; i++) {
			if (!spliceSites[i].isDonor())
				result+= "<";
			result+= spliceSites[i].getPos();
			if (spliceSites[i].isDonor())
				result+= ">";
			result+= "---";
		}
		result+= getEnd();
		
		return result;
	}

	/**
	 * Inspects whether there is a homology with one of the given species
	 * @deprecated not recomended
	 * @param rr
	 * @return
	 */
	public boolean hasHomologyWith(Species[] hSpec) {

		if (homologies== null) {
			return false;
		}
		
		for (int i = 0; i < hSpec.length; i++) 
			if (homologies.get(hSpec[i])== null|| ((Vector) homologies.get(hSpec[i])).size()< 1)
				return false;
			
		return true; 
	}

	/**
	 * @return Returns the confidence.
	 */
	public String getConfidence() {
		return confidence;
	}
	/**
	 * @return Returns the type.
	 */
	public int getType() {
		return type;
	}

	String chromosome = null;
	Species species = null;

	/**
	 * @return
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @return
	 */
	public Species getSpecies() {
		return species;
	}

	/**
	 * @param string
	 */
	public void setChromosome(String string) {
		chromosome= string;
	}

	/**
	 * @param string
	 */
	public void setSpecies(Species newSpecies) {
		species= newSpecies;
	}

	public SpliceSite[] getSpliceSites() {
		return spliceSites;
	}

	private static int uniqueID= 0;
	
	public static String getUniqueID() {
		
		String nb= Integer.toString(uniqueID++);
		for (int i = 0; i < 13- nb.length(); i++) 
			nb= "0"+ nb;
		
		return nb;
	}

	
}

