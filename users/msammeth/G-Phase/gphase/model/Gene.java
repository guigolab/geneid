package gphase.model;

import gphase.Constants;
import gphase.algo.AlignmentGenerator;
import gphase.graph.Tuple;
import gphase.tools.Array;

import java.awt.List;
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
	SpliceSite[] tssSites= null, tesSites= null;
	AbstractSite[] abstractSites= null; // ordered according to position
	ASMultiVariation[] asComplexes= null;
	AbstractSite[] sites= null;
	TU[] tu= null;
	Exon[] exons= null;
	Intron[] introns= null;
	ASVariation[] vars= null;
	int varCC= -1;
	
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
	
	/**
	 * Returns Map(non-red event x Vector[](alt tIDs))
	 * @return
	 */
	public HashMap getVariantGroups() {
		ASVariation[] varsNorm= getASVariations(ASMultiVariation.FILTER_NONE);
		ASVariationGroup[] vars= null;
		if (varsNorm!= null) {
			vars= new ASVariationGroup[varsNorm.length];
			for (int i = 0; i < vars.length; i++) 
				vars[i]= new ASVariationGroup(varsNorm[i]);
		}
		HashMap groupsHash= new HashMap();
		Comparator compi= new ASVariation.StructureComparator();
		for (int j = 0; vars!= null&& j < vars.length; j++) {
			Vector[] v= (Vector[]) groupsHash.remove(vars[j]);
			if (v== null) {
				v= new Vector[2];
				v[0]= new Vector();
				v[1]= new Vector();
			}
			v[0]= (Vector) gphase.tools.Arrays.addUnique(v[0], vars[j].getTranscript1().getTranscriptID());
			v[1]= (Vector) gphase.tools.Arrays.addUnique(v[1], vars[j].getTranscript2().getTranscriptID());
			groupsHash.put(vars[j], v);
		}
		
		return groupsHash;
	}
	

	/**
	 * @return
	 */
	public VariantGroup_SpliceChain[] getVariantGroups_spliceChains() {
		ASVariation[] vars= getASVariations(ASMultiVariation.FILTER_NONE);
		HashMap variantMap= new HashMap(); // SC x Variant Group
		for (int j = 0; vars!= null&& j < vars.length; j++) {
			SpliceSite[] sc1= vars[j].getSpliceChain1();
			SpliceSite[] sc2= vars[j].getSpliceChain2();
			
			if (sc1== null) {	// if empty, invert
				sc1= new SpliceSite[sc2.length];
				for (int i = sc2.length- 1; i >= 0; --i) 
					sc1[i- sc2.length+ 1]= sc2[i];
			}
			Array sa1= new Array(sc1);
			VariantGroup_SpliceChain var= (VariantGroup_SpliceChain) variantMap.get(sa1); 
			if (var== null) {
				var= new VariantGroup_SpliceChain(sc1);
				var.addTranscriptID(vars[j].getTranscript1().getTranscriptID());
				var.addASevent(vars[j]);
				variantMap.put(sa1, var);
			} else {
				var.addTranscriptID(vars[j].getTranscript1().getTranscriptID());
				var.addASevent(vars[j]);
			}
			
			if (sc2== null) {	// if empty, invert
				sc2= new SpliceSite[sc1.length];
				for (int i = sc1.length- 1; i >= 0; --i) 
					sc2[i- sc1.length+ 1]= sc1[i];
			}
			Array sa2= new Array(sc2);
			var= (VariantGroup_SpliceChain) variantMap.get(sa2); 
			if (var== null) {
				var= new VariantGroup_SpliceChain(sc2);
				var.addTranscriptID(vars[j].getTranscript2().getTranscriptID());
				var.addASevent(vars[j]);
				variantMap.put(sa2, var);
			} else {
				var.addTranscriptID(vars[j].getTranscript2().getTranscriptID());
				var.addASevent(vars[j]);
			}
		}

			// result
		Object[] vals= variantMap.values().toArray();
		return (VariantGroup_SpliceChain[]) gphase.tools.Arrays.toField(vals);
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
	

	
	public void initTU() {
		Vector codV= new Vector();
		Vector ncV= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].isCoding())
				codV.add(transcripts[i]);
			else
				ncV.add(transcripts[i]);
		}
		
			// cluster first coding ones
		for (int i = 0; i < codV.size(); i++) {
			Transcript tr= ((Transcript) codV.elementAt(i));
			Translation t= tr.getTranslations()[0];
			int j;
			for (j = 0; tu!= null&& j < tu.length; j++) {
				if (tu[j].overlaps(t)) {	// overlapping CDS
					SpliceSite[] ss= tr.getSpliceChain();
					int k;
					for (k = 0; k < ss.length; k++)	// mind 1 overlapping SS 
						if (tu[j].contains(ss[k]))
							break;
					if (k< ss.length)
						break;
				}
			}
			if (tu== null)
				tu= new TU[] {new TU(this, tr)};
			else {
				if(j< tu.length) 
					tu[j].addTranscript(tr);
				else
					tu= (TU[]) gphase.tools.Arrays.add(tu, new TU(this, tr));
			}
		}
		
			// then join nc transcripts
		for (int i = 0; i < ncV.size(); i++) {
			Transcript tr= ((Transcript) ncV.elementAt(i));
			for (int j = 0; tu!= null&& j < tu.length; j++) {
				SpliceSite[] ss= tr.getSpliceChain();
				int k;
				for (k = 0; k < ss.length; k++)	// mind 1 overlapping SS 
					if (tu[j].contains(ss[k]))
						break;
				if (k< ss.length)
					tu[j].addTranscript(tr);
			}
		}
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
	
	public SpliceSite[] getSpliceSites(int spliceType) {

		if (spliceSites== null)
			return null;
		
		Vector v= new Vector(spliceSites.length);
		for (int i = 0; i < spliceSites.length; i++) {
			if (spliceType== SpliceSite.ALTERNATE_SS&& !spliceSites[i].isConstitutive()) 
				v.add(spliceSites[i]);
			else if (spliceType== SpliceSite.CONSTITUTIVE_SS&& spliceSites[i].isConstitutive()) 
				v.add(spliceSites[i]);
			else if (spliceType!= SpliceSite.CONSTITUTIVE_SS&& spliceType!= SpliceSite.ALTERNATE_SS)
				v.add(spliceSites[i]);
		}
		
		return (SpliceSite[]) gphase.tools.Arrays.toField(v);
	}
	
	public SpliceSite[] getSpliceSites(int spliceType, int regionType) {
		SpliceSite[] sites= getSpliceSites(spliceType);
		DirectedRegion region= getRegion(regionType);
		return (SpliceSite[]) DirectedRegion.contained(region, sites);
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
		DirectedRegion reg;
		if (!isProteinCoding())
			System.currentTimeMillis();
		if (isForward()) {
			int x= getMinCDSStart();
			if (x== 0)
				return null;
			reg= new DirectedRegion(getStart(), x- 1, getStrand());	// utr starts outside CDS			
		} else {
			int x= getMaxCDSStart();
			if (x== 0)
				return null;
			reg= new DirectedRegion(x- 1, getEnd(), getStrand());
		}
		reg.setChromosome(getChromosome());
		return reg;
		
	}
	
	public DirectedRegion getRealCDS() {
		
		DirectedRegion reg;
		if (isForward())
			reg= new DirectedRegion(getMaxCDSStart(), getMinCDSEnd(), getStrand());
		else 
			reg= new DirectedRegion(getMaxCDSEnd(), getMinCDSStart(), getStrand());
		reg.setChromosome(getChromosome());
		return reg;

	}
	
	public DirectedRegion getMaxCDS() {
		
		DirectedRegion reg;
		if (isForward())
			reg= new DirectedRegion(getMinCDSStart(), getMaxCDSEnd(), getStrand());
		else 
			reg= new DirectedRegion(getMinCDSEnd(), getMaxCDSStart(), getStrand());
		reg.setChromosome(getChromosome());
		return reg;

	} 
	
	public boolean isRealCDS(DirectedRegion reg) {
		if (getRealCDS().contains(reg))
			return true;
		return false;
	}
	
	public boolean isReal5UTR(DirectedRegion reg) {
		if (getReal5UTR().contains(reg))
			return true;
		return false;
	}
	
	public DirectedRegion getReal3UTR() {
		
		DirectedRegion reg;
		if (isForward()) {
			int x= getMaxCDSEnd();
			if (x== 0)
				return null;
			reg= new DirectedRegion(x+1, getEnd(), getStrand());	// utr starts outside CDS			
		} else {
			int x= getMinCDSEnd();
			if (x== 0)
				return null;
			reg= new DirectedRegion(getStart(), x+ 1, getStrand());	// neg strand -(-1)
		}
		
		reg.setChromosome(getChromosome());
		return reg;
		
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

	public AbstractSite addAbstractSite(AbstractSite as) {
		if (abstractSites== null) {
			abstractSites= new AbstractSite[] {as};
			return as;
		}
		
		int p= Arrays.binarySearch(abstractSites, as, new AbstractSite.PositionComparator());
		if (p>= 0) 
			return abstractSites[p];
		else {
			abstractSites= (AbstractSite[]) gphase.tools.Arrays.insert(abstractSites, as, p);
			return as;
		}
	}
	
	public int[] getTranslationInitSites() {
		int[] met= new int[0];
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].isNonCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.isOpenEnded5()|| !getTranscripts()[i].isATGStart())	// exclude truncated 
				continue;
			met= gphase.tools.Arrays.addUnique(met, tln.get5PrimeEdge());
		}
		return met;
	}
	
	public TranslationLocus[] getTranslationLoci() {
		int[] met= getTranslationInitSites();
		Vector v= new Vector();
		for (int i = 0; i < met.length; i++) 
			v.add(getTranscripts(met[i]));
		Transcript[] nc= getNonCodingTranscripts();
		if (nc!= null&& nc.length> 0)
			v.add(nc);
		
		TranslationLocus[] tls= new TranslationLocus[v.size()];
		for (int i = 0; i < tls.length; i++) 
			tls[i]= new TranslationLocus(
					this, (Transcript[]) v.elementAt(i));
		return tls;
	}
	
	public Transcript[] getTranscripts(int tlnInit) {
		Vector v= new Vector();
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].isNonCoding())
				continue;
			Translation tln= getTranscripts()[i].getTranslations()[0];
			if (tln.isOpenEnded5())	// exclude truncated 
				continue;
			if (tln.get5PrimeEdge()== tlnInit)
				v.add(getTranscripts()[i]);
		}
		return ((Transcript[]) gphase.tools.Arrays.toField(v));
	}
	
	public int countCodingSpliceForms() {
		int cnt= 0;
		for (int i = 0; i < getTranscripts().length; i++) 
			if (getTranscripts()[i].isCoding())
				++cnt;
		return cnt;
	}
	
	public Transcript[] getNonCodingTranscripts() {
		Vector v= new Vector();
		for (int i = 0; i < getTranscripts().length; i++) 
			if (getTranscripts()[i].isNonCoding())
				v.add(getTranscripts()[i]);
		return ((Transcript[]) gphase.tools.Arrays.toField(v));
	}
	
	/**
	 * @param ss
	 * @return
	 */
	public boolean addSpliceSite(SpliceSite ss) {
		
//		if (1== 1)
//			return false;
//		
//		if (getSite(ss)!= null)
//			return false;
		
			// insert new splice site
		if(spliceSites== null) 
			spliceSites= new SpliceSite[] {ss};
		else {
			int p= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
			if (p< 0)
				spliceSites= (SpliceSite[]) gphase.tools.Arrays.insert(spliceSites, ss, p);
			else {
				for (int i = 0; i < ss.getExons().length; i++) 
					ss.getExons()[i].replaceSite(ss, spliceSites[p]);
				
				for (int i = 0; i < ss.getTranscripts().length; i++) 
					spliceSites[p].addTranscript(ss.getTranscripts()[i]);
				
			}
		}
		return true;
	}
	
	public boolean replaceSite(SpliceSite oldSite, SpliceSite newSite) {
		for (int i = 0; spliceSites!= null&& i < spliceSites.length; i++) {
			if (spliceSites[i]== oldSite) {
				spliceSites[i]= newSite;
				return true;
			}
		}
		return false;
	}
	
	public AbstractSite getTSS(AbstractSite as) {
		Comparator compi= new AbstractSite.PositionComparator();	// make static
		if (tssSites== null) {	// init
			HashMap<Integer,AbstractSite> map= new HashMap<Integer,AbstractSite>();
			for (int i = 0; i < transcripts.length; i++) {
				Integer pos= new Integer(transcripts[i].getTSSPos());
				AbstractSite site= map.remove(pos);
				if (site== null)
					site= new AbstractSite(pos.intValue());
				site.addTranscript(transcripts[i]);
				map.put(pos, site);
			}
			tssSites= (AbstractSite[]) gphase.tools.Arrays.convertTo(AbstractSite.class, map.values().toArray());
			Arrays.sort(tssSites, compi);
		}
		
		int p= Arrays.binarySearch(tssSites, as, compi);
		assert(p>= 0);
		return tssSites[p];
	}

	public SpliceSite getSite(SpliceSite as) {
		int p= -1;
		if (as.isSpliceSite()) {
			if (getSpliceSites()== null)
				return null;
			p= Arrays.binarySearch(getSpliceSites(), as, SpliceSite.getDefaultPositionTypeComparator());
			if (p>= 0)
				as= getSpliceSites()[p];
		} else if (as.isTSS()) {
			p= Arrays.binarySearch(getTSSsites(), as, SpliceSite.getDefaultPositionTypeComparator());
			if (p>= 0)
				as= getTSSsites()[p];
		} else if (as.isTES()) {
			p= Arrays.binarySearch(getTESsites(), as, SpliceSite.getDefaultPositionTypeComparator());
			if (p>= 0)
				as= getTESsites()[p];
		} else
			System.err.println("ERROR: unknown type "+as.getType());
		
		return as;
	}
	public SpliceSite[] getTESsites() {
		if (tesSites== null) {	// init
			HashMap<Integer,SpliceSite> map= new HashMap<Integer,SpliceSite>();
			for (int i = 0; i < transcripts.length; i++) {
				Integer pos= new Integer(transcripts[i].getTESPos());
				SpliceSite site= map.get(pos);
				if (site== null)
					site= new SpliceSite(pos.intValue(), SpliceSite.TYPE_TES);
				site.addTranscript(transcripts[i]);
				map.put(pos, site);
			}
			tesSites= (SpliceSite[]) gphase.tools.Arrays.convertTo(SpliceSite.class, map.values().toArray());
			Arrays.sort(tesSites, SpliceSite.getDefaultPositionTypeComparator());
		}
		return tesSites;
	}
	
	public SpliceSite[] getTSSsites() {
		if (tssSites== null) {	// init
			HashMap<Integer,SpliceSite> map= new HashMap<Integer,SpliceSite>();
			for (int i = 0; i < transcripts.length; i++) {
				Integer pos= new Integer(transcripts[i].getTSSPos());
				SpliceSite site= map.get(pos);
				if (site== null)
					site= new SpliceSite(pos.intValue(), SpliceSite.TYPE_TSS);
				site.addTranscript(transcripts[i]);
				map.put(pos, site);
			}
			tssSites= (SpliceSite[]) gphase.tools.Arrays.convertTo(SpliceSite.class, map.values().toArray());
			Arrays.sort(tssSites, SpliceSite.getDefaultPositionTypeComparator());
		}
		return tssSites;
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
	
	public boolean isConstitutive(DirectedRegion intron, int eventFilter) {
		DirectedRegion[] intrns= getIntrons();
		for (int i = 0; i < intrns.length; i++) {
			if (!intrns[i].overlaps(intron))
				continue;
			if (intron.get5PrimeEdge()!= intrns[i].get5PrimeEdge()||
					intron.get3PrimeEdge()!= intrns[i].get3PrimeEdge())
				return false;
		}
		return true;
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

	public Gene(String newGeneID) {
		geneID= newGeneID;
		setID("gene");
		//setStrand(getStrand());	// hae?
	}
	
	public Gene(Species spec, String newGeneID) {
		this (newGeneID);
		setSpecies(spec);
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
	
	public String toStringSSPattern() {
		String s= "";
		for (int i = 0; getSpliceSites()!= null&& i < getSpliceSites().length; i++) {
			if (getSpliceSites()[i].isDonor())
				s+= "^";
			else
				s+= "-";
		}
		return s;
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
		for (int i = 0; trans.getExons()!= null&& i < trans.getExons().length; i++) {
			Exon ex= trans.getExons()[i]; 
			ex.removeTranscript(trans);
			if (ex.getDonor()!= null) {
				ex.getDonor().removeTranscript(trans);
				if (ex.getTranscripts().length== 0) 	// remove exon/ss
					ex.getDonor().removeExon(ex);
				if (ex.getDonor().getExons()== null|| ex.getDonor().getExons().length< 1)
					removeSpliceSite(ex.getDonor());
			} 
			if (ex.getAcceptor()!= null) {
				ex.getAcceptor().removeTranscript(trans);
				if (ex.getTranscripts().length== 0) 	// remove exon/ss
					ex.getAcceptor().removeExon(ex);
				if (ex.getAcceptor().getExons()== null|| ex.getAcceptor().getExons().length< 1)
					removeSpliceSite(ex.getAcceptor());
			}
		}
		
			// remove from tu
		for (int i = 0; tu!= null&& i < tu.length; i++) {
			tu[i].removeTranscript(trans);
			if (tu[i].getTranscripts().length< 1)
				tu= (TU[]) gphase.tools.Arrays.remove(tu, tu[i]);
			if (tu.length< 1)
				tu= null;
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
	
	public void repairAlignmentErrors() {
		SpliceSite[] ss= getSpliceSites();
		for (int i = 0; i < ss.length; i++) {
			
		}
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

	public Intron[] getIntrons() {
		if (introns == null) {
			Transcript[] trpts= getTranscripts();
			Vector<Intron> intronV= new Vector<Intron>();
			Comparator compi= new AbstractRegion.PositionComparator();	// position
			for (int i = 0; trpts!= null&& i < trpts.length; i++) {
				Intron[] intrs= trpts[i].getIntrons();
				Vector<Intron> tmpV= new Vector<Intron>();
				for (int j = 0; intrs!= null&& j < intrs.length; j++) {
					int k;
					for (k = 0; k < intronV.size(); k++) {
						if (compi.compare(intronV.elementAt(k), intrs[j])== 0) {
							intronV.elementAt(k).addTranscript(intrs[j].getTranscripts()[0]);
							break;
						}
					}
					if (k== intronV.size())
						tmpV.add(intrs[j]);
				}
				for (int j = 0; j < tmpV.size(); j++) 
					intronV.add(tmpV.elementAt(j));
			}
			introns= (Intron[]) gphase.tools.Arrays.toField(intronV);
		}

		return introns;
	}
	
	public Vector<DirectedRegion> getAdjacentStructures(DirectedRegion reg, boolean prime5, boolean exon) {
		DirectedRegion[] regs= null;
		if (exon)
			regs= getExons();
		else
			regs= getIntrons();
		
		Vector<DirectedRegion> resReg= new Vector<DirectedRegion>();
		for (int i = 0; regs!= null&& i < regs.length; i++) {
			if (prime5&& regs[i].get3PrimeEdge()+ 1== reg.get5PrimeEdge())
				resReg.add(regs[i]);
			else if ((!prime5)&& regs[i].get5PrimeEdge()- 1== reg.get3PrimeEdge())
				resReg.add(regs[i]);
		}
		
		return resReg;
	}
	/**
	 * @return
	 */
	public Exon[] getExons() {

//		Vector v= new Vector();
//		for (int i = 0; i < transcripts.length; i++) 
//			v= (Vector) gphase.tools.Arrays.addUnique(v, transcripts[i].getExons());
//			
//		Exon[] exons= new Exon[v.size()];
//		for (int i = 0; i < v.size(); i++) 
//			exons[i]= (Exon) v.elementAt(i);
//		
		return exons;
	}
	
	public void merge(Gene anotherGene) {
		
		//spliceSites= (SpliceSite[]) gphase.tools.Arrays.addAll(spliceSites, anotherGene.getSpliceSites());
		for (int i = 0; i < anotherGene.getTranscripts().length; i++) 
			addTranscript(anotherGene.getTranscripts()[i]);
		// done in addTranscript()
//		for (int i = 0; i < anotherGene.getTranscripts().length; i++)  
//			updateBoundaries(anotherGene.getTranscripts()[i]);
	}
	
	public DirectedRegion[] getExonicRegions() {
	
		DirectedRegion[] superExons= DirectedRegion.unite((DirectedRegion[]) getExons());
		// if (superExons.length== getExons().length)
		//	System.err.println("No AS");
		return superExons;
	}
	
	public String getExonicRegionsSplicedSequence() {
		
		DirectedRegion[] superExons= getExonicRegions();
		
		String result= "";
		for (int i = 0; i < superExons.length; i++) {
			String tmp= Graph.readSequence(getSpecies(), getChromosome(), isForward(),
					Math.abs(superExons[i].getStart()), Math.abs(superExons[i].getEnd()));
			if (!isForward())
				tmp= gphase.tools.Arrays.reverseComplement(tmp);
			result+= tmp;
		}
		return result;
	}
	
	public int[] getExonicRegionsSSCoords() {
		SpliceSite[] ss= getSpliceSites();
		Comparator compi= new SpliceSite.PositionComparator();
		Arrays.sort(ss, compi);
		
		int[] result= new int[ss.length];
		DirectedRegion[] sExons= getExonicRegions();
		int j= 0;
		int pos= 0;
		for (int i = 0; j< ss.length&& i < sExons.length; i++) {
			while (j< ss.length&& sExons[i].contains(ss[j].getPos())) {	// check for complained SSs
				int aPos= pos+ ss[j].getPos()- sExons[i].get5PrimeEdge();
					// forward/rev does not matter in this case since sequences are already reversed
				result[j++]= aPos;
			}
			pos+= sExons[i].getLength();
		}
	
		return result;
	}
	
	/**
	 * 
	 * @param seqs	return by parameter !!!
	 * @return
	 * @deprecated too mem-intensive
	 */
	int align_qalign(String[] seqs) {
//		
//		QAlign qalign= new QAlign();
//		try {
//			qalign.setAll(
//			    0,	// weighting tree
//			    1,	// output console
//			    1,	// simultaneous alignment
//			    seqs,
//			    CostTable.DNARNA,
//			    5,	// minimal epsilon
//			    null,
//			    null,
//			    null);
//			qalign.run();
//			seqs[0]= qalign.getSimultaneousLayout()[0];
//			seqs[1]= qalign.getSimultaneousLayout()[1];
//		} catch (CancelException e) {
//			e.printStackTrace();
//		}
//		
//		return qalign.getCost();
		return -1;
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
		AbstractSite o= new AbstractSite(pos);
		//o.setPos(pos);
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
	
	public ASVariation[] getASVariations(int codingCode, HashMap refTrptIDs) {
		if (vars== null|| (varCC!= codingCode)) {
			ASMultiVariation[] asm= getASMultiVariations(refTrptIDs);
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
				else if (codingCode== ASMultiVariation.FILTER_STRUCTURALLY)
					as= asm[i].getASVariationsStructurallyFiltered();
				else if (codingCode== ASMultiVariation.FILTER_CONTAINED_IN_CDS) {
					as= asm[i].getASVariationsAll();
					for (int j = 0; j < as.length; j++) {
						if (as[j].isContainedCDS())
							resVec.add(as[j]);
					}
					continue;
				}
				for (int j = 0; as!= null&& j < as.length; j++) 
					resVec.add(as[j]);
			}
			vars= (ASVariation[]) gphase.tools.Arrays.toField(resVec);
			varCC= codingCode;
		}			
		return vars;
	}

	public static ASVariation[] getASVariations(Transcript[] t, HashMap refTrptIDs) {
		
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		Vector result= new Vector();

		for (int x = 0; x < t.length; x++) {
			for (int y = x+1; y < t.length; y++) {
				if (refTrptIDs!= null&& refTrptIDs.size()> 0&&
						refTrptIDs.get(t[x].getTranscriptID())== null&& 
						refTrptIDs.get(t[y].getTranscriptID())== null)
					continue;

				SpliceSite[][] spliceChains= new SpliceSite[2][];
				spliceChains[0]= t[x].getSpliceChainComplete();
				spliceChains[1]= t[y].getSpliceChainComplete();
				for (int i = 0; i < spliceChains.length; i++) { 						
					Arrays.sort(spliceChains[i], c);	// TODO remove, efficiency, just to make sure..
				}
				
					// find common sites
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
				for (int i = 0; i < tokens.size(); i++) {
					pos= (int[]) tokens.elementAt(i);
								
					SpliceSite[][] ass= new SpliceSite[pos.length][];
					for (int j = 0; j < ass.length; j++) {
						int len= pos[j]- posOld[j]- 1;
						if (len< 0)
							break;
						ass[j]= new SpliceSite[len];
						for (int k = 0; k < ass[j].length; k++) { 
							ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
						}
					}
					for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
						if (ass[j]!= null&& ass[j].length> 0) {
							ASVariation var= new ASVariation(t[x], t[y], ass[0], ass[1]);
							SpliceSite flank5;
							if (i== 0)
								flank5= null;
							else
								flank5= spliceChains[j][posOld[j]];
							SpliceSite flank3= spliceChains[j][pos[j]];
							var.setAnchors(flank5, flank3);
							result.add(var);
							break;
						}
					posOld= pos;
				}
				
					// last
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
						ASVariation var= new ASVariation(t[x], t[y], ass[0], ass[1]);
						SpliceSite flank5;
						if (tokens.size()== 0)
							flank5= null;
						else
							flank5= spliceChains[j][posOld[j]];
						SpliceSite flank3= null;
						var.setAnchors(flank5, flank3);
						result.add(var);
						break;
					}
				

			}
		}
		
		ASVariation[] vars= (ASVariation[]) gphase.tools.Arrays.toField(result);
		
		return vars;
	}
	
	/**
	 * gets non-redundant set of coding exons
	 * @param completelyCoding
	 * @return
	 */
	public Exon[] getCodingExons(boolean completelyCoding) {
		Exon[] ex= getExons();
		Vector<Exon> resEx= new Vector<Exon>();
		for (int i = 0; i < ex.length; i++) {
			if (completelyCoding&& ex[i].isCodingSomewhere5Prime()&& 
					ex[i].isCodingSomewhere3Prime())
				resEx.add(ex[i]);
			if ((!completelyCoding)&& ex[i].overlapsCDS())
				resEx.add(ex[i]);
		}
		return (Exon[]) gphase.tools.Arrays.toField(resEx);
	}
	
	public int getTranscriptNbFromSource(String source) {
		int nb= 0;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranscriptID().contains(source))
				++nb;
		}
		return nb;
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
							//if (ex1[k].get)
						}
					}
					if (m== ex2.length&& !constitutive)
						v.add(ex1);
				}
			}
		}
		return null;
	}
	
	public ASMultiVariation[] getASMultiVariations(HashMap refTrptIDs) {
			
			asComplexes= null;
			if (asComplexes == null) {
				if (transcripts== null|| transcripts.length< 2) 
					return null;			
				
				SpliceSite[][] spliceChains= new SpliceSite[transcripts.length][];
				for (int i = 0; i < transcripts.length; i++) 
					spliceChains[i]= transcripts[i].getSpliceChainComplete();
				
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
							
							if (refTrptIDs!= null&& refTrptIDs.size()> 0&& 
									refTrptIDs.get(transcripts[i].getTranscriptID())== null&& 
									refTrptIDs.get(transcripts[j].getTranscriptID())== null)
								continue;
							
							SpliceSite[][] ss2= new SpliceSite[2][];
							ss2[0]= cluster[i];
							ss2[1]= cluster[j];
							if (ss2[0].length< 1&& ss2[1].length< 1)	// skip when both are empty
								continue;
							
		
							Transcript[] tt= new Transcript[2];
							tt[0]= transcripts[i];
							tt[1]= transcripts[j];
							
							AbstractSite[] trimSites= new AbstractSite[2];
							//ss2= ASVariation.trim(ss2, tt, trimSites);	// trim left in first, last at right
							
							//Vector cc2= tokenizeASClusters(ss2, false, false);
							Vector cc2= tokenizeASEvents(ss2, tt, false, false);
							for (int k = 0; k < cc2.size(); k++) {			// maybe more than one AS variation for a pair
	//							SpliceSite[][] pwEvent= (SpliceSite[][]) cc2.elementAt(k);
	//							ASVariation as2= new ASVariation(transcripts[i], transcripts[j],
	//									pwEvent[0], pwEvent[1]);
	//							as2Events.add(as2);
								as2Events.add(cc2.elementAt(k));
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

	public static ASMultiVariation[] getASMultiVariations(Transcript[] t, HashMap refTrptIDs) {
			
			ASMultiVariation[] asComplexes= null;
			if (t== null|| t.length< 2) 
				return null;			
			
			SpliceSite[][] spliceChains= new SpliceSite[t.length][];
			for (int i = 0; i < t.length; i++) 
				spliceChains[i]= t[i].getSpliceChainComplete();
			
			Vector spliceClusters= tokenizeASClusters(spliceChains, false, false);	// get splice clusters across all sequences
																			//TODO warning! TSS/TES flanked events included !!
			
				// determine pw splice vars
			Vector asComp= new Vector(spliceClusters.size());
			int[] rightEdge= new int[t.length];
			for (int x = 0; x < spliceClusters.size(); x++) {	// iterate complexes
				
				SpliceSite[][] cluster= (SpliceSite[][]) spliceClusters.elementAt(x);
				
				Vector as2Events= new Vector();
				for (int i = 0; i < t.length; i++) {		// iterate pw combinations for a complex 
					for (int j = (i+1); j < t.length; j++) {
						
						if (refTrptIDs!= null&& refTrptIDs.size()> 0&& 
								refTrptIDs.get(t[i].getTranscriptID())== null&& 
								refTrptIDs.get(t[j].getTranscriptID())== null)
							continue;
						
						SpliceSite[][] ss2= new SpliceSite[2][];
						ss2[0]= cluster[i];
						ss2[1]= cluster[j];
						if (ss2[0].length< 1&& ss2[1].length< 1)	// skip when both are empty
							continue;
						
	
						Transcript[] tt= new Transcript[2];
						tt[0]= t[i];
						tt[1]= t[j];
						
						AbstractSite[] trimSites= new AbstractSite[2];
						//ss2= ASVariation.trim(ss2, tt, trimSites);	// trim left in first, last at right
						
						//Vector cc2= tokenizeASClusters(ss2, false, false);
						Vector cc2= tokenizeASEvents(ss2, tt, false, false);
						for (int k = 0; k < cc2.size(); k++) {			// maybe more than one AS variation for a pair
//							SpliceSite[][] pwEvent= (SpliceSite[][]) cc2.elementAt(k);
//							ASVariation as2= new ASVariation(transcripts[i], transcripts[j],
//									pwEvent[0], pwEvent[1]);
//							as2Events.add(as2);
							as2Events.add(cc2.elementAt(k));
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
			
			return asComplexes;
		}

	public ASMultiVariation[] getASMultiVariations2() {
		
		if (transcripts== null|| transcripts.length< 2) 
			return null;			
		
		SpliceSite[][] spliceChains= new SpliceSite[transcripts.length][];
		for (int i = 0; i < transcripts.length; i++) 
			spliceChains[i]= transcripts[i].getSpliceChain();
		
		Vector spliceClusters= tokenizeASClusters2(spliceChains, false, false);	// get splice clusters across all sequences
																		//TODO warning! TSS/TES flanked events included !!
		ASMultiVariation[] multiVars= new ASMultiVariation[spliceClusters.size()];
		for (int x = 0; x < spliceClusters.size(); x++) {	// iterate complexes
			//SpliceSite[][] parts= (SpliceSite[][]) spliceClusters.elementAt(x);
			HashMap trptHash= (HashMap) spliceClusters.elementAt(x);
			multiVars[x]= new ASMultiVariation(trptHash);
			System.currentTimeMillis();
		}
	
		return multiVars;
	}
	
	public int[] getASMultiVariationsDebug() {
		
		asComplexes= null;
		int[] result= new int[3];
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
			result[0]++;
			for (int x = 0; x < spliceClusters.size(); x++) {	// iterate complexes
				result[1]++;
				
				SpliceSite[][] cluster= (SpliceSite[][]) spliceClusters.elementAt(x);
				
				Vector as2Events= new Vector();
				for (int i = 0; i < transcripts.length; i++) {		// iterate pw combinations for a complex 
					for (int j = (i+1); j < transcripts.length; j++) {
						result[2]++;
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

		return result;
	}
	
	HashMap addTrptHash(HashMap trptHash, SpliceSite[] schain, Transcript trpt) {
		Object[] o= trptHash.keySet().toArray();
		int i;
		for (i = 0; i < o.length; i++) {
			SpliceSite[] keyChain= (SpliceSite[]) o[i];
			if (keyChain.length!= schain.length)
				continue;
			int j;
			for (j = 0; j < keyChain.length; j++) {
				if (keyChain[j].getPos()!= schain[j].getPos())
					break;
			}
			if (j>= keyChain.length) 
				break;
		}
		Vector v= new Vector();
		if (i< o.length) {
			v= (Vector) trptHash.remove(o[i]);
		}
		v.add(trpt);
		trptHash.put(schain, v);
		return trptHash;
	}
	
	Vector tokenizeASEvents(SpliceSite[][] spliceChains, Transcript[] tt, AbstractSite[] trimSites, boolean omitStart, boolean omitEnd) {
		
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
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) 
					break;
			
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (i== 0)
					flank5= addAbstractSite(trimSites[0]);
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= null;
				flank3= spliceChains[j][pos[j]];
				var.setAnchors(flank5, flank3);
				result.add(var);
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
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) 
					break;
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (tokens.size()== 0)
					flank5= addAbstractSite(trimSites[0]);
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= addAbstractSite(trimSites[1]);
				var.setAnchors(flank5, flank3);
				result.add(var);
			}
	
		}
		
		return result;
	}

	static Vector tokenizeASEvents(SpliceSite[][] spliceChains, Transcript[] tt, boolean omitStart, boolean omitEnd) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			AbstractSite s= spliceChains[0][i];
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
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) 
					break;
			
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (i== 0)
					flank5= null;
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= spliceChains[j][pos[j]];
				var.setAnchors(flank5, flank3);
				result.add(var);
				
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
			int j;
			for (j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) 
					break;
			if (j< ass.length) {
				ASVariation var= new ASVariation(tt[0], tt[1], ass[0], ass[1]);
				AbstractSite flank5= null;
				if (tokens.size()== 0)
					flank5= null;
				else
					flank5= spliceChains[j][posOld[j]];
				AbstractSite flank3= null;
				var.setAnchors(flank5, flank3);
				result.add(var);
				
			}

		}
		
		return result;
	}

	static Vector getASVariations(SpliceSite[][] spliceChains) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		for (int i = 0; i < spliceChains.length; i++) 	// just to make sure
			Arrays.sort(spliceChains[i], c);
		
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
						
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) {
					result.add(ass);
					break;
				}
			posOld= pos;
		}
		
			// last
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
		
		
		return result;
	}

	static Vector tokenizeASClusters(SpliceSite[][] spliceChains, boolean omitStart, boolean omitEnd) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionEqualSSTypeComparator();	// here! SpliceSite.PositionComparator
		for (int i = 0; i < spliceChains.length; i++) 	// just to make sure
			Arrays.sort(spliceChains[i], c);
		
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
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
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

	Vector tokenizeASClusters(AbstractSite[][] spliceChains, boolean omitStart, boolean omitEnd) {
		
			// find splice sites common in ALL transcripts
		SpliceSite.PositionComparator c= new SpliceSite.PositionComparator();	// here! SpliceSite.PositionComparator
		Vector tokens= new Vector();
		for (int i = 0; i < spliceChains[0].length; i++) {
			AbstractSite s= spliceChains[0][i];
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
			AbstractSite[][] ass= new AbstractSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
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
			AbstractSite[][] ass= new AbstractSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				ass[j]= new AbstractSite[pos[j]- posOld[j]- 1];
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

	Vector tokenizeASClusters2(SpliceSite[][] spliceChains, boolean omitStart, boolean omitEnd) {
		
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
			HashMap trptHash= new HashMap();
			pos= (int[]) tokens.elementAt(i);
			if (omitStart&& i== 0) {	// omit splice chains containing start
				posOld= pos;
				continue;
			}
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			Transcript[] trpt= new Transcript[pos.length];
			for (int j = 0; j < ass.length; j++) {
				int len= pos[j]- posOld[j]- 1;
				if (len< 0)
					break;
				ass[j]= new SpliceSite[len];
				for (int k = 0; k < ass[j].length; k++) { 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];
				}
				addTrptHash(trptHash, ass[j], transcripts[j]);
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs
				if (ass[j]!= null&& ass[j].length> 0) {
					//result.add(ass);
					result.add(trptHash);
					break;
				}
			posOld= pos;
		}
		
		if (!omitEnd) {
			HashMap trptHash= new HashMap();
			pos= new int[spliceChains.length];		// anchor after sequence end
			for (int i = 0; i < pos.length; i++) 
				pos[i]= spliceChains[i].length;
			SpliceSite[][] ass= new SpliceSite[pos.length][];
			for (int j = 0; j < ass.length; j++) {
				ass[j]= new SpliceSite[pos[j]- posOld[j]- 1];
				for (int k = 0; k < (pos[j]- posOld[j]- 1); k++) 
					ass[j][k]= spliceChains[j][posOld[j]+ k+ 1];				
				addTrptHash(trptHash, ass[j], transcripts[j]);
			}
			for (int j = 0; j < ass.length; j++) 	// at least one transcript needs to provide alternative SSs 
				if (ass[j].length> 0) {
					//result.add(ass);
					result.add(trptHash);
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
	 * called during clustering process
	 * @param newExon
	 * @return
	 */
	public boolean addExon(Exon newExon) {
		Comparator compi= new Exon.EqualComparator();	//Exon.IdentityComparator(); no, mess with dynamic ss init
		int x= -1;	// insertion point if empty
		if (exons!= null) {
			x= Arrays.binarySearch(exons, newExon, compi);
			if (x>= 0) {
				for (int i = 0; i < newExon.getTranscripts().length; i++) 
					newExon.getTranscripts()[i].replaceExon(newExon, exons[x]);
				return false;	// get other exon
			} else
				exons= (Exon[]) gphase.tools.Arrays.insert(
						exons, newExon, x);
		} else 
			exons= new Exon[] {newExon};
		
		addSpliceSite(newExon.getAcceptor());
		addSpliceSite(newExon.getDonor());
		return true;

	}
	
	/**
	 * 
	 * @param newExon
	 * @return exon already known or <code>null</code>
	 */
	public Exon getExon(Exon newExon) {
		if (exons== null|| exons.length< 1)
			return null;
		Comparator compi= new Exon.PositionSSComparator();	// before identity comparator
		int x= Arrays.binarySearch(exons, newExon, compi);
		if (x>= 0) 
			return exons[x];
		else
			return null;	// not found	
	}
	
	/**
	 * 
	 * @param newExon
	 * @return exon already known or <code>null</code>
	 */
	public boolean removeExon(Exon newExon) {
		if (exons== null|| exons.length< 1)
			return false;
		for (int i = 0; i < exons.length; i++) {
			if (exons[i]== newExon) {
				exons= (Exon[]) gphase.tools.Arrays.remove(exons, i);	// get other exon
				return true;
			}
		}

		return false;	// not found	

	}
	
	public void markAlternativeSpliceSites() {
		SpliceSite[] sites= getSpliceSites();
		Comparator compi= new SpliceSite.PositionTypeComparator();
		for (int i = 0; i < sites.length; i++) {
			if (sites[i].getTranscripts().length== getTranscripts().length) {
				sites[i].setModality(SpliceSite.CONSTITUTIVE);
				continue;
			}
			for (int j = 0; j < getTranscripts().length; j++) {
				SpliceSite[] ss= getTranscripts()[j].getSpliceSitesAll();	// not schain
				int p= Arrays.binarySearch(ss, sites[i], compi);
				if (p< 0) {
					p= -(p+ 1);
					if (p!= 0&& p!= ss.length) {
						sites[i].setModality(SpliceSite.ALTERNATIVE);
						break;
					}
				}
			}
		}
	}
	
	/**
	 * @param transcripts
	 */
	public boolean addTranscript(Transcript newTranscript) {

			// search transcript
//		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) 
//			if (transcripts[i].getStableID().equalsIgnoreCase(newTranscript.getStableID()))
//				return false;

		newTranscript.gene= this;
		updateBoundaries(newTranscript);
		
		if(transcripts== null) 
			transcripts= new Transcript[] {newTranscript};
		else {
			Transcript[] nTranscripts= new Transcript[transcripts.length+ 1];
			for (int i= 0; i < transcripts.length; i++) 
				nTranscripts[i]= transcripts[i];
			nTranscripts[nTranscripts.length- 1]= newTranscript;
			newTranscript.setGene(this);
			transcripts= nTranscripts;
		}
		
		for (int i = 0;newTranscript.getExons()!= null && 
					i < newTranscript.getExons().length; i++) {
			addExon(newTranscript.getExons()[i]);
		}
		//sites= null;	// are to be re-inited
		
		return true;
	}
	
	public Transcript getNameTranscript() {
		int minStart= Integer.MAX_VALUE;
		for (int i = 0; i < transcripts.length; i++) {
			if (transcripts[i].get5PrimeEdge()< minStart)
				minStart= transcripts[i].get5PrimeEdge();
		}
		Transcript t= null;
		for (int i = 0; i < transcripts.length; i++) {
			if ((transcripts[i].get5PrimeEdge()== minStart) &&(t== null ||
				(transcripts[i].getTranscriptID().compareTo(t.getTranscriptID())< 0)))
					t= transcripts[i];
		}
		return t;
	}
	
	public void updateBoundaries(DirectedRegion reg) {
		// update boundaries
		if (strand== 0)
			setStrand(reg.getStrand());
		if (reg.getStart()!= 0&& (start== 0|| Math.abs(reg.getStart())< Math.abs(getStart())))
			setStart(reg.getStart());
		if (reg.getEnd()!= 0&& (end== 0|| Math.abs(reg.getEnd())> Math.abs(getEnd())))
			setEnd(reg.getEnd());
		if (chromosome== null)
			setChromosome(reg.getChromosome());
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
	
	public boolean hasVariation(DirectedRegion intron, 
			String varCode, int filterCode) {
		SpliceSite donor= new SpliceSite(this, intron.get5PrimeEdge()- 1, true);
		SpliceSite acceptor= new SpliceSite(this, intron.get3PrimeEdge()+ 1, false);
		ASVariation[] vars= getASVariations(filterCode);
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			if (!vars[i].toString().equals(varCode))
				continue;
			Comparator compi= new SpliceSite.PositionComparator();
			SpliceSite[] sc1= vars[i].getSpliceChain1();
			if (sc1!= null) {
				int p1= Arrays.binarySearch(sc1, donor, compi);
				int p2= Arrays.binarySearch(sc1, acceptor, compi);
				if (p1>= 0&& p2>= 0&& p2==p1+1)
					return true;
			}
			SpliceSite[] sc2= vars[i].getSpliceChain2();
			if (sc2!= null) {
				int p1= Arrays.binarySearch(sc2, donor, compi);
				int p2= Arrays.binarySearch(sc2, acceptor, compi);
				if (p1>= 0&& p2>= 0&& p2==p1+1)
					return true;
			}
		}
		return false;
		
	}

	/**
	 * @return Returns the confidence.
	 */
	public String getConfidence() {
		return confidence;
	}
	
	public Exon[] getExons(int region) {
		Vector v= new Vector();
		Comparator compi= new DirectedRegion.PositionComparator();
		for (int i = 0; i < transcripts.length; i++) 
			v= (Vector) gphase.tools.Arrays.addUnique(v, transcripts[i].getExons(), compi);
			
		DirectedRegion reg= getRegion(region);
		if (reg== null)
			return null;
		for (int i = 0; i < v.size(); i++) 
			if (!reg.contains((Exon) v.elementAt(i))) {
				v.remove(i--);
				break;
			}
		
		return  (Exon[]) gphase.tools.Arrays.toField(v);
	}
	
	public DirectedRegion[] getIntrons(int region) {
		Intron[] introns= getIntrons();
		for (int i = 0; i < introns.length; i++) {
			if (DirectedRegion.checkOverlap(introns[i].getTranscript().))
		}
		Vector v= new Vector();
		Comparator compi= new DirectedRegion.PositionComparator();
		for (int i = 0; i < transcripts.length; i++) 
			v= (Vector) gphase.tools.Arrays.addUnique(v, transcripts[i].getIntrons(), compi);
			
		DirectedRegion reg= getRegion(region);
		if (reg== null)
			return null;
		for (int i = 0; i < v.size(); i++) 
			if (!reg.contains((DirectedRegion) v.elementAt(i))) {
				v.remove(i--);
				break;
			}
		
		return  (DirectedRegion[]) gphase.tools.Arrays.toField(v);
	}
	
	
	/**
	 * just returning global regions, eg real/max/transcript utr/cds
	 * intronic/exonic arrays delegated to submethods..
	 * @param regionID
	 * @return
	 */ 
	public DirectedRegion getRegion(int regionID) {
		
		Object o= null;
		if (regionID== REGION_REAL_5UTR)
			o= getReal5UTR();
		else if (regionID== REGION_REAL_CDS)
			o= getRealCDS();
		else if (regionID== REGION_REAL_3UTR)
			o= getReal3UTR();
		else if (regionID== REGION_MAX_CDS)
			o= getMaxCDS();
		else if (regionID== REGION_COMPLETE_GENE)
			o= this;	// new DirectedRegion(getStart(), getEnd(), getStrand());
		
		return (DirectedRegion) o;
	}
	

	/**
	 * @return Returns the type.
	 */
	public int getType() {
		return type;
	}

	String chromosome = null;
	public static final String[] REGION_ID= 
	{"REGION_COMPLETE_GENE", "REGION_REAL_5UTR", "REGION_REAL_CDS", "REGION_REAL_3UTR",
		"REGION_MAX_5UTR", "REGION_MAX_CDS", "REGION_MAX_3UTR", "REGION_TRANSCRIPT_5UTR",
		"REGION_TRANSCRIPT_CDS", "REGION_TRANSCRIPT_3UTR"};
	public static final int REGION_COMPLETE_GENE= 0;
	public static final int REGION_REAL_5UTR= 1;
	public static final int REGION_REAL_CDS= 2;
	public static final int REGION_REAL_3UTR= 3;
	public static final int REGION_MAX_5UTR= 4;
	public static final int REGION_MAX_CDS= 5;
	public static final int REGION_MAX_3UTR= 6;
	public static final int REGION_TRANSCRIPT_5UTR= 7;
	public static final int REGION_TRANSCRIPT_CDS= 8;
	public static final int REGION_TRANSCRIPT_3UTR= 9;

	/**
	 * @return
	 */
	public String getChromosome() {
		return chromosome;
	}

	/**
	 * @param string
	 */
	public void setChromosome(String string) {
		String stringU= string.toUpperCase();
		if (stringU.startsWith("SCAFFOLD")|| stringU.startsWith("REFTIG")
				|| stringU.startsWith("CONTIG")|| stringU.startsWith("CHR")
				|| stringU.startsWith("GROUP"))
			chromosome= string;
		else {
			if (stringU.equals("MT"))	// correct Ensembl to GRIB jargon
				chromosome= "M";
			else {						// add cher
				if (string.startsWith("0"))
					string= string.substring(1, string.length());
				chromosome= "chr"+ string;
			}
		}
	}

	public SpliceSite[] getSpliceSites() {
		return spliceSites;	// evtl lazy extractions from exons
//		spliceSites= new SpliceSite[0];
//		Comparator compi= new SpliceSite.PositionTypeComparator();
//		for (int i = 0; i < getTranscripts().length; i++) {
//			SpliceSite[] ss= getTranscripts()[i].getSpliceChain();
//			for (int j = 0; j < ss.length; j++) {
//				int p= Arrays.binarySearch(spliceSites, ss[j], compi);
//				if (p< 0)
//					spliceSites= (SpliceSite[]) gphase.tools.Arrays.insert(spliceSites, ss[j], p);
//				else
//					spliceSites[p].addTranscript(getTranscripts()[i]);
//			}
//		}
//		return spliceSites;
	}

	public HashMap getSeparationIndex() {
		ASMultiVariation[] mvars= getASMultiVariations2();
		HashMap sepHash= new HashMap();
		for (int i = 0; mvars!= null&& i < mvars.length; i++) {
			for (int j = (i+1); j < mvars.length; j++) {
					// only for not overlapping bubbles, cannot overlap for the moment..
				//if (mvars[i].getRegion().overlaps(mvars[j].getRegion()))
				//	continue;
				
				Tuple[] interSep= Tuple.intersect(mvars[i].getSeparationTuples(), mvars[j].getSeparationTuples());
				
				if (interSep== null|| interSep.length< 1)
					continue;
				Integer index= new Integer(interSep.length);
				Vector v= (Vector) sepHash.remove(index);
				if (v== null)
					v= new Vector();
				v.add(new ASMultiVariation[] {mvars[i], mvars[j]});
				sepHash.put(index, v);
			}
		}
		return sepHash;
	}

	private static int uniqueID= 0;
	
	public static String getUniqueID() {
		
		String nb= Integer.toString(uniqueID++);
		for (int i = 0; i < 13- nb.length(); i++) 
			nb= "0"+ nb;
		
		return nb;
	}

	
}

