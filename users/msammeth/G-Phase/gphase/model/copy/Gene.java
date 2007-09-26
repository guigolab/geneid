package gphase.model.copy;

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

	
	SpliceSite[] spliceSites= null; // ordered according to position
	Exon[] exons= null;
	boolean construct= true;	// used for hashing during construction
	
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
	public SpliceSite addSpliceSite(SpliceSite ss) {
		
//		if (1== 1)
//			return false;
//		
//		if (getSite(ss)!= null)
//			return false;
		
			// insert new splice site
		if(spliceSites== null) { 
			spliceSites= new SpliceSite[] {ss};
			return ss;
		} else {
			int p= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());	// must positionType, since there can be 2 at the same pos
			boolean replaceNew= true;
			if (p< 0) {	// check another time for "variant" site
				int q;
				if (ss.getType()== SpliceSite.TYPE_TSS) {
					ss.type= SpliceSite.TYPE_ACCEPTOR;
					q= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
					if (q>= 0) 	// found
						p= q;	// replace new ss
					else
						ss.type= SpliceSite.TYPE_TSS;
				}
				if (ss.getType()== SpliceSite.TYPE_TES) {
					ss.type= SpliceSite.TYPE_DONOR;
					q= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
					if (q>= 0) 	// found
						p= q;
					else
						ss.type= SpliceSite.TYPE_TES;
				}
				if (ss.getType()== SpliceSite.TYPE_ACCEPTOR) {
					ss.type= SpliceSite.TYPE_TSS;
					q= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
					if (q>= 0) {	// found
						p= q;
						replaceNew= false;
					} else
						ss.type= SpliceSite.TYPE_ACCEPTOR;
				}
				if (ss.getType()== SpliceSite.TYPE_DONOR) {
					ss.type= SpliceSite.TYPE_TES;
					q= Arrays.binarySearch(spliceSites, ss, SpliceSite.getDefaultPositionTypeComparator());
					if (q>= 0) {	// found
						p= q;
						replaceNew= false;
					} else
						ss.type= SpliceSite.TYPE_DONOR;
				}
			}
			
				// merge
			if (p< 0)
				spliceSites= (SpliceSite[]) gphase.tools.Arrays.insert(spliceSites, ss, p);
			else {
				if (isConstruct()) {
					if (replaceNew) {
						Iterator iter= ss.getExonHash().values().iterator();
						while (iter.hasNext())
							((Exon) iter.next()).replaceSite(ss, spliceSites[p]);
						
						iter= ss.getTransHash().values().iterator();
						while (iter.hasNext())
							spliceSites[p].addTranscript(((Transcript) iter.next()));
					} else {
						Iterator iter= spliceSites[p].getExonHash().values().iterator();
						while (iter.hasNext())
							((Exon) iter.next()).replaceSite(spliceSites[p], ss);
						
						iter= spliceSites[p].getTransHash().values().iterator();
						while (iter.hasNext())
							ss.addTranscript(((Transcript) iter.next()));
						spliceSites[p]= ss;
					}
				} else {
					if (replaceNew) {
						for (int i = 0; i < ss.getExons().length; i++) 
							ss.getExons()[i].replaceSite(ss, spliceSites[p]);
						for (int i = 0; i < ss.getTranscripts().length; i++) 
							spliceSites[p].addTranscript(ss.getTranscripts()[i]);
					} else {
						for (int i = 0; i < spliceSites[p].getExons().length; i++) 
							spliceSites[p].getExons()[i].replaceSite(spliceSites[p], ss);
						for (int i = 0; i < spliceSites[p].getTranscripts().length; i++) 
							ss.addTranscript(spliceSites[p].getTranscripts()[i]);
						spliceSites[p]= ss;
					}
				}
				if (replaceNew)
					return spliceSites[p];
				else
					return ss;
			}
			return ss;
		}
		
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
		
	String geneID= null;
		
	
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
	
	public void repairAlignmentErrors() {
		SpliceSite[] ss= getSpliceSites();
		for (int i = 0; i < ss.length; i++) {
			
		}
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
		SpliceSite flank5= addSpliceSite(newExon.getAcceptor());
		SpliceSite flank3= addSpliceSite(newExon.getDonor());

		if (exons!= null) {
			int x= Arrays.binarySearch(exons, newExon, Exon.getDefaultPositionComparator());	//Exon.IdentityComparator(); no, mess with dynamic ss init
			if (x>= 0) {
				for (int i = 0; i < newExon.getTranscripts().length; i++) 
					newExon.getTranscripts()[i].replaceExon(newExon, exons[x]);
				return false;	// get other exon
			} else
				exons= (Exon[]) gphase.tools.Arrays.insert(
						exons, newExon, x);
		} else 
			exons= new Exon[] {newExon};
		
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

	private static int uniqueID= 0;
	
	public static String getUniqueID() {
		
		String nb= Integer.toString(uniqueID++);
		for (int i = 0; i < 13- nb.length(); i++) 
			nb= "0"+ nb;
		
		return nb;
	}

	public boolean isConstruct() {
		return construct;
	}

	public void setConstruct(boolean construct) {
		this.construct = construct;
	}

	
}

