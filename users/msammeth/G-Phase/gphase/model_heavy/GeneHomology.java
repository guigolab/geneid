/*
 * Created on Nov 8, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model_heavy;

import gphase.Constants;
import gphase.tools.Arrays;

import java.io.Serializable;
import java.util.Collection;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Vector;

import qalign.algo.CancelException;
import qalign.algo.CostTable;
import qalign.algo.msa.QAlign;

/**
 * 
 * 
 * @author msammeth
 */
public class GeneHomology implements Serializable {

	public static class StableIDComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			GeneHomology gh1= (GeneHomology) arg0;
			GeneHomology gh2= (GeneHomology) arg1;
			
			return gh1.getStableIDs().compareTo(gh2.getStableIDs());
		}
	}
	
	final static String[] METHODS= {	// field "type" in table "method_link" (ensembl.compara)
			"ENSEMBL_ORTHOLOGUES",
			"ENSEMBL_PARALOGUES",
			"ENSEMBL_HOMOLOGUES",	// new
			"FAMILY",
			"SYNTENY",
			"MLAGAN",
			"TRANSLATED_BLAT",
			"PHUSION_BLASTN",
			"PHUSION_BLASTN_TIGHT",
			"BLASTZ_GROUP",
			"BLASTZ_GROUP_TIGHT",
			"BLASTZ_NET",
			"BLASTZ_NET_TIGHT",
			"BLASTZ_NET_RECIP_NET"
	};
	final static String[] TYPES_old= {	// field "description" in table "homology" (ensembl.compara)
			"UBRH", 
			"MBRH", 			// multiple best reciprocal hit
			"RHS", 				// reziprocal hit syntheny
			"DWGA",				// derived from whole genome alignment 
			"YoungParalogues"	
	};
	final static String[] TYPES= {
		"ortholog_one2one",
		"apparent_ortholog_one2one",
		"ortholog_one2many",
		"between_species_paralog",
		"ortholog_many2many",
		"within_species_paralog"
	};
	final static String[] SUBTYPES= {
		"ENSEMBL_ORTHOLOGUES",
		"ENSEMBL_PARALOGUES",
		"SYNTHENY"
	};
	final static String[] SUBTYPES_old= {
			"DUP1.1",
			"DUP1.2",
			"DUP1.3",
			"DUP1.4",
			"DUP1.6",
			"DUP1.7",
			"SYN",
			"complex",
	};
	
	public final static GeneHomology[] toGeneHomologyArray(Collection c) {
		
		GeneHomology[] result= new GeneHomology[c.size()];
		Iterator iter= c.iterator();
		int ctr= 0;
		while(iter.hasNext())
			result[ctr++]= (GeneHomology) iter.next();
		
		return result;
	}
	
	Gene gene1= null;
	Gene gene2= null;
	SpliceSiteHomology[] homologSpliceSites= null;
	
		// DB: ensembl_compara, homology
	int method= Constants.NOINIT;		// homology.description & homology.subtype
	double dn;
	double ds;
	double n;
	double s;
	double lnl;
	double thresholdOnDS;

		// DB: ensembl_compara, homology
	int type= Constants.NOINIT;
	int subtype= Constants.NOINIT;	

	
		// DB: ensembl_compara, homology_member 
		// pairwise homologies, see docu !!
	int g1PercCov;
	int g1PercId;
	int g1PercPos;

	int g2PercCov;
	int g2PercId;
	int g2PercPos;
	
	SpliceSiteHomology[] homologSplices= null;

	/**
	 * 
	 */
	public GeneHomology(Gene newGene1, Gene newGene2) {
		gene1= newGene1;
		gene2= newGene2;
	}

	public boolean orthologSpliceSitesIdentified() {
		return (homologSplices!= null);
	}
	
	public void identifyOrthologSpliceSites() {

		if (gene1== null|| gene2== null)
			return;
		
		SpliceSite[] ss1= gene1.getSpliceSites();
		SpliceSite[] ss2= gene2.getSpliceSites();
		SpliceSiteHomology[][] ssRelations= new SpliceSiteHomology[ss1.length][ss2.length];
		
		for (int i = 0; i < ss1.length; i++) {
			for (int j = 0; j < ss2.length; j++) {
					
				String sstring1= Graph.readSequence(ss1[i]);	// get sequences
				String sstring2= Graph.readSequence(ss2[j]);
					
				QAlign qal= new QAlign(); 	// align
				try {
					qal.setSequences(new String[] {sstring1, sstring2});
					qal.setCostTable(CostTable.DNA);
					qal.setSimultaneousAlignment(false);
					qal.run();
				} catch (CancelException e) {
					; // :)
				}
				
				SpliceSiteHomology ssh= new SpliceSiteHomology(ss1[i], ss2[j]);
				ssh.setAlignment(qal.getProgressiveLayout());
				ssh.setCost(qal.getCost());
				ssRelations[i][j]= ssh; 
			}
		}
				
			// find (U?)BRH
		Vector ubrh= new Vector(Math.min(ss1.length, ss2.length));
		for (int i = 0; i < ssRelations.length; i++) {
			float m= Float.MAX_VALUE;
			int mPos= -1;
			for (int j = 0; j < ssRelations[i].length; j++)		// find max in line 
				if (ssRelations[i][j].getCost()< m) {		// only bigger allowed, ok?!?
					m= ssRelations[i][j].getCost();
					mPos= j;
				}
				
			float mc= Float.MAX_VALUE;
			int mcPos= -1;
			for (int j = 0; j < ssRelations.length; j++) {		// equals max in column?!
				if (ssRelations[j][mPos].getCost()< m) {		// only bigger allowed, ok?!?
					mc= ssRelations[j][mPos].getCost();
					mcPos= j;
				}
			}
			
			if (mcPos== i) 	// best reciprocal hit
				ubrh.add(ssRelations[i][mPos]);			
		}
		
		this.homologSpliceSites= (SpliceSiteHomology[]) Arrays.toField(ubrh);
	}
	
	public Gene getOtherGene(Gene query) {
		
		if (query== null)
			return null;
		
		if (query== gene1)
			return gene2;
		if (query== gene2)
			return gene1;
		
		return null;		
	}
	
	/**
	 * 
	 * @param descrpt method descriptor
	 * @return <code>false</code> if method is already set or the method
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setMethod(String descrpt) {
		
		if (descrpt== null|| method!= Constants.NOINIT)
			return false;
		
		int i= Constants.findIgnoreCase(descrpt, METHODS);

		if (i< 0) {
			System.err.println("Unknown method "+ descrpt);
			return false;
		}

		method= i;
		return true;
	}
	
	/**
	 * 
	 * @param descrpt type descriptor
	 * @return <code>false</code> if type is already set or the type
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setType(String descrpt) {
		
		if (descrpt== null|| type!= Constants.NOINIT)
			return false;
		
		int i= Constants.findIgnoreCase(descrpt, TYPES);

		if (i< 0) {
			System.err.println("Unknown type "+ descrpt);
			return false;
		}

		type= i;
		return true;
	}
	
	/**
	 * 
	 * @param descrpt subtype descriptor
	 * @return <code>false</code> if subtype is already set or the subtype
	 * descriptor is not valid, <code>true</code> otherwise.  
	 */
	public boolean setSubtype(String descrpt) {
		
		if (descrpt== null|| type!= Constants.NOINIT)
			return false;
		
		int i= Constants.findIgnoreCase(descrpt, SUBTYPES);

		if (i< 0) {
			System.err.println("Unknown subtype "+ descrpt);
			return false;
		}

		subtype= i;
		return true;
	}
	
	/**
	 * @return Returns the dn.
	 */
	public double getDn() {
		return dn;
	}
	/**
	 * @param dn The dn to set.
	 */
	public void setDn(double dn) {
		this.dn = dn;
	}

	/**
	 * @return Returns the ds.
	 */
	public double getDs() {
		return ds;
	}
	/**
	 * @param ds The ds to set.
	 */
	public void setDs(double ds) {
		this.ds = ds;
	}
		/**
	 * @return Returns the lnl.
	 */
	public double getLnl() {
		return lnl;
	}
	/**
	 * @param lnl The lnl to set.
	 */
	public void setLnl(double lnl) {
		this.lnl = lnl;
	}
	/**
	 * @return Returns the n.
	 */
	public double getN() {
		return n;
	}
	/**
	 * @param n The n to set.
	 */
	public void setN(double n) {
		this.n = n;
	}
	/**
	 * @return Returns the s.
	 */
	public double getS() {
		return s;
	}
	/**
	 * @param s The s to set.
	 */
	public void setS(double s) {
		this.s = s;
	}
	/**
	 * @return Returns the thresholdOnDS.
	 */
	public double getThresholdOnDS() {
		return thresholdOnDS;
	}
	/**
	 * @param thresholdOnDS The thresholdOnDS to set.
	 */
	public void setThresholdOnDS(double thresholdOnDS) {
		this.thresholdOnDS = thresholdOnDS;
	}
		/**
		 * @return Returns the type.
		 */
		public int getType() {
			return type;
		}
		/**
		 * @param type The type to set.
		 */
		public void setType(int type) {
			this.type = type;
		}
	/**
	 * @return Returns the gene1.
	 */
	public Gene getGene1() {
		return gene1;
	}
	/**
	 * @param gene1 The gene1 to set.
	 */
	public void setGene1(Gene gene1) {
		this.gene1 = gene1;
	}
	/**
	 * @return Returns the gene2.
	 */
	public Gene getGene2() {
		return gene2;
	}
	/**
	 * @param gene2 The gene2 to set.
	 */
	public void setGene2(Gene gene2) {
		this.gene2 = gene2;
	}
	
	/**
	 * @return concatenated stableID Strings of both members
	 */
	public String getStableIDs() {
		return (getGene1().getStableID()+ getGene2().getStableID());
	}
		/**
		 * @return Returns the g1PercCov.
		 */
		public int getG1PercCov() {
			return g1PercCov;
		}
		/**
		 * @param percCov The g1PercCov to set.
		 */
		public void setG1PercCov(int percCov) {
			g1PercCov = percCov;
		}
	/**
	 * @return Returns the g1PercId.
	 */
	public int getG1PercId() {
		return g1PercId;
	}
	/**
	 * @param percId The g1PercId to set.
	 */
	public void setG1PercId(int percId) {
		g1PercId = percId;
	}
	/**
	 * @return Returns the g1PercPos.
	 */
	public int getG1PercPos() {
		return g1PercPos;
	}
	/**
	 * @param percPos The g1PercPos to set.
	 */
	public void setG1PercPos(int percPos) {
		g1PercPos = percPos;
	}
	/**
	 * @return Returns the g2PercCov.
	 */
	public int getG2PercCov() {
		return g2PercCov;
	}
	/**
	 * @param percCov The g2PercCov to set.
	 */
	public void setG2PercCov(int percCov) {
		g2PercCov = percCov;
	}
	/**
	 * @return Returns the g2PercId.
	 */
	public int getG2PercId() {
		return g2PercId;
	}
	/**
	 * @param percId The g2PercId to set.
	 */
	public void setG2PercId(int percId) {
		g2PercId = percId;
	}
	/**
	 * @return Returns the g2PercPos.
	 */
	public int getG2PercPos() {
		return g2PercPos;
	}
	/**
	 * @param percPos The g2PercPos to set.
	 */
	public void setG2PercPos(int percPos) {
		g2PercPos = percPos;
	}
}
