/*
 * Created on Mar 31, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.Constants;
import gphase.tools.Arrays;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * 
 * 
 * @author micha
 */
public class Species implements Serializable {
	
	public static final String ENS_ID= "ENS";
	public static final String ENS_IMCB_ID= "SIN";
	public static final String GENOSCOPE_ID= "GS";
	public static final String FLYBASE_ID= "C";
	
	public static final String[] SP_NAMES_ENS_PFX= new String[] {
			ENS_ID+			"", 
			ENS_ID+			"PTR",
			ENS_ID+			"MUS", 
			ENS_ID+			"RNO",
			ENS_ID+			"CAF",
			ENS_ID+			"BTA",
			ENS_ID+			"MOD",
			
			ENS_ID+			"GAL",
			ENS_ID+			"XET",
			ENS_ID+			"DAR",
			ENS_IMCB_ID+	"FRU",	
			GENOSCOPE_ID+ 	"TEN", // GSTEN0039005408, HOXA1, HOXA10, HOXBb1, ...
			
			FLYBASE_ID, // CR for genes on forward strand, CG for genes on reverse strand
			ENS_ID+			"ANG",
			ENS_ID+			"APM",
			
			"",		// W09B6.4, Y59E9AR.9, ZK994.1, ZK994.3, ZK994.4, ... (no rule)
			""		// RDN58-1, snRNA, 15S-RNA, YAL016W, ... (no plan)
	};
	
	public static final boolean isValidEnsemblPrefix(String pfx) {
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) 
			if (pfx.equalsIgnoreCase(SP_NAMES_ENS_PFX[i]))
				return true;
		return false;
	}
	
	public static final String[] SP_NAMES_BINOMIAL= new String[] {
			"homo_sapiens", 	// mammals
			"pan_troglodytes",
			"mus_musculus", 
			"rattus_norvegicus",
			"canis_familiaris",
			"bos_taurus",
			"monodelphis_domestica",
			
			"gallus_gallus",	// other chordates
			"xenopus_tropicalis",
			"danio_rerio",
			"fugu_rubripes",
			"tetraodon_nigroviridis",
			
			"drosophila_melanogaster",	// insects
			"anopheles_gambiae",
			"apis_mellifera",
			
			"caenorhabditis_elegans",
			"saccharomyces_cerevisiae"
	};
	
	// public static final String[] SP_NAMES_ABBREV= new String[] {"hsapiens", "mmusculus", "rnorvegicus"};
	public static final String[] SP_NAMES_COMMON= new String[] {
			"human",	// mammals 
			"chimp",
			"mouse", 
			"rat",
			"dog",
			"cow",
			"opossum",

			"chicken",	// other chordates
			"frog",
			"zebrafish",
			"fugu",
			"tetraodon",
			
			"fruitfly",
			"mosquito",
			"bee",
			"worm",
			"yeast"
	};

	public static Species[] createByBinomialName(String[] binomialNames) {
		
		if (binomialNames== null)
			return null;
		
		Species[] specs= new Species[binomialNames.length];
		for (int i = 0; i < specs.length; i++) 
			specs[i]= new Species(binomialNames[i]);
		
		return specs;
	}
	
	EncodeRegion[] encodeRegions= null;
	
	int spNumber= -1;

	
	public static HashMap basePairMapping= new HashMap(5);
	static {
		basePairMapping.put(new Character('G'), new Character('C'));
		basePairMapping.put(new Character('C'), new Character('G'));
		basePairMapping.put(new Character('A'), new Character('T'));
		basePairMapping.put(new Character('T'), new Character('A'));
		basePairMapping.put(new Character('N'), new Character('N'));
	}

	public boolean addEncodeRegion(EncodeRegion region) {
		
		if (region== null)
			return false;
		
		if (encodeRegions== null) {
			encodeRegions= new EncodeRegion[]{region};
			return true;
		}
		
		for (int i = 0; i < encodeRegions.length; i++) 
			if (encodeRegions[i].equals(region))
				return false;
		
		EncodeRegion[] newRegions= new EncodeRegion[encodeRegions.length+ 1];
		for (int i = 0; i < encodeRegions.length; i++) 
			newRegions[i]= encodeRegions[i];
		newRegions[encodeRegions.length]= region;
		encodeRegions= newRegions;
		return true;
	}
	
	/**
	 * @return Returns the ensemblPrefix.
	 */
	public String getEnsemblPrefix() {
		return SP_NAMES_ENS_PFX[spNumber];
	}

	public String getBinomialName() {
		return SP_NAMES_BINOMIAL[spNumber];
	}
	
	public String getCommonName() {
		return SP_NAMES_COMMON[spNumber];
	}
	
	public String getNameFirst() {
		return SP_NAMES_BINOMIAL[spNumber].substring(0, SP_NAMES_BINOMIAL[spNumber].indexOf('_'));
	}	
	
	static final String abbreviateBinomialName(String binName) {
		return binName.charAt(0)+ binName.substring(binName.indexOf("_")+ 1);		
	}
	
	public String getNameAbbrev() {
		return abbreviateBinomialName(SP_NAMES_BINOMIAL[spNumber]);
	}
	
	
	int buildVersion= -1; 	// assuming a common build src (like NCBI, RGSC, ...)
	HashMap geneHash= null;

	public static final String getCommonNameForPrefix(String ensemblPFX) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) 
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPFX))
				return SP_NAMES_COMMON[i];
		
		return null;
	}
	
	public static final String getBinomialForCommonName(String commonName) {
		
		for (int i = 0; i < SP_NAMES_COMMON.length; i++) 
			if (SP_NAMES_COMMON[i].equalsIgnoreCase(commonName))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}
	
	public static final String getBinomialForFirstName(String firstName) {
		
		firstName= firstName.toLowerCase();
		for (int i = 0; i < SP_NAMES_BINOMIAL.length; i++) 
			if (SP_NAMES_BINOMIAL[i].substring(0,SP_NAMES_BINOMIAL[i].indexOf('_')).equals(firstName))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}	
	
	public static final String getBinomialForEnsemblPfx(String ensemblPfx) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) 
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPfx))
				return SP_NAMES_BINOMIAL[i];
		
		return null;
	}		
	
	public static final String getBinomialForSomeName(String someName) {
		
		String bin= getBinomialForCommonName(someName);
		if (bin!= null)
			return bin;
		
		bin= getBinomialForFirstName(someName);
		if (bin!= null)
			return bin;
		
		bin= getBinomialForEnsemblPfx(someName);
		if (bin!= null)
			return bin;

		return null;
	}	
	
	public static final String getAbbrevNameForBionomial(String ensemblPFX) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) {
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPFX))
				return abbreviateBinomialName(SP_NAMES_BINOMIAL[i]);
		}
		
		return null;
	}
	
	public boolean remove(Gene g, boolean silent) {
		
		if (!silent)
			System.out.print("Removing gene "+ g.getStableID()+":");
		Object o= geneHash.remove(g.getStableID());		
		if (o!= null) {
			
			g= (Gene) o;
			GeneHomology[]  hgies= g.getHomologies();
			for (int i = 0; hgies!= null&& i < hgies.length; i++) {
				Gene gg= hgies[i].getOtherGene(g);
				gg.removeHomology(g);
			}
			
			if (!silent)
				System.out.println(" ok.");
			return true;
		}
		if (!silent)
			System.out.println(" failed!");
		return false;
	}

	public static int getSpeciesNumber(String someName) {
		
		return getSpeciesNumberForBinomialName(getBinomialForSomeName(someName));
	}
	
	public static int getSpeciesNumberForBinomialName(String binName) {
		int i;
		for (i = 0; i < SP_NAMES_BINOMIAL.length; i++) 
			if (SP_NAMES_BINOMIAL[i]== binName)
				break;
		if (i>= SP_NAMES_BINOMIAL.length) {
			System.err.println("Unknown name "+ binName);
			return -1;
		}
		
		return i;
	}
	
	public Species(String someName) {

		int i= getSpeciesNumber(someName); 
			// else
		spNumber= i;
	}
	
	/**
	 * Filters out genes that have at least <code>minSpec</code> orthologs in different species.
	 * 
	 * @deprecated no longer used, filtering already done by DB query
	 * @param minSpec how many species
	 * @param unique  allow only one homology per species
	 * @return
	 */
	public Gene[] getASEvolutiveGenes(int minSpec, boolean unique) {
		if (geneHash== null)
			return null;
	
		Iterator iter= geneHash.values().iterator();
		Gene g;
		Vector o= new Vector();
		while (iter.hasNext()) {			// iterate over all genes of species
			g= (Gene) iter.next();
			if (g.getHomologies().length< (minSpec- 1))// filter for min species conservation
				continue;
			int ref= g.getTranscripts().length;
			GeneHomology[] homol= g.getHomologies();
			Gene[] h= new Gene[homol.length];
			for (int i = 0; i < h.length; i++) 
				h[i]= homol[i].getOtherGene(g);
			for (int i = 0; i < h.length; i++) 	// get genes w different nb of transcript
				if (h[i].getTranscripts().length!= ref) {
					o.add(g);
					break;						// evolutive condition -> changing transcripts in time
				}
			if (unique)							// unique (ortholog check)
				for (int i = 0; i < h.length; i++) 
					for (int j = (i+1); j < h.length; j++) 
						if (h[i].getSpecies()== h[j].getSpecies())
							o.remove(g);
		}

		Gene[] result= new Gene[o.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) o.elementAt(i);
		return result;	
	}
	
	public Gene[] getASGenes() {
		
		if (geneHash== null)
			return null;
	
		Iterator iter= geneHash.values().iterator();
		Gene g;
		Vector o= new Vector();
		while (iter.hasNext()) {
			g= (Gene) iter.next();
			if (g.getTranscripts().length> 1)
				o.add(g);
		}

		Gene[] result= new Gene[o.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) o.elementAt(i);
		return result;	
	}
	
	
	public static char[] invertSequence(char[] seq) {
		
		char[] invertedSeq= new char[seq.length];
		for (int i = 0; i < invertedSeq.length; i++) {
			
			char c= seq[i];
			boolean lowerCase= false;
			if (Character.isLowerCase(c)) {
				c= Character.toUpperCase(c);
				lowerCase= true;
			}
			invertedSeq[i]= ((Character) basePairMapping.get(new Character(c))).charValue();
			if (lowerCase)
				invertedSeq[i]= Character.toLowerCase(invertedSeq[i]);
		}
		
		return invertedSeq;
	}
	
	/**
	 * @deprecated no longer used, see Graph.getSequenceDirectory
	 */
	public String getSequenceDirectory() {
	
		String seqDirName= Character.toUpperCase(getBinomialName().charAt(0))+ "."
				+ getBinomialName().substring(getBinomialName().indexOf('_')+ 1);
		File speciesGenome= new File(Constants.DATA_DIR+ File.separator
				+ Constants.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		Pattern patty= Pattern.compile("(\\D+)(\\d+)(\\D*)");
		Matcher matty;
		int highestVersion= -1;
		String goldenPath= null;
		String[] list= speciesGenome.list();
		for (int i = 0; i < list.length; i++) {
			matty= patty.matcher(list[i]);
			if (!matty.matches())
				continue;
			if (!matty.group(1).equals("golden_path_"))
				continue;
			int ver= Integer.parseInt(matty.group(2));
			if (ver> highestVersion) {
				highestVersion= ver;
				goldenPath= list[i];
			}
		}
		
		return speciesGenome.getAbsolutePath()+ File.separator
				+ goldenPath+ File.separator+ Constants.CHROMOSOME_DIR;
	}
	
	public Species(String newBinomialName, int nbGenes) {
		
		this(newBinomialName);
		setEstimatedGeneNb(nbGenes);
	}
	
	/**
	 * Initilizes the <code>HashMap</code> with the specified capacity and the
	 * specified loading factor for the expected number of genes. 
	 * 
	 * @param nbGenes
	 * @return <code>true</code> when the <code>HashMap</code> was successfully
	 * initialized, <code>false</code> when it had already been filled with
	 * values beforeahead.
	 */
	public boolean setEstimatedGeneNb(int nbGenes, int loadFactor) {
		
		if (geneHash!= null&& geneHash.values().size()> 0)
			return false;
			
			// else
		geneHash= new HashMap(nbGenes, loadFactor);
		return true;
	}
	
	/**
	 * Initilizes the <code>HashMap</code> with the capacity and a default
	 * loading factor of <code>nbGenes/100</code> for the expected number 
	 * of genes.
	 * 
	 * @param nbGenes
	 * @return <code>true</code> when the <code>HashMap</code> was successfully
	 * initialized, <code>false</code> when it had already been filled with
	 * values beforeahead.
	 */
	public boolean setEstimatedGeneNb(int nbGenes) {
	
		return setEstimatedGeneNb(nbGenes, nbGenes/100);
	}

	
	public Gene getGene(String stableID) {

		return (Gene) geneHash.get(stableID);		// get gene
	}	
	public Gene[] getGenes() {
	
		if (geneHash== null)
			return null;
	
		Object[] o= geneHash.values().toArray();
		return Gene.toGeneArray(o);	
	}
	
	
	public Iterator getGeneIterator() {
		return geneHash.values().iterator();
	}


	public boolean addGene(Gene newGene) {
		
			// get speciesHash
		if (geneHash== null) 
			geneHash= new HashMap();
		
			// 
		if (geneHash.get(newGene.getStableID())!= null)
			return false;	// gene exists
		geneHash.put(newGene.getStableID(), newGene);	// else
		return true;
	}

	public int getTranscriptCount() {
		if (geneHash== null) 
			return 0;
		Gene[] ge= getGenes();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			cnt+= ge[i].getTranscriptCount();
		return cnt;
	}
	public Transcript[] getTranscripts(String speciesID) {
		
		Gene[] genes= getGenes();
		Vector transV= new Vector(genes.length);
		int nctr= 0;
		for (int i = 0; i < genes.length; i++) {			// for all genes
			Transcript[] trans= genes[i].getTranscripts();
			if (trans== null) {
				nctr++;
				continue;
			}
			for (int j = 0; j < trans.length; j++) 			// add all transcripts
				transV.add(trans[j]);
		}
		if (nctr> 0)
			System.err.println(nctr+ " genes wo transcipts!");
		
		Transcript[] result= new Transcript[transV.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Transcript) transV.elementAt(i);
		return result;
	}

	static final long serialVersionUID = 5861506554499528191L;

	/**
	 * More effective than <code>getGenes().length</code> since it does not 
	 * convert <code>Object[]</code> to <code>Gene[]</code>.
	 * @return 
	 */
	public int getGeneNb() {
	
		if (geneHash== null|| geneHash.values()== null)
			return 0;
		
		return geneHash.values().size();
	}
	
	public SpliceSite[] getSpliceSites(int regionType, int ssType) {
		int perc= 0;
		Gene[] ge= getGenes();
		Vector v= new Vector();
		for (int i = 0; i < ge.length; i++) {
			SpliceSite[] ssites= ge[i].getSpliceSites(ssType, regionType);
			for (int j = 0; ssites!= null&& j < ssites.length; j++) 
				v.add(ssites[j]);
		}
		return (SpliceSite[]) Arrays.toField(v);
	}
	/**
	 * @return Returns the buildVersion.
	 */
	public int getBuildVersion() {
		return buildVersion;
	}
	public Exon[] getExons() {

		Vector v= new Vector();
		Iterator iter= geneHash.values().iterator();
		while(iter.hasNext()) {
			Gene ge= (Gene) iter.next();
			Exon[] ex= ge.getExons();
			for (int i = 0; i < ex.length; i++) {
				v.add(ex[i]);
			}
		}
		return (Exon[]) Arrays.toField(v);
	}

	/**
	 * @param buildVersion The buildVersion to set.
	 */
	public void setBuildVersion(int buildVersion) {
		this.buildVersion = buildVersion;
	}

	public String toString() {
		return getCommonName();
	}
	public EncodeRegion[] getEncodeRegions() {
		return encodeRegions;
	}
	
	public int countNonProteinCodingLoci() {
		Gene[] ge= getGenes();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			if (!ge[i].isProteinCoding())
				++cnt;
		return cnt;
	}
	
	public void initTU() {
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) {
			ge[i].initTU();
		}
	}
	
	public void recluster() {
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[][] newTranscripts= ge[i].recluster();
			
			if (newTranscripts== null|| newTranscripts.length== 0) {
				remove(ge[i], true);
				continue;
			}
			
			if (newTranscripts.length< 2)	// nothing happened, still 1 cluster
				continue;
			
				// create new clusters
			for (int j = 0; j < newTranscripts.length; j++) {
				Gene nuGene= new Gene(this, Gene.getUniqueID());
				nuGene.setChromosome(ge[i].getChromosome());
				nuGene.setStrand(ge[i].getStrand());
				for (int k = 0; k < newTranscripts[j].length; k++) 	{
					Transcript t= new Transcript(nuGene, newTranscripts[j][k].getTranscriptID());
					t.setStrand(newTranscripts[j][k].getStrand());
					t.setStart(newTranscripts[j][k].getStart());
					t.setEnd(newTranscripts[j][k].getEnd());
					t.setHUGO(newTranscripts[j][k].getHUGO());
					//Transcript t= newTranscripts.removeAll
					nuGene.addTranscript(t);
					Exon[] ex= newTranscripts[j][k].getExons();
					for (int m = 0; m < ex.length; m++) {
						Exon e= new Exon(t, ex[m].getExonID(), ex[m].getStart(), ex[m].getEnd());
						e.setStartCDS(ex[m].getStartCDS());
						e.setEndCDS(ex[m].getEndCDS());
						e.setFrame(ex[m].getFrame());
						t.addExon(e);
					}
					Translation[] trans= newTranscripts[j][k].getTranslations();
					for (int n = 0; trans!= null&& n < trans.length; n++) {
						Translation tra= new Translation(t);
						tra.setStart(newTranscripts[j][k].getTranslations()[n].getStart());	// bugfix 21.6.
						tra.setEnd(newTranscripts[j][k].getTranslations()[n].getEnd());
						t.addTranslation(tra);
					}
				}
				geneHash.put(nuGene.getGeneID(), nuGene);
			}			
			remove(ge[i], true);
		}
	}

	public ASVariation[][] getASVariations(int filter) {
		
		ASVariation[][] asClasses= null;
		if (asClasses == null) {
			Gene[] ge= getGenes();
			int asVariations= 0;
			Comparator compi= new ASVariation.SpliceStringComparator();
			Vector asClassesV= new Vector();
			for (int i = 0; i < ge.length; i++) {	// all genes
				ASMultiVariation[] as= ge[i].getASMultiVariations();
				for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
					ASVariation[] asvars= null;
					if (filter== ASMultiVariation.FILTER_NONE)	// all
						asvars= as[j].getASVariationsAll();
					else if (filter== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
						asvars= as[j].getASVariationsHierarchicallyFiltered();
					else if (filter== ASMultiVariation.FILTER_CODING_REDUNDANT)
						asvars= as[j].getASVariationsClusteredCoding();
					else if (filter== ASMultiVariation.FILTER_STRUCTURALLY)
						asvars= as[j].getASVariationsStructurallyFiltered();
					asVariations+= asvars.length;
					for (int k = 0; k < asvars.length; k++) {	// get pw variations
						int m;
						for (m = 0; m < asClassesV.size(); m++)	// compare against already existing classes 
							if (compi.compare(asvars[k],
									((Vector) asClassesV.elementAt(m)).elementAt(0))== 0) {
								((Vector) asClassesV.elementAt(m)).add(asvars[k]);
								break;
							}
						if (m>= asClassesV.size()) {
							Vector newClass= new Vector();
							newClass.add(asvars[k]);
							asClassesV.add(newClass);
						}
					}
				}
			}
			
			try {
				asClasses= (ASVariation[][]) gphase.tools.Arrays.toField(asClassesV);
			} catch (ClassCastException e) {
				System.err.println("No AS Classes!");
			}
		}
	
		return asClasses;
	}
}
