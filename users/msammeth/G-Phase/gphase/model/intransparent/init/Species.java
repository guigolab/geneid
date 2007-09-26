/*
 * Created on Mar 31, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.intransparent.init;

import gphase.Constants;
import gphase.io.gtf.GTFChrReader;
import gphase.tools.Arrays;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.io.PrintStream;
import java.io.Serializable;
import java.lang.reflect.Method;
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
	
	String genomeVersion= null;
	public static final String ENS_ID= "ENS";
	public static final String ENS_IMCB_ID= "SIN";
	public static final String GENOSCOPE_ID= "GS";
	public static final String FLYBASE_ID= "C";
	
	public static final String[] ASSEMBLY_PFX= new String[] {
		"ENSEMBL", "HG", "hg", "NCBI", "golden_path_", "gp", "WASHUC", "JGI", "Zv", "FUGU"
	};
	// org=
	public static final String[][][] SP_GENOME_VERSIONS= new String[][][] {
		new String[][] {new String[] {"GP200405", "NCBI35", "HG17"}, new String[] {"GP200603","NCBI36", "HG18"}},
		new String[][] {new String[] {"GP200603"}},
		new String[][] {new String[] {"GP200603", "MM8", "NCBI36"}},
		new String[][] {new String[] {"GP200411", "RGSC43"}},
		new String[][] {new String[] {"GP200505", "CANFAM2"}},
		new String[][] {new String[] {"GP200503", "BTAU2"}, new String[] {"BTAU31"}},
		new String[][] {new String[] {}},	// opossum
		
		new String[][] {new String[] {"GP200605", "WASHUC2"}},
		new String[][] {new String[] {"GP200508", "JGI41"}},
		new String[][] {new String[] {"GP200603", "ZV6"}},
		new String[][] {new String[] {}},	// fugu
		new String[][] {new String[] {}},	// tetraodon
		
		new String[][] {new String[] {"BDGP43"}},
		new String[][] {new String[] {}},	// mosquito
		new String[][] {new String[] {"GP200501", "APIMEL2"}},
		new String[][] {new String[] {"GP200607", "WS160"}},
		
		new String[][] {new String[] {}},	// yeast
		new String[][] {new String[] {}},	// seasquirt
		new String[][] {new String[] {}},	// seqsquirt2
		new String[][] {new String[] {}},	// platypus
		new String[][] {new String[] {}}	// cress
		
	};

	public static final String[][][] SP_ANNOTATION_VERSIONS= new String[][][] {
		new String[][] {new String[] {"GENCODE200510"}, new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"PANTRO21", "ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42"}, new String[] {"ENSEMBL43"}},	// cow
		new String[][] {new String[] {}},	// opossum

		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {}},	// fugu
		new String[][] {new String[] {}},	// tetraodon
		
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {}},	// mosquito
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		new String[][] {new String[] {"ENSEMBL42", "ENSEMBL43"}},
		
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}},
		new String[][] {new String[] {}}

	};
	public static final String[] SP_UCSC_CGI_STRINGS= new String[] {
		//&clade=vertebrate&org=Mouse&db=mm8
		// org=

			"Homo_sapiens;db=hg17",	// &db=hg18 newest, but for gencode..
			"Chimp",	// db=panTro
			"Mouse",	// db=mm8
			"Rat",		// &db=rn4
			"Dog",		// &db=canFam2
			"Cow",		//&db=bosTau2
			"Opossum",	//&db=monDom4
			
			"Chicken",	//&db=galGal3
			"X.+tropicalis",	//&db=xenTro2
			"Zebrafish",	//&db=danRer4
			"Fugu",			//&db=fr1
			"Tetraodon",	//&db=tetNig1
			
			"D.+melanogaster",	// db=dm2
			"A.+gambiae",	//&db=anoGam1
			"A.+mellifera",  //&db=anoGam1	// not in ensembl
			"C.+elegans",	//&db=ce2
			
			"S.+cerevisiae",	//&db=sacCer1
			"",
			"",
			"",
			""
			
	};
	
	public static final String[][] SP_UCSC_GENOME_VERSIONS= new String[][] {
		//&clade=vertebrate&org=Mouse&db=mm8
		// org=

		new String[] {"hg17", "hg18"},	// &db=hg18 newest, but for gencode..
		new String[] {"panTro"},
		new String[] {"mm8"},
		new String[] {"rn4"},
		new String[] {"canFam2"},
		new String[] {"bosTau2"},
		new String[] {"monDom4"},
		
		new String[] {"galGal3"},
		new String[] {"xenTro2"},
		new String[] {"danRer4"},	// zfish
		new String[] {"fr1"},	// fugu
		new String[] {"tetNig1"},	// tetraodon
		
		new String[] {"dm2"},	// droso
		new String[] {"anoGam1"},	// mosquito
		new String[] {"apiMel2"},  	// honeybee
		new String[] {"ce2"},	// worm
		
		new String[] {"sacCer1"},	// yeast
		new String[] {"ci2"},	// ciona
		new String[] {""},
		new String[] {""},
		new String[] {""}
		
		
};
	
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
			
			"WO",		// W09B6.4, Y59E9AR.9, ZK994.1, ZK994.3, ZK994.4, ... (no rule)
			"RDN"		// RDN58-1, snRNA, 15S-RNA, YAL016W, ... (no plan)
			
			// species missing
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
			
			"gallus_gallus",	// vertebrates
			"xenopus_tropicalis",
			"danio_rerio",
			"takifugu_rubripes",
			"tetraodon_nigroviridis",
			
			"drosophila_melanogaster",	// insects
			"anopheles_gambiae",
			"apis_mellifera",	// not in ensembl
			
			"caenorhabditis_elegans",	// deuterostome?
			"saccharomyces_cerevisiae",
			
			"ciona_intestinalis",	// chordata?
			"ciona_savignyi",
			"ornithorhynchus_anatinus",	// platypus, vertebrate?
			
			"arabidopsis_thaliana"
	};
	
	public static final String[] SP_NAMES_COMMON = new String[] { 
		"human", 
		"chimp", 
		"mouse", 
		"rat", 
		"dog", 
		"cow", 
		"opossum", 
		
		"chicken", "frog", "zebrafish", "fugu", "tetraodon", 
		
		"fruitfly", "mosquito", "honeybee", "worm", 
		
		"yeast", "seasquirt", "seqsquirt2", "platypus", "cress" };
	public static final String[][] SPECIFIC_ANNOTATIONS = new String[][] { 
		new String[] {"Known", "RefSeq", "EnsEmbl", "MGC"}, 
		new String[] {"chimp"}, 
		new String[] {"Known", "RefSeq", "EnsEmbl", "MGC"}, 
		new String[] {"rat"}, 
		new String[] {"dog"}, 
		new String[] {"cow"}, 
		new String[] {"opossum"}, 
		new String[] {"chicken"}, 
		new String[] {"frog"}, 
		new String[] {"zebrafish"}, 
		new String[] {"fugu"}, 
		new String[] {"tetraodon"}, 
		new String[] {"RefSeq", "FlyBase"}, 
		new String[] {"mosquito"}, 
		new String[] {"honeybee"}, 
		new String[] {"RefSeq", "WormBase"}, 
		new String[] {"yeast"}, 
		new String[] {"seasquirt"}, 
		new String[] {"seqsquirt2"}, 
		new String[] {"platypus"}, 
		new String[] {"cress"} 
	};
	public static final String[] SP_NAMES_ANDRE = new String[] { "human", "mouse", "rat", "dog", "cow", "chicken", "tetraodon"};
	// public static final String[] SP_NAMES_ABBREV= new String[] {"hsapiens", "mmusculus", "rnorvegicus"};
	public static final String[] SP_NAMES_METAZOA= new String[] {
			"human",	// mammals 
			"chimp",
			"mouse", 
			"rat",
			"dog",
			"cow",

			"chicken",	// other chordates
			"frog",
			"zebrafish",
			
			"fruitfly",
			"honeybee",		// not in ensembl
			"worm"

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
	
	public int spNumber= -1;

	
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

	public String getAbbreviatedName() {
		String s= SP_NAMES_BINOMIAL[spNumber];
		int p= s.indexOf("_");
		s= Character.toUpperCase(s.charAt(0))+ "."+ s.substring(p+1, s.length());
		return s;
	}
	
	public String getCommonName() {
		return SP_NAMES_COMMON[spNumber];
	}
	
	public String getNameFirst() {
		return SP_NAMES_BINOMIAL[spNumber].substring(0, SP_NAMES_BINOMIAL[spNumber].indexOf('_'));
	}	
	
	public static final String getAbbrevNameForBinomial(String binName) {
		return binName.charAt(0)+ binName.substring(binName.indexOf("_")+ 1);		
	}
	
	public String getNameAbbrev() {
		return getAbbrevNameForBinomial(SP_NAMES_BINOMIAL[spNumber]);
	}
	
	

	HashMap geneHash= null;

	public static final int getSpeNrForPrefix(String ensemblPFX) {
		
		ensemblPFX= ensemblPFX.toUpperCase();
		
		for (int i = SP_NAMES_ENS_PFX.length- 1; i >= 0 ; --i)	// iterate backw, human matches everything 
			if (ensemblPFX.startsWith(SP_NAMES_ENS_PFX[i]))
				return i;
		
		return -1;
	}
	
	public static final String getCommonNameForPrefix(String ensemblPFX) {
		
		int speNr= getSpeNrForPrefix(ensemblPFX);
		if (speNr>= 0)
			return SP_NAMES_COMMON[speNr];
		
		return null;
	}
	
	public static final String getGenomeVersionForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return SP_GENOME_VERSIONS[i][j][0];
		return null;
	}
	
	public static final int getGenomeVerNbForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return j;
		return -1;
	}
	
	public static final int getGenomeVerNb(String genomeVersion) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(genomeVersion))
						return j;
		return -1;
	}

	public static final int getSpeciesNumberForAnnotation(String build) {
		for (int i = 0; i < SP_ANNOTATION_VERSIONS.length; i++) 
			for (int j = 0; j < SP_ANNOTATION_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_ANNOTATION_VERSIONS[i][j].length; k++) 
					if (SP_ANNOTATION_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return i;
		return -1;
	}
	
	public static final int getSpeciesNumberForGenomeVersion(String genomeVer) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(genomeVer))
						return i;
		return -1;
	}

	public static final String getCommonNameForGenomeVersion(String build) {
		for (int i = 0; i < SP_GENOME_VERSIONS.length; i++) 
			for (int j = 0; j < SP_GENOME_VERSIONS[i].length; j++) 
				for (int k = 0; k < SP_GENOME_VERSIONS[i][j].length; k++) 
					if (SP_GENOME_VERSIONS[i][j][k].equalsIgnoreCase(build))
						return SP_NAMES_COMMON[i];
		return null;
	}
	
	public static final String decodeEnsemblPfx(String ensemblID) {
		Pattern patti= Pattern.compile("^(\\D{4,})\\d{11}.*");
		Matcher m= patti.matcher(ensemblID);
		if (m.matches())
			return m.group(1);
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
	
	public static final String getAbbrevNameForPrefix(String ensemblPFX) {
		
		for (int i = 0; i < SP_NAMES_ENS_PFX.length; i++) {
			if (SP_NAMES_ENS_PFX[i].equalsIgnoreCase(ensemblPFX))
				return getAbbrevNameForBinomial(SP_NAMES_BINOMIAL[i]);
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
	
	public static int getAnnotationNumber(String someName, String annoName) {
		String[] annotations= SPECIFIC_ANNOTATIONS[getSpeciesNumber(someName)];
		for (int i = 0; i < annotations.length; i++) {
			if (annotations[i].equals(annoName))
				return i;
		}
		return -1;
	}
	
	public static int getSpeciesNumberForBinomialName(String binName) {
		int i;
		for (i = 0; i < SP_NAMES_BINOMIAL.length; i++) 
			if (SP_NAMES_BINOMIAL[i].equals(binName))
				break;
		if (i>= SP_NAMES_BINOMIAL.length) {
			//System.err.println("Unknown name "+ binName);
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
	
	public String getDefaultGenomeVersion() {
		return Species.getDefaultGenomeVersion(getSpNumber());
	}
	
	public static String getDefaultGenomeVersion(int speNb) {
		return SP_GENOME_VERSIONS[speNb][0][0];
	}
	
	
	public Iterator getGeneIterator() {
		return geneHash.values().iterator();
	}

	public String[] getChromosomes() {
		Gene[] ge= getGenes();
		HashMap chrMap= new HashMap();
		for (int i = 0; i < ge.length; i++) {
			Object o= chrMap.get(ge[i].getChromosome());
			if (o== null)
				chrMap.put(ge[i].getChromosome(), ge[i].getChromosome());
		}
		
		return (String[]) Arrays.toField(chrMap.keySet());
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
	public String getAnnotationVersion() {
		return annotationVersion;
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
	public void setAnnotationVersion(String buildVersion) {
		this.annotationVersion = buildVersion;
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
	
	public void filter() {
		String cname= SP_NAMES_COMMON[getSpNumber()];
		Gene[] ge= getGenes();
		Vector hNames= new Vector();
		Vector otherNames= new Vector();
		int hCtr= 0;
		Vector remChrV= new Vector();
		int b4GenNb= ge.length, b4TrptNb= getTranscriptCount();
		for (int i = 0; i < ge.length; i++) {			
			String chr= ge[i].getChromosome();
			
			if (cname.equalsIgnoreCase("cow")) {
			if (ge[i].getChromosome().equals("30"))
				ge[i].setChromosome("X");	// ENSEMBL is funny
			}
			
				// accept: integers
			boolean isInt= true;
			try {  Integer.parseInt(chr);}
			catch (NumberFormatException e) {isInt= false;}
			if (isInt)
				continue;
			
				// accept: integers+ char
			String pfx= chr.substring(0, chr.length()- 1);
			isInt= true;
			try {  Integer.parseInt(pfx);}
			catch (NumberFormatException e) {isInt= false;}
			if (isInt) {
				if (chr.charAt(chr.length()-1)!= 'h')	// Droso, heterochromatin..
					continue;
				else {
					++hCtr;
					remChrV= Arrays.addUnique(remChrV, chr, String.CASE_INSENSITIVE_ORDER);
					hNames.add(ge[i].getGeneID());
					remove(ge[i], true);
					continue;
				}
			}
			
				// accept: roman numbers
			String[] romanNb= new String[] {
				"I", "II", "III", "IV", "V",
				"VI", "VII", "VIII", "IX", "X",
				"XI", "XII", "XIII", "XIV", "XV",
				"XVI", "XVII", "XVIII", "XIX", "XX",
				"XXI", "XXII", "XXIII", "XXIV", "XXV",
				"XXVI", "XXVII", "XXVIII", "XXIX", "XXX"
			};
			int j;
			for (j = 0; j < romanNb.length; j++) 
				if (chr.equalsIgnoreCase(romanNb[j]))
					break;
			if (j< romanNb.length)
				continue;
			
			// accept: extra chromosomes
			String[] extraChr= new String[] {
				"W", "X", "Y", "Z" 
			};
			for (j = 0; j < extraChr.length; j++) {
				if (chr.equalsIgnoreCase(extraChr[j]))
					break;
				if (chr.equalsIgnoreCase(extraChr[j]+"h")) {
					++hCtr;
					remChrV= Arrays.addUnique(remChrV, chr, String.CASE_INSENSITIVE_ORDER);
					hNames.add(ge[i].getGeneID());
					remove(ge[i], true);
					continue;
				}
			}
			if (j< extraChr.length)
				continue;
			
			
				// remove
			String[] noIDs= new String[] {
					"M", "Mt", "MtDNA",	// mitochondrial
					"U", "Un", "Unkn"	// un-...
			};
			for (j = 0; j < noIDs.length; j++) 				
				if (chr.equalsIgnoreCase(noIDs[j]))
					break;
			if (j< noIDs.length|| (chr.indexOf("_random")== (chr.length()- 7))) {
				++hCtr;
				remChrV= Arrays.addUnique(remChrV, chr, String.CASE_INSENSITIVE_ORDER);
				hNames.add(ge[i].getGeneID());
				remove(ge[i], true);
				continue;
			}

				// rest
			if (chr.startsWith("scaffold_")|| chr.startsWith("reftig_")|| 
					chr.startsWith("Contig")|| chr.startsWith("contig"))
				continue;
			else {
				++hCtr;
				remChrV= Arrays.addUnique(remChrV, chr, String.CASE_INSENSITIVE_ORDER);
				hNames.add(ge[i].getGeneID());
				remove(ge[i], true);
				continue;
			}
//			if (cname.equalsIgnoreCase("fruitfly")) {
//				String chr= ge[i].getChromosome();
//				if (Character.toLowerCase(chr.charAt(chr.length()- 1))== 'h') {
//					++hCtr;
//					hNames.add(ge[i].getGeneID());
//					remove(ge[i], true);
//				} else
//					otherNames.add(ge[i].getGeneID());
//			}
//			
//			
//			if (cname.equalsIgnoreCase("human")|| cname.equalsIgnoreCase("mouse")
//					|| cname.equalsIgnoreCase("cow")|| cname.equalsIgnoreCase("zebrafish")) {
//				String chr= ge[i].getChromosome();
//				try {
//					Integer.parseInt(chr);
//				} catch (NumberFormatException e) {
//						// zebrafish doesnt have X,Y, but igual
//					String[] otherChr= new String[] {"X", "Y"};	
//					int j;
//					for (j = 0; j < otherChr.length; j++) 
//						if (otherChr[j].equalsIgnoreCase(chr))
//							break;
//					if (j== otherChr.length) {
//						++hCtr;
//						hNames.add(ge[i].getGeneID());
//						remove(ge[i], true);
//					} else
//						otherNames.add(ge[i].getGeneID());
//				}
//			}
		}

		if (remChrV.size()> 0)
			System.out.print("Dont know what is");
		for (int i = 0; i < remChrV.size(); i++) 
			 System.out.print(" chr"+remChrV.elementAt(i));
		if (remChrV.size()> 0)
			System.out.println(", removing it");
		
		int ctrDup= 0;
		for (int i = 0; i < hNames.size(); i++) {
			int j;
			for (j = 0; j < otherNames.size(); j++) 
				if (((String) hNames.elementAt(i)).equalsIgnoreCase((String) otherNames.elementAt(j)))
					break;
			if (j< otherNames.size())
				++ctrDup;
		}		
		
		String percGen= Float.toString(((float) hCtr* 100f)/ b4GenNb);
		percGen= percGen.substring(0, percGen.indexOf('.'));
		if (hCtr> 0)
			System.out.println(SP_NAMES_COMMON[getSpNumber()]+" removed "+hCtr+" genes ("+percGen+"%) in strange DNA (heterochromatin, fragments, etc.), "+ctrDup+" of them with duplicates in chromosomal data.");
	}

	
	/**
	 * only splits genes, after transcript removal
	 * cannot unite them
	 */
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

	public int getSpNumber() {
		return spNumber;
	}

	/**
	 */
	public HashMap filterNonGTAGIntrons() {
		
		String[] chroms= getChromosomes();
		HashMap chrHash= new HashMap(chroms.length);
		for (int i = 0; i < chroms.length; i++) {
			int[] ratio= new int[2];
			ratio[0]= 0;
			ratio[1]= 0;
			chrHash.put(chroms[i], ratio);
		}
		
		int ctrTrpt= 0, ctrGene= 0;
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] trans= ge[i].getTranscripts();
			for (int j = 0; j < trans.length; j++) {
				if (trans[j].getNonGTAGIntron(chrHash)) {
					ge[i].removeTranscript(trans[j]);
					++ctrTrpt;
				}
			}
			if (ge[i].getTranscriptCount()< 1) {	// remove gene if there are no more transcripts
				remove(ge[i], true);
				++ctrGene;
			}
		}
		recluster();
		
		
		Object[] keys= chrHash.keySet().toArray();
		java.util.Arrays.sort(keys);
		int nonGTAG= 0, total= 0;
		for (int i = 0; i < keys.length; i++) {
			int[] ratio= (int[]) chrHash.get(keys[i]);
			String rat= Float.toString((ratio[0]* 100f)/ ratio[1]);
			rat= rat.substring(0, rat.indexOf('.')+ 2);
			if (keys.length< 50)
				System.out.print(keys[i]+": "+ratio[0]+"/"+ratio[1]+" ("+rat+"%)  ");
			nonGTAG+= ratio[0]; 
			total+= ratio[1];
		}
		String rat= Float.toString((nonGTAG* 100f)/ total);
		rat= rat.substring(0, rat.indexOf('.')+ 2);
		System.out.println("\n=> total "+nonGTAG+"/"+total+" ("+rat+"%)");
		System.out.println(SP_NAMES_COMMON[getSpNumber()]+" removed "+ctrTrpt+" transcripts, by this "+ctrGene+" complete genes" +
				" due to non-GT/AG introns.");
		return chrHash;
	}
	String annotationVersion = null;

	public String getGenomeVersion() {
		return genomeVersion;
	}

	public void setGenomeVersion(String genomeVersion) {
		this.genomeVersion = genomeVersion;
	}

	public static final String getAnnotation(String speName, String genomeVer, String annoName, String[] keywords) {
		for (int i = 0; keywords!= null&& i < keywords.length; i++) 
			keywords[i]= keywords[i].toUpperCase();
		
		speName= speName.toUpperCase();
		annoName= annoName.toUpperCase();
		String[] list= new File(Constants.SUBDIR_ANNOTATION).list();		
		Vector v= new Vector();
		for (int i = 0; i < list.length; i++) {
			String s= list[i].toUpperCase();
			if (s.contains(speName)&& s.contains(annoName)) {
				if (genomeVer!= null) {
					genomeVer= genomeVer.toUpperCase();
					if (!s.contains(genomeVer))
						continue;
				}
				if (keywords!= null) {
					int x;
					for (x = 0; x < keywords.length; x++) 
						if (s.contains(keywords[x]))
							break;
					if (x< keywords.length)
						v.add(new File(Constants.SUBDIR_ANNOTATION+ File.separator+ list[i]).getAbsolutePath());
				} else if (!new gphase.tools.File(list[i]).getExtension().contains("_"))
					v.add(new File(Constants.SUBDIR_ANNOTATION+ File.separator+ list[i]).getAbsolutePath());
			}
		}
		
		if (v.size()== 0) {
			System.err.println("No annotation found for "+speName+", "+annoName);
			return null;
		}
		if (v.size()== 1) 
			return (String) v.elementAt(0);
	
		System.out.println("Ambiguous information, press <CR> for fetching newest relase (according to 200x date).");
		try {
			System.in.read();
		} catch (IOException e) {
			e.printStackTrace();
		}
		int max= -1;
		int maxP= -1;
		for (int i = 0; i < v.size(); i++) {
			String s= (String) v.elementAt(i);
			int p= s.indexOf("200");// 200x
			if (p< 0) {
				System.err.println("No date available for "+s);
				continue;
			}
			int x= Integer.parseInt(s.substring(p+3, p+6));
			if (x> max) {
				max= x;
				maxP= i;
			}
		}
		
		return (String) v.elementAt(maxP);
	}
}
