/*
 * Created on Mar 8, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.Constants;
import gphase.NMDSimulator;

import gphase.algo.AlignmentGenerator;
import gphase.algo.AlignmentWrapper;
import gphase.db.EnsemblDBAdaptor;
import gphase.ext.ClustalWrapper;
import gphase.graph.SpliceGraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.RandomAccessFile;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Date;
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
public class Graph implements Serializable {
	
	static final long serialVersionUID = 6025857021729053707L;
	HashMap speciesHash= null;	// maps EnsemblPfx to Species
	public static final String GRAPH_ENCODE_REF_GENES= "check_enc.fasta";

	static HashMap buildMap= new HashMap(3);	// maps species builds (standard source) to date
	/*
	 * The Ensembl database names have the format:
	 * 
	 * <species><database type><release number><data version>
	 * 
	 * As an example, the name of the Human Core database for the Ensembl 32 release would be:
	 * homo_sapiens_core_32_35e
	 * A letter suffix to the data version indicates a change in data without a change in assembly; e.g. a new gene build.
	 * 
	 * --
	 * 
	 * Basically there are 3 identifiers for a genome version:
	 * NCBI build (e.g., 35), trivial name (hg17), date (200411)
	 * (which then are mapped to a Ensembl-version in ensemblDB)
	 */
//	static {
//		
//		HashMap map;
//		for (int i = 0; i < EnsemblDBAdaptor.SPECIES.length; i++) {
//			if (EnsemblDBAdaptor.SPECIES[i].contains("homo_sapiens")) {
//				map= new HashMap(1);
//				map.put(new Integer(35), new Integer(200405));	// NCBI build 35 = golden path hg17
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			}else if (EnsemblDBAdaptor.SPECIES[i].contains("mus_musculus")) {
//				map= new HashMap(2);
//				map.put(new Integer(33), new Integer(200405));	// NCBI build 33 = golden path mm5 (200405)
//				map.put(new Integer(34), new Integer(200503));	// NCBI build 34 = golden path mm6
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("rattus_norvegicus")) {
//				map= new HashMap(1);
//				map.put(new Integer(34), new Integer(200411));	// RGSC build 3.45 (no golden path)
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("canis_familiaris")) {
//				map= new HashMap(1);
//				map.put(new Integer(1), new Integer(200407));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("gallus_gallus")) {
//				map= new HashMap(1);
//				map.put(new Integer(1), new Integer(200402));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			} else if (EnsemblDBAdaptor.SPECIES[i].contains("pan_troglodytes")) {
//				map= new HashMap(1);
//				map.put(new Integer(3), new Integer(200311));	// UCSC
//				buildMap.put(EnsemblDBAdaptor.SPECIES[i], map);
//			}
//		}
//		
//	}

	/**
	 * @deprecated no longer used
	 */
	// assuming that homologs are always all included or none
	public Graph(Gene[] newGenes) {
		
		speciesHash= new HashMap();
		for (int i = 0; i < newGenes.length; i++) {
			
			if (speciesHash.get(newGenes[i].getSpecies().getEnsemblPrefix())== null) {
				addSpecies(newGenes[i].getSpecies());
				newGenes[i].getSpecies().geneHash= new HashMap();	// reinit geneMap
				newGenes[i].getSpecies().addGene(newGenes[i]);
			} else 
				newGenes[i].getSpecies().addGene(newGenes[i]);
			
			for (int j = 0; j < newGenes[i].getHomologies().length; j++) {
				if (speciesHash.get(newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]).getSpecies().getEnsemblPrefix())== null) {
					addSpecies(newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]).getSpecies());
					newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]).getSpecies().geneHash= new HashMap();	// reinit geneMap
					newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]).getSpecies().addGene(newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]));
				} else
					newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]).getSpecies().addGene(newGenes[i].getHomologies()[j].getOtherGene(newGenes[i]));
			} 
		}
		
	}
	
	public void recluster() {
		Iterator iter= speciesHash.values().iterator();
		while (iter.hasNext()) {
			Species spe= (Species) iter.next();
			spe.recluster();
		}
	}
	
	public void initTU() {
		Iterator iter= speciesHash.values().iterator();
		while (iter.hasNext()) {
			Species spe= (Species) iter.next();
			spe.initTU();
		}
	}
	
	/**
	 * Removes a gene and its homologies, but homologs only if they do not have any further homologies in the graph.
	 * @deprecated not in use.
	 * 
	 * @param rgene
	 * @return
	 */
	public Gene[] removePurge(Gene rgene) {
		
		Vector removed= new Vector();
		removed.add(rgene);
		
			// destroy the gene
		rgene.getSpecies().remove(rgene, true);	// eliminate from species
		if (rgene.getSpecies().getGeneNb()< 1)
			speciesHash.remove(decodeStableID(rgene.getStableID())[0]);
		
			// investigate homologs
		Gene[] homologs= rgene.getHomologGenes();
		for (int i = 0; i < homologs.length; i++) {
			homologs[i].removeHomology(rgene);	// eliminate gene from homolog list
			if (homologs[i].getHomologies().length< 1) {	// eliminate gene				
				homologs[i].getSpecies().remove(homologs[i], true);	// eliminate from species
				removed.add(homologs[i]);
			}
			if (homologs[i].getSpecies().getGeneNb()< 1) 
				speciesHash.remove(decodeStableID(homologs[i].getStableID())[0]);
		}
		
		Gene[] result= new Gene[removed.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (Gene) removed.elementAt(i);
		
		return result;
	}
	
	/**
	 * Removes a gene and kills also all its homologies as well as their further homologies.
	 * 
	 * @param rgene
	 * @return
	 */
	public Gene[] removeKill(Gene rgene) {
		
		Vector removed= new Vector();
		removed.add(rgene);
		
			// destroy the gene
		rgene.getSpecies().remove(rgene, true);	// eliminate from species
		if (rgene.getSpecies().getGeneNb()< 1)
			speciesHash.remove(decodeStableID(rgene.getStableID())[0]);
		
			// investigate homologs
		Gene[] homologs= rgene.getHomologGenes();
		for (int i = 0; homologs!= null&& i < homologs.length; i++) {
			homologs[i].removeHomology(homologs[i]);	// eliminate gene from homolog list
			homologs[i].getSpecies().remove(homologs[i], true);	// eliminate from species
			
			removed.add(homologs[i]);

			if (homologs[i].getSpecies().getGeneNb()< 1) 
				speciesHash.remove(decodeStableID(homologs[i].getStableID())[0]);
		}
		
		return Gene.toGeneArray(removed);
	}
	
	public static String[] decodeStableID(String stableID) {

		Pattern patty= Pattern.compile("(\\D*)(\\d*)");	// speciesID, stableID
		Matcher matty= patty.matcher(stableID);
		if (!matty.matches())
			System.out.println(stableID);
		
		return new String[] {matty.group(1).substring(0,matty.group(1).length()- 1), matty.group(2)};
	}		
		
	public static String decodeStableID(String stableID, boolean number) {
	
		String[] both= decodeStableID(stableID);
		
		if (number) 
			return both[1];
		else
			return both[0];
	}

	public Graph() {
		speciesHash= new HashMap();
	}

	public int countNonCodingLoci() {
		Iterator iter= speciesHash.values().iterator();
		int cnt= 0;
		while (iter.hasNext()) {
			Species spe= (Species) iter.next();
			cnt+= spe.countNonProteinCodingLoci();
		}
		return cnt;
	}
	public Graph(int nbSpecies) {
		speciesHash= new HashMap(nbSpecies);
	}
	
	public Graph(Species[] newSpecies) {
		this(newSpecies.length);
		addSpecies(newSpecies);
	}


	public static void main(String[] args) {

			// get Graph
		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
		Graph g= adaptor.getGraphAllHomologs();
		g.getSpeciesByName("mus_musculus").setBuildVersion(34);
		GraphHandler.writeOut(g, GraphHandler.getGraphAbsPath());
		
			// read in sequence
//		Exon e= g.getGene("ENSG00000198763").getExon("ENSE00001435686");
//		Gene g0= g.getGene("ENSMUST00000049656");
//		String seq= g.readSequence(g0);  
//		System.out.println(seq);
	}
	
	public Gene getGene(String stableGeneID, String speciesID) {
	
		Species spec= (Species) speciesHash.get(speciesID);		// get species
		if (spec== null)
			return null;

		return spec.getGene(stableGeneID);	
	}

	public Gene[] getGenes() {
		
		Iterator iter= speciesHash.values().iterator();
		Gene[] result= new Gene[0];
		while (iter.hasNext()) {
			Gene[] tmp= ((Species) iter.next()).getGenes();
			Gene[] old= result;
			result= new Gene[old.length+ tmp.length];
			for (int i = 0; i < old.length; i++) 
				result[i]= old[i];
			for (int i = 0; i < tmp.length; i++) 
				result[old.length+ i]= tmp[i];
		}
		
		return result;	
	}
	
	public Exon[] getExons(int region) {
		Vector v= new Vector();
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) 
			v= (Vector) gphase.tools.Arrays.addAll(v, ge[i].getExons(region));
		return (Exon[]) gphase.tools.Arrays.toField(v);
	}

	public Species getSpeciesByEnsemblPrefix(String speciesPfx) {
		
		if (Species.isValidEnsemblPrefix(speciesPfx))
			return (Species) speciesHash.get(speciesPfx);
		
		return getSpeciesByName(Species.SP_NAMES_BINOMIAL[11]);	// tetraodon, genoscan
	}
	
	
	/**
	 * Inefficient - traverses hash via loop. 
	 * Species names are to be given in the format "mus_musculus". 
	 * (Warning, contains an ineficient iteration over all species.)
	 * @param speciesName
	 * @return
	 */
	public Species getSpeciesByName(String binomialName) {
		
		if (speciesHash== null)
			return null;
		
		Iterator iter= speciesHash.values().iterator();
		Species spec;
		while (iter.hasNext()) {		// iterate hashmap
			spec= (Species) iter.next();
			String s= spec.getBinomialName();
			if (s.equalsIgnoreCase(binomialName))	// "mus_musculus"
				return spec;
		}
		
		return null;
	}	
	
	public int getTranscriptCount() {
		if (speciesHash== null) 
			return 0;
		Species[] ge= getSpecies();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) 
			cnt+= ge[i].getTranscriptCount();
		return cnt;
		
	}
	
	public int getClusterCountWithAS() {
		Gene[] ge= getGenes();
		int cnt= 0;
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as= ge[i].getASVariations(ASMultiVariation.FILTER_CODING_REDUNDANT);
			if (as!= null&& as.length> 0)
				++cnt;
		}
		
		return cnt;
	}
	
	public Gene[] getClustersWithAS() {
		Gene[] ge= getGenes();
		Vector geV= new Vector();
		for (int i = 0; i < ge.length; i++) {
			ASVariation[] as= ge[i].getASVariations(ASMultiVariation.FILTER_CODING_REDUNDANT);
			if (as!= null&& as.length> 0)
				geV.add(ge[i]);
		}
		
		return (Gene[]) gphase.tools.Arrays.toField(geV);
	}
	
	public Species getSpeciesByGeneID(String stableGeneID) {
		
			String[] idDec= decodeStableID(stableGeneID);		
			
			return getSpeciesByEnsemblPrefix(idDec[0]);		// get species
	}	
	
	
	public String countGenesTranscriptsExons() {
		
		String result= "[G,T,E]: ";
		Iterator specIter= speciesHash.values().iterator();
		while(specIter.hasNext()) {
			Iterator geneIter= ((Species) specIter.next()).geneHash.values().iterator();
			int geneCount= 0, transCount= 0, exonCount= 0;
			while(geneIter.hasNext()) {
				Gene gene= (Gene) geneIter.next();
				if (gene.getTranscripts()!= null)
					transCount+= gene.getTranscripts().length;
				if (gene.getExons()!= null)
					exonCount+= gene.getExons().length;
				++geneCount;
			}
			result+= "("+geneCount+","+transCount+","+exonCount+") ";
		}
		
		return result;
	}
	
	
	public boolean addGene(Gene newGene) {
		
			// get speciesHash
		Species spec= (Species) speciesHash.get(newGene.getSpecies().getCommonName()); 
		return spec.addGene(newGene);
	}
	
	/**
	 * 
	 * @param newSpecies
	 * @return <code>true</code> if new <code>Species[]</code> has been added successfully
	 */
	public boolean addSpecies(Species newSpecies) {
		
		if (newSpecies== null)
			return false;
		
		if (speciesHash== null)
			speciesHash= new HashMap();
		
		String key= newSpecies.getCommonName();
		if (speciesHash.get(key)!= null)
			return false;	// return false if key already occupied
		
		speciesHash.put(key, newSpecies);
		return true;
	}
	public Exon[] getExons() {

		Vector v= new Vector();
		Iterator iter= speciesHash.values().iterator();
		while(iter.hasNext()) {
			Species sp= (Species) iter.next();
			Exon[] ex= sp.getExons();
			for (int i = 0; i < ex.length; i++) {
				v.add(ex[i]);
			}
		}
		return (Exon[]) gphase.tools.Arrays.toField(v);
	}

	/**
	 * @return
	 */
	public Species[] getSpecies() {

		Collection c= speciesHash.values();
		Species[] sp= new Species[c.size()];
		Iterator it= c.iterator();
		int i= 0;
		while (it.hasNext())
			sp[i++]= (Species) it.next();
		return sp;
	}
	
	public SpliceSite[][] getSpliceSites(int regionType) {
		Species[] spe= getSpecies();
		Vector v1= new Vector();
		Vector v2= new Vector();
		getASVariations(ASMultiVariation.FILTER_NONE);	// init AS marks
		for (int i = 0; i < spe.length; i++) {
			SpliceSite[] ssites= spe[i].getSpliceSites(regionType, SpliceSite.ALTERNATE_SS);
		for (int j = 0; ssites!= null&& j < ssites.length; j++) 
				v1.add(ssites[j]);
			ssites= spe[i].getSpliceSites(regionType, SpliceSite.CONSTITUTIVE_SS);
			for (int j = 0; ssites!= null&& j < ssites.length; j++) 
				v2.add(ssites[j]);
		}
		SpliceSite[][] result= new SpliceSite[2][];
		result[0]= (SpliceSite[]) gphase.tools.Arrays.toField(v1);
		result[1]= (SpliceSite[]) gphase.tools.Arrays.toField(v2);
		return result;
	}
	
	public SpliceSite[] getSpliceSites(int spliceType, int regionType) {
		Vector v= new Vector();
		Gene[] ge= getGenes();
		for (int i = 0; i < ge.length; i++) {
			SpliceSite[] ss= ge[i].getSpliceSites(spliceType, regionType);
			v= (Vector) gphase.tools.Arrays.addAll(v, ss);
		}
		
		return (SpliceSite[]) gphase.tools.Arrays.toField(v);
	}

	public static String readSequence(Gene g) {
		
		String seq= readSequence(
				g.getSpecies(),
				g.getChromosome(),
				g.isForward(),
				g.getStart(),
				g.getEnd()
		);
		
		if (!g.isForward()) 
			seq= invertSequence(seq);
		return seq;
	}
	
	public static String readSequence(SpliceSite s) {
		String seq= readSequence(
				s.getGene().getSpecies(),
				s.getGene().getChromosome(),
				s.getGene().isForward(),
				s.getPos()- SpliceSite.SPLICE_SITE_FLANK_SIZE,
				s.getPos()+ SpliceSite.SPLICE_SITE_FLANK_SIZE
		);
		
		if (!s.getGene().isForward()) 
			seq= invertSequence(seq);
		return seq;
	}
	
	public static String readSequence(Exon e) {
		
		String seq= readSequence(
				e.getGene().getSpecies(),
				e.getGene().getChromosome(),
				e.getGene().isForward(),
				e.getStart(),
				e.getEnd()
		);
		
		if (!e.getGene().isForward()) 
			seq= invertSequence(seq);
		return seq;
	}

	public static String readSequence(Transcript t) {
		
		String seq= readSequence(
				t.getGene().getSpecies(),
				t.getGene().getChromosome(),
				t.getGene().isForward(),
				t.getStart(),
				t.getEnd()
		);
		
		if (!t.getGene().isForward()) 
			seq= invertSequence(seq);
		return seq;
	}
	
	public static String readSequence(DirectedRegion reg) {
		
		String seq= readSequence(
				reg.getSpecies(),
				reg.getChromosome(),
				reg.isForward(),
				reg.getStart(),
				reg.getEnd()
		);
		
			// neg strand already reversed in subroutine
		return seq;
	}

	public static String invertSequence(String seq) {
		
		seq= reverseSequence(seq);
		seq= complementarySequence(seq); 
		return seq;
	}
	
	public static String reverseSequence(String seq) {
		
		StringBuffer sb= new StringBuffer(seq.length());
		for (int i = (seq.length()-1); i >= 0; i--) 
			sb.append(seq.charAt(i));
		return sb.toString();
	}
	
	public static String complementarySequence(String seq) {
		
		StringBuffer sb= new StringBuffer(seq.length());
		for (int i = 0; i < seq.length(); i++) {
			
			char c= seq.charAt(i);
			boolean wasLow= Character.isLowerCase(c);
			c= Constants.NA_COMPL_IUPAC[Character.toUpperCase(c)- 65];
			if (wasLow)
				c= Character.toLowerCase(c);
			sb.append(c);
		}
		return sb.toString();
	}
	
	/**
	 * 
	 * @param speRealName
	 * @param chromosome
	 * @param forwardStrand
	 * @param start 1st position to be read
	 * @param end	1st position not to be read
	 * @return
	 */
	public static String readSequence(Species spe, String chromosome, boolean forwardStrand, int start, int end) {
			
			if (start< 0) {	// neg strand genes
				start= -start;
				end= -end;
			}
			if (start> end) {
				int h= start;
				start= end;
				end= h;
			}
				
		
			start--;
			if (chromosome.equals("MT"))	// correct Ensembl to GRIB jargon
				chromosome= "M";
			byte[] seq= new byte[end- start];
			try {
//				System.out.println(getSequenceDirectory(speRealName)+ File.separator+ "chr"+ chromosome+ Constants.CHROMOSOME_EXT);
//				System.out.println(start+"-"+end);
				RandomAccessFile raf= new RandomAccessFile(new File(
						getSequenceDirectory(spe)+ 
						File.separator+ "chr"+ chromosome+ Constants.CHROMOSOME_EXT),
						"r");
				String read= raf.readLine();
				while (!read.startsWith(">"))
					read= raf.readLine();
				int offset= read.length();
				read= "";
				while (read.trim().length()< 1)
					read= raf.readLine();
				int line= read.length();
				
				raf.seek(offset+1+start+ (start/line));	// fpointer is different from reading point!
				int pos= 0;
				int nextN= (line- (start%line));				// read (end of) first line
				while (pos+ nextN< seq.length) {		// full lines
					raf.readFully(seq,pos,nextN);
					raf.skipBytes(1);
					pos+= nextN;
					nextN= line;
				}
				raf.readFully(seq,pos,seq.length- pos);	// read start of last line
				raf.close();
			} catch (Exception e) {
				System.err.println("Problems reading sequence "+ e.getMessage());
				e.printStackTrace();
			}
			
			String s= new String(seq);
			if (!forwardStrand) 
				s= gphase.tools.Arrays.reverseComplement(s);
			return s;
		}

	public static String getSequenceDirectory(Species spe) {

		// Species spe= getSpeciesByName(realName);
		String realName= spe.getBinomialName();
//		HashMap hm= (HashMap) buildMap.get(realName);
//		Object o= hm.get(new Integer(spe.getBuildVersion()));
//		int dateID= ((Integer) ((HashMap) buildMap.get(realName))
//						.get(new Integer(spe.getBuildVersion()))).intValue();		// extract ID (date of build)
		
		
		String seqDirName= Character.toUpperCase(realName.charAt(0))+ "."
				+ realName.substring(realName.indexOf('_')+ 1);
		File speciesGenome= new File(Constants.DATA_DIR+ File.separator
				+ Constants.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		String[] list= speciesGenome.list();
		int i;
		for (i = 0; i < list.length; i++) {
			if (list[i].indexOf(new Integer(spe.getBuildVersion()).toString())>= 0)	// dateID
				break;
		}
		if (i> list.length- 1) {
			System.err.println("Build not found: "+ realName+" ("+spe.getBuildVersion()+")");
		}
		
		
		return speciesGenome.getAbsolutePath()+ File.separator
				+ list[i]+ File.separator+ Constants.CHROMOSOME_DIR;
	}

	public static String getSequenceDirectory_old(String realName) {
	
		String seqDirName= Character.toUpperCase(realName.charAt(0))+ "."
				+ realName.substring(realName.indexOf('_')+ 1);
		File speciesGenome= new File(Constants.DATA_DIR+ File.separator
				+ Constants.SEQUENCES_SUBDIR+ File.separator
				+ seqDirName);
		Pattern patty= Pattern.compile("^(\\D+)(\\d+).*");	// assuming that eg mouse "...mm5" is not relevant
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

	
	public void init() {
		//initExonSplicePattern();
		initExonHomology();
//		filter(g);
//		System.out.println("Graph filtered ---");
//		
//		g.init();
//		System.out.println("Graph inited ---");
	}
	
	void initExonSplicePattern() {
//		for (int i = 0; i < getSpecies().length; i++) {
//			for (int j = 0; j < getSpecies()[i].getGenes().length; j++) {
//				Gene g= getSpecies()[i].getGenes()[j];
//				switch (g.getTranscripts().length) {
//					case 0: System.err.println("Not transcribed gene "+ g.getStableID()
//							+ " ("+ getSpecies()[i].getRealName()+")");
//							break;
//					
//					case 1: for (int k = 0; k < g.getExons().length; k++) 
//								g.getExons()[k].addIdentica(Exon.TYPE_CO);
//							break;
//				
//					default: {
//						int nbTranscripts= g.getTranscripts().length;
//						
//					}
//				}
//			}
//		}
	}
	
	
	/**
	 *deprecated no longer in use
	 *
	 */
	public void initExonHomology() {
			
			for (int i = 0; i < getSpecies().length; i++) 			// species pair
				for (int j = i; j < getSpecies().length; j++) {		// incl same species (paralogs)
					Gene[] genesI= getSpecies()[i].getGenes();
	
					for (int k = 0; k < genesI.length; k++) {			// for all homolg gene pairs
						Gene[] hGenesIJ= genesI[k].getHomologGenes(getSpecies()[j]);
						for (int h = 0; h < hGenesIJ.length; h++) {
							
								// perform alignments
							for (int e = 0; e < genesI[k].getExons().length; e++) {	// for each exon pair
								for (int f = 0; f < hGenesIJ[h].getExons().length; f++) {
	
									String[] seqs= new String[] {
											readSequence(genesI[k].getExons()[e]), readSequence(hGenesIJ[h].getExons()[f])};
									String[] names= new String[] {genesI[k].getExons()[e].getExonID(), hGenesIJ[h].getExons()[f].getExonID()};
									
									String delme= AlignmentGenerator.writeOutTemp(names, seqs);
									ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(delme));
									new File(delme).delete();					// remove output
									//new File(cw.getOutputFName()).delete();	// readin result immediately
	
										// create hit
									PWHit hit= new PWHit(genesI[k].getExons()[e], hGenesIJ[h].getExons()[f]);
									hit.setScore(cw.getScore());
									hit.setAlignment(cw.getLayout());
									genesI[k].getExons()[e].addHit(hGenesIJ[h], hit);
									hGenesIJ[h].getExons()[f].addHit(genesI[k], hit);
								}
							}	// exon pair
							
								// assign relations
							for (int e = 0; e < genesI[k].getExons().length; e++) {	// for each exon pair
								
							}
							for (int f = 0; f < hGenesIJ[h].getExons().length; f++) {
							
						}
					}	// gene pair
					
				}
			}	// species pairs
		
		
		
			Gene[] as= getSpecies()[0].getGenes();	 
			for (int i = 0; i < as.length; i++) {
				Gene[] hGenes= as[i].getHomologGenes();
				System.out.println(i+": initing exon homology of "+as[i].getStableID());
	
					// gfx progress
				int sum= 0;
				for (int k = 0; k < hGenes.length; k++) {
					for (int j = 0; j < hGenes[k].getExons().length; j++) {				// align pw
						sum+= as[i].getHomologies().length- k;
					}				
				}
				char[] c= new char[sum];
				Arrays.fill(c, '.');
				System.out.println(c);
	
					// pw alignment of exons
				for (int k = 0; k < as[i].getHomologies().length; k++) {
					
					for (int j = 0; j < hGenes[k].getExons().length; j++) {				// align pw
	
						for (int ii = 0; ii < as[i].getExons().length; ii++) {				// with base gene 					
							
							String[] seqs= new String[] {readSequence(as[i].getExons()[ii]), readSequence(hGenes[k].getExons()[j])};
							String[] names= new String[] {as[i].getExons()[ii].getExonID(), hGenes[k].getExons()[j].getExonID()};
							
							PWHit hit= new PWHit(as[i].getExons()[ii], hGenes[k].getExons()[j]);
							ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
							hit.setScore(cw.getScore());
							hit.setAlignment(cw.getLayout());
							as[i].getExons()[ii].addHit(hGenes[k], hit);
							hGenes[k].getExons()[j].addHit(as[i], hit);
						}
						
						System.out.print('*');		// gfx output
						System.out.flush();
						
						for (int kk = (k+1); kk < as[i].getHomologies().length; kk++) {
							
							for (int ii = 0; ii < hGenes[kk].getExons().length; ii++) {				// with other homologs 					
								
								String[] seqs= new String[] {readSequence(hGenes[kk].getExons()[ii]), readSequence(hGenes[k].getExons()[j])};
								String[] names= new String[] {hGenes[kk].getExons()[ii].getExonID(), hGenes[k].getExons()[j].getExonID()};
								
								PWHit hit= new PWHit(hGenes[kk].getExons()[ii], hGenes[k].getExons()[j]);
								ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
								hit.setScore(cw.getScore());
								hit.setAlignment(cw.getLayout());
								hGenes[kk].getExons()[ii].addHit(hGenes[k], hit);
								hGenes[k].getExons()[j].addHit(hGenes[kk], hit);
							}
							System.out.print('*');
							System.out.flush();
						}
						
					}
				}
				System.out.println();
			}
			
	
			// output
			for (int i = 0; i < as.length; i++) {
				Gene[] hGenes= as[i].getHomologGenes();
				for (int j = 0; j < as[i].getExons().length; j++) { 
					for (int k = 0; k < as[i].getHomologies().length; k++) {
						// System.out.print(i+ " - "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologs()[k].getStableID()+ " : ");
						PWHit[] bestHits= getBRH(as[i].getExons()[j], hGenes[k], true);
						if (bestHits.length== 1) {							// add UBRH
							Exon e1= (Exon) bestHits[0].getObject1();
							Exon e2= (Exon) bestHits[0].getObject2();
							e1.addHomolog(e2.getGene(), e2);
							e2.addHomolog(e1.getGene(), e1);
						
						} else if (bestHits.length> 1)
							System.err.println("MBRHs ("+bestHits.length+") -"+i+ "- "+as[i].getExons()[j].getExonID()+" x "+ hGenes[k].getStableID());
						else 	// < 1, one exon not found
							System.err.println("noBRHs -"+i+ "- "+as[i].getExons()[j].getExonID()+" ("+as[i].getExons().length+") x "+ 
									hGenes[k].getStableID()+ " ("+ hGenes[k].getExons().length+")");
							
	//					System.out.print(getBRH(as[i].getExons()[j], as[i].getHomologs()[k], true).length);
	//					System.out.println(" BRH / "+ as[i].getExons()[j].getHits(as[i].getHomologs()[k]).length+ " hits.");
					}
				}
			}
			
			
		}

	/**
	 *deprecated no longer in use
	 *
	 */
	public void initSpliceSiteHomology() {
			
			for (int i = 0; i < getSpecies().length; i++) 			// species pair
				for (int j = i; j < getSpecies().length; j++) {		// incl same species (paralogs)
					Gene[] genesI= getSpecies()[i].getGenes();

					for (int k = 0; k < genesI.length; k++) {			// for all homolg gene pairs
						Gene[] hGenesIJ= genesI[k].getHomologGenes(getSpecies()[j]);
						if (hGenesIJ== null)
							continue;
						for (int h = 0; h < hGenesIJ.length; h++) {
							
								// perform alignments
							for (int e = 0; genesI[k].getSpliceSites()!= null&& e < genesI[k].getSpliceSites().length; e++) {	// for each exon pair
								for (int f = 0; hGenesIJ[h].getSpliceSites()!= null&& f < hGenesIJ[h].getSpliceSites().length; f++) {

									SpliceSite e1= genesI[k].getSpliceSites()[e];
									SpliceSite e2= hGenesIJ[h].getSpliceSites()[f];
									String[] seqs= new String[] {
											readSequence(e1.getGene().getSpecies(), e1.getGene().getChromosome(), e1.getGene().isForward(), Math.abs(e1.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e1.getPos())+ SpliceSite.DELTA_RANGE), 
											readSequence(e2.getGene().getSpecies(), e2.getGene().getChromosome(), e2.getGene().isForward(), Math.abs(e2.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e2.getPos())+ SpliceSite.DELTA_RANGE)};
									String[] names= new String[] {e1.getExons()[0].getExonID()+(e1.isDonor()?"^":"-"), e2.getExons()[0].getExonID()+(e2.isDonor()?"^":"-")};
									
									String delme= AlignmentGenerator.writeOutTemp(names, seqs);
									ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(delme));
									new File(delme).delete();					// remove output
									//new File(cw.getOutputFName()).delete();	// readin result immediately

										// create hit
									PWHit hit= new PWHit(genesI[k].getExons()[e], hGenesIJ[h].getExons()[f]);
									hit.setScore(cw.getScore());
									hit.setAlignment(cw.getLayout());
									e1.addHit(hGenesIJ[h], hit);
									e2.addHit(genesI[k], hit);
								}
							}	// ss pair
							
					}	// gene pair
					
				}
			}	// species pairs

	
			// output
			for (int i = 0; i < getSpecies().length; i++) 			// species pair
				for (int j = i; j < getSpecies().length; j++) {		// incl same species (paralogs)
					Gene[] genesI= getSpecies()[i].getGenes();

					for (int k = 0; k < genesI.length; k++) {			// for all homolg gene pairs
						Gene[] hGenesIJ= genesI[k].getHomologGenes(getSpecies()[j]);
						if (hGenesIJ== null)
							continue;
						for (int h = 0; h < hGenesIJ.length; h++) {
							for (int e = 0; genesI[k].getSpliceSites()!= null&& e < genesI[k].getSpliceSites().length; e++) {	// for each exon pair
								PWHit[] bestHits= getBRH(genesI[k].getSpliceSites()[e], hGenesIJ[h], true);
								if (bestHits.length== 1) {							// add UBRH
									SpliceSite e1= (SpliceSite) bestHits[0].getObject1();
									SpliceSite e2= (SpliceSite) bestHits[0].getObject2();
									e1.addHomolog(e2.getGene(), e2);
									e2.addHomolog(e1.getGene(), e1);
								
								} else if (bestHits.length> 1)
									System.err.println("MBRHs ("+bestHits.length+") -"+i+ "- "+genesI[k].getSpliceSites()[e].getExons()[0].getExonID()+" x "+ 
											hGenesIJ[h].getStableID());
								else 	// < 1, one exon not found
									System.err.println("noBRHs -"+i+ "- "+genesI[k].getSpliceSites()[e].getExons()[0].getExonID()+" ("+genesI[k].getSpliceSites().length+") x "+ 
											hGenesIJ[h].getStableID()+ " ("+ hGenesIJ[h].getSpliceSites().length+")");

								
							}			
					}
				}
			}
	}
		
		

	void initSpliceSiteHomology_alt() {
			Gene[] as= getSpecies()[0].getGenes();	 
			for (int i = 0; i < as.length; i++) {
				Gene[] hGenes= as[i].getHomologGenes();
				System.out.println(i+": initing exon homology of "+as[i].getStableID());
	
					// gfx progress
				int sum= 0;
				for (int k = 0; k < hGenes.length; k++) {
					for (int j = 0; j < hGenes[k].getExons().length; j++) {				// align pw
						sum+= as[i].getHomologies().length- k;
					}				
				}
				char[] c= new char[sum];
				Arrays.fill(c, '.');
				System.out.println(c);
	
					// pw alignment of exons
				for (int k = 0; k < as[i].getHomologies().length; k++) {
					
					for (int j = 0; j < hGenes[k].getSpliceSites().length; j++) {				// align pw
	
						for (int ii = 0; ii < as[i].getSpliceSites().length; ii++) {				// with base gene 					
							
							SpliceSite e1= as[i].getSpliceSites()[ii];
							SpliceSite e2= hGenes[k].getSpliceSites()[j];
							String[] seqs= new String[] {
									readSequence(e1.getGene().getSpecies(), e1.getGene().getChromosome(), e1.getGene().isForward(), Math.abs(e1.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e1.getPos())+ SpliceSite.DELTA_RANGE), 
									readSequence(e2.getGene().getSpecies(), e2.getGene().getChromosome(), e2.getGene().isForward(), Math.abs(e2.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e2.getPos())+ SpliceSite.DELTA_RANGE)};
							String[] names= new String[] {e1.getExons()[0].getExonID()+(e1.isDonor()?"^":"-"), e2.getExons()[0].getExonID()+(e2.isDonor()?"^":"-")};
							
							PWHit hit= new PWHit(as[i].getExons()[ii], hGenes[k].getExons()[j]);
							ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
							hit.setScore(cw.getScore());
							hit.setAlignment(cw.getLayout());
							as[i].getSpliceSites()[ii].addHit(hGenes[k], hit);
							hGenes[k].getSpliceSites()[j].addHit(as[i], hit);
						}
						
						System.out.print('*');		// gfx output
						System.out.flush();
						
						for (int kk = (k+1); kk < as[i].getHomologies().length; kk++) {
							
							for (int ii = 0; ii < hGenes[kk].getSpliceSites().length; ii++) {				// with other homologs 					
								
								SpliceSite e1= hGenes[k].getSpliceSites()[j];
								SpliceSite e2= hGenes[kk].getSpliceSites()[ii];
								String[] seqs= new String[] {
										readSequence(e1.getGene().getSpecies(), e1.getGene().getChromosome(), e1.getGene().isForward(), Math.abs(e1.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e1.getPos())+ SpliceSite.DELTA_RANGE), 
										readSequence(e2.getGene().getSpecies(), e2.getGene().getChromosome(), e2.getGene().isForward(), Math.abs(e2.getPos())- SpliceSite.DELTA_RANGE, Math.abs(e2.getPos())+ SpliceSite.DELTA_RANGE)};
								String[] names= new String[] {e1.getExons()[0].getExonID()+(e1.isDonor()?"^":"-"), e2.getExons()[0].getExonID()+(e2.isDonor()?"^":"-")};
								
								PWHit hit= new PWHit(hGenes[kk].getExons()[ii], hGenes[k].getExons()[j]);
								ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
								hit.setScore(cw.getScore());
								hit.setAlignment(cw.getLayout());
								hGenes[kk].getSpliceSites()[ii].addHit(hGenes[k], hit);
								hGenes[k].getSpliceSites()[j].addHit(hGenes[kk], hit);
							}
							System.out.print('*');
							System.out.flush();
						}
						
					}
				}
				System.out.println();
			}
			
	
			// output
			for (int i = 0; i < as.length; i++) {
				Gene[] hGenes= as[i].getHomologGenes();
				for (int j = 0; j < as[i].getSpliceSites().length; j++) { 
					for (int k = 0; k < as[i].getHomologies().length; k++) {
						// System.out.print(i+ " - "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologs()[k].getStableID()+ " : ");
						PWHit[] bestHits= getBRH(as[i].getSpliceSites()[j], hGenes[k], true);
						if (bestHits.length== 1) {							// add UBRH
							Exon e1= (Exon) bestHits[0].getObject1();
							Exon e2= (Exon) bestHits[0].getObject2();
							e1.addHomolog(e2.getGene(), e2);
							e2.addHomolog(e1.getGene(), e1);
						
						} else if (bestHits.length> 1)
							System.err.println("MBRHs ("+bestHits.length+") -"+i+ "- "+as[i].getExons()[j].getExonID()+" x "+ hGenes[k].getStableID());
						else 	// < 1, one exon not found
							System.err.println("noBRHs -"+i+ "- "+as[i].getExons()[j].getExonID()+" ("+as[i].getExons().length+") x "+ 
									hGenes[k].getStableID()+ " ("+ hGenes[k].getExons().length+")");
							
	//					System.out.print(getBRH(as[i].getExons()[j], as[i].getHomologs()[k], true).length);
	//					System.out.println(" BRH / "+ as[i].getExons()[j].getHits(as[i].getHomologs()[k]).length+ " hits.");
					}
				}
			}
			
			
		}

	/**
	 * @deprecated intransparent loop structure
	 *
	 */
	void initExonHomology_old() {
			
/*		Gene[] as= getSpecies()[0].getGenes();	 
		for (int i = 0; i < as.length; i++) {
			System.out.println(i+": initing exon homology of "+as[i].getStableID());

				// gfx progress
			int sum= 0;
			for (int k = 0; k < as[i].getHomologies().length; k++) {
				for (int j = 0; j < as[i].getHomologies()[k].getExons().length; j++) {				// align pw
					sum+= as[i].getHomologies().length- k;
				}				
			}
			char[] c= new char[sum];
			Arrays.fill(c, '.');
			System.out.println(c);

				// pw alignment of exons
			for (int k = 0; k < as[i].getHomologies().length; k++) {
				
				for (int j = 0; j < as[i].getHomologies()[k].getExons().length; j++) {				// align pw

					for (int ii = 0; ii < as[i].exons.length; ii++) {				// with base gene 					
						
						String[] seqs= new String[] {readSequence(as[i].exons[ii]), readSequence(as[i].getHomologies()[k].getExons()[j])};
						String[] names= new String[] {as[i].exons[ii].getStableID(), as[i].getHomologies()[k].getExons()[j].getStableID()};
						
						PWHit hit= new PWHit(as[i].exons[ii], as[i].getHomologies()[k].getExons()[j]);
						ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
						hit.setScore(cw.getScore());
						hit.setAlignment(cw.getLayout());
						as[i].exons[ii].addHit(as[i].getHomologies()[k], hit);
						as[i].getHomologies()[k].getExons()[j].addHit(as[i], hit);
					}
					
					System.out.print('*');		// gfx output
					System.out.flush();
					
					for (int kk = (k+1); kk < as[i].getHomologies().length; kk++) {
						
						for (int ii = 0; ii < as[i].getHomologies()[kk].exons.length; ii++) {				// with other homologs 					
							
							String[] seqs= new String[] {readSequence(as[i].getHomologies()[kk].exons[ii]), readSequence(as[i].getHomologies()[k].getExons()[j])};
							String[] names= new String[] {as[i].getHomologies()[kk].exons[ii].getStableID(), as[i].getHomologies()[k].getExons()[j].getStableID()};
							
							PWHit hit= new PWHit(as[i].getHomologies()[kk].exons[ii], as[i].getHomologies()[k].getExons()[j]);
							ClustalWrapper cw= ((ClustalWrapper)AlignmentGenerator.alignClustal(AlignmentGenerator.writeOutTemp(names, seqs)));
							hit.setScore(cw.getScore());
							hit.setAlignment(cw.getLayout());
							as[i].getHomologies()[kk].exons[ii].addHit(as[i].getHomologies()[k], hit);
							as[i].getHomologies()[k].getExons()[j].addHit(as[i].getHomologies()[kk], hit);
						}
						System.out.print('*');
						System.out.flush();
					}
					
				}
			}
			System.out.println();
		}
		

		// output
		for (int i = 0; i < as.length; i++) {
			for (int j = 0; j < as[i].getExons().length; j++) { 
				for (int k = 0; k < as[i].getHomologies().length; k++) {
					// System.out.print(i+ " - "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologs()[k].getStableID()+ " : ");
					PWHit[] bestHits= getBRH(as[i].getExons()[j], as[i].getHomologies()[k], true);
					if (bestHits.length== 1) {							// add UBRH
						Exon e1= (Exon) bestHits[0].getObject1();
						Exon e2= (Exon) bestHits[0].getObject2();
						e1.addHomolog(e2.getGene(), e2);
						e2.addHomolog(e1.getGene(), e1);
					
					} else if (bestHits.length> 1)
						System.err.println("MBRHs ("+bestHits.length+") -"+i+ "- "+as[i].getExons()[j].getStableID()+" x "+ as[i].getHomologies()[k].getStableID());
					else 	// < 1, one exon not found
						System.err.println("noBRHs -"+i+ "- "+as[i].getExons()[j].getStableID()+" ("+as[i].getExons().length+") x "+ 
								as[i].getHomologies()[k].getStableID()+ " ("+ as[i].getHomologies()[k].getExons().length+")");
						
//					System.out.print(getBRH(as[i].getExons()[j], as[i].getHomologs()[k], true).length);
//					System.out.println(" BRH / "+ as[i].getExons()[j].getHits(as[i].getHomologs()[k]).length+ " hits.");
				}
			}
		}
*/		
		
	}
	
	PWHit[] getBRH(Exon e, Gene g, boolean max) {
	
		Vector result= new Vector(e.getHits().length);
		
		PWHit[] hits= e.getHits(g);
		int best= 0;
		for (int i = 0; i < hits.length; i++) {
			if (max&& hits[i].getScore()>= best) {			// for scores, look for max
				best= hits[i].getScore();
				PWHit[] counterHits= ((Exon) hits[i].getOtherObject(e)).getHits(e.getGene());
				int j;
				for (j = 0; j < counterHits.length; j++) 
					if (counterHits[j].getScore()> best)
						break;
				if (j>= counterHits.length)		// BRH spotted
					result.add(hits[i]);
			} else if ((!max)&& hits[i].getScore()<= best) {	// for costs, look for min
				best= hits[i].getScore();
				PWHit[] counterHits= ((Exon) hits[i].getOtherObject(e)).getHits();
				int j;
				for (j = 0; j < counterHits.length; j++) 
					if (counterHits[j].getScore()< best)
						break;
				if (j>= counterHits.length)		// BRH spotted
					result.add(hits[i]);
			} 
		}
		
		PWHit[] res= new PWHit[result.size()];
		for (int i = 0; i < res.length; i++) 
			res[i]= (PWHit) result.elementAt(i);
		
		return res;
	}

	PWHit[] getBRH(SpliceSite e, Gene g, boolean max) {

		Vector result= new Vector(e.getHits().length);
		
		PWHit[] hits= e.getHits(g);
		int best= 0;
		for (int i = 0; i < hits.length; i++) {
			if (max&& hits[i].getScore()>= best) {			// for scores, look for max
				best= hits[i].getScore();
				PWHit[] counterHits= ((Exon) hits[i].getOtherObject(e)).getHits(e.getGene());
				int j;
				for (j = 0; j < counterHits.length; j++) 
					if (counterHits[j].getScore()> best)
						break;
				if (j>= counterHits.length)		// BRH spotted
					result.add(hits[i]);
			} else if ((!max)&& hits[i].getScore()<= best) {	// for costs, look for min
				best= hits[i].getScore();
				PWHit[] counterHits= ((Exon) hits[i].getOtherObject(e)).getHits();
				int j;
				for (j = 0; j < counterHits.length; j++) 
					if (counterHits[j].getScore()< best)
						break;
				if (j>= counterHits.length)		// BRH spotted
					result.add(hits[i]);
			} 
		}
		
		PWHit[] res= new PWHit[result.size()];
		for (int i = 0; i < res.length; i++) 
			res[i]= (PWHit) result.elementAt(i);
		
		return res;
	}


	/**
	 * filters off not alternatively spliced genes.
	 *
	 */
	public void filterSingleTranscriptGenes() {
		System.out.println("filter single transcript genes");
		int rgen= 0;
		Species[] spec= getSpecies();
		for (int i = 0; i < spec.length; i++) {
			Gene[] ge= spec[i].getGenes();	// make static, method always generates newly -> bad when removing
			for (int j = 0; j < ge.length; j++) 
				if (ge[j].getTranscripts().length< 2) {
					spec[i].remove(ge[j], true);
					++rgen;
				}
		}
		System.out.println("==> "+rgen+"nonAS genes filtered out.");
	}
	
	public ASVariation[][] search(int[][] pattern) {
		
		ASVariation[][] vars= getASVariations();
		Vector result= new Vector();
		for (int i = 0; i < vars.length; i++) {
			if (vars[i][0].contains(pattern))
				result.add(vars[i]);
		}
				
		ASVariation[][] found= null;
		try {
			found= (ASVariation[][]) gphase.tools.Arrays.toField(result);
		} catch (ClassCastException e) {
			System.err.println("Pattern not found!");
		}
		return found;
	}
	
	/**
	 * filters out Ensembl crap (too small introns, splice donor/acceptors at the same position, ..)
	 * @param graph
	 * @return
	 */
	public void filterNonsense() {
			System.out.println("filter nonsense");
			
			int rgen= 0; 
			int rtran= 0; 
			Species[] spec= getSpecies();
			for (int i = 0; i < spec.length; i++) {
				Gene[] ge= spec[i].getGenes();	// make static, method always generates newly -> bad when removing
				for (int j = 0; j < ge.length; j++) {
					boolean removeGene= false;
					Transcript[] tra= ge[j].getTranscripts();
					for (int k = 0; !removeGene && k < tra.length; k++) {

							// check for removing transcript
						boolean removeTranscript= false;
						for (int m = 0; !removeGene && m < tra[k].getExons().length; m++) 
							if (((m< tra[k].getExons().length- 1)&& 
									(tra[k].getExons()[m+1].getStart()- tra[k].getExons()[m].getEnd()< Constants.MIN_INTRON_LENGTH))
								|| tra[k].getExons()[m].getEnd()- tra[k].getExons()[m].getStart()<= Constants.MIN_EXON_LENGTH) {
								removeTranscript= true;
								break;
							}
						if (removeTranscript) {
							ge[j].removeTranscript(tra[k]);
							++rtran;
							if (ge[j].getTranscripts().length< 1) {
								removeGene= true;
								break;
							}
							continue;	// next transcript
						}

						// check for removing complete gene
						for (int m = 0; !removeGene && m < tra[k].getExons().length; m++) 
							for (int n = m+1; !removeGene && n < tra[k].getExons().length; n++) 
								if ((tra[k].getExons()[m].getDonor()!= null&&
									tra[k].getExons()[n].getAcceptor()!= null&&
									tra[k].getExons()[m].getDonor().getPos()==
										tra[k].getExons()[n].getAcceptor().getPos())||
									(tra[k].getExons()[m].getAcceptor()!= null&&
									tra[k].getExons()[n].getDonor()!= null&&
									tra[k].getExons()[m].getAcceptor().getPos()==
										tra[k].getExons()[n].getDonor().getPos())) {
												
										removeGene= true;
										break;
								}
							

					}
					
					if (removeGene) {
						spec[i].remove(ge[j], true);
						++rgen;
					}
				}
			}
		
			System.out.println("==> Nonsense removed "+rgen+" genes and "+rtran+ " transcripts.");
		
		}
	
	/**
	 */
	public void filterNonCodingTranscripts() {
		System.out.println("filter non-coding transcripts");
		Iterator iter= speciesHash.values().iterator();
		while (iter.hasNext()) {
			Species spec= (Species) iter.next();
			Gene[] ge= spec.getGenes();
			for (int i = 0; i < ge.length; i++) {
				Transcript[] trans= ge[i].getTranscripts();
				for (int j = 0; j < trans.length; j++) {
					if (trans[j].isNonCoding())
						ge[i].removeTranscript(trans[j]);
				}
				if (ge[i].getTranscriptCount()< 1)	// remove gene if there are no more transcripts
					spec.remove(ge[i], true);
			}
		}
		recluster();
	}

	/**
	 */
	public void filterNMDTranscripts() {
		System.out.println("filter NMD transcripts");
		Iterator iter= speciesHash.values().iterator();
		NMDSimulator nmd;
		while (iter.hasNext()) {
			Species spec= (Species) iter.next();
			Gene[] ge= spec.getGenes();
			for (int i = 0; i < ge.length; i++) {
				Transcript[] trans= ge[i].getTranscripts();
				for (int j = 0; j < trans.length; j++) {					
					nmd= new NMDSimulator(trans[j]);
					Translation tln;
					if (trans[j].getTranslations()== null|| trans[j].getTranslations()[0]== null)
						tln= trans[j].findHavanaORF();	// *
					else
						tln= trans[j].getTranslations()[0];
					if (tln!= null&& nmd.isNMD(tln))	// (tln== null|| nmd.isNMD(tln)) for removing trpts wo ORF
						ge[i].removeTranscript(trans[j]);
				}
				if (ge[i].getTranscriptCount()< 1)	// remove gene if there are no more transcripts
					spec.remove(ge[i], true);
			}
		}
		recluster();
	}
	
	public void filterCodingTranscripts() {
		System.out.println("filter coding transcripts");
		Iterator iter= speciesHash.values().iterator();
		while (iter.hasNext()) {
			Species spec= (Species) iter.next();
			Gene[] ge= spec.getGenes();
			for (int i = 0; i < ge.length; i++) {
				Transcript[] trans= ge[i].getTranscripts();
				for (int j = 0; j < trans.length; j++) {
					if (!trans[j].isNonCoding())
						ge[i].removeTranscript(trans[j]);
				}
				if (ge[i].getTranscriptCount()< 1)	// remove gene if there are no more transcripts
					spec.remove(ge[i], true);
			}
		}
		recluster();
	}
	
	/**
	 * Filters to leave only genes in that are within one of the specified regions.
	 * All regions have to be from the same species.
	 * 
	 * @param regions
	 */
	public GeneHomology[] filterForRegions(Region[] regions) {
		
		if (regions== null|| regions.length< 1)
			return null;
		
		Species[] spec= getSpecies();
		Gene[] tgenes= regions[0].getSpecies().getGenes();
		if (tgenes== null|| tgenes.length< 1)
			return null;
		
		for (int i = 0; i < regions.length; i++)	// check 
			if (regions[i].getSpecies()!= tgenes[0].getSpecies())
				System.err.println("Regions with not homolg species "+ 
						tgenes[0].getSpecies()+ " <> "+ regions[i].getSpecies());

		Vector xHomol= new Vector();
		for (int j = 0; j < tgenes.length; j++) {
			int i;
			for (i = 0; i < regions.length; i++) 
				if (regions[i].contains(tgenes[j]))
					break;
			
			if (i>= regions.length) 	// not found
				removeKill(tgenes[j]);
			else {
				GeneHomology[] h= tgenes[j].getHomologies();
				for (int k = 0; k < h.length; k++) 
					xHomol.add(h[k]);
			}
		}
		
		return GeneHomology.toGeneHomologyArray(xHomol);
	}

	
	public ASVariation[][][] getASEvents() {
		
		if (asAbstractClasses == null) {
			Vector abstractV = new Vector();
			for (int i = 0; i < getASVariations().length; i++) {
				int j;
				for (j = 0; j < abstractV.size(); j++) {
					if (((ASVariation[]) ((Vector) abstractV.elementAt(j)).elementAt(0))[0].toStringASEvents().equals(
							getASVariations()[i][0].toStringASEvents()))	// merge ?						
					((Vector) abstractV.elementAt(j)).add(getASVariations()[i]);
					break;
				}
				if (j>= abstractV.size()) {		// add new super class
					Vector v= new Vector();
					v.add(getASVariations()[i]);
					abstractV.add(v);
				}
			}
			
			asAbstractClasses= (ASVariation[][][]) gphase.tools.Arrays.toField(abstractV);
		}

		return asAbstractClasses;
	}
	
	/**
	 * 
	 * @param newSpecies
	 * @return <code>true</code> if new <code>Species[]</code> has been added successfully
	 */
	public boolean addSpecies(Species[] newSpecies) {
		
		if (newSpecies== null)
			return false;
		
		boolean result= true;
		for (int i = 0; i < newSpecies.length; i++)  
			result&= addSpecies(newSpecies[i]);
		
		return result;
	}

	public ASMultiVariation[][] getASMultiVariations(int k) {
		Gene[] ge= getGenes();
		int x= 0; 
		Vector varsV= new Vector();
		for (x = 0; x < ge.length; x++)  {
			SpliceGraph gr= new SpliceGraph(ge[x].getTranscripts());
			gr.init();
			ASMultiVariation[] vars= gr.getMultiVariations(k);
			
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				int j;
				for (j = 0; j < varsV.size(); j++) {
					Vector tmpV= (Vector) varsV.elementAt(j);
					ASMultiVariation v= (ASMultiVariation) tmpV.elementAt(0);
					if (v.toString().equals(vars[i].toString())) {
						tmpV.add(vars[i]);
						break;
					}
				}
				if (j== varsV.size()) {	// not found
					Vector tmpV= new Vector();
					tmpV.add(vars[i]);
					varsV.add(tmpV);
				}
			}
		}
		
		return (ASMultiVariation[][]) gphase.tools.Arrays.toField(varsV);
	}
	
	public ASMultiVariation[][] getASMultiVariations_alternativ() {
		
		Gene[] ge= getGenes();
		Vector asMV= new Vector();
		for (int i = 0; i < ge.length; i++) {	// all genes
			ASMultiVariation[] as= ge[i].getASMultiVariations();
			asMV.add(as);
		}
		
		ASMultiVariation[][] multiV= null;
		try {
			multiV= (ASMultiVariation[][]) gphase.tools.Arrays.toField(asMV);
		} catch (ClassCastException e) {
			System.err.println("No AS Classes!");
		}
		
		return multiV;
	}

	public ASMultiVariation[][] getASMultiVariations2() {
		
		Gene[] ge= getGenes();
		Vector mvClasses= new Vector();
		for (int i = 0; i < ge.length; i++) {	// all genes
			ASMultiVariation[] as= ge[i].getASMultiVariations2();
			if (as== null)
				continue;
			for (int j = 0; j < as.length; j++) {	// sort MVars into classes
				int k;
				for (k = 0; k < mvClasses.size(); k++) {
					Vector tmpClass= (Vector) mvClasses.elementAt(k);
					if (((ASMultiVariation) tmpClass.elementAt(0)).toString().equals(as[j].toString())) {
						tmpClass.add(as[j]);
						break;
					}
				}
				if (k>= mvClasses.size()) {
					Vector v= new Vector();
					v.add(as[j]);
					mvClasses.add(v);
				}
			}
		}
		
		ASMultiVariation[][] multiV= null;
		try {
			multiV= (ASMultiVariation[][]) gphase.tools.Arrays.toField(mvClasses);
		} catch (ClassCastException e) {
			System.err.println("No multi variations!");
		}
		
		return multiV;
	}
	
	public ASVariation[][] getASVariations(int filter) {
		
		asClasses= null;
		System.out.println("Retrieving pw AS variations, filtered "+ ASMultiVariation.FILTER_TO_STRING[filter]);
		if (asClasses == null) {
			Gene[] ge= getGenes();
			asVariations= 0;
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

	public int countASVariations() {
		Gene[] ge= getGenes();
		int nb= 0;
		for (int i = 0; i < ge.length; i++) {
			ASMultiVariation[] as= ge[i].getASMultiVariations();
			for (int j = 0; as!= null&& j < as.length; j++) {
				nb+= as[j].getASVariationsHierarchicallyFiltered().length;
			}
		}		
	
		return nb;
	}

	ASVariation[][] asClasses;
	ASVariation[][][] asAbstractClasses;
	
	int asVariations;
}
