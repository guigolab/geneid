/*
 * Created on Mar 31, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase;

import java.io.File;
import java.util.Arrays;
import java.util.Date;

import qalign.OSChecker;

/**
 * 
 * 
 * @author micha
 */
public class Constants {

	public static final String getDateString() {
		return "["+ new Date(System.currentTimeMillis()).toString()+ "]";
	}
	
	static final String[] trees= {"", "human", "(human, mouse)", "(human (mouse rat))"};
	
	public static String DATA_DIR;	// where there is a lot of space..
	static {
		if (OSChecker.isSunOS())
			DATA_DIR= "/vol/cluster-data/micha/gphase";
		else if (OSChecker.isLinux())
			DATA_DIR= "/home/msammeth";
	}
	public final static String SEQUENCES_SUBDIR= "genomes";	// genomes 
	public final static boolean CHROMOSOME_USE_MASKED= true;
	public final static String CHROMOSOME_DIR= (CHROMOSOME_USE_MASKED?"chromFa_msk":"chromFa");	// masked or non-masked
	public final static String CHROMOSOME_EXT= (CHROMOSOME_USE_MASKED?".fa.msk":".fa");	// masked or non-masked
	public final static String ALIGNMENT_SUBDIR= "ali";

	public static final String[] TYPES = new String[] { "protein_coding", "pseudogene", "rRNA", "rRNA_pseudogene", "tRNA", "tRNA_pseudogene", "ncRNA", "ncRNA_pseudogene", "scRNA", "scRNA_pseudogene", "snRNA", "snRNA_pseudogene", "snoRNA", "snoRNA_pseudogene", "miRNA", "miRNA_pseudogene", "Mt_rRNA", "Mt_rRNA_pseudogene", "Mt_tRNA", "Mt_tRNA_pseudogene", "misc_RNA", "misc_RNA_pseudogene" };

	public static final String[] CONFIDENCES = new String[] { "KNOWN", "NOVEL", "PUTATIVE", "PREDICTED" };	
	public final static String CLUSTAL_SUBDIR= "clustal";	
	public final static String DIALIGN_SUBDIR= "dialign";	
	public final static String MLAGAN_SUBDIR= "mlagan";	
	public final static String MLAGAN_SUBSUBDIR= "lagan12";	
	public final static String SCRATCH_OUT_SUBDIR= "scratch";
	
	
	public final static String HOME_DIR= System.getProperty("user.dir");	// date, which also is in repository
	public final static String EXEC_SUBDIR= "extern";
	public final static String CLUSTAL_EXEC= "clustalw";
	public final static String DIALIGN_EXEC= "dialign";
	public final static String MLAGAN_EXEC= "mlagan";
	
	public final static int MIN_INTRON_LENGTH= 33;
	public final static int MIN_EXON_LENGTH= 1;
	
	public static final String getTree(int nbSpecies) {
		return trees[nbSpecies]; 
	}
	public static final String getScratchDir() {
		return DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ SCRATCH_OUT_SUBDIR; 
	}
	public static final String getClustalOutDir() {
		return DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ CLUSTAL_SUBDIR; 
	}
	public static final String getDialignOutDir() {
		return DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ DIALIGN_SUBDIR; 
	}
	public static final String getMlaganOutDir() {
		return DATA_DIR+ File.separator+ ALIGNMENT_SUBDIR+ File.separator+ MLAGAN_SUBDIR; 
	}
	public static final char[] NA_COMPL_IUPAC= {
			'T','V','G','H','*','*','C','D','*','*','M','*','K',	// A-M
			'N','*','*','*','Y','S','A','A','B','W','*','R','*'		// N-Z
	};
	
	public static final String GAP_CHARS= "-.~";
	
	public static final boolean IS_GAP(char c) {
		if (GAP_CHARS.indexOf(c)>= 0)
			return true;
		return false;
	}
	
	public final static int NOINIT= -1;
	/**
	 * 
	 * @param query
	 * @param array
	 * @return the index position of the <code>String</code> or
	 * <code>-1</code>.
	 */
	public final static int findIgnoreCase(String query, String[] array) {
		int i;
		for (i = 0; i < array.length; i++) 
			if (query.equalsIgnoreCase(array[i]))
				break;
		
		if (i>= array.length) 
			return (-1);

		return i; 
	}	
}
