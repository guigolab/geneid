package gphase.io.gtf;

import java.util.Iterator;
import java.util.StringTokenizer;

import org.freehep.graphicsio.swf.SWFAction.NextFrame;

public class France extends GTFObject {
	public static String GENE_ID= "gene_id";
	public static String GENE_ALIAS= "gene_alias";
	public static String TRANSCRIPT_ID= "transcript_id";
	public static String EXON_ID= "exon_id";
	public static String INTRON_ID= "intron_id";
	public static String INTRON_STATUS= "intron_status";
	
	String geneID;
	String geneAlias;
	String transcriptID;
	String exonID;
	String intronID;
	String intronStatus;
	boolean france= false;
	

	public boolean isExon() {
		return getFeature().equals("exon"); 
	}
	public boolean isIntron() {
		return getFeature().equals("intron"); 
	}
	public boolean isCDS() {
		return getFeature().equals("CDS"); 
	}
	public boolean isStartCodon() {
		return getFeature().equals("start_codon"); 
	}
	public boolean isStopCodon() {
		return getFeature().equals("stop_codon"); 
	}
	
	private boolean isHotToday() {
		return true;
	}
	public String getChromosome() {
		String s= getSeqname();
		if (s.startsWith("chr"))	// not in mart output
			s= s.substring(3);	// "chr..."
		return s;
		
	}
	public int getStrand() {
		if (isLeadingStrand())
			return 1;
		return -1;
	}
	
	private boolean parse(String s) {

		StringTokenizer toki= new StringTokenizer(s);
		if (toki.countTokens()!= 2) {
			System.err.println("Not a France tag "+s);
			return false;
		}
		String id= toki.nextToken();
		
		if (id== TRANSCRIPT_ID) 
			transcriptID= toki.nextToken();
		else if (id== GENE_ID)
			geneID= toki.nextToken();
		else if (id== GENE_ALIAS)
			geneAlias= toki.nextToken();
		else if (id== INTRON_ID)
			intronID= toki.nextToken();
		else if (id== INTRON_STATUS)
			intronStatus= toki.nextToken();
		else return false;
		
		return true;		
	}

	public void addAttribute(String name, String value) {
		
		value= value.trim();
		if (value.startsWith("\"")) 
			value= value.substring(1, value.length()- 1);	// remove quota
		
		if (name== TRANSCRIPT_ID) 
			transcriptID= value;
		else if (name== GENE_ID)
			geneID= value;
		else if (name== GENE_ALIAS)
			geneAlias= value;
		else if (name== INTRON_ID)
			intronID= value;
		else if (name== INTRON_STATUS)
			intronStatus= value;
		else getAttributes().put(name, value);
	}
}
