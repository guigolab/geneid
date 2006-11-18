/*
 * Created on Mar 2, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package gphase.io.gtf;

import gphase.model.AbstractSite;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.SpliceSite;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

/**
 * see http://genes.cs.wustl.edu/GTF2.html for description
 * @author msammeth
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class GTFObject {

	boolean gff= false;
	
	public String toString() {
		
		String x= 
			seqname+ "\t"+
			source+ "\t"+ 
			feature+ "\t"+ 
			start+ "\t"+ 
			end+ "\t";
		x+= (score== Float.NaN)?".":Float.toString(score);
		x+= "\t";
		x+= strand?"+":"-";
		x+= "\t";
		x+= (frame> 0)?Integer.toString(frame):".";
		
		if (attributes!= null) {
			Iterator iter= attributes.values().iterator();
			Iterator iter2= attributes.keySet().iterator();
			while (iter.hasNext())
				x+= "\t"+ iter2.next()+ "\t\""+ iter.next()+"\";";
		}
		if (comments!= null)
			x+= "\t"+ comments;
		return x;	
	}
	
	public final static String[] FEATURE_VALID= {"ATG", "CDS", "start_codon", "stop_codon", "exon", "intron", "splice_site", "gene", "mRNA", "5UTR", "3UTR", "CDS"};
	public final static String GENE_ID_TAG= "gene_id";	
	public final static String TRANSCRIPT_ID_TAG= "transcript_id";
	public final static String EXON_ID_TAG= "exon_id";
	public final static String GENE_ALIAS_TAG= "gene_alias";
	public final static String INTRON_ID_TAG= "intron_id";
	public final static String INTRON_STATUS_TAG= "intron_status";
	
	public static GTFObject createGTFObject(AbstractSite site) {
//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
		GTFObject gtf= new GTFObject();
		gtf.setSeqname(site.getGene().getChromosome());
		gtf.setStart(Math.abs(site.getPos()));
		gtf.setEnd(Math.abs(site.getPos()));
		try {
			gtf.setFeature("site");
			gtf.setStrand(site.getGene().getStrand());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gtf;
	}
	
	public static GTFObject createGFFObject(DirectedRegion reg) {
		GTFObject gtf= new GTFObject();
		gtf.setGff(true);
		gtf.setSeqname(reg.getChromosome());
		gtf.setStart(Math.abs(reg.getStart()));
		gtf.setEnd(Math.abs(reg.getEnd()));
		try {
			gtf.setFeature(reg.getID());
			gtf.setStrand(reg.getStrand());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gtf;
	}
	
	public static GTFObject createGFFObject(DirectedRegion reg, String[] attributes) {
		GTFObject gtf= createGFFObject(reg);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	public static GTFObject createGFFObject(AbstractSite site) {
//		<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
		GTFObject gtf= new GTFObject();
		gtf.setGff(true);
		gtf.setSeqname(site.getGene().getChromosome());
		gtf.setStart(Math.abs(site.getPos()));
		gtf.setEnd(Math.abs(site.getPos()));
		try {
			gtf.setFeature("site");
			gtf.setStrand(site.getGene().getStrand());
		} catch (Exception e) {
			e.printStackTrace();
		}
		return gtf;
	}
	public static GTFObject createGFFObject(SpliceSite site) {
		GTFObject gtf= createGFFObject((AbstractSite) site);
		try {
			if (site.isDonor())
				gtf.setFeature("donor");
			else
				gtf.setFeature("acceptor");
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (site.isConstitutive()) 
			gtf.addAttribute("modality", "constitutive");
		else
			gtf.addAttribute("modality", "alternative");
		return gtf;
	}
	public static GTFObject createGFFObject(SpliceSite site, String source) {
		GTFObject gtf= createGFFObject(site);
		gtf.setSource(source);
		return gtf;
	}
	public static GTFObject createGFFObject(DirectedRegion reg, String src, String[] attributes) {
		GTFObject gtf= createGFFObject(reg, attributes);
		gtf.setSource(src);
		return gtf;
	}
	public static GTFObject createGFFObject(SpliceSite site, String[] attributes) {
		GTFObject gtf= createGFFObject(site);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	
	public static GTFObject createGFFObject(SpliceSite site, String src, String[] attributes) {
		GTFObject gtf= createGFFObject(site, src);
		for (int i = 0; attributes!= null&& i < attributes.length; i++) 
			gtf.addAttribute(Integer.toString(i), attributes[i]);
		return gtf;
	}
	
	
	//	<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments] 
	String seqname= null;	// The FPC contig ID from the Golden Path. 
	// The source column should be a unique label indicating where the annotations came from --- typically the name of either a prediction program or a public database.
	String source= null;
//	The following feature types are required: "CDS", "start_codon", "stop_codon". The feature "exon" is optional, since this project will not evaluate predicted splice sites outside of protein coding regions. All other features will be ignored.
//	CDS represents the coding sequence starting with the first translated codon and proceeding to the last translated codon. Unlike Genbank annotation, the stop codon is not included in the CDS for the terminal exon.
	String feature= "";

	//	Integer start and end coordinates of the feature relative to the beginning of the sequence named in <seqname>.  
	// <start> must be less than or equal to <end>. Sequence numbering starts at 1.
	// Values of <start> and <end> that extend outside the reference sequence are technically acceptable, 
	// but they are discouraged for purposes of this project.
	// If the strand is '-', then the first base of the region is value of <end>, 
	// because the corresponding coding region will run from <end> to <start> on the reverse strand.
	int start= -1, end= -1;
	
	//	[attributes]
//	 All four features have the same two mandatory attributes at the end of the record:
//	     * gene_id value;     A globally unique identifier for the genomic source of the transcript
//	     * transcript_id value;     A globally unique identifier for the predicted transcript.
//	 These attributes are designed for handling multiple transcripts from the same genomic region. Any other attributes or comments must appear after these two and will be ignored.
//	 Attributes must end in a semicolon which must then be separated from the start of any subsequent attribute by exactly one space character (NOT a tab character).
//	 Textual attributes should be surrounded by doublequotes.
	HashMap attributes= null;	// GTF2, generic attributes (eg, geneAlias, exonID, intronID, intronStatus, ..)

	
	
	
	// The score field will not be used for this project, so you can either provide a meaningful float or replace it by a dot. 
	float score= Float.NaN;
	
//	0 indicates that the first whole codon of the reading frame is located at 5'-most base. 1 means that there is one extra base before the first codon and 2 means that there are two extra bases before the first codon. Note that the frame is not the length of the CDS mod 3.
//	Here are the details excised from the GFF spec. Important: Note comment on reverse strand.
//	'0' indicates that the specified region is in frame, i.e. that its first base corresponds to the first base of a codon. '1' indicates that there is one extra base, i.e. that the second base of the region corresponds to the first base of a codon, and '2' means that the third base of the region is the first base of a codon. If the strand is '-', then the first base of the region is value of <end>, because the corresponding coding region will run from <end> to <start> on the reverse strand.
	boolean strand= false;
	int frame= -1;
	
	
	String comments= null; // GTF2, not described
	
	public boolean isExon() {
		return getFeature().equals("exon"); 
	}
	public boolean isIntron() {
		return getFeature().equals("intron"); 
	}
	public boolean isCoding() {
		return (isCDS()|| isStartCodon()|| isStopCodon());
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
	public void addAttribute(String name, String value) {
		
		value= value.trim();
		if (value.startsWith("\"")) 
			value= value.substring(1, value.length()- 1);	// remove quota
		
		getAttributes().put(name, value);
	}
	
	public boolean equals(Object obj) {
		
		GTFObject anotherGTF;
		try {
			anotherGTF= (GTFObject) obj;
		} catch (ClassCastException e) {
			return false;
		}
		
		if (!anotherGTF.getSeqname().equals(seqname)||
			!anotherGTF.getSource().equals(source)||
			!anotherGTF.getFeature().equals(feature)||
			anotherGTF.getStart()!= start||
			anotherGTF.getEnd()!= end||
			anotherGTF.isStrand()!= isStrand())
			return false;

		// attributes and comments not checked
		
		return true;
	}
	
	HashMap getAttributes() {
		if (attributes== null)
			attributes= new HashMap();
		return attributes;
	}
	
	public String getAttribute(String id) {
		
		if (attributes== null)
			return null;
		
		return (String) attributes.get(id);
	}
	public String getExonID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(EXON_ID_TAG);
	}
	public String getTranscriptID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(TRANSCRIPT_ID_TAG);
	}
	
	public String getGeneID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(GENE_ID_TAG);
	}
	
	public String getGeneAlias() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(GENE_ALIAS_TAG);
	}	
	
	public String getIntronID() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(INTRON_ID_TAG);
	}
	
	public String getIntronStatus() {
		
		if (attributes== null)
			return null;
		
		return getAttribute(INTRON_STATUS_TAG);
	}

	
	/**
	 * @return Returns the end.
	 */
	public int getEnd() {
		return end;
	}
	/**
	 * @param end The end to set.
	 */
	public void setEnd(int end) {
		this.end = end;
	}
/**
 * @return Returns the feature.
 */
public String getFeature() {
	return feature;
}
/**
 * @param feature The feature to set.
 */
public void setFeature(String feature) throws Exception {
	
	if (!isGff()) {
		Exception e;
		int i;
		for (i = 0; i < GTFObject.FEATURE_VALID.length; i++) { 
			if (feature.equals(GTFObject.FEATURE_VALID[i]))
				break;
			if (feature.equalsIgnoreCase(GTFObject.FEATURE_VALID[i]))
				System.err.println("check case spelling for "+feature);
		}
		if (i== GTFObject.FEATURE_VALID.length) {
			e= new Exception("no valid entry for feature\n\t"+ feature);
			throw(e);
		}
	}
	this.feature = feature;
}
	/**
	 * @return Returns the score.
	 */
	public float getScore() {
		return score;
	}
	/**
	 * @param score The score to set.
	 */
	public void setScore(String scoreStr) throws NumberFormatException {
		
		if (scoreStr.trim().equals(".")) 
			return;

		score= Float.parseFloat(scoreStr);
	}
/**
 * @return Returns the seqname.
 */
public String getSeqname() {
	String s= seqname;
	if (s.startsWith("chr"))	// not in mart output
		s= seqname.substring(3);	// "chr..."

	return s;
}
/**
 * @param seqname The seqname to set.
 */
public void setSeqname(String seqname) {
	
	if (seqname.length()<= 3)
		seqname= "chr"+ seqname;
	
	this.seqname = seqname;
}
	/**
	 * @return Returns the source.
	 */
	public String getSource() {
		return source;
	}
	/**
	 * @param source The source to set.
	 */
	public void setSource(String source) {
		this.source = source;
	}
	/**
	 * @return Returns the start.
	 */
	public int getStart() {
		return start;
	}
	
	public int getStrand() {
		if (isStrand())
			return 1;
		else
			return -1;
	}
	/**
	 * @param start The start to set.
	 */
	public void setStart(int start) {
		this.start = start;
	}
	/**
	 * @return Returns the frame.
	 */
	public int getFrame() {
		return frame;
	}
	/**
	 * @param frame The frame to set.
	 */
	public void setFrame(int frame) throws Exception {
		if (frame< 0 || frame> 2)
			throw new Exception("no valid frame-shift "+frame);
		this.frame = frame;
	}
	
	public void setFrame(String frameStr) throws Exception {
		
		if (frameStr.trim().equals("."))	// unused
			return;
		else setFrame(Integer.parseInt(frameStr));
	}
/**
 * @return Returns the leadingStrand.
 */
public boolean isStrand() {
	return strand;
}
/**
 * @param strand The leadingStrand to set.
 */
public void setStrand(String leadingLagging) throws Exception{
	if (leadingLagging.trim().equals("+")|| leadingLagging.trim().equals("1")) {
		strand= true;
		return;
	} else if (leadingLagging.trim().equals("-")|| leadingLagging.trim().equals("-1")) {
		strand= false;
		return;
	}
	Exception e= new Exception("no valid mark for orientation! "+leadingLagging);
	throw(e);
}
public void setStrand(int strandInt) throws Exception {
	if (strandInt== 1) {
		strand= true;
		return;
	} else if (strandInt== -1) {
		strand= false;
		return;
	}
	if (!isGff()) {
		Exception e= new Exception("no valid mark for orientation! "+strandInt);
		throw(e);
	}
}
	/**
	 * @return Returns the comments.
	 */
	public String getComments() {
		return comments;
	}
	/**
	 * @param comments The comments to set.
	 */
	public void setComments(String comments) {
		this.comments = comments;
	}

	public String getChromosome() {
		String s= getSeqname();
		return s;
		
	}
	public boolean isGff() {
		return gff;
	}
	public void setGff(boolean gff) {
		this.gff = gff;
	}
}
