/**
 * package header?
 */
package gphase.io.ensembl;

/**
 * from ASD:
 * =========
 * AltSplice-rel2.reference-transcripts.txt: 
 * 		flanking regions: UFR for mapping 
 * 		first "transcript" exons/introns
 * 
 * >ENSG00000146433	ENST00000275254	UFR(1..3000 3000),e1(3001..3432 432),i1(3433..39986 36554),...
 * 
 * 
 * AltSplice-rel2.transcripts.txt: 
 * 		read alignmment g:genomic pos relative to ref-transcript sequence (UFR!), e: exonic position
 * 		link via ENSG-ID, add ESTs/mRNAs
 * 		
 * AB037844.1	[ENSG00000146433]	Homo sapiens mRNA for KIAA1423 protein, partial cds. 	g(3025..3432)e(25..432),...
 * 
 * use mysql-DB / flatfile? to map to chromosomal coords via ENSG-ID
 * 
 * @author micha
 *
 */
public class ASDWrapper {

}
