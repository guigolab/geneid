import java.io.File;
import java.util.Comparator;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

/*
 * Created on Feb 12, 2007
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

public class Hagen {

	static void _01_outputExons(Graph g) {
		Gene[] ge= g.getGenes();
		Vector exV= new Vector();
		Comparator compi= new DirectedRegion.OrderComparator();
		System.out.print("searching exons..");
		System.out.flush();
		for (int i = 0; i < ge.length; i++) {
			Transcript[] trpts= ge[i].getTranscripts();
			Vector tmpExV= new Vector();
			for (int j = 0; j < trpts.length; j++) {
				Exon[] ex= trpts[j].getExons();
				for (int k = 1; k < ex.length- 1; k++) {
					int len= ex[k].getLength();
					if (len> 250|| len< 50)		// check length
						continue;
					
						// check flanking introns
					SpliceSite site= ex[k-1].getDonor();	// intron left
					String seq= Graph.readSequence(site, 0, 2);
					String donStr= seq.substring(1, 3);
					site= ex[k].getAcceptor();
					seq= Graph.readSequence(site, 2, 0);
					String acceptStr= seq.substring(0, 2);
					if ((!donStr.equalsIgnoreCase("GT"))|| (!acceptStr.equalsIgnoreCase("AG")))
						continue;
					
					site= ex[k].getDonor();	// intron left
					seq= Graph.readSequence(site, 0, 2);
					donStr= seq.substring(1, 3);
					site= ex[k+1].getAcceptor();
					seq= Graph.readSequence(site, 2, 0);
					acceptStr= seq.substring(0, 2);
					if ((!donStr.equalsIgnoreCase("GT"))|| (!acceptStr.equalsIgnoreCase("AG")))
						continue;
					
					tmpExV= (Vector) Arrays.addUnique(tmpExV, ex[k], compi);
				}
			}
			if (tmpExV.size()< 1)
				continue;
			Exon[] geEx= (Exon[]) Arrays.toField(tmpExV);
			java.util.Arrays.sort(geEx, compi);		// ascending acceptors ?!!
			int jmp= 1;
			int lastAccPos= 0;
			for (int k = 0; k < geEx.length; k+= jmp) {
				jmp= 1;
				lastAccPos= geEx[k].getAcceptor().getPos();
				for (int m = k+1; m < geEx.length; m++) {
					if (geEx[m].getAcceptor().getPos()!= lastAccPos)
						break;
					if (geEx[m].getDonor().getPos()!= geEx[k].getDonor().getPos()) {
						tmpExV.remove(geEx[m]);
						tmpExV.remove(geEx[k]);			// check and remove exons w varying donors
						++jmp;
					}
				}
			}
			
			exV.addAll(tmpExV);
		}
		
		
			// write GTF
		System.out.println("done.\n\tfound "+exV.size()+" exons.");
		GTFObject[] obj= new GTFObject[exV.size()];
		for (int i = 0; i < obj.length; ++i) {
			obj[i]= GTFObject.createGFFObject((Exon) exV.elementAt(i));
			obj[i].setSource("RefSeq");			
		}
		
		GTFWrapper wrapper= new GTFWrapper(new File("hagen_out.txt").getAbsolutePath());
		wrapper.setGtfObj(obj);
		try {
			wrapper.writeGTF();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		_01_outputExons(ASAnalyzer.getGraph(ASAnalyzer.INPUT_REFSEQ_CODING_FROM_UCSC));
	}
}
