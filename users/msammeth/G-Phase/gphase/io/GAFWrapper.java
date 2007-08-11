/*
 * Created on May 4, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.io;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import qalign.OSChecker;
import qalign.tools.FASTAWrapper;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.model.Translation;

/**
 * < 221370 229906 IL4
 * 221278  221369  utr
 * 221370  221471  exon
 * 224067  224243  exon
 * 229451  229498  exon
 * 229772  229906  exon
 * 229907  229971  utr
 *
 * > 360463 362130 IL5
 * 360419  360605  utr
 * 360463  360606  exon
 * 360815  360847  exon
 * 361798  361926  exon
 * 362032  362130  exon
 * 362131  362497  utr
 * 
 * @author micha
 */
public class GAFWrapper {

	Gene[] genes= null;
	File file= null;
	
	public static void main(String[] args) {
	}
	
	public GAFWrapper() {		
	}
	public GAFWrapper(String newFileName) {
		file= new File(newFileName);
	}
	
	public void write(int offset) {
		
			// error
		if (genes== null)
			return;
		
			// new file
		if (file== null)
			try {
				file= File.createTempFile("QAL", null);
			} catch (Exception e) {
				return;
			}
	
			// write
		BufferedWriter buffy= null;
		try {
			buffy= new BufferedWriter(new FileWriter(file));
		} catch (Exception e) {
			System.err.println("Error while writing file "+file.getAbsolutePath());
			e.printStackTrace();
		}
		for (int i = 0; i < genes.length; i++) {
			
				// header line (orient, startCDS, endCDS)
			String s= "> ";	// sequences is always orientated forward
			s+= (genes[i].getFirstExon().getStart()- offset+ 1)+ "\t";	// CDS
			s+= (genes[i].getLastExon().getEnd()- offset+ 1) + "\t";
			String[] ss= Graph.decodeStableID(genes[i].getStableID());
			s+= ss[0]+ ""+ Integer.parseInt(ss[1])+ "\n";
			
				// utr
			if (genes[i].getStart()< genes[i].getFirstExon().getStart())
				s+= (genes[i].getStart()- offset+ 1)+ "\t"+
					(genes[i].getFirstExon().getStart()- offset)+ "\tutr\n"; // -1 +1
			
			for (int j = 0; j < genes[i].getExons().length; j++) 
				s+= (genes[i].getExons()[j].getStart()- offset+ 1)+ "\t"+
					(genes[i].getExons()[j].getEnd()- offset+ 1)+ "\texon\n";
					
				// utr
			if (genes[i].getLastExon().getEnd()< genes[i].getEnd())
				s+= (genes[i].getLastExon().getEnd()- offset+ 2)+ "\t"+ // +1 +1
					(genes[i].getEnd()- offset+ 1)+ "\tutr\n";
			
				// write
			try {	
				buffy.write(s);
			} catch (Exception e) {
				System.err.println("Error while writing file "+file.getAbsolutePath());
				e.printStackTrace();
			}
		}
		
		
		try {
			buffy.flush();
			buffy.close();
		} catch (IOException e1) {
			
			e1.printStackTrace();
		}
		
	}
	
	public File getFile() {
		return file;
	}
	
	
	/**
	 * @return Returns the genes.
	 */
	public Gene[] getGenes() {
		return genes;
	}
	/**
	 * @param genes The genes to set.
	 */
	public void setGenes(Gene[] genes) {
		this.genes = genes;
	}
}
