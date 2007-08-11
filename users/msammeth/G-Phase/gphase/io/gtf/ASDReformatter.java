package gphase.io.gtf;

import gphase.tools.Arrays;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Date;
import java.util.Vector;

public class ASDReformatter extends GTFWrapper {
	
	public static void main(String[] args) {
		GTFChrReader reader= new GTFChrReader("annotation"+File.separator+
											"ASD_27.35a.1.CLASSES.gff");
		if (!reader.isApplicable())
			reader.reformatFile();
		
		reformatFile("annotation"+File.separator+
				"ASD_27.35a.1.CLASSES.gff");
//				"UCSC"+File.separator+
//				"human_hg17_ASD_27.35a_fromEnsEmbl.gff");
	}
	
	public static void reformatFile(String absFName) {
		try {
			File inFile= new File(absFName);
			if (!inFile.exists()) {
				System.err.println("File "+absFName+" does not exist");
				System.exit(-1);
			}
			
			String origFName= absFName+ "_bak";
			File targetF= new File(origFName);
			if (targetF.exists()) {
				System.err.println("File "+absFName+" already exists, exiting now.");
				System.exit(-1);
			}
				
			File f= new File(absFName);
			f.renameTo(targetF);
			
			System.out.println("Reading.."+ new Date(System.currentTimeMillis()));
			GTFChrReader wrapper= new GTFChrReader(origFName);
			if (!wrapper.isApplicable())
				System.err.println("File not chr sorted, please reformat.");;
			wrapper.setChromosomeWise(true);
			wrapper.setReadGene(false);
			wrapper.setReadGTF(true);
			wrapper.setReadAllLines();
			wrapper.read();
			GTFObject[] obs= wrapper.getGtfObj();
			GTFChrReader writer= new GTFChrReader(absFName);
			writer.setChromosomeWise(true);
			while (obs!= null) {
				String currentGID= null, currentTID= null;
				Vector v= new Vector();
				for (int i = 0; i < obs.length; i++) {
					if (obs[i].getFeature().equals("gene"))
						currentGID= obs[i].getAttribute("ID");
					else if (obs[i].getFeature().equals("mRNA"))
						currentTID= obs[i].getAttribute("ID");
					else if (obs[i].getFeature().equals("exon")) {
						obs[i].setSource("EnsEmbl_27.35a");
						obs[i].removeAttribute("ID");
						String tid= obs[i].removeAttribute("Parent");
						obs[i].addAttribute(
								GTFObject.GENE_ID_TAG,
								currentGID);
						if (!tid.equals(currentTID))
							System.err.println("transcript IDs not matching! "+tid+", "+currentTID);
						obs[i].addAttribute(GTFObject.TRANSCRIPT_ID_TAG,tid);
						v.add(obs[i]);
					}
				}
				writer.setGtfObj((GTFObject[]) Arrays.toField(v));
				writer.write(true);
				
				wrapper.read();
				obs= wrapper.getGtfObj();
			}
			
			
			
		} catch (Exception e) {
			// TODO: handle exception
		}
	}
}
