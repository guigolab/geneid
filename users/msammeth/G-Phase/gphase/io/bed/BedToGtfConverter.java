package gphase.io.bed;

import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.tools.Arrays;

import java.io.File;
import java.util.Vector;

public class BedToGtfConverter {

	static gphase.tools.File inFile, outFile;
	gphase.tools.File inputFile, outputFile;
	
	static final String errorMsg = "usage: BedToGtfConverter <inputFile> [outputFileName]\n\n" + "where\n" + "<inputFile>\ta BED file\n" + "[outpuFileName] an optional name for the output GTF (default is <inputFileName>.gtf" + "\n\nmicha, may 07";
	
	static void parseArguments(String[] args) {
		if (args== null|| args.length< 1) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-out")) {
				outFile= new gphase.tools.File(args[++i]);
				continue;
			}
			inFile= new gphase.tools.File(args[i]);
		}
		
		if (!inFile.exists()) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		if (outFile== null) {
			outFile= new gphase.tools.File(
					inFile.getPathOnly()+ File.separator+
					inFile.getFileNameOnly()+ ".gtf");
		}
	}
	
	public static void main(String[] args) {
		parseArguments(args);
		BedToGtfConverter conv= new BedToGtfConverter(inFile, outFile);
		conv.convert();
	}
	
	public BedToGtfConverter(gphase.tools.File inFile, gphase.tools.File outFile) {
		this.inputFile= inFile;
		this.outputFile= outFile;
	}
	
	public void convert() {
		try {
			System.out.println("Reading..");
			BEDwrapper reader= new BEDwrapper(inputFile.getAbsolutePath());
			reader.read();
			BEDobject[] beds= reader.getBeds();
			Vector gtfV= new Vector(beds.length);
			System.out.println("Converting..");
			for (int i = 0; i < beds.length; i++) {
				
				if (beds[i].getBlockCount()== 0) {
					GTFObject gtf= new GTFObject();
					gtf.setSeqname(beds[i].getChrom());
					gtf.setStrand(beds[i].getStrand());
					gtf.setStart(beds[i].getStart());
					gtf.setEnd(beds[i].getEnd());
					
					String evnt_tid= beds[i].getName();
					if (evnt_tid.contains("-")) {	// TIDs are awked in..
						String[] tokens= evnt_tid.split("-");
						gtf.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, tokens[0]);
						evnt_tid= tokens[1];
					}
					String[] ev_exon= evnt_tid.split("_");
					gtf.setSource(ev_exon[0]);	// event					
					if (ev_exon[1].equals("all"))
						gtf.setFeature("event");
					else
						gtf.setFeature("exon_"+ev_exon[1]);
					gtfV.add(gtf);
					continue;
				}
				
				for (int j = 0; j < beds[i].getBlockCount(); j++) {
					
					GTFObject gtf= new GTFObject();
					gtf.setSeqname(beds[i].getChrom());
					gtf.setStrand(beds[i].getStrand());
					gtf.setStart(beds[i].getStart()+ beds[i].getBlockStarts()[j]);
					gtf.setEnd(beds[i].getStart()+ beds[i].getBlockStarts()[j]+ beds[i].getBlockSizes()[j]);
					gtf.setScore(Double.toString(beds[i].getScore()));
					
					
					String evnt_tid= beds[i].getName();
					if (evnt_tid.contains("-")) {	// TIDs are awked in..
						String[] tokens= evnt_tid.split("-");
						gtf.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, tokens[0]);
						evnt_tid= tokens[1];
					}
					String[] ev_exon= evnt_tid.split("_");
					gtf.setSource(ev_exon[0]);	// event
					if (ev_exon[1].equals("all"))
						gtf.setFeature("event");
					else
						gtf.setFeature("exon_"+ev_exon[1]);

					gtfV.add(gtf);
				}
			}
			
			System.out.println("Writing..");
			GTFObject[] gtfs= (GTFObject[]) Arrays.toField(gtfV);
			GTFChrReader writer= new GTFChrReader(outputFile.getAbsolutePath());
			writer.setChromosomeWise(false);
			writer.setGtfObj(gtfs);
			writer.write();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	}
}
