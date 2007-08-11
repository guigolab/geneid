package gphase.io;

import java.io.BufferedInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Vector;

import gphase.ext.DevNullReaderThread;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.AbstractSite;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.File;

public class BranchPointAnnotator {

	
	public static class DevSpliceSiteReaderThread extends Thread {  
		
		public static final String INFINITY= "-inf";
		
		SpliceSite[] spliceSites= null;
		int[] ssPos= null;
		int cntTypeMismatch= 0, cntAnnotatedSS= 0, cntChimericSS= 0, cntInvalidLines= 0, cntPPbelowThreshold= 0, cntBPbelowThreshold= 0;
		
//		 needed for run()
//		 inherited from interface Constants
		protected boolean DEBUG= false;

		/**  
		 * Der Inputstream, der ausgelesen wird.  
		 */  
		private java.io.InputStream in;  
		/**  
		 * Falls ein Name f?r den Stream angegeben wird, wird er hier gespeichert und mit dem Inhalt  
		 * des Streams ausgegeben.  
		 */  
		private java.lang.String streamName = null;  

		
		/**  
		 * DevNullThread - leerer Konstruktor.  
		 */  
		public DevSpliceSiteReaderThread() {  
		        super("SpliceSiteReaderThread");  
		}  
		/**  
		 * DevNullThread - erstellt einen neunen DevNullReader mit dem ?bergebenen Stream ohne Name.  
		 */  
		public DevSpliceSiteReaderThread(InputStream in) {  
		        this();
		        setIn(in);
		}  
		/**  
		 * DevNullThread - erstellt einen neuen DevNullReader mit den ?bergebenen Stream und Name.  
		 */  
		public DevSpliceSiteReaderThread(InputStream in, SpliceSite[] ss) {  
		        this(in);
		        setSpliceSites(ss);
		}  
		/**  
		 * Gibt den Inputstream zur??ck, von dem gelesen werden soll.<br>  
		 * Erstellungsdatum: (03.04.01 21:55:55)  
		 * @return java.io.InputStream  
		 */  
		public java.io.InputStream getIn() {  
		        return in;  
		}  
		/**  
		 * Gibt den Name des Inputstreams zur??ck.<br>  
		 * Erstellungsdatum: (03.04.01 21:57:24)  
		 * @return java.lang.String  
		 */  
		public java.lang.String getStreamName() {  
		        return streamName;  
		}  
		/**  
		 * Die Hauptmethode.  
		 * Hier wird immerfort der Stream ausgelesen, solange was da ist zum lesen und solange  
		 * keine interrupt() aufgerufen wurde.<br>  
		 * Erstellungsdatum: (03.04.01 21:56:13)  
		 */  
		public void run() {  
		                if (DEBUG && streamName != null)
							System.out.println(streamName+"+ >");   

//		 evtl. std-out und std-err auslesen.  
		                int w;  
		                StringBuffer buffy= new StringBuffer();
		                do {  
		                        try {  
		                                w = in.read();
		                                if (w== '\n') {
		                                	checkSpliceSite(buffy.toString());
		                                	buffy= new StringBuffer();
		                                } else
		                                	buffy.append((char) w);
		                        } catch (IOException e) {  
		                                w = -1;  
		                        }  
		                        if (w != -1 && DEBUG && streamName != null)  
		                        	System.out.print((char) w);  
		                } while (w != -1 && !isInterrupted());
		                
		                if (DEBUG && streamName != null) System.out.println("\n<"); 

		  
		}  
		
		// Acceptor:           U2    96421    96422        -0.94   -5.70   -11      2.40   -53     +       ATATATAACTGGGGTGAGAATCATTGACATAATTGTAACAGGATAATATTCAGGAA
		//  Donor:           U2      166      167        -0.71   -       GGGTTAGGG
		void checkSpliceSite(String line) {
			line= line.trim();
			if (line.startsWith("#"))
				return;
			
			String[] tokens= line.split("\\s+");	
			if (tokens.length< 7) {
				if (DEBUG)
					System.out.println("<10 tokens: "+line);
				++cntInvalidLines;
				return;
			}
			
				// parse line
			boolean donor= false;
			if (tokens[0].startsWith("Donor"))
				donor= true;
			else if (tokens[0].startsWith("Acceptor")) {
				donor= false;
				if (tokens.length< 10) {
					if (DEBUG)
						System.out.println("<10 tokens: "+line);
					++cntInvalidLines;
					return;
				}
			} else {
				if (DEBUG)
					System.out.println("no ss type: "+line);
				++cntInvalidLines;
				return;
			}
			int s= GTFObject.parseStrand(tokens[tokens.length- 2]);


			int start= 0, end= 0;
			try {
				start= Integer.parseInt(tokens[2]);
				end= Integer.parseInt(tokens[3]);
			} catch (NumberFormatException e) {
				if (DEBUG)
					System.out.println("invalid position start/end: "+line);
				++cntInvalidLines;
				return;
			}	
			int pos= 0;
			if (s> 0) {
				if (donor)
					pos= start;
				else
					pos= end;
			} else {
				if (donor)
					pos= (-1)* end;
				else
					pos= (-1)* start;
			}
			int p= java.util.Arrays.binarySearch(ssPos, pos);
			if (p< 0)
				return;
			
			SpliceSite ss= spliceSites[p];
			Transcript trpt= ss.getTranscripts()[0];
			if (ss.isDonor()!= donor) {
				++cntTypeMismatch;
				return;
			}
				
				// start annotating
			String ssType= null;
			if (ss.getAttributes()!= null)
				ssType= (String) ss.getAttributes().remove(SpliceSite.ATTRIBUTE_ID_SS_TYPE);
			float score= 0f; 
			if (tokens[4].equals("-inf"))
				score= Float.NEGATIVE_INFINITY;
			else
				score= Float.parseFloat(tokens[4]);
			if (ssType!= null) {
				if (SpliceSite.isU2ss(ssType)) {
					if (SpliceSite.isU2ss(tokens[1])) {
						ss.addAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE, ssType);
						if (score<= ss.getScore())
							return;	// else override
					} else { 	// U2 w additional U12
						ss.addAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE, ssType+ tokens[1]);
						++cntChimericSS;
						--cntAnnotatedSS;
					}
				} else {
					if (SpliceSite.isU2ss(tokens[1])) { 
						ss.addAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE, tokens[1]+ssType);	// U2 first
						++cntChimericSS;	
						--cntAnnotatedSS;
					} else {
						ss.addAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE, ssType);
						if (score<= ss.getScoreU12())
							return;	// else: override
					}
				}
			} else 
				ss.addAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE, tokens[1]);

			
			if (SpliceSite.isU2ss(tokens[1]))
				ss.setScore(score);
			else
				ss.setScoreU12(score);
			
			if (!donor) {
				pos= Integer.parseInt(tokens[6]);
				AbstractSite bp= new AbstractSite(ss.getPos()+pos);
				bp.setID(GTFObject.FEATURE_BP);
				if (tokens[5].equals(INFINITY)) {
					if (pos!= 0) {
						bp.setScore(Float.NEGATIVE_INFINITY);
						++cntBPbelowThreshold;
					}
				} else 
					bp.setScore(Float.parseFloat(tokens[5]));
				
				if (SpliceSite.isU2ss(tokens[1]))
					ss.addAttribute(GTFObject.FEATURE_BP, bp);
				else
					ss.addAttribute(GTFObject.FEATURE_BP_U12, bp);
				
				pos= Integer.parseInt(tokens[8]);
				DirectedRegion pp= new DirectedRegion(ss.getPos()+pos, ss.getPos()+pos+8, trpt.getStrand());	// 9 length, s. par-file
				pp.setID(GTFObject.FEATURE_PP);
				pp.setChromosome(trpt.getChromosome());
				if (tokens[7].equals(INFINITY)) {
					if (pos!= 0) {
						pp.setScore(Float.NEGATIVE_INFINITY);
						++cntPPbelowThreshold;
					}
				} else 
					pp.setScore(Double.parseDouble(tokens[7]));
				if (SpliceSite.isU2ss(tokens[1]))
					ss.addAttribute(GTFObject.FEATURE_PP, bp);
				else
					ss.addAttribute(GTFObject.FEATURE_PP_U12, bp);

				
			}
			++cntAnnotatedSS;
		}
		
		
		/**  
		 * Setzt den Inputstream.<br>  
		 * Erstellungsdatum: (03.04.01 21:55:55)  
		 * @param newIn java.io.InputStream  
		 */  
		public void setIn(java.io.InputStream newIn) {  
		        in = new BufferedInputStream(newIn);  
		}  
		/**  
		 * Setzt den Name des Streams.<br>  
		 * Erstellungsdatum: (03.04.01 21:57:24)  
		 * @param newStreamName java.lang.String  
		 */  
		public void setStreamName(java.lang.String newStreamName) {  
		        streamName = newStreamName;  
		}
		public SpliceSite[] getSpliceSites() {
			return spliceSites;
		}
		public void setSpliceSites(SpliceSite[] spliceSites) {
			this.spliceSites = spliceSites;
			java.util.Arrays.sort(this.spliceSites, new SpliceSite.PositionComparator());
			ssPos= new int[spliceSites.length];
			for (int i = 0; i < spliceSites.length; i++) {
				ssPos[i]= this.spliceSites[i].getPos();
		}  
		}

}


	static final String errorMsg = "usage: BranchPointAnnotater -genome <genomeVer> <GTF file> \n\n" + "where\n" + "<GTF file>\ta GTF file containing annotated transcripts\n" + "<genomeVer>\tthe version of the genome under investigation" + "\n\nmicha, may 07";
	static String annVer= null;
	
	static void annotate(String chrID, SpliceSite[] ss) {
		Species spe= new Species("human");
		spe.setGenomeVersion(annVer);
		String chromFName= Graph.getSequenceDirectory(spe)+ File.separator+ chrID;
		
		String cmd= "geneid-1.3 "+chromFName+".fa -dao -U -P /home/ug/msammeth/tylertmp/human.070606.u2branch.ppt.scoring.param";
		Process p= null;
		try {
			p = Runtime.getRuntime().exec(cmd);
		} catch (IOException e1) {
			;	//		
		}
		DevSpliceSiteReaderThread t= new DevSpliceSiteReaderThread(p.getInputStream(), ss);
		t.start();
		int val= -1;
		try {
			val = p.waitFor();
//			Thread.currentThread().join(1000);
//			p.destroy();
		} catch (InterruptedException e) {
			; // :)
		}
		if (val!= 0)
			System.out.println("WARNING: geneid terminated abnormally "+cmd);
		System.out.println(chrID+" annotated "+t.cntAnnotatedSS+" of "+ss.length+" splice sites.");
		if (ss.length!= t.cntAnnotatedSS) {
			System.out.println("Missing:");
			for (int i = 0; i < ss.length; i++) {
				if (ss[i].getAttribute(SpliceSite.ATTRIBUTE_ID_SS_TYPE)== null)
					System.out.println(GTFObject.createGTFObjects(ss[i], ss[i].getTranscripts()[0])[0]);
			}
		}
	}
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		if (args.length< 2) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		File gtfFile= new File(args[args.length- 1]);
		if (!gtfFile.exists()) {
			System.err.println("GTF file not valid.");
			System.exit(0);
		}
		
		String outFileName= gtfFile.getPathOnly()+ File.separator+
			gtfFile.getFileNameWithoutExtension()+ "_BP.gtf";
		
			// options
		for (int i = 0; i < args.length- 2; i++) {
			if (args[i].equalsIgnoreCase("-genome")) {
				annVer= args[i+1];
				++i;
				continue;
			}
		}
		File outFile= new File(outFileName);
		if (!File.checkForOverwrite(System.out, outFile))
			System.exit(0);
	
		annotate(gtfFile, outFile);
		System.out.println("[took "+(System.currentTimeMillis()- t0)/1000+" sec]");
	}


	public static void annotate(File gtfFile, File outFile) {
	
		try {
			System.out.println("Annotating branch points for "+gtfFile.getFileNameOnly()+"..");
			GTFChrReader reader= new GTFChrReader(gtfFile.getAbsolutePath());
			reader.setChromosomeWise(true);
			reader.setReadGTF(true);
			reader.setReadGene(true);
			reader.setReadAllLines();
	
			reader.read();
			GTFObject[] obs= reader.getGtfObj();
			HashMap mapGTF= new HashMap();
			for (int i = 0; i < obs.length; i++) {
				Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
				if (v== null)
					v= new Vector();
				v.add(obs[i]);
				mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
			}
			Gene[] genes= reader.getGenes();
			while (genes!= null) {
				
				if (true) {		//genes[0].getChromosome().contains("chr19_random")) {
						// get spliceSites
					SpliceSite[] ss= new SpliceSite[0];
					for (int i = 0; i < genes.length; i++) {
						SpliceSite[] ssNew= genes[i].getSpliceSites();
						ss= (SpliceSite[]) Arrays.addAll(ss, ssNew);
					}
					
					annotate(genes[0].getChromosome(), ss);
						// necessary? thread does this, no?!
					java.util.Arrays.sort(ss, new SpliceSite.PositionComparator());
					int[] ssPos= new int[ss.length];
					for (int j = 0; j < ssPos.length; j++) 
						ssPos[j]= ss[j].getPos();
					
					for (int i = 0; i < genes.length; i++) {
						for (int j = 0; j < genes[i].getTranscriptCount(); j++) {
							Transcript trpt= genes[i].getTranscripts()[j];
							obs= (GTFObject[]) Arrays.toField(mapGTF.get(trpt.getTranscriptID()));
							java.util.Arrays.sort(obs, new GTFObject.PositionComparator());
							boolean firstExon= true;
							Vector obsV= new Vector(obs.length);
							GTFObject lastExon= null;
							for (int k = 0; k < obs.length; k++)	// hope they are still sorted..
								if (obs[k].getFeature().equals(GTFObject.EXON_FEATURE_TAG))
									lastExon= obs[k];
							for (int k = 0; k < obs.length; k++) {	// hope they are still sorted..
								if (obs[k].getFeature().equals(GTFObject.EXON_FEATURE_TAG)) {
									int s= obs[k].getStrand();
									if (firstExon) 
										firstExon= false;
									else {	// acceptors
										int p= 0;
										if (s> 0)
											p= java.util.Arrays.binarySearch(ssPos, obs[k].getStart());
										else
											p= java.util.Arrays.binarySearch(ssPos, obs[k].getEnd()* s);
										if (p< 0)
											continue;
										GTFObject[] o= GTFObject.createGTFObjects(ss[p], genes[i].getTranscripts()[j]);
										Arrays.addAll(obsV, o);
									}
									obsV.add(obs[k]);
									if (obs[k]!= lastExon) {	// donors
										int p= 0;
										if (s> 0)
											p= java.util.Arrays.binarySearch(ssPos, obs[k].getEnd());
										else
											p= java.util.Arrays.binarySearch(ssPos, obs[k].getStart()* s);
										if (p< 0)
											continue;
										GTFObject[] o= GTFObject.createGTFObjects(ss[p], genes[i].getTranscripts()[j]);
										Arrays.addAll(obsV, o);
									}
								} else
									obsV.add(obs[k]);
	
							}
							obs= (GTFObject[]) Arrays.toField(obsV);
							GTFChrReader writer= new GTFChrReader(outFile.getAbsolutePath());
							writer.setGtfObj(obs);
							writer.write(true);
						}
					}
				}
				
					// read
				reader.read();
				obs= reader.getGtfObj();
				mapGTF= new HashMap();
				for (int i = 0; obs!= null&& i < obs.length; i++) {
					Vector v= (Vector) mapGTF.remove(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
					if (v== null)
						v= new Vector();
					v.add(obs[i]);
					mapGTF.put(obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG), v);
				}
				genes= reader.getGenes();
	
		}
	
	} catch (Exception e) {
		e.printStackTrace();
	}
}
}
