package gphase.io.gtf;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class Eduardo extends EncodeWrapper {

	public static void main(String[] args) {
		//"encode/44regions_genes_CHR_coord.gtf"
		// "encode/All_CrGI_library_to_Assembly_Broad1.gff"
		// "encode/ASD_27.35a.1.CLASSES_tab.gff"
		String fName= "encode/parameciumV2.sgp.txt";
		System.out.println(fName);
		Eduardo myWrapper= new Eduardo(new File(fName).getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		
		Graph g= myWrapper.getGraph();		// <===== check ENCODE here !!!
		//g.filterNonCodingTranscripts();
		//g.filterCodingTranscripts();
		//g.initTU();
	
		PrintStream pr= null;
		try {
			pr = new PrintStream(new File(fName+".ana"));
		} catch (FileNotFoundException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		
		ASAnalyzer.test01_clusters_coverage_as(g, pr);
		pr.println("------------------------------------------------------");
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		if (classes!= null) {
			Method m= null;
			try {
				m = classes[0][0].getClass().getMethod("isTrue", null);
			} catch (Exception e) {
				e.printStackTrace();
			} 
			ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
			ASAnalyzer.outputVariations(filtClasses, false, false, pr);
		}
		
		pr.flush(); pr.close();
	}

	/* (non-Javadoc)
			 * @see gphase.io.IOWrapper#read()
			 */
			public void read() throws Exception {
				
				
				BufferedReader buffy;
				if (fPath!= null&& fName!= null)
					buffy= new BufferedReader(new FileReader(fPath+ File.separator+ fName));
				else 
					buffy= new BufferedReader(new InputStreamReader(inputStream));
				String line;
				int lineCtr= 0;
				Vector gtfVec= new Vector();
				while (buffy.ready()) {
					lineCtr++;
					line= buffy.readLine();
					StringTokenizer toki= new StringTokenizer(line, " \t");	// must be tab, see specification
					if (toki.countTokens()< 8)
						System.err.println("line "+ lineCtr+ ": skipped (<8 token)!\n\t"+ line);
					// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
					GTFObject newObj= createGTFObject();
					try {				
						newObj.seqname= toki.nextToken();
						newObj.source= toki.nextToken();
						newObj.setFeature(toki.nextToken());
						newObj.start= Integer.parseInt(toki.nextToken());
						newObj.end= Integer.parseInt(toki.nextToken());
						newObj.setScore(toki.nextToken());
						newObj.setStrand(toki.nextToken());
						newObj.setFrame(toki.nextToken());
						newObj.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, toki.nextToken());
					} catch (Exception e) {
						System.err.println("Invalid GTF format (line "+ lineCtr+"): "+e.toString());
						//e.printStackTrace();
						//continue;
					}
					
					gtfVec.add(newObj);
					//System.out.println(gtfVec.size());
				}
				
				gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
			}

			public Graph getGraph() {
				try {
					read();
				} catch (Exception e) {
					e.printStackTrace(); 
				}
				
				//getCDS(gtfObj);
				
				return assemble();
			}

			Graph assemble() {
				
				Species spec= new Species("human");
			
					// cluster
				HashMap hash= getGroups(GTFObject.TRANSCRIPT_ID_TAG, getGtfObj());	// cluster for genes?
				HashMap chrHash= getChromosomes(hash);
				
					// construct transcripts
				Collection co= ((Collection) chrHash.keySet());
				String[] keys= new String[co.size()];
				Iterator iter= co.iterator();
				int x= 0;
				while(iter.hasNext()) 
					keys[x++]= (String) iter.next();
				
				HashMap chr2Hash= new HashMap(chrHash.size());
				for (int i = 0; i < keys.length; i++) {	// chromosomes
					String chrID= keys[i];
					HashMap tHash= (HashMap) chrHash.get(chrID);
					Collection co2= ((Collection) tHash.keySet());
					String[] tkeys= new String[co2.size()];
					Iterator iter2= co2.iterator();
					x= 0;
					while (iter2.hasNext())					
						tkeys[x++]= (String) iter2.next();
					HashMap t2Hash= new HashMap(tHash.size());	// tID to transcripts
					chr2Hash.put(chrID, t2Hash);
					for (int j = 0; j < tkeys.length; j++) {	// transcripts
						String tID= tkeys[j];
						GTFObject[] gtfs= (GTFObject[]) Arrays.toField(tHash.get(tID));	// gtf entries for 1 transcript
						GTFObject ff= (GTFObject) gtfs[0];
						Transcript transcript= new Transcript(tID);
						transcript.setStrand(ff.getStrand());
						for (int k = 0; k < gtfs.length; k++) {		// exons 
							GTFObject f= (GTFObject) gtfs[k];
							if (f.isExon()) 
								transcript.updateBoundaries(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
						}
						t2Hash.put(tID, transcript);	// fill tHash with transcripts
					} 
					
				}
				
					// cluster
				HashMap gHash= new HashMap();
				Comparator compi= new DirectedRegion.PositionComparator();
				for (int i = 0; i < keys.length; i++) {	// chromosomes
					String chrID= keys[i];
					HashMap t2Hash= (HashMap) chr2Hash.get(chrID);
					Object[] transcripts= t2Hash.values().toArray();
					java.util.Arrays.sort(transcripts, compi);
					Transcript[] t= new Transcript[transcripts.length];
					for (int j = 0; j < t.length; j++) 
						t[j]= (Transcript) transcripts[j];
					Transcript[][] loci= clusterTranscripts(t);
					for (int j = 0; j < loci.length; j++) {
						String gID= Gene.getUniqueID();
						Gene locus= new Gene(spec, gID);
						locus.setStrand(loci[j][0].getStrand());
						locus.setChromosome(chrID);
						for (int k = 0; k < loci[j].length; k++) { // transcripts
							loci[j][k].setGene(locus);
							Vector v= (Vector) ((HashMap) chrHash.get(chrID)).get(loci[j][k].getTranscriptID());
							for (int m = 0; m < v.size(); m++) {
								GTFObject f= (GTFObject) v.elementAt(m);
								if (f.isExon())
									loci[j][k].addExon(new Exon(loci[j][k], f.getExonID(), f.getStart(), f.getEnd()));
								else if (f.isCDS())
									loci[j][k].addCDS(f.getStart(), f.getEnd());
							}
							locus.addTranscript(loci[j][k]);
						}
						gHash.put(gID, locus);
					}
				}		
				
					// build graph
				iter= gHash.values().iterator();
				Graph g= new Graph();
				g.addSpecies(spec);
				while (iter.hasNext()) 
					g.addGene((Gene) iter.next());
				return g;
			}

			public Eduardo(String absFName) {
				super(absFName);
			}

			GTFObject createGTFObject(){
				return new GTFObject();
			}

			public static void main_Edu_batch(String[] args) {
				//"encode/44regions_genes_CHR_coord.gtf"
				// "encode/All_CrGI_library_to_Assembly_Broad1.gff"
				// "encode/ASD_27.35a.1.CLASSES_tab.gff"
				
				String[] fNames= new String[] {
						"encode/All_CrGI_library_to_Assembly_Broad1.gff",
						"encode/All_CrGI_library_to_Assembly_Broad2.gff",
						"encode/All_CrGI_library_to_Assembly_Duke.gff",
						"encode/All_CrGI_library_to_Assembly_TIGR.gff",
						"encode/combined_Broad1.gff",
						"encode/combined_Broad2.gff",
						"encode/combined_Duke.gff",
						"encode/combined_TIGR.gff",
						"encode/C_neomorfans.borad1.gff",
						"encode/Cneo_CH99.broad2_exon.gff",
						"encode/geneid_predictions_duke_exon.gff",
						"encode/cryptococcus_neoformans_JEC21.TIGR.gff"
				};
				
				for (int i = 0; i < fNames.length; i++) {
					System.gc();
					System.out.println(fNames[i]);
					Eduardo myWrapper= new Eduardo(new File(fNames[i]).getAbsolutePath()); // testGTF.gtf
					try {
						myWrapper.read();
					} catch (Exception e) {
						e.printStackTrace(); 
					}
					
					Graph g= myWrapper.getGraph();		// <===== check ENCODE here !!!
					//g.filterNonCodingTranscripts();
					//g.filterCodingTranscripts();
					//g.initTU();
				
					PrintStream pr= null;
					try {
						pr = new PrintStream(new File(fNames[i]+".ana"));
					} catch (FileNotFoundException e1) {
						// TODO Auto-generated catch block
						e1.printStackTrace();
					}
					
					ASAnalyzer.test01_clusters_coverage_as(g, pr);
					pr.println("------------------------------------------------------");
					ASVariation[][] classes= g.getASVariations(1);
					classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
					if (classes!= null) {
						Method m= null;
						try {
							m = classes[0][0].getClass().getMethod("isTrue", null);
						} catch (Exception e) {
							e.printStackTrace();
						} 
						ASVariation[][] filtClasses= (ASVariation[][]) Arrays.filter(classes, m);
						ASAnalyzer.outputVariations(filtClasses, false, false, pr);
					}
					
					pr.flush(); pr.close();
				}
			}

}
