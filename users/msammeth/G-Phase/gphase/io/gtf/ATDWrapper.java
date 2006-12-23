package gphase.io.gtf;

import gphase.algo.ASAnalyzer;
import gphase.db.EnsemblDBAdaptor;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.ENCODE;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Vector;

public class ATDWrapper extends EncodeWrapper {

	void getCDS(GTFObject[] gtfObj) {
		EnsemblDBAdaptor ens= new EnsemblDBAdaptor();
		for (int i = 0; i < gtfObj.length; i++) {
			if (gtfObj[i].getFeature()!= "gene")
				continue;
			String genID= gtfObj[i].getAttribute("ID");	// ID=ENSG00000146433
			ens.re
		}
	}
	
	public static Graph filterForENCODERegions(Graph g) {
		Gene[] ge= g.getGenes();
		DefaultRegion[] regs= ENCODE.getEncodeRegions();
//		for (int i = 0; i < regs.length; i++) 
//			regs[i]= ENCODE.convertToChromosomalCoord(regs[i]);
		for (int i = 0; i < ge.length; i++) {
			int j;
			for (j = 0; j < regs.length; j++) 
				if (regs[j].overlaps(ge[i]))
					break;
			if (j>= regs.length)
				g.removeKill(ge[i]);
		}
		
		return g;
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
		spec.setBuildVersion(17);
		
			// genes
		GTFObject[] gtfs= getGtfObj();
		Vector geneV= new Vector();
		for (int i = 0; i < gtfs.length; i++) {
			if (gtfs[i].getFeature().equalsIgnoreCase("gene")) {
				Gene tmpG= new Gene(spec, gtfs[i].getAttribute("ID"));
				tmpG.setChromosome(gtfs[i].getSeqname());
				tmpG.setStrand(gtfs[i].getStrand());
				tmpG.setStart(gtfs[i].getStart());
				tmpG.setEnd(gtfs[i].getEnd());
				spec.addGene(tmpG);
			}
		}
		
		HashMap tHash= new HashMap(spec.getGeneNb());
		for (int i = 0; i < gtfs.length; i++) {
			if (gtfs[i].getFeature().equalsIgnoreCase("mRNA")) {
				Gene tmpG= spec.getGene(gtfs[i].getAttribute("Parent"));
				Transcript tmp= new Transcript(tmpG, gtfs[i].getAttribute("ID"));
				//tmp.setChromosome(gtfs[i].getSeqname());
				tmp.setStrand(gtfs[i].getStrand());
				tmp.setStart(gtfs[i].getStart());
				tmp.setEnd(gtfs[i].getEnd());
				tmpG.addTranscript(tmp);
				tHash.put(tmp.getTranscriptID(), tmp);
			}
		}

		for (int i = 0; i < gtfs.length; i++) {
			if (!gtfs[i].getFeature().equalsIgnoreCase("mRNA")&& 
					!gtfs[i].getFeature().equalsIgnoreCase("gene")) {
				Transcript tmpT= (Transcript) tHash.get(gtfs[i].getAttribute("Parent"));
				Exon tmp= new Exon(tmpT, gtfs[i].getAttribute("ID"), gtfs[i].getStart(), gtfs[i].getEnd());
				tmpT.addExon(tmp);
			}
		}
		
			// build graph
		Graph g= new Graph();
		g.addSpecies(spec);
		return g;
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
				} catch (Exception e) {
					System.err.println("Invalid GTF format: line "+ lineCtr);
					//e.printStackTrace();
					//continue;
				}
				
					// optional attributes
				int smc= line.indexOf(';');		// GTF2
				if (smc>= 0) {
					String ss= toki.nextToken();
	//				toki= new StringTokenizer(line, " \t");	// must be tab, see specification
	//				for (int i = 0; i < 8; i++) 
	//					ss= toki.nextToken();
	//				String h= line.substring(0, smc);			// last ';'
	//				h= line.substring(0, h.lastIndexOf(' '));	// two ' ' tokens before
	//				h= line.substring(0, h.lastIndexOf(' '));
					String h= line.substring(line.indexOf(ss), line.length()).trim();	// skip that part
					
					toki= new StringTokenizer(h, ";");		// attributes
					while (toki.hasMoreTokens()) {
						h= toki.nextToken().trim();
						int sep= h.indexOf('=');
						if (sep < 0) {						// comments
							String s= h;
							while (toki.hasMoreTokens())
								s+= " "+ toki.nextToken();
							newObj.setComments(s);
						}
						if (sep>= 0) {
							String id= h.substring(0, sep);
							String val= h.substring(sep+ 1, h.length());
							newObj.addAttribute(id, val);
						}
					}
				}
				
//				if (newObj.getFeature().equalsIgnoreCase("mRNA")|| newObj.getFeature().equalsIgnoreCase("gene"))
//					continue;
				gtfVec.add(newObj);
				//System.out.println(gtfVec.size());
			}
			
			gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
		}

	public static void main(String[] args) {
			//"encode/44regions_genes_CHR_coord.gtf"
			// "encode/RefSeqGenes_fromUCSC.gtf"
			// "encode/ASD_27.35a.1.CLASSES_tab.gff"
			
			String fName= "encode/ASD_27.35a.1.CLASSES_tab.gff" ;
			ATDWrapper myWrapper= new ATDWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
			try {
				myWrapper.read();
			} catch (Exception e) {
				e.printStackTrace(); 
			}
			boolean encode= false;
			if (fName.startsWith("encode/44regions_genes_CHR_coord"))
				encode= true;
			
			Graph g= myWrapper.getGraph();		// <===== check ENCODE here !!!
			//g= filterForENCODERegions(g);
			//g.filterNonCodingTranscripts();
			//g.filterCodingTranscripts();
			g.initTU();
	
			PrintStream pr= null;
			try {
				pr = new PrintStream(new File("test04_asd.txt"));
			} catch (FileNotFoundException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			//ASAnalyzer.test01_clusters_coverage_as(g, pr);
			ASAnalyzer.test04_determineVariations_rev(g, pr);
		}

	public ATDWrapper(String absFName) {
		super(absFName);
	}

}
