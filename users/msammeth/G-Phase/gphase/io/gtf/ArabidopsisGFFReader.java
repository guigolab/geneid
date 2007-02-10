/*
 * Created on Feb 8, 2007
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.io.gtf;

import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Vector;

public class ArabidopsisGFFReader extends EncodeWrapper {
	
	public ArabidopsisGFFReader(String absFName) {
		super(absFName);
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
			if (toki.countTokens()< 1)
				continue;	// empty lines at file concat points
			// <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
			GTFObject newObj= createGTFObject();
			try {				
				String s= toki.nextToken();
				int p= s.indexOf(":TIGR5:");				
				newObj.seqname= s.substring(p+7,p+8);	// chromosome
				newObj.source= toki.nextToken();
				toki.nextToken();	// empty '.' !!!
				
				s= toki.nextToken();
				newObj.start= Integer.parseInt(s);
				s= toki.nextToken();
				newObj.end= Integer.parseInt(s);
				newObj.setScore(toki.nextToken());
				newObj.setStrand(toki.nextToken());
				newObj.setFrame(toki.nextToken());
				
					// gene_id=...; transcript_id=...; exon_id=...
				while (toki.hasMoreTokens()) {
					s= toki.nextToken().trim();
					p= s.indexOf('=');
					if (p>= 0) {
						String key= s.substring(0, p);
						String val=	s.substring(p+1, s.length()-1);	// trailing ';'
						if (key.equalsIgnoreCase("gene_id")) {
							newObj.addAttribute(GTFObject.GENE_ID_TAG, val);
							if (newObj.getFeature()== null|| newObj.getFeature().equals(""))
								newObj.setFeature("gene");
						} 
						if (key.equalsIgnoreCase("transcript_id")) {
							newObj.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, val);
							if (newObj.getFeature()== null|| newObj.getFeature().equals(""))
								newObj.setFeature("mRNA");
						} 
						if (key.equalsIgnoreCase("exon_id")) {
							newObj.addAttribute(GTFObject.EXON_ID_TAG, val);
							if (newObj.getFeature()== null|| newObj.getFeature().equals(""))
								newObj.setFeature("exon");
						}
					}
				}
				
			} catch (Exception e) {
				System.err.println("Invalid GTF format (line "+ lineCtr+"): "+e.toString());
				//e.printStackTrace();
				//continue;
			}
			
			if (newObj.getFeature()!= null&& (!newObj.getFeature().equals("")))
				gtfVec.add(newObj);
			//System.out.println(gtfVec.size());
		}
		
		gtfObj= (GTFObject[]) Arrays.toField(gtfVec);
	}

	Graph assemble() {
		
		Species spec= new Species("cress");
	
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
						transcript.setBoundaries(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
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

	public Graph getGraph() {
		try {
			read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		
		//getCDS(gtfObj);
		
		return assemble();
	}

	GTFObject createGTFObject(){
		return new GTFObject();
	}
	
	

}
