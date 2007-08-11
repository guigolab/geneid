package gphase.io;

import gphase.Analyzer;
import gphase.Constants;
import gphase.algo.ASAnalyzer;
import gphase.algo.AlgoHandler;
import gphase.db.MapTable;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.DoubleVector;
import gphase.tools.Formatter;
import gphase.tools.IntVector;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class DomainWrapper extends DefaultIOWrapper {
	
	final static String ID_START_PROTEIN= ">mRNA_NCBI_ID:";
	final static String ID_NO_HITS= "No hits above threshold.";
	
	public static void main(String[] args) {
		
		DomainWrapper myDomainWrapper= new DomainWrapper("douglas"+File.separator+"H.sapiens.All.mRNAs.with.CDS.pfam.parsed");
		PrintStream p= System.out;
		try {
			p= new PrintStream("."+File.separator+"ovlDomains.txt");
			myDomainWrapper.checkOverlappingDomains(p);
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (1== 1)
			System.exit(0);
		
		
		
		System.out.print("Loading domains..");
		System.out.flush();
		DouglasDomainWrapper wrapper= new DouglasDomainWrapper(
				"douglas"+File.separator+"Dmela20070107.mif.seq.pfam.parsed");
		try {
			wrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}		
		HashMap map= wrapper.getMap();
		Object[] o= map.keySet().toArray();
		String[] someIDs= new String[o.length];
		int cntDomains= 0;
		for (int i = 0; i < o.length; i++) { 
//			System.out.print(o[i]+"\t");
			Vector v= (Vector) map.get(o[i]); 
			cntDomains+= v.size(); 
//			for (int j = 0; j < v.size(); j++) 
//				System.out.print(v.elementAt(j)+",");
//			System.out.println();
			someIDs[i]= (String) o[i];
		}
		System.out.println("found "+someIDs.length+" genes, "+cntDomains+" domains.\n");
		
		
			// get basic transcripts
		String absFName= Constants.getLatestUCSCAnnotation("fruitfly", "RefSeq", null);
		System.out.println("Searching for reference transcripts in "+absFName+".");
		GTFChrReader reader= new GTFChrReader(absFName);
		reader.setChromosomeWise(false);
		reader.setFiltSomeIDs(someIDs);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		int cntMultiTrpt= 0;
		int totNTcovByORFs= 0;
		for (int i = 0; i < genes.length; i++) {
			totNTcovByORFs+= genes[i].getTranscripts()[0].getTranslations()[0].getLength();
			if (genes[i].getTranscriptCount()> 1)
				++cntMultiTrpt;
		}
		System.out.println("Found "+genes.length+" genes ("+cntMultiTrpt+"multi-transcript).\n");		
//		System.out.print("Not found: ");
//		String[] notFound= reader.getFiltSomeIDsNotFound();
//		for (int i = 0; i < notFound.length; i++) 
//			System.out.print(notFound[i]+" ");
//		System.out.println();
		
		
		
			// map AA -> NT
		System.out.print("Mapping domain info to genome..");
		System.out.flush();
		Vector vRegs= new Vector();
		HashMap mapProtTID= new HashMap(someIDs.length);
		int totNTcovByDom= 0;
		for (int i = 0; i < someIDs.length; i++) {
			int x;
			for (x = 0; x < genes.length; x++) {
				assert(genes[x].getTranscripts().length== 1);
				String[] protIDs= genes[x].getTranscripts()[0].getTranslations()[0].getProteinIDsAll();
				int k;
				for (k = 0; k < protIDs.length; k++) 
					if (protIDs[k].equals(someIDs[i]))
						break;
				if (k< protIDs.length)
					break;
			}
			if (x== genes.length) 
				continue;
			mapProtTID.put(someIDs[i], genes[x].getTranscripts()[0].getTranscriptID());
			Vector v= (Vector) map.get(someIDs[i]);
			for (int j = 0; j < v.size(); j++) {
				DefaultRegion reg= (DefaultRegion) v.elementAt(j);
				Transcript trpt= genes[x].getTranscripts()[0];				
				int start= trpt.getTranslations()[0].getGenomicPosition((reg.getStart()- 1)*3);		// TODO check for 0-based
				int end= trpt.getTranslations()[0].getGenomicPosition((reg.getEnd())*3);		// -1 +1
				DirectedRegion genReg= new DirectedRegion(start, end, trpt.getStrand());
				genReg.setChromosome(trpt.getChromosome());
				genReg.setID(reg.getID());
				genReg.setScore(reg.getScore());
				vRegs.add(genReg);
				totNTcovByDom+= genReg.getLength();
			}
		}
		System.out.println("mapped "+vRegs.size()+" domains.");
		
		
		
			// get all transcripts in the regions
		System.out.print("Getting all transcripts in locus of domain..");
		System.out.flush();
		reader= new GTFChrReader(Constants.getLatestUCSCAnnotation("fruitfly", "RefSeq", null));
		DirectedRegion[] diregs= (DirectedRegion[]) Arrays.toField(vRegs);
		reader.setFiltRegs(diregs);
		reader.setSilent(false);
		reader.setChromosomeWise(false);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		genes= reader.getGenes();
		int cntDomSingle= 0, cntSingle= 0;
		for (int i = 0; i < someIDs.length; i++) {	// another time, map to protIDs
			String tid= (String) mapProtTID.get(someIDs[i]);
			int x;
			for (x = 0; x < genes.length; x++) {
				int j;
				for (j = 0; j < genes[x].getTranscripts().length; j++) {
					if (genes[x].getTranscripts()[j].getTranscriptID().equals(tid))
						break;
				}
				if (j< genes[x].getTranscripts().length)
					break;
			}
			if (x== genes.length) 
				continue;
			
			int a= 0;
			for (int j = 0; j < genes[x].getTranscripts().length; j++) 
				if (genes[x].getTranscripts()[j].isCoding())
					++a;
			int b= 0;
			for (int j = 0; j < genes[x].getExons().length; j++) 
				if (genes[x].getExons()[j].isCoding())
					++b;
			
			if (a< 2|| b< 2) {
				++cntSingle;
				cntDomSingle+= ((Vector) map.get(someIDs[i])).size();
			}
		}
		System.out.println("found "+(genes.length- cntSingle)+" genes ("+cntDomSingle+
				" domains) with multiple transcripts/exons.\n");
		
			// retrieve events
		System.out.print("Retrieving events..");
		System.out.flush();
		HashMap localMap= new HashMap();
		Vector vASNDom= new Vector();
		int cntASgenes= 0, cntASEvents= 0, cntASinDom= 0;
		HashMap mapLandscape= new HashMap(), mapLandscapeNonDom= new HashMap(), mapLandscapeAll= new HashMap();
		for (int i = 0; i < genes.length; i++) { 
			ASVariation[][] vars= Analyzer.getASVariations(genes[i], ASMultiVariation.FILTER_HIERARCHICALLY, ASVariation.TYPE_ALL);
			int oldCntASEvents= cntASEvents;
			for (int x = 0; vars!= null&& x < vars.length; x++) {
				for (int xx = 0; xx < vars[x].length; xx++) {
					if (!vars[x][xx].isContainedCDS())
						continue;
					++cntASEvents;
					DirectedRegion asReg= vars[x][xx].getRegion();
					int y;
					for (y = 0; y < diregs.length; y++) 
						if (asReg.overlaps(diregs[y])) {
							++cntASinDom;
							Vector v= (Vector) localMap.get(diregs[y]);
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							localMap.put(diregs[y], v);
							
							v= (Vector) mapLandscape.get(vars[x][xx].toString());
							if (v== null)
								v= new Vector();
							v.add(vars[x][xx]);
							mapLandscape.put(vars[x][xx].toString(), v);
							break;
						}
					if (y== diregs.length) {
						vASNDom.add(vars[x]);
						
						Vector v= (Vector) mapLandscapeNonDom.get(vars[x][xx].toString());
						if (v== null)
							v= new Vector();
						v.add(vars[x][xx]);
						mapLandscapeNonDom.put(vars[x][xx].toString(), v);
					}
					Vector v= (Vector) mapLandscapeAll.get(vars[x][xx].toString());
					if (v== null)
						v= new Vector();
					v.add(vars[x][xx]);
					mapLandscapeAll.put(vars[x][xx].toString(), v);
				}
			}
			if (cntASEvents!= oldCntASEvents)
				++cntASgenes;
		}
		System.out.println("found "+cntASgenes+" genes with AS.\n\n");
		
		
		System.out.println("Domain\tevents");
		o= localMap.keySet().toArray();
		int cntDomAS= 0;
		for (int i = 0; i < o.length; i++) {
			System.out.print(o[i]+"\t");
			Vector v= (Vector) localMap.get(o[i]);
			if (v== null|| v.size()== 0) 
				continue;
			++cntDomAS;
			for (int j = 0; j < v.size(); j++) 
				System.out.print(v.elementAt(j)+",");
			System.out.println();
		}
		System.out.println(cntDomAS+" domains containing at least 1 AS event.\n");
		
		System.out.println("Landscape Domains");
		Vector vv= new Vector();
		Object[] keys= mapLandscape.keySet().toArray();
		for (int i = 0; i < keys.length; i++) 
			vv.add(mapLandscape.get(keys[i]));
		ASVariation[][] vars= (ASVariation[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		ASAnalyzer.outputVariations(vars, false, false, System.out);
		System.out.println("\n");
		
		System.out.println("Landscape NonDomains");
		vv= new Vector();
		keys= mapLandscapeNonDom.keySet().toArray();
		for (int i = 0; i < keys.length; i++) 
			vv.add(mapLandscapeNonDom.get(keys[i]));
		vars= (ASVariation[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		ASAnalyzer.outputVariations(vars, false, false, System.out);
		System.out.println("\n");
		
		System.out.println("Landscape All");
		vv= new Vector();
		keys= mapLandscapeAll.keySet().toArray();
		for (int i = 0; i < keys.length; i++) 
			vv.add(mapLandscapeAll.get(keys[i]));
		vars= (ASVariation[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		ASAnalyzer.outputVariations(vars, false, false, System.out);
		System.out.println("\n");

		System.out.println();
		float percGenSingle= ((genes.length- cntSingle)*100f)/ genes.length;
		float percGenAS= (cntASgenes*100f)/ (genes.length- cntSingle);
		System.out.println(Formatter.fprint(percGenSingle, 2)+"% multi transcript genes, of them "+
				Formatter.fprint(percGenAS, 2)+"% with AS.");
		
		float percDomSingle= ((cntDomains- cntDomSingle)*100f)/ cntDomains;  
		float percDomAS= (cntDomAS*100f)/ (cntDomains- cntDomSingle);  
		System.out.println(Formatter.fprint(percDomSingle, 2)+"% domains in loci with multiple transcripts, of them "+
				Formatter.fprint(percDomAS, 2)+"% with AS.");
		
		float percASinDom= (cntASinDom*100f)/ cntASEvents;  
		float percNTcov= (totNTcovByDom*100f)/ totNTcovByORFs;  
		System.out.println(Formatter.fprint(percASinDom, 2)+"% of AS events are in domains, which occupy a fraction of NT of "+
				Formatter.fprint(percNTcov, 2)+"%");
		
	}
	
	public DomainWrapper(String absFName) {
		super(absFName);
	}
	
	HashMap map, mapConflict; 
	
	public boolean isApplicable() throws Exception {
		// TODO Auto-generated method stub
		return false;
	}

	public void read() throws Exception {
		
		map= new HashMap();
		BufferedReader buffy= new BufferedReader(new FileReader(fPath+ 
				File.separator+ fName));
		String line= buffy.readLine();
		while (buffy.ready()&& !line.startsWith(ID_START_PROTEIN))
			line= buffy.readLine();
		DoubleVector remScoreV= new DoubleVector();
		HashMap mapNames= new HashMap();
		File ovlFile= new File("domains_overlap.txt");
		if (ovlFile.exists())
			ovlFile.delete();
		int cntTrpts= 0, cntDomains= 0, cntDomLeft= 0;
		while (buffy.ready()) {
			++cntTrpts;
			String protID= line.substring(ID_START_PROTEIN.length());
			Vector v= new Vector();
			line= buffy.readLine();
			mapConflict= new HashMap();
			while (buffy.ready()&& !line.startsWith(ID_START_PROTEIN)) {
				if (line.startsWith(ID_NO_HITS)) {
					DirectedRegion dummy= new DirectedRegion();
					dummy.addAttribute(DomainToGenomeMapper.DOMAIN_ANNOTATION_ATTRIBUTE_ID, DomainToGenomeMapper.DOMAIN_ANNOTATION_DOM_NONE);
					v.add(dummy);
				} else {
					++cntDomains;
					String[] tokens= line.split("\t");
					int start= Integer.parseInt(tokens[1]);
					int end= Integer.parseInt(tokens[2]);
					double score= Double.parseDouble(tokens[3]);
					DefaultRegion reg= new DefaultRegion(start, end);
					reg.setID(tokens[0]);
					if (mapNames.get(tokens[0])== null)
						mapNames.put(tokens[0], tokens[0]);
					reg.setScore(score);
					reg.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, protID);
					reg.addAttribute(DomainToGenomeMapper.DOMAIN_ANNOTATION_ATTRIBUTE_ID, DomainToGenomeMapper.DOMAIN_ANNOTATION_DOM_ANNOTATED);
					boolean add= true;	// greedily add
					for (int i = 0; i < v.size(); i++) {
						DefaultRegion reg2= (DefaultRegion) v.elementAt(i);
						if (reg.overlaps(reg2)) {
							Vector vConfl= (Vector) mapConflict.remove(reg);	// collect conflicting ones
							Vector vConfl2= (Vector) mapConflict.remove(reg2);	// collect conflicting ones
							if (vConfl== null)
								vConfl= new Vector();
							for (int j = 0; vConfl2!= null&& j < vConfl2.size(); j++) 	// joing
								vConfl= Arrays.addUnique(vConfl, vConfl2.elementAt(j));
							vConfl= Arrays.addUnique(vConfl, reg);
							mapConflict.put(reg, vConfl);
							mapConflict.put(reg2, vConfl);
							if (reg.getScore()< reg2.getScore())
								add= true;
							else
								add= false;
						}
					}
					if (true) {	// always add
						for (int i = 0; i < v.size(); i++) {
							DefaultRegion reg2= (DefaultRegion) v.elementAt(i);
							if (reg.overlaps(reg2)) {
								remScoreV.add(reg2.getScore());
								//v.remove(i--);	// always add
							}
						}
						v.add(reg);
					} else {
						remScoreV.add(reg.getScore());
					}
				}
				line= buffy.readLine();
			}
			if (line.startsWith(ID_NO_HITS)&& buffy.ready())
				line= buffy.readLine();
			if (v.size()> 0) 
				map.put(protID, v);
			cntDomLeft+= v.size();
			
			BufferedWriter writerOvl= new BufferedWriter(new FileWriter(ovlFile.getAbsolutePath()));
			Iterator iter= mapConflict.values().iterator();
			while (iter.hasNext()) {
				Vector vv= (Vector) iter.next();
				for (int i = 0;vv!= null&& i < vv.size(); i++) {
					DefaultRegion reg2= (DefaultRegion) vv.elementAt(i);
					writerOvl.write(protID+"\t"+reg2.getID()+"\t"+reg2.getStart()+"\t"+reg2.getEnd()+"\t"+reg2.getScore()+"\n");
				}
			}
			writerOvl.flush(); writerOvl.close();

		}
	
		Distribution dist= new Distribution(remScoreV.toDoubleArray());
		System.out.println("Read "+(cntDomains)+" domains ("+
				mapNames.size()+" different types, "+cntTrpts+" transcripts), of which I kept "+
				remScoreV.size()+" domains although they overlap (med score "+Formatter.fprint(dist.getMedian(), 5)+
				", details in "+ovlFile.getAbsolutePath()+"), left "+
				map.size()+" domains.");
		
	
	}

	public void checkOverlappingDomains(PrintStream outP) throws Exception {
		
		map= new HashMap();
		BufferedReader buffy= new BufferedReader(new FileReader(fPath+ 
				File.separator+ fName));
		String line= buffy.readLine();
		while (buffy.ready()&& !line.startsWith(ID_START_PROTEIN))
			line= buffy.readLine();
		DoubleVector remScoreV= new DoubleVector();
		while (buffy.ready()) {
			String protID= line.substring(ID_START_PROTEIN.length());
			Vector v= new Vector();
			line= buffy.readLine();
			while (buffy.ready()&& !line.startsWith(ID_START_PROTEIN)&& 
					!line.startsWith(ID_NO_HITS)) {
				String[] tokens= line.split("\t");
				int start= Integer.parseInt(tokens[1]);
				int end= Integer.parseInt(tokens[2]);
				double score= Double.parseDouble(tokens[3]);
				DefaultRegion reg= new DefaultRegion(start, end);
				reg.setID(tokens[0]);
				reg.setScore(score);
				v.add(reg);
				line= buffy.readLine();
			}
			if (line.startsWith(ID_NO_HITS)&& buffy.ready())
				line= buffy.readLine();
			if (v.size()> 0) {
				Vector allOvlV= new Vector();
				for (int i = 0; i < v.size(); i++) {
					Vector ovlV= new Vector();
					for (int j = i+1; j < v.size(); j++) {
						DefaultRegion reg1= (DefaultRegion) v.elementAt(j-1);	// assuming order on the domains
						DefaultRegion reg2= (DefaultRegion) v.elementAt(j);
						if (reg1.overlaps(reg2))
							ovlV.add(reg2);
						else
							break;
					}
					if (ovlV.size()> 0) {
						ovlV.insertElementAt(v.elementAt(i), 0);	// 1st overlapping
						allOvlV.add(ovlV);
						i+= v.size()- 1;
					} 
				}
				if (allOvlV.size()> 0)
					map.put(protID, allOvlV);
			}
		}

		outP.println("Overlapping domains:");
		Object[] keys= map.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			Vector allOvlV= (Vector) map.get(keys[i]);
			outP.println("\n"+keys[i]+":");
			for (int j = 0; j < allOvlV.size(); j++) {
				Vector v= (Vector) allOvlV.elementAt(j);
				for (int k = 0; k < v.size(); k++) {
					outP.println(v.elementAt(k));
				}
				outP.println();
			}
		}
		outP.flush(); outP.close();
	}

	public void write() throws Exception {
		// TODO Auto-generated method stub

	}

	public HashMap getMap() {
		return map;
	}
	
	public HashMap getMapConflict() {
		if (mapConflict== null)
			return null;
		Object[] values= mapConflict.values().toArray();
		HashMap mp= new HashMap();
		for (int i = 0; i < values.length; i++)	// make non-redundant 
			if (mp.get(values[i])== null)
				mp.put(values[i], null);
		values= mp.keySet().toArray();
		mapConflict= new HashMap();
		for (int i = 0; i < values.length; i++) {
			Vector regV= (Vector) values[i];
			mapConflict.put(
					((DefaultRegion) regV.elementAt(0)).getAttribute(GTFObject.TRANSCRIPT_ID_TAG),
					regV);
		}
		return mapConflict;
	}

}
