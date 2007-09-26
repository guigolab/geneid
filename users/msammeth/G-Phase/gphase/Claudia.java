package gphase;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import com.sun.org.apache.xerces.internal.dom.DOMMessageFormatter;

import gphase.io.DomainToGenomeMapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.ASEventold;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.ASVariationWithRegions;
import gphase.model.AbstractRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.File;
import gphase.tools.Formatter;

public class Claudia {
	public static final String ID_EVENT_TAG= "event_id";
	public static final String ID_EXON_ROLE= "exon_role_id";
	public static final String ID_EXON_ROLE_A3= "exon_A3";
	public static final String ID_EXON_ROLE_A= "exon_A";
	public static final String ID_EXON_ROLE_A5= "exon_A5";
	public static final String ID_EXON_ROLE_B= "exon_B";
	public static final String ID_EVENT_FEATURE= "event";
	
	public static int cntExtrem= 0, cntMultiB= 0, cntInvFlanks= 0, cntGTFevents= 0, cntFailedEvents= 0, cntMappedEvents= 0, cntDomainEvents= 0, cntNoDomainEvents= 0; 

	
	public static void main(String[] args) {
	
	}
	
	static DirectedRegion[] getClaudiaExons(File clauFile, String chrID) {
		GTFChrReader clauReader= new GTFChrReader(clauFile.getAbsolutePath());
		clauReader.setChromosomeWise(true);
		clauReader.setReadGene(false);
		clauReader.setReadGTF(true);
		clauReader.setFiltChrIDs(new String[] {chrID});
		clauReader.setReadFeatures(new String[] {
				ID_EXON_ROLE_A3, 
				ID_EXON_ROLE_A5,
				ID_EXON_ROLE_B});
		try {
			clauReader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		GTFObject[] obs= clauReader.getGtfObj();
		if (obs== null)
			return new DirectedRegion[0];
		
			// convert to regions
		DirectedRegion[] regs= new DirectedRegion[obs.length];
		for (int i = 0; i < regs.length; i++) {
			regs[i]= new DirectedRegion(obs[i].getStart(), obs[i].getEnd(), obs[i].getStrand());
			regs[i].addAttribute(GTFObject.TRANSCRIPT_ID_TAG, obs[i].getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
			regs[i].addAttribute(ID_EXON_ROLE, obs[i].getFeature());
			regs[i].addAttribute(ID_EVENT_TAG, obs[i].getSource());
		}
		java.util.Arrays.sort(regs, new AbstractRegion.PositionComparator());
		
		return regs;
	}
	
	static ASVariation[] getCompleteEvents(File inFile, File claudiaFile) {
		long t0= System.currentTimeMillis();
		
		System.out.println("Retrieving AS events from "+inFile.getFileNameOnly()+".");
		GTFChrReader reader= new GTFChrReader(inFile.getAbsolutePath());
		reader.setChromosomeWise(true);
		reader.setReadFeatures(new String[]{"exon", "CDS", "domain"});
		
		
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		
		HashMap<String,Vector<ASVariation>> eventMap= new HashMap<String,Vector<ASVariation>>();
		while (genes!= null&& genes.length> 0) {
			DirectedRegion[] clauExons= getClaudiaExons(claudiaFile, genes[0].getChromosome());
			for (int i = 0; i < genes.length; i++) {
				if (!genes[i].isProteinCoding())
					continue;
				
				ASVariation[] varsSimple= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
				if (varsSimple== null|| varsSimple.length== 0)
					continue;
				
				
				
					// get genomic domain regions
				HashMap<String,Vector<DirectedRegion>> domMap= new HashMap<String,Vector<DirectedRegion>>();
				Transcript[] trpts= genes[i].getTranscripts();
				for (int j = 0; j < trpts.length; j++) {
					Object[] keys= trpts[j].getAttributes().keySet().toArray();
					Vector<DirectedRegion> trptDomV= new Vector<DirectedRegion>();
					for (int m = 0; m < keys.length; m++) {
						DirectedRegion[] diregV= (DirectedRegion[]) trpts[j].getAttribute(keys[m]);							
						int minStart= Integer.MAX_VALUE, maxEnd= Integer.MIN_VALUE;
						String lastID= null;
						for (int k = 0; k < diregV.length; k++) {
							if (diregV[k].get5PrimeEdge()< minStart)
								minStart= diregV[k].get5PrimeEdge();
							if (diregV[k].get3PrimeEdge()> maxEnd)
								maxEnd= diregV[k].get3PrimeEdge();
						}
						DirectedRegion reg= new DirectedRegion();
						reg.setStrand(trpts[j].getStrand());
						reg.setChromosome(trpts[j].getChromosome());
						reg.setStart(minStart);		// TODO check for 0-based
						reg.setEnd(maxEnd);		// -1 +1
						reg.setScore(diregV[0].getScore());
						reg.setID(diregV[0].getID());
						trptDomV.add(reg);
					}
					if (trptDomV.size()> 0) {
						domMap.put(trpts[j].getTranscriptID(),trptDomV);
					}
				}

				// cnt events
				Vector<ASVariationWithRegions> varsDomV= new Vector<ASVariationWithRegions>(); 
				ASVariationWithRegions[] vars= new ASVariationWithRegions[varsSimple.length];
				for (int j = 0; j < vars.length; j++) 
					vars[j]= new ASVariationWithRegions(varsSimple[j], 
							domMap.get(varsSimple[j].getTranscript1().getTranscriptID()),
							domMap.get(varsSimple[j].getTranscript2().getTranscriptID()));
				Comparator domOvlCompi= new ASVariationDomain.StructureComparator();
				for (int j = 0; vars!= null&& j < vars.length; j++) {
					DirectedRegion reg= vars[j].getRegion();
					Vector<DirectedRegion> vRegs= (Vector<DirectedRegion>) domMap.get(vars[j].getTranscript1().getTranscriptID());
					int k;
					for (k = 0; vRegs!= null&& k < vRegs.size(); k++) {
						if (vRegs.elementAt(k).overlaps(reg))
							break;
					}
					if (vRegs!= null&& k< vRegs.size()) {
						varsDomV.add(vars[j]);
						continue;
					}
					vRegs= (Vector<DirectedRegion>) domMap.get(vars[j].getTranscript2().getTranscriptID());
					for (k = 0; vRegs!= null&& k < vRegs.size(); k++) {
						if (vRegs.elementAt(k).overlaps(reg))
							break;
					}
					if (vRegs!= null&& k< vRegs.size()) 
						varsDomV.add(vars[j]);
				}
				if (varsDomV.size()== 0)
					continue;
				vars= (ASVariationWithRegions[]) ASVariation.removeRedundancy((ASVariation[]) Arrays.toField(varsDomV), domOvlCompi);
				for (int j = 0; vars!= null&& j < vars.length; j++) {
					Vector<ASVariation> v= eventMap.remove(vars[j].toString());
					if (v== null)
						v= new Vector<ASVariation>();
					v.add(vars[j]);
					eventMap.put(vars[j].toString(), v);
				}
			}
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
			System.gc();
		}
		
			// output
		ASVariationWithRegions[][] allVars= new ASVariationWithRegions[eventMap.size()][];
		Object[] keys= eventMap.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			allVars[i]= (ASVariationWithRegions[]) Arrays.toField(eventMap.get(keys[i]));
			java.util.Arrays.sort(allVars[i]);
		}
		allVars= (ASVariationWithRegions[][]) Arrays.sort2DFieldRev(allVars);
		
		for (int i = 0; i < allVars.length; i++) 
			System.out.println(allVars[i][0]+"\t"+allVars[i].length);
		return allVars;

	}
	
	HashMap mapTrpts= null, mapDom= null;
	Gene[] genes= null;
	GTFObject[] obs= null;
	ASVariationWithRegions[][][] clauEvents= null;
	
	
	public Claudia(Gene[] newGenes, GTFObject[] newObs) {
		genes= newGenes;
		obs= newObs;
	}
	
	public HashMap getDomainMap() {
		if (mapDom == null) {
			mapDom= new HashMap();
			for (int i = 0; i < genes.length; i++) {
				Transcript[] trpts= genes[i].getTranscripts();
				for (int j = 0; j < trpts.length; j++) {
					if (trpts[j].getAttributes()== null|| trpts[j].getAttributes().size()== 0)
						continue;
					Object[] keys= trpts[j].getAttributes().keySet().toArray();
					Vector<DirectedRegion> trptDomV= new Vector<DirectedRegion>();
					for (int m = 0; m < keys.length; m++) {
						DirectedRegion[] diregV= (DirectedRegion[]) trpts[j].getAttribute(keys[m]);							
						int minStart= Integer.MAX_VALUE, maxEnd= Integer.MIN_VALUE;
						String lastID= null;
						for (int k = 0; k < diregV.length; k++) {
							if (diregV[k].get5PrimeEdge()< minStart)
								minStart= diregV[k].get5PrimeEdge();
							if (diregV[k].get3PrimeEdge()> maxEnd)
								maxEnd= diregV[k].get3PrimeEdge();
						}
						DirectedRegion reg= new DirectedRegion();
						reg.setStrand(trpts[j].getStrand());
						reg.setChromosome(trpts[j].getChromosome());
						reg.setStart(minStart);		// TODO check for 0-based
						reg.setEnd(maxEnd);		// -1 +1
						reg.setScore(diregV[0].getScore());
						reg.setID(diregV[0].getID());
						trptDomV.add(reg);
					}
					if (trptDomV.size()> 0) {
						mapDom.put(trpts[j].getTranscriptID(),trptDomV);
					}
				}
			}
		}

		return mapDom;

	}

	public HashMap getMapTrpts() {
		if (mapTrpts == null) {
			mapTrpts = new HashMap();
			for (int i = 0; genes!= null&& i < genes.length; i++) 
				for (int j = 0; j < genes[i].getTranscriptCount(); j++) 
					mapTrpts.put(genes[i].getTranscripts()[j].getTranscriptID(), genes[i].getTranscripts()[j]);
		}

		return mapTrpts;
	}
	
	public ASVariationWithRegions[] getClauEvents() {
		if (clauEvents== null) {
			HashMap mapTrpt= getMapTrpts(), map= new HashMap(), mapRef= new HashMap();
	
			// collect domain groups, events
			for (int i = 0; i < obs.length; i++) {
				String id= null;
				id= obs[i].getSource();
				int p1= id.indexOf('.');
				int p2= id.lastIndexOf('.');
				if (p1== p2) {
					mapRef.put(id, obs[i].getTranscriptID());	// ref transcript
					continue;	// complete event
				}
				Vector v= (Vector) map.get(id);
				if (v== null)
					v= new Vector();
				v.add(obs[i]);
				map.put(id, v);	// evID
			}
			Vector evV= new Vector(map.size());
			
				// construct complete domains / events
			Object[] keys= map.keySet().toArray();
			cntGTFevents+= keys.length;
			BufferedWriter writer= null;
			try {
				writer= new BufferedWriter(new FileWriter("missed_tIDs_all.txt", true));
			} catch (IOException e) {
				e.printStackTrace();
			}
			for (int i = 0; i < keys.length; i++) {
				Vector v= (Vector) map.remove(keys[i]);
				GTFObject obj= (GTFObject) v.elementAt(0);
				String evID= (String) keys[i];
				int p= evID.lastIndexOf('.');
				
				String tID1= (String) mapRef.get(evID.substring(0,p));	// TODO: ensure that event is read before
				String tID2= obj.getAttribute(GTFObject.TRANSCRIPT_ID_TAG);
				Transcript t1= (Transcript) mapTrpt.get(tID1);
				Transcript t2= (Transcript) mapTrpt.get(tID2);
				if (t1== null|| t2== null) {
					if (tID1== null|| tID2== null) {
						System.out.println("incomplete info for ev "+evID+", transcripts "+tID1+" x "+tID2);
					} else {
						try {
							if (/*tID1.startsWith("NM")&&*/ t1== null) {
								System.out.println("Not found: "+tID1);
								writer.write(tID1+"\n");
							} 
							if (/*tID2.startsWith("NM")&&*/ t2== null) {
								System.out.println("Not found: "+tID2);
								writer.write(tID2+"\n");
							}
							
						} catch (Exception e) {
							e.printStackTrace();
						}
					}
						
					continue;
				}
				++cntMappedEvents;
				
					// get regions
				Vector regA5v= new Vector(), regA3v= new Vector(), regBv= new Vector();
				for (int j = 0; j < v.size(); j++) {
					GTFObject o= (GTFObject) v.elementAt(j);
					DirectedRegion reg= new DirectedRegion(o.getStart(), o.getEnd(), o.getStrand());
					reg.setChromosome(o.getChromosome());
					reg.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
					if (o.getFeature().equals(Claudia.ID_EXON_ROLE_A5))
						regA5v.add(reg);
					else if (o.getFeature().equals(Claudia.ID_EXON_ROLE_A3)) 
						regA3v.add(reg);
					else if (o.getFeature().equals(Claudia.ID_EXON_ROLE_B)) {
						regBv.add(reg);
					}
				}
				
				if (regA5v.size()== 0|| regA3v.size()== 0) {
					++cntExtrem;
					continue;
				}
				if (regBv.size()== 0) {
					System.out.println("WARNING: no B-exon found: "+evID);
					continue;
				}
				if (regBv.size()> 1) {
					++cntMultiB;
				}
								
				
					// build splice chains
				Vector sc1V= new Vector(), sc2V= new Vector();
				if (getClauEvents2(t1, t2, regA5v, regBv, regA3v, sc1V, sc2V)< 0) {
					++cntFailedEvents;
					continue;
				}
				
					// construct event
				SpliceSite[] sChain1= (SpliceSite[]) Arrays.toField(sc1V);
				if (sChain1== null)
					sChain1= new SpliceSite[0];
				SpliceSite[] sChain2= (SpliceSite[]) Arrays.toField(sc2V);
				if (sChain2== null)
					sChain2= new SpliceSite[0];
				ASVariation baseVar= new ASVariation(t1, t2, sChain1, sChain2);
				baseVar.addAttribute(Claudia.ID_EVENT_TAG, evID);
				
				t1= baseVar.getTranscript1();
				t2= baseVar.getTranscript2();
				ASVariationWithRegions var= new ASVariationWithRegions(baseVar, 
						(Vector) getDomainMap().get(t1.getTranscriptID()),
						(Vector) getDomainMap().get(t2.getTranscriptID()));
				if (var.getReg1().length< 1&& var.getReg2().length< 1) {
					++cntNoDomainEvents;
					continue;
				}
				++cntDomainEvents;
				evV.add(var);
			}
			
			clauEvents= (ASVariationWithRegions[]) Arrays.toField(evV);
			Comparator domOvlCompi= new ASVariationDomain.IdentityComparator();
			clauEvents= (ASVariationWithRegions[]) ASVariation.removeRedundancy(clauEvents, domOvlCompi);
			
			try {
				writer.flush();
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
		
		return clauEvents;
	}

	/**
	 * Only reads B-regions and overlaps them with ASTAlaVista events.
	 * @return
	 */
	public ASVariationWithRegions[][][] getClauEvents2() {
		if (clauEvents== null) {
			HashMap mapTrpt= getMapTrpts(), map= new HashMap(), mapRef= new HashMap();

			// collect events
			for (int i = 0; i < obs.length; i++) {
				String id= null;
				id= obs[i].getSource();
				int p1= id.indexOf('.');
				int p2= id.lastIndexOf('.');
				if (p1== p2) {
					mapRef.put(id, obs[i].getTranscriptID());	// ref transcript
					continue;	// complete event
				}
				Vector v= (Vector) map.get(id);
				if (v== null)
					v= new Vector();
				v.add(obs[i]);
				map.put(id, v);	// evID
			}
			Vector evV= new Vector(map.size());

				// prepare claudia events
			Object[] keys= map.keySet().toArray();
			HashMap bRegHash= new HashMap();
			for (int i = 0; i < keys.length; i++) {
				Vector v= (Vector) map.remove(keys[i]);
				Vector regA5v= new Vector(), regA3v= new Vector(), regBv= new Vector();
				for (int j = 0; j < v.size(); j++) {
					GTFObject o= (GTFObject) v.elementAt(j);
					DirectedRegion reg= new DirectedRegion(o.getStart(), o.getEnd(), o.getStrand());	// strand works?!
					reg.setChromosome(o.getChromosome());
					reg.addAttribute(GTFObject.TRANSCRIPT_ID_TAG, reg.getAttribute(GTFObject.TRANSCRIPT_ID_TAG));
					if (o.getFeature().equals(Claudia.ID_EXON_ROLE_A5)|| o.getFeature().equals(Claudia.ID_EXON_ROLE_A))
						regA5v.add(reg);
					else if (o.getFeature().equals(Claudia.ID_EXON_ROLE_A3)) 
						regA3v.add(reg);
					else if (o.getFeature().equals(Claudia.ID_EXON_ROLE_B)) {
						regBv.add(reg);
					}
				}
				
				if (regBv== null|| regBv.size()== 0) {
					System.out.println("WARNING: no B-probes found, "+keys[i]);
					continue;
				}
					
				DirectedRegion r= null;
				int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
				for (int j = 0; j < regBv.size(); j++) {
					r= (DirectedRegion) regBv.elementAt(j);
					if (r.get5PrimeEdge()< min)
						min= r.get5PrimeEdge();
					if (r.get3PrimeEdge()> max)
						max= r.get3PrimeEdge();
				}
				DirectedRegion bReg= new DirectedRegion(min, max, r.getStrand());
				bReg.setChromosome(r.getChromosome());
				bReg.addAttribute(Claudia.ID_EVENT_TAG, keys[i]);
				bRegHash.put(bReg, Arrays.toField(regBv));
			}
			DirectedRegion[] bRegs= (DirectedRegion[]) Arrays.toField(bRegHash.keySet());
			java.util.Arrays.sort(bRegs, new DirectedRegion.EndComparator());	// according to asc end positions, asc start
			DirectedRegion[] bRegStarts= new DirectedRegion[bRegs.length];
			for (int i = 0; i < bRegStarts.length; i++) 
				bRegStarts[i]= bRegs[i];
			java.util.Arrays.sort(bRegStarts, new DirectedRegion.StartComparator());
			
				// get events
			HashMap varsMap= new HashMap();
			for (int i = 0; i < genes.length; i++) {
				
				if (genes[i].getTranscriptCount()> 1000) {
					System.out.println("WARNING: gene on chr "+genes[i].getChromosome()+" with "+genes[i].getTranscriptCount()+" transcritps " +
							"("+genes[i].getTranscripts()[0]+"..)");
					continue;
				}
				
				ASVariation[] singleVars= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
				if (singleVars== null)
					continue;
				ASVariation[][] varsS= (ASVariation[][]) ASMultiVariation.clusterIdenticalEvents(singleVars);
				
				ASVariationWithRegions[][] vars= new ASVariationWithRegions[varsS.length][]; 					
				Comparator compi= new DirectedRegion.EndsBeforeComparator();
				Comparator compi2= new DirectedRegion.StartComparator();
				for (int j = 0; j < varsS.length; j++) {

					vars[j]= new ASVariationWithRegions[varsS[j].length];
						// convert to domain events
					for (int k = 0; k < varsS[j].length; k++) {
						vars[j][k]= new ASVariationWithRegions(varsS[j][k], 
								(Vector) getDomainMap().get(varsS[j][k].getTranscript1().getTranscriptID()),
								(Vector) getDomainMap().get(varsS[j][k].getTranscript2().getTranscriptID()));

					}
					
						// map claudia to astalavista
					DirectedRegion varReg= vars[j][0].getRegion();
					int p= java.util.Arrays.binarySearch(bRegs, varReg, compi);
					if (p< 0)
						p= -(p+1);
					
						// ends
					for (int k = p; k < bRegs.length; k++) {
						if (bRegs[k].get3PrimeEdge()> varReg.get3PrimeEdge())
							break;
						if (!varReg.overlaps(bRegs[k]))
							continue;
						DirectedRegion[] regsB= (DirectedRegion[]) bRegHash.get(bRegs[k]);
						DirectedRegion[] regsVar= vars[j][0].getExonicRegions();
						DirectedRegion[] both= DirectedRegion.intersect(regsB, regsVar);
						if (both!= null&& both.length> 0) {	// TODO
//							int a= 0;	// check for all b-probes are covered by the event
//							for (a = 0; a < regsB.length; a++) {
//								int b= 0;
//								for (b = 0; b < both.length; b++) 
//									if (both[b].overlaps(regsB[a]))
//										break;
//								if (b== both.length)
//									break;
//							}
//							if (a< regsB.length)	// actually they have to be covered in the same transcript variant
//								continue;	
							for (int m = 0; m < vars[j].length; m++) {
								Vector v= null;
								if (vars[j][m].getAttributes()!= null)
									v= (Vector) vars[j][m].getAttributes().remove(Claudia.ID_EVENT_TAG);
								if (v== null)
									v= new Vector();
								v= Arrays.addUnique(v, bRegs[k].getAttribute(Claudia.ID_EVENT_TAG));
								vars[j][m].addAttribute(Claudia.ID_EVENT_TAG, v);
							}
						}
					}

						// starts, exhaustive search
					p= java.util.Arrays.binarySearch(bRegs, varReg, compi2);
					if (p< 0)
						p= -(p+1);
					for (int k = p; k < bRegStarts.length; k++) {
						if (bRegStarts[k].get5PrimeEdge()> varReg.get3PrimeEdge())
							break;
						if (!varReg.overlaps(bRegStarts[k]))
							continue;
						DirectedRegion[] regsB= (DirectedRegion[]) bRegHash.get(bRegStarts[k]);
						DirectedRegion[] regsVar= vars[j][0].getExonicRegions();
						DirectedRegion[] both= DirectedRegion.intersect(regsB, regsVar);
						if (both!= null&& both.length> 0) {
							for (int m = 0; m < vars[j].length; m++) {
								Vector v= null;
								if (vars[j][m].getAttributes()!= null)
									v= (Vector) vars[j][m].getAttributes().remove(Claudia.ID_EVENT_TAG);
								if (v== null)
									v= new Vector();
								v= Arrays.addUnique(v, bRegStarts[k].getAttribute(Claudia.ID_EVENT_TAG));
								vars[j][m].addAttribute(Claudia.ID_EVENT_TAG, v);
							}
						}
					}

						// structural clustering
					Vector v= (Vector) varsMap.remove(vars[j][0].toString());
					if (v== null)
						v= new Vector();
					v.add(vars[j]);
					varsMap.put(vars[j][0].toString(), v);
				}
			}
			clauEvents= new ASVariationWithRegions[varsMap.size()][][];
			keys= varsMap.keySet().toArray();
			for (int i = 0; i < keys.length; i++) 
				clauEvents[i]= (ASVariationWithRegions[][]) Arrays.toField(varsMap.get(keys[i]));
			
		}
		
		return clauEvents;
	}
	
	int getClauEvents1(Transcript t1, Transcript t2, Vector regA5v, Vector regBv, Vector regA3v,
			Vector sc1V, Vector sc2V) {
		// find exons
		//(regA5== null|| reg.get5PrimeEdge()> regA5.get3PrimeEdge()))
		//(regA3== null|| reg.get3PrimeEdge()< regA3.get5PrimeEdge()))
		Exon a5_1= null, a5_2= null; 
		for (int j = regA5v.size()- 1; j >=0; --j) {	// hope that they are sorted
			a5_1= t1.getExon((DirectedRegion) regA5v.elementAt(j));
			a5_2= t2.getExon((DirectedRegion) regA5v.elementAt(j));
			if (a5_1!= null&& a5_2!= null)
				break;
		}
		Exon a3_1= null, a3_2= null; 
		for (int j = 0; j < regA3v.size(); ++j) {	// hope that they are sorted
			a3_1= t1.getExon((DirectedRegion) regA3v.elementAt(j));
			a3_2= t2.getExon((DirectedRegion) regA3v.elementAt(j));
			if (a5_1!= null&& a5_2!= null)
				break;
		}
		if (a5_1== null|| a3_1== null|| a5_2== null|| a3_2== null) {
			++cntInvFlanks;
			return -1;
		}
		Vector b_1v= new Vector(), b_2v= new Vector(); 
		for (int j = 0; j < regBv.size(); ++j) {	// hope that they are sorted
			Exon e= t1.getExon((DirectedRegion) regBv.elementAt(j));
			if (e!= null)
				b_1v.add(e);
			e= t2.getExon((DirectedRegion) regBv.elementAt(j));
			if (e!= null)
				b_2v.add(e);
		}
		if (!(b_1v.size()== 0^ b_2v.size()== 0))
			System.out.println("WARNING: variable part found in both transcripts, "+
					t1.getTranscriptID()+" x "+t2.getTranscriptID());
		Exon[] b_1= (Exon[]) Arrays.toField(b_1v);
		if (b_1== null)
			b_1= new Exon[0];
		Exon[] b_2= (Exon[]) Arrays.toField(b_2v);
		if (b_2== null)
			b_2= new Exon[0];
	
		
		int c1= 0, c2= 0;
			// start with AD in the conserved flank probe, not captured by B
		if (a5_1!= a5_2&& a5_1!= a3_1&& a5_2!= a3_2 &&	// exclude IR
				a5_1.getDonor().getPos()!= a5_2.getDonor().getPos()&&
				((b_1.length> 0&& a5_1!= b_1[0])|| (b_2.length> 0&& a5_2!= b_2[0]))) {
			sc1V.add(a5_1.getDonor());	// starts with AD
			sc2V.add(a5_2.getDonor());
		}
		while (true) {
			if ((b_1.length> 0&& a5_1== b_1[c1]&& a5_1!= a3_1)||
					(b_2.length> 0&& a5_2== b_2[c2]&& a5_2!= a3_2)) {
				sc1V.add(a5_1.getDonor());	// starts with AD
				sc2V.add(a5_2.getDonor());
			} else if ((b_1.length> 0&& a3_1== b_1[c1]&& a5_1!= a3_1)||
					(b_2.length> 0&& a3_2== b_2[c2]&& a5_2!= a3_2)) {
				sc1V.add(a3_1.getAcceptor());	// ends with AA
				sc2V.add(a3_2.getAcceptor());	
			} else {
				if (b_1.length> 0) {
					Exon refLeft= (c1== 0)?a5_1:b_1[c1-1];
					Exon refRight= (c1== b_1.length- 1)?a3_1:b_1[c1+1];
					if (b_1[c1]== refLeft&& b_1[c1]== refRight) {	// a5= b= a3
						sc1V.add(refLeft.getDonor());		// intron ret
						sc1V.add(refRight.getAcceptor());
					} else {								// skipped exon
						sc1V.add(b_1[c1].getAcceptor());
						sc1V.add(b_1[c1].getDonor());
					}
				} else if (b_2.length> 0) {
					Exon refLeft= (c2== 0)?a5_2:b_2[c2-1];
					Exon refRight= (c2== b_2.length- 1)?a3_2:b_2[c2+1];
					if (b_2[c2]== refLeft&& b_2[c2]== refRight) {	// a5= b= a3
						sc2V.add(refLeft.getDonor());		// intron ret
						sc2V.add(refRight.getAcceptor());
					} else {								// skipped exon
						sc2V.add(b_2[c2].getAcceptor());
						sc2V.add(b_2[c2].getDonor());
					}
				}
			}
			if (++c1== b_1.length|| ++c2== b_2.length)
				break;
		}
		// ends with AA in the conserved flank probe, not captured by B
		if (a3_1!= a3_2&& a5_1!= a3_1&& a5_2!= a3_2 && 
				a3_1.getAcceptor().getPos()!= a3_2.getAcceptor().getPos()&&
				((b_1.length> 0&& a3_1!= b_1[b_1.length- 1])|| (b_2.length> 0&& a3_2!= b_2[b_2.length- 1]))) {
			sc1V.add(a3_1.getAcceptor());	// starts with AD
			sc2V.add(a3_2.getAcceptor());
		}
		
		return 0;
	
	}

	/**
	 * find conserved A5 flank, A3 flank and all variable exons in between
	 * @param t1
	 * @param t2
	 * @param regA5v
	 * @param regBv
	 * @param regA3v
	 * @param sc1V
	 * @param sc2V
	 * @return
	 */
	int getClauEvents2(Transcript t1, Transcript t2, Vector regA5v, Vector regBv, Vector regA3v,
			Vector sc1V, Vector sc2V) {
		// find exons
		//(regA5== null|| reg.get5PrimeEdge()> regA5.get3PrimeEdge()))
		//(regA3== null|| reg.get3PrimeEdge()< regA3.get5PrimeEdge()))
		Exon a5_1= null, a5_2= null; 
		for (int j = regA5v.size()- 1; j >=0; --j) {	// hope that they are sorted
			a5_1= t1.getExon((DirectedRegion) regA5v.elementAt(j));
			a5_2= t2.getExon((DirectedRegion) regA5v.elementAt(j));
			if (a5_1!= null&& a5_2!= null)
				break;
		}
		Exon a3_1= null, a3_2= null; 
		for (int j = 0; j < regA3v.size(); ++j) {	// hope that they are sorted
			a3_1= t1.getExon((DirectedRegion) regA3v.elementAt(j));
			a3_2= t2.getExon((DirectedRegion) regA3v.elementAt(j));
			if (a5_1!= null&& a5_2!= null)
				break;
		}
		Comparator compi= new Exon.PositionSSComparator();
		if (a5_1== null|| a3_1== null|| a5_2== null|| a3_2== null|| (compi.compare(a5_1,a3_1)== 0&& compi.compare(a5_2,a3_2)== 0)) {
			++cntInvFlanks;
			System.out.println("WARNING: at least one not valid flank found.");
			return -1;
		}
		Vector b_1v= new Vector(), b_2v= new Vector(); 
		Exon[] e= t1.getExons(a5_1, a3_1);
		for (int i = 0;e!= null&& i < e.length; i++) 
			b_1v.add(e[i]);
		e= t2.getExons(a5_2, a3_2);
		for (int i = 0;e!= null&& i < e.length; i++) 
			b_2v.add(e[i]);
		if (b_1v.size()> 0&& b_2v.size()> 0) {	// correct for missed flanks
			for (int i = 0; i < b_1v.size()&& i< b_2v.size(); i++) {
				Exon reg1= (Exon) b_1v.elementAt(i);
				Exon reg2= (Exon) b_2v.elementAt(i);
				if (reg1.overlaps(reg2)) {
					a5_1= reg1; a5_2= reg2;
					b_1v.remove(i); b_2v.remove(i--);
				}
			}
			for (int i = 0; i < b_1v.size()&& i< b_2v.size(); i++) {
				Exon reg1= (Exon) b_1v.elementAt(b_1v.size()- 1- i);
				Exon reg2= (Exon) b_2v.elementAt(b_2v.size()- 1- i);
				if (reg1.overlaps(reg2)) {
					a3_1= reg1; a3_2= reg2;
					b_1v.remove(b_1v.size()- 1- i); b_2v.remove(b_2v.size()- 1- i--);
				}
			}
		}
//		if (b_1v.size()> 0&& b_2v.size()> 0)
//			System.out.println("WARNING: variable part found in both transcripts, "+
//					t1.getTranscriptID()+" x "+t2.getTranscriptID());
		Exon[] b_1= (Exon[]) Arrays.toField(b_1v);
		if (b_1== null)
			b_1= new Exon[0];
		Exon[] b_2= (Exon[]) Arrays.toField(b_2v);
		if (b_2== null)
			b_2= new Exon[0];

		
		int c1= 0, c2= 0;
			// start with AD in the conserved flank probe, not captured by B
		if (compi.compare(a5_1,a5_2)!= 0&& compi.compare(a5_1,a3_1)!= 0&& compi.compare(a5_2,a3_2)!= 0 &&	// exclude IR
				a5_1.getDonor().getPos()!= a5_2.getDonor().getPos()) {
			sc1V.add(a5_1.getDonor());	// starts with AD
			sc2V.add(a5_2.getDonor());
		}
		while (c1< b_1.length|| c2< b_2.length) {
			if (b_1.length> 0&& c1< b_1.length) {
				sc1V.add(b_1[c1].getAcceptor());	// skipped exon
				sc1V.add(b_1[c1].getDonor());
				++c1;
			} 
			if (b_2.length> 0&& c2< b_2.length) {
				sc2V.add(b_2[c2].getAcceptor());
				sc2V.add(b_2[c2].getDonor());
				++c2;
			}
		}
		// ends with AA in the conserved flank probe, not captured by B
		if (compi.compare(a3_1,a3_2)!= 0&& compi.compare(a5_1,a3_1)!= 0&& compi.compare(a5_2,a3_2)!= 0 && 
				a3_1.getAcceptor().getPos()!= a3_2.getAcceptor().getPos()) {
			sc1V.add(a3_1.getAcceptor());	// ends with AA
			sc2V.add(a3_2.getAcceptor());
		}
		// pure intron retention events
		if (b_1.length== 0&& b_2.length== 0) {
			if (compi.compare(a5_1,a3_1)== 0&& compi.compare(a5_2,a3_2)!= 0) {
				sc2V.add(a5_2.getDonor());
				sc2V.add(a3_2.getAcceptor());
			}
			else if (compi.compare(a5_2,a3_2)== 0&& compi.compare(a5_1,a3_1)!= 0) {
				sc1V.add(a5_1.getDonor());
				sc1V.add(a3_1.getAcceptor());
			}
		}
		
			// check
		if ((sc1V.size()== 0&& sc2V.size()== 0)|| ((sc1V.size()+sc2V.size())%2!= 0)) {
			System.out.println("WARNING: invalid event, giving up.");
			return -1;
		}
		for (int i = 0; i < sc1V.size(); i++) 
			if (sc1V.elementAt(i)== null) {
				System.out.println("WARNING: invalid event, got to give up.");
				return -1;
			}
		for (int i = 0; i < sc2V.size(); i++) 
			if (sc2V.elementAt(i)== null) {
				System.out.println("WARNING: invalid event, I have to give up.");
				return -1;
			}
		
		return 0;
	
	}
	
	public ASVariationWithRegions[] getASEvents() {
		
	}
	
}
