package gphase;

import java.awt.Color;
import java.awt.Component;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import javax.imageio.ImageIO;

import org.freehep.util.export.ExportFileType;

import gphase.gui.Circle;
import gphase.gui.CopyOfSpliceOSigner;
import gphase.gui.SpliceOSigner;
import gphase.io.DomainWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.ASVariationWithRegions;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.DoubleVector;
import gphase.tools.File;
import gphase.tools.Formatter;

public class StructuratorDomains extends AStaLaVista {
	
	private static boolean checkGene(Gene g) {
		for (int i = 0; i < g.getTranscripts().length; i++) {
			if (g.getTranscripts()[i].getTranscriptID().startsWith("NM_004396"))
				return true;
		}
		return false;
	}
	
	
	public static ASVariationWithRegions[][] getEvents(File inFile) {
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
				while (genes!= null) {
					for (int i = 0; i < genes.length; i++) {
						if (!genes[i].isProteinCoding())
							continue;
	//					if (!checkGene(genes[i]))
	//						continue;
	
						ASVariation[] varsSimple= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
						if (varsSimple== null|| varsSimple.length== 0)
							continue;
						
	
							// get genomic domain regions
						HashMap<String,Vector<DirectedRegion>> domMap= new HashMap<String,Vector<DirectedRegion>>();
						Transcript[] trpts= genes[i].getTranscripts();
						for (int j = 0; j < trpts.length; j++) {
							Vector<DirectedRegion> trptDomV= new Vector<DirectedRegion>();
							if (trpts[j].getAttributes()!= null) {
								Object[] keys= trpts[j].getAttributes().keySet().toArray();
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
							}
							for (int k = 0; k < trptDomV.size(); k++) {
								for (int m = k+1; m < trptDomV.size(); m++) {
									if (((DirectedRegion) trptDomV.elementAt(m)).overlaps((DirectedRegion) trptDomV.elementAt(k)))
										System.out.println("WARNING: transcript "+trpts[j].getTranscriptID()+" has overlapping domains.");
								}
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
						Comparator domOvlCompi= new ASVariationDomain.IdentityComparator();
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
							//if (vars[j].is_contained_in_CDS())		// ***
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


	public static ASVariationWithRegions[][][] getEventsNew(File inFile) {
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
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) {
					if (!genes[i].isProteinCoding())
						continue;
//					if (!checkGene(genes[i]))
//						continue;

					ASVariation[] varsSimple= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
					if (varsSimple== null|| varsSimple.length== 0)
						continue;
					ASVariation[][] varsSimpleCl= ASMultiVariation.clusterIdenticalEvents(varsSimple);

						// get genomic domain regions
					HashMap<String,Vector<DirectedRegion>> domMap= new HashMap<String,Vector<DirectedRegion>>();
					Transcript[] trpts= genes[i].getTranscripts();
					for (int j = 0; j < trpts.length; j++) {
						Vector<DirectedRegion> trptDomV= new Vector<DirectedRegion>();
						if (trpts[j].getAttributes()!= null) {
							Object[] keys= trpts[j].getAttributes().keySet().toArray();
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
						}
						for (int k = 0; k < trptDomV.size(); k++) {
							for (int m = k+1; m < trptDomV.size(); m++) {
								if (((DirectedRegion) trptDomV.elementAt(m)).overlaps((DirectedRegion) trptDomV.elementAt(k)))
									System.out.println("WARNING: transcript "+trpts[j].getTranscriptID()+" has overlapping domains.");
							}
						}
						
						if (trptDomV.size()> 0) {
							domMap.put(trpts[j].getTranscriptID(),trptDomV);
						}
					}
	
					// cnt events
					ASVariationWithRegions[][] vars= new ASVariationWithRegions[varsSimpleCl.length][];
					for (int j = 0; j < vars.length; j++) {
						vars[j]= new ASVariationWithRegions[varsSimpleCl[j].length];
						for (int k = 0; k < vars[j].length; k++) {
							vars[j][k]= new ASVariationWithRegions(varsSimpleCl[j][k], 
									domMap.get(varsSimpleCl[j][k].getTranscript1().getTranscriptID()),
									domMap.get(varsSimpleCl[j][k].getTranscript2().getTranscriptID()));
						}
					}
					Comparator domOvlCompi= new ASVariationDomain.IdentityComparator();
					Vector domClV= new Vector();
					for (int j = 0; vars!= null&& j < vars.length; j++) {
						Vector<ASVariationWithRegions> varsDomV= new Vector<ASVariationWithRegions>(); 
						for (int k = 0; k < vars[j].length; k++) {
							DirectedRegion reg= vars[j][k].getRegion();
							Vector<DirectedRegion> vRegs= (Vector<DirectedRegion>) domMap.get(vars[j][k].getTranscript1().getTranscriptID());
							int m;
							for (m = 0; vRegs!= null&& m < vRegs.size(); m++) {
								if (vRegs.elementAt(m).overlaps(reg))
									break;
							}
							if (vRegs!= null&& m< vRegs.size()) {
								varsDomV.add(vars[j][k]);
								continue;
							}
							vRegs= (Vector<DirectedRegion>) domMap.get(vars[j][k].getTranscript2().getTranscriptID());
							for (m = 0; vRegs!= null&& m < vRegs.size(); m++) {
								if (vRegs.elementAt(m).overlaps(reg))
									break;
							}
							if (vRegs!= null&& k< vRegs.size())
							//if (vars[j].is_contained_in_CDS())		// ***
								varsDomV.add(vars[j][k]);
						}
						if (varsDomV.size()!= 0)
							domClV.add(varsDomV);
					}
					
					if (domClV.size()== 0)
						continue;
					vars= (ASVariationWithRegions[][]) Arrays.toField(domClV); 
					//ASVariation.removeRedundancy((ASVariation[]) Arrays.toField(varsDomV), domOvlCompi);
					for (int j = 0; vars!= null&& j < vars.length; j++) {
						Vector v= eventMap.remove(vars[j][0].toString());
						if (v== null)
							v= new Vector<ASVariation>();
						v.add(vars[j]);
						eventMap.put(vars[j][0].toString(), v);
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
			ASVariationWithRegions[][][] allVars= new ASVariationWithRegions[eventMap.size()][][];
			Object[] keys= eventMap.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				allVars[i]= (ASVariationWithRegions[][]) Arrays.toField(eventMap.get(keys[i]));
				//java.util.Arrays.sort(allVars[i]);
			}
			allVars= (ASVariationWithRegions[][][]) Arrays.sortNDFieldRev((Object) allVars);
			
			for (int i = 0; i < allVars.length; i++) 
				System.out.println(allVars[i][0][0]+"\t"+allVars[i].length);
			return allVars;
		}
	
	
		public static void main(String[] args) {
			parseArguments(args);
			
			ASVariationWithRegions[][][] vars= getEventsNew(new File(fName));
			path= "."+File.separator+"tst"+File.separator;
			String[] dir= new File(path).list();
			for (int i = 0; i < dir.length; i++) {
				new File(path+dir[i]).delete();
			}
			if (outputASTA) {
				try {
					PrintStream p= new PrintStream(path+"landscape.asta");
					//writeASTA(vars, p);
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			if (html) {
				writeHTML(vars, path);
			}
			
			SpliceOSigner.writeOutDomainColorMap();
			//writeLandscape(vars, baseDir);
			//writeEvents(vars, baseDir);
			//writePictures(vars, baseDir);
		}


		public static ASVariation[][] getEvents() {
				long t0= System.currentTimeMillis();
				String input= "H.sapiens.All.mRNAs.with.CDS.pfam.parsed";
				String speName= "human";
				String annoName= "RefSeq";
				String[] keywords= new String[] {"_xref"}; 
				
				System.out.print("Loading domains..");
				System.out.flush();
				DomainWrapper wrapper= new DomainWrapper(
						"douglas"+File.separator+input);
				try {
					wrapper.read();
				} catch (Exception e) {
					e.printStackTrace();
				}		
				HashMap map= wrapper.getMap(); 
				Object[] o= map.keySet().toArray();
				String[] someIDs= new String[o.length];
				int cntDomains= 0;
				DoubleVector scores= new DoubleVector();
				for (int i = 0; i < o.length; i++) { 
		//			out.print(o[i]+"\t");
					Vector v= (Vector) map.get(o[i]); 
					for (int j = 0; j < v.size(); j++) { 
		//				out.print(v.elementAt(j)+",");
						if (v.elementAt(j)== null)
							continue;
						DefaultRegion reg= (DefaultRegion) v.elementAt(j);
						boolean dontAdd= false;
						for (int k = j+1; k < v.size(); k++) {
							if (v.elementAt(k)== null)
								continue;
							DefaultRegion reg2= (DefaultRegion) v.elementAt(k);
							if (reg.overlaps(reg2)) {
								//System.err.println("Overlapping domains in "+o[i]+": "+reg+",\t"+v.elementAt(k));
								if (reg.getScore()> reg2.getScore())
									v.setElementAt(null, k);
								else {
									v.setElementAt(null, j);
									dontAdd= true;
								}
							}
						}
						if (dontAdd)
							continue;
						scores.add(reg.getScore());
						++cntDomains; 
					}
					for (int j = 0; j < v.size(); j++) {
						if (v.elementAt(j)== null)
							v.remove(j--);					
					}
		//			out.println();
					someIDs[i]= (String) o[i];
				}		
				Distribution dist= new Distribution(scores.toDoubleArray());
				double median= dist.getMedian();
				double threshold= -1d;
				
					// get basic transcripts
				String absFName= Species.getAnnotation(speName, null, annoName, keywords); 
					//Constants.getLatestUCSCAnnotation(speName, annoName, keywords);
				System.out.println("Searching for reference transcripts in "+absFName+".");
				GTFChrReader reader= new GTFChrReader(absFName);
				reader.setChromosomeWise(true);
				try {
					reader.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				Gene[] genes= reader.getGenes();
				
				HashMap<String,Vector<ASVariation>> eventMap= new HashMap<String,Vector<ASVariation>>();
				while (genes!= null) {
					for (int i = 0; i < genes.length; i++) {
						if (!genes[i].isProteinCoding())
							continue;
		
						ASVariation[] vars= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
						if (vars== null|| vars.length== 0)
							continue;
		
							// get genomic domain regions
						Transcript[] trpts= genes[i].getTranscripts();
						for (int j = 0; j < trpts.length; j++) {
							Vector v= (Vector) map.remove(trpts[j].getTranscriptID());
							Vector diregV= new Vector();
							for (int k = 0; v!= null&& k < v.size(); k++) {
								DefaultRegion regD= (DefaultRegion) v.elementAt(k);
								if (threshold< 0d&& regD.getScore()< threshold)
									continue;
								Translation tln= trpts[j].getTranslations()[0];
								DirectedRegion reg= new DirectedRegion();
								reg.setStrand(trpts[j].getStrand());
								reg.setChromosome(trpts[j].getChromosome());
								reg.setStart(tln.getGenomicPosition((regD.getStart()- 1)*3));		// TODO check for 0-based
								reg.setEnd(tln.getGenomicPosition((regD.getEnd())*3));		// -1 +1
								reg.setScore(regD.getScore());
								reg.setID(regD.getID());
								diregV.add(reg);
							}
							map.put(trpts[j].getTranscriptID(),diregV);
						}
		
						// cnt events
						Vector<ASVariation> varsDomV= new Vector<ASVariation>(); 
						Comparator domOvlCompi= new ASVariationDomain.IdentityComparator();
						for (int j = 0; vars!= null&& j < vars.length; j++) {
							DirectedRegion reg= vars[j].getRegion();
							Vector<DirectedRegion> vRegs= (Vector<DirectedRegion>) map.get(vars[j].getTranscript1().getTranscriptID());
							int k;
							for (k = 0; k < vRegs.size(); k++) {
								if (vRegs.elementAt(k).overlaps(reg))
									break;
							}
							if (k< vRegs.size()) {
								varsDomV.add(vars[j]);
								continue;
							}
							vRegs= (Vector<DirectedRegion>) map.get(vars[j].getTranscript2().getTranscriptID());
							for (k = 0; k < vRegs.size(); k++) {
								if (vRegs.elementAt(k).overlaps(reg))
									break;
							}
							if (k< vRegs.size()) 
								varsDomV.add(vars[j]);
						}
						if (varsDomV.size()== 0)
							continue;
						vars= ASVariation.removeRedundancy((ASVariation[]) Arrays.toField(varsDomV), domOvlCompi);
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
				ASVariation[][] allVars= new ASVariation[eventMap.size()][];
				Object[] keys= eventMap.keySet().toArray();
				for (int i = 0; i < keys.length; i++) 
					allVars[i]= (ASVariation[]) Arrays.toField(eventMap.get(keys[i]));
				allVars= (ASVariation[][]) Arrays.sort2DFieldRev(allVars);
				
				for (int i = 0; i < allVars.length; i++) 
					System.out.println(allVars[i][0]+"\t"+allVars[i].length);
				return allVars;
			}


		static String writeDomainStructPic(ASVariationWithRegions var) {
		  Color[] cols= SpliceOSigner.getDomainColors(var.getRegionIDs());
		  String fName= convertDomName(var.toStringRegionCode(),cols);
		  fName+= ".gif";
			try {
			      File f= new File(path+File.separator+fName);
				  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			      SpliceOSigner component= new SpliceOSigner(var);
			      component.setSize(component.getPreferredSize());
			      if (!f.exists())
			    	  t.exportToFile(f,component,component.getParent(),null,null);
			} catch (Exception e) {
				e.printStackTrace();
			}
			return fName;
		}
		
		public static String writeDotPic(Color c) {
			try {
				String imgID= Integer.toHexString(c.getRGB())+".gif";
				File f= new File(path+File.separator+imgID);
				ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
				Component component= new Circle(c,8);
			    component.setSize(component.getPreferredSize());
			      if (!f.exists())
			    	  t.exportToFile(f,component,component.getParent(),null,null);
			     return imgID;
			} catch (Exception e) {
				// TODO: handle exception
			}
			return null;
		}
		

		protected static String convertDomName(String in, Color[] cols) {
			StringBuffer sb= new StringBuffer(in);
				// kill spaces
			for (int i = 0; i < sb.length(); i++) 
				if (sb.charAt(i)== ' ')
					sb.deleteCharAt(i--);
			for (int i = 0; cols!= null&& i < cols.length; i++) 
				sb.append("C"+Integer.toHexString(cols[i].getRGB()));
			String out= sb.toString();
			
			Integer name= (Integer) fileNameMap.get(out);
			if (name== null) {
				name= new Integer(++eventTypeNr);
				fileNameMap.put(out, name);
			}
			return name.toString();
			
			
//			out= out.replace('^', 'd');
//			out= out.replace('-', 'a');
//			out= out.replace('[', 'o');
//			out= out.replace(']', 'c');
//			out= out.replace(',', 'I');			
//			
//			if (out.length()> 245) {
//				System.out.println("WARNING: shorted filename "+out);
//				out= out.substring(0, 245);	// max fname len
//			}
//			
//			return out;
		}


		public static void writeHTML(ASVariationWithRegions[][] vars, String fName) {
			
			PrintStream p;
			try {
				writePiePicture(vars, fName);
				writePictures(vars, fName);
				
				p= new PrintStream(fName+ "landscape.html");
				writeLandscapeHTML(vars, p);
				p.flush(); p.close();
				
				for (int i = 0; vars!= null&& i < vars.length; i++) {
					String varStr= vars[i][0].toString();
					String varFStr= convertFName(varStr)+ ".html";	//varStr.replace(' ', '_')+ ;
					
					p= new PrintStream(fName+ varFStr);
					StructuratorDomains.writeEventsHTML(vars[i], p);
					p.flush(); p.close();
				}
				
			} catch (Throwable e) {
				e.printStackTrace();
			}
		}


		public static void writeHTML(ASVariationWithRegions[][][] vars, String fName) {
			
			PrintStream p;
			try {
				writeHTMLlandscapePie(vars, fName);
				writePictures(vars, fName);
				
				p= new PrintStream(fName+ "landscape.html");
				writeHTMLlandscape(vars, p);
				p.flush(); p.close();
				
				for (int i = 0; vars!= null&& i < vars.length; i++) {
					String varStr= vars[i][0][0].toString();
					String varFStr= convertFName(varStr)+ ".html";	//varStr.replace(' ', '_')+ ;
					
					p= new PrintStream(fName+ varFStr);
					StructuratorDomains.writeEventsHTML(vars[i], p);
					p.flush(); p.close();
				}
				
			} catch (Throwable e) {
				e.printStackTrace();
			}
		}


								protected static void writeEventsHTML(ASVariationWithRegions[] vars, PrintStream p) {
										//		p.println("<HTML>");
										//		p.println("<HEAD>");
										//		p.println("<TITLE>Event Details</TITLE>");
										//		include(p, STYLE_FILE);
										//		p.println("</HEAD>");
										//		p.println("<BODY>");
										
												
												include(p, HEADER_FILE);
												headline(p, "Event Details", null);
												String varStr= vars[0].toString();
												String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
												
												varStr= varStr.replaceAll("\\s", "");
												p.println("<DIV class=\"userspace\">");
												p.println("<a href=\"landscape.html\"><IMG class=\"pnt\" src=\"/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
														"<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>");
												p.println("<DIV class=\"userspace\" align=\"center\">");
												p.println("<CENTER><img src=\""+ varFStr+".gif\"><br>" +
														"<FONT size=\"2\" face=\"Arial,Lucida Sans\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ varStr+ "</b></FONT></CENTER><br><br>");
												p.println("</DIV>");
												p.println("<TABLE border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
												// org=D.+melanogaster
												// org=Homo_sapiens&db=hg17
												//&clade=vertebrate&org=Mouse&db=mm8
												int ps= speStr.indexOf(";db=");	
												if (ps< 0)
													ps= speStr.length();
												String spec= speStr.substring(0, ps);
												ps= spec.indexOf("+");	
												if (ps>= 0)
													spec= spec.substring(0, ps)+ " "+ spec.substring(ps+1, spec.length());
												ps= spec.indexOf("_");	
												if (ps>= 0)
													spec= spec.substring(0, ps)+ " "+spec.substring(ps+1, spec.length());
												p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Nr</b></FONT></TD>" +
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Schema</b></FONT>" +
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Domains</b></FONT>" +
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Coordinates</b></FONT>" +
																"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Transcript ID</b> <i>Chromosome</i>: alt. Splice Sites</FONT></TD>" +
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
														"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Genome Browser</b></FONT>" +
																"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">");
												if (speStr!= null)
													p.print(spec+"</FONT></TD></TR>");
												p.println("</FONT></TD></TR>");
												for (int i = 0; i < vars.length; i++) {
													String colStr= "";
													if (i%2 == 0)
														colStr= TABLE_EVEN_COLOR;
													else
														colStr= TABLE_ODD_COLOR;
													
														// event nr
													p.println("<TR bgcolor=\""+colStr+"\">");			
													p.println("\t<TD valign=\"middle\" align=\"center\">"
															+ (i+1)+ "</TD>");
													p.println("\t<TD></TD>"); 
													
														// picture
													String domPicName= writeDomainStructPic(vars[i]);
													p.println("\t<TD valign=\"middle\" align=\"center\"><img src=\""+ domPicName+"\">"
															+ "</TD>");
													p.println("\t<TD></TD>"); 
													
														// domains
													DirectedRegion[] dom1= vars[i].getReg1();
													DirectedRegion[] dom2= vars[i].getReg2();
													HashMap<String,Vector<DirectedRegion>> domGrpHash= new HashMap<String,Vector<DirectedRegion>>();
													p.print("\t<TD valign=\"middle\" align=\"left\">");
													for (int j = 0; j < dom1.length; j++) {
														String id= ASVariationWithRegions.stripOrderSuffix(dom1[j].getID());
														Vector<DirectedRegion> v= domGrpHash.get(id);
														if (v== null)
															v= new Vector<DirectedRegion>();
														v.add(dom1[j]);
														domGrpHash.put(id, v);
													}
													for (int j = 0; j < dom2.length; j++) {	// spacer between transcripts
														String id= ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
														Vector<DirectedRegion> v= domGrpHash.get(id);
														if (v== null)
															continue;
														v.add(null);
														domGrpHash.put(id, v);
													}
													for (int j = 0; j < dom2.length; j++) {
														String id= ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
														Vector<DirectedRegion> v= domGrpHash.get(id);
														if (v== null) {
															v= new Vector<DirectedRegion>();
															v.add(null);	// spacer
														}
														v.add(dom2[j]);
														domGrpHash.put(id, v);
													}
													
													Object[] domVal= domGrpHash.keySet().toArray();
													for (int j = 0; j < domVal.length; j++) {
														String id= (String) domVal[j];
														Vector<DirectedRegion> v= (Vector<DirectedRegion>) domGrpHash.get(id);
														Color c= SpliceOSigner.getDomainColor(v.elementAt(v.size()- 1).getID());	// to be stripped, take original
														String dotFName= writeDotPic(c);
														p.print("\t<img src=\""+ dotFName+"\">&nbsp;<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+id+"</b> ");
														if (v.elementAt(0)!= null)
															p.print(vars[i].getTranscript1().getTranscriptID()+" ");
														for (int k = 0; k < v.size(); k++) {
															if (v.elementAt(k)== null) {
																p.print("<br>"+vars[i].getTranscript2().getTranscriptID()+" ");
																continue;
															}
															p.print(Formatter.fprint(v.elementAt(k).getScore(),2)+" ");
														}
														p.print("<br>");
													}
													p.println("\t</FONT></TD><TD></TD>"); 
													
													
													
														// coordinates
													String pos1Str= "";
													SpliceSite[] su= vars[i].getSpliceUniverse();
													for (int j = 0; j < vars[i].getSpliceChain1().length; j++)
														if (vars[i].getSpliceChain1()[j].getPos()> 0)
															pos1Str+= vars[i].getSpliceChain1()[j].getPos()+ ", ";
														else
															pos1Str+= Math.abs(vars[i].getSpliceChain1()[j].getPos())+ ", ";
													if (pos1Str.length()> 2)
														pos1Str= pos1Str.substring(0, pos1Str.length()- 2);
													String pos2Str= "";
													for (int j = 0; j < vars[i].getSpliceChain2().length; j++) 
														if (vars[i].getSpliceChain2()[j].getPos()> 0)
															pos2Str+= vars[i].getSpliceChain2()[j].getPos()+ ", ";
														else
															pos2Str+= Math.abs(vars[i].getSpliceChain2()[j].getPos())+ ", ";
													if (pos2Str.length()> 2)
														pos2Str= pos2Str.substring(0, pos2Str.length()- 2);
													
														// swap according to rules
													String s= vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
															+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><b>"
															+ vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
															+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</FONT></TD>";
								
								//					if (vars[i].getSpliceChain2().length== 0|| 
								//							(vars[i].getSpliceChain1().length> 0&& vars[i].getSpliceChain1()[0].getPos()< vars[i].getSpliceChain2()[0].getPos()))
								//						s= vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
								//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><b>"
								//							+ vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
								//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</FONT></TD>";
								//					else
								//						s= vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
								//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"<br><b>"
								//							+ vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
								//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"</FONT></TD>";
													
													p.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+ s);
													p.println("\t<TD></TD>"); 
										
														// linkout
													SpliceSite[] flanks= vars[i].getFlankingSpliceSites();
													int begin, end;
													if (flanks[0]== null)
														begin= Math.max(vars[i].getTranscript1().get5PrimeEdge(), vars[i].getTranscript2().get5PrimeEdge());
													else
														begin= flanks[0].getPos();
													if (flanks[1]== null)
														end= Math.min(vars[i].getTranscript1().get3PrimeEdge(), vars[i].getTranscript2().get3PrimeEdge());
													else
														end= flanks[1].getPos();
													if (begin< 0) {
														int h= -begin;
														begin= -end;
														end= h;
													}
													String ucscBase= UCSC_GB_CGI +"org="+ speStr+";";
													String refName= vars[i].getGene().getChromosome()+"_"+vars[i].getGene().getNameTranscript();
													if (annStr!= null)
														ucscBase+= "db="+ annStr+ ";";
													if (annStr.equalsIgnoreCase("hg18"))
														ucscBase+= "&hgt.customText=http://genome.imim.es/~msammeth/tstBED/"+refName+".bed;";
													// ucscBase+= UCSC_STANDARD_PAR+"knownGene=hide;mrna=pack;";
													// http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg17&hgt.customText=http://genome.imim.es/g
													String ucscEventFlanks= ucscBase+ "position="
													+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
													+ "hgFind.matches="+vars[i].getTranscript1()+","+vars[i].getTranscript2()+",";
													
													begin= su[0].getPos();
													end= su[su.length- 1].getPos();
													if (begin< 0) {
														int h= -begin;
														begin= -end;
														end= h;
													}
													String ucscEventVar= ucscBase+ "position="
													+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
													+ "hgFind.matches="+vars[i].getTranscript1()+","+vars[i].getTranscript2()+",";
													
													if (annStr!= null)
														p.println("\t<TD valign=\"middle\" align=\"center\" width=\"120\">"
																+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""+ ucscEventVar+"\">Show&nbsp;Alternative&nbsp;Parts>></a><br>" +
																	"<a href=\""+ ucscEventFlanks+"\">Show&nbsp;Complete&nbsp;Event>></a><br></FONT></TD>");
													else
														p.println("<TD valign=\"middle\" align=\"center\" width=\"120\"></TD>");	// empty
													p.println("\t<TD></TD>"); 
													
													p.println("</TR>");
												}
												p.println("</TABLE>");
												p.println("</div><!-- closing main -->");
												
												include(p, TRAILER_FILE);
											}


								protected static void writeEventsHTML(ASVariationWithRegions[][] vars,
			PrintStream p) {
		//		p.println("<HTML>");
		//		p.println("<HEAD>");
		//		p.println("<TITLE>Event Details</TITLE>");
		//		include(p, STYLE_FILE);
		//		p.println("</HEAD>");
		//		p.println("<BODY>");

		include(p, HEADER_FILE);
		headline(p, "Event Details", null);
		String varStr = vars[0][0].toString();
		String varFStr = convertFName(varStr); //varStr.replace(' ', '_');

		varStr = varStr.replaceAll("\\s", "");
		p.println("<DIV class=\"userspace\">");
		p
				.println("<a href=\"landscape.html\"><IMG class=\"pnt\" src=\"/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">"
						+ "<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>");
		p.println("<DIV class=\"userspace\" align=\"center\">");
		p
				.println("<CENTER><img src=\""
						+ varFStr
						+ ".gif\"><br>"
						+ "<FONT size=\"2\" face=\"Arial,Lucida Sans\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"
						+ varStr + "</b></FONT></CENTER><br><br>");
		p.println("</DIV>");
		p
				.println("<TABLE border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
		// org=D.+melanogaster
		// org=Homo_sapiens&db=hg17
		//&clade=vertebrate&org=Mouse&db=mm8
		int ps = speStr.indexOf(";db=");
		if (ps < 0)
			ps = speStr.length();
		String spec = speStr.substring(0, ps);
		ps = spec.indexOf("+");
		if (ps >= 0)
			spec = spec.substring(0, ps) + " "
					+ spec.substring(ps + 1, spec.length());
		ps = spec.indexOf("_");
		if (ps >= 0)
			spec = spec.substring(0, ps) + " "
					+ spec.substring(ps + 1, spec.length());
		p
				.print("<TR><TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Nr</b></FONT></TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"right\">&nbsp&nbsp</TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Schema</b></FONT>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"right\">&nbsp&nbsp</TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Domains</b></FONT>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"right\">&nbsp&nbsp</TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Coordinates</b></FONT>"
						+ "<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Transcript ID</b> <i>Chromosome</i>: alt. Splice Sites</FONT></TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"right\">&nbsp&nbsp</TD>"
						+ "<TD bgcolor=\""
						+ TABLE_HEADER_COLOR
						+ "\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Genome Browser</b></FONT>"
						+ "<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">");
		if (speStr != null)
			p.print(spec + "</FONT></TD></TR>");
		p.println("</FONT></TD></TR>");
		for (int i = 0; i < vars.length; i++) {
			String colStr = "";
			if (i % 2 == 0)
				colStr = TABLE_EVEN_COLOR;
			else
				colStr = TABLE_ODD_COLOR;

			// event nr
			p.println("<TR bgcolor=\"" + colStr + "\">");
			p.println("\t<TD valign=\"middle\" align=\"center\">" + (i + 1)
					+ "</TD>");
			p.println("\t<TD></TD>");

			// picture
			String domPicName = writeDomainStructPic(vars[i][0]);
			p.println("\t<TD valign=\"middle\" align=\"center\"><img src=\""
					+ domPicName + "\">" + "</TD>");
			p.println("\t<TD></TD>");

			// domains
			DirectedRegion[] dom1 = vars[i][0].getReg1();
			DirectedRegion[] dom2 = vars[i][0].getReg2();
			HashMap<String, Vector<DirectedRegion>> domGrpHash = new HashMap<String, Vector<DirectedRegion>>();
			p.print("\t<TD valign=\"middle\" align=\"left\">");
			for (int j = 0; j < dom1.length; j++) {
				String id = ASVariationWithRegions.stripOrderSuffix(dom1[j].getID());
				Vector<DirectedRegion> v = domGrpHash.get(id);
				if (v == null)
					v = new Vector<DirectedRegion>();
				v.add(dom1[j]);
				domGrpHash.put(id, v);
			}
			for (int j = 0; j < dom2.length; j++) { // spacer between transcripts
				String id = ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
				Vector<DirectedRegion> v = domGrpHash.get(id);
				if (v == null)
					continue;
				v.add(null);
				domGrpHash.put(id, v);
			}
			for (int j = 0; j < dom2.length; j++) {
				String id = ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
				Vector<DirectedRegion> v = domGrpHash.get(id);
				if (v == null) {
					v = new Vector<DirectedRegion>();
					v.add(null); // spacer
				}
				v.add(dom2[j]);
				domGrpHash.put(id, v);
			}

			Object[] domVal = domGrpHash.keySet().toArray();
			for (int j = 0; j < domVal.length; j++) {
				String id = (String) domVal[j];
				Vector<DirectedRegion> v = (Vector<DirectedRegion>) domGrpHash
						.get(id);
				Color c = SpliceOSigner.getDomainColor(v
						.elementAt(v.size() - 1).getID()); // to be stripped, take original
				String dotFName = writeDotPic(c);
				p
						.print("\t<img src=\""
								+ dotFName
								+ "\">&nbsp;<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"
								+ id + "</b> ");
				if (v.elementAt(0) != null)
					p
							.print(vars[i][0].getTranscript1()
									.getTranscriptID()
									+ " ");
				for (int k = 0; k < v.size(); k++) {
					if (v.elementAt(k) == null) {
						p.print("<br>"
								+ vars[i][0].getTranscript2().getTranscriptID()
								+ " ");
						continue;
					}
					p.print(Formatter.fprint(v.elementAt(k).getScore(), 2)
							+ " ");
				}
				p.print("<br>");
			}
			p.println("\t</FONT></TD><TD></TD>");

			// coordinates
			String pos1Str = "";
			SpliceSite[] su = vars[i][0].getSpliceUniverse();
			for (int j = 0; j < vars[i][0].getSpliceChain1().length; j++)
				if (vars[i][0].getSpliceChain1()[j].getPos() > 0)
					pos1Str += vars[i][0].getSpliceChain1()[j].getPos() + ", ";
				else
					pos1Str += Math.abs(vars[i][0].getSpliceChain1()[j]
							.getPos())
							+ ", ";
			if (pos1Str.length() > 2)
				pos1Str = pos1Str.substring(0, pos1Str.length() - 2);
			String pos2Str = "";
			for (int j = 0; j < vars[i][0].getSpliceChain2().length; j++)
				if (vars[i][0].getSpliceChain2()[j].getPos() > 0)
					pos2Str += vars[i][0].getSpliceChain2()[j].getPos() + ", ";
				else
					pos2Str += Math.abs(vars[i][0].getSpliceChain2()[j]
							.getPos())
							+ ", ";
			if (pos2Str.length() > 2)
				pos2Str = pos2Str.substring(0, pos2Str.length() - 2);

			// swap according to rules
			String s = vars[i][0].getTranscript1().getTranscriptID()
					+ "</b> <i>" + vars[i][0].getGene().getChromosome()
					+ "</i>: " + pos1Str + "<br><b>"
					+ vars[i][0].getTranscript2().getTranscriptID()
					+ "</b> <i>" + vars[i][0].getGene().getChromosome()
					+ "</i>: " + pos2Str + "</FONT></TD>";

			//					if (vars[i].getSpliceChain2().length== 0|| 
			//							(vars[i].getSpliceChain1().length> 0&& vars[i].getSpliceChain1()[0].getPos()< vars[i].getSpliceChain2()[0].getPos()))
			//						s= vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
			//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><b>"
			//							+ vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
			//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</FONT></TD>";
			//					else
			//						s= vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
			//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"<br><b>"
			//							+ vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
			//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"</FONT></TD>";

			p
					.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"
							+ s);
			p.println("\t<TD></TD>");

			// linkout
			SpliceSite[] flanks = vars[i][0].getFlankingSpliceSites();
			int begin, end;
			if (flanks[0] == null)
				begin = Math.max(vars[i][0].getTranscript1().get5PrimeEdge(),
						vars[i][0].getTranscript2().get5PrimeEdge());
			else
				begin = flanks[0].getPos();
			if (flanks[1] == null)
				end = Math.min(vars[i][0].getTranscript1().get3PrimeEdge(),
						vars[i][0].getTranscript2().get3PrimeEdge());
			else
				end = flanks[1].getPos();
			if (begin < 0) {
				int h = -begin;
				begin = -end;
				end = h;
			}
			String ucscBase = UCSC_GB_CGI + "org=" + speStr + ";";
			String refName = vars[i][0].getGene().getChromosome() + "_"
					+ vars[i][0].getGene().getNameTranscript();
			if (annStr != null)
				ucscBase += "db=" + annStr + ";";
			if (annStr.equalsIgnoreCase("hg18"))
				ucscBase += "&hgt.customText=http://genome.imim.es/~msammeth/tstBED/"
						+ refName + ".bed;";
			// ucscBase+= UCSC_STANDARD_PAR+"knownGene=hide;mrna=pack;";
			// http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg17&hgt.customText=http://genome.imim.es/g
			String ucscEventFlanks = ucscBase + "position="
					+ vars[i][0].getGene().getChromosome() + ":"
					+ (begin - UCSC_FLANK) + "-" + (end + UCSC_FLANK) + ";"
					+ "hgFind.matches=" + vars[i][0].getTranscript1() + ","
					+ vars[i][0].getTranscript2() + ",";

			begin = su[0].getPos();
			end = su[su.length - 1].getPos();
			if (begin < 0) {
				int h = -begin;
				begin = -end;
				end = h;
			}
			String ucscEventVar = ucscBase + "position="
					+ vars[i][0].getGene().getChromosome() + ":"
					+ (begin - UCSC_FLANK) + "-" + (end + UCSC_FLANK) + ";"
					+ "hgFind.matches=" + vars[i][0].getTranscript1() + ","
					+ vars[i][0].getTranscript2() + ",";

			if (annStr != null)
				p
						.println("\t<TD valign=\"middle\" align=\"center\" width=\"120\">"
								+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""
								+ ucscEventVar
								+ "\">Show&nbsp;Alternative&nbsp;Parts>></a><br>"
								+ "<a href=\""
								+ ucscEventFlanks
								+ "\">Show&nbsp;Complete&nbsp;Event>></a><br></FONT></TD>");
			else
				p
						.println("<TD valign=\"middle\" align=\"center\" width=\"120\"></TD>"); // empty
			p.println("\t<TD></TD>");

			p.println("</TR>");
		}
		p.println("</TABLE>");
		p.println("</div><!-- closing main -->");

		include(p, TRAILER_FILE);
	}

}
