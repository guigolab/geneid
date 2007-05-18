package gphase;


import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.text.Collator;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import org.apache.commons.collections.BidiMap;
import org.apache.commons.collections.bidimap.DualHashBidiMap;

import prefuse.render.NullRenderer;

import gphase.Structurator;
import gphase.algo.ASAnalyzer;
import gphase.db.MapTable;
import gphase.ext.DevNullReaderThread;
import gphase.ext.DevPipeReaderThread;
import gphase.io.TabDelimitedFormatWrapper;
import gphase.io.gtf.Andre;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.io.gtf.ProgressiveIOWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.tools.Distribution;
import gphase.tools.Formatter;
import gphase.tools.IntVector;


public class ProteinAnalyzer {
	static String[] species= null;
	static int[] testAnnot= null;
	static int[] icount= null;
	
	public static void main(String[] args) {
		_00_analyzeDouglas();
		//_00_checkHubNoHubOverlap();
		//_01_analyzeDouglasOverlapGenewise();
		if (1== 1)
			System.exit(0);
		
		
		
		species= new String[] {"fruitfly"};
		//icount= new int[] {12};
		testAnnot= new int[] {1};
		
		String mName= "_01_analyzeASEvents2";
		Class[] sig= new Class[] {Gene.class, PrintStream.class, HashMap.class};
		Method[] m= null;
		try {
			m= new Method[] {ProteinAnalyzer.class.getMethod(mName, sig)};
		} catch (Exception e) {
			e.printStackTrace();
		} 
		if (m== null) 
			System.err.println("Method "+mName+" not found.");
		else
			_00_mainLoopChromosome(m);
	}
	static String getDouglasAbbreviation(String someName) {
		String binName= Species.getBinomialForSomeName(someName);
		int p= binName.indexOf("_");
		String dougName= binName.substring(0,1).toUpperCase()+
					binName.substring(p+1, p+5).toLowerCase();	// eg, Mmusc
		return dougName;
	}
	
	static void _00_analyzeDouglas() {
		final String[] DOUGLAS_SPECIES= new String[] {
			"human", "mouse", "fruitfly", "worm"
		};
		
		String subdir= "douglas";
		String[] list= new File(subdir).list();
		for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
//			if (!DOUGLAS_SPECIES[i].equals("worm"))
//				continue;
			String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
			for (int j = 0; j < list.length; ++j) { 
				if (!list[j].startsWith(speName)|| new File(subdir+ File.separator+ list[j]).isDirectory())
					continue;
				
				System.out.println(list[j]);
				TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
						subdir+ File.separator+ list[j]);
				try {
					wrapper.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				String[] ids= wrapper.getColumn(0);
				System.out.println("Read "+ids.length+" prot ids.");
				
				File dir= new File(subdir+File.separator+list[j]+"_landscape");
				if (dir.exists())
					dir.delete();
				dir.mkdir();
				
				String comSpeName= DOUGLAS_SPECIES[i];
				String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
				String[] notFoundIDs= ids;	//nonRefSeqIDs;
				Vector v= new Vector();
				for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
							k < specAnnos.length; k++) {
					String[] specIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], specAnnos[k], notFoundIDs);
					notFoundIDs= MapTable.getNotFoundIDs();
					
					String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
					EncodeWrapper gtf= new EncodeWrapper(fName);
					try {
						PrintStream err= System.err;
						System.setErr(new PrintStream("delme.tmp"));
						gtf.read(specIDs);
						new File("delme.tmp").delete();
						System.setErr(err);
					} catch (Exception e) {
						e.printStackTrace();
					}
					GTFObject[] obj= gtf.getGtfObj();
					v= (Vector) gphase.tools.Arrays.addAll(v, obj);
					int x= 0;
					if (specIDs!= null)
						x= specIDs.length;
					int y= 0;
					if (obj!= null)
						y= obj.length;
					System.out.println("Retrieved "+x+" new gene ids, "+y+" new gtf objects. Total "+v.size());
				}
				
					// standard annot, eg REFSEQ, not many AS in Droso..
				String[] RefSeqIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], "RefSeq", notFoundIDs);
				String fName= Constants.getLatestUCSCAnnotation(comSpeName, "RefSeq", null);
				EncodeWrapper gtf= new EncodeWrapper(fName);
				try {
					PrintStream err= System.err;
					System.setErr(new PrintStream("delme.tmp"));
					gtf.read(RefSeqIDs);
					new File("delme.tmp").delete();
					System.setErr(err);
				} catch (Exception e) {
					e.printStackTrace();
				}
				GTFObject[] obj= gtf.getGtfObj();
				v= (Vector) gphase.tools.Arrays.addAll(v, gtf.getGtfObj());
				int x= 0;
				if (RefSeqIDs!= null)
					x= RefSeqIDs.length;
				int y= 0;
				if (obj!= null)
					y= obj.length;
				System.out.println("Retrieved "+x+" gene ids, "+y+" new gtf objects. Total "+v.size());
				
	
				
				
//				GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
//				System.out.println("analyzing "+objs.length+" gtf objects..");
//				Graph g= EncodeWrapper.assemble(false, objs);
//				
//				ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
//				gphase.tools.Arrays.sort2DFieldRev(vars);
//				Structurator.writeHTML(vars, dir+File.separator+list[j]);
			}
			
		}
	}
	static void _00_analyzeGO() {
		final String[] species= new String[] {
			"human", "mouse", "rat", "fruitfly", "worm"
		};
		
		String subdir= "go";
		String[] list= new File(subdir).list();
		for (int i = 0; i < species.length; i++) {
			if (!species[i].equals("worm"))
				continue;
			System.out.println(species[i]);
			
			String binName= Species.getBinomialForCommonName(species[i]);	// "Homo Sapiens"
			String abbName= Species.getAbbrevNameForBinomial(binName);	// H_sapiens	
			abbName= Character.toUpperCase(abbName.charAt(0))+"_"+abbName.substring(1);
			for (int j = 0; j < list.length; ++j) { 
				if (!list[j].startsWith(abbName)|| new File(subdir+ File.separator+ list[j]).isDirectory())
					continue;
				
				System.out.println(list[j]);
				TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
						subdir+ File.separator+ list[j]);
				try {
					wrapper.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				String[] ids= wrapper.getColumn(1);	// gene symbol
				System.out.println("Read "+ids.length+" prot ids.");
				
				File dir= new File(subdir+File.separator+list[j]+"_landscape");
				if (dir.exists())
					dir.delete();
				dir.mkdir();
				
				String comSpeName= species[i];
				String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
				String[] notFoundIDs= ids;	//nonRefSeqIDs;
				Vector v= new Vector();
				for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
							k < specAnnos.length; k++) {
					String[] specIDs= MapTable.getPrimaryGeneID(species[i], specAnnos[k], notFoundIDs);
					notFoundIDs= MapTable.getNotFoundIDs();
					
					String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
					EncodeWrapper gtf= new EncodeWrapper(fName);
					try {
						PrintStream err= System.err;
						System.setErr(new PrintStream("delme.tmp"));
						gtf.read(specIDs);
						new File("delme.tmp").delete();
						System.setErr(err);
					} catch (Exception e) {
						e.printStackTrace();
					}
					GTFObject[] obj= gtf.getGtfObj();
					v= (Vector) gphase.tools.Arrays.addAll(v, obj);
					int x= 0;
					if (specIDs!= null)
						x= specIDs.length;
					int y= 0;
					if (obj!= null)
						y= obj.length;
					System.out.println("Retrieved "+x+" new gene ids, "+y+" new gtf objects. Total "+v.size());
				}
				
					// standard annot, eg REFSEQ, not many AS in Droso..
				String[] RefSeqIDs= MapTable.getPrimaryGeneID(species[i], "RefSeq", notFoundIDs);
				String fName= Constants.getLatestUCSCAnnotation(comSpeName, "RefSeq");
				EncodeWrapper gtf= new EncodeWrapper(fName);
				try {
					PrintStream err= System.err;
					System.setErr(new PrintStream("delme.tmp"));
					gtf.read(RefSeqIDs);
					new File("delme.tmp").delete();
					System.setErr(err);
				} catch (Exception e) {
					e.printStackTrace();
				}
				GTFObject[] obj= gtf.getGtfObj();
				v= (Vector) gphase.tools.Arrays.addAll(v, gtf.getGtfObj());
				int x= 0;
				if (RefSeqIDs!= null)
					x= RefSeqIDs.length;
				int y= 0;
				if (obj!= null)
					y= obj.length;
				System.out.println("Retrieved "+x+" gene ids, "+y+" new gtf objects. Total "+v.size());
				
				GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
				System.out.println("analyzing "+objs.length+" gtf objects..");
				Graph g= EncodeWrapper.assemble(false, objs);
				
				ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
				gphase.tools.Arrays.sort2DFieldRev(vars);
				PrintStream p= new PrintStream(dir+File.separator+"landscape");
				ASAnalyzer.outputVariations(vars, false, false, p);				
				Structurator.writeHTML(vars, dir+File.separator+list[j]);
			}
			
		}
	}
	static void _00_analyzeDouglasOverlap() {
			final String[] DOUGLAS_SPECIES= new String[] {
				"human", "mouse", "fruitfly", "worm"
			};
			
			String subdir= "douglas";
			String[] list= new File(subdir).list();
			for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
//				if (!DOUGLAS_SPECIES[i].equals("mouse"))
//					continue;
				String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
				for (int j = 0; j < list.length; ++j) { 
					if (!list[j].startsWith(speName)|| new File(subdir+ File.separator+ list[j]).isDirectory())
						continue;
					
					System.out.println(list[j]);
					TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
							subdir+ File.separator+ list[j]);
					try {
						wrapper.read();
					} catch (Exception e) {
						e.printStackTrace();
					}
					String[] ids= wrapper.getColumn(0);
					System.out.println("Read "+ids.length+" prot ids.");
					for (int k = 0; k < ids.length; k++) {	// trim version numbers
						int p=  ids[k].lastIndexOf('.');
						if (p>= 0)
							ids[k]= ids[k].substring(0, p);
					}
					
					
					String comSpeName= DOUGLAS_SPECIES[i];
					String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
					String[] notFoundIDs= ids;	//nonRefSeqIDs;
					Vector v= new Vector();
					for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
								k < specAnnos.length; k++) {
						String[] specIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], specAnnos[k], notFoundIDs);
						notFoundIDs= MapTable.getNotFoundIDs();
						
						String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
						EncodeWrapper gtf= new EncodeWrapper(fName);
						try {
	//						PrintStream err= System.err;
	//						System.setErr(new PrintStream("delme.tmp"));
							gtf.read(specIDs);
	//						new File("delme.tmp").delete();
	//						System.setErr(err);
						} catch (Exception e) {
							e.printStackTrace();
						}
						GTFObject[] obj= gtf.getGtfObj();
						v= (Vector) gphase.tools.Arrays.addAll(v, obj);
						int x= 0;
						if (specIDs!= null)
							x= specIDs.length;
						int y= 0;
						if (obj!= null)
							y= obj.length;
						System.out.println("Retrieved "+x+" new gene ids, "+y+" new gtf objects. Total "+v.size());
					}
					if (notFoundIDs!= null&& notFoundIDs.length> 0) {
						System.out.print("Not found: ");
						for (int k = 0; k < notFoundIDs.length; k++) 
							System.out.print(notFoundIDs[k]+" ");
						System.out.println();
					}
	
						// get reference gene set
					GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
					System.out.println("analyzing "+objs.length+" gtf objects..");
					Gene[] ge= ProgressiveIOWrapper.assemble(objs);
					HashMap chrHash= new HashMap();
					for (int k = 0; k < ge.length; k++) {
						Vector geV= (Vector) chrHash.get(ge[k].getChromosome());
						if (geV== null) {
							geV= new Vector();
							chrHash.put(ge[k].getChromosome(), geV);
						}
						geV.add(ge[k]);
					}
					
						// overlap with reference annotations
					for (int k = 0; specAnnos!= null&& k < specAnnos.length; k++) {
						ProgressiveIOWrapper gtf= new ProgressiveIOWrapper(Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null));
						try {
							gtf.read();
						} catch (Exception e) {
							e.printStackTrace();
						}
						Species sp= gtf.assemble();
						ge= sp.getGenes();
						System.out.println("Analyzing reference annotation "+specAnnos[k]+", "+ge.length+" genes.");
						
						File dir= new File(subdir+File.separator+comSpeName+"_"+specAnnos[k]+"_landscape");
						if (dir.exists())
							dir.delete();
						dir.mkdir();
	
						System.out.print("retrieving variations..");
						System.out.flush();
						ASVariation[][] vars= sp.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
						gphase.tools.Arrays.sort2DFieldRev(vars);
						Structurator.writeHTML(vars, dir+File.separator+list[j]);
						System.out.println("done");
	
						Vector geneV= new Vector(ge.length/ 2);
						for (int m = 0; m < ge.length; m++) {
							Vector geV= (Vector) chrHash.get(ge[m].getChromosome());
							if (geV== null) {
								sp.remove(ge[m], true);
								continue;
							}
							int n= 0;
							for (n = 0; n < geV.size(); n++) 
								if (ge[m].overlaps((Gene) geV.elementAt(n)))
									break;
							if (n== geV.size())
								sp.remove(ge[m], true);
						}
						ge= sp.getGenes();
						System.out.println("Found "+ge.length+" genes.");
						
						dir= new File(subdir+File.separator+list[j]+"_"+specAnnos[k]+"_landscape");
						if (dir.exists())
							dir.delete();
						dir.mkdir();
						
						vars= sp.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
						gphase.tools.Arrays.sort2DFieldRev(vars);
						Structurator.writeHTML(vars, dir+File.separator+list[j]);
					}
					
				}
				
			}
		}
	static void _01_analyzeDouglasOverlapGenewise() {
			final String[] DOUGLAS_SPECIES= new String[] {
				"human", "mouse", "fruitfly", "worm"
			};
			
			String subdir= "douglas";
			String[] list= new File(subdir).list();
			for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
	//			if (!DOUGLAS_SPECIES[i].equals("mouse"))
	//				continue;
				String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
				for (int j = 0; j < list.length; ++j) { 
					if (!(list[j].startsWith(speName)&& list[j].contains(".interaction.count"))|| new File(subdir+ File.separator+ list[j]).isDirectory())
						continue;
					
					System.out.println(list[j]);
					TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
							subdir+ File.separator+ list[j]);
					try {
						wrapper.read();
					} catch (Exception e) {
						e.printStackTrace();
					}
					String[] all_ids= wrapper.getColumn(0);
					String[] strActions= wrapper.getColumn(1);
					HashMap map= new HashMap();	//# interact - V(protIDs)
					for (int k = 0; k < strActions.length; k++) {
						String key= strActions[k].trim();
						Vector v= (Vector) map.get(key);
						if (v== null) {
							v= new Vector();
							map.put(key, v);
						}
						v.add(all_ids[k]);
					}
					
					System.out.println("Read "+all_ids.length+" prot ids.");
					for (int k = 0; k < all_ids.length; k++) {	// trim version numbers
						int p=  all_ids[k].lastIndexOf('.');
						if (p>= 0)
							all_ids[k]= all_ids[k].substring(0, p);
					}
					
						// iterate all sets
					Object[] keys= map.keySet().toArray();
					Arrays.sort(keys);
					for (int x = keys.length- 1; x >= 0; --x) {
						if (!keys[x].equals("11"))
							continue;		// REMOVE !!!
						String comSpeName= DOUGLAS_SPECIES[i];
						String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
						
						String[] notFoundIDs= (String[]) gphase.tools.Arrays.toField((Vector) map.get(keys[x]));	//nonRefSeqIDs;
						System.out.println("Proteins with "+keys[x]+" interactions: "+notFoundIDs.length+".");
						Vector v= new Vector();
						for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
									k < specAnnos.length; k++) {
							String[] specIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], specAnnos[k], notFoundIDs);
							notFoundIDs= MapTable.getNotFoundIDs();
							
							String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
							EncodeWrapper gtf= new EncodeWrapper(fName);
							try {
	//							PrintStream err= System.err;
	//							System.setErr(new PrintStream("delme.tmp"));
								gtf.read(specIDs);
	//							new File("delme.tmp").delete();
	//							System.setErr(err);
							} catch (Exception e) {
								e.printStackTrace();
							}
							GTFObject[] obj= gtf.getGtfObj();
							v= (Vector) gphase.tools.Arrays.addAll(v, obj);
							int z= 0;
							if (specIDs!= null)
								z= specIDs.length;
							int y= 0;
							if (obj!= null)
								y= obj.length;
							System.out.println("Retrieved "+z+" new gene ids, "+y+" new gtf objects. Total "+v.size());
						}
						if (notFoundIDs!= null&& notFoundIDs.length> 0) {
							System.out.print("Not found: ");
							for (int k = 0; k < notFoundIDs.length; k++) 
								System.out.print(notFoundIDs[k]+" ");
							System.out.println();
						}
	
							// get reference gene set
						GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
						if (objs== null) {
							System.out.println(keys[x]+" INTERACTIONS: skipped, no gtf-obj.");
							continue;
						}
						System.out.println("analyzing "+objs.length+" gtf objects..");
						Gene[] ge= ProgressiveIOWrapper.assemble(objs);
						if (ge== null) {
							System.out.println(keys[x]+" INTERACTIONS: skipped, no genes.");
							continue;
						}
						HashMap chrHash= new HashMap();
						for (int k = 0; k < ge.length; k++) {
							Vector geV= (Vector) chrHash.get(ge[k].getChromosome());
							if (geV== null) {
								geV= new Vector();
								chrHash.put(ge[k].getChromosome(), geV);
							}
							geV.add(ge[k]);
						}
						
							// overlap with reference annotations
						for (int k = 0; specAnnos!= null&& k < specAnnos.length; k++) {
							if (!specAnnos[k].contains("MGC"))
								continue;						// REMOVE !!!
							ProgressiveIOWrapper gtf= new ProgressiveIOWrapper(Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k]));
							try {
								gtf.read();
							} catch (Exception e) {
								e.printStackTrace();
							}
							Species sp= gtf.assemble();
							ge= sp.getGenes();
							System.out.println("Analyzing reference annotation "+specAnnos[k]+", "+ge.length+" genes.");
							
							Vector geneV= new Vector(ge.length/ 2);
							for (int m = 0; m < ge.length; m++) {
								Vector geV= (Vector) chrHash.get(ge[m].getChromosome());
								if (geV== null) {
									sp.remove(ge[m], true);
									continue;
								}
								int n= 0;
								for (n = 0; n < geV.size(); n++) 
									if (ge[m].overlaps((Gene) geV.elementAt(n)))
										break;
								if (n== geV.size())
									sp.remove(ge[m], true);
							}
							ge= sp.getGenes();
							System.out.println("Found "+ge.length+" genes.");
							
							ASVariation[][] vars= sp.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
							IntVector degreeV= new IntVector();
							geneV= new Vector();
							for (int m = 0; vars!= null&& m < vars.length; m++) {		
								int deg= vars[m][0].getDegree();
								for (int n = 0; n < vars[m].length; n++) { 
									degreeV.add(deg);
									geneV= (Vector) gphase.tools.Arrays.addUnique(geneV, vars[m][n].getTranscript1().getGene());
								}
							}
	
							System.out.println(keys[x]+" INTERACTIONS");
	
								// degree
							Distribution dist= new Distribution(degreeV.toIntArray());
							System.out.println("DEGREE: med "+dist.getMedian()+", mea "+dist.getMean()+", std "+dist.getStandardDeviation());
							gphase.tools.Arrays.sort2DFieldRev(vars);
							
								// transcripts
							IntVector nbTrans= new IntVector();
							IntVector nbExons= new IntVector();
							for (int m = 0; m < geneV.size(); m++) {
								Gene gen= (Gene) geneV.elementAt(m);
								nbTrans.add(gen.getTranscriptCount());
								nbExons.add(gen.getExons().length);
							}
							dist= new Distribution(nbTrans.toIntArray());
							System.out.println("#TRPT  med "+dist.getMedian()+", mea "+dist.getMean()+", std "+dist.getStandardDeviation());
							dist= new Distribution(nbExons.toIntArray());
							System.out.println("#EXON  med "+dist.getMedian()+", mea "+dist.getMean()+", std "+dist.getStandardDeviation());
							
							//Structurator.writeHTML(vars, dir+File.separator+list[j]);
						}
						
					}
					
				}
				
			}
		}
	static void _00_mainLoopChromosome(Method[] target) {
			final String[] DOUGLAS_SPECIES= new String[] {
				"human", "mouse", "fruitfly", "worm"};
			String subdir= "douglas";
			String[] list= new File(subdir).list();
			
			for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
				
					// start at species
				if (species!= null) {
					int j;
					for (j = 0; j < species.length; j++) 
						if (species[j].equalsIgnoreCase(DOUGLAS_SPECIES[i]))
							break;
					if (j== species.length)
						continue;
				}
				
				String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
				for (int j = 0; j < list.length; ++j) { 
					if (!(list[j].startsWith(speName)&& list[j].contains(".interaction.count"))|| new File(subdir+ File.separator+ list[j]).isDirectory())
						continue;
					
					System.out.println("\n\n=== "+list[j]+" ===");
					TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
							subdir+ File.separator+ list[j]);
					try {
						wrapper.read();
					} catch (Exception e) {
						e.printStackTrace();
					}
					String[] all_ids= wrapper.getColumn(0);
					String[] strActions= wrapper.getColumn(1);
					HashMap map= new HashMap();	//# interact - V(protIDs)
					for (int k = 0; k < strActions.length; k++) {
						String key= strActions[k].trim();
						Vector v= (Vector) map.get(key);
						if (v== null) {
							v= new Vector();
							map.put(key, v);
						}
						v.add(all_ids[k]);
					}
					
					System.out.println("Read "+all_ids.length+" prot ids.\n\n");
					for (int k = 0; k < all_ids.length; k++) {	// trim version numbers
						int p=  all_ids[k].lastIndexOf('.');
						if (p>= 0)
							all_ids[k]= all_ids[k].substring(0, p);
					}
					
						// iterate all protein w x connections
					Object[] keys= map.keySet().toArray();
					int[] ikeys= new int[keys.length];
					for (int k = 0; k < ikeys.length; k++) 
						ikeys[k]= Integer.parseInt((String) keys[k]);
					Arrays.sort(ikeys);
					for (int x = 0; x < ikeys.length; ++x) {
//						if (!keys[x].equals("11"))
//							continue;		// REMOVE !!!
						String comSpeName= DOUGLAS_SPECIES[i];
						String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
						

							// (2) GET GTF-OBJECT for GENE-IDs
						String[] notFoundIDs= (String[]) gphase.tools.Arrays.toField((Vector) map.get(Integer.toString(ikeys[x])));	//nonRefSeqIDs;
							// check interaction count
						if (icount!= null) {
							int k;
							for (k = 0; k < icount.length; k++) {
								if (ikeys[x] == icount[i])
									break;
							}
							if (k== icount.length)
								continue;
						}
						System.out.println("\n=> Proteins with "+ikeys[x]+" interactions: "+notFoundIDs.length+".\n");
						Vector v= new Vector();
						GTFChrReader gtf= null;
						for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
									k < specAnnos.length; k++) {
							if (k> 1)	// !!! REMOVE !!!
								continue;
							String[] specIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], specAnnos[k], notFoundIDs);
							notFoundIDs= MapTable.getNotFoundIDs();
							
							String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
							gtf= new GTFChrReader(fName);	// filters DNA automatically
							//gtf.setFiltTrptIDs(specIDs);
							try {
	//							PrintStream err= System.err;
	//							System.setErr(new PrintStream("delme.tmp"));
								if (specIDs!= null)
									gtf.read(specIDs);	// geneIDs, because primaryGeneID, may contain (rarely) >1 transcript
	//							new File("delme.tmp").delete();
	//							System.setErr(err);
							} catch (Exception e) {
								e.printStackTrace();
							}
							GTFObject[] obj= gtf.getGtfObj();
							v= (Vector) gphase.tools.Arrays.addAll(v, obj); 
							int z= 0;
							if (specIDs!= null)
								z= specIDs.length;
							int y= 0;
							if (obj!= null)
								y= obj.length;
							System.out.println("\tRetrieved "+z+" new gene ids, "+y+" new gtf objects. Total "+v.size()+"\n");
						}
						
						
						if (notFoundIDs!= null&& notFoundIDs.length> 0) {
							System.out.print("\tNot found: ");
							for (int k = 0; k < notFoundIDs.length; k++) 
								System.out.print(notFoundIDs[k]+" ");
							System.out.println();
						}
	
							// get reference gene set
						GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
						if (objs== null) {
							System.out.println(ikeys[x]+" INTERACTIONS: skipped, no gtf-obj.");
							continue;
						}
						System.out.print("\tanalyzing "+objs.length+" gtf objects..");
						Gene[] ge= gtf.assemble(objs);	//ProgressiveIOWrapper.assemble(objs);
						System.out.println(" found "+ge.length+" genes.");
						if (ge== null) {
							System.out.println(ikeys[x]+" INTERACTIONS: skipped, no genes.");
							continue;
						}
						HashMap localMap= new HashMap();
						for (int q = 0; q < ge.length; q++) {	// chromosome genes
							for (int m = 0; m < target.length; m++) {	// methods
								try {
									target[m].invoke(null, new Object[]{ge[q], null, localMap});
								} catch (Exception e) {
									e.printStackTrace();
								}
							} 
						}
						
							// output
						System.out.println("\ntranscripts encoding protein only");
						for (int m = 0; m < target.length; m++) {	// methods
							try {
								File f= new File(subdir+File.separator+
										DOUGLAS_SPECIES[i]+"_"+"prot-only_"+target[m].getName());
								f.createNewFile();
								PrintStream p= new PrintStream(f);
								target[m].invoke(null, new Object[]{null, p, localMap});
								p.flush();
								p.close();
							} catch (Exception e) {
								e.printStackTrace();
							}
						}

						HashMap chrHash= new HashMap();
						for (int k = 0; k < ge.length; k++) {
							Vector geV= (Vector) chrHash.get(ge[k].getChromosome());
							if (geV== null) {
								geV= new Vector();
								chrHash.put(ge[k].getChromosome(), geV);
							}
							geV.add(ge[k]);
						}
						
						
						
							// GET GENES, OVERLAP 
						Object[] o= chrHash.keySet().toArray();
						String[] chrIDs= new String[o.length];
						for (int k = 0; k < chrIDs.length; k++) 
							chrIDs[k]= (String) o[k];
						Arrays.sort(o);
						
						
							// 3. overlap with reference annotations
							// ...
						System.out.println("\ncomplete clusters");
						for (int k = 0; o.length> 0&& specAnnos!= null&& k < specAnnos.length; k++) {
							if (testAnnot!= null) {
								int m;
								for (m = 0; m < testAnnot.length; m++) {
									if (testAnnot[m]== k)
										break;
								}
								if (m== testAnnot.length)
									continue;
							} else if (k> 2)
								break;						// REMOVE !!!
							System.out.print(specAnnos[k]+" ");
							System.out.flush();
							int chrPtr= 0;
							gtf= new GTFChrReader(
									Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], new String[] {".gtf_xref"}));
							gtf.setSilent(true);
//							gtf.skipToChromosome((String) o[chrPtr++]);
							Vector w= new Vector();
//							while (!gtf.isEOF()) {
//								Gene[] gene= gtf.readNextChromosome(ge);
//								w= (Vector) gphase.tools.Arrays.addAll(v, gene);
//							}
							
							gtf.sweepToChromosome((String) o[chrPtr++]);
							gtf.setSilent(true);
							gtf.setFiltChrIDs(chrIDs);
							gtf.setChromosomeWise(true);
							Gene[] gene= null;
							try {
								gtf.read();
							} catch (Exception e1) {
								e1.printStackTrace();
							}
							Vector geV= null;
							if (gene!= null&& gene.length> 0)
								geV= (Vector) chrHash.get(gene[0].getChromosome());
							
							while (gene!= null) {
								for (int q = 0; geV!= null&& q < gene.length; q++) {
									int n= 0;
									for (n = 0; n < geV.size(); n++) 
										if (gene[q].overlaps((Gene) geV.elementAt(n)))
											break;
									if (n< geV.size())  
										w.add(gene[q]);
								}
								
								if (chrPtr== o.length)
									break;	// all chr analyzed
								
								gtf.sweepToChromosome((String) o[chrPtr++]);
								try {
									gtf.read();
								} catch (Exception e) {
									e.printStackTrace();
								}	// next
								gene= gtf.getGenes();
								if (gene!= null&& gene.length> 0)
									geV= (Vector) chrHash.get(gene[0].getChromosome());
							}
							System.out.println("Found "+w.size()+" genes.");

								// analyze target
							localMap= new HashMap();
							for (int q = 0; q < w.size(); q++) {	// chromosome genes
								for (int m = 0; m < target.length; m++) {	// methods
									try {
										target[m].invoke(null, new Object[]{w.elementAt(q), null, localMap});
									} catch (Exception e) {
										e.printStackTrace();
									}
								} 
							}
							
								// output
							for (int m = 0; m < target.length; m++) {	// methods
								try {
									File f= new File(subdir+File.separator+
											DOUGLAS_SPECIES[i]+"_"+specAnnos[k]+"_"+"cluster_"+target[m].getName());
									f.createNewFile();
									PrintStream p= new PrintStream(f);
									target[m].invoke(null, new Object[]{null, p, localMap});
									p.flush();
									p.close();
								} catch (Exception e) {
									e.printStackTrace();
								}
							}


								
						}	// end ref annotations
					}	// end interaction#
				}	// files (input)
			}	// species
				
		}
	

	static void _00_mainLoopSpecies(Method[] target) {
		final String[] DOUGLAS_SPECIES= new String[] {
			"human", "mouse", "fruitfly", "worm"};
		String subdir= "douglas";
		String[] list= new File(subdir).list();
		
		for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
//			if (!DOUGLAS_SPECIES[i].equals("mouse"))
//				continue;
			String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
			for (int j = 0; j < list.length; ++j) { 
				if (!(list[j].startsWith(speName)&& list[j].contains(".interaction.count"))|| new File(subdir+ File.separator+ list[j]).isDirectory())
					continue;
				
				System.out.println(list[j]);
				TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(
						subdir+ File.separator+ list[j]);
				try {
					wrapper.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				String[] all_ids= wrapper.getColumn(0);
				String[] strActions= wrapper.getColumn(1);
				HashMap map= new HashMap();	//# interact - V(protIDs)
				for (int k = 0; k < strActions.length; k++) {
					String key= strActions[k].trim();
					Vector v= (Vector) map.get(key);
					if (v== null) {
						v= new Vector();
						map.put(key, v);
					}
					v.add(all_ids[k]);
				}
				
				System.out.println("Read "+all_ids.length+" prot ids.");
				for (int k = 0; k < all_ids.length; k++) {	// trim version numbers
					int p=  all_ids[k].lastIndexOf('.');
					if (p>= 0)
						all_ids[k]= all_ids[k].substring(0, p);
				}
				
					// iterate all sets
				Object[] keys= map.keySet().toArray();
				Arrays.sort(keys);
				for (int x = keys.length- 1; x >= 0; --x) {
//					if (!keys[x].equals("11"))
//						continue;		// REMOVE !!!
					String comSpeName= DOUGLAS_SPECIES[i];
					String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(comSpeName)];
					
					String[] notFoundIDs= (String[]) gphase.tools.Arrays.toField((Vector) map.get(keys[x]));	//nonRefSeqIDs;
					System.out.println("Proteins with "+keys[x]+" interactions: "+notFoundIDs.length+".");
					Vector v= new Vector();
					for (int k = 0; specAnnos!= null&& notFoundIDs!= null&& notFoundIDs.length> 0&&
								k < specAnnos.length; k++) {
						String[] specIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], specAnnos[k], notFoundIDs);
						notFoundIDs= MapTable.getNotFoundIDs();
						
						String fName= Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k], null);
						EncodeWrapper gtf= new EncodeWrapper(fName);
						try {
//							PrintStream err= System.err;
//							System.setErr(new PrintStream("delme.tmp"));
							gtf.read(specIDs);
//							new File("delme.tmp").delete();
//							System.setErr(err);
						} catch (Exception e) {
							e.printStackTrace();
						}
						GTFObject[] obj= gtf.getGtfObj();
						v= (Vector) gphase.tools.Arrays.addAll(v, obj);
						int z= 0;
						if (specIDs!= null)
							z= specIDs.length;
						int y= 0;
						if (obj!= null)
							y= obj.length;
						System.out.println("Retrieved "+z+" new gene ids, "+y+" new gtf objects. Total "+v.size());
					}
					
					
					if (notFoundIDs!= null&& notFoundIDs.length> 0) {
						System.out.print("Not found: ");
						for (int k = 0; k < notFoundIDs.length; k++) 
							System.out.print(notFoundIDs[k]+" ");
						System.out.println();
					}

						// get reference gene set
					GTFObject[] objs= (GTFObject[]) gphase.tools.Arrays.toField(v);
					if (objs== null) {
						System.out.println(keys[x]+" INTERACTIONS: skipped, no gtf-obj.");
						continue;
					}
					System.out.println("analyzing "+objs.length+" gtf objects..");
					Gene[] ge= ProgressiveIOWrapper.assemble(objs);
					if (ge== null) {
						System.out.println(keys[x]+" INTERACTIONS: skipped, no genes.");
						continue;
					}
					HashMap chrHash= new HashMap();
					for (int k = 0; k < ge.length; k++) {
						Vector geV= (Vector) chrHash.get(ge[k].getChromosome());
						if (geV== null) {
							geV= new Vector();
							chrHash.put(ge[k].getChromosome(), geV);
						}
						geV.add(ge[k]);
					}
					
						// overlap with reference annotations
					for (int k = 0; specAnnos!= null&& k < specAnnos.length; k++) {
//						if (!specAnnos[k].contains("MGC"))
//							continue;						// REMOVE !!!
						
						EncodeWrapper gtf= new EncodeWrapper(Constants.getLatestUCSCAnnotation(comSpeName, specAnnos[k]));
						Species sp= gtf.getGraph(false).getSpecies()[0];
						
						ge= sp.getGenes();
						System.out.println("Analyzing reference annotation "+specAnnos[k]+", "+ge.length+" genes.");
						
						Vector geneV= new Vector(ge.length/ 2);
						for (int m = 0; m < ge.length; m++) {
							Vector geV= (Vector) chrHash.get(ge[m].getChromosome());
							if (geV== null) {
								sp.remove(ge[m], true);
								continue;
							}
							int n= 0;
							for (n = 0; n < geV.size(); n++) 
								if (ge[m].overlaps((Gene) geV.elementAt(n)))
									break;
							if (n== geV.size())
								sp.remove(ge[m], true);
						}
						ge= sp.getGenes();
						System.out.println("Found "+ge.length+" genes.");
						
						for (int m = 0; m < target.length; m++) {
							try {
								PrintStream p= new PrintStream(subdir+ File.separator+
										target[m].getName()+"_"+comSpeName+"_"+specAnnos[k]);
								target[m].invoke(null, new Object[] {sp, p, null});
								p.flush();
								p.close();
							} catch (Exception e) {
								e.printStackTrace();
							} 
						}
					}
					
				}
				
			}
			
		}
			
	}
	
	static public void _01_analyzeASEvents(Gene ge, PrintStream p, HashMap map) {
		final String ID_DEGREE= "degree";
		final String ID_CPLXEV= "cplxEv";
		final String ID_SPLEV= "splEv";
		final String ID_SHORTSC= "shortSc";
		final String ID_LONGSC= "longSc";
		final String ID_NB_TRPTS= "nbTrpts";
		final String ID_NB_EXONS= "nbExons";
	
		final String ID_DEGREE_CDS= "deg_cds";
		final String ID_CPLXEV_CDS= "cplxEv_cds";
		final String ID_SPLEV_CDS= "splEv_cds";
		final String ID_SHORTSC_CDS= "shortSc_cds";
		final String ID_LONGSC_CDS= "longSc_cds";
		final String ID_NB_TRPTS_CDS= "nbTrpts_cds";
		final String ID_NB_EXONS_CDS= "nbExons_cds";
		
		
		final String ID_DEGREE_5UTR= "deg_5utr";
		final String ID_CPLXEV_5UTR= "cplxEv_5utr";
		final String ID_SPLEV_5UTR= "splEv_5utr";
		final String ID_SHORTSC_5UTR= "shortSc_5utr";
		final String ID_LONGSC_5UTR= "longSc_5utr";
		
		IntVector degreeV= (IntVector) map.get(ID_DEGREE);
		if (degreeV== null) {
			degreeV= new IntVector();
			map.put(ID_DEGREE, degreeV);
		}
		Integer cplxEvents= (Integer) map.get(ID_CPLXEV);
		if (cplxEvents== null) 
			cplxEvents= new Integer(0);
		Integer splEvents= (Integer) map.get(ID_SPLEV);
		if (splEvents== null) 
			splEvents= new Integer(0);
		Integer shortSc= (Integer) map.get(ID_SHORTSC);
		if (shortSc== null) 
			shortSc= new Integer(0);
		Integer longSc= (Integer) map.get(ID_LONGSC);
		if (longSc== null) 
			longSc= new Integer(0);
		IntVector nbTrptV= (IntVector) map.get(ID_NB_TRPTS);
		if (nbTrptV== null) {
			nbTrptV= new IntVector();
			map.put(ID_NB_TRPTS, nbTrptV);
		}
		IntVector nbExonV= (IntVector) map.get(ID_NB_EXONS);
		if (nbExonV== null) {
			nbExonV= new IntVector();
			map.put(ID_NB_EXONS, nbExonV);
		}
		
			// CDS
		IntVector degreeCDSV= (IntVector) map.get(ID_DEGREE_CDS);
		if (degreeCDSV== null) {
			degreeCDSV= new IntVector();
			map.put(ID_DEGREE_CDS, degreeCDSV);
		}
		Integer cplxEventsCDS= (Integer) map.get(ID_CPLXEV_CDS);
		if (cplxEventsCDS== null) 
			cplxEventsCDS= new Integer(0);
		Integer splEventsCDS= (Integer) map.get(ID_SPLEV_CDS);
		if (splEventsCDS== null) 
			splEventsCDS= new Integer(0);
		Integer shortScCDS= (Integer) map.get(ID_SHORTSC_CDS);
		if (shortScCDS== null) 
			shortScCDS= new Integer(0);
		Integer longScCDS= (Integer) map.get(ID_LONGSC_CDS);
		if (longScCDS== null) 
			longScCDS= new Integer(0);
		IntVector nbTrptCDSV= (IntVector) map.get(ID_NB_TRPTS_CDS);
		if (nbTrptCDSV== null) {
			nbTrptCDSV= new IntVector();
			map.put(ID_NB_TRPTS_CDS, nbTrptCDSV);
		}
		IntVector nbExonCDSV= (IntVector) map.get(ID_NB_EXONS_CDS);
		if (nbExonCDSV== null) {
			nbExonCDSV= new IntVector();
			map.put(ID_NB_EXONS_CDS, nbExonCDSV);
		}
		
		// 5UTR
		IntVector degree5UTRV= (IntVector) map.get(ID_DEGREE_5UTR);
		if (degree5UTRV== null) {
			degree5UTRV= new IntVector();
			map.put(ID_DEGREE_5UTR, degree5UTRV);
		}
		Integer cplxEvents5UTR= (Integer) map.get(ID_CPLXEV_5UTR);
		if (cplxEvents5UTR== null) 
			cplxEvents5UTR= new Integer(0);
		Integer splEvents5UTR= (Integer) map.get(ID_SPLEV_5UTR);
		if (splEvents5UTR== null) 
			splEvents5UTR= new Integer(0);
		Integer shortSc5UTR= (Integer) map.get(ID_SHORTSC_5UTR);
		if (shortSc5UTR== null) 
			shortSc5UTR= new Integer(0);
		Integer longSc5UTR= (Integer) map.get(ID_LONGSC_5UTR);
		if (longSc5UTR== null) 
			longSc5UTR= new Integer(0);
	
			// collect
		if (ge!= null) {
			ASVariation[] vars= ge.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			Vector geneV= new Vector();
			for (int m = 0; vars!= null&& m < vars.length; m++) {
				int deg= vars[m].getDegree();
				degreeV.add(deg);
				if (deg> 2) {
					cplxEvents= new Integer(cplxEvents.intValue()+ 1);
					map.put(ID_CPLXEV, cplxEvents);
					int[] sym= vars[m].getBalance();
					shortSc= new Integer(shortSc.intValue()+ sym[0]);
					map.put(ID_SHORTSC, shortSc);
					longSc= new Integer(longSc.intValue()+ sym[1]);
					map.put(ID_LONGSC, longSc);
					if (vars[m].is_affecting_CDS()) {
						cplxEventsCDS= new Integer(cplxEventsCDS.intValue()+ 1);
						map.put(ID_CPLXEV_CDS, cplxEventsCDS);
						shortScCDS= new Integer(shortScCDS.intValue()+ sym[0]);
						map.put(ID_SHORTSC_CDS, shortScCDS);
						longScCDS= new Integer(longScCDS.intValue()+ sym[1]);
						map.put(ID_LONGSC_CDS, longScCDS);
					}
					if (vars[m].is_affecting_5UTR()) {
						cplxEvents5UTR= new Integer(cplxEvents5UTR.intValue()+ 1);
						map.put(ID_CPLXEV_5UTR, cplxEvents5UTR);
						shortSc5UTR= new Integer(shortSc5UTR.intValue()+ sym[0]);
						map.put(ID_SHORTSC_5UTR, shortSc5UTR);
						longSc5UTR= new Integer(longSc5UTR.intValue()+ sym[1]);
						map.put(ID_LONGSC_5UTR, longSc5UTR);
					}
						
				} else {
					splEvents= new Integer(splEvents.intValue()+ 1);
					map.put(ID_SPLEV, splEvents);
					if (vars[m].is_affecting_CDS()) { 
						splEventsCDS= new Integer(splEventsCDS.intValue()+ 1);
						map.put(ID_SPLEV_CDS, splEventsCDS);
					}
					if (vars[m].is_affecting_5UTR()) { 
						splEvents5UTR= new Integer(splEvents5UTR.intValue()+ 1);
						map.put(ID_SPLEV_5UTR, splEvents5UTR);
					}
				}
			}
	
			gphase.tools.Arrays.sort2DFieldRev(vars);
			
			nbTrptV.add(ge.getTranscriptCount());
			int cnt= 0;
			for (int i = 0; i < ge.getTranscripts().length; i++) {
				if (ge.getTranscripts()[i].isCoding())
					++cnt;
			}
			nbTrptCDSV.add(cnt);
			
			nbExonV.add(ge.getExons().length);
			cnt= 0;
			for (int i = 0; i < ge.getExons().length; i++) {
				if (ge.getExons()[i].isCodingCompletely())
					++cnt;
			}
			nbExonCDSV.add(cnt);
			
		}
		
			// degree
		if (p!= null) {
			System.out.println("#sVar\t#cVar\tsVar/cVar\tshortSC\tlongSC\tshort/long\tDEGREE med\tmea\tstd\tGENE\tTRPT/G\tEXON/G");
			p.println("#sVar\t#cVar\tsVar/cVar\tshortSC\tlongSC\tshort/long\tDEGREE med\tmea\tstd\tGENE\tTRPT/G\tEXON/G");
			System.out.print(splEvents+"\t"+cplxEvents+"\t"+Formatter.fprint(((float) splEvents/cplxEvents), 2)+"\t");
			p.print(splEvents+"\t"+cplxEvents+"\t"+Formatter.fprint(((float) splEvents/cplxEvents), 2)+"\t");
			System.out.print(shortSc+"\t"+longSc+"\t"+Formatter.fprint(((float) shortSc/longSc), 2)+"\t");
			p.print(shortSc+"\t"+longSc+"\t"+Formatter.fprint(((float) shortSc/longSc), 2)+"\t");
			Distribution dist= new Distribution(degreeV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			dist= new Distribution(nbTrptV.toIntArray());
			System.out.print(dist.getMedian()+"\t");
			p.print(dist.getMedian()+"\t");
			dist= new Distribution(nbExonV.toIntArray());
			System.out.println(dist.getMedian());
			p.println(dist.getMedian());
			
			System.out.print(splEventsCDS+"\t"+cplxEventsCDS+"\t"+Formatter.fprint(((float) splEventsCDS/cplxEventsCDS), 2)+"\t");
			p.print(splEventsCDS+"\t"+cplxEventsCDS+"\t"+Formatter.fprint(((float) splEventsCDS/cplxEventsCDS), 2)+"\t");
			System.out.print(shortScCDS+"\t"+longScCDS+"\t"+Formatter.fprint(((float) shortScCDS/longScCDS), 2)+"\t");
			p.print(shortScCDS+"\t"+longScCDS+"\t"+Formatter.fprint(((float) shortScCDS/longScCDS), 2)+"\t");
			dist= new Distribution(degreeCDSV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			dist= new Distribution(nbTrptCDSV.toIntArray());
			System.out.print(dist.getMedian()+"\t");
			p.print(dist.getMedian()+"\t");
			dist= new Distribution(nbExonCDSV.toIntArray());
			System.out.println(dist.getMedian());
			p.println(dist.getMedian());
	
			System.out.print(splEvents5UTR+"\t"+cplxEvents5UTR+"\t"+Formatter.fprint(((float) splEvents5UTR/cplxEvents5UTR), 2)+"\t");
			p.print(splEvents5UTR+"\t"+cplxEvents5UTR+"\t"+Formatter.fprint(((float) splEvents5UTR/cplxEvents5UTR), 2)+"\t");
			System.out.print(shortSc5UTR+"\t"+longSc5UTR+"\t"+Formatter.fprint(((float) shortSc5UTR/longSc5UTR), 2)+"\t");
			p.print(shortSc5UTR+"\t"+longSc5UTR+"\t"+Formatter.fprint(((float) shortSc5UTR/longSc5UTR), 2)+"\t");
			dist= new Distribution(degree5UTRV.toIntArray());
			System.out.println(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			p.println(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			
			//Structurator.writeHTML(vars, dir+File.separator+list[j]);
		}
		
	}
	static public void _01_analyzeASEvents2(Gene ge, PrintStream p, HashMap map) {
		final String ID_NB_ALL_GENES= "allGenes";
		final String ID_NB_STE= "nbSingleExonOrTranscript";
		final String ID_NB_AS_GENES= "nbGenesAS";
		final String ID_NB_TRPTS= "nbTrpts";	// wo STE
		final String ID_NB_EXONS= "nbExons";
		
		final String ID_EV_GENE= "eventsPerGene";
		final String ID_SPL_EV_GENE= "nbSplEvGene";
		final String ID_CPLX_EV_GENE= "nbCplEvGene";
		final String ID_CPLX_DEGREE= "cplxDegree";
		final String ID_SHORT_LEN= "shortLen";
		final String ID_LONG_LEN= "longLen";

		
		Integer nbAllGenes= (Integer) map.get(ID_NB_ALL_GENES);
		if (nbAllGenes== null) 
			nbAllGenes= new Integer(0);
		IntVector cplxDegreeV= (IntVector) map.get(ID_CPLX_DEGREE);
		if (cplxDegreeV== null) {
			cplxDegreeV= new IntVector();
			map.put(ID_CPLX_DEGREE, cplxDegreeV);
		}
		IntVector shortLenV= (IntVector) map.get(ID_SHORT_LEN);
		if (shortLenV== null) {
			shortLenV= new IntVector();
			map.put(ID_SHORT_LEN, shortLenV);
		}
		IntVector longLenV= (IntVector) map.get(ID_LONG_LEN);
		if (longLenV== null) {
			longLenV= new IntVector();
			map.put(ID_LONG_LEN, longLenV);
		}
		Integer nbSTE= (Integer) map.get(ID_NB_STE);
		if (nbSTE== null) 
			nbSTE= new Integer(0);
		
		IntVector nbTrptV= (IntVector) map.get(ID_NB_TRPTS);
		if (nbTrptV== null) {
			nbTrptV= new IntVector();
			map.put(ID_NB_TRPTS, nbTrptV);
		}
		IntVector nbExonV= (IntVector) map.get(ID_NB_EXONS);
		if (nbExonV== null) {
			nbExonV= new IntVector();
			map.put(ID_NB_EXONS, nbExonV);
		}
		Integer nbASgenes= (Integer) map.get(ID_NB_AS_GENES);
		if (nbASgenes== null) 
			nbASgenes= new Integer(0);
		
		IntVector eventsGeneV= (IntVector) map.get(ID_EV_GENE);
		if (eventsGeneV== null) {
			eventsGeneV= new IntVector();
			map.put(ID_EV_GENE, eventsGeneV);
		}
		IntVector cplxEvGeneV= (IntVector) map.get(ID_CPLX_EV_GENE);
		if (cplxEvGeneV== null) {
			cplxEvGeneV= new IntVector();
			map.put(ID_CPLX_EV_GENE, cplxEvGeneV);
		}
		IntVector splEvGeneV= (IntVector) map.get(ID_SPL_EV_GENE);
		if (splEvGeneV== null) {
			splEvGeneV= new IntVector();
			map.put(ID_SPL_EV_GENE, splEvGeneV);
		}
		

			// collect
		if (ge!= null) {
				// gene attrribs
			nbAllGenes= new Integer(nbAllGenes.intValue()+ 1);
			map.put(ID_NB_ALL_GENES, nbAllGenes);
			
			int cnt= 0;
			for (int i = 0; i < ge.getTranscripts().length; i++) {
				if (ge.getTranscripts()[i].isCoding())
					++cnt;
			}
			int cnt2= 0;
			for (int i = 0; i < ge.getExons().length; i++) {
				if (ge.getExons()[i].isCoding())
					++cnt2;
			}
			
			if (ge.getTranscripts().length< 2|| ge.getExons().length< 2) {
				nbSTE= new Integer(nbSTE.intValue()+ 1);
				map.put(ID_NB_STE, nbSTE);
			} else {
				nbTrptV.add(cnt);
				nbExonV.add(cnt2);
			}
			
				// AS vars
			ASVariation[] vars= ge.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			int cntCplx= 0, cntSpl= 0;
			for (int m = 0; vars!= null&& m < vars.length; m++) {				
				if (!vars[m].is_affecting_CDS()) 
					continue;
				int deg= vars[m].getDegree();
				if (deg> 2) {
					cplxDegreeV.add(deg);
					++cntCplx;
				} else {
					++cntSpl;
				}
			}
			int sum= cntSpl+ cntCplx;
			if (sum> 0) {
				cplxEvGeneV.add(cntCplx);
				splEvGeneV.add(cntSpl);
				eventsGeneV.add(vars.length);
				nbASgenes= new Integer(nbASgenes.intValue()+ 1);
				map.put(ID_NB_AS_GENES, nbASgenes);
			}
	
		}
		
			// output
		if (p!= null) {
			System.out.println("#allG\t#STEG\t[%]\t#ASG\t[%]\tcTrp/G\tcEx/G\t" +
					"ev/G\ttotEv\tsplEv/G\ttotSpl\tcplxEv/G\ttotCplx\tcplxDeg\tmean\tstdev");
			p.println("#allG\t#STEG\t[%]\t#ASG\t[%]\tcTrp/G\tcEx/G\t" +
					"ev/G\ttotEv\tsplEv/G\ttotSpl\tcplxEv/G\ttotCplx\tcplxDeg\tmean\tstdev");
			
			System.out.print(nbAllGenes+"\t"+nbSTE+"\t"+Formatter.fprint(((float) nbSTE* 100/nbAllGenes), 2)+"\t");
			p.print(nbAllGenes+"\t"+nbSTE+"\t"+Formatter.fprint(((float) nbSTE* 100/nbAllGenes), 2)+"\t");
			System.out.print(nbASgenes+"\t"+Formatter.fprint(((float) nbASgenes* 100/(nbAllGenes- nbSTE)), 2)+"\t");
			p.print(nbASgenes+"\t"+Formatter.fprint(((float) nbASgenes* 100/(nbAllGenes- nbSTE)), 2)+"\t");
			Distribution dist= new Distribution(nbTrptV.toIntArray());
			System.out.print(dist.getMedian()+"\t");
			p.print(dist.getMedian()+"\t");
			dist= new Distribution(nbExonV.toIntArray());
			System.out.print(dist.getMedian()+"\t");
			p.print(dist.getMedian()+"\t");
			dist= new Distribution(eventsGeneV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			dist= new Distribution(splEvGeneV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			dist= new Distribution(cplxEvGeneV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getTotal()+"\t");
			dist= new Distribution(cplxDegreeV.toIntArray());
			System.out.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			p.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
			
			System.out.println();
			p.println();
		}
		
	}
	static public void _01_analyzeASEvents(Species sp, PrintStream p, HashMap map) {
		
		ASVariation[][] vars= sp.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);		
		IntVector degreeV= new IntVector();
		Vector geneV= new Vector();
		int splEvents= 0, cplxEvents= 0, shortSc= 0, longSc= 0;
		for (int m = 0; vars!= null&& m < vars.length; m++) {		
			int deg= vars[m][0].getDegree();
			for (int n = 0; n < vars[m].length; n++) { 
				degreeV.add(deg);
				geneV= (Vector) gphase.tools.Arrays.addUnique(geneV, vars[m][n].getTranscript1().getGene());
			}
			if (deg> 2)
				cplxEvents+= vars[m].length;
			else
				splEvents+= vars[m].length;
			int[] sym= vars[m][0].getBalance();
			shortSc+= sym[0];
			longSc+= sym[1];
		}

			// transcripts
		IntVector nbTrptV= new IntVector();
		IntVector nbExonV= new IntVector();
		for (int m = 0; m < geneV.size(); m++) {
			Gene gen= (Gene) geneV.elementAt(m);
			nbTrptV.add(gen.getTranscriptCount());
			nbExonV.add(gen.getExons().length);
		}
		
		System.out.println("#sVar\t#cVar\tsVar/cVar\tshortSC\tlongSC\tshort/long\tDEGREE med\tmea\tstd\tGENE\tTRPT/G\tEXON/G");
		p.println("#sVar\t#cVar\tsVar/cVar\tshortSC\tlongSC\tshort/long\tDEGREE med\tmea\tstd\tGENE\tTRPT/G\tEXON/G");
		System.out.print(splEvents+"\t"+cplxEvents+"\t"+Formatter.fprint(((float) splEvents/cplxEvents), 2)+"\t");
		p.print(splEvents+"\t"+cplxEvents+"\t"+Formatter.fprint(((float) splEvents/cplxEvents), 2)+"\t");
		System.out.print(shortSc+"\t"+longSc+"\t"+Formatter.fprint(((float) shortSc/longSc), 2)+"\t");
		p.print(shortSc+"\t"+longSc+"\t"+Formatter.fprint(((float) shortSc/longSc), 2)+"\t");
		Distribution dist= new Distribution(degreeV.toIntArray());
		System.out.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
		p.print(dist.getMedian()+"\t"+dist.getMean()+"\t"+dist.getStandardDeviation()+"\t");
		dist= new Distribution(nbTrptV.toIntArray());
		System.out.print(dist.getMedian()+"\t");
		p.print(dist.getMedian()+"\t");
		dist= new Distribution(nbExonV.toIntArray());
		System.out.println(dist.getMedian());
		p.println(dist.getMedian());
		//Structurator.writeHTML(vars, dir+File.separator+list[j]);
	}
		
		
	static void _00_checkHubNoHubOverlap() {
		final String[] DOUGLAS_SPECIES= new String[] {
			"human", "mouse", "fruitfly", "worm"
		};
		
		String subdir= "douglas";
		String[] list= new File(subdir).list();
		for (int i = 0; i < DOUGLAS_SPECIES.length; i++) {
			if (!DOUGLAS_SPECIES[i].equals("fruitfly"))
				continue;
			String speName= getDouglasAbbreviation(DOUGLAS_SPECIES[i]);
			String hubName= null;
			String noHubName= null;
			for (int j = 0; j < list.length; j++) { 
				if (!list[j].startsWith(speName)|| new File(subdir+ File.separator+ list[j]).isDirectory())
					continue;
				if (list[j].contains(".nonhubs"))
					noHubName= list[j];
				else if (list[j].contains(".hubs"))
					hubName= list[j];
			}

			TabDelimitedFormatWrapper wrapper= new TabDelimitedFormatWrapper(subdir+ File.separator+ hubName);
			try {
				wrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			String[] hubIDs= wrapper.getColumn(0);
			
			wrapper= new TabDelimitedFormatWrapper(subdir+ File.separator+ noHubName);
			try {
				wrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			String[] noHubIDs= wrapper.getColumn(0);
			
			Arrays.sort(noHubIDs);
			Vector v= new Vector();
			for (int j = 0; j < hubIDs.length; j++) {
				int p= Arrays.binarySearch(noHubIDs, hubIDs[j]);
				if (p>= 0) 
					v.add(hubIDs[j]);
			}
			
			System.out.println(speName+" overlaps by "+v.size()+" on proteinID level, "
					+((float) v.size()/hubIDs.length)+((float) v.size()/noHubIDs.length));
			
			String[] hubGenIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], hubIDs);
			String[] noHubGenIDs= MapTable.getPrimaryGeneID(DOUGLAS_SPECIES[i], noHubIDs);
			Arrays.sort(noHubGenIDs);
			v= new Vector();
			for (int j = 0; j < hubGenIDs.length; j++) {
				int p= Arrays.binarySearch(noHubGenIDs, hubGenIDs[j]);
				if (p>= 0) 
					v.add(hubGenIDs[j]);
			}

			System.out.println(speName+" overlaps by "+v.size()+" on gene level, "
					+((float) v.size()/hubGenIDs.length)+((float) v.size()/noHubGenIDs.length));
			
		}
	}
}
