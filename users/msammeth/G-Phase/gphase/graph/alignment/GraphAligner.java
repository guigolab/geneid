package gphase.graph.alignment;

import java.awt.Point;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Vector;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// import prefuse.data.Graph;

import gphase.Constants;
import gphase.algo.ASAnalyzer;
import gphase.db.EnsemblDBAdaptor;
import gphase.db.OrthologousGeneWrapper;
import gphase.graph.SpliceGraph;
import gphase.graph.SpliceNode;
import gphase.graph.SplicePath;
import gphase.graph.gui.GraphView;
import gphase.graph.gui.PFGraph;
import gphase.io.BackwardFileReader;
import gphase.io.gtf.Andre;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.tools.DoubleVector;
import gphase.tools.Formatter;

public class GraphAligner {

	static boolean inclFlanks= false; 
	static String outDir= "andre";
	
	public static void alignENSEMBLHomologGenesAndre(String startSpe1, String startSpe2, String startID1, String startID2) {
		
		if (startSpe1== null&& startID1!= null)
			startSpe1= Species.getCommonNameForPrefix(Species.decodeEnsemblPfx(startID1));
		if (startSpe2== null&& startID2!= null)
			startSpe2= Species.getCommonNameForPrefix(Species.decodeEnsemblPfx(startID2));
		if (startSpe1!= null|| startSpe2!= null)
			System.out.println("skipping to spe1 "+startSpe1+" x spe2 "+startSpe2+": "+
					startID1+ " x "+startID2);
		
		
		String[] speIDs= Species.SP_NAMES_ANDRE;
		for (int i = 0; i < speIDs.length; i++) {
			if (startSpe1!= null&& !startSpe1.equals(speIDs[i]))
				continue;
			if (startSpe1!= null) {
				if (!startSpe1.equals(speIDs[i]))
					continue;
				startSpe1= null;
			}
			Species spe1= new Species(speIDs[i]);
			Andre andre1= new Andre(Constants.DATA_DIR+ File.separator+ Constants.SEQUENCES_SUBDIR 
					+ File.separator+ "caipirinha"+ File.separator+ spe1.getAbbreviatedName()+".gff");
			for (int j = (i+1); j < speIDs.length; j++) {
				if (startSpe2!= null) {
					if (!startSpe2.equals(speIDs[j]))
						continue;
					startSpe2= null;
				}
				Species spe2= new Species(speIDs[j]);
				Andre andre2= new Andre(Constants.DATA_DIR+ File.separator+ Constants.SEQUENCES_SUBDIR 
						+ File.separator+ "caipirinha"+ File.separator+ spe2.getAbbreviatedName()+".gff");
				OrthologousGeneWrapper ortho= new OrthologousGeneWrapper(spe1, spe2);
				
				File f= new File(outDir+File.separator+
							spe1.getAbbreviatedName()+"-"+spe2.getAbbreviatedName()+".ali");
				if (f.exists()) {
					try {
						BackwardFileReader breader= new BackwardFileReader(f, "rw");
						String line= breader.readLineBackward();
						while (!line.startsWith("ENS"))
							line= breader.readLineBackward();
						String[] tokens= line.split("\t");
						startID1= tokens[0];
						startID2= tokens[1];
						startSpe1= Species.getCommonNameForPrefix(Species.decodeEnsemblPfx(tokens[0]));
						startSpe2= Species.getCommonNameForPrefix(Species.decodeEnsemblPfx(tokens[1]));
						System.out.println("skipping to spe1 "+startSpe1+" x spe2 "+startSpe2+": "+
								startID1+ " x "+startID2);
						byte[] b= new byte[(int) (f.length()- breader.getFilePointer())];	// +1 -> but goes to the CR to add
						for (int k = 0; k < b.length; k++) 
							b[k]= ' ';
						breader.writeChar('\n');
						breader.write(b);
						breader.writeChar('\n');
						breader.close();
						if (!startSpe2.equals(speIDs[j]))
							continue;
					} catch (Exception e) {
						e.printStackTrace();
					}
					
				}
				
					// output log
				PrintStream p= null; 
				try {
					if (f.exists())
						p= new PrintStream(new FileOutputStream(f, true));	// append
					else
						p= new PrintStream(f.getAbsolutePath());
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
				
				GeneHomology homol= ortho.getNextGeneHomology();
				while (homol!= null) {
					System.gc();
					if ((startID1!= null&& !homol.getGene1().getGeneID().equals(startID1))||
							(startID2!= null&& !homol.getGene2().getGeneID().equals(startID2))) {
						homol= ortho.getNextGeneHomology();
						continue;
					}
					if (startID1!= null|| startID2!= null) {
						startID1= null;
						startID2= null;
					}
					String iter= homol.getGene1().getGeneID()+"\t"+homol.getGene2().getGeneID()+ "\t";
					System.out.print(iter);
					System.out.flush();
					p.print(iter);
					p.flush();
//					if (!iter.startsWith("ENSG00000003393\tENSMUSG00000026024")) {
//						homol= ortho.getNextGeneHomology();
//						continue;	// REMOVE
//					}
//					System.out.println();
					
					long t0= System.currentTimeMillis();
					andre1.read(homol.getGene1().getGeneID());
					Gene ge1= andre1.assembleGene(spe1);
					andre2.read(homol.getGene2().getGeneID());
					Gene ge2= andre2.assembleGene(spe1);
					if (ge1.getSpliceSites()== null|| ge1.getSpliceSites().length< 1||
							ge2.getSpliceSites()== null|| ge2.getSpliceSites().length< 1) {
						System.out.println();
						p.println();
						homol= ortho.getNextGeneHomology();
						continue;
					}

					long t1= System.currentTimeMillis();
					SpliceGraph g1= new SpliceGraph(ge1.getTranscripts());
					g1.init(inclFlanks);
					SpliceGraph g2= new SpliceGraph(ge2.getTranscripts());
					g2.init(inclFlanks);
					long t2= System.currentTimeMillis();
					
					String s= "\t";
					try {
						Mapping[] maps= align2_DP(g1, g2, -1, p);
						s+= "\tcost: "+maps[0].getCost();
					} catch (OutOfMemoryError e) {
						System.err.println(e.getMessage());
						p.println();
					}
					
					long t3= System.currentTimeMillis();
					
						// measures
					s+= "\t"+				
						TimeUnit.SECONDS.convert(t1-t0, TimeUnit.MILLISECONDS)+ " "+
						TimeUnit.SECONDS.convert(t2-t1, TimeUnit.MILLISECONDS)+ " "+
						TimeUnit.SECONDS.convert(t3-t2, TimeUnit.MILLISECONDS)+ " = "+
						TimeUnit.SECONDS.convert(t3-t0, TimeUnit.MILLISECONDS)+"sec";
					System.out.println(s);
					//p.println(s);
					
					
					homol= ortho.getNextGeneHomology();
					//break;	// REMOVE
				}
				p.flush(); p.close();
			}
		}
	}
	
	static void alignHomologs(Gene[] hGenes) {
		int cntFail= 0;
		for (int i = 0; i < hGenes.length; i++) {
			for (int j = i+1; j < hGenes.length; j++) {
				SpliceGraph g1= new SpliceGraph(hGenes[i].getTranscripts());
				g1.init(false);
				SpliceGraph g2= new SpliceGraph(hGenes[j].getTranscripts());
				g2.init(false);
				try {
					Mapping[] maps= align(g1, g2, -1);
					if (maps!= null)
						System.out.println(maps[0].toString());
				} catch (OutOfMemoryError e) {
					++cntFail; // :)
					e.printStackTrace();
				}
			}
		}
	}
	
	public static void main(String[] args) {
		
		String cost= "Len";
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-flanks"))
				inclFlanks= true;
			else if (args[i].equalsIgnoreCase("-outDir")&& i+1< args.length)
				outDir= args[i+1];
			else if (args[i].equalsIgnoreCase("-cost")&& i+1< args.length)
				cost= args[i+1];
		}
		
		System.out.print("Transcript flanks (tss/tes) ");
		if (!inclFlanks)
			System.out.print("NOT ");
		System.out.println("included.");
		System.out.println("Cost scheme "+cost);
		
		File f= new File(outDir);
		if (!f.exists())
			f.mkdir();
		System.out.println("Output directory "+f.getAbsolutePath());
		
		
		Mapping.setCostID(cost);
		//alignENSEMBLHomologGenes();
		alignENSEMBLHomologGenesAndre(null, null, null, null);
	}
	
	public static void test() {
		Species spec= new Species("human");
		Gene ge1= new Gene(spec, "ge1");
		ge1.setChromosome("1");
		ge1.setStrand(1);
		Transcript t11= new Transcript(ge1, "t11");
		t11.setStrand(1);
		t11.addExon(new Exon(t11, "e11_1", 1, 2));
		t11.addExon(new Exon(t11, "e11_2", 3, 4));
		t11.addExon(new Exon(t11, "e11_3", 500, 600));
		t11.addExon(new Exon(t11, "e11_4", 800, 900));
		ge1.addTranscript(t11);
		Transcript t12= new Transcript(ge1, "t12");
		t12.setStrand(1);
		t12.addExon(new Exon(t12, "e12_1", 1, 2));
		t12.addExon(new Exon(t12, "e12_2", 3, 4));
		t12.addExon(new Exon(t12, "e12_3", 500, 700));
		t12.addExon(new Exon(t12, "e12_4", 800, 900));
		ge1.addTranscript(t12);
		SpliceGraph g1= new SpliceGraph(new Transcript[] {t11,t12});
		g1.init();
		g1.getBubbles();
		
		Gene ge2= new Gene(spec, "ge2");
		ge2.setChromosome("2");
		ge2.setStrand(1);
		Transcript t21= new Transcript(ge2, "t21");
		t21.setStrand(1);
		t21.addExon(new Exon(t21, "e21_1", 1, 2));
		t21.addExon(new Exon(t21, "e21_2", 300, 400));
		t21.addExon(new Exon(t21, "e21_3", 600, 700));
		ge2.addTranscript(t21);
		Transcript t22= new Transcript(ge2, "t22");
		t22.setStrand(1);
		t22.addExon(new Exon(t22, "e22_1", 1, 2));
		t22.addExon(new Exon(t22, "e22_2", 300, 500));
		t22.addExon(new Exon(t22, "e22_3", 600, 700));
		ge2.addTranscript(t22);
		SpliceGraph g2= new SpliceGraph(new Transcript[] {t21,t22});
		g2.init();
		g2.getBubbles();
		
		Mapping[] maps= align(g1, g2);
		for (int i = 0; i < maps.length; i++) {
			System.out.println(maps[i]);
		}
		//GraphView.demo(ge1.getGeneID(), new PFGraph(g1), ge2.getGeneID(), new PFGraph(g2));
	}
	
	SpliceGraph g1, g2;
	
	public GraphAligner(SpliceGraph newG1, SpliceGraph newG2) {
		
	}
	
	public static Mapping[] align(SpliceGraph g1, SpliceGraph g2, int maxAli) {
		
		Comparator compi= new SpliceNode.PositionTypeComparator();
		SpliceNode[] listI= g1.getNodeList();
		Arrays.sort(listI, compi);
		SpliceNode[] listJ= g2.getNodeList();
		Arrays.sort(listJ, compi);
		MappingPriorityQueue q= new MappingPriorityQueue(Math.max(listI.length, listJ.length), 
				new Mapping.PriorityComparator());
		
		if (listI.length> 0) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(listI[0], null);
			q.add(m);
		}
		if (listJ.length> 0) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(null, listJ[0]);
			q.add(m); 
		}
		if (listI.length> 0&& listJ.length> 0) {	// both unaligned
			Mapping m= new Mapping(g1, g2);			
			m.addMapping(listI[0], null);	
			m.addMapping(null, listJ[0]);
			m.maxRel= 4;
			q.add(m); 
		}
		if ((listI.length> 0&& listJ.length> 0)&& isAligneable(listI[0], listJ[0])) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(listI[0], listJ[0]);
			q.add(m);
		}
		Mapping map= (Mapping) q.poll();
		double ulCost= getNaiveCost(g1, g2, listI, listJ); 	// Double.MAX_VALUE;
		Vector optMaps= new Vector();
		//int maxAli= 5;
		while (map!= null&& map.getCost()<= ulCost) {
			
			if (maxAli> 0&& optMaps.size()== maxAli)
				break;
	
				// generate possibilities for next i, next j
			int nextI= -1, nextJ= -1;
			if (map.getMaxI()!= null)
				nextI= Arrays.binarySearch(listI, map.getMaxI(), compi);
			if (map.getMaxJ()!= null)
				nextJ= Arrays.binarySearch(listJ, map.getMaxJ(), compi);
	
			String s= Formatter.fprint(map.getCost(), 2)+" ";
			if (map.maxRel== 1)
				s+= nextI+",-\t->\t";
			else if (map.maxRel== 2)
				s+= "-,"+nextJ+"\t->\t";
			else if (map.maxRel== 3)
				s+= nextI+","+nextJ+"\t->\t";
			else if (map.maxRel== 4)
				s+= nextI+"|"+nextJ+"\t->\t";
			
	
			if (((nextI+ 1< listI.length)&& (nextJ+ 1< listJ.length))&& 
					isAligneable(listI[nextI+ 1], listJ[nextJ+ 1])) {	// both aligned
				Mapping m= null;
				try {
					m= (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(listI[nextI+ 1], listJ[nextJ+ 1]);
				System.out.println(s+(nextI+1)+","+(nextJ+1)+" "+Formatter.fprint(m.getCost(), 2));
				q.add(m);
			}
			if (nextI+ 1< listI.length&& nextJ+ 1< listJ.length&& 
					map.getMaxRel()!= 1&& map.getMaxRel()!= 2) {	// both unaligned
				Mapping m= null;
				try {
					m= (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(listI[nextI+ 1], null);
				m.addMapping(null, listJ[nextJ+ 1]);
				m.maxRel= 4;
				System.out.println(s+(nextI+1)+"|"+(nextJ+1)+" "+Formatter.fprint(m.getCost(), 2));
				q.add(m);
			}
			if (nextI+ 1< listI.length&& map.getMaxRel()!= 2) {	// one gapped (allow only chains of gaps)
				Mapping m= null;
				try {
					m = (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(listI[nextI+ 1], null);
				System.out.println(s+(nextI+1)+",-"+" "+Formatter.fprint(m.getCost(), 2));
				q.add(m);
			}
			if (nextJ+ 1< listJ.length&& map.getMaxRel()!= 1) {	// other gapped (allow only chains of gaps)
				Mapping m= null;
				try {
					m= (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(null, listJ[nextJ+ 1]);
				System.out.println(s+"-,"+(nextJ+1)+" "+Formatter.fprint(m.getCost(), 2));
				q.add(m);
			}
			
				// update cheapest path
			if ((nextI+ 1>= listI.length)&& (nextJ+ 1>= listJ.length)) {
				if (map.getCost()<= ulCost) {
					if (maxAli< 0&& map.getCost()< ulCost) {
						ulCost= map.getCost();
					}
					optMaps.add(map);
				}
			}
			
			map= (Mapping) q.poll();	// next
		}
		
		return (Mapping[]) gphase.tools.Arrays.toField(optMaps);
	}

	public static Mapping[] align2(SpliceGraph g1, SpliceGraph g2, int maxAli) {
			
			Mapping.setCostID("LenExp");
			Comparator compi= new SpliceNode.PositionTypeComparator();
			SpliceNode[] listI= g1.getNodeList();
			Arrays.sort(listI, compi);
			SpliceNode[] listJ= g2.getNodeList();
			Arrays.sort(listJ, compi);
			MappingPriorityQueue q= new MappingPriorityQueue(Math.max(listI.length, listJ.length), 
					new Mapping.PriorityExtensionComparator());
			
				// exhaustively align roots
			for (int i = 0; i < g1.getRoots().length; i++) {
				for (int j = 0; j < g2.getRoots().length; j++) {
					Mapping m= new Mapping(g1, g2);
					m.addMapping(g1.getRoots()[i], g2.getRoots()[j]);
					q.add(m);
				}
			}
			
			Mapping map= (Mapping) q.poll();
			double ulCost= getNaiveCost2(g1, g2, listI, listJ); 	// Double.MAX_VALUE;
			Vector optMaps= new Vector();
			//int maxAli= 5;
			while (map!= null&& map.getCost()<= ulCost) {
				
				if (maxAli> 0&& optMaps.size()== maxAli)
					break;
	
					// generate possibilities for next i, next j
				int nextI= -1, nextJ= -1;
				if (map.getMaxI()!= null)
					nextI= Arrays.binarySearch(listI, map.getMaxI(), compi);
				if (map.getMaxJ()!= null)
					nextJ= Arrays.binarySearch(listJ, map.getMaxJ(), compi);
		
				
				SpliceNode[] nextPair= map.getNextAlignPair(listI, listJ, nextI, nextJ);
				if (nextPair== null) {
					map= (Mapping) q.poll();
					continue;
				} else 
					q.add(map);
				
				int toI= Arrays.binarySearch(listI, nextPair[0], compi);
				int toJ= Arrays.binarySearch(listJ, nextPair[1], compi);
				System.out.println(nextI+","+nextJ+" "+
						Formatter.fprint(map.getCost(), 2)+ "\t->\t"+
						toI+","+toJ+" "+isAligneable(nextPair[0], nextPair[1])+" "+
						Formatter.fprint(map.getCost()+ map.getCost(nextPair[0], nextPair[1]), 2)+ " "+
						(map.getCost()+ map.getCost(nextPair[0], nextPair[1])< ulCost)+
						" "+map.offIJ);
				
	
				if (isAligneable(nextPair[0], nextPair[1])&&
						map.getCost()+ map.getCost(nextPair[0], nextPair[1])< ulCost) {
					Mapping m= null;
					try {
						m= (Mapping) map.clone();
					} catch (CloneNotSupportedException e) {
						e.printStackTrace();
					}
					m.addMapping(nextPair[0], nextPair[1]);
					q.add(m);
				}
				
					// align all resting SSs
	//			for (int i = nextI+ 1; i < listI.length; i++) {
	//				for (int j = 0; j < listJ.length; j++) {
	//					if (isAligneable(listI[i], listJ[j])&&
	//							map.getCost()+ map.getCost(listI[i], listJ[j])< ulCost) {
	//						Mapping m= null;
	//						try {
	//							m= (Mapping) map.clone();
	//						} catch (CloneNotSupportedException e) {
	//							e.printStackTrace();
	//						}
	//						m.addMapping(listI[i], listJ[j]);
	//						q.add(m);
	//					}
	//				}
	//			}
				
					// update cheapest path
				if ((nextI+ 1>= listI.length)&& (nextJ+ 1>= listJ.length)) {
					if (map.getCost()<= ulCost) {
						if (maxAli< 0&& map.getCost()< ulCost) {
							ulCost= map.getCost();
						}
						optMaps.add(map);
					}
				}
				
				map= (Mapping) q.poll();	// next
			}
			
			return (Mapping[]) gphase.tools.Arrays.toField(optMaps);
		}

	public static Mapping[] align2_DP(SpliceGraph g1, SpliceGraph g2, int maxAli, PrintStream pr) {
		
		Comparator compi= new SpliceNode.PositionTypeComparator();
		SpliceNode[] listI= g1.getNodeList();
		Arrays.sort(listI, compi);
		SpliceNode[] listJ= g2.getNodeList();
		Arrays.sort(listJ, compi);
		MappingPriorityQueue q= new MappingPriorityQueue(Math.max(listI.length, listJ.length), 
				new Mapping.PriorityExtensionComparator());
		
			// init matrix
		double[][] matrix= new double[listI.length][];
		for (int i = 0; i < matrix.length; i++) {
			matrix[i]= new double[listJ.length];
			for (int j = 0; j < matrix[i].length; j++) 
				matrix[i][j]= -1d;
		}
		Point[][] pointer= new Point[listI.length][listJ.length];
		for (int i = 0; i < matrix.length; i++) {
			pointer[i]= new Point[listJ.length];
			for (int j = 0; j < matrix[i].length; j++) 
				pointer[i][j]= null;
		}

		// init exhaustively align roots
		for (int i = 0; i < g1.getRoots().length; i++) {
			int rI= Arrays.binarySearch(listI, g1.getRoots()[i], compi);
			for (int j = 0; j < g2.getRoots().length; j++) {
				int rJ= Arrays.binarySearch(listJ, g2.getRoots()[j], compi);
				matrix[rI][rJ]= 0d;
			}
		}
		

			// calc matrix
		double ul= getNaiveCost2_DP(g1, g2, listI, listJ);
		Vector sinkV= new Vector();
		for (int i = 0; i < listI.length; i++) {
			for (int j = 0; j < listJ.length; j++) {
				if (isAligneable(listI[i], listJ[j])) {
						// iterate back matrix
					for (int k = 0; k < i; k++) {
						for (int l = 0; l < j; l++) {
							if (matrix[k][l]< 0)
								continue;
							double addCost= getCost_exonicLenExpGrowing(g1, g2, listI[k], listJ[l], listI[i], listJ[j]);
							if (addCost< 0)
								continue;	
							double newCost= matrix[k][l]+ addCost;
							if (newCost> ul) 
								continue;	// lipman
							if (matrix[i][j]< 0) {
								matrix[i][j]= newCost;
								pointer[i][j]= new Point(k, l);
							} else if (newCost<= matrix[i][j]) {
								matrix[i][j]= newCost;
								pointer[i][j]= new Point(k, l);
							}							
						}
					}
					if (matrix[i][j]>= 0&& listI[i].getOutDegree()< 1&& listJ[j].getOutDegree()< 1)
						sinkV.add(new Point(i,j));
				} else 
					matrix[i][j]= -1d;
			}
		}
		
			// retrieve best solution
		Point p= null;
		double min= Double.MAX_VALUE;
		for (int i = 0; i < sinkV.size(); i++) {
			Point tmP= (Point) sinkV.elementAt(i);
			if (matrix[tmP.x][tmP.y]< min) {
				min= matrix[tmP.x][tmP.y];
				p= tmP;
			}
		}
		
			// trace back 
		Vector v= new Vector();
		while (p!= null)  {
			v.add(p);
			p= pointer[p.x][p.y];
		}
		Mapping m= new Mapping(g1, g2);
		for (int i = v.size()- 1; i >= 0; --i) {
			p= (Point) v.elementAt(i);
			m.addMapping(listI[p.x], listJ[p.y]);
		}
		m.cost= min;
		
		String s= m.toStringAndre();
		System.out.println(s);
		pr.println(s);

		System.out.println(m.toString(listI, listJ));
		
		return new Mapping[] {m};
	}

	/**
		 * Performs all naive gapless alignment between the two splicechains.
		 * 
		 * @param g1
		 * @param g2
		 * @return
		 */
		public static strictfp double getNaiveCost(SpliceGraph g1, SpliceGraph g2, SpliceNode[] listI, SpliceNode[] listJ) {
			
			int diff= Math.abs(listI.length- listJ.length);
			boolean rev= false;
			if (listI.length> listJ.length) {	// shorter in listI
				SpliceNode[] h= listI;
				listI= listJ;
				listJ= h;
				rev= true;
			}
			double ulCost= Double.MAX_VALUE;
			Mapping ulMap= null;
			for (int i = 0; i < diff; i++) {
				Mapping m= null;
				if (rev) 
					m= new Mapping(g2, g1);		// for finding correctly the pathes afterwards
				else
					m= new Mapping(g1, g2);
				
					// heading
				for (int j = 0; j < i; j++) 
					m.addMapping(null, listJ[j]);
	
					// align
				int iAli= 0;
				for (int j = 0; iAli < listI.length&& (j+i)< listJ.length; j++) {
					if (isAligneable(listI[iAli], listJ[j+i])) 
						m.addMapping(listI[iAli++], listJ[j+i]);
					else
						m.addMapping(null, listJ[j+i]);	// leave longer unaligned
				}
				
					// trainling
				if (iAli< listI.length) {	// longer ran empty
					for (int j = iAli; j < listI.length; j++) 
						m.addMapping(listI[j], null);
				
				} else {					// shorter list ran empty
					int rest= (diff-i);
					for (int j = 0; j < rest; j++) 
						m.addMapping(null, listJ[listI.length+ j]);
				}
				
				if (m.getCost()< ulCost) {
					ulCost= m.getCost();
					ulMap= m;
				}
			}
			
	//		if (rev)
	//			System.out.println("\nul\n"+ulMap.toStringInverse());
	//		else
	//			System.out.println("\nul\n"+ulMap.toString());
	
			return ulCost;
		}

	/**
	 * Performs all naive gapless alignment between the two splicechains.
	 * 
	 * @param g1
	 * @param g2
	 * @return
	 */
	public static strictfp double getNaiveCost2(SpliceGraph g1, SpliceGraph g2, SpliceNode[] listI, SpliceNode[] listJ) {
		
		int diff= Math.abs(listI.length- listJ.length);
		boolean rev= false;
		if (listI.length> listJ.length) {	// shorter in listI
			SpliceNode[] h= listI;
			listI= listJ;
			listJ= h;
			rev= true;
		}
		double ulCost= Double.MAX_VALUE;
		Mapping ulMap= null;
		for (int i = 0; i < diff; i++) {
			Mapping m= null;
			if (rev) 
				m= new Mapping(g2, g1);		// for finding correctly the pathes afterwards
			else
				m= new Mapping(g1, g2);
			
				// align
			int iAli= 0;
			for (int j = 0; iAli < listI.length&& (j+i)< listJ.length; j++) {
				if (isAligneable(listI[iAli], listJ[j+i])) 
					m.addMapping(listI[iAli++], listJ[j+i]);
			}
			if ((i+iAli-1)>= listJ.length)
				continue;	// bigger list ran empty
			
			if (m.getCost()< ulCost) {
				ulCost= m.getCost();
				ulMap= m;
			}
		}
		
		System.out.println("ul: "+ ulCost);
		if (rev)
			System.out.println(ulMap.toStringInverse(listI, listJ));
		else
			System.out.println(ulMap.toString(listI, listJ));
	
		return ulCost;
	}

	/**
		 * Performs all naive gapless alignment between the two splicechains.
		 * 
		 * @param g1
		 * @param g2
		 * @return
		 */
		public static strictfp double getNaiveCost2_DP(SpliceGraph g1, SpliceGraph g2, SpliceNode[] listI, SpliceNode[] listJ) {
			
			int diff= Math.abs(listI.length- listJ.length);
			boolean rev= false;
			if (listI.length> listJ.length) {	// shorter in listI
				SpliceNode[] h= listI;
				listI= listJ;
				listJ= h;
				rev= true;
			}
			double ulCost= Double.MAX_VALUE;
			Vector ulMap= null;
			DoubleVector ulScores= null;
			for (int i = 0; i < diff; i++) {
				
					// align
				int iAli= 0;
				double cost= 0d;
				SpliceNode[] last= new SpliceNode[] {null, null};			
				Vector v= new Vector();
				DoubleVector scoreV= new DoubleVector();
				for (int j = 0; iAli < listI.length&& (j+i)< listJ.length; j++) {
					if (isAligneable(listI[iAli], listJ[j+i])) {
						double addCost;
						if (rev)
							addCost= getCost_exonicLenExpGrowing(g2, g1, last[0], last[1], listI[iAli], listJ[j+i]);
						else
							addCost= getCost_exonicLenExpGrowing(g1, g2, last[0], last[1], listI[iAli], listJ[j+i]);
						if (addCost< 0)
							continue;
						
						cost+= addCost;
						scoreV.add(addCost);
						
						Point p;
						if (rev)
							p= new Point(j+i, iAli);
						else
							p= new Point(iAli, j+i);
						v.add(p);
						last= new SpliceNode[] {listI[iAli], listJ[j+i]};
						++iAli;
					}
				}
				
				
				if (iAli< listI.length|| (i+iAli-1)>= listJ.length)
					continue;	// bigger list ran empty
				
				if (cost< ulCost) {
					ulCost= cost;
					ulMap= v;
					ulScores= scoreV;
				}
			}
			
//			System.out.println("ul: "+ ulCost);
//			for (int i = 0; i < ulMap.size(); i++) {
//				Point p= (Point) ulMap.elementAt(i);
//				System.out.print("["+p.x+","+p.y+"], ");
//			}
//			System.out.println();
//			for (int i = 0; i < ulMap.size(); i++) 
//				System.out.print(Formatter.fprint(ulScores.get(i), 4)+", ");
//			System.out.println();
	//		if (rev)
	//			System.out.println(ulMap.toStringInverse(listI, listJ));
	//		else
	//			System.out.println(ulMap.toString(listI, listJ));
	
			return ulCost;
		}

	/**
	 * Performs all naive gapless alignment between the two splicechains.
	 * 
	 * @param g1
	 * @param g2
	 * @return
	 */
	public static strictfp double getNaiveCost2_DP_Score(SpliceGraph g1, SpliceGraph g2, SpliceNode[] listI, SpliceNode[] listJ) {
		
		int diff= Math.abs(listI.length- listJ.length);
		boolean rev= false;
		if (listI.length> listJ.length) {	// shorter in listI
			SpliceNode[] h= listI;
			listI= listJ;
			listJ= h;
			rev= true;
		}
		double ulCost= Double.MAX_VALUE;
		Vector ulMap= null;
		DoubleVector ulScores= null;
		for (int i = 0; i < diff; i++) {
			
				// align
			int iAli= 0;
			double cost= 0d;
			SpliceNode[] last= new SpliceNode[] {null, null};			
			Vector v= new Vector();
			DoubleVector scoreV= new DoubleVector();
			for (int j = 0; iAli < listI.length&& (j+i)< listJ.length; j++) {
				if (isAligneable(listI[iAli], listJ[j+i])) {
					double addCost;
					if (rev)
						addCost= getCost_exonicLenExp(g2, g1, last[0], last[1], listI[iAli], listJ[j+i]);
					else
						addCost= getCost_exonicLenExp(g1, g2, last[0], last[1], listI[iAli], listJ[j+i]);
					if (addCost< 0)
						continue;
					
					cost+= addCost;
					scoreV.add(addCost);
					
					Point p;
					if (rev)
						p= new Point(j+i, iAli);
					else
						p= new Point(iAli, j+i);
					v.add(p);
					last= new SpliceNode[] {listI[iAli], listJ[j+i]};
					++iAli;
				}
			}
			
			
			if (iAli< listI.length|| (i+iAli-1)>= listJ.length)
				continue;	// bigger list ran empty
			
			if (cost> llCost) {
				llCost= cost;
				ulMap= v;
				ulScores= scoreV;
			}
		}
		
		System.out.println("ul: "+ ulCost);
		for (int i = 0; i < ulMap.size(); i++) {
			Point p= (Point) ulMap.elementAt(i);
			System.out.print("["+p.x+","+p.y+"], ");
		}
		System.out.println();
		for (int i = 0; i < ulMap.size(); i++) 
			System.out.print(Formatter.fprint(ulScores.get(i), 4)+", ");
		System.out.println();
//		if (rev)
//			System.out.println(ulMap.toStringInverse(listI, listJ));
//		else
//			System.out.println(ulMap.toString(listI, listJ));

		return ulCost;
	}
	
	static boolean isAligneable(SpliceNode s1, SpliceNode s2) {
		
		if (s1== null|| s2== null)
			return true;		
		
		if ((s1.isDonor()== true&& s2.isDonor()== true)||
				(s1.isAcceptor()== true&&  s2.isAcceptor()== true)||
				(s1.isTSS()== true&&  s2.isTSS()== true)||
				(s1.isTES()== true&&  s2.isTES()== true))
			return true;
		return false;
	}

	public static void alignENSEMBLHomologGenes() {
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			Graph g= adaptor.getGraphAllHomologs(EnsemblDBAdaptor.SPECIES_ISMB);
	//		g.filterNonsense();
	//		GraphHandler.writeOut(g, "graph_filt-nons.oos");
	//		g.filterNonCodingTranscripts();
	//		GraphHandler.writeOut(g, "graph_filt-nons_filt-nc.oos");
	//
			//EnsemblDBAdaptor.removeNotAllHomologGenes(g);
			// "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish", "fruitfly", "mosquito"
			// 5e "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish"
			String[] speNames= new String[] {"human", "mouse", "dog"};
			EnsemblDBAdaptor.removeNotAllHomologGenes(g, speNames);
			//GraphHandler.writeOut(g, "graph_filt-nons_filt-nc_remAll.oos");
	
			int ctr= 0;
			for (int i = 0; i < speNames.length; i++) {
				Species sp= g.getSpeciesByName(Species.getBinomialForCommonName(speNames[i]));
				Gene[] ge= sp.getGenes();
				for (int j = 0; j < ge.length; j++) {
					//System.out.println(ge[j].getGeneID());
					++ctr;
				}
			}
			System.out.println("==\nHomolog gene nb: "+ctr+"/"+speNames.length+"="+(ctr/speNames.length));
			
			Species sp= g.getSpeciesByName(Species.getBinomialForCommonName(speNames[0]));
			Gene[] ge= sp.getGenes();
			Vector v= new Vector(ctr/speNames.length);
			for (int j = 0; j < ge.length; j++) 
				v.add(ge[j]);
			
			alignHomologs((Gene[]) gphase.tools.Arrays.toField(v));
			
		}

	static strictfp double getCost_exonicLenExpScore(SpliceGraph g1, SpliceGraph g2, SpliceNode fromI, SpliceNode fromJ, SpliceNode nI, SpliceNode nJ) {
		
		if (fromI== null&& fromJ== null)	// for aligning roots
			return 0d; 
		
			// get new pathes (since last aligned node or root)
		SplicePath[] pathesI= g1.findPathes(fromI, nI);
		SplicePath[] pathesJ= g2.findPathes(fromJ, nJ);
		
			// at least one with NO valid pathes (eg, src of both graphs)
		if ((pathesI== null|| pathesI.length== 0)|| (pathesJ== null|| pathesJ.length== 0))
			return -1d;	// not alignable
	
			// two with valid pathes, match table
		SplicePath[] lessPathes, morePathes;
		if (pathesI.length< pathesJ.length) {
			lessPathes= pathesI;
			morePathes= pathesJ;
		} else {
			lessPathes= pathesJ;
			morePathes= pathesI;
		}

			// get cheapest combination for relations, 
			// do not penalize non-present chains,
			// allow non-injectional mappings
		double sum= 0d;
		for (int i = 0; i < lessPathes.length; i++) {
			double max= 0d;
			for (int j = 0; j < morePathes.length; j++) {
				int vw= Math.abs(lessPathes[i].getExonicLength());
				int tu= Math.abs(morePathes[j].getExonicLength());
//				int vwExonEdges= lessPathes[i].getExonicPieces();
//				int tuExonEdges= morePathes[i].getExonicPieces();
//				int maxEEdges= Math.max(vwExonEdges, tuExonEdges);
				double lenRel;
				if (vw> 0&& tu> 0)
					lenRel= Math.min(((double) tu/ vw), ((double) vw/ tu));
				else if (vw== 0&& tu== 0)
					lenRel= 0d;	// no exons
				else {
					assert (vw== 0^ tu== 0);
					lenRel= 0d;	// here 1 is ok (div by 0)
				}
				if (lenRel> max)
					max= lenRel;
			}
			sum+= max;
		}
		
		return sum;
	}

	static strictfp double getCost_exonicLenExp(SpliceGraph g1, SpliceGraph g2, SpliceNode fromI, SpliceNode fromJ, SpliceNode nI, SpliceNode nJ) {
			
			if (fromI== null&& fromJ== null)	// for aligning roots
				return 0d; 
			
				// get new pathes (since last aligned node or root)
			SplicePath[] pathesI= g1.findPathes(fromI, nI);
			SplicePath[] pathesJ= g2.findPathes(fromJ, nJ);
			
				// at least one with NO valid pathes (eg, src of both graphs)
			if ((pathesI== null|| pathesI.length== 0)|| (pathesJ== null|| pathesJ.length== 0))
				return -1d;	// not alignable
		
				// two with valid pathes, match table
			SplicePath[] lessPathes, morePathes;
			if (pathesI.length< pathesJ.length) {
				lessPathes= pathesI;
				morePathes= pathesJ;
			} else {
				lessPathes= pathesJ;
				morePathes= pathesI;
			}
	
				// get cheapest combination for relations, 
				// do not penalize non-present chains,
				// allow non-injectional mappings
			double sum= 0d;
			for (int i = 0; i < lessPathes.length; i++) {
				double min= 1d;
				for (int j = 0; j < morePathes.length; j++) {
					int vw= Math.abs(lessPathes[i].getExonicLength());
					int tu= Math.abs(morePathes[j].getExonicLength());
					int vwExonEdges= lessPathes[i].getExonicPieces();
					int tuExonEdges= morePathes[i].getExonicPieces();
					int nrEdges= Math.min(vwExonEdges, tuExonEdges);
					double lenRel;
					if (vw> 0&& tu> 0) {
						lenRel= 1d-Math.min(((double) tu/ vw), ((double) vw/ tu));
						
					} else if (vw== 0&& tu== 0)
						lenRel= 0d;	// no exons
					else {
						assert (vw== 0^ tu== 0);
						lenRel= 1d;	// here 1 is ok (div by 0)
					}
					if (lenRel< min)
						min= lenRel;
				}
				sum+= min;
			}
			
			return sum;
		}

	static strictfp double getCost_exonicLenExpGrowing(SpliceGraph g1, SpliceGraph g2, SpliceNode fromI, SpliceNode fromJ, SpliceNode nI, SpliceNode nJ) {
		
		if (fromI== null&& fromJ== null)	// for aligning roots
			return 0d; 
		
			// get new pathes (since last aligned node or root)
		SplicePath[] pathesI= g1.findPathes(fromI, nI);
		SplicePath[] pathesJ= g2.findPathes(fromJ, nJ);
		
			// at least one with NO valid pathes (eg, src of both graphs)
		if ((pathesI== null|| pathesI.length== 0)|| (pathesJ== null|| pathesJ.length== 0))
			return -1d;	// not alignable
	
			// two with valid pathes, match table
		SplicePath[] lessPathes, morePathes;
		if (pathesI.length< pathesJ.length) {
			lessPathes= pathesI;
			morePathes= pathesJ;
		} else {
			lessPathes= pathesJ;
			morePathes= pathesI;
		}
	
			// get cheapest combination for relations, 
			// do not penalize non-present chains,
			// allow non-injectional mappings
		double sum= 0d;
		for (int i = 0; i < lessPathes.length; i++) {
			double min= 1d;
			for (int j = 0; j < morePathes.length; j++) {
				int[] vw= lessPathes[i].getExonicLengthes();
				int[] tu= morePathes[j].getExonicLengthes();
				double lenRel;
				if (vw!= null&& vw.length> 0&& tu!= null&& tu.length> 0) {
					int totVW= 0;
					for (int k = 0; k < vw.length; k++) 
						totVW+= vw[k];
					int totTU= 0;
					for (int k = 0; k < tu.length; k++) 
						totTU+= tu[k];
					lenRel= Math.min(((double) totTU/ totVW), ((double) totVW/ totTU));
					for (int k = 0; k < vw.length; k++) {
						for (int l = 0; l < tu.length; l++) {
							if (Math.min(((double) vw[k]/tu[l]), ((double) tu[l]/vw[k]))
									> lenRel)	// similarity rises when decomposing 
								return -1d;	// not valid
						}
					}
					lenRel= 1d- lenRel;
				} else if ((vw== null|| vw.length== 0)&& (tu== null|| tu.length== 0))
					lenRel= 0d;	// no exons
				else {
					//assert (vw== 0^ tu== 0);
					lenRel= 1d;	// here 1 is ok (div by 0)
				}
				if (lenRel< min)
					min= lenRel;
			}
			sum+= min;
		}
		
		return sum;
	}
}
