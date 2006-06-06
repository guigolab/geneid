/*
 * Created on Mar 1, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.algo;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Calendar;
import java.util.Collection;
import java.util.Comparator;
import java.util.Date;
import java.util.Vector;

import javax.swing.JFrame;
import javax.swing.JWindow;

import com.p6spy.engine.logging.appender.StdoutLogger;
import com.sun.org.apache.xalan.internal.xsltc.compiler.util.FilterGenerator;

import gphase.db.EnsemblDBAdaptor;
import gphase.gui.Paparazzi;
import gphase.gui.SpliceOSigner;
import gphase.io.gtf.EncodeWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.regex.Automaton;
import gphase.regex.ChessMaster;
import gphase.regex.RegExp;
import gphase.tools.Arrays;
import gphase.tools.Time;

/**
 * 
 * 
 * @author msammeth
 */
public class ASAnalyzer {

	Graph graph;
	public ASAnalyzer(Graph newGraph) {
		this.graph= newGraph;
	}

	public static void analyze1_transcripts_clusters(Graph g) {
		System.out.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
		ASAnalyzer.filterSingleExonClusters(g);
		int x= g.getGenes().length;
		int y= g.getTranscriptCount();
		System.out.println("clusters: "+x+" transcripts: "+y+"\t"+ ((float) y/ (float) x));
		int z= g.getClustersWithAS();
		System.out.println("clusters w as: "+z+"\t"+((float) z/ (float) x));
	}
	
	public static void filterSingleExonClusters(Graph g) {
		Species[] spe= g.getSpecies();
		for (int i = 0; i < spe.length; i++) {
			Gene[] ge= spe[i].getGenes();
			for (int j = 0; j < ge.length; j++) {
				if (ge[j].getExons().length== 1)
					g.removeKill(ge[j]);
			}
		}
	}
	
	public static void getAllASVariation(Graph g) {
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		outputVariations(classes, false, false, System.out);
		int cnt= 0;
		try {
			PrintStream printer= new PrintStream("out005_encode_vars_all_non-redundant");
			for (int i = 0; i < classes.length; i++) {
				cnt+= classes[i].length;
				for (int j = 0; j < classes[i].length; j++) {
					classes[i][j].outputDetail(printer);
				}
				printer.print("\n\n");
			}
			printer.flush();
			printer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(cnt);
	}
	public static void getVincenco(Graph g) {
	
		ASVariation[][] classes= g.getASVariations(1);
		ASVariation[] events= getVariation("( // 1^2=)", classes);
		int cnt1= 0, cnt2= 0, cnt3= 0, cnt4= 0;
		Comparator compi= new SpliceSite.PositionComparator();
		for (int i = 0; i < events.length; i++) {
			if ((events[i].getTranscript1().getSpliceChain().length> 0&&
					compi.compare(events[i].getBorderingSpliceSites()[0],
					events[i].getTranscript1().getSpliceChain()[0])== 0)||
				(events[i].getTranscript2().getSpliceChain().length> 0&& 
					compi.compare(events[i].getBorderingSpliceSites()[0],
						events[i].getTranscript2().getSpliceChain()[0])== 0)) {
				cnt1++;
				if (events[i].getLengthDiff(false)< 150) {
					cnt3++;
					if (events[i].isSemiProteinCoding()) {
						cnt4++;
						events[i].outputDetail(System.out);
					}
				}
			}
			if (events[i].getLengthDiff(false)< 150) {
				cnt2++;
			}
		} 
		
		System.out.println(cnt1+","+cnt2+","+cnt3+","+cnt4);
	}
	
	
	public static void lengthVariation(Graph g) {
			ASVariation[][] classes= g.getASVariations(1);
			classes= filter(classes, "isNotAtAllCoding");
			classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
			//writeLabels(classes);
			classes= sortLabels(classes);
			
			
			
				// for plots
//			ASVariation[][] vRev= new ASVariation[classes[0].length][];
//			for (int i = 0; i < vRev.length; i++) {
//				vRev[i]= new ASVariation[classes.length];
//				for (int j = 0; j < vRev[i].length; j++) {
//					if (classes[j]== null|| classes[j].length<= i) 
//						vRev[i][j]= null;
//					else
//						vRev[i][j]= classes[j][i];
//				}
//			}
//			classes= vRev;
//			
//			//ASVariation[] events= getVariation("( // 1^2=)", classes);
//			for (int x = 0; x < classes[0].length; x++) 
//				System.out.print(classes[0][x]+ ",");
//			System.out.println();
//			for (int x = 0; x < classes.length; x++) {
//				ASVariation[] events= classes[x];
//				for (int i = 0; i < events.length; i++) {
//					if (events[i]== null) {
//						System.out.print(",");
//						continue;
//					}
//					int[] a= events[i].getLength(true);
//					int diffA= Math.abs(a[0]- a[1]);
//	//				events[i].outputDetail(System.out);
//					System.out.print(diffA+ ",");
//				}
//				System.out.println();
//			}
		}

	public static void lengthVariationModulo(Graph g) {
		ASVariation[][] classes= g.getASVariations(1);
		classes= filter(classes, "isProteinCoding");
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		
		//ASVariation[] events= getVariation("( // 1^2=)", classes);
		int cnt0= 0, cnt1= 0, cnt2= 0;
		int scnt0= 0, scnt1= 0, scnt2= 0, stot= 0;
		int tot= 0;
		for (int x = 0; x < classes.length; x++) {
			cnt1= 0; cnt2= 0; cnt0= 0;
			tot= 0;
			ASVariation[] events= classes[x];
			System.out.println(events[0].toString());
			for (int i = 0; i < events.length; i++) {
				int[] a= events[i].getLength(true);
				int[] b= events[i].getLength(false);
				int diffA= Math.abs(a[0]- a[1]);
//				events[i].outputDetail(System.out);
//				System.out.println(a[0]+","+a[1]+": "+diffA+"("+(diffA%3)+")");
				if (diffA%3== 0)
					cnt0++;
				if (diffA%3== 1)
					cnt1++;
				if (diffA%3== 2)
					cnt2++;
			}
			tot+= events.length;
			System.out.println("0: "+cnt0+"\t1: "+cnt1+"\t2: "+cnt2+"\t("+tot+"): "+((double) cnt0/ tot));
			scnt1+= cnt1;
			scnt0+= cnt0;
			scnt2+= cnt2;
			stot+= tot;
		}
		System.out.println("===> 0: "+scnt0+"\t1: "+scnt1+"\t2: "+scnt2+"\t("+stot+"): "+((double) scnt0/ stot));
	}
	
	public static void outputCheckFile() {
		// "encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/44regions_genes_CHR_coord.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		String fName= "encode/44regions_genes_CHR_coord.gtf";
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		boolean encode= false;
		if (fName.equals("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		Graph g= myWrapper.getGraph(encode);
		
		PrintStream p= null;
		Calendar c= Calendar.getInstance();
			// check c.MONTH for bug?!
		String outName= "consistency_"+(c.get(c.MONTH)+ 1)+"_"+c.get(c.DAY_OF_MONTH)+"_"+c.get(c.HOUR_OF_DAY)+"_"+c.get(c.MINUTE);
		String sffx= ".chk";
		try {
			p= new PrintStream(outName+"_all"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//determineVariations(g, p);
		determineVariationsDetails(g, "isTrue", p);
		p.flush();
		p.close();
		
		try {
			p= new PrintStream(outName+"_cod"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isProteinCoding", p);
		p.flush();
		p.close();

		try {
			p= new PrintStream(outName+"_part"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isPartiallyCoding", p);
		p.flush();
		p.close();

		try {
			p= new PrintStream(outName+"_non"+sffx);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		determineVariationsDetails(g, "isNotAtAllCoding", p);
		p.flush();
		p.close();
}
	
	public static void outputASStatistics(Graph g) {
		g.getASVariations(ASMultiVariation.FILTER_NONE);
		int[] res= g.getSpliceSites(AbstractRegion.REGION_5UTR);
		System.out.println("5UTR: "+res[0]+" "+res[1]+"\t"+((float) res[0]/(float) res[1]));
		res= g.getSpliceSites(AbstractRegion.REGION_CDS);
		System.out.println("CDS: "+res[0]+" "+res[1]+"\t"+((float) res[0]/(float) res[1]));
		res= g.getSpliceSites(AbstractRegion.REGION_3UTR);
		System.out.println("3UTR: "+res[0]+" "+res[1]+"\t"+((float) res[0]/(float) res[1]));
	}
	/** @deprecated ?
	 * 
	 * 
	 * @param g
	 * @return
	 */
	public static double getPercASSpliceSites(Graph g) {
		Gene[] ge= g.getGenes();
		int altPS= 0;
		int totPS= 0;
		Vector v= new Vector();
		for (int i = 0; i < ge.length; i++) {
			SpliceSite[] ps= ge[i].getSpliceSites();
			if (ps== null)
				continue;
			totPS+= ps.length;
			for (int j = 0; j < ps.length; j++) {
				v.add(ps[j]);
				if (!ps[j].isConstitutive()) {
					++altPS;
				}
			}
		}
		BufferedWriter buffy= null;
		try {
			buffy = new BufferedWriter(new FileWriter("ss.txt"));
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		for (int i = 0; i < v.size(); i++) {
			SpliceSite ss= (SpliceSite) v.elementAt(i);
			System.out.println(ss.getGene().getChromosome()+" "+ss);
			try {
				buffy.write(ss.getGene().getChromosome()+" "+ss+"\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		try {
			buffy.flush();
			buffy.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println(altPS+","+totPS);
		return (double) altPS/ totPS;
	}
	
	public static void diff(ASVariation[][] classes1, ASVariation[][] classes2) {
		int cnt= 0;
		Comparator compi= null; //new ASVariation.IdentityComparator();
		for (int i = 0; i < classes1.length; i++) {
			String id= classes1[i][0].toString();
				// check missing in 2nd file
			int j;
			for (j = 0; j < classes2.length; j++) 
				if (classes2[j][0].toString().equals(id))
					break;
			if (j== classes2.length) {
				System.err.println("Missing in 2nd file "+id);
				for (int k = 0; k < classes2.length; k++) 
					classes1[i][k].outputDetail(System.out);
				cnt+= classes1[i].length;
				continue;
			}
			
				// check in 1st file
			for (int k = 0; k < classes1[i].length; k++) {
				int m;
				for (m = 0; m < classes2[j].length; m++) 
					if (compi.compare(classes1[i][k], classes2[j][m])== 0)
						break;
				if (m== classes2.length) {
					System.err.println("Missing in 2nd file ");
					classes1[i][k].outputDetail(System.out);
					++cnt;
				}
			}
			// check in 2nd file
			for (int k = 0; k < classes2[j].length; k++) {
				int m;
				for (m = 0; m < classes1[i].length; m++) 
					if (compi.compare(classes1[i][k], classes2[j][m])== 0)
						break;
				if (m== classes2.length) {
					System.err.println("Missing in 1st file ");
					classes1[i][k].outputDetail(System.out);
					++cnt;
				}
			}
		}
		
		//check missing in 1st file 
		for (int i = 0; i < classes2.length; i++) {
			String id= classes2[i][0].toString();
			int j;
			for (j = 0; j < classes1.length; j++) 
				if (classes1[j][0].toString().equals(id))
					break;
			if (j== classes1.length) {
				System.err.println("Missing in 1st file "+id);
				for (int k = 0; k < classes1[j].length; k++) 
					classes1[i][k].outputDetail(System.out);
				cnt+= classes1[i].length;
				continue;
			}
		}
			

	}
	
	public static void analyzeGene(Graph g, String geneID, String methodName, int filter) {
		Gene ge= g.getGene(geneID, "human");
		System.out.println("--- Report for Gene "+geneID+" ("+methodName+"): ---");
		ASVariation[] as= ge.getASVariations(filter);
		Method m= null;
		try {
			m = as[0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		as= filter(new ASVariation[][] {as}, m)[0];
		
		
		for (int i = 0; i < as.length; i++) {
			as[i].outputDetail(System.out);
			System.out.println();
		}
	}

	static void addVirtualGene1(Graph g) {
		Gene newGene= new Gene(g.getSpeciesByEnsemblPrefix("ENS"), "ENSG00001000001");
		Transcript trans1= new Transcript(newGene, "ENST00000100001");
		newGene.addTranscript(trans1);
		Transcript trans2= new Transcript(newGene, "ENST00000100002");
		newGene.addTranscript(trans2);
		Exon e= new Exon(trans1, "ENSE00000100001", 10,40);
		trans1.addExon(e);
		trans2.addExon(e);
		e= new Exon(trans1, "ENSE00000100001", 50,150);
		trans1.addExon(e);
		e= new Exon(trans1, "ENSE00000100002", 100,250);
		trans2.addExon(e);
		e= new Exon(trans1, "ENSE00000100003", 200,250);
		trans1.addExon(e);
		
		g.addGene(newGene);
	}
	
	public static ASVariation[] getVariation(String s, ASVariation[][] classes) {
		for (int i = 0; i < classes.length; i++) {
			if (classes[i][0].toString().equals(s))
				return classes[i];
		}
		return null;
	}
	
	public static ASVariation[] getAlternateDonors(ASVariation[][] classes) {
		Vector v= new Vector();
		for (int i = 0; i < classes.length; i++) {
			if ((classes[i][0].toString().startsWith("(1^")&& 
					classes[i][0].toString().contains("// 2^"))
					|| (classes[i][0].toString().startsWith("(2^")&& 
							classes[i][0].toString().contains("// 1^")))
				v.add(classes[i]);
		}
		return (ASVariation[]) Arrays.toField(v);
	}
	

	public static ASVariation[][] getSpecificVariations(Graph g, PrintStream p) {
		
		ASVariation[][] classes= g.getASVariations(0);
		Comparator compi= new ASVariation.SpliceChainUniqueCodingComparator();
		Comparator c= new ASVariation.SpliceChainCodingHierarchyFilter();
		ASVariation[][] newClasses= new ASVariation[classes.length][];
		for (int i = 0; i < classes.length; i++) {
			Vector v= new Vector();
			for (int j = 0; j < classes[i].length; j++) {
				for (int k = (j+1); k < classes[i].length; k++) {
					if (compi.compare(classes[i][j], classes[i][k])== 0) {	 
						v.add(classes[i][j]);	// keep one
					} else {
						for (int m = 0; m < v.size(); m++) {				
							if ((c.compare(v.elementAt(m), classes[i][k])== 0)
									|| (c.compare(v.elementAt(m), classes[i][j])== 0))
								v.remove(m);	// remove both
						}
					}
				}
			}
			newClasses[i]= (ASVariation[]) Arrays.toField(v);
		}
		
		return outputVariations(newClasses, true, false, System.out);
	}
	/**
	 * @param g
	 * @param respectCoding
	 * @return
	 */
	public static ASVariation[][] determineVariations(Graph g, PrintStream p) {
		
		ASVariation[][] classes= g.getASVariations(1);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);
		
		try {
			Method m = classes[0][0].getClass().getMethod("isTrue", null);
			ASVariation[][] filtClasses= filter(classes, m);
			outputVariations(filtClasses, true, false, p);
			
			m = classes[0][0].getClass().getMethod("isProteinCoding", null);
			filtClasses= filter(classes, m);
			outputVariations(filtClasses, true, false, p);
	
			m = classes[0][0].getClass().getMethod("isPartiallyCoding", null);
			filtClasses= filter(classes, m);
			outputVariations(filtClasses, true, false, p);
	
			m = classes[0][0].getClass().getMethod("isNotAtAllCoding", null);
			filtClasses= filter(classes, m);
			return outputVariations(filtClasses, true, false, p);
		} catch (Exception e) {
			e.printStackTrace();
		} 
		return null;
	}

	public static void determineVariationsDetails(Graph g, String methodName, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		for (int i = 0; i < classes.length; i++) {
			for (int j = 0; j < classes[i].length; j++) {
				classes[i][j].outputDetail(p);
			}
		}
	}

	public static ASVariation[][] determineVariations(Graph g, String methodName, PrintStream p) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		return outputVariations(classes, true, false, p);
	}

	public static ASVariation[][] determineVariations(Graph g, String methodName) {
		ASVariation[][] classes= g.getASVariations(1);
		Method m= null;
		try {
			m = classes[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		classes= filter(classes, m);
		classes= (ASVariation[][]) Arrays.sort2DFieldRev(classes);	
		return outputVariations(classes, true, false, System.out);
	}

	static String[] readLabels() {
		String fName= "labelSort.txt";
		Vector v= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			while (buffy.ready()) {
				String s= buffy.readLine().trim();
				if (s.length()> 0)
					v.add(s);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return (String[]) Arrays.toField(v);
	}
	
	
	static void writeLabels(ASVariation[][] classes) {
		String fName= "labelSort.txt";
		Vector v= new Vector();
		try {
			BufferedWriter buffy= new BufferedWriter(new FileWriter(fName));
			for (int i = 0; i < classes.length; i++) {
				buffy.write(classes[i][0].toString()+"\n");
			}
			buffy.flush();
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	public static ASVariation[][] sortLabels(ASVariation[][] classes) {
		String[] labels= readLabels();
		Vector v= new Vector(classes.length);
		for (int i = 0; i < labels.length; i++) {
			int j;
			for (j = 0; j < classes.length; j++) { 
				if (classes[j][0].toString().equals(labels[i])) {
					v.add(classes[j]);
					break;
				}
			}
			if (j== classes.length) {
				v.add(null);
			}
			
		}
			
		return (ASVariation[][]) Arrays.toField(v);

	}
	
	/**
	 * @param classes
	 * @param sortConsistently
	 * @param writeFile
	 * @param p
	 * @return
	 */
	public static ASVariation[][] outputVariations(ASVariation[][] classes, boolean sortConsistently, boolean writeFile, PrintStream p) {
		
		if (classes== null|| p== null)
			return classes;
		
			// sort
		int nullObj= 0;
		if (sortConsistently) {
			String[] labels= readLabels();
			Vector v= new Vector(classes.length);
			for (int i = 0; i < labels.length; i++) {
				int j;
				for (j = 0; j < classes.length; j++) { 
					if (classes[j]== null)
						continue;
					if (classes[j][0].toString().equals(labels[i])) {
						v.add(classes[j]);
						break;
					}
				}
				if (j== classes.length) {
					v.add(null);
					nullObj++;
				}
				
			}
				
			classes= (ASVariation[][]) Arrays.toField(v);
		}
		
		int threshold= 5;
		int totNum= 0, thrNum= 0;
		for (int i = 0; i < classes.length; i++) { 
			if (classes[i]== null)
				continue;
			totNum+= classes[i].length;
			if (classes[i]!= null&& classes[i].length>= threshold)
				thrNum+= classes[i].length;
		}
		p.println("ASVariations: "+ totNum+ ", structural distinct classes: "+ (classes.length- nullObj));
		int cnt1= 0, cnt2= 0;
		// for excluding low length population
//		int tot= 0;
//		for (int i = 0; i < classes.length; i++) {
//			if (classes[i]== null) 
//				continue;
//			int rel= 0;
//			for (int j = 0; j < classes[i].length; j++) 
//				if (classes[i][j].getLengthDiff(true)> 10)
//					tot++;
//		}
		for (int i = 0; i < classes.length; i++) {
			if (classes[i]== null) {
				p.println();
				continue;
			}
			// exlcuding short length population
//			int rel= 0;
//			for (int j = 0; j < classes[i].length; j++) {
//				if (classes[i][j].getLengthDiff(true)> 10)
//					rel++;
//			}
//			p.println(rel+ " "+ ((float) rel/ tot)+" "+classes[i][0]);
			p.println(classes[i].length+ " "+ ((float) classes[i].length/ totNum)+ " "+classes[i][0]);

			if (classes[i].length< threshold) {
				++cnt1;
				cnt2+= classes[i].length;
			}
		}
		//p.println("low represented ("+threshold+") events: "+ cnt1+ " ("+cnt2+")\n");
		p.println("\n\n");
		
		//GraphHandler.writeOut(g);
		if (writeFile)
			writeLabels(classes);
		
		return classes;
	}

	public static ASVariation[][] debug(Graph g, int howMany, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;
		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		for (int i = classes.length-1; i > classes.length-(10+1); i--) {
			System.out.println(classes[i][0]);
			
			SpliceSite[] sc= classes[i][0].getSpliceChain1();
			System.out.print(classes[i][0].getGene().getChromosome()+"\t");
			for (int j = 0; j < sc.length; j++) {
				System.out.print(sc[j]+" ");
			}
			System.out.println();
			
			sc= classes[i][0].getSpliceChain2();
			System.out.print(classes[i][0].getGene().getChromosome()+"\t");
			for (int j = 0; j < sc.length; j++) {
				System.out.print(sc[j]+" ");
			}
			System.out.println("\n\n");
		}
		return classes;
	}

	public static ASVariation[][] lowRepresented(Graph g, int howMany, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;
		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		for (int i = 0; i < classes.length; i++) {
			
			if (classes[i].length> howMany)
				continue;
			
			System.out.println("#"+classes[i].length+"\t"+classes[i][0]);
			for (int k = 0; k < classes[i].length; k++) {
				
				SpliceSite[] sc= classes[i][k].getSpliceChain1();
				System.out.print(classes[i][k].getGene().getChromosome()+"\t");
				for (int j = 0; j < sc.length; j++) {
					System.out.print(sc[j]+" ");
				}
				System.out.println();
				
				sc= classes[i][k].getSpliceChain2();
				System.out.print(classes[i][k].getGene().getChromosome()+"\t");
				for (int j = 0; j < sc.length; j++) {
					System.out.print(sc[j]+" ");
				}
				System.out.println();
			}
			System.out.println("\n\n");
		}
		return classes;
	}
	
	static Graph identifyOrthologSpliceSites(Graph g) {
		
		for (int i = 0; i < g.getSpecies().length; i++) {
			Gene[] ge= g.getSpecies()[i].getGenes();
			for (int j = 0; j < ge.length; j++) {
				GeneHomology[] h= ge[j].getHomologies();	// work on homologies
				for (int k = 0; k < h.length; k++) {
					if (!h[k].orthologSpliceSitesIdentified())
						h[k].identifyOrthologSpliceSites();	// only identify if necessary
				}
			}
		}
		
		return g;
	}
	
	static ASVariation[] findClass(Graph g, String s, boolean respectCoding) {
		int x= 0;
		if (respectCoding)
			x= 1;

		ASVariation[][] classes= g.getASVariations(x);
		classes= (ASVariation[][]) Arrays.sort2DField(classes);
		int totNum= 0;
		for (int i = 0; i < classes.length; i++) {
			if (classes[i][0].toString().equals(s))
				return classes[i];
		}
		return null;
	}
	
	static void filterEvents(File rawFile) {
		try {
		
			BufferedReader buffy= new BufferedReader(new FileReader(rawFile));
			BufferedWriter buffy2= new BufferedWriter(new FileWriter(rawFile.getPath()+ "_filtered"));
			while (buffy.ready()) {
				String s= buffy.readLine();
				if (s.length()< 80)
					buffy2.write(s+ "\n");				
			}
			buffy2.flush();
			buffy2.close();
		} catch (Exception e) {
			System.err.println("Problems:");
			e.printStackTrace();
		}
	}
	
	static ASVariation[][] match(Graph g, String r1, String r2, boolean respectCoding) {
		
		RegExp reg1= new RegExp(r1);
		RegExp reg2= new RegExp(r2);
		Automaton a1= reg1.toAutomaton();
		Automaton a2= reg2.toAutomaton();
		
		ASVariation[][] vars= g.getASVariations(1);
		int hits= 0;
		Vector hitClasses= new Vector();
		for (int i = 0; i < vars.length; i++) {
			a1.init(); a2.init();
			ChessMaster c= new ChessMaster(
					new Automaton[] {a1, a2},
					new SpliceSite[][] {vars[i][0].getSpliceChain1(), vars[i][0].getSpliceChain2()}
			);
			for (int j = 0; j < vars[i].length; j++) {
				if (vars[i][j].getGene().getStableID().equals("ENSG00000186777"))
					System.currentTimeMillis();
			}
			if (c.match()) {
				hits+= vars[i].length;
				hitClasses.add(vars[i]);
			}
		}
		
		return (ASVariation[][]) Arrays.toField(hitClasses);
	}
	
	/**
	 * 
	 * @param vars
	 * @param regions relative positions (1-based) of the splice sites including
	 * the regions for the length plot. Negative values are counted from the end. 
	 */
	public static void lengthPlot(Graph g, boolean respectCoding) {
		
		ASVariation[][] vars= g.getASVariations(1);
		int base= 0;
		for (int i = 0; i < vars.length; i++) 
			base+= vars[i].length;
		long t0= System.currentTimeMillis();
		ASVariation[][] hits= null;
		int filterVal= -1;
		
			// donors
		hits= match(g, "^D1$", "^D2(3,4)*$", true);
		//hits= filter5UTR(hits);
		int cnt= countHits(hits);
		System.out.println("\nDonor type 1 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^D2$", "^D1(3,4)*$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDonor type 2 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^,A,4$", "^,A,(1,2)*3$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nAcceptor type 1 ("+cnt+")");
		outputLength(hits);
		
		hits= match(g, "^,A,3$", "^,A,(1,2)*4$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nAcceptor type 2 ("+cnt+")");
		outputLength(hits);

		hits= match(g, "^,D,1,6$", "^,D,2(3,4)*5$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDon/Acc type 1 ("+cnt+")");
		outputLength(hits);

		hits= match(g, "^,D,2,5$", "^,D,1(3,4)*6$", true);
		//hits= filter5UTR(hits);
		cnt= countHits(hits);
		System.out.println("\nDon/Acc type 1 ("+cnt+")");
		outputLength(hits);
	}
	private static void outputLength(ASVariation[][] hits) {
		for (int i = 0; i < hits.length; i++) {
			for (int j = 0; j < hits[i].length; j++) {
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= (sUniverse[k+1].getPos()- sUniverse[k].getPos())+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)% 3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
			}
		}
	}
	private static ASVariation[][] filter(ASVariation[][] hits, String methodName) {
		Method m= null;
		try {
			m = hits[0][0].getClass().getMethod(methodName, null);
		} catch (SecurityException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (NoSuchMethodException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		return filter(hits, m);
	}
	private static ASVariation[][] filter(ASVariation[][] hits, Method target) {
		
		Vector result= new Vector(hits.length);
		for (int i = 0; i < hits.length; i++) {
			Vector v= new Vector(hits[i].length);
			for (int j = 0; j < hits[i].length; j++) {
				Object o= null;
				try {
					o= target.invoke(hits[i][j], null);
				} catch (IllegalArgumentException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (IllegalAccessException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} catch (InvocationTargetException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				if (((Boolean) o).booleanValue())
					v.add(hits[i][j]);
			}
			if (v.size()> 0)
				result.add(v);
		}
		
		return (ASVariation[][]) Arrays.toField(result);
	}
	
	private static int countHits(ASVariation[][] vars) {
		int sum= 0;
		for (int i = 0; i < vars.length; i++) 
			sum+= vars[i].length;
		return sum;
	}
	
	/**
	 * 
	 * @param vars
	 * @param regions relative positions (1-based) of the splice sites including
	 * the regions for the length plot. Negative values are counted from the end. 
	 */
	static void lengthPlot_old(Graph g) {
		
		ASVariation[][] vars= g.getASVariations(1);
		int base= 0;
		for (int i = 0; i < vars.length; i++) 
			base+= vars[i].length;
		long t0= System.currentTimeMillis();
		ASVariation[][] hits= null;
		int filterVal= -1;
		
			// donors
		hits= match(g, "^D1$", "^D2(3,4)*$", true);
		int hitcount= 0;
		int filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("Donor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		int nul= 0, eins= 0, zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= (sUniverse[k+1].getPos()- sUniverse[k].getPos())+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)% 3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");
		
		hits= match(g, "^D2$", "^D1(3,4)*$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		int eq= 0, diff11= 0, diff12= 0, diff21= 0, diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[1].getPos()- sUniverse[0].getPos();
				int skip= 0;
				for (int k = 2; k < sUniverse.length; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				if (s== null)
					s= hits[i][j].getTranscript2().getSuccSpliceSite(sUniverse[sUniverse.length- 1]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+"\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
		
			// acceptors
		hits= match(g, "^,A,4$", "^,A,(1,2)*3$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nAcceptor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		nul= 0; eins= 0; zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos();
				int skip= 0;
				for (int k = 0; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos() + 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getPredSpliceSite(sUniverse[0]);
				if (s== null)
					s= hits[i][j].getTranscript2().getPredSpliceSite(sUniverse[0]);
				int b= (s.getPos()- sUniverse[sUniverse.length- 1].getPos())+ 1;
				int x= (alt+skip)%3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+ alt+"\t"+skip+"\t"+
						x+ "\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");

		hits= match(g, "^,A,3$", "^,A,(1,2)*4$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nAcceptor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		eq= 0; diff11= 0; diff12= 0; diff21= 0; diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos();
				int skip= 0;
				for (int k = 0; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				int a= sUniverse[2].getPos()- sUniverse[1].getPos()+ 1;
				SpliceSite s= hits[i][j].getTranscript1().getPredSpliceSite(sUniverse[0]);
				if (s== null)
					s= hits[i][j].getTranscript2().getPredSpliceSite(sUniverse[0]);
				int b= (sUniverse[0].getPos()- s.getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+ "\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
		
		hits= match(g, "^,D,1,6$", "^,D,2(3,4)*5$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor/Acceptor T1 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		nul= 0; eins= 0; zwei= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= (sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos())+
						(sUniverse[1].getPos()- sUniverse[0].getPos());
				
				int skip= 0;	
				for (int k = 2; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				
				int a= (sUniverse[2].getPos()- sUniverse[1].getPos())+ 1;
				int b= (sUniverse[sUniverse.length- 2].getPos()- sUniverse[sUniverse.length- 3].getPos())+ 1;
				int x= (alt+skip)%3;
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+"\t"+
						+x+"\t"+a+"\t"+b);
				if (x== 0)
					nul++;
				else if(x== 1)
					eins++;
				else if(x== 2)
					zwei++;
			}
		}
		System.out.println("--- 0:"+nul+", 1:"+eins+", 2:"+zwei+"("+hitcount+")");

		
		hits= match(g, "^,D,2,5$", "^,D,1(3,4)*6$", true);
		hitcount= 0;
		filterHits= 0;
		for (int i = 0; i < hits.length; i++) { 
			if (hits[i].length< filterVal)
				continue;
			hitcount+= hits[i].length;
			++filterHits;
		}
		System.out.println("\nDonor/Acceptor T2 ("+((float) hitcount/ base)+ ","+ hitcount+","+filterHits+")");
		eq= 0; diff11= 0; diff12= 0; diff21= 0; diff22= 0;
		hitcount= 0;
		for (int i = 0; i < hits.length; i++) {
			if (hits[i].length< filterVal)
				continue;
			int[][] pos= new int[hits[i][0].getSpliceUniverse().length/ 2][];
			for (int j = 0; j < hits[i].length; j++) {
				if (hits[i][j].isProteinCoding())
					continue;
				hitcount++;
				SpliceSite[] sUniverse= hits[i][j].getSpliceUniverse();
				int alt= (sUniverse[sUniverse.length-1].getPos()- sUniverse[sUniverse.length-2].getPos())+
						(sUniverse[1].getPos()- sUniverse[0].getPos());
				
				int skip= 0;	
				for (int k = 2; k < sUniverse.length- 2; k= k+2) 
					skip+= sUniverse[k+1].getPos()- sUniverse[k].getPos()+ 1;
				
				int a= (sUniverse[2].getPos()- sUniverse[1].getPos())+ 1;
				int b= (sUniverse[sUniverse.length- 2].getPos()- sUniverse[sUniverse.length- 3].getPos())+ 1;
				int x= (alt%3), y= (skip%3); 
				System.out.println(hits[i][j].getGene().getStableID()+"\t"+alt+"\t"+skip+"\t"+
						x+"\t"+y+"\t"+a+"\t"+b);
				if (x== y)
					eq++;
				else if(Math.abs(x-y)== 1) {
					if (x> y)
						diff11++;
					else 
						diff12++;
				} else if(Math.abs(x-y)== 2) {
					if (x> y)
						diff21++;
					else
						diff22++;
				}
			}
		}
		System.out.println("--- eq:"+eq+", 1>:"+diff11+", 1<"+diff12+", 2>"+diff21+", 2<"+diff22+"("+hitcount+")");
	}
	
	public static void main(String[] args) {
		
//		filterEvents(new File("human_graph_as_variations_non-redundant.log"));

		outputCheckFile();
		
//		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
//		//Graph g= null;
//		Graph g= adaptor.getGraphAllGenes(new Species("chimp"));
//		
//		lengthPlot(g);
//		//determineVariations(g);
//		if (1== 1)
//			System.exit(0);
//		
//		determineVariations(g);
//		
//		ASVariation[] var= findClass(g, "(1^3= // 2^4=)");
//
//		for (int i = 0;  var!= null&& i < var.length; i++) {
//			Gene ge= var[i].getGene();
//			SpliceOSigner painter= new SpliceOSigner(ge);
//			JFrame frame= new JFrame();
//			frame.getContentPane().add(painter);
//			frame.pack();
//			frame.setVisible(true);
//		}
//		
//		if (1== 1)
//			System.exit(0);
//		
//		ASVariation[][][] clazzes= g.getASEvents();
//		try {
//			Paparazzi.writeGIF(null, new File("test.gif"));
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//		if (1== 1)
//			System.exit(0);
//		clazzes= (ASVariation[][][]) Arrays.sort2DField(clazzes);
//		System.out.println("asSuperClasses "+ clazzes.length);		
//		for (int i = 0; i < clazzes.length; i++) {
//			System.out.println(clazzes[i].length+ " "+clazzes[i][0][0].toASEventString());
//
//			for (int j = 0; j < clazzes[i].length; j++) 
//				System.out.println("\t"+ clazzes[i][j].length+ " "+clazzes[i][j][0]);
//		}
//
//		//GraphHandler.writeOut(g);
//		
//		
//		
//		
//			//System.out.print(ge[i].getStableID()+"\t");
////			ASMultiVariation[] as= ge[i].getASComplexes();
////			boolean breakk= false;
////			for (int j = 0; as!= null&& !breakk&& j < as.length; j++) 
////				for (int k = 0; !breakk&& k < as[j].getASVariations().length; k++) 
////					for (int l = 0; !breakk&& l < as[j].getASVariations()[k].getASEvents().length; l++) 
////						if (as[j].getASVariations()[k].getASEvents()[l].isIntronRetention())
////							breakk= true;
////			if (ge[i].getStableID().equals("ENSG00000107643")) {
////				JFrame myFrame= new JFrame(); 
////				myFrame.getContentPane().add(new TranscriptPainter(ge[i]));
////				myFrame.pack();
////				myFrame.setVisible(true);
////			}				
////		}
////		System.out.println(asevents+ "AS Variation");
//		
	}
}
