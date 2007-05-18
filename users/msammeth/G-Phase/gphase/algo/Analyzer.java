package gphase.algo;

import java.awt.Point;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import org.apache.commons.collections.BidiMap;
import org.apache.commons.collections.bidimap.DualHashBidiMap;
import org.freehep.graphicsio.swf.SWFAction.StopSounds;

import sun.security.krb5.internal.crypto.c;

import gphase.Constants;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.Translation;
import gphase.tools.Arrays;
import gphase.tools.Distribution;
import gphase.tools.DoubleVector;
import gphase.tools.Formatter;
import gphase.tools.IntVector;
import gphase.tools.KMP;
import gphase.tools.Sequence;

public class Analyzer {

	static String[] ese= null, ess= null;
	static final String _04_ACCUMULATIVE_COUNT= "accumulativeCount";
	
	public static class VectorSizeComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			Vector v1= (Vector) o1;
			Vector v2= (Vector) o2;
			if (v1.size()< v2.size())
				return -1;
			if (v1.size()> v2.size())
				return 1;
			return 0;
		}
	}

	public static void main(String[] args) {
		
			// parse
		String startSpe= null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-startSpe")&& (i+1< args.length))
				startSpe= args[i+1];
		}
		
		startSpe= "frog";
		
			// get method things
		//_00_checkConsistency();
		String mName= "_04_codonExtendTruncation";
		String outDir= "ASanalysis";
		Class[] sig= new Class[] {Gene.class, PrintStream.class, HashMap.class};
		Method[] m= null;
		try {
			m= new Method[] {
					Analyzer.class.getMethod("_01_generalStatistics", sig),
					Analyzer.class.getMethod(mName, sig)};
		} catch (Exception e) {
			e.printStackTrace();
		}
		if (m== null) 
			System.err.println("Target method "+mName+" not found.");
		else
			AlgoHandler._00_mainLoopEnsemblChrom(Species.SP_NAMES_METAZOA, startSpe, m);
			//AlgoHandler._00_mainLoopHuman(m);

			
	}
	
	public static void outputLandscape(HashMap landscapeMap, PrintStream p) {
		DualHashBidiMap map= new DualHashBidiMap(landscapeMap);
		Object[] values= map.values().toArray();
		Comparator compi= new VectorSizeComparator();
		java.util.Arrays.sort(values, compi);
		BidiMap revMap= map.inverseBidiMap();
		for (int i = values.length- 1; i >= 0; --i) 
			p.println(((Vector) values[i]).size()+ "\t"+ revMap.get(values[i]));
	}
	
	public static String[] getESE() {
		if (ese == null) {
			ese= getMotifs(ASAnalyzer.INPUT_ESE_HAGEN_PENTAMERS);
		}
	
		return ese;
	}

	public static String[] getESS() {
		if (ess == null) {
			ess= getMotifs(ASAnalyzer.INPUT_ESS_HAGEN_PENTAMERS);
		}

		return ess; 
	}
	
	
	public static void _00_checkConsistency() {
		GTFChrReader reader= new GTFChrReader(Constants.getLatestUCSCAnnotation("worm", "RefSeq", null));
//		DualHashBidiMap varMap= getASVariations(reader, ASMultiVariation.FILTER_NONE, -1);
//		BidiMap iMap= varMap.inverseBidiMap();
//		Object[] o= iMap.keySet().toArray();
//		java.util.Arrays.sort(o);
		Vector v= getASVariationCoordinates(reader, ASMultiVariation.FILTER_NONE, -1);
		try {
			PrintStream p= new PrintStream("delme_GTFReader");
			//for (int i = o.length- 1; i >= 0; --i)
			Object[] o= v.toArray();
			java.util.Arrays.sort(o);
			for(int i= 0; i< v.size(); ++i) 
				p.println(o[i]);
			p.flush(); p.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		EncodeWrapper oldReader= new EncodeWrapper(Constants.getLatestUCSCAnnotation("worm", "RefSeq", null));
		Graph g= oldReader.getGraph(false);
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_NONE);
		Arrays.sort2DFieldRev(vars);
		try {
			PrintStream p= new PrintStream("delme_oldReader");
			v= new Vector();
			for (int i = 0; i < vars.length; i++) 
				for (int j = 0; j < vars[i].length; j++) 
					v.add(vars[i][j].toStringCoordinates());
			Object[] o= v.toArray();
			java.util.Arrays.sort(o);
			for (int i = 0; i < o.length; i++) 
				p.println(o[i]);
			p.flush(); p.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}
	

	/**
	 * new version
	 * @return
	 */
	public static ASVariation[][] getASVariations(GTFChrReader reader, int filtType, Method[] filterMethods) {
		
		HashMap classesHash= new HashMap();
		int asVariations= 0;
		reader.setChromosomeWise(true);
		Gene[] ge= null;
		try {
			reader.read();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		ge= reader.getGenes();
		while (ge!= null) {
			for (int i = 0; ge!= null&& i < ge.length; i++) {
				ASMultiVariation[] as= ge[i].getASMultiVariations();
				for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
					ASVariation[] asvars= null;
					if (filtType== ASMultiVariation.FILTER_NONE)	// all
						asvars= as[j].getASVariationsAll();
					else if (filtType== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
						asvars= as[j].getASVariationsHierarchicallyFiltered();
					else if (filtType== ASMultiVariation.FILTER_CODING_REDUNDANT)
						asvars= as[j].getASVariationsClusteredCoding();
					else if (filtType== ASMultiVariation.FILTER_STRUCTURALLY)
						asvars= as[j].getASVariationsStructurallyFiltered();
					
					asVariations+= asvars.length;
					for (int k = 0; k < asvars.length; k++) {	// get pw variations
						String key= asvars[k].toString();
						Vector v= (Vector) classesHash.get(key);
						if (v== null) {
							v= new Vector();
							classesHash.put(key, v);
						}
						v.add(asvars[k]);
					}
				}
			}
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			ge= reader.getGenes();
			System.gc();
		}
		
		Vector vv= new Vector();
		Object[] o= classesHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) 
			vv.add(classesHash.get(o[i]));
		
		ASVariation[][] vars= (ASVariation[][]) Arrays.toField(vv);
		for (int i = 0; filterMethods!= null&& i < filterMethods.length; i++) {
			vars= (ASVariation[][]) Arrays.filter(vars, filterMethods[i]);
		}
		Arrays.sort2DFieldRev(vars);
		
		return vars;
	}

	public static String[][] getASVariations(GTFChrReader reader, int filtType, int regCode) {
		
		HashMap classesHash= new HashMap();
		int asVariations= 0;
		reader.setChromosomeWise(true);
		Gene[] ge= null;
		try {
			reader.read();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		ge= reader.getGenes();
		while (ge!= null) {
			for (int i = 0; ge!= null&& i < ge.length; i++) {
				ASMultiVariation[] as= ge[i].getASMultiVariations();
				for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
					ASVariation[] asvars= null;
					if (filtType== ASMultiVariation.FILTER_NONE)	// all
						asvars= as[j].getASVariationsAll();
					else if (filtType== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
						asvars= as[j].getASVariationsHierarchicallyFiltered();
					else if (filtType== ASMultiVariation.FILTER_CODING_REDUNDANT)
						asvars= as[j].getASVariationsClusteredCoding();
					else if (filtType== ASMultiVariation.FILTER_STRUCTURALLY)
						asvars= as[j].getASVariationsStructurallyFiltered();
					
					asVariations+= asvars.length;
					for (int k = 0; k < asvars.length; k++) {	// get pw variations
						String key= asvars[k].toString();
						Vector v= (Vector) classesHash.get(key);
						if (v== null) {
							v= new Vector();
							classesHash.put(key, v);
						}
						v.add(asvars[k].toString());
					}
				}
			}
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			ge= reader.getGenes();
			System.gc();
		}
		
		Vector vv= new Vector();
		Object[] o= classesHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) 
			vv.add(classesHash.get(o[i]));
		
		String[][] vars= (String[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		
		return vars;
	}

	public static ASVariation[][] getASVariations(Gene[] ge, int filtType, int regCode) {
		HashMap varMap= new HashMap();
		for (int i = 0; i < ge.length; i++) {
			ASVariation[][] vars= getASVariations(ge[i], filtType, regCode);
			for (int j = 0; vars!= null&& j < vars.length; j++) {
				String key= vars[j][0].toString();
				Vector v= (Vector) varMap.get(key);
				if (v== null)
					v= new Vector();
				v= (Vector) Arrays.addAll(v, vars[j]);
				varMap.put(key, v);
			}
		}
		
		Vector vv= new Vector();
		Object[] keys= varMap.keySet().toArray();
		for (int i = 0; i < keys.length; i++) 
			vv.add(varMap.get(keys[i]));
		ASVariation[][] vars= (ASVariation[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		
		return vars;
	}
	public static ASVariation[][] getASVariations(Gene ge, int filtType, int regCode) {
		
		HashMap classesHash= new HashMap();
		int asVariations= 0;
		ASMultiVariation[] as= ge.getASMultiVariations();
		for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
			ASVariation[] asvars= null;
			if (filtType== ASMultiVariation.FILTER_NONE)	// all
				asvars= as[j].getASVariationsAll();
			else if (filtType== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
				asvars= as[j].getASVariationsHierarchicallyFiltered();
			else if (filtType== ASMultiVariation.FILTER_CODING_REDUNDANT)
				asvars= as[j].getASVariationsClusteredCoding();
			else if (filtType== ASMultiVariation.FILTER_STRUCTURALLY)
				asvars= as[j].getASVariationsStructurallyFiltered();
			
			asVariations+= asvars.length;
			for (int k = 0; k < asvars.length; k++) {	// get pw variations
				String key= asvars[k].toString();
				Vector v= (Vector) classesHash.get(key);
				if (v== null) {
					v= new Vector();
					classesHash.put(key, v);
				}
				v.add(asvars[k]);
			}
		}

		
		Vector vv= new Vector();
		Object[] o= classesHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) {
			if (classesHash.get(o[i])== null)
				System.currentTimeMillis();
			vv.add(classesHash.get(o[i]));
		}
		
		ASVariation[][] vars= (ASVariation[][]) Arrays.toField(vv);
		Arrays.sort2DFieldRev(vars);
		
		return vars;
	}

	public static Vector getASVariationCoordinates(GTFChrReader reader, int filtType, int regCode) {
		
		Comparator compi= new ASVariation.SpliceStringComparator();
		Vector v= new Vector();
		int asVariations= 0;
	
		System.out.println("analyzing "+reader.getFileName()+"..");
		Gene[] ge= null;
		reader.read();
		ge= reader.getGenes();
		while (ge!= null) {
			for (int i = 0; ge!= null&& i < ge.length; i++) {
				ASMultiVariation[] as= ge[i].getASMultiVariations();
				for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
					ASVariation[] asvars= null;
					if (filtType== ASMultiVariation.FILTER_NONE)	// all
						asvars= as[j].getASVariationsAll();
					else if (filtType== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
						asvars= as[j].getASVariationsHierarchicallyFiltered();
					else if (filtType== ASMultiVariation.FILTER_CODING_REDUNDANT)
						asvars= as[j].getASVariationsClusteredCoding();
					else if (filtType== ASMultiVariation.FILTER_STRUCTURALLY)
						asvars= as[j].getASVariationsStructurallyFiltered();
					
					asVariations+= asvars.length;
					for (int k = 0; k < asvars.length; k++) 
						v.add(asvars[k].toStringElza());		//toStringElza());
				}
			}
			ge= reader.read();
			System.gc();
		}
		System.out.println("done.");
		
		return v;
	}

	public static String[] getASVariationCoordinates(GTFChrReader reader, String code, int filtType, int regCode) {
		
		Comparator compi= new ASVariation.SpliceStringComparator();
		Vector v= new Vector();
		int asVariations= 0;

		System.out.println("analyzing "+reader.getFileName()+"..");
		Gene[] ge= null;
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		while (ge!= null) {
			for (int i = 0; ge!= null&& i < ge.length; i++) {
				ASMultiVariation[] as= ge[i].getASMultiVariations();
				for (int j = 0; as!= null&& j < as.length; j++) {	// get complexes of variations 
					ASVariation[] asvars= null;
					if (filtType== ASMultiVariation.FILTER_NONE)	// all
						asvars= as[j].getASVariationsAll();
					else if (filtType== ASMultiVariation.FILTER_HIERARCHICALLY) 	// hierarchical
						asvars= as[j].getASVariationsHierarchicallyFiltered();
					else if (filtType== ASMultiVariation.FILTER_CODING_REDUNDANT)
						asvars= as[j].getASVariationsClusteredCoding();
					else if (filtType== ASMultiVariation.FILTER_STRUCTURALLY)
						asvars= as[j].getASVariationsStructurallyFiltered();
					
					asVariations+= asvars.length;
					for (int k = 0; k < asvars.length; k++) 
						if (code== null|| asvars[k].toString().equals(code))
							v.add(asvars[k].toStringElza());		//toStringElza());
				}
			}
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			ge= reader.getGenes();
			System.gc();
		}
		System.out.println("done.");
		String[] result= (String[]) Arrays.toField(v);
		
		return result;
	}

	public static void _03c_exonIntron(GTFChrReader reader, PrintStream p) {
		
		String[] ese= ASAnalyzer.getMotifs(ASAnalyzer.INPUT_ESE_HAGEN_PENTAMERS);
		String[] ess= ASAnalyzer.getMotifs(ASAnalyzer.INPUT_ESS_HAGEN_PENTAMERS);
	
		IntVector inSizesV= new IntVector(), eseDensity80V= new IntVector(), 
			eseDensity40V= new IntVector(), eseDensity15V= new IntVector(), 
			essDensityV= new IntVector(), vStart= new IntVector(), 
			vEnd= new IntVector(), vPTrans= new IntVector(), vPGene= new IntVector();
		
		Gene[] ge= reader.read();
		float avgRemRate= 0f;
		for (int i = 0; ge!= null&& i < ge.length; i++) {
			
			float fracRem= _03c_exonIntron(ge[i], ese, ess, inSizesV, eseDensity80V, 
					eseDensity40V, eseDensity15V, essDensityV, vStart, vEnd, vPTrans, vPGene);
			
			if (fracRem> 0.5f) {
				System.out.println("Removed "+Formatter.fprint(fracRem*100f, 2)+"% non-GTAG exons.");
				p.println("Removed "+Formatter.fprint(fracRem*100f, 2)+" non-GTAG exons.");
			}
			avgRemRate+= fracRem;
			
			ge= reader.read();
		}
		
		System.out.println("Averagely removed "+Formatter.fprint(avgRemRate* 100f, 2)+"% of the exons.");
		System.out.println();
		System.out.println("[Exons]\tmean\tmedian\tstd dev");
		p.println("[Exons]\tmean\tmedian\tstd dev");
		
			// size and sequence
		Distribution dist= new Distribution(inSizesV.toIntArray());
		String[] statStr= dist.toStatString();
		System.out.println("len\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("len\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(eseDensity80V.toIntArray());
		statStr= dist.toStatString();
		System.out.println("eseD80\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("eseD80\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(eseDensity40V.toIntArray());
		statStr= dist.toStatString();
		System.out.println("eseD40\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("eseD40\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(eseDensity15V.toIntArray());
		statStr= dist.toStatString();
		System.out.println("eseD15\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("eseD15\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(essDensityV.toIntArray());
		statStr= dist.toStatString();
		System.out.println("essDens\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("essDens\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(vStart.toIntArray());
		statStr= dist.toStatString();
		System.out.println("start\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("start\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(vEnd.toIntArray());
		statStr= dist.toStatString();
		System.out.println("end\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("end\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(vPTrans.toIntArray());
		statStr= dist.toStatString();
		System.out.println("pTrpt\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("pTrpt\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		dist= new Distribution(vPGene.toIntArray());
		statStr= dist.toStatString();
		System.out.println("pGene\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		p.println("pGene\t"+statStr[0]+"\t"+statStr[1]+"\t"+statStr[2]);
		
		System.out.println();
		p.println();		
	}

	public static void _00_mainLoop(Method[] target, String outDir) {

		File f= new File(outDir);
		if (f.exists()) {
			String[] dir= f.list();
			if (dir.length> 0) {
				System.out.println("Confirm overwrite "+f.getName());
				int x= 0;
				try {
					x = System.in.read();
				} catch (IOException e) {
					e.printStackTrace();
				}
				if (x== 'y')
					f.delete();
				else
					return;
			}
		} else
			f.mkdir();
		
		for (int x = 0; x < target.length; x++) {
			for (int i = 0; i < Species.SP_NAMES_METAZOA.length; i++) {
//				if (!Species.SP_NAMES_METAZOA.equals("worm"))
//					continue;
				System.out.println(Species.SP_NAMES_METAZOA[i]);
				String iname= "graph"+ File.separator+ Species.SP_NAMES_METAZOA[i];
				Species dummySpec= new Species(Species.SP_NAMES_METAZOA[i]);
				File inF= new File("graph"+File.separator+Species.SP_NAMES_METAZOA[i]+"_42_filtDNA");
				if (inF.exists()) {
					System.err.println("File not found: "+inF.getAbsolutePath());
					continue;
				}
				
				try {
					GTFChrReader reader= new GTFChrReader(f.getAbsolutePath());
					PrintStream p= new PrintStream(f.getAbsoluteFile()+ File.separator+
							Species.SP_NAMES_METAZOA+ "_"+ target[x].getName());
					target[x].invoke(null, new Object[] {reader, p});
					p.flush(); p.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}		
		}

	}
	
	
	
	public static void _03c_exonIntron(Gene ge, PrintStream p, HashMap map)  {
		String inSizeID= "inSizes";
		String gcContentID= "gcContent";
		String inNumberID= "inNumber";
		String nonGTAGintrons_ID= "nonGTAGintrons";
		String exSizeID= "exSizes";
		String eseDensity80ID= "eseDensity80";
		String eseDensity40ID= "eseDensity40";
		String eseDensity15ID= "eseDensity15";
		String essDensityID= "essDensity";
		String startID= "start";
		String endID= "end";
		String pTransID= "pTrans";
		String pGeneID= "pGene";
		String nonGTAGexons_ID= "nonGTAGexons";
		String exNumberID= "exNumber";

			// init
		String[] ese= getESE(), ess= getESS();
		IntVector inSizesV= (IntVector) map.get(inSizeID);
		if (inSizesV== null) {
			inSizesV= new IntVector();
			map.put(inSizeID, inSizesV);
		}
		DoubleVector gcContentV= (DoubleVector) map.get(gcContentID);
		if (gcContentV== null) {
			gcContentV= new DoubleVector();
			map.put(gcContentID, gcContentV);
		}
		Integer nonGTAGintrons= (Integer) map.get(nonGTAGintrons_ID);
		if (nonGTAGintrons== null) 
			nonGTAGintrons= new Integer(0);
		Integer intronNumber= (Integer) map.get(inNumberID);
		if (intronNumber== null) 
			intronNumber= new Integer(0);
		
		IntVector exSizesV= (IntVector) map.get(exSizeID);
		if (exSizesV== null) {
			exSizesV= new IntVector();
			map.put(exSizeID, exSizesV);
		}
		IntVector eseDensity80V= (IntVector) map.get(eseDensity80ID);
		if (eseDensity80V== null) {
			eseDensity80V= new IntVector();
			map.put(eseDensity80ID, eseDensity80V);
		}
		IntVector eseDensity40V= (IntVector) map.get(eseDensity40ID);
		if (eseDensity40V== null) {
			eseDensity40V= new IntVector();
			map.put(eseDensity40ID, eseDensity40V);
		}
		IntVector eseDensity15V= (IntVector) map.get(eseDensity15ID);
		if (eseDensity15V== null) {
			eseDensity15V= new IntVector();
			map.put(eseDensity15ID, eseDensity15V);
		}
		IntVector essDensityV= (IntVector) map.get(essDensityID);
		if (essDensityV== null) {
			essDensityV= new IntVector();
			map.put(essDensityID, essDensityV);
		}
		IntVector vStart= (IntVector) map.get(startID);
		if (vStart== null) {
			vStart= new IntVector();
			map.put(startID, vStart);
		}
		IntVector vEnd= (IntVector) map.get(endID);
		if (vEnd== null) {
			vEnd= new IntVector();
			map.put(endID, vEnd);
		}
		IntVector vPTrans= (IntVector) map.get(pTransID);
		if (vPTrans== null) {
			vPTrans= new IntVector();
			map.put(pTransID, vPTrans);
		}
		IntVector vPGene= (IntVector) map.get(pGeneID);
		if (vPGene== null) {
			vPGene= new IntVector();
			map.put(pGeneID, vPGene);
		}
		Integer nonGTAGexons= (Integer) map.get(nonGTAGexons_ID);
		if (nonGTAGexons== null) 
			nonGTAGexons= new Integer(0);
		Integer exonNumber= (Integer) map.get(exNumberID);
		if (exonNumber== null) 
			exonNumber= new Integer(0);
		
		if (ge!= null) {
			final int FLANKSIZE_STADLER= 80;
			final int FLANKSIZE_EDU= 40;
			final int FLANKSIZE_GIL= 15; 
			int cntNonGTAG= 0;
			
				// exons
			Exon[] regs= ge.getExons();
			int ccnt= 0;
			for (int k = 0; regs!= null&& k < regs.length; k++) {
				Exon reg= regs[k];
				int off5= 0, off3= 0;
				if (reg.getAcceptor()!= null&& reg.getAcceptor().isAcceptor()) 	// not TSS exon
					off5= 2;
				if (reg.getDonor()!= null&& reg.getDonor().isDonor()) 	// not TSS exon
					off3= 2;					
				int start= Math.abs(reg.get5PrimeEdge()- off5);
				int end= Math.abs(reg.get3PrimeEdge()+ off3);
				if (start> end) {
					int h= start; start= end; end= h;
				}
				
					// read seq
				String seq= null;
				try {
					seq= Graph.readSequence(reg.getSpecies(), reg.getChromosome(), reg.isForward(),
							start, end);
				} catch (Exception ex) {
					ex.printStackTrace();
					continue;
				}
				if (off5> 0&& (!seq.substring(0, 2).toUpperCase().equals("AG"))) {
					++cntNonGTAG;
					continue;
				}
				if (off3> 0&& (!seq.substring(seq.length()- 2, seq.length()).toUpperCase().equals("GT"))) {
					++cntNonGTAG;
					continue;
				}
				
				
				seq= seq.substring(off5, seq.length()- off3);
				++ccnt;	// new valid exon
				vStart.add(regs[k].get5PrimeEdge()- ge.get5PrimeEdge());
				vEnd.add(ge.get3PrimeEdge()- regs[k].get3PrimeEdge());
				int len= reg.getLength();
				
				int cnt= 0;
				for (int x = 0; x < ess.length; ++x) {	// unbeschr f. ESS
					KMP kmp= new KMP(seq, ess[x], false);
					cnt+= kmp.countMatches();
				}
				essDensityV.add(cnt);
				
					// ESE stadler
				cnt= 0;
				if (len<= (2* FLANKSIZE_STADLER)) {	// flank f. ESE
					for (int x = 0; x < ese.length; ++x) {
						KMP kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity80V.add(cnt); 
				
				} else {	// search the 2 flanks
					String seq0= seq.substring(0, FLANKSIZE_STADLER);
					seq= seq.substring(seq.length()- FLANKSIZE_STADLER, seq.length());
					for (int x = 0; x < ese.length; x++) {
						KMP kmp= new KMP(seq0, ese[x], false);
						cnt+= kmp.countMatches();
						kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity80V.add(cnt); 
				}
				//System.out.print(cnt+"\t");
				
				// ESE edu
				cnt= 0;
				if (len<= (2* FLANKSIZE_EDU)) {	// flank f. ESE
					for (int x = 0; x < ese.length; ++x) {
						KMP kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity40V.add(cnt); 
				
				} else {	// search the 2 flanks
					String seq0= seq.substring(0, FLANKSIZE_EDU);
					seq= seq.substring(seq.length()- FLANKSIZE_EDU, seq.length());
					for (int x = 0; x < ese.length; x++) {
						KMP kmp= new KMP(seq0, ese[x], false);
						cnt+= kmp.countMatches();
						kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity40V.add(cnt); 
				}
				//System.out.print(cnt+"\t");
				
				// ESE gil
				cnt= 0;
				if (len<= (2* FLANKSIZE_GIL)) {	// flank f. ESE
					for (int x = 0; x < ese.length; ++x) {
						KMP kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity15V.add(cnt); 
				
				} else {	// search the 2 flanks
					String seq0= seq.substring(0, FLANKSIZE_GIL);
					seq= seq.substring(seq.length()- FLANKSIZE_GIL, seq.length());
					for (int x = 0; x < ese.length; x++) {
						KMP kmp= new KMP(seq0, ese[x], false);
						cnt+= kmp.countMatches();
						kmp= new KMP(seq, ese[x], false);
						cnt+= kmp.countMatches();
					}
					eseDensity15V.add(cnt); 				
				}
				//System.out.println(cnt);
				
			}	// iter regions
			// save
			nonGTAGexons= new Integer(nonGTAGexons.intValue()+ cntNonGTAG);
			map.put(nonGTAGexons_ID, nonGTAGexons);
			exonNumber= new Integer(exonNumber.intValue()+ regs.length);
			map.put(exNumberID, exonNumber);
			
			
				// introns
			DirectedRegion[] intrRegs= ge.getIntrons();
			int cntSkipIntrons= 0;
			for (int i = 0; intrRegs!= null&& i < intrRegs.length; i++) {
				inSizesV.add(intrRegs[i].getLength());
				String seq= Graph.readSequence(ge.getSpecies(), ge.getChromosome(), ge.isForward(), 
						intrRegs[i].getStart(), intrRegs[i].getEnd()).toUpperCase();
				if (seq.length()< 4|| 
						(!seq.substring(0,2).equals("GT"))|| (!seq.substring(seq.length()- 2).equals("AG"))) {
					++cntSkipIntrons;
					continue;
				}
				int gc= Sequence.countGC(seq);
				gcContentV.add(((double) gc)/ seq.length());
			} 
			if (intrRegs!= null) {
				intronNumber= new Integer(intronNumber.intValue()+ intrRegs.length);
				map.put(inNumberID, intronNumber);
				nonGTAGintrons= new Integer(nonGTAGintrons.intValue()+ cntSkipIntrons);
				map.put(nonGTAGintrons_ID, nonGTAGintrons);
			}
			
			
				// gene / transcript level
			Transcript[] trpts= ge.getTranscripts();
			for (int j = 0; j < trpts.length; j++) 
				vPTrans.add(trpts[j].getExons().length);
			vPGene.add(regs.length);
	
		
		} 

		if (p!= null) {
			Distribution distPge= new Distribution(vPGene.toIntArray()); 
			Distribution distPtrpt= new Distribution(vPTrans.toIntArray()); 
			Distribution distInsize= new Distribution(inSizesV.toIntArray());
			Distribution distStart= new Distribution(vStart.toIntArray()); 
			Distribution distEnd= new Distribution(vEnd.toIntArray()); 
			Distribution distESE80= new Distribution(eseDensity80V.toIntArray()); 
			Distribution distESE40= new Distribution(eseDensity40V.toIntArray()); 
			Distribution distESE15= new Distribution(eseDensity15V.toIntArray()); 
			Distribution distESS= new Distribution(essDensityV.toIntArray());
			Distribution distIntrSize= new Distribution(inSizesV.toIntArray());
			Distribution distGCcoDistribution= new Distribution(gcContentV.toDoubleArray());
			
			System.out.println("%nonGTAG\tpTrpt\tpGene\tinsize\tstart\tend\tese80\tese40\tese15\tess"+
					"\tintrSize\tgc%\tnonGTAGintr");			
			System.out.println(
					Formatter.fprint((nonGTAGexons.intValue()* 10d)/ exonNumber.intValue(), 2)+ "%\t"+
					distPtrpt.getMedian()+ "\t"+
					distPge.getMedian()+ "\t"+
					distInsize.getMedian()+ "\t"+
					distStart.getMedian()+ "\t"+
					distEnd.getMedian()+ "\t"+
					distESE80.getMedian()+ "\t"+
					distESE40.getMedian()+ "\t"+
					distESE15.getMedian()+ "\t"+
					distESS.getMedian()+ "\t"+
					distIntrSize.getMedian()+ "\t"+
					distGCcoDistribution.getMedian()+ "\t"+
					Formatter.fprint((nonGTAGintrons.intValue()* 10d)/ intronNumber.intValue(), 2)+ "%"
					);
			if (p!= null) {
				p.println("%nonGTAG\tpTrpt\tpGene\tinsize\tstart\tend\tese80\tese40\tese15\tess"+
							"\tintrSize\tgc%\tnonGTAGintr");						
				p.println(
						Formatter.fprint((nonGTAGexons.intValue()* 10d)/ exonNumber.intValue(), 2)+ "%\t"+
						distPtrpt.getMedian()+ "\t"+ 
						distPge.getMedian()+ "\t"+
						distInsize.getMedian()+ "\t"+
						distStart.getMedian()+ "\t"+
						distEnd.getMedian()+ "\t"+
						distESE80.getMedian()+ "\t"+
						distESE40.getMedian()+ "\t"+
						distESE15.getMedian()+ "\t"+
						distESS.getMedian()+ "\t"+
						distIntrSize.getMedian()+ "\t"+
						distGCcoDistribution.getMedian()+ "\t"+
						Formatter.fprint((nonGTAGintrons.intValue()* 10d)/ intronNumber.intValue(), 2)+ "%"
						);
			}

		}
		
		
	}

	public static String[] getMotifs(String fName) {
		Vector v= null;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			v= new Vector();
			while (buffy.ready()) {
				String line= buffy.readLine().trim();
				v.add(line);
			}
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		//System.out.println("read "+ v.size()+ " motifs from "+fName);
		return ((String[]) Arrays.toField(v));
	}

	/**
	 * based on exons
	 * @param g
	 * @param coding
	 */
	public static void _01_generalStatistics(Gene ge, PrintStream p, HashMap map) {
		
			final String cntSEGene_ID= "cntSEGene";	// single exon genes 
			final String cntSTGene_ID= "cntSTGene";	// single trpt genes 
			final String cntSTEGene_ID= "cntSTEGene";	// genes wo AS 
			final String cntPlusGene_ID= "cntPlusGene";	// single exon genes 
			final String cntMinusGene_ID= "cntMinusGene";	// single exon genes 
			final String trptGene_ID= "trptGene"; 
			final String exonGene_ID= "exonGene"; 
			final String exonTrpt_ID= "exonTrpt"; 

			Integer cntPlusGene= (Integer) map.get(cntPlusGene_ID);			
			if (cntPlusGene== null) {
				cntPlusGene= new Integer(0);
				map.put(cntPlusGene_ID, cntPlusGene);
			}
			Integer cntMinusGene= (Integer) map.get(cntMinusGene_ID);			
			if (cntMinusGene== null) {
				cntMinusGene= new Integer(0);
				map.put(cntMinusGene_ID, cntMinusGene);
			}
			Integer cntSTGene= (Integer) map.get(cntSTGene_ID);			
			if (cntSTGene== null) {
				cntSTGene= new Integer(0);
				map.put(cntSTGene_ID, cntSTGene);
			}
			Integer cntSEGene= (Integer) map.get(cntSEGene_ID);			
			if (cntSEGene== null) {
				cntSEGene= new Integer(0);
				map.put(cntSEGene_ID, cntSEGene);
			}
			Integer cntSTEGene= (Integer) map.get(cntSTEGene_ID);			
			if (cntSTEGene== null) {
				cntSTEGene= new Integer(0);
				map.put(cntSTEGene_ID, cntSTEGene);
			}
			IntVector trptGeneV= (IntVector) map.get(trptGene_ID);
			if (trptGeneV== null) {
				trptGeneV= new IntVector();
				map.put(trptGene_ID, trptGeneV);
			}
			IntVector exonGeneV= (IntVector) map.get(exonGene_ID);
			if (exonGeneV== null) {
				exonGeneV= new IntVector();
				map.put(exonGene_ID, exonGeneV);
			}
			IntVector exonTrptV= (IntVector) map.get(exonTrpt_ID);
			if (exonTrptV== null) {
				exonTrptV= new IntVector();
				map.put(exonTrpt_ID, exonTrptV);
			}

			
			
			
			if (ge!= null) {
				trptGeneV.add(ge.getTranscriptCount());
				exonGeneV.add(ge.getExons().length);
				for (int i = 0; i < ge.getTranscripts().length; i++) 
					exonTrptV.add(ge.getTranscripts()[i].getExons().length);
				if (ge.getExons().length== 1) {
					cntSEGene= new Integer(cntSEGene.intValue()+ 1);
					map.put(cntSEGene_ID, cntSEGene);
					if (ge.getTranscriptCount()== 1) {
						cntSTEGene= new Integer(cntSTEGene.intValue()+ 1);
						map.put(cntSTEGene_ID, cntSTEGene);
					}
				}
				if (ge.getTranscriptCount()== 1) {
					cntSTGene= new Integer(cntSTGene.intValue()+ 1);
					map.put(cntSTGene_ID, cntSTGene);
				}
				if (ge.isForward()) {
					cntPlusGene= new Integer(cntPlusGene.intValue()+ 1);
					map.put(cntPlusGene_ID, cntPlusGene);
				} else {
					cntMinusGene= new Integer(cntMinusGene.intValue()+ 1);
					map.put(cntMinusGene_ID, cntMinusGene);
				}
				
			}
	
			
			
			if (p!= null) {
				Distribution distTrptGene= new Distribution(trptGeneV.toIntArray());
				Distribution distExonGene= new Distribution(exonGeneV.toIntArray());
				Distribution distExonTrpt= new Distribution(exonTrptV.toIntArray());
				HashMap histo;
				Object[] keys;
				
				int[] exonCnts= exonGeneV.toIntArray();
				int cntExon= 0;
				for (int i = 0; i < exonCnts.length; i++) 
					cntExon+= exonCnts[i];
				
				System.out.println("#genes\t+\t-\tST\tSE\tSTE\t#transcripts\t#exons");
				p.println("#genes\t+\t-\tST\tSE\tSTE\t#transcripts\t#exons");
				System.out.println(trptGeneV.length+"\t"+cntPlusGene+"\t"+cntMinusGene+"\t"+cntSTGene+"\t"+cntSEGene+"\t"+cntSTEGene+"\t"+
						exonTrptV.length+"\t"+cntExon);
				p.println(trptGeneV.length+"\t"+cntPlusGene+"\t"+cntMinusGene+"\t"+cntSTGene+"\t"+cntSEGene+"\t"+cntSTEGene+"\t"+
						exonTrptV.length+"\t"+cntExon);
				System.out.println("med(trpt/gene)\tmed(exons/gene)\tmed(exons/trpt)"+
						"\tmea(trpt/gene)\tmea(exons/gene)\tmea(exons/trpt)");
				p.println("med(trpt/gene)\tmed(exons/gene)\tmed(exons/trpt)"+
						"\tmea(trpt/gene)\tmea(exons/gene)\tmea(exons/trpt)");
				System.out.println(distTrptGene.getMedian()+"\t"+distExonGene.getMedian()+"\t"+distExonTrpt.getMedian()+"\t"+
						distTrptGene.getMean()+"\t"+distExonGene.getMean()+"\t"+distExonTrpt.getMean());
				p.println(distTrptGene.getMedian()+"\t"+distExonGene.getMedian()+"\t"+distExonTrpt.getMedian()+"\t"+
						distTrptGene.getMean()+"\t"+distExonGene.getMean()+"\t"+distExonTrpt.getMean());
			}
			
	
		}

	/**
	 * based on exons
	 * @param g
	 * @param coding
	 */
	public static void _01_lengthDistribution(Gene ge, PrintStream p, HashMap map) {
		
			final String cntSEGene_ID= "cntSEGene";	// single exon genes 
			final String cntSTGene_ID= "cntSTGene";	// single exon genes 
			final String cntPlusGene_ID= "cntPlusGene";	// single exon genes 
			final String cntMinusGene_ID= "cntMinusGene";	// single exon genes 
			final String trptGene_ID= "trptGene"; 
			final String exonGene_ID= "exonGene"; 
			final String exonTrpt_ID= "exonTrpt"; 
			final String ID_GENOMIC_LEN= "genmoicLen";
			final String ID_EXONIC_LEN= "exonicLen";
			final String ID_INTRONIC_LEN= "intronicLen";
			final String ID_GENOMIC_LEN_CDS= "genmoicLenCDS";
			final String ID_EXONIC_LEN_CDS= "exonicLenCDS";
			final String ID_INTRONIC_LEN_CDS= "intronicLenCDS";
			final String ID_GENOMIC_LEN_5UTR= "genmoicLen5UTR";
			final String ID_EXONIC_LEN_5UTR= "exonicLen5UTR";
			final String ID_INTRONIC_LEN_5UTR= "intronicLen5UTR";
			
				// length distributions
			IntVector genLenV= (IntVector) map.get(ID_GENOMIC_LEN);
			if (genLenV== null) {
				genLenV= new IntVector();
				map.put(ID_GENOMIC_LEN, genLenV);
			}
			IntVector exLenV= (IntVector) map.get(ID_EXONIC_LEN);
			if (exLenV== null) {
				exLenV= new IntVector();
				map.put(ID_EXONIC_LEN, exLenV);
			}
			IntVector inLenV= (IntVector) map.get(ID_INTRONIC_LEN);
			if (inLenV== null) {
				inLenV= new IntVector();
				map.put(ID_INTRONIC_LEN, inLenV);
			}
			IntVector genLenCDSV= (IntVector) map.get(ID_GENOMIC_LEN_CDS);
			if (genLenCDSV== null) {
				genLenCDSV= new IntVector();
				map.put(ID_GENOMIC_LEN_CDS, genLenCDSV);
			}
			IntVector exLenCDSV= (IntVector) map.get(ID_EXONIC_LEN_CDS);
			if (exLenCDSV== null) {
				exLenCDSV= new IntVector();
				map.put(ID_EXONIC_LEN_CDS, exLenCDSV);
			}
			IntVector inLenCDSV= (IntVector) map.get(ID_INTRONIC_LEN_CDS);
			if (inLenCDSV== null) {
				inLenCDSV= new IntVector();
				map.put(ID_INTRONIC_LEN_CDS, inLenCDSV);
			}
			IntVector genLen5UTRV= (IntVector) map.get(ID_GENOMIC_LEN_5UTR);
			if (genLen5UTRV== null) {
				genLen5UTRV= new IntVector();
				map.put(ID_GENOMIC_LEN_5UTR, genLen5UTRV);
			}
			IntVector exLen5UTRV= (IntVector) map.get(ID_EXONIC_LEN_5UTR);
			if (exLen5UTRV== null) {
				exLen5UTRV= new IntVector();
				map.put(ID_EXONIC_LEN_5UTR, exLen5UTRV);
			}
			IntVector inLen5UTRV= (IntVector) map.get(ID_INTRONIC_LEN_5UTR);
			if (inLen5UTRV== null) {
				inLen5UTRV= new IntVector();
				map.put(ID_INTRONIC_LEN_5UTR, inLen5UTRV);
			}
			
			
			if (ge!= null) {
				genLenV.add(ge.getLength());
				for (int i = 0; i < ge.getTranscripts().length; i++) {
					exLenV.add(ge.getTranscripts()[i].getExonicLength());
					inLenV.add(ge.getTranscripts()[i].getIntronicLength());
					
				}
				
			}

			
			
			if (p!= null) {
				Distribution distTrptGene= new Distribution(trptGeneV.toIntArray());
				Distribution distExonGene= new Distribution(exonGeneV.toIntArray());
				Distribution distExonTrpt= new Distribution(exonTrptV.toIntArray());
				HashMap histo;
				Object[] keys;
				
				int[] exonCnts= exonGeneV.toIntArray();
				int cntExon= 0;
				for (int i = 0; i < exonCnts.length; i++) 
					cntExon+= exonCnts[i];
				
				System.out.println("#genes\t+\t-\tST\tSE\t#transcripts\t#exons\t#STgene\tmed(trpt/gene)\tmed(exons/gene)\tmed(exons/trpt)"+
						"\tmea(trpt/gene)\tmea(exons/gene)\tmea(exons/trpt)");
				p.println("#genes\t+\t-\tST\tSE\t#transcripts\t#exons\t#STgene\tmed(trpt/gene)\tmed(exons/gene)\tmed(exons/trpt)");
				System.out.println(trptGeneV.length+"\t"+cntPlusGene+"\t"+cntMinusGene+"\t"+cntSTGene+"\t"+cntSEGene+"\t"+exonTrptV.length+"\t"+cntExon+"\t"+
						distTrptGene.getMedian()+"\t"+distExonGene.getMedian()+"\t"+distExonTrpt.getMedian()+"\t"+
						distTrptGene.getMean()+"\t"+distExonGene.getMean()+"\t"+distExonTrpt.getMean());
			}
			

		}

	/**
	 * based on exons
	 * @param g
	 * @param coding
	 */
	public static void _04_codonExtendTruncation(Gene ge, PrintStream p, HashMap map) {
		
			//Graph g, boolean coding, boolean accumulate);
		 	final int DELTA= 100;	//10; 5= max to have a full codon in utr
			final int EXONIC_POS= 1;
			String[] stops= Translation.STOP_CODONS;
			final String accCodonV_ID= "accCodonV"; 
			final String donCodonV_ID= "donCodonV"; 
			
			IntVector accCodonV= (IntVector) map.get(accCodonV_ID);
			if (accCodonV== null) {
				accCodonV= new IntVector();
				map.put(accCodonV_ID, accCodonV);
			}
			IntVector donCodonV= (IntVector) map.get(donCodonV_ID);
			if (donCodonV== null) {
				donCodonV= new IntVector();
				map.put(donCodonV_ID, donCodonV);
			}
			
			
			
			if (ge!= null) {
				Exon[] ex= ge.getExons();
				for (int i = 0; i < ex.length; i++) {
					if (ex[i].getAcceptor()!= null&& ex[i].isCodingSomewhere5Prime()) {
						String seq= Graph.readSequence(ge.getSpecies(), ge.getChromosome(), ge.isForward(),
								ex[i].get5PrimeEdge()+ EXONIC_POS, ex[i].get5PrimeEdge()- DELTA);
						int frame= ex[i].getFrame();
						int x= seq.length()- (2- frame)- 3;	// EXONIC_POS contained in seq.length()
						int cod;
						for (cod= 1; x >= 0; x-= 3, ++cod) {
							String codon= seq.substring(x, x+3);
							int j;
							for (j = 0; j < stops.length; j++) 
								if (codon.toUpperCase().equals(stops[j]))
									break;
							if (j< stops.length) {
								accCodonV.add(cod);
								break;				
							}
						}
						
							
					}
					if (ex[i].getDonor()!= null&& ex[i].isCoding3Prime()) {
						String seq= Graph.readSequence(ge.getSpecies(), ge.getChromosome(), ge.isForward(),
								ex[i].get3PrimeEdge()- EXONIC_POS, ex[i].get3PrimeEdge()+ DELTA);
						int frame= ex[i].get3PrimeFrame();
						int x= (frame+ EXONIC_POS)% 3;
						int cod;
						for (cod= 1; (x+ 3) < seq.length(); x+= 3, ++cod) {
							String codon= seq.substring(x, x+3);
							int j;
							for (j = 0; j < stops.length; j++) 
								if (codon.toUpperCase().equals(stops[j]))
									break;
							if (j< stops.length) {
								donCodonV.add(cod);
								break;				
							}
						}
					}
				}
			}
	
			
			
			if (p!= null) {
				Distribution dist= new Distribution(accCodonV.toIntArray());
				HashMap histo= dist.getHistogram();
				Distribution dist2= new Distribution(donCodonV.toIntArray());
				HashMap histo2= dist2.getHistogram();
				System.out.println("- acceptors -\t\t\t\t- donors");
				System.out.println("#codon\ttrunc\taccum\t\t#codon\ttrunc\taccum");
				p.println("- acceptors -\t\t\t\t- donors");
				p.println("#codon\ttrunc\taccum\t\t#codon\ttrunc\taccum");
				
				Object[] keys= null;
				if (histo.keySet().size()> histo2.keySet().size())
					keys= histo.keySet().toArray();
				else
					keys= histo2.keySet().toArray();
				
				java.util.Arrays.sort(keys);
				int sumAcc= 0, sumDon= 0;
				for (int i = 0; i < keys.length; i++) {
					Integer valAcc= (Integer) histo.get(keys[i]);
					if (valAcc== null) {
						System.out.print(keys[i]+ "\t0\t"+ sumAcc+"\t");
						p.print(keys[i]+ "\t0\t"+ sumAcc+"\t");
					} else {
						sumAcc+= valAcc.intValue();
						System.out.print(keys[i]+ "\t"+valAcc.intValue()+"\t"+ sumAcc+"\t\t");
						p.print(keys[i]+ "\t"+valAcc.intValue()+"\t"+ sumAcc+"\t\t");
					}
					Integer valDon= (Integer) histo2.get(keys[i]);
					if (valDon== null) {
						System.out.println(keys[i]+ "\t0\t"+ sumDon+"\t");
						p.println(keys[i]+ "\t0\t"+ sumDon+"\t");
					} else {
						sumDon+= valDon.intValue();
						System.out.println(keys[i]+ "\t"+valDon.intValue()+"\t"+ sumDon);
						p.println(keys[i]+ "\t"+valDon.intValue()+"\t"+ sumDon);
					}
					
				}
				
			}
			
	
		}
}
