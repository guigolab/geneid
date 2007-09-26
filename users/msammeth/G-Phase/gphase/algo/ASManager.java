package gphase.algo;

import java.io.PrintStream;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import gphase.NMDSimulator2;
import gphase.io.gtf.GTFObject;
import gphase.model_heavy.ASMultiVariation;
import gphase.model_heavy.ASVariation;
import gphase.model_heavy.ASVariationWithRegions;
import gphase.model_heavy.DirectedRegion;
import gphase.model_heavy.Gene;
import gphase.model_heavy.Transcript;
import gphase.model_heavy.Translation;
import gphase.tools.Arrays;
import gphase.tools.Formatter;

public class ASManager {
	
	public static int maxTranscripts= 500;
	
	static public class AbundancyComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			ASVariation[][] v1= (ASVariation[][]) o1;
			ASVariation[][] v2= (ASVariation[][]) o2;
			
			if (v1.length> v2.length)
				return -1;
			if (v2.length> v1.length)
				return 1;
			
			if (v1[0][0].toStringStructureCode().length()< v2[0][0].toStringStructureCode().length())
				return -1;
			if (v2[0][0].toStringStructureCode().length()< v1[0][0].toStringStructureCode().length())
				return 1;
			return 0;
		}
	}
	
	static public void sortLandscape(ASVariation[][][] vars) {
		if (vars== null)
			return;
		Comparator compi= new AbundancyComparator();
		java.util.Arrays.sort(vars, compi);
	}
	
	static public void printLandscape(ASVariation[][][] vars, PrintStream p) {
		printLandscape(vars, p, null);
	}

	static public void printStats(ASVariation[][][] vars, PrintStream p) {
		int evCnt= 0, clsCnt= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			evCnt+= vars[i].length;
			++clsCnt;
		}
			
		p.print(evCnt+" events in "+clsCnt+" structural different classes found.\n");
	}
	// 1*2^ , 0
	static public void printLandscape(ASVariation[][][] vars, PrintStream p, String detailID) {
		int cntAllEv= Arrays.countFieldsRek(vars, 0, 0, 1);
		//vars= (ASVariation[][][]) Arrays.sortNDFieldRev(vars);
		sortLandscape(vars);
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			p.print(vars[i].length+"\t"+Formatter.fprint(vars[i].length* 100d/cntAllEv, 2)+
					"\t"+vars[i][0][0].toStringStructureCode()+"\n");
			for (int j = 0; detailID!= null&& vars[i][0][0].toStringStructureCode().equals(detailID)&&
					j < vars[i].length; j++) {
				p.print("\t"+vars[i][j][0].toStringUCSCtpair()+"\n");
			}
				
		}
		
	}
	
	static public Gene[] filterNonCodingTranscripts(Gene[] g) {
		Vector gV= new Vector(g.length);
		for (int i = 0; i < g.length; i++) {
			Vector remV= new Vector();
			for (int j = 0; j < g[i].getTranscripts().length; j++) 
				if (!g[i].getTranscripts()[j].isCoding())
					remV.add(g[i].getTranscripts()[j]);
			
			for (int j = 0; j < remV.size(); j++) 
				g[i].removeTranscript((Transcript) remV.elementAt(j));
			if (g[i].getTranscripts().length> 0)
				gV.add(g[i]);
		}
		return (Gene[]) Arrays.toField(gV);
	}
	
	static public Gene[] filterCodingTranscripts(Gene[] g) {
		Vector gV= new Vector(g.length);
		for (int i = 0; i < g.length; i++) {
			Vector remV= new Vector();
			for (int j = 0; j < g[i].getTranscripts().length; j++) 
				if (g[i].getTranscripts()[j].isCoding())
					remV.add(g[i].getTranscripts()[j]);
			
			for (int j = 0; j < remV.size(); j++) 
				g[i].removeTranscript((Transcript) remV.elementAt(j));
			if (g[i].getTranscripts().length> 0)
				gV.add(g[i]);
		}
		return (Gene[]) Arrays.toField(gV);
	}
	
	static public Object[] filterNMDTranscripts(Object[] o, Method[] m) {
		if (o== null|| m== null)
			return null;
		
			// check methods
		Vector v= new Vector(m.length);
		for (int i = 0; i < m.length; i++) {
			if ((!m[i].getReturnType().equals(Boolean.class))|| (m[i].getParameterTypes().length!= 0))
				System.err.println("Invalid filter method "+m[i]);
			v.add(m[i]);
		}
		m= (Method[]) Arrays.toField(v);
		
		v= new Vector();
		for (int i = 0; i < o.length; i++) {
			for (int j = 0; j < m.length; j++) {
				try {
					Boolean b= (Boolean) m[j].invoke(o[i], null);
					if (b.booleanValue())
						v.add(o[i]);
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
			}
		}
		
		return v.toArray();
	}
	
	static public Gene[] filterNMDTranscripts(Gene[] g) {
		Vector gV= new Vector(g.length);
		for (int i = 0; i < g.length; i++) {
			Vector remV= new Vector();
			for (int j = 0; j < g[i].getTranscripts().length; j++) {
				NMDSimulator2 nmd= new NMDSimulator2(g[i].getTranscripts()[j]);
				Translation tln;
				if (g[i].getTranscripts()[j].getTranslations()== null|| 
						g[i].getTranscripts()[j].getTranslations()[0]== null)
					tln= g[i].getTranscripts()[j].findHavanaORF();	// *
				else
					tln= g[i].getTranscripts()[j].getTranslations()[0];
				if (tln!= null&& nmd.isNMD(tln)) 	// (tln== null|| nmd.isNMD(tln)) for removing trpts wo ORF
					remV.add(g[i].getTranscripts()[j]);
			}
			
			for (int j = 0; j < remV.size(); j++) 
				g[i].removeTranscript((Transcript) remV.elementAt(j));
			if (g[i].getTranscripts().length> 0)
				gV.add(g[i]);
		}
		return (Gene[]) Arrays.toField(gV);
	}
	
	static public ASVariation[][][] getASVariations(Gene[] g) {
		return getASVariations(g, ASMultiVariation.FILTER_NONE, null);
	}
	
	static public ASVariation[][][] getASVariations(Gene[] g, HashMap trptIDs) {
		return getASVariations(g, ASMultiVariation.FILTER_NONE, trptIDs);
	}

	static public ASVariation[][][] getASVariations(Transcript[] t, HashMap refTrptIDs) {
		
		if (t.length> maxTranscripts) {
			System.out.println("WARNING: skipped gene with "+t.length+" transcripts: "+DirectedRegion.getUnion(t));
			return null;
		}
		ASVariation[] vars= Gene.getASVariations(t, refTrptIDs);
		ASVariation[][] redVars= ASMultiVariation.clusterIdenticalEvents(vars);
		
		return clusterStructuralEqualEvents(redVars);
		
	}
	
	static public ASVariation[][][] clusterStructuralEqualEvents(ASVariation[][] redVars) {
		HashMap mapStructures= new HashMap();
		for (int j = 0; redVars!= null&& j < redVars.length; j++) {
			Vector v= (Vector) mapStructures.remove(redVars[j][0].toStringStructureCode());
			if (v== null)
				v= new Vector();
			v.add(redVars[j]);
			mapStructures.put(redVars[j][0].toStringStructureCode(), v);		
		}
		return (ASVariation[][][]) Arrays.toField(mapStructures.values());
	}

	static public ASVariation[][][] getASVariations(Gene[] g, int filterCode, HashMap refTrptIDs) {
		
		HashMap mapStructures= new HashMap();
		for (int i = 0; i < g.length; i++) {
			if (g[i].getTranscriptCount()> maxTranscripts) {
				System.out.println("WARNING: skipped gene with "+g[i].getTranscriptCount()+" transcripts: "+g[i].getNameTranscript());
				continue;
			}
			ASVariation[] vars= g[i].getASVariations(filterCode, refTrptIDs);
			ASVariation[][] redVars= ASMultiVariation.clusterIdenticalEvents(vars);
			for (int j = 0; redVars!= null&& j < redVars.length; j++) {
				Vector v= (Vector) mapStructures.remove(redVars[j][0].toStringStructureCode());
				if (v== null)
					v= new Vector();
				v.add(redVars[j]);
				mapStructures.put(redVars[j][0].toStringStructureCode(), v);
			}
		}
		
		return (ASVariation[][][]) Arrays.toField(mapStructures.values());
	}
	
	/**
	 * counts events in dimension 2
	 * @param vars
	 * @return
	 */
	static public int countDifferentEvents(ASVariation[][][] vars) {
		return Arrays.countFieldsRek(vars, 0, 0, 1);
	}
	
	static public ASVariation[][][] filterNonGTAGintrons(ASVariation[][][] vars) {
		Vector result= new Vector();
		for (int i = 0; i < vars.length; i++) {
			Vector keepV= new Vector();
			for (int j = 0; j < vars[i].length; j++) 
				if (vars[i][j][0].hasOnlyGTAGintrons())
					keepV.add(vars[i][j]);
			if (keepV.size()> 0)
				result.add(keepV);
		}
		return (ASVariation[][][]) Arrays.toField(result);
	}

	static public ASVariation[][][] filter(ASVariation[][][] vars, Method m) {
		return filter(vars, m, 2);
	}
	/**
	 * filter such that all variations for that Method m returns <code>true</code> are in the
	 * returned set. 
	 * @param vars
	 * @param m
	 * @param level
	 * @return
	 */
	static public ASVariation[][][] filter(ASVariation[][][] vars, Method m, int level) {
		Vector result= new Vector();
		try {
			for (int i = 0; i < vars.length; i++) {
				Vector keepV= new Vector();
				if (level== 1) {	// only check first level
					for (int j = 0; j < vars[i].length; j++) 
						if (((Boolean) m.invoke(null, vars[i][j][0])).booleanValue())
							keepV.add(vars[i][j]);
				} else {	// level= 2
					if (level!= 2) {
						System.out.println("Error: invalid filter level "+level);
						return vars;
					}
					for (int j = 0; j < vars[i].length; j++) {
						Vector keepV2= new Vector();
						for (int k = 0; k < vars[i][j].length; k++) 
							if (((Boolean) m.invoke(null, vars[i][j][k])).booleanValue())
								keepV2.add(vars[i][j][k]);
						if (keepV2.size()> 0)
							keepV.add(keepV2);
					}
				}
				if (keepV.size()> 0)
					result.add(keepV);
			}
		} catch (Exception e) {
			System.out.println("Error during filtering: "+e.getMessage());
			return vars;
		}
		return (ASVariation[][][]) Arrays.toField(result);
	}
}
