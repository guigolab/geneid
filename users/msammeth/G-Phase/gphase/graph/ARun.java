package gphase.graph;

import gphase.algo.ASAnalyzer;
import gphase.model.ASMultiVariation;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.SpliceSite;

import java.io.PrintStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class ARun {
	public static void main(String[] args) {
		outputSeparationIndex();
	}
	
	public static void outputSeparationIndex() {
		PrintStream p= null;
		p= System.out;
		
		Graph gMod= ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE);
		Gene[] ge= gMod.getGenes();
		HashMap sepHash= new HashMap();
		for (int i = 0; i < ge.length; i++) {
			HashMap tmpSep= ge[i].getSeparationIndex();
			Object[] keys= tmpSep.keySet().toArray();
			for (int j = 0; j < keys.length; j++) {
				Vector v= (Vector) sepHash.remove(keys[j]);
				if (v== null)
					v= new Vector();
				v.addAll((Vector) tmpSep.get(keys[j]));
				sepHash.put(keys[j], v);
			}
		}
		
		Object[] keys= sepHash.keySet().toArray();
		Arrays.sort(keys);
		for (int i = keys.length- 1; i >= 0; --i) {
			Vector v= (Vector) sepHash.get(keys[i]);
			for (int j = 0; j < v.size(); j++) {
				ASMultiVariation[] mvars= (ASMultiVariation[]) v.elementAt(j);
				String mv0= gphase.tools.Arrays.toString(mvars[0].getTranscriptsFromHash());
				String mv1= gphase.tools.Arrays.toString(mvars[1].getTranscriptsFromHash());				
				SpliceSite site0= mvars[0].getMinSpliceSite();
				SpliceSite site1= mvars[1].getMinSpliceSite();
				p.println(keys[i]+ "\t"+ mvars[0]+ mv0+ site0.getPos()+ 
						"\n"+ mvars[1]+ mv1+ site1.getPos()+ "\n");
			}
		}
	}
}
