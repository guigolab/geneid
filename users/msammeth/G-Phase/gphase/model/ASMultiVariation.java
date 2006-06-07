/*
 * Created on Mar 2, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.tools.Arrays;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Vector;

/**
 * 
 * 
 * @author msammeth
 */
public class ASMultiVariation implements Serializable {

	public static final int FILTER_NONE= 0;
	public static final int FILTER_HIERARCHICALLY= 1;
	public static final int FILTER_CODING_REDUNDANT= 2;
	
	ASVariation[] asVariations= null;
	
	public ASMultiVariation(ASVariation[] newASVariations) {
		this.asVariations= newASVariations;
	}
	
	ASVariation[] removeRedundancy(ASVariation[] vars) {
		
		Vector v= new Vector(vars.length);
		for (int i = 0; i < vars.length; i++) 
			v.add(vars[i]);
		
		
		Comparator compi= null;
		compi= new ASVariation.CodingHierarchyFilter();
		for (int i = 0; i < v.size(); i++) 
			for (int j = i+1; j < v.size(); j++) {
				int c= compi.compare(v.elementAt(i), v.elementAt(j));
				if (c== 0|| c== 1)
					v.remove(j--);
				if (c== 2) {
					v.remove(i--);
					break;
				}
				// case of (-1), keep both (different in structure!)
			}
		
		return (ASVariation[]) Arrays.toField(v);
	}

	ASVariation[] removeRedundancyCoding(ASVariation[] vars) {
		
		Vector v= new Vector(vars.length);
		for (int i = 0; i < vars.length; i++) 
			v.add(vars[i]);
		
		
		Comparator compi= null;
		compi= new ASVariation.CodingComparator();
		for (int i = 0; i < v.size(); i++) 
			for (int j = i+1; j < v.size(); j++) {
				if (compi.compare(v.elementAt(i), v.elementAt(j))== 0)
					v.remove(j--);
			}
		
		return (ASVariation[]) Arrays.toField(v);
	}
	
	public ASVariation[] getASVariationsAll() {
		return asVariations;
	}
	
		public ASVariation[] getASVariationsHierarchicallyFiltered() {
		return removeRedundancy(asVariations);
	}
		
	public ASVariation[] getASVariations(Transcript transcriptID_1, Transcript transcriptID_2) {
		if (asVariations== null)
			return null;
		
		Vector v= new Vector();
		for (int i = 0; i < asVariations.length; i++) 
			if (asVariations[i].containsTranscripts(transcriptID_1, transcriptID_2))
				v.add(asVariations[i]);
		
		if (v.size()< 1)
			return null;
		ASVariation[] result= new ASVariation[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (ASVariation) v.elementAt(i);
		
		return result;
	}

	public ASVariation[] getASVariationsClusteredCoding() {
		return removeRedundancyCoding(asVariations);
	}


	public String toString() {
		String result= getDegree()+ ": ( ";
		for (int i = 0; i < asVariations.length; i++) {
			result+= asVariations[i]+"; ";
		}
		result= result.substring(0,result.length()- 2)+ " )";
		return result;
	}

	public int getDegree() {
		int deg= 0;
		for (int i = 0; i < asVariations.length; i++) 
			deg+= asVariations[i].getDegree();
		
		return deg;
	}
}
