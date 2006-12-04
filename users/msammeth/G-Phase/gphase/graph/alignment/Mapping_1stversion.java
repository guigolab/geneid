package gphase.graph.alignment;

import gphase.graph.SpliceEdge;
import gphase.graph.SpliceNode;
import gphase.model.AbstractSite;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class Mapping_1stversion {

	public static class PriorityComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			if (((Mapping_1stversion) arg0).getCost()< ((Mapping_1stversion) arg1).getCost())
				return -1;
			if (((Mapping_1stversion) arg0).getCost()> ((Mapping_1stversion) arg1).getCost())
				return 1;
			return 0;
		}
	}
	
	int cost= 0;
	HashMap mapTableI= new HashMap(), mapTableJ= new HashMap();
	SpliceNode maxI= null, maxJ= null;
	int maxRel= -1;
	
	public Mapping_1stversion() {
	}
	
	protected Object clone() throws CloneNotSupportedException {
		
		Mapping_1stversion map= new Mapping_1stversion();
		map.cost= getCost();
		map.mapTableI= (HashMap) mapTableI.clone();
		map.mapTableJ= (HashMap) mapTableJ.clone();
		map.maxI= getMaxI();
		map.maxJ= getMaxJ();
		map.maxRel= maxRel;
		
		return map;
	}
	
	public void addMapping(SpliceNode nI, SpliceNode nJ) {
		
		if (nI!= null) {
			maxI= nI;
			maxRel= 1;
		}
		if (nJ!= null) {
			maxJ= nJ;
			maxRel= 2;
		}
		
		if (nI!= null&& nJ!= null) {
			mapTableI.put(nI, nJ);
			mapTableJ.put(nJ, nI);
			maxRel= 3;
		}
		
		cost+= getCost_weight(nI, nJ);
	}
	
	int getCost_edge(SpliceNode nI, SpliceNode nJ) {		
		if (nI== null) 
			return nJ.getInDegree();
		if (nJ== null)
			return nI.getInDegree();
		
			// matched node, compare hits of (mis-)aligned splices
		int ctr= 0;
		SpliceEdge[]  edgesI= nI.getInEdges();
		SpliceEdge[]  edgesJ= nJ.getInEdges();
		for (int i = 0; edgesI!= null&& i < edgesI.length; i++) {	// check predecessors of i
			Object o= mapTableI.get(edgesI[i].getTail());
			if (o== null) {
				++ctr;
				//continue;
			} else {
				int j;
				for (j = 0; j < edgesJ.length; j++) 
					if (edgesJ[j].getTail()== o)	// aligned with a predecesseor of j
						break;
				if (j>= edgesJ.length)
					++ctr;
			}
		}
		for (int i = 0; edgesJ!= null&& i < edgesJ.length; i++) {	// check predecessors of j
			Object o= mapTableJ.get(edgesJ[i].getTail());
			if (o== null) {
				++ctr;
				//continue;
			} else {
				int j;
				for (j = 0; j < edgesI.length; j++) 
					if (edgesI[j].getTail()== o)	// aligned with a predecesseor of i
						break;
				if (j>= edgesI.length)
					++ctr;
			}
		}
		
		return ctr;
	}

	int getCost_weight(SpliceNode nI, SpliceNode nJ) {		
		int fac= 100; 	// factor, to not make it double
		if (nI== null) 
			return nJ.getInDegree()* fac;
		if (nJ== null)
			return nI.getInDegree()* fac;
		
			// matched node, compare hits of (mis-)aligned splices
		int ctr= 0;
		SpliceEdge[]  edgesI= nI.getInEdges();
		SpliceEdge[]  edgesJ= nJ.getInEdges();
		for (int i = 0; edgesI!= null&& i < edgesI.length; i++) {	// check predecessors of i
			Object o= mapTableI.get(edgesI[i].getTail());
			if (o== null) {
				ctr+= fac;
				//continue;
			} else {
				int j;
				for (j = 0; j < edgesJ.length; j++) 
					if (edgesJ[j].getTail()== o)	// aligned with a predecesseor of j
						break;
				if (j>= edgesJ.length) {
					ctr+= fac;
				} else {	// penalty for incongruences in edge length
					int distI= nI.getSite().getPos()- edgesI[i].getTail().getSite().getPos();
					int distJ= nJ.getSite().getPos()- edgesJ[j].getTail().getSite().getPos();
					if (distI<= 0|| distJ<= 0)
						System.err.println("dist< 0");
					double w= 1d- Math.min((double) distI/distJ, (double) distJ/distI);
					ctr+= w* fac;
				}
			}
		}
		for (int i = 0; edgesJ!= null&& i < edgesJ.length; i++) {	// check predecessors of j
			Object o= mapTableJ.get(edgesJ[i].getTail());
			if (o== null) {
				ctr+= fac;
			} else {
				int j;
				for (j = 0; j < edgesI.length; j++) 
					if (edgesI[j].getTail()== o)	// aligned with a predecesseor of i
						break;
				if (j>= edgesI.length) {
					ctr+= fac;
				} else {	// penalty for incongruences in edge length
					int distI= nI.getSite().getPos()- edgesI[j].getTail().getSite().getPos();
					int distJ= nJ.getSite().getPos()- edgesJ[i].getTail().getSite().getPos();
					if (distI<= 0|| distJ<= 0)
						System.err.println("dist< 0");
					double w= 1d- Math.min((double) distI/distJ, (double) distJ/distI);
					ctr+= w* fac;
				}
			}
		}
		
		return ctr;
	}

	int getCost_ancestorPath(SpliceNode nI, SpliceNode nJ) {		
		if (nI== null)
			return nJ.getFromList().size();
		if (nJ== null)
			return nI.getFromList().size();
		
			// matched node, compare hits of closed bubbles
		int ctr= 0;
		Set keysI= nI.getFromList().keySet();
		Set keysJ= nJ.getFromList().keySet();
		Iterator iter= keysI.iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (mapTableI.get(o)== null) {
				ctr+= ((Vector) nI.getFromList().get(o)).size();	// all pathes missed
			} else {	// src mapped
				int pI= ((Vector) nI.getFromList().get(o)).size();	// distance= diff in path#
				int pJ= 0;
				if (nJ.getFromList().get(mapTableI.get(o))!= null)
					pJ= ((Vector) nJ.getFromList().get(mapTableI.get(o))).size();
				ctr+= Math.abs(pI- pJ);
			}
		}
		iter= keysJ.iterator();
		while (iter.hasNext()) {
			Object o= iter.next();
			if (mapTableJ.get(o)== null) {
				ctr+= ((Vector) nJ.getFromList().get(o)).size();	// all pathes missed
			}	// differences in mapped src already added
		}
		
		return ctr;
	}

	public int getCost() {
		return cost;
	}

	public HashMap getMapTableI() {
		return mapTableI;
	}

	public HashMap getMapTableJ() {
		return mapTableJ;
	}

	public SpliceNode getMaxI() {
		return maxI;
	}

	public SpliceNode getMaxJ() {
		return maxJ;
	}
	
	public String toString() {
		Object[] o= mapTableI.keySet().toArray();
		SpliceNode[] nodes= new SpliceNode[o.length];
		for (int i = 0; i < nodes.length; i++) 
			nodes[i]= (SpliceNode) o[i];
		Arrays.sort(nodes, new SpliceNode.PositionComparator());
		
		String result= "";
		for (int i = 0; i < nodes.length; i++) 
			result+= "["+ nodes[i].getSite().getPos()+ ","+ 
				((SpliceNode) mapTableI.get(nodes[i])).getSite().getPos()+ "], ";
		
		if (result.length()> 2)
			result= result.substring(0, result.length()- 2);
		result+= ": "+ getCost();
		
		return result;
	}

	public int getMaxRel() {
		return maxRel;
	}
}
