package gphase.graph.alignment;

import gphase.graph.SpliceBubble;
import gphase.graph.SpliceEdge;
import gphase.graph.SpliceGraph;
import gphase.graph.SpliceNode;
import gphase.graph.SplicePath;
import gphase.model.AbstractSite;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

public class Mapping {

	public static class PriorityComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			if (((Mapping) arg0).getCost()< ((Mapping) arg1).getCost())
				return -1;
			if (((Mapping) arg0).getCost()> ((Mapping) arg1).getCost())
				return 1;
			return 0;
		}
	}
	
	double cost= 0;
	HashMap mapTableI= new HashMap(), mapTableJ= new HashMap();
	SpliceNode maxI= null, maxJ= null;
	int maxRel= -1;
	SpliceGraph g1, g2;
	
	
	public Mapping(SpliceGraph newG1, SpliceGraph newG2) {
		this.g1= newG1;
		this.g2= newG2;
	}
	
	protected Object clone() throws CloneNotSupportedException {
		
		Mapping map= new Mapping(g1,g2);
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
			if (nJ== null&& nI.getOutDegree()== 0)
				cost+= getCost_edge(nI, nJ); // getCost_pathes(nI, nJ);
		}
		if (nJ!= null) {
			maxJ= nJ;
			maxRel= 2;
			if (nI== null&& nJ.getOutDegree()== 0)
				cost+= getCost_edge(nI, nJ); // getCost_pathes(nI, nJ);
		}
		
		if (nI!= null&& nJ!= null) {
			mapTableI.put(nI, nJ);
			mapTableJ.put(nJ, nI);
			maxRel= 3;
			cost+= getCost_edge(nI, nJ); // getCost_pathes(nI, nJ);
		}
		
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

	SplicePath[] getShortestPathes(SplicePath[] pathes) {
			// find pathes sharing at least one transcript
		Vector remV= new Vector();
		for (int i = 0; i < pathes.length; i++) {
			for (int j = 0; j < pathes.length; j++) {
				if (i== j)
					continue;
				if (SpliceBubble.intersects(pathes[i].getTranscripts(), pathes[j].getTranscripts())) {
					if (pathes[i].edgeLength()> pathes[j].edgeLength())
						remV.add(pathes[i]);
					else
						remV.add(pathes[j]);
				}
			}
		}
		
			// remove longer ones
		Vector v= new Vector();
		for (int i = 0; i < pathes.length; i++) {
			int j;
			for (j = 0; j < remV.size(); j++) 
				if (pathes[i]== remV.elementAt(j))
					break;
			if (j== remV.size())
				v.add(pathes[i]);
		}
		
		
		return (SplicePath[]) gphase.tools.Arrays.toField(v);
	}
	double getCost(SplicePath path1, SplicePath path2) {
		SpliceEdge[] edges= path1.getEdges();
		int exLen1= 0;
		int inNb1= 0;
		for (int i = 0; i < edges.length; i++) {
			if (edges[i].isExonic())
				exLen1+= edges[i].getLength();
			else
				++inNb1;
		}
		
		edges= path1.getEdges();
		int exLen2= 0;
		int inNb2= 0;
		for (int i = 0; i < edges.length; i++) {
			if (edges[i].isExonic())
				exLen2+= edges[i].getLength();
			else
				++inNb2;
		}
		
		double exID, inID;
		if (exLen1== 0|| exLen2== 0)
			exID= 1d;
		else
			exID= 1- Math.min((double) exLen1/ exLen2, (double) exLen2/ exLen1);
		
		if (inNb1== 0|| inNb2== 0)
			inID= 1d;
		else
			inID= 1- Math.min((double) inNb1/ inNb2, (double) inNb2/ inNb1);
		
		return (exID+ inID);
	}
	
	double getCost_pathes(SpliceNode nI, SpliceNode nJ) {
		
			// get sources and targets of delimited pathes
		SpliceNode[] srcI= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableI.keySet());
		gphase.tools.Arrays.addAll(srcI, g1.getRoots());
		SpliceNode[] srcJ= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableJ.keySet());
		gphase.tools.Arrays.addAll(srcJ, g2.getRoots());
		SpliceNode[] tgtI= (nI!= null)? new SpliceNode[] {nI}:g1.getLeafs();
		SpliceNode[] tgtJ= (nJ!= null)? new SpliceNode[] {nJ}:g2.getLeafs();
		
		
			// get pathes in both graphs
		SplicePath[] pathes1= g1.getPathes(srcI, tgtI);
		pathes1= getShortestPathes(pathes1);
		SplicePath[] pathes2= g2.getPathes(srcJ, tgtJ);
		pathes2= getShortestPathes(pathes2);
		
			// match pathes
		if (pathes1.length< pathes2.length) {	// swap, outer loop over longer path
			SplicePath[] h= pathes1;
			pathes1= pathes2;
			pathes2= h;
		}
		double totCosts= 0d;
		for (int i = 0; i < pathes1.length; i++) {
			double min= Double.MAX_VALUE;
			for (int j = 0; j < pathes2.length; j++) {
				double cost= getCost(pathes1[i], pathes2[j]);
				if (cost< min)
					min= cost;
			}
			totCosts+= cost;
		}
		
		return totCosts;
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

	public double getCost() {
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
