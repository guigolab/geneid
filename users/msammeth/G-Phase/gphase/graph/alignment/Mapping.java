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

	boolean freeShift= false;
	static String costID= "Len";
	
	public static class PriorityComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			if (((Mapping) arg0).getCost()< ((Mapping) arg1).getCost())
				return -1;
			if (((Mapping) arg0).getCost()> ((Mapping) arg1).getCost())
				return 1;
			return 0;
		}
	}
	
	public static class PriorityExtensionComparator extends PriorityComparator {
		public int compare(Object arg0, Object arg1) {
			if (((Mapping) arg0).offIJ< ((Mapping) arg1).offIJ)
				return -1;
			if (((Mapping) arg0).offIJ> ((Mapping) arg1).offIJ)
				return 1;
			return super.compare(arg0, arg1);
		}
	}

	double cost= 0;
	HashMap mapTableI= new HashMap(), mapTableJ= new HashMap();
	SpliceNode maxI= null, maxJ= null;
	int maxRel= -1;
	SpliceGraph g1, g2;
	int offI= 0, offJ= 0, offIJ= 0;
	
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
	
	public double getCost(SpliceNode nI, SpliceNode nJ) {
		if (costID.contains("LenExp"))
			return getCost_exonicLenExp(nI, nJ);	
		if (costID.contains("LenPiece"))
			return getCost_exonicLenPiece(nI, nJ);	//getCost_exonicLenPiece(nI, nJ);
		if (costID.contains("Len"))	// after lenpiece
			return getCost_exonicLen(nI, nJ);	//getCost_exonicLen
		return getCost_edge(nI, nJ);
	}
	
	public void addMapping(SpliceNode nI, SpliceNode nJ) {
		
		cost+= getCost(nI, nJ); // getCost_pathes(nI, nJ);
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

	int getCost_exonicEdge(SpliceNode nI, SpliceNode nJ) {

		if ((nI!= null&& nI.isAcceptor())|| (nJ!= null&& nJ.isAcceptor()))
			return 0;	// do not penalize intronic edges
		
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

	double getCost_exonicLen_edges(SpliceNode nI, SpliceNode nJ) {
		double d1= getCost_exonicLen(nI, nJ);
		int i2= getCost_exonicEdge(nI, nJ);
		return (d1+ i2);
	}
	
	strictfp double getCost_exonicLen(SpliceNode nI, SpliceNode nJ) {
		
		assert(!(nI== null&& nJ== null));
			//	allow freedom in the last node iff tss
			// for tss handled by 2 sets of null pathes
		if ((nI== null&& nJ.isTES())|| (nJ== null&& nI.isTES())) {	// end node ??
			return 0d; 	// not aligned node, no penalty
		}
		
			// get new pathes (since last aligned node or root)
		SplicePath[] pathesI= null, pathesJ= null;
		if (nI!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableI.keySet());
			SpliceNode[] roots= g1.getRoots(); 
			pathesI= fpath(roots, src, nI, g1);
		}		
		
		if (nJ!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableJ.keySet());
			SpliceNode[] roots= g2.getRoots(); 
			pathesJ= fpath(roots, src, nJ, g2);
		}
		
			// both no valid pathes (eg, src of both graphs)
		if ((pathesI== null|| pathesI.length== 0)&& (pathesJ== null|| pathesJ.length== 0))
			return 0d;	// no penalty
	
			// only one with valid pathes
		if (pathesI== null|| pathesI.length== 0) {
			if (!freeShift)
				return (double) pathesJ.length;
			int cnt= 0;
			for (int i = 0; i < pathesJ.length; i++) 
				if (pathesJ[i].getSrc().getInDegree()> 0)
					++cnt;		// not good, alignment will go over the corner
			return (double) cnt;
		}
		if (pathesJ== null|| pathesJ.length== 0) { 
			if (!freeShift)
				return (double) pathesI.length;
			int cnt= 0;
			for (int i = 0; i < pathesJ.length; i++) 
				if (pathesI[i].getSrc().getInDegree()> 0)
					++cnt;		// only non-border not aligned exonic areas
			return (double) cnt;
		}
		
			// two with valid pathes, match table
		double[][] exLenDiff= new double[pathesI.length][pathesJ.length];
		for (int i = 0; i < pathesI.length; i++) {
			for (int j = 0; j < pathesJ.length; j++) {
				int vw= Math.abs(pathesI[i].getExonicLength());
				int tu= Math.abs(pathesJ[j].getExonicLength());
				if (vw> 0&& tu> 0)
					exLenDiff[i][j]= 1d-Math.min(((double) tu/ vw), ((double) vw/ tu));
				else if (vw== 0&& tu== 0)
					exLenDiff[i][j]= 0d;
				else {
					assert (vw== 0^ tu== 0);
					exLenDiff[i][j]= 1d;
				}
			}
		}
		
			// add up minima
		double sum= 0d;
		for (int i = 0; i < exLenDiff.length; i++) {
			double min= 2d;
			for (int j = 0; j < exLenDiff[i].length; j++) 
				if (exLenDiff[i][j]< min)
					min= exLenDiff[i][j];
			sum+= min/2;
		}
		for (int i = 0; i < exLenDiff[0].length; i++) {
			double min= 2d;
			for (int j = 0; j < exLenDiff.length; j++) 
				if (exLenDiff[j][i]< min)
					min= exLenDiff[j][i];
			sum+= min/2;
		}
		
		return sum;
	}

	strictfp double getCost_exonicLenExp(SpliceNode nI, SpliceNode nJ) {
		
		if (nI.isTSS()|| nJ.isTSS())
			return 0d; 
		
			// get new pathes (since last aligned node or root)
		SplicePath[] pathesI= null, pathesJ= null;
		if (nI!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableI.keySet());
			SpliceNode[] roots= g1.getRoots(); 
			pathesI= fpath(roots, src, nI, g1);
		}		
		
		if (nJ!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableJ.keySet());
			SpliceNode[] roots= g2.getRoots(); 
			pathesJ= fpath(roots, src, nJ, g2);
		}
		
			// at least one with NO valid pathes (eg, src of both graphs)
		if ((pathesI== null|| pathesI.length== 0)|| (pathesJ== null|| pathesJ.length== 0))
			return 1d;	// not alignable
	
			// two with valid pathes, match table
		double[][] exLenDiff= new double[pathesI.length][pathesJ.length];
		for (int i = 0; i < pathesI.length; i++) {
			for (int j = 0; j < pathesJ.length; j++) {
				int vw= Math.abs(pathesI[i].getExonicLength());
				int tu= Math.abs(pathesJ[j].getExonicLength());
				if (vw> 0&& tu> 0)
					exLenDiff[i][j]= 1d-Math.min(((double) tu/ vw), ((double) vw/ tu));
				else if (vw== 0&& tu== 0)
					exLenDiff[i][j]= 0d;
				else {
					assert (vw== 0^ tu== 0);
					exLenDiff[i][j]= 1d;
				}
			}
		}
		
			// add up minima
		double sum= 0d;
		for (int i = 0; i < exLenDiff.length; i++) {
			double min= 2d;
			for (int j = 0; j < exLenDiff[i].length; j++) 
				if (exLenDiff[i][j]< min)
					min= exLenDiff[i][j];
			sum+= min;
		}
		
		return sum;
	}

	strictfp double getCost_exonicLenPiece(SpliceNode nI, SpliceNode nJ) {
		
		assert(!(nI== null&& nJ== null));
		if ((nI== null&& nJ.getOutDegree()> 0)|| (nJ== null&& nI.getOutDegree()> 0)) {	// end node ??
			return 0d; 	// not aligned node, no penalty
		}
		
			// get new pathes (since last aligned node or root)
		SplicePath[] pathesI= null, pathesJ= null;
		if (nI!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableI.keySet());
			SpliceNode[] roots= g1.getRoots(); 
			pathesI= fpath(roots, src, nI, g1);
		}		
		
		if (nJ!= null) {
			SpliceNode[] src= (SpliceNode[]) gphase.tools.Arrays.toField(mapTableJ.keySet());
			SpliceNode[] roots= g2.getRoots(); 
			pathesJ= fpath(roots, src, nJ, g2);
		}
		
			// both no valid pathes (eg, src of both graphs)
		if ((pathesI== null|| pathesI.length== 0)&& (pathesJ== null|| pathesJ.length== 0))
			return 0d;	// no penalty

			// only one with valid pathes
		if (pathesI== null|| pathesI.length== 0) 
			return (double) pathesJ.length* 2;	// 2 for exLen && piece
		if (pathesJ== null|| pathesJ.length== 0) 
			return (double) pathesI.length* 2;
		
			// two with valid pathes, match table
		double[][] exLenDiff= new double[pathesI.length][pathesJ.length];
		for (int i = 0; i < pathesI.length; i++) {
			for (int j = 0; j < pathesJ.length; j++) {
				int vw= pathesI[i].getExonicLength();
				int tu= pathesJ[j].getExonicLength();
				if (vw> 0&& tu> 0)
					exLenDiff[i][j]= 1d-Math.min(((double) tu/ vw), ((double) vw/ tu));
				else if (vw== 0&& tu== 0)
					exLenDiff[i][j]= 0d;
				else {
					assert (vw== 0^ tu== 0);
					exLenDiff[i][j]= 1d;
				}
				
				int exPiecVW= pathesI[i].getExonicPieces(); 
				int exPiecTU= pathesJ[j].getExonicPieces();
				if (exPiecVW> 0&& exPiecTU> 0)
					exLenDiff[i][j]+= 1d-Math.min(((double) exPiecTU/ exPiecVW), ((double) exPiecVW/ exPiecTU));
				else if (exPiecVW== 0&& exPiecTU== 0)
					exLenDiff[i][j]= 0d;
				else {
					assert (exPiecVW== 0^ exPiecTU== 0);
					exLenDiff[i][j]= 1d;	
				}
			}
		}
		
			// add up minima
		double sum= 0d;
		for (int i = 0; i < exLenDiff.length; i++) {
			double min= 2d;
			for (int j = 0; j < exLenDiff[i].length; j++) 
				if (exLenDiff[i][j]< min)
					min= exLenDiff[i][j];
			sum+= min/2;
		}
		for (int i = 0; i < exLenDiff[0].length; i++) {
			double min= 2d;
			for (int j = 0; j < exLenDiff.length; j++) 
				if (exLenDiff[j][i]< min)
					min= exLenDiff[j][i];
			sum+= min/2;
		}
		
		return sum;
	}
	
	SplicePath[] fpath(SpliceNode[] roots, SpliceNode[] src, SpliceNode n, SpliceGraph g) {
		
		SpliceEdge[] edge= n.getInEdges();
		if (edge== null)
			return null;
		
		if (src!= null)
			Arrays.sort(src, new SpliceNode.PositionComparator());
		Vector v= new Vector();
		for (int i = 0; edge!= null&& i < edge.length; i++) {
			SplicePath[] pathes= null;
			if (src!= null)
				for (int j = src.length- 1; j >= 0; --j) {
					pathes= g.findPathes(src[j], edge[i]);
					if (pathes!= null)
						break;
				}
			if (pathes!= null) {
				for (int j = 0; j < pathes.length; j++) 
					v.add(pathes[j]);
			} else {	// no aligned src found for edge, get from roots || src== null
				Vector vv= new Vector();
				for (int k = 0; k < roots.length; k++) {
					pathes= g.findPathes(roots[k], edge[i]);
					if (pathes!= null)
						for (int j = 0; j < pathes.length; j++) 
							vv.add(pathes[j]);
				}
				pathes= (SplicePath[]) gphase.tools.Arrays.toField(vv);
				// assert(pathes!= null);	// can happen: tgt== root, no inedge, no path.
				for (int j = 0; pathes!= null&& j < pathes.length; j++) 
					v.add(pathes[j]);
			}
		}

		return (SplicePath[]) gphase.tools.Arrays.toField(v);
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
		for (int i = 0; edges!= null&& i < edges.length; i++) {
			if (edges[i].isExonic())
				exLen1+= edges[i].getLength();
			else
				++inNb1;
		}
		
		edges= path2.getEdges();
		int exLen2= 0;
		int inNb2= 0;
		for (int i = 0; edges!= null&& i < edges.length; i++) {
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

	public SpliceNode[] getNextAlignPair(SpliceNode[] listI, SpliceNode[] listJ, int idxI, int idxJ) {
		
		while(true) {
			if (offI< offIJ&& idxI+offI+1< listI.length&& idxJ+offIJ< listJ.length)
				return new SpliceNode[] {listI[idxI+(offI++)], listJ[idxJ+offIJ]};
			if (offJ< offIJ&& idxJ+offJ+1< listJ.length&& idxI+offIJ< listI.length)
				return new SpliceNode[] {listI[idxI+offIJ], listJ[idxJ+(offJ++)]};
			++offIJ;
			offI= 1; offJ= 1;
			if (idxI+ offIJ< listI.length&& idxJ+ offIJ< listJ.length)
				return new SpliceNode[] {listI[idxI+offIJ], listJ[idxJ+offIJ]};
			if (idxI+ offIJ> listI.length&& idxJ+ offIJ> listJ.length)
				return null;
		}
	}
	
	public SpliceNode getMaxI() {
		return maxI;
	}

	public SpliceNode getMaxJ() {
		return maxJ;
	}
	
	public String toStringInverse() {
		Object[] o= mapTableJ.keySet().toArray();
		SpliceNode[] nodes= new SpliceNode[o.length];
		for (int i = 0; i < nodes.length; i++) 
			nodes[i]= (SpliceNode) o[i];
		Arrays.sort(nodes, new SpliceNode.PositionComparator());
		
		String result= "";
		for (int i = 0; i < nodes.length; i++) 
			result+= "["+ nodes[i].getSite().getPos()+ ","+ 
				((SpliceNode) mapTableJ.get(nodes[i])).getSite().getPos()+ "], ";
		
		if (result.length()> 2)
			result= result.substring(0, result.length()- 2);
		result+= ": "+ getCost();
		
		return result;
	}

	public String toStringInverse(SpliceNode[] listI, SpliceNode[] listJ) {
		Object[] o= mapTableJ.keySet().toArray();
		SpliceNode[] nodes= new SpliceNode[o.length];
		for (int i = 0; i < nodes.length; i++) 
			nodes[i]= (SpliceNode) o[i];
		Comparator compi= new SpliceNode.PositionComparator();
		Arrays.sort(nodes, compi);
		
		String result= "";
		for (int i = 0; i < nodes.length; i++) {
			int p1= Arrays.binarySearch(listJ, nodes[i], compi);
			int p2= Arrays.binarySearch(listI, mapTableJ.get(nodes[i]), compi);
			result+= "["+ p1+ ","+ p2 + "], ";
		}
		
		if (result.length()> 2)
			result= result.substring(0, result.length()- 2);
		result+= ": "+ getCost();
		
		return result;
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

	public String toString(SpliceNode[] listI, SpliceNode[] listJ) {
		Object[] o= mapTableI.keySet().toArray();
		SpliceNode[] nodes= new SpliceNode[o.length];
		for (int i = 0; i < nodes.length; i++) 
			nodes[i]= (SpliceNode) o[i];
		Comparator compi= new SpliceNode.PositionComparator();
		Arrays.sort(nodes, compi);
		
		String result= "";
		for (int i = 0; i < nodes.length; i++) { 
			int p1= Arrays.binarySearch(listI, nodes[i], compi);
			int p2= Arrays.binarySearch(listJ, mapTableI.get(nodes[i]), compi);
			result+= "["+ p1+ ","+ p2+ "], ";
		}
		
		if (result.length()> 2)
			result= result.substring(0, result.length()- 2);
		result+= ": "+ getCost();
		
		return result;
	}

	//gene id spec1       \t      geneid spec2      \t        chrpos spec1     \t      chrpos spec2
	public String toStringAndre() {
		Object[] o= mapTableI.keySet().toArray();
		SpliceNode[] nodes= new SpliceNode[o.length];
		for (int i = 0; i < nodes.length; i++) 
			nodes[i]= (SpliceNode) o[i];
		Arrays.sort(nodes, new SpliceNode.PositionComparator());
		
		String result1= "";
		String result2= "";
		for (int i = 0; i < nodes.length; i++) {
			result1+= Math.abs(nodes[i].getSite().getPos())+ " "; 
			result2+= Math.abs(((SpliceNode) mapTableI.get(nodes[i])).getSite().getPos())+ " ";
		}
		if (result1.length()> 0)
			result1= result1.substring(0, result1.length()-1);
		if (result2.length()> 0)
			result2= result2.substring(0, result2.length()-1);
		
		return (result1+ "\t"+ result2);
	}


	public int getMaxRel() {
		return maxRel;
	}

	public boolean isFreeShift() {
		return freeShift;
	}

	public void setFreeShift(boolean freeShift) {
		this.freeShift = freeShift;
	}

	public static String getCostID() {
		return costID;
	}

	public static void setCostID(String newCostID) {
		costID = newCostID;
	}
}
