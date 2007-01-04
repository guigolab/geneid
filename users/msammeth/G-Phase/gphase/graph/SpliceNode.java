package gphase.graph;

import java.util.Comparator;
import java.util.HashMap;

import gphase.model.AbstractSite;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class SpliceNode {
	
	public static class PositionTypeComparator extends PositionComparator {
		public int compare(Object arg0, Object arg1) {
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;
			
			arg0= ((SpliceNode) arg0).getSite();
			arg1= ((SpliceNode) arg1).getSite();
			
				// same position, distinguish type
			if ((!(arg0 instanceof SpliceSite))&& arg1 instanceof SpliceSite)
				return -1;
			if ((!(arg1 instanceof SpliceSite))&& arg0 instanceof SpliceSite)
				return 1;
			if ((!(arg1 instanceof SpliceSite))&& (!(arg0 instanceof SpliceSite)))
				return 0;
			
			SpliceSite s0= (SpliceSite) arg0;
			SpliceSite s1= (SpliceSite) arg1;
			if (s0.isAcceptor()&& s1.isDonor())
				return -1;
			if (s1.isAcceptor()&& s0.isDonor())
				return 1;
			return 0;
		}
	}
	
	public static class PositionComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			SpliceNode node0= (SpliceNode) arg0;
			SpliceNode node1= (SpliceNode) arg1;
			if (node0.getSite().getPos()< node1.getSite().getPos())
				return -1;
			if (node0.getSite().getPos()> node1.getSite().getPos())
				return 1;
			return 0;
		}
	}
	
	AbstractSite site= null;
	SpliceEdge[] inEdges= null;
	SpliceEdge[] outEdges= null;
	HashMap fromList= null;	// maps ss pos to #pathes
	
	public SpliceNode (AbstractSite newSS) {
		this.site= newSS;
	}
	
	public String toString() {
		return Integer.toString(getSite().getPos());
	}
	
	public int getOutDegree(boolean max) {
		return getDegree(outEdges, max);
	}
	
	public int getOutDegree() {
		return getOutDegree(false);
	}
	
	public int getInDegree() {
		return getInDegree(false);
	}
	
	public int getInDegree(boolean max) {
		return getDegree(inEdges, max);
	}
	int getDegree(SpliceEdge[] edges, boolean max) {
		if (edges== null)
			return 0;

		int x= 0;
		if (max)
			for (int i = 0; i < edges.length; i++) 
				x+= edges[i].getTranscripts().length;
		else
			x= edges.length;
		
		return x;
	}
	

	public boolean isDonor() {
		try {
			return ((SpliceSite) site).isDonor();
		} catch (Exception e) {
			return false;
		}
	}
	
	public boolean isAcceptor() {
		try {
			return !((SpliceSite) site).isDonor();
		} catch (Exception e) {
			return false;
		}
	}
	
	public boolean isTSS() {
		return (getInDegree()== 0);
	}
	
	public boolean isTES() {
		return (getOutDegree()== 0);
	}
	
	public void addOutEdge(SpliceEdge newOutEdge) {
		outEdges= addEdge(outEdges, newOutEdge);
	}
	
	public void addInEdge(SpliceEdge newInEdge) {
		inEdges= addEdge(inEdges, newInEdge);
	}
	
	public void removeInEdge(SpliceEdge remEdge) {
		inEdges= (SpliceEdge[]) Arrays.remove(inEdges, remEdge);
	}
	
	public void removeOutEdge(SpliceEdge remEdge) {
		outEdges= (SpliceEdge[]) Arrays.remove(outEdges, remEdge);
	}
	
	SpliceEdge[] addEdge(SpliceEdge[] edges, SpliceEdge newEdge) {
		
		if (edges== null) 
			return new SpliceEdge[] {newEdge};
		
			// check for non-redundancy
//		for (int i = 0; i < edges.length; i++) 
//			if (edges[i].equals(newEdge)) {
//				Transcript[] trans= newEdge.getTranscripts();
//				for (int j = 0; j < trans.length; j++) 
//					edges[i].addTranscript(trans[j]);
//				return edges;
//			}
		
			// add
//		SpliceEdge[] newEdges= new SpliceEdge[edges.length+ 1];
//		for (int i = 0; i < edges.length; i++) 
//			newEdges[i]= edges[i];
//		newEdges[newEdges.length- 1]= newEdge;
		edges= (SpliceEdge[]) Arrays.extendField(edges, newEdge);
		return edges;
	}

	public AbstractSite getSite() {
		return site;
	}

	public HashMap getFromList() {
		if (fromList == null) {
			fromList = new HashMap();	// mainly for initing tss/tes JIT..
		}
		return fromList;
	}

	public SpliceEdge[] getInEdges() {
		return inEdges;
	}

	public SpliceEdge[] getOutEdges() {
		return outEdges;
	}

	public void setFromList(HashMap fromList) {
		this.fromList = fromList;
	}
	
	
}
