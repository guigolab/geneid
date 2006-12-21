package gphase.graph;

import gphase.model.ASMultiVariation;
import gphase.model.AbstractSite;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.model.ASVariation.SpliceChainComparator;
import gphase.model.AbstractSite.PositionComparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Stack;
import java.util.Vector;

import javax.swing.tree.TreePath;

import com.sun.org.apache.bcel.internal.generic.IndexedInstruction;

public class SpliceGraph {
	public static class NodeOrderComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			AbstractSite a0= ((SpliceNode) arg0).getSite(); 
			AbstractSite a1= ((SpliceNode) arg1).getSite(); 
			Comparator compi= new AbstractSite.PositionComparator();
			int x= compi.compare(a0, a1);
			if (x!= 0)
				return x;
			if (!(a0 instanceof SpliceSite) && a1 instanceof SpliceSite)
				return -1;
			if (a0 instanceof SpliceSite && !(a1 instanceof SpliceSite))
				return 1;
			return 0;
		}
	}
	
	Transcript[] transcripts= null;
	HashMap nodeList= null;	// SpliceNodes / Vector{SplicePathes..}
	SpliceBubble[] bubbles= null;
	
	public SpliceGraph(Transcript[] newTranscripts) {
		this.transcripts= newTranscripts;
		nodeList= new HashMap();	// optimize by providing default size ?!
	}
	
	public void init() {
			// add splice sites
		for (int i = 0; i < transcripts.length; i++) {
			SpliceSite[] sss= transcripts[i].getSpliceChain();
			if (sss== null|| sss.length< 1)
				continue;
			SpliceNode node= new SpliceNode(sss[0]);
			addNode(node);
			for (int j = 1; j < sss.length; j++) {
				node= new SpliceNode(sss[j]);
				node= addNode(node);				
				createEdge(getNode(sss[j-1]), node, new Transcript[] {transcripts[i]});
			}
		}
		
			// add border anchors
		Vector tssV= new Vector();
		for (int i = 0; i < transcripts.length; i++) {
			for (int j = i+1; j < transcripts.length; j++) {
					// find trim points
				int[] t= trim(new Transcript[] {transcripts[i], transcripts[j]});
				int end5= t[0];
				int end3= t[1];

				if (end3== 0|| end5== 0|| end3== end5) {
					if (end3== 0^ end5== 0)
						System.err.println("Assertion failed: only one valid trimming site, cannot be!");
					continue;
				}
				
					// create sites here, for having both transcripts
				AbstractSite tss= new AbstractSite(end5);
				tss.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
				AbstractSite tes= new AbstractSite(end3);
				tes.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
				createEdge(transcripts[i], tss, tes);
				createEdge(transcripts[j], tss, tes);
				
				SpliceNode tmpTss= (SpliceNode) nodeList.get(tss);
				int k;
				for (k = 0; k < tssV.size(); k++) 
					if (tmpTss== tssV.elementAt(k))
						break;
				if (k== tssV.size())
					tssV.add(tmpTss);	// cannot check for forbidden connections at that point
									// due to transitivity (3 or more transcripts)
				
				tmpTss= (SpliceNode) nodeList.get(tes);
				for (k = 0; k < tssV.size(); k++) 
					if (tmpTss== tssV.elementAt(k))
						break;
				if (k== tssV.size())
					tssV.add(tmpTss);

//				if (end5!= 0) {	// overlap found
//					SpliceSite chk1= transcripts[i].getSpliceSite(end5);	// forbid anchor -> acceptor
//					SpliceSite chk2= transcripts[j].getSpliceSite(end5);	// no valid event
//					if ((chk1== null|| chk1.isDonor())&& (chk2== null|| chk2.isDonor())) {
//						AbstractSite tss= new AbstractSite(end5);
//						tss.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
//						SpliceNode tssN= new SpliceNode(tss);
//						tssN= addNode(tssN);
//							// connect to schains
//						createEdge(transcripts[i], tss, compi, tssN, true);
//						createEdge(transcripts[j], tss, compi, tssN, true);
//					}
//				}
//
//				int end3= t[1];
//				if (end3!= 0) {
//					SpliceSite chk1= transcripts[i].getSpliceSite(end3);	// forbid donor -> anchor
//					SpliceSite chk2= transcripts[j].getSpliceSite(end3);
//					if ((chk1== null|| chk1.isAcceptor())&& (chk2== null|| chk2.isAcceptor())) {
//						AbstractSite tes= new AbstractSite(end3);
//						tes.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
//						SpliceNode tesN= new SpliceNode(tes);
//						tesN= addNode(tesN);
//						createEdge(transcripts[i], tes, compi, tesN, false);
//						createEdge(transcripts[j], tes, compi, tesN, false);
//					}
//				}
			}
		}
			// forbid AS -> acceptor, donor -> AS
		contractInvalidEdges(tssV);
			// forbid 2 AS adjacent in same exon 
		eliminateRedundantTSS(tssV);
	}
	
	public void init_stable() {
		// add splice sites
	for (int i = 0; i < transcripts.length; i++) {
		SpliceSite[] sss= transcripts[i].getSpliceChain();
		if (sss== null|| sss.length< 1)
			continue;
		SpliceNode node= new SpliceNode(sss[0]);
		addNode(node);
		for (int j = 1; j < sss.length; j++) {
			node= new SpliceNode(sss[j]);
			node= addNode(node);				
			createEdge(getNode(sss[j-1]), node, new Transcript[] {transcripts[i]});
		}
	}
	
		// add border anchors
	Vector tssV= new Vector();
	for (int i = 0; i < transcripts.length; i++) {
		for (int j = i+1; j < transcripts.length; j++) {
				// find trim points
			int[] t= trim(new Transcript[] {transcripts[i], transcripts[j]});
			int end5= t[0];
			int end3= t[1];

			if (end3== 0|| end5== 0|| end3== end5) {
				if (end3== 0^ end5== 0)
					System.err.println("Assertion failed: only one valid trimming site, cannot be!");
				continue;
			}
			
				// create sites here, for having both transcripts
			AbstractSite tss= new AbstractSite(end5);
			tss.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
			AbstractSite tes= new AbstractSite(end3);
			tes.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
			createEdge(transcripts[i], tss, tes);
			createEdge(transcripts[j], tss, tes);
			
			SpliceNode tmpTss= (SpliceNode) nodeList.get(tss);
			int k;
			for (k = 0; k < tssV.size(); k++) 
				if (tmpTss== tssV.elementAt(k))
					break;
			if (k== tssV.size())
				tssV.add(tmpTss);	// cannot check for forbidden connections at that point
								// due to transitivity (3 or more transcripts)
			
			tmpTss= (SpliceNode) nodeList.get(tes);
			for (k = 0; k < tssV.size(); k++) 
				if (tmpTss== tssV.elementAt(k))
					break;
			if (k== tssV.size())
				tssV.add(tmpTss);

//			if (end5!= 0) {	// overlap found
//				SpliceSite chk1= transcripts[i].getSpliceSite(end5);	// forbid anchor -> acceptor
//				SpliceSite chk2= transcripts[j].getSpliceSite(end5);	// no valid event
//				if ((chk1== null|| chk1.isDonor())&& (chk2== null|| chk2.isDonor())) {
//					AbstractSite tss= new AbstractSite(end5);
//					tss.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
//					SpliceNode tssN= new SpliceNode(tss);
//					tssN= addNode(tssN);
//						// connect to schains
//					createEdge(transcripts[i], tss, compi, tssN, true);
//					createEdge(transcripts[j], tss, compi, tssN, true);
//				}
//			}
//
//			int end3= t[1];
//			if (end3!= 0) {
//				SpliceSite chk1= transcripts[i].getSpliceSite(end3);	// forbid donor -> anchor
//				SpliceSite chk2= transcripts[j].getSpliceSite(end3);
//				if ((chk1== null|| chk1.isAcceptor())&& (chk2== null|| chk2.isAcceptor())) {
//					AbstractSite tes= new AbstractSite(end3);
//					tes.addTranscripts(new Transcript[] {transcripts[i], transcripts[j]});
//					SpliceNode tesN= new SpliceNode(tes);
//					tesN= addNode(tesN);
//					createEdge(transcripts[i], tes, compi, tesN, false);
//					createEdge(transcripts[j], tes, compi, tesN, false);
//				}
//			}
		}
	}
		// forbid AS -> acceptor, donor -> AS
	contractInvalidEdges(tssV);
}	/**
	 * Removes redundant TSS/TES, which are further away from the first real SS
	 * and (by this) representing a smaller transcript partition set with
	 * no additional variation added.
	 * @param tssV
	 */
	void eliminateRedundantTSS(Vector tssV) {
		Vector chkedNodesV= new Vector(tssV.size());
		Vector remNodesV= new Vector();
		for (int i = 0; i < tssV.size(); i++) {
				// chk tss/tes
			SpliceNode tmpNode= (SpliceNode) tssV.elementAt(i);
			boolean tss= false;
			if (tmpNode.getInEdges()== null|| tmpNode.getInEdges().length< 1)
				tss= true;
			else if (tmpNode.getOutEdges()== null|| tmpNode.getOutEdges().length< 1)
				tss= false;
			else
				System.err.println("Neither tss nor tes!");
			
				// iterate all target not of all out/in edges
			SpliceEdge[] tmpEdges= tss?tmpNode.getOutEdges():tmpNode.getInEdges();
			for (int j = 0; j < tmpEdges.length; j++) {
				SpliceNode targetNode= tss?tmpEdges[j].getHead():tmpEdges[j].getTail();
				int k;
				for (k = 0; k < chkedNodesV.size(); k++) 
					if (chkedNodesV.elementAt(k)== targetNode)
						break;
				if (k< chkedNodesV.size())		// already purified..
					continue;
				
					// check for other abstract sites connected to this node
				chkedNodesV.add(targetNode);
				SpliceEdge[] targetEdges= tss?targetNode.getInEdges():targetNode.getOutEdges();
				SpliceNode refNode= tmpNode;
				for (k = 0; k < targetEdges.length; k++) {
					SpliceNode altSrcNode= tss?targetEdges[k].getTail():targetEdges[k].getHead();
					if (altSrcNode== refNode)
						continue;
					if (!(altSrcNode.getSite() instanceof SpliceSite)) {
						if (tss) {
							if (refNode.getSite().getPos()< altSrcNode.getSite().getPos()) {
								if (SpliceBubble.contained(refNode.getSite().getTranscripts(), altSrcNode.getSite().getTranscripts())) {
									remNodesV.add(refNode);
									refNode= altSrcNode;
								}
							} else if (SpliceBubble.contained(altSrcNode.getSite().getTranscripts(), refNode.getSite().getTranscripts())) {
								remNodesV.add(altSrcNode);
							}
						} else {	// tes
							if (refNode.getSite().getPos()> altSrcNode.getSite().getPos()) {
								if (SpliceBubble.contained(refNode.getSite().getTranscripts(), altSrcNode.getSite().getTranscripts())) {
									remNodesV.add(refNode);
									refNode= altSrcNode;
								}
							} else if (SpliceBubble.contained(altSrcNode.getSite().getTranscripts(), refNode.getSite().getTranscripts())) {
								remNodesV.add(altSrcNode);
							}
						}
					}
				}
			}
		}
		
			// delete
		for (int i = 0; i < remNodesV.size(); i++) {
			tssV.remove(remNodesV.elementAt(i));	// necessary?
			removeNode((SpliceNode) remNodesV.elementAt(i));
		}
	}
	
	public void removeNode(SpliceNode delNode) {
		SpliceEdge[] edges= delNode.getInEdges();
		for (int j = 0; edges!= null&& j < edges.length; j++) 
			edges[j].getTail().removeOutEdge(edges[j]);
		edges= delNode.getOutEdges();
		for (int j = 0; edges!= null&& j < edges.length; j++) 
			edges[j].getHead().removeInEdge(edges[j]);
		
		Object o= nodeList.remove(delNode.getSite());
		System.currentTimeMillis();

	}
	
	void contractInvalidEdges(Vector tssV) {
		
		for (int i = 0; i < tssV.size(); i++) {
			SpliceNode tmpNode= (SpliceNode) tssV.elementAt(i);
			boolean tss= false;
			if (tmpNode.getInEdges()== null|| tmpNode.getInEdges().length< 1)
				tss= true;
			else if (tmpNode.getOutEdges()== null|| tmpNode.getOutEdges().length< 1)
				tss= false;
			else
				System.err.println("Neither tss nor tes!");
			
			SpliceEdge[] edges= null;
			if (tss)
				edges= tmpNode.getOutEdges();
			else 
				edges= tmpNode.getInEdges();
			
			for (int j = 0; j < edges.length; j++) {	// check all links from tss/ to tes: edges[j]
				AbstractSite site= null;
				if (tss)
					site= edges[j].getHead().getSite();
				else
					site= edges[j].getTail().getSite();

					// validity check: allowed links
				if (tss&& ((!(site instanceof SpliceSite))|| ((SpliceSite) site).isDonor())||
						!tss&& ((!(site instanceof SpliceSite))|| ((SpliceSite) site).isAcceptor()))
					continue;
				
					// not allowed link: edges[j]
				SpliceNode accN= null;
				if (tss) 
					accN= edges[j].getHead();	// other node of delete link
				else
					accN= edges[j].getTail();
				for (int k = 0; k < edges[j].getTranscripts().length; k++) {
					int m;
					SpliceEdge[] edges2= null;
					if (tss)
						edges2= accN.getOutEdges();	// get other links for triangle re-linking
					else
						edges2= accN.getInEdges();
					if (edges2== null)
						continue;
					for (m = 0; m < edges2.length; m++) {
						int n;
						for (n = 0; n < edges2[m].getTranscripts().length; n++) {
							if (edges2[m].getTranscripts()[n]== edges[j].getTranscripts()[k])
								break;
						}
						if (n< edges2[m].getTranscripts().length)	// transcript found
							break;
					}
					if (m< edges2.length) {
						SpliceEdge[] edges3= null;	// compensatory outgoing links from triangle node to be shortcut
													// ... do NOT remove links, if there is another tss / tes with
													// different transcript set (!)
						if (tss)
							edges3= accN.getOutEdges();
						else
							edges3= accN.getInEdges();
						int x;
						for (x = 0; x < edges3.length; x++) {
							int y;
							for (y = 0; y < edges3[x].getTranscripts().length; y++) 
								if (edges3[x].getTranscripts()[y]== edges[j].getTranscripts()[k])
									break;
							if (y< edges3[x].getTranscripts().length)
								break;
						}
						if (x>= edges3.length)	// do only delete, if transcript doesnt support to continue schain
							edges2[m].removeTranscript(edges[j].getTranscripts()[k]);
						
						if (edges2[m].getTranscripts()== null|| edges2[m].getTranscripts().length< 1) {
							edges2[m].getTail().removeOutEdge(edges2[m]);
							if (edges2[m].getTail().getOutEdges()== null||
									edges2[m].getTail().getOutEdges().length < 1)
								nodeList.remove(edges2[m].getTail());
							edges2[m].getHead().removeInEdge(edges2[m]);
							if (edges2[m].getHead().getInEdges()== null||
									edges2[m].getHead().getInEdges().length < 1)
								nodeList.remove(edges2[m].getHead());
						}
						if (tss)
							createEdge(tmpNode, edges2[m].getHead(), new Transcript[] {edges[j].getTranscripts()[k]});
						else
							createEdge(edges2[m].getTail(), tmpNode, new Transcript[] {edges[j].getTranscripts()[k]});
					} else {
						System.err.println("Transcript for contraction not found!");
					}
				}
				edges[j].getTail().removeOutEdge(edges[j]);
				if (edges[j].getTail().getOutEdges()== null||
						edges[j].getTail().getOutEdges().length < 1)
					nodeList.remove(edges[j].getTail());
				edges[j].getHead().removeInEdge(edges[j]);
				if (edges[j].getHead().getInEdges()== null||
						edges[j].getHead().getInEdges().length < 1)
					nodeList.remove(edges[j].getHead());
				
			}
		}
	}
	void createEdge(Transcript trans, AbstractSite trimEnd5, AbstractSite trimEnd3) {
		
			// make nodes, check before whether sides are valid
		SpliceNode tssN= new SpliceNode(trimEnd5);
		SpliceNode tssN2= addNode(tssN);
		if ((tssN!= tssN2)&& (tssN2.getInDegree(false)== 0))
			tssN= tssN2;					// distinguish tss / tes nodes, do not merge them
		SpliceNode tesN= new SpliceNode(trimEnd3);
		SpliceNode tesN2= addNode(tesN);	// distinguish tss / tes nodes, do not merge them
		if ((tesN!= tesN2)&& (tesN2.getOutDegree(false)== 0))
			tesN= tesN2;

			// try to get ss to connect to
		SpliceSite[] sc= trans.getSpliceChain();
		PositionComparator compi= new AbstractSite.PositionComparator();
		int p5= Arrays.binarySearch(sc, trimEnd5, compi);
		if (p5< 0) 
			p5= -(p5+1);
//		else if (sc[p5].isAcceptor())	// do not connect tss -> acceptor
//			++p5;						// try next ss
		int p3= Arrays.binarySearch(sc, trimEnd3, compi);
		if (p3< 0) 
			p3= -(p3+1)- 1;// one before for end connections
//		else if (sc[p3].isDonor())	// do not connect donor -> tes
//			--p3;				// try prev ss
		
		
		if (p5>= 0&& p5< sc.length&& p5<= p3) {
			createEdge(tssN, getNode(sc[p5]), new Transcript[] {trans});
		} else if (p3< 0|| p3>= sc.length|| p3< p5) {	// no valid ss for both trimmed ends, connect them (for events with empty sc)
			createEdge(tssN, tesN, new Transcript[] {trans});
		}
		if (p3>= 0&& p3< sc.length&& p3>= p5) {
			createEdge(getNode(sc[p3]), tesN, new Transcript[] {trans});
		}
	}

	boolean createEdge_old(Transcript trans, AbstractSite trimEnd, Comparator compi, SpliceNode node, boolean tss) {
		SpliceSite[] sc= trans.getSpliceChain();
		int p= Arrays.binarySearch(sc, trimEnd, compi);
		if (p< 0) {
			p= -(p+1);
			if (!tss)
				--p;	// one before for end connections
		}
		try {
			SpliceNode tmpNode= (SpliceNode) getNode(sc[p]);
			if (tss) {
				if (((SpliceSite) tmpNode.getSite()).isAcceptor())	// do not connect tss->acceptor
					tmpNode= (SpliceNode) getNode(sc[p+1]);	// try next ss (if any..)
				createEdge(node, tmpNode, new Transcript[] {trans});
			} else {
				if (((SpliceSite) tmpNode.getSite()).isDonor())		// do not connect donor -> tes
					tmpNode= (SpliceNode) getNode(sc[p+1]);	// try next ss (if any..)
				createEdge(tmpNode, node, new Transcript[] {trans});
			}
		} catch (ArrayIndexOutOfBoundsException e) {
			return false; // :)				
		}
		return true;
	}
	
	SpliceEdge createEdge(SpliceNode tail, SpliceNode head, Transcript[] trans) {
		SpliceNode node= (SpliceNode) nodeList.get(tail.getSite());
		if (node!= null) {
			SpliceEdge[] edges= node.getOutEdges();
			if (edges!= null) 
				for (int j = 0; j < edges.length; j++) {
					if (edges[j].getHead()== head) {
						edges[j].addTranscripts(trans);
						return edges[j];
					}
				}
		}
		
		return new SpliceEdge(tail, head, trans);
	}
	
	public SpliceNode addNode(SpliceNode newNode) {
		AbstractSite key= newNode.getSite();
		SpliceNode node= (SpliceNode) nodeList.get(key);
		if (node!= null) {	// check also for same content (tss / ss)
			node.getSite().addTranscripts(newNode.getSite().getTranscripts());
			return node;
		}
		
			// else
		nodeList.put(key, newNode);
		return newNode;
	}
	
	public SpliceNode getNode(SpliceSite ss) {
		return (SpliceNode) nodeList.get(ss);
	}
	
	/**
	 * @return node list sorted according to genomic positions
	 */
	public SpliceNode[] getNodeList() {
		
		Collection c= nodeList.values();
		SpliceNode[] nodes= new SpliceNode[c.size()];
		nodes= (SpliceNode[]) c.toArray(nodes);
		
		Arrays.sort(nodes, new NodeOrderComparator());
		return nodes;
	}
	
	public SpliceNode[] getRoots() {
		SpliceNode[] nodes= getNodeList();
		Vector v= new Vector();
		for (int i = 0; i < nodes.length; i++) 
			if (nodes[i].getInDegree()== 0)
				v.add(nodes[i]);
		
		return (SpliceNode[]) gphase.tools.Arrays.toField(v);
	}
	
	public SpliceNode[] getLeafs() {
		SpliceNode[] nodes= getNodeList();
		Vector v= new Vector();
		for (int i = 0; i < nodes.length; i++) 
			if (nodes[i].getOutDegree()== 0)
				v.add(nodes[i]);
		
		return (SpliceNode[]) gphase.tools.Arrays.toField(v);
	}
	
	
	public SplicePath[] getPathes(SpliceNode[] src, SpliceNode[] tgt) {
		Vector v= new Vector();
		for (int i = 0; i < src.length; i++) {
			v= getPathesRek(src[i], tgt, new SplicePath(), v);
		}
		return (SplicePath[]) gphase.tools.Arrays.toField(v);
	}
	
	Vector getPathesRek(SpliceNode nd, SpliceNode[] tgt, SplicePath ndPath, Vector result) {
		
			// abort
		for (int i = 0; i < tgt.length; i++) 
			if (nd== tgt[i]) {
				result.add(ndPath);
				return result;
			}
		
			// rekursion
		SpliceEdge[] out= nd.getOutEdges();
		for (int i = 0; out!= null&& i < out.length; i++) 
			getPathesRek(out[i].getHead(), tgt, ndPath.exendPath(out[i]), result);
		
		
		return result;
	}
	
	public SpliceEdge[] getEdgeList() {
		// assuming that all edges have a src node
		SpliceNode[] n= getNodeList();
		Vector eV= new Vector(n.length);
		for (int i = 0; i < n.length; i++) 
			for (int j = 0; j < n[i].getOutDegree(); j++) 
				eV.add(n[i].getOutEdges()[j]);
			
		SpliceEdge[] e= (SpliceEdge[]) gphase.tools.Arrays.toField(eV);
		return e;
	}
	
	void checkRedundancy(SpliceBubble blob) {
		
		Vector rmV= new Vector();
		int x= 0;
		boolean addBlob= true;
		for (x = 0; x < bubbles.length; x++) 
			if (bubbles[x].getSource().getSite().getPos()>= blob.getSource().getSite().getPos())
				break;
		
		Transcript[][] nuParts= blob.getTranscriptPartitions();
		for (int i = x; i < bubbles.length; i++) {	// all bubbles after p have same start and more little ends
			
			if (bubbles[i].getSink().getSite().getPos()> blob.getSink().getSite().getPos())
				break;
			
				// in contained, check transcript set
			Transcript[][] tmpPart= bubbles[i].getTranscriptPartitions();
			if (SpliceBubble.contained(nuParts, tmpPart))
				addBlob= false; // blob has wider borders and equal (or kleiner) transcript set, discard blob 
			else if (tmpPart.length< nuParts.length&&	// NEW: to be REALLY contained !!
					SpliceBubble.contained(tmpPart, nuParts)) {	// blob is larger, but contains other transcript set
				rmV.add(bubbles[i]);
				blob.addContainedBubble(bubbles[i]);
			}			
		}
		
			// remove
		for (int i = 0; i < rmV.size(); i++) 
			bubbles= (SpliceBubble[]) gphase.tools.Arrays.remove(bubbles, rmV.elementAt(i));
		
			// add blob
		if (!addBlob)
			return;
		if (bubbles== null)
			bubbles= new SpliceBubble[0];
		
		int p= Arrays.binarySearch(bubbles, blob, new SpliceBubble.PositionComparator());
		if (p>= 0) {
			if (bubbles[p].getSource().getSite() instanceof SpliceSite== blob.getSource().getSite() instanceof SpliceSite&&
					bubbles[p].getSink().getSite() instanceof SpliceSite== blob.getSink().getSite() instanceof SpliceSite)
				System.err.println("assertion failed !");
		}
		bubbles= (SpliceBubble[]) gphase.tools.Arrays.insert(bubbles, blob, p);
	
	}

	void checkRedundancy_old(SpliceBubble blob) {
		
		Vector rmV= new Vector();
		int x= 0;
		boolean addBlob= true;
		for (x = 0; x < bubbles.length; x++) 
			if (bubbles[x].getSource().getSite().getPos()>= blob.getSource().getSite().getPos())
				break;
		
		Transcript[][] nuParts= blob.getTranscriptPartitions();
		System.currentTimeMillis();
		for (int i = x; i < bubbles.length; i++) {	// all bubbles after p have same start and more little ends
			
			if (bubbles[i].getSink().getSite().getPos()> blob.getSink().getSite().getPos())
				break;
			
				// in contained, check transcript set
			Transcript[][] tmpPart= bubbles[i].getTranscriptPartitions();
			if (SpliceBubble.contained(nuParts, tmpPart))
				addBlob= false; // blob has wider borders and equal (or kleiner) transcript set, discard blob 
			else if (tmpPart.length< nuParts.length&&	// NEW: to be REALLY contained !!
					SpliceBubble.contained(tmpPart, nuParts)) {	// blob is larger, but contains other transcript set
				rmV.add(bubbles[i]);
				blob.addContainedBubble(bubbles[i]);
			}			
		}
		
			// remove
		for (int i = 0; i < rmV.size(); i++) 
			bubbles= (SpliceBubble[]) gphase.tools.Arrays.remove(bubbles, rmV.elementAt(i));
		
			// add blob
		if (!addBlob)
			return;
		if (bubbles== null)
			bubbles= new SpliceBubble[0];
		
		int p= Arrays.binarySearch(bubbles, blob, new SpliceBubble.PositionComparator());
		if (p>= 0) {
			if (bubbles[p].getSource().getSite() instanceof SpliceSite== blob.getSource().getSite() instanceof SpliceSite&&
					bubbles[p].getSink().getSite() instanceof SpliceSite== blob.getSink().getSite() instanceof SpliceSite)
				System.err.println("assertion failed !");
		}
		bubbles= (SpliceBubble[]) gphase.tools.Arrays.insert(bubbles, blob, p);

	}

	void checkRedundancy_veryold(SpliceBubble blob) {
		
		int p= Arrays.binarySearch(bubbles, blob, new SpliceBubble.PositionComparator());
		if (p< 0)
			p= -(p+1);
		else
			++p;
		Vector rmV= new Vector();
		for (int i = p; i < bubbles.length; i++) {	// all bubbles after p have same start and more little ends
													// or bigger start (and same or more little ends)
			if (blob.contains(bubbles[i]))
				rmV.add(bubbles[i]);
		}
		
			// remove
		for (int i = 0; i < rmV.size(); i++) 
			bubbles= (SpliceBubble[]) gphase.tools.Arrays.remove(bubbles, rmV.elementAt(i));		
	}
	
	/**
	 * 
	 * @param k (-1) for max extended variations
	 * @return
	 */
	public ASMultiVariation[] getMultiVariations_stable(int k) {
		
		SpliceBubble[] blobs= getBubbles();
		if (blobs== null)
			return null;
		Vector multiVars= new Vector(blobs.length);
		for (int i = 0; i < blobs.length; i++) {
			SpliceNode[][] pathes= blobs[i].getPathes();
			SpliceSite[][] sc= new SpliceSite[pathes.length][];
			HashMap transHash= new HashMap();
			for (int j = 0; j < sc.length; j++) {
				if (pathes[j]== null)
					sc[j]= new SpliceSite[0];
				else {
					sc[j]= new SpliceSite[pathes[j].length];
					for (int x = 0; x < pathes[j].length; x++) 
						sc[j][x]= (SpliceSite) pathes[j][x].getSite();
				}
				
				transHash.put(sc[j], blobs[i].getTransHash().get(pathes[j]));
			}
			
				// generate mvariations
			if (k< 0) {
				multiVars.add(new ASMultiVariation(sc, transHash));
			} else {	// generate all k-mers
				// this I will do for you, Sylvain
			}
			
		}
		return (ASMultiVariation[]) gphase.tools.Arrays.toField(multiVars);
	}
	
	/**
	 * 
	 * @param k (-1) for max extended variations
	 * @return
	 */
	public ASMultiVariation[] getMultiVariations(int k) {
		
		SpliceBubble[] blobs= getBubbles();
		if (blobs== null)
			return null;
		Vector multiVars= new Vector(blobs.length);
		HashMap usedVars= new HashMap();
		Comparator compi= new ASMultiVariation.PositionComparator();
		for (int i = 0; i < blobs.length; i++) {
			extractVars(k, blobs[i], usedVars, multiVars);
		}
			// filter redundancy, TODO very inefficient
		Vector remVector= new Vector(multiVars.size());
		for (int j = 0; j < multiVars.size(); j++) 
			for (int x = j+1; x < multiVars.size(); x++) 
				if (compi.compare(multiVars.elementAt(j), multiVars.elementAt(x))== 0)
					remVector.add(multiVars.elementAt(x));
		multiVars.removeAll(remVector);
		
		return (ASMultiVariation[]) gphase.tools.Arrays.toField(multiVars);
	}
	
	void extractVars(int k, SpliceBubble blob, HashMap procH, Vector multiVars) {
		
			// check whether this bubble has already been processed
		if (procH.get(blob)!= null)
			return;
		procH.put(blob, blob);
		
		SpliceNode[][] pathes= blob.getPathes();
		SpliceSite[][] sc= new SpliceSite[pathes.length][];
		HashMap transHash= new HashMap();
		for (int j = 0; j < sc.length; j++) {
			if (pathes[j]== null)
				sc[j]= new SpliceSite[0];
			else {
				sc[j]= new SpliceSite[pathes[j].length];
				for (int x = 0; x < pathes[j].length; x++) 
					sc[j][x]= (SpliceSite) pathes[j][x].getSite();
			}
			
			transHash.put(sc[j], blob.getTransHash().get(pathes[j]));
		}
		
			// generate mvariations
		if (k< 2) {	// output max variations
			multiVars.add(new ASMultiVariation(sc, transHash));	// TODO children !!!?!!
		} else { // generate all k-mers
				// check for contained subbubbles
			SpliceBubble[] sbubl= blob.getChildren();
			Vector transNot= null;
			if (sbubl!= null&& sbubl.length> 0) {
				transNot= new Vector(sbubl.length);
				for (int j = 0; j < sbubl.length; j++) {
					transNot.add(sbubl[j].getTranscriptPartitions());
					extractVars(k, sbubl[j], procH, multiVars);
				}
			}
			
			reKmer(k, sc, transHash, 0, new int[0], multiVars, transNot);
		}
		

	}
	private void reKmer(int maxDepth, SpliceSite[][] sc, HashMap transHash, int minVal, int[] vals, Vector result, Vector transNot) {
		
		if (vals.length== maxDepth) {
			SpliceSite[][] ss= new SpliceSite[vals.length][];
			HashMap hh= new HashMap(3);
			for (int i = 0; i < ss.length; i++) { 
				ss[i]= sc[vals[i]];
				hh.put(ss[i], transHash.get(ss[i]));
			}
			
				// check whether combination is fully contained in any sub-bubble
			Transcript[][] thisT= (Transcript[][]) gphase.tools.Arrays.toField(hh.values().toArray());
			for (int i = 0; transNot!= null&& i < transNot.size(); i++) {
				Transcript[][] subublT= (Transcript[][]) transNot.elementAt(i);
				if (SpliceBubble.contained(thisT, subublT))
					return;	// skip
			}
			
				// else add
		    ASMultiVariation x= new ASMultiVariation(ss, hh);
			result.add(x);
			return;
		}
		
		for (int i = minVal; i <= (sc.length- maxDepth+ 1); i++) {
			int[] nuVals= new int[vals.length+ 1];
			for (int j = 0; j < vals.length; j++) 
				nuVals[j]= vals[j];
			nuVals[vals.length]= i;
			reKmer(maxDepth, sc, transHash, i+1, nuVals, result, transNot);
		}
	}
	
	public SpliceBubble[] getBubbles() {
		
		if (bubbles== null) {
			//System.err.println(transcripts[0].getTranscriptID());
			SpliceBubble[] allBubs= findBubbles();
			bubbles= cascadeBubbles(allBubs);
		}
		return bubbles;
	}
	
	public SpliceBubble[] getBubbles_stable() {
		
		if (bubbles== null) {
			SpliceNode[] nodes= getNodeList();
			for (int i = 0; i < nodes.length; i++) {
				
				SpliceNode tmpNode= nodes[i];
				HashMap map= tmpNode.getFromList();
				
					// opening bubbles: add node to fromList
				SpliceEdge[] outs= tmpNode.getOutEdges();
				for (int j = 0; outs!= null&& j < outs.length; j++) {	// outs.length> 1?
					HashMap newMap= new HashMap(map.size());
					SpliceNode[] fromNodes= (SpliceNode[]) gphase.tools.Arrays.toField(map.keySet().toArray());
					for (int k = 0; fromNodes!= null&& k < fromNodes.length; k++) {	// remove orphan fromNodes
						Vector pathes= (Vector) map.get(fromNodes[k]);
						Vector newPathes= new Vector();
						for (int h = 0; h < pathes.size(); h++) {
							SplicePath p= ((SplicePath) pathes.elementAt(h)).createPath(outs[j]);
							if (p!= null)
								newPathes.add(p);
						}
						if (newPathes.size()> 0)	// fromNode no longer valid, remove
							newMap.put(fromNodes[k], newPathes);
					}
					SplicePath sp= new SplicePath(outs[j]);
					Vector tmpV= new Vector();
					tmpV.add(sp);
					newMap.put(tmpNode, tmpV);	// add this node to path
					
						// merge, if target node has already a fromList
					if (outs[j].getHead().getFromList()!= null) {
						HashMap mList= outs[j].getHead().getFromList();
						SpliceNode[] mKeys= (SpliceNode[]) gphase.tools.Arrays.toField(mList.keySet().toArray());
						if (mKeys!= null) { 
							Arrays.sort(mKeys, new SpliceNode.PositionComparator());	// iterate backwards
							for (int k = mKeys.length- 1; k >= 0; --k) {	// merge pathes
								Vector mV= (Vector) mList.get(mKeys[k]);
								Vector v= (Vector) newMap.get(mKeys[k]);
								if (v== null) 
									newMap.put(mKeys[k], mV);
								else {	// here bubbles are closed !!
									mKeys[k].getOutEdges();
									outs[j].getHead().getInEdges();
									
									for (int h = 0; h < mV.size(); h++) 
										v.add(mV.elementAt(h));
									SpliceBubble blob= new SpliceBubble(mKeys[k], outs[j].getHead(), 
											(SplicePath[]) gphase.tools.Arrays.toField(v));
									if (bubbles== null) 
										bubbles= new SpliceBubble[] {blob};
									else 
										checkRedundancy(blob);	// check for REAL contained
								}								
							}
						}
					}
					outs[j].getHead().setFromList(newMap);	// set merged fromList
				}
			}	// end for all nodes				
		}
		return bubbles;
	}
	
	SpliceBubble[] cascadeBubbles(SpliceBubble[] bubs) {
		if (bubs== null)
			return null;
		
			// establish bubble hierachy
		Arrays.sort(bubs, new SpliceBubble.SizeComparator());
		int lastSize= -1;
		int pos= -1;
		for (int i = 0; i < bubs.length; i++) {
			if (lastSize== -1|| bubs[i].getSize()> lastSize) {
				lastSize= bubs[i].getSize();
				pos= i;
			}
			for (int j = pos; j < bubs.length; j++) {
				if (i== j)
					continue;
				if (bubs[j].containsGeometrically(bubs[i])) {
					SpliceBubble[] c1= bubs[i].getChildren();	// remove double represented children
					SpliceBubble[] c2= bubs[j].getChildren();
					for (int k = 0; c1!= null&& k < c1.length; k++) 
						for (int m = 0; c2!= null&& m < c2.length; m++) 
							if (c1[k]== c2[m]) {
								bubs[j].removeChild(c2[m]);
								c2[m].removeParent(bubs[j]);
							}
					bubs[i].addParent(bubs[j]);
					bubs[j].addChild(bubs[i]);
					if (bubs[i]== bubs[j])
						System.out.println("error");
				}
			}
		}
		
			// collect leaves
		Vector bubV= new Vector();
		for (int i = 0; i < bubs.length; i++) 
			if (!bubs[i].hasChildren())
				bubV.add(bubs[i]);
		SpliceBubble[] leafBubs= (SpliceBubble[]) gphase.tools.Arrays.toField(bubV);
	
			// create intersections
		Vector interBubV= new Vector();
		for (int i = 0; i < leafBubs.length; i++) {
			for (int j = i+1; j < leafBubs.length; j++) {	// try every 2 hierarchies once
				intersectBubbles_bottomUp(leafBubs[i], leafBubs[j], interBubV);
			}
		}

		// find top-level bubbles
		bubV= new Vector();
		for (int i = 0; i < bubs.length; i++) 
			if (!bubs[i].hasParents())
				bubV.add(bubs[i]);
		SpliceBubble[] topBubs= (SpliceBubble[]) gphase.tools.Arrays.toField(bubV);
		return topBubs;
	}
	
	void intersectBubbles_bottomUp(SpliceBubble bub, Vector chkBub) {
		if (!bub.hasParents())
			return;
		for (int i = 0; i < bub.getParents().length; i++) {
			intersectBubbles_bottomUp(bub, bub.getParents()[i], chkBub);
		}
	}
	
	void intersectBubbles_bottomUp(SpliceBubble bub, SpliceBubble blob, Vector chkBub) {
		if (!blob.hasParents())
			return;
		
		intersect(bub, blob, chkBub);
		
		for (int i = 0; blob.getParents()!= null&& i < blob.getParents().length; i++) 
			intersectBubbles_bottomUp(bub, blob.getParents()[i], chkBub);
		for (int i = 0; bub.getParents()!= null&& i < bub.getParents().length; i++) 
			intersectBubbles_bottomUp(bub.getParents()[i], blob, chkBub);
		
//		Transcript[][] bubT= bub.getTranscriptPartitions();
//		SpliceBubble[][] bubAnc= bub.getAncestors();
//		for (int i = 0; i < bubAnc.length; i++) {	// all pathes
//			Vector leafV= new Vector(bubAnc.length);
//			for (int j = 0; j < bubAnc[i].length; j++) {
//				Transcript[][] bubAncT= bubAnc[i][j].getTranscriptPartitions();
//				for (int k = (j+1); k < bubAnc[i].length; k++) {	// intersect upward
//					if (SpliceBubble.containedByTranscript(bubT, bubAncT)) {
//						SpliceBubble.intersect(bubAnc[i][j], bubAnc[i][k]);		
//					} else {	// hang in as sister of parent
//						for (int m = 0; bubAnc[i][j].getParents()!= null&& m < bubAnc[i][j].getParents().length; m++) 
//							bubAnc[i][j].getParents()[m].removeChild(bubAnc[i][j]);
//						bubAnc[i][j].setParents(bubAnc[i][k].getParents());
//					}
//					
//				}
//			}
//		}
		
	}

	void intersectBubbles_topDown(SpliceBubble[] sisters, Vector chkBubV) {
		if (sisters== null)
			return;
		
		for (int i = 0; i < sisters.length; i++) {
			for (int j = i+1; j < sisters.length; j++) {
				if (SpliceGraph.intersect(sisters[i], sisters[j], chkBubV)) {
					for (int k = 0; k < sisters[i].getChildren().length; k++) {
						SpliceGraph.intersect(sisters[j], sisters[i].getChildren()[k], chkBubV);
						for (int m = 0; m < sisters[j].getChildren().length; m++) 
							SpliceGraph.intersect(sisters[i].getChildren()[k], sisters[j].getChildren()[m], chkBubV);
					}
					for (int k = 0; k < sisters[j].getChildren().length; k++) 
						SpliceGraph.intersect(sisters[i], sisters[j].getChildren()[k], chkBubV);
				} 
				intersectBubbles(sisters[i].getChildren(), chkBubV);		
				intersectBubbles(sisters[j].getChildren(), chkBubV);		// TODO inefficient, iterated multiple times
			}
		}
	}
	
	SpliceBubble[] findBubbles() {
		Vector bubV= new Vector();
		SpliceNode[] nodes= getNodeList();
		for (int i = 0; i < nodes.length; i++) {
			
			SpliceNode tmpNode= nodes[i];
			HashMap map= tmpNode.getFromList();
			
				// opening bubbles: add node to fromList
			SpliceEdge[] outs= tmpNode.getOutEdges();
			for (int j = 0; outs!= null&& j < outs.length; j++) {	// outs.length> 1?
					
				SpliceNode head= outs[j].getHead(); 
				
					// extend existing pathes, create new path from this node
				HashMap newMap= new HashMap(map.size());
				SpliceNode[] fromNodes= (SpliceNode[]) gphase.tools.Arrays.toField(map.keySet().toArray());
				for (int k = 0; fromNodes!= null&& k < fromNodes.length; k++) {	// remove orphan fromNodes
					Vector pathes= (Vector) map.get(fromNodes[k]);
					Vector newPathes= new Vector();
					for (int h = 0; h < pathes.size(); h++) {
						SplicePath p= ((SplicePath) pathes.elementAt(h)).exendPath(outs[j]);
						if (p!= null)
							newPathes.add(p);
					}
					if (newPathes.size()> 0)	// fromNode no longer valid, remove
						newMap.put(fromNodes[k], newPathes);
				}
				SplicePath sp= new SplicePath(outs[j]);
				Vector tmpV= new Vector();
				tmpV.add(sp);
				newMap.put(tmpNode, tmpV);	// add this node to path
				
					// merge, if target node has already a fromList
				if (head.getFromList()!= null) {	
					HashMap mList= head.getFromList();		// nodes with pathes origins to target node
					SpliceNode[] mKeys= (SpliceNode[]) gphase.tools.Arrays.toField(mList.keySet().toArray());
					Vector partV= new Vector();	// partitions
					if (mKeys!= null) { 
						Arrays.sort(mKeys, new SpliceNode.PositionComparator());
						for (int k = mKeys.length- 1; k >= 0; --k) {	// iterate backwards over from nodes already stored in target
							Vector mV= (Vector) mList.get(mKeys[k]);		// pathes from a node already leading to target node
							Vector v= (Vector) newMap.get(mKeys[k]);		// new pathes from same node to target node
							if (v== null) 
								newMap.put(mKeys[k], mV);
							else {
								for (int h = 0; h < mV.size(); h++)		// merge pathes 
									v.add(mV.elementAt(h));

									// condition for creating bubbles (maximality of pathes)
//								if (head.getInDegree()== ((Integer) tab.get(head)).intValue()) {
										// collect data
									SplicePath[] pathes= (SplicePath[]) gphase.tools.Arrays.toField(v);
									Transcript[][] pathT= new Transcript[pathes.length][];
									for (int n = 0; n < pathT.length; n++) 
										pathT[n]= pathes[n].getTranscripts();
									
										// check bubble for minimal boundaries (wrt start)
									if (partV.size()> 0) {
										int n;
										for (n = 0; n < partV.size(); n++) 
											if (SpliceBubble.contained(pathT, (Transcript[][]) partV.elementAt(n)))
												break;
										if (n< partV.size())
											continue; // skip bubble (smaller bubble w same transcript set already found)
									}
									
									partV.add(pathT);
									SpliceBubble blob= new SpliceBubble(mKeys[k], head, pathes);
										// remove partial bubbles
									boolean add= true;
									for (int m = 0; m < bubV.size(); m++) {
										SpliceBubble tmpBub= (SpliceBubble) bubV.elementAt(m);
										if (tmpBub.getSource().getSite().getPos()== blob.getSource().getSite().getPos()&& 
												tmpBub.getSink().getSite().getPos()== blob.getSink().getSite().getPos()) {
											if (tmpBub.getSource()== blob.getSource()&& tmpBub.getSink()== blob.getSink()) 
												bubV.remove(m); 
											else {	// AbstractSite vc SS
												if (SpliceBubble.contained(blob.getTranscriptPartitions(), tmpBub.getTranscriptPartitions())) 
													bubV.remove(m);
												else
													add= false;
											}
											break;
										}
									}
									if (add)
										bubV.add(blob);
//								}
							}								
						}
					}
				}
				head.setFromList(newMap);	// set merged fromList
			}
		}	// end for all nodes				

		if (bubV.size()> 0)
			return (SpliceBubble[]) gphase.tools.Arrays.toField(bubV);
		return null;
	}

	public SpliceBubble[] getBubbles_old(boolean nonredundant) {
		
		if (bubbles== null) {
			SpliceNode[] nodes= getNodeList();
			for (int i = 0; i < nodes.length; i++) {
				int inDegree= nodes[i].getInDegree(false);
					// closing bubbles
				if (inDegree> 1) {
					HashMap[] fromLists= new HashMap[inDegree];
					nodes[i].setFromList(new HashMap());
					for (int j = 0; j < inDegree; j++) {
						fromLists[j]= nodes[i].getInEdges()[j].getTail().getFromList();
						Iterator iter= fromLists[j].keySet().iterator();
						while (iter.hasNext()) {
							Object key= iter.next();
							nodes[i].getFromList().put(key, new Integer(0));	// put with collisions, create superset
						}
					}
					if (nodes[i].getFromList().size()> 0) {
						//int[] pos= gphase.tools.Arrays.toPrimitive((Integer[]) gphase.tools.Arrays.toField(nodes[i].getFromList().keySet().toArray()));
						AbstractSite[] sites= (AbstractSite[]) gphase.tools.Arrays.toField(nodes[i].getFromList().keySet().toArray()); 
						Arrays.sort(sites, new AbstractSite.PositionComparator());
						//for (int j = pos.length- 1; j >= 0; --j) {
						for (int j = sites.length- 1; j >= 0; --j) {
							Object key= sites[j];
							int val= 0;
							for (int x = 0; x < fromLists.length; x++) {	// unite lists, add up
								Integer addI= ((Integer) fromLists[x].get(key));
								if (addI== null)
									continue;
								int add= addI.intValue();
								if (val> 0) {	// always holds add> 0 (otherwise not in map)
									SpliceNode src= (SpliceNode) nodeList.get(key);
									SpliceBubble blob= new SpliceBubble(src, nodes[i]);	
									if (bubbles== null)	// add to bubble set
										bubbles= new SpliceBubble[] {blob};
									else {
										int p= Arrays.binarySearch(bubbles, blob, new SpliceBubble.PositionComparator());
										if (p< 0)
											bubbles= (SpliceBubble[]) gphase.tools.Arrays.insert(bubbles, blob, p);
										if (nonredundant)
											checkRedundancy(blob);
									}
									
								}
								val+= add;
							}
								// eventually remove src, when all outgoing transcripts cross through that point
							Transcript[] srcTrans= ((SpliceNode) nodeList.get(key)).getSite().getTranscripts();
							Transcript[] snkTrans= nodes[i].getSite().getTranscripts();
							int x;
							for (x = 0; x < srcTrans.length; x++) {
								int y;
								for (y = 0; y < snkTrans.length; y++) 
									if (snkTrans[y]== srcTrans[x])
										break;
								if (y== snkTrans.length)	// not found
									break;
							}
							if (x< srcTrans.length) 
								nodes[i].getFromList().put(key, new Integer(val));
							else
								nodes[i].getFromList().remove(key);
						}
		//				for (int j = 0; j < fromLists.length; j++) {	// unite lists, add up
		//					Iterator iter= fromLists[j].keySet().iterator();
		//					while (iter.hasNext()) {
		//						Object key= iter.next();
		//						int count= ((Integer) fromList.remove(key)).intValue();
		//						int add= ((Integer) fromLists[j].get(key)).intValue();
		//						int ncount= count+ add;
		//						if (count> 0) {	// always holds add> 0 (otherwise not in map)
		//							SpliceNode src= (SpliceNode) nodeList.get(key);
		//							SpliceBubble blob= new SpliceBubble(src, nodes[i]);	
		//							if (bubbles== null)	// add to bubble set
		//								bubbles= new SpliceBubble[] {blob};
		//							else {
		//								int x= Arrays.binarySearch(bubbles, blob, new SpliceBubble.PositionComparator());
		//								if (x< 0)
		//									bubbles= (SpliceBubble[]) gphase.tools.Arrays.insert(bubbles, blob, x);
		//								checkRedundancy(bubbles, blob);
		//							}
		//							if (src.getOutDegree(true)> ncount) // remove from list, if all outgoing sc closed
		//								fromList.put(key, new Integer(count+ add));	
		//						}
		//					}
		//				}
					}
				} else {	// init now done in SpliceNode.getFromList()
//					if (nodes[i].getInEdges()== null)	// tss
//						nodes[i].setFromList(new HashMap());
//					else
					if (nodes[i].getInEdges()!= null&& nodes[i].getInEdges().length> 0)
						nodes[i].setFromList((HashMap) nodes[i].getInEdges()[0].getTail().getFromList().clone());	// copy fromList from parent
				}
				
					// opening bubbles: add node to fromList
				if (nodes[i].getOutDegree(false)> 1) 
					nodes[i].getFromList().put(nodes[i].getSite(), new Integer(1));	// 1 connection between node and itself
			}
		}
		
		return bubbles;
	}

	public SplicePath findPath(SpliceNode pSrc, SpliceNode pSnk, Transcript[] pPart) {
		
		SpliceNode tmpNode= pSrc;
		SplicePath path= null;
		while (tmpNode!= pSnk) {
			SpliceEdge[] oEdges= tmpNode.getOutEdges();
			Vector edgeV= new Vector();
			for (int j = 0; j < oEdges.length; j++) 
				if (SpliceBubble.contained(pPart, oEdges[j].getTranscripts()))
					edgeV.add(oEdges[j]);
				// chk if snk is an abstract site
			SpliceEdge tmpEdge= null;
			int i;
			SpliceEdge[] oldVecData= (SpliceEdge[]) gphase.tools.Arrays.toField(edgeV);
			for (i = 0; i < oldVecData.length; i++) {
				tmpEdge= oldVecData[i];
				if (!(tmpEdge.getHead().getSite() instanceof SpliceSite)) {
					if (tmpEdge.getHead()== pSnk)
						break;	// sink reached
					else
						edgeV.remove(tmpEdge);	// delete abstract sites with are not sink (no out edges anyway)
				}
			}
			if (i== oldVecData.length) {	// no abstract site as sink found
				if (edgeV.size()!= 1)
					System.err.println("Assertion failed: ambigous path");
				else
					tmpEdge= (SpliceEdge) edgeV.elementAt(0);
			}
			
			if (path== null) 
				path= new SplicePath(tmpEdge);
			else
				path= path.exendPath(tmpEdge);
			tmpNode= tmpEdge.getHead();
		}
		
		if (tmpNode== pSnk)
			return path;
		return null;
	}
	
	public boolean intersect(SpliceBubble bub0, SpliceBubble bub1, Vector chkBubV) {
		
		if ((bub0== bub1)|| bub0.isIBubble()|| bub1.isIBubble())
			return false;
		
			// same anchors
		if (bub0.getSource().getSite().getPos()== bub1.getSource().getSite().getPos()&&
				bub0.getSink().getSite().getPos()== bub1.getSink().getSite().getPos())
			return false;
		
			// no intersection
		if ((bub0.getSink().getSite().getPos()< bub1.getSource().getSite().getPos())||
				(bub1.getSink().getSite().getPos()< bub0.getSource().getSite().getPos()))
			return false;
		
			// find intersection nodes
		SpliceNode iSrc, iSnk;
		iSrc= (Math.abs(bub0.getSource().getSite().getPos())> Math.abs(bub1.getSource().getSite().getPos()))?
				bub0.getSource():bub1.getSource();	// maxPos of srcs
		iSnk= (Math.abs(bub0.getSink().getSite().getPos())< Math.abs(bub1.getSink().getSite().getPos()))?
				bub0.getSink():bub1.getSink();	// maxPos of srcs
				
				
			// find intersection partitions
		Transcript[][] part0= bub0.getTranscriptPartitions();
		Transcript[][] part1= bub1.getTranscriptPartitions();
		Vector v= new Vector();	// collects transcripts found in the intersection
		for (int i = 0; i < part0.length; i++) {
			for (int j = 0; j < part1.length; j++) {
				Transcript[][] iPart= SpliceBubble.intersect(part0[i], part1[j]);	// [0] sobre tiny, [1] sobre wide, [2] intersect
				if (iPart!= null&& iPart[2]!= null) 
						v.add(iPart[2]);	//v= (Vector) Arrays.addUnique(v, iPart[2]);	// intersection
			}
		}
		
			// get pathes
		Transcript[][] iPart= (Transcript[][]) gphase.tools.Arrays.toField(v);
		if (iPart== null|| iPart.length< 1)
			return false;
		SplicePath[] iPathes= new SplicePath[iPart.length];
		for (int i = 0; i < iPart.length; i++) {
			iPathes[i]= findPath(iSrc, iSnk, iPart[i]);
		}
			// no path
		if (iPathes== null|| iPathes[0]== null|| iPathes[0].getNodeV().size()== 0)
			return false;
			// redundancy filter pathes
		Vector iPathesV= new Vector();
		for (int i = 0; i < iPathes.length; i++) {
			int j;
			for (j = i+1; j < iPathes.length; j++) {
				Transcript[] t1= iPathes[i].getTranscripts();
				Transcript[] t2= iPathes[j].getTranscripts();
				if (t1.length!= t2.length)
					continue;
				
				Transcript[][] it= SpliceBubble.intersect(t1, t2);
				if (it== null||((it[0]== null|| it[0].length== 0)&& (it[1]== null|| it[1].length== 0)))
					break;	// identical path found
			}
			if (j== iPathes.length)
				iPathesV.add(iPathes[i]);
		}
		iPathes= (SplicePath[]) gphase.tools.Arrays.toField(iPathesV);
		
			// create bubble
		SpliceBubble interBub= new SpliceBubble(iSrc, iSnk, iPathes);
		interBub.setIBubble(true);
		int x;	// redundancy check
		for (x = 0; x < chkBubV.size(); x++) 
			if (interBub.equals(chkBubV.elementAt(x))) {
				interBub= (SpliceBubble) chkBubV.elementAt(x);
				break;
			}
		if (x== chkBubV.size()) {
			chkBubV.add(interBub);
			return false;
		}
		
			// insert into hierarchy
//		if (iPathes.length> 0) {
////			SpliceBubble[] c= bub0.getChildren();
////			for (int i = 0; c!= null&& i < c.length; i++) {
////				if (interBub== c[i])
////					continue;
////				if (interBub.contains(c[i])) {
////					try {
////						c[i].removeParent(bub0);
////					} catch (Exception e) {;}
////					interBub.addChild(c[i]);
////					c[i].addParent(interBub);
////					bub0.removeChild(c[i]);
////				}
////			}
//			if (interBub!= bub0) {
//				bub0.addChild(interBub);
//				interBub.addParent(bub0);
//			}
//			
////			c= bub1.getChildren();
////			for (int i = 0; c!= null&& i < c.length; i++)  {
////				if (interBub== c[i])
////					continue;
////				if (interBub.contains(c[i])) {
////					try {
////						c[i].removeParent(bub1);
////					} catch (Exception e) {;}
////					interBub.addChild(c[i]);
////					c[i].addParent(interBub);
////					bub1.removeChild(c[i]);
////				}
////			}
//			if (interBub!= bub1) {
//				bub1.addChild(interBub);
//				interBub.addParent(bub1);
//			}
//		} else 
//			return false;
		
		return true;
	}

	/**
		 * trims splice chains from both sides to the first splice site
		 * that is covered by exonic positions in all other splice chains.
		 *  
		 * @param vars
		 * @return
		 */
		public static int[] trim(Transcript[] trans) {
			
			SpliceSite[][] vars= new SpliceSite[trans.length][];
			for (int i = 0; i < vars.length; i++) 
				vars[i]= trans[i].getSpliceChain();
			
				// create splice universe
			Vector suVec= new Vector();
			for (int i = 0; i < vars.length; i++) {
				for (int j = 0; j < vars[i].length; j++) {
					int k;
					for (k = 0; k < suVec.size(); k++) 
						if (((SpliceSite) suVec.elementAt(k)).getPos()== vars[i][j].getPos())
							break;
					if (k== suVec.size())
						suVec.add(vars[i][j]);
				}
			}
			Vector tsVec= new Vector();	// add tss, tes
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].get5PrimeEdge())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].get5PrimeEdge()));
			}
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].get3PrimeEdge())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].get3PrimeEdge()));
			}
			int[] su= new int[suVec.size()+ tsVec.size()];	// convert to int[]
			for (int i = 0; i < suVec.size(); i++) 
				su[i]= ((SpliceSite) suVec.elementAt(i)).getPos();
			for (int i = 0; i < tsVec.size(); i++) 
				su[i+suVec.size()]= ((Integer) tsVec.elementAt(i)).intValue();
			
			Comparator compi= new AbstractSite.PositionComparator();
			Arrays.sort(su);
			int min= su[0];
			int max= su[su.length- 1];
			
				// trim end5
			while (true) {
				// check for exonic area
				int i;
				for (i = 0; i < vars.length; i++) {
					if (Transcript.getSpliceSiteByPos(vars[i], min)== null) {	
	//					if (vars[i]== null)	// if not containing edge site ----> bullshit, intron retention in 3'UTR
	//						break;
						SpliceSite pred= Transcript.getPredSpliceSite(vars[i], min);
						SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], min);
						if (pred== null) {
							if (succ== null) {
								if (!trans[i].contains(min))
									break;
							} else {
								if (trans[i].isUpstream(min)|| !succ.isDonor())
									break;
							}
						} else {
							if (succ== null) {
								if (!pred.isAcceptor()|| trans[i].isDownstream(min))
									break;
							} else {
								if (!pred.isAcceptor()|| !succ.isDonor())
									break;
							}
						}
					}
				}
				if (i== vars.length)	// finish: all exonic
					break;
						
				int s= Transcript.getSuccPos(su, min);		// containing the splice site, accepted, next
				if (s== 0) {
					min= 0;
					break;
				}
				min= s;
			}
					
			// trim end3
			while (true) {
				// check for exonic area
				int i;
				for (i = 0; i < vars.length; i++) {
					if (Transcript.getSpliceSiteByPos(vars[i], max)== null) {
	//					if (vars[i].length== 0)	// cannot succeed in containing ----> bullshit, intron retention in 3'UTR
	//						break;
						SpliceSite pred= Transcript.getPredSpliceSite(vars[i], max);
						SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], max);
						if (pred== null) {
							if (succ== null) {
								if (!trans[i].contains(max))
									break;
							} else {
								if (trans[i].isUpstream(max)|| !succ.isDonor())
									break;
							}
						} else {
							if (succ== null) {
								if (!pred.isAcceptor()|| trans[i].isDownstream(max))
									break;
							} else {
								if (!pred.isAcceptor()|| !succ.isDonor())
									break;
							}
						}
					}
				}
				if (i== vars.length)
					break;
						
				int s= Transcript.getPredPos(su, max);		// containing the splice site, accepted, next
				if (s== 0) {
					max= 0;
					break;
				}					
				max= s;
			}
						
			return new int[] {min, max};
		}
		
}
