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
import java.util.Iterator;
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
		System.currentTimeMillis();
		for (int i = x; i < bubbles.length; i++) {	// all bubbles after p have same start and more little ends
			
			if (bubbles[i].getSink().getSite().getPos()> blob.getSink().getSite().getPos())
				break;
			
				// in contained, check transcript set
			Transcript[][] tmpPart= bubbles[i].getTranscriptPartitions();
			if (SpliceBubble.contained(nuParts, tmpPart))
				addBlob= false; // blob has wider borders and equal transcript set, discard blob 
			else if (SpliceBubble.contained(tmpPart, nuParts)) {
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
	public ASMultiVariation[] getMultiVariations(int k) {
		
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
	
	public SpliceBubble[] getBubbles() {
		
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
