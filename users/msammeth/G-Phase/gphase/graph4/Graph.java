package gphase.graph4;

import java.io.FileWriter;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import com.sun.org.apache.xml.internal.utils.NodeVector;

import sun.reflect.ReflectionFactory.GetReflectionFactoryAction;

import gphase.io.gtf.GTFChrReader;
import gphase.model.ASEvent;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.sgraph.Bubble;
import gphase.tools.DoubleVector;

public class Graph {
	public static void main(String[] args) {
		_070808_test();
	}
	
	static void _070808_test() {
		// 
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_CDS.gtf		
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC0703.gtf
		gphase.tools.File file= new gphase.tools.File("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf");
		Species species= new Species("human");
		species.setGenomeVersion("hg18");
		boolean output= false; 
		
		long t0= System.currentTimeMillis();
		GTFChrReader reader= new GTFChrReader(file.getAbsolutePath());
		try {
			reader.read();
		} catch (Exception e1) { 
			e1.printStackTrace();
		}
		if (output)
			System.out.println("read: "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
		
		Gene[] g= reader.getGenes();
		int cnt= 0;
		HashMap<String, Integer> evMap= new HashMap<String, Integer>();
		long cumulGC= 0l, cumulEV= 0l;
		while (g!= null) {
			
			for (int i = 0; i < g.length; i++) {
				
				if (g[i].getTranscriptCount()== 1)
					continue;
				
				Graph gr= new Graph(g[i]);
				if (output) {
					System.out.print(">  t= "+g[i].getTranscriptCount());
				}
				long t1= System.currentTimeMillis();
				gr.init();
				long t2= System.currentTimeMillis();
				long dA= (t2-t1);
				cumulGC+= dA;
				if (output) {
					System.out.print(", n= "+gr.nodeHash.size()+", e="+gr.edgeHash.size()+". aufbau "+(dA/1000)+" sec, ");
					System.out.flush();
				}
				int oldSize= evMap.size();
				gr.getEventsByPathes(2, evMap);
				cnt+= evMap.size()- oldSize;
				long dB= (System.currentTimeMillis()- t2);
				cumulEV+= dB;
				if (output) 
					System.out.println("events "+(evMap.size()-oldSize)+", extract "+(dB/1000)+" sec.");
				
				g[i]= null;
//				System.gc();
//				Graph g2= new Graph(g[i]);
//				g2.constructGraph();
//				Vector ev2= g2.extractASevents(2);			
//				System.currentTimeMillis();
			}

			System.gc();
			long tx= System.currentTimeMillis();
			try {
				reader.read();
			} catch (Exception e1) {
				e1.printStackTrace();
			}
			if (output)
				System.out.println("read: "+((System.currentTimeMillis()- tx)/ 1000)+" sec.");
			g= reader.getGenes();
			System.gc();
			Thread.currentThread().yield();
		}
		
		System.out.println("found "+cnt+" events.");
		Object[] o= evMap.keySet().toArray();
		int[] size= new int[o.length];
		for (int i = 0; i < o.length; i++) 
			size[i]= evMap.get(o[i]).intValue();
		Vector v= new Vector();
		v.add(o);
		gphase.tools.Arrays.synchroneousSort(size, v);
		for (int i = 0; i < o.length; i++) {
			System.out.println(size[i]+"\t"+o[i]);
		}
		System.out.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec, construct "+cumulGC+" msec, extract "+cumulEV+" msec.");
	}
	
	public boolean isRoot(Node v) {
		if (v.equals(root))
			return true;
		return false;
	}

	public boolean isLeaf(Node v) {
		if (v.equals(leaf))
			return true;
		return false;
	}
	
	public static class TranscriptByNameComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			return ((Transcript) o1).getTranscriptID().compareTo(((Transcript) o2).getTranscriptID());
		}
	}
	public static boolean isNull(long[] sig) {
		for (int i = 0; i < sig.length; i++) 
			if (sig[i]!= 0l)
				return false;
		return true;
	}
	
	public static long[] intersect(long[] sig1, long[] sig2) {
		assert(sig1.length== sig2.length);
		long[] reSig= new long[sig1.length];
		for (int i = 0; i < sig2.length; i++) 
			reSig[i]= sig1[i]& sig2[i];
		return reSig;
	}
	
	public static long[] unite(long[] sig1, long[] sig2) {
		assert(sig1.length== sig2.length);
		long[] reSig= new long[sig1.length];
		for (int i = 0; i < sig2.length; i++) 
			reSig[i]= sig1[i]| sig2[i];
		return reSig;
	}
	
	public static long[] without(long[] a, long[] b) {
		long[] c= new long[a.length];
		for (int i = 0; i < c.length; i++) 
			c[i]= a[i]^ b[i];
		return c;
	}
	
	
	public static long[] unite(Vector<long[]> a) {
		if (a.size()== 1)
			return a.elementAt(0);
		long[] inter= unite(a.elementAt(0), a.elementAt(1));
		for (int i = 2; i < a.size(); i++) 
			inter= unite(inter, a.elementAt(i));
		return inter;
	}
	
	public static int intersect(Vector<long[]> ps1, Vector<long[]> ps2, Vector<long[]> interSet, boolean onlyFullIntersect) {
		int cntFullIntersect= 0;
		for (int k = 0; k < ps1.size(); k++) {
			for (int h = 0; h < ps2.size(); h++) {
				long[] sect= intersect(ps1.elementAt(k), ps2.elementAt(h));
				if (isNull(sect))
					continue;
				int fis= 0;
				if (equalSet(sect,ps1.elementAt(k)))
					fis= 1;
				if (interSet!= null) {
					if  (!onlyFullIntersect)
						interSet.add(sect);
					else if (fis== 1)
						interSet.add(ps1.elementAt(k));
				}
				cntFullIntersect+= fis;
				if (onlyFullIntersect&& fis== 1)
					break;
			}
		}
		return cntFullIntersect;
	}
	
	public static boolean equalSet(long[] a, long[] b) {
		for (int i = 0; i < b.length; i++) {
			if (a[i]!= b[i])
				return false;
		}
		return true;
	}
	
	
	Vector<Path[]> generateTuples(int k, Vector<Vector<Path>> buckets) {
	
			int[] idx= new int[k];
			for (int i = 0; i < idx.length; i++) 
				idx[i]= i;
			
			int[] combi= new int[k];
	//		int permutNb= 1;
	//		for (int i = 0; i < combi.length; i++) 
	//			;
			Vector<Path[]> tuples= new Vector<Path[]>();
			while (idx[0]< (buckets.size()-k+ 1)) {
				
					// now all combinations
				for (int j = 0; j < combi.length; j++) 
					combi[j]= 0;
				while (true) {
					
					Path[] p= new Path[k];
					for (int j = 0; j < combi.length; j++) 
						p[j]= buckets.elementAt(idx[j]).elementAt(combi[j]);
					tuples.add(p);
					
					int negPtr= combi.length-1;
					while (negPtr>= 0) {
						if (combi[negPtr]< (buckets.elementAt(idx[negPtr]).size()-1)) {
							++combi[negPtr];
							break;
						} 
							// else
						combi[negPtr]= 0;
						--negPtr;
					}
					if (negPtr< 0)
						break;
				}
				
				//
				
				int negPtr= idx.length-1;
				while (negPtr>= 0) {
					if (idx[negPtr]< (buckets.size()-1)) {
						int c= ++idx[negPtr];
						for (int i = negPtr+1; i < idx.length; i++) 
							idx[negPtr]= ++c;
						break;
					} 
						// else
					--negPtr;
				}
			}
	
			return tuples;
		}

	Vector<long[][]> generateTuples(int k, long[][][] buckets) {	//Vector<Vector<Path>> buckets

		int[] idx= new int[k];
		for (int i = 0; i < idx.length; i++) 
			idx[i]= i;
		
		int[] combi= new int[k];
//		int permutNb= 1;
//		for (int i = 0; i < combi.length; i++) 
//			;
		Vector<long[][]> tuples= new Vector<long[][]>();
		while (idx[0]< (buckets.length-k+1)) {
			
				// now all combinations
			for (int j = 0; j < combi.length; j++) 
				combi[j]= 0;
			while (true) {
				
				long[][] p= new long[k][];
				for (int j = 0; j < combi.length; j++) 
					p[j]= buckets[idx[j]][combi[j]];
				tuples.add(p);
				
				int negPtr= combi.length-1;
				while (negPtr>= 0) {
					if (combi[negPtr]< (buckets[idx[negPtr]].length-1)) {
						++combi[negPtr];
						break;
					} 
						// else
					combi[negPtr]= 0;
					--negPtr;
				}
				if (negPtr< 0)
					break;
			}
			
			//
			
			int negPtr= idx.length-1;
			while (negPtr>= 0) {
				if (idx[negPtr]< (buckets.length-1)) {
					int c= ++idx[negPtr];
					for (int i = negPtr+1; i < idx.length; i++) 
						idx[negPtr]= ++c;
					break;
				} 
					// else
				--negPtr;
			}
		}

		return tuples;
	}

	Vector<Path[]> generateTuples(int k, int x, int[] a, Vector<Path>buckets, Vector<Path[]> tuples) {
		if (k== a.length) {
			
			HashMap<Integer,Integer> bubbleHash= new HashMap<Integer, Integer>();
			int start= buckets.elementAt(0).getSecondNode().getSite().getPos();
			int end= buckets.elementAt(0).getSecondLastNode().getSite().getPos();
			boolean hasSS= false, overlap= true;
			for (int z= 0; z < a.length; z++) {
					// check for bubble
				int[] bids= buckets.elementAt(a[z]).getBubbleIDs();
				for (int i = 0; i < bids.length; i++) {
					Integer id= new Integer(bids[i]);
					Integer val= bubbleHash.get(id);
					if (val== null)
						val= new Integer(0);
					val= new Integer(val.intValue()+1);
					bubbleHash.put(id, val);
				}
				
					// every has to overlap with each other
				if (z!= 0&& isRoot(buckets.elementAt(0).getSource())&&
						isLeaf(buckets.elementAt(0).getSink())) {
					Node sec= buckets.elementAt(a[z]).getSecondNode();
					Node secLast= buckets.elementAt(a[z]).getSecondLastNode();
					if (end< sec.getSite().getPos()|| start> secLast.getSite().getPos()) {
						overlap= false;
						break;
					}
				}
				
				if (buckets.elementAt(a[z]).hasSpliceSite())
					hasSS= true;
			}
			if ((!hasSS)|| (!overlap))
				return tuples;
			
			boolean bubble= false;
			Object[] o= bubbleHash.values().toArray();
			for (int i = 0; i < o.length; i++) {
				if (((Integer) o[i]).intValue()== a.length) {
					bubble= true;
					break;
				}
			}
			if (bubble)
				return tuples;

				// add
			Path[] tuple= new Path[a.length];
			for (int i = 0; i < a.length; i++) 
				tuple[i]= buckets.elementAt(a[i]);				
			tuples.add(tuple);
			return tuples;
		}
		
		for (int i = x; i < buckets.size()-(a.length- k- 1); i++) {
			a[k]= i;
			generateTuples(k+1,i+1,a,buckets,tuples);
		}
		return tuples;
	}
	
	static TranscriptByNameComparator defaultTranscriptByNameComparator= new TranscriptByNameComparator();
	static Node.PositionTypeComparator defaultNodeByPositionTypeComparator= new Node.PositionTypeComparator();
	
	Gene gene;
	Transcript[] trpts;
	int taSize;
	HashMap<SpliceSite, Node> nodeHash= new HashMap<SpliceSite, Node>();
	HashMap<Integer, Edge> edgeHash= new HashMap<Integer, Edge>();
	Node root, leaf;
	Node[] nodesInGenomicOrder= null;
	
	public Graph(Gene g) {
		this.gene= g;
		trpts= g.getTranscripts();
		Arrays.sort(trpts, defaultTranscriptByNameComparator);
		taSize= trpts.length/ 64;
		if (trpts.length%64!= 0)
			++taSize;
	}

	long[] encodeTset(Transcript[] t) {
		long[] taVector= new long[taSize];
		for (int i = 0; i < taVector.length; i++) 
			taVector[i]= 0l;
		for (int i = 0; i < t.length; i++) {
			int p= Arrays.binarySearch(trpts, t[i], defaultTranscriptByNameComparator);
			assert(p>= 0);
			int cnt= 0;
			while(p>= 64) {
				p-= 64;
				++cnt;
			}
			taVector[cnt]|= gphase.tools.Math.pow(2,p);
		}
		return taVector;
	}
	
	
	
	public static int getTranscriptNb(long[] c) {
		int cnt= 0;
		for (int i = 0; i < c.length; i++) {
			int base= 0;
			long val= 1;
			int sum= 0;
			while(base< 64) {
				
				if ((c[i]&val)!= 0l) {
					++cnt;
					sum+= val;
					if (sum== c[i])
						break;
				}
				++base;
				val*= 2;
			}
		}
		return cnt;
	}
	
	public Node getNode(SpliceSite ss) {
		return nodeHash.get(ss);
	}
	
	public Edge getEdge(Node v, Node w) {
		return edgeHash.get(new Integer(v.getSite().getPos()+w.getSite().getPos()));
	}
	
	public Node createNode(SpliceSite ss) {
		Node n= getNode(ss);
		if (n== null) {
			n= new Node(ss, encodeTset(ss.getTranscripts()));
			nodeHash.put(ss,n);
		}
		return n;
	}
	
	public Edge createEdge(Node v, Node w, long[] newTset) {
		Edge e= getEdge(v,w);
		if (e== null) {
			e= new Edge(v,w);
			e.setTranscripts(newTset);
			edgeHash.put(new Integer(e.getTail().getSite().getPos()+e.getHead().getSite().getPos()),e);
		} else 
			e.setTranscripts(unite(e.getTranscripts(),newTset));
		return e;
	}

	public Edge addEdge(Node v, Node w, long[] newTset) {
		Edge e= new Edge(v,w);
		e.setTranscripts(newTset);
		edgeHash.put(new Integer(e.getTail().getSite().getPos()+e.getHead().getSite().getPos()),e);
		return e;
	}
	
	public void init() {
		constructGraph();
		contractGraph(root);
	}
	
	void constructGraph() {
		Transcript[] t= gene.getTranscripts();
		HashMap<Node, Node> rootHash= new HashMap<Node, Node>(), leafHash= new HashMap<Node, Node>();
		for (int i = 0; i < t.length; i++) {
			SpliceSite[] ss= t[i].getSpliceSitesAll();
			for (int j = 1; j < ss.length; j++) {
				Node v= createNode(ss[j-1]);
				if (j== 1&& rootHash.get(v)== null)
					rootHash.put(v,v);
				Node w= createNode(ss[j]);
				if (j== ss.length- 1&& leafHash.get(v)== null)
					leafHash.put(w,w);
				createEdge(v,w,encodeTset(new Transcript[] {t[i]}));
			}
		}
		
		SpliceSite ss= new SpliceSite(Integer.MIN_VALUE, SpliceSite.TYPE_NOT_INITED);
		ss.setTranscripts(trpts);
		root= createNode(ss);
		Object[] o= rootHash.keySet().toArray();
		
			// collect transcripts where o[i] is first ss
		for (int i = 0; i < o.length; i++) { 
			Vector<Transcript> v= new Vector<Transcript>();
			for (int j = 0; j < trpts.length; j++) {
				if (trpts[j].getSpliceSitesAll()[0]== ((Node) o[i]).getSite())
					v.add(trpts[j]);
			}
			Transcript[]  tt= new Transcript[v.size()];
			for (int j = 0; j < tt.length; j++) 
				tt[j]= v.elementAt(j);
			createEdge(root, (Node) o[i], encodeTset(tt));
		}
		
		ss= new SpliceSite(Integer.MAX_VALUE, SpliceSite.TYPE_NOT_INITED);
		ss.setTranscripts(trpts);
		leaf= createNode(ss);
		o= leafHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) { 
			Vector<Transcript> v= new Vector<Transcript>();
			for (int j = 0; j < trpts.length; j++) {
				if (trpts[j].getSpliceSitesAll()[trpts[j].getSpliceSitesAll().length-1]== ((Node) o[i]).getSite())
					v.add(trpts[j]);
			}
			Transcript[]  tt= new Transcript[v.size()];
			for (int j = 0; j < tt.length; j++) 
				tt[j]= v.elementAt(j);
			createEdge((Node) o[i], leaf, encodeTset(tt));
		}
	}
	
	//DFS
	void contractGraph_new_notWork(Node n) {
		if (n.isProcessed())
			return;
		n.setProcessed(true);
		
		Set<Edge> outV= n.getOutEdges();
		Set<Edge> inV= n.getInEdges();
		Path p= null;
		Vector<Edge> edgeV= new Vector<Edge>();
		while (outV.size()== 1&& (p== null|| inV.size()== 1)) {
			Edge e= outV.iterator().next();
			if (p== null)
				p= new Path(n);
			else
				p.addNode(n);
			edgeV.add(e);
			n= e.getHead();
			outV= n.getOutEdges();
			inV= n.getInEdges();			
		}
		
		if (p!= null&& edgeV.size()> 1) {
			p.addNode(edgeV.elementAt(edgeV.size()- 1).getHead());	// complete
			Vector ns= p.getNodesAndEdges();
			
			for (int i = 0; i < edgeV.size(); i++) 
				removeEdge(edgeV.elementAt(i));
			
			Edge e= null;
			if (p.getSource().getInEdges().size()!= 1) {
				for (int i = 1; i < ns.size()- 1; i++) 
					removeNode((Node) ns.elementAt(i));
				e= createEdge(p.getSource(), p.getSink());
				p.removeSource();
			} else { 
				e= createEdge(((Edge) p.getSource().getInEdges().toArray()[0]).getTail(), p.getSink());
				for (int i = 0; i < ns.size()- 1; i++) 
					removeNode((Node) ns.elementAt(i));
			}
			p.removeSink();
			e.setPath(p);
		}
		
		if (outV.size()> 1) {
			Object[] o= outV.toArray();
			for (int i = 0; i < o.length; i++) {
				contractGraph(((Edge) o[i]).getHead());
			}
		} else 
			contractGraph(n);
	}

	//DFS
	void contractGraph(Node n) {
		if (n.isProcessed())
			return;
		n.setProcessed(true);
		
		Vector<Edge> outV= n.getOutEdges();
		Vector<Edge> inV= n.getInEdges();
		Path p= null;
		Vector<Edge> edgeV= new Vector<Edge>();
		Vector<Node> nodeV= new Vector<Node>();
		while (outV.size()== 1&& inV.size()< 2) {
			Edge f= null;
			Iterator<Edge> it= inV.iterator();
			Edge e= outV.iterator().next();
			if (it.hasNext())
				f= it.next();
			if (p== null) {
				if (f!= null) {
					edgeV.add(f);
					nodeV.add(f.getTail());
				}
				p= new Path();
				if (f== null)
					p.setSourceEdge(e);
				else 
					p.setSourceEdge(f);
				p.setSinkEdge(e);
				p.setTranscripts(e.getTranscripts());
			} else
				p.setSinkEdge(e);
			
			edgeV.add(e);
			nodeV.add(n);
			n= e.getHead();
			outV= n.getOutEdges();
			inV= n.getInEdges();			
		}
		
			// contract
		if (edgeV.size()> 1) {
			for (int i = 1; i < nodeV.size(); i++) 
				removeNode((Node) nodeV.elementAt(i));	// not the first one
			
			for (int i = 0; i < edgeV.size(); i++)  {
				removeEdge(edgeV.elementAt(i));
			}
			
			Edge e= addEdge(p.getSourceNode(), p.getSinkNode(), p.getTranscripts());	// cannot use create, will return same edge for alt. ways
			e.setPath(p);
		}

		
		Object[] o= outV.toArray();
		for (int i = 0; i < o.length; i++) {
			contractGraph(((Edge) o[i]).getHead());
		}
	}
	
	void initPartitions(int n, HashMap<String, Integer> evMap) {		
			Node[] nodes= getNodesInGenomicOrder();
			Vector<Node> nodesWithOutPSS= new Vector<Node>(nodes.length/ 2);
			for (int i = 0; i < nodes.length; i++) {
				Set<Edge> inEdges= nodes[i].getInEdges();
				Set<Edge> outEdges= nodes[i].getOutEdges();
				if (inEdges.size()>= n) {
					PartitionSet inPartition= nodes[i].getInPartitionSet();
					for (int j = nodesWithOutPSS.size()- 1; j >= 0; --j) {
						PartitionSet outPartition= nodesWithOutPSS.elementAt(j).getOutPartitionSet();
						Vector<long[]> interPartVec= new Vector<long[]>(Math.max(outPartition.getPartitions().size(), inPartition.getPartitions().size()));
						int cntFullIS= intersect(outPartition.getPartitions(), inPartition.getPartitions(), interPartVec, false);
	
						// bubble:: and here will be the problem of all !!!
						if (interPartVec.size()>= n) {
								// retrieve event
							long[][][] splitPart= outPartition.getSplitPartitions();
							Vector<long[][]> tuples= generateTuples(n, splitPart);
							for (int k = 0; k < tuples.size(); k++) {
								
									// check AS event
								SpliceSite[][] ss= new SpliceSite[n][];
								Transcript[][] tt= new Transcript[n][];
								for (int h = 0; h < tt.length; h++) {
									tt[h]= decodeTset(tuples.elementAt(k)[h]);
									ss[h]= tt[h][0].getSpliceSitesBetween(nodesWithOutPSS.elementAt(j).getSite(),nodes[i].getSite());
								}
								for (int h = 0; h < tt.length; h++) {
									tt[h]= decodeTset(tuples.elementAt(k)[h]);
									ss[h]= tt[h][0].getSpliceSitesBetween(nodesWithOutPSS.elementAt(j).getSite(),nodes[i].getSite());
								}
								ASEvent ev= new ASEvent(tt,ss);
								boolean event= true;
								if (isRoot(nodesWithOutPSS.elementAt(j))|| isLeaf(nodes[i]))
									event= ev.isASevent();
								if (event) {
									Integer nr= evMap.get(ev.toString());
									if (nr== null)
										nr= new Integer(0);
									evMap.put(ev.toString(),new Integer(nr.intValue()+1));
	//								try {
	//									FileWriter writer= new FileWriter("new.asta", true);
	//									writer.write(ev.toStringASTA()+"\n");
	//									writer.flush(); writer.close();
	//								} catch (Exception e) {
	//									e.printStackTrace();
	//								}
								}
								
							}
							
								// merge partitions for this and aterior outPartitions (outer Bubbles)
							for (int k = j; k >= 0; --k) {
								PartitionSet ops= nodesWithOutPSS.elementAt(k).getOutPartitionSet();
								Vector<long[]> fullISset= new Vector<long[]>(Math.min(ops.getPartitions().size(),inPartition.getPartitions().size()));
								int cntFIS= intersect(ops.getPartitions(),interPartVec, fullISset, true);
								// (cntFIS== ops.getPartitions().size())  	// could be removed, but cannot happen
								if (cntFIS> 1) { 	// have to be merged
									ops.mergePartitions(fullISset);
									if (ops.getPartitions().size()< n) {
										nodesWithOutPSS.remove(k--);
										--j;
										continue;
									}
								}
								
							}
	
							// outPartition fully satisfied by Bubble, remove
							if (cntFullIS== outPartition.getPartitions().size()) {
								nodesWithOutPSS.remove(j--);
								break;
							}
						} // if bubble
					} // nodesWithOutPSS
					
				}
				
				if (outEdges.size()>= n) {
					nodesWithOutPSS.add(nodes[i]);
					
					// split pathes. Not necessary for bubble retrieval, but for collecting pathes without iterating the graph
					for (int j = 0; j < nodesWithOutPSS.size(); j++) {
						nodesWithOutPSS.elementAt(j).getOutPartitionSet().splitPathes(nodes[i].getInPartitionSet());
					}
				}
			}
		}

	public void getEventsByPathes(int n, HashMap<String, Integer> evMap) {		
		Node[] nodes= getNodesInGenomicOrder();
		Vector<Node> nodesWithOutPSS= new Vector<Node>(nodes.length/ 2);
		for (int i = 0; i < nodes.length; i++) {
			Vector<Edge> inEdges= nodes[i].getInEdges();
			Vector<Edge> outEdges= nodes[i].getOutEdges();
			if (inEdges.size()>= n) {
				
				Object[] fNodes= nodes[i].getFromNodes().toArray();
				Arrays.sort(fNodes,defaultNodeByPositionTypeComparator);
				for (int j = fNodes.length- 1; j >= 0; --j) {
					Node fromNode= (Node) fNodes[j];
					HashMap<Edge, Vector<Path>> fromPmap= nodes[i].getFromNodeMap().get(fNodes[j]);
					// bubble
					if (fromPmap.size()>= n) {

							// retrieve event
						Vector<Vector<Path>> splitPart= new Vector(nodes[i].getFromNodeMap().get(fNodes[j]).values());
						Vector<Path[]> tuples= generateTuples(n, splitPart);
						for (int k = 0; k < tuples.size(); k++) {
							
								// check AS event
							SpliceSite[][] ss= new SpliceSite[n][];
							Transcript[][] tt= new Transcript[n][];
							for (int h = 0; h < tt.length; h++) {
								tt[h]= decodeTset(tuples.elementAt(k)[h].getTranscripts());
								ss[h]= tt[h][0].getSpliceSitesBetween(
										tuples.elementAt(k)[h].getSourceNode().getSite(),nodes[i].getSite());
							}
							
							ASEvent ev= new ASEvent(tt,ss);
							boolean event= true;
							if (isRoot(((Node) fNodes[j]))|| isLeaf(nodes[i]))
								event= ev.isASevent();
							if (event) {
								Integer nr= evMap.get(ev.toString());
								if (nr== null)
									nr= new Integer(0);
								evMap.put(ev.toString(),new Integer(nr.intValue()+1));
								try {
									FileWriter writer= new FileWriter("new.asta", true);
									writer.write(ev.toStringASTA()+"\n");
									writer.flush(); writer.close();
								} catch (Exception e) {
									e.printStackTrace();
								}
							}
						}	// for all tuples
						
						// kill pathes
						nodes[i].getFromNodeMap().remove(fNodes[j]);
						if (((Node) fNodes[j]).getOutEdges().size()!= splitPart.size()) {
							Vector<Path> newPartition= new Vector<Path>();
							for (int k = 0; k < splitPart.size(); k++) 
								for (int m = 0; m < splitPart.elementAt(k).size(); m++) 
									newPartition.add(splitPart.elementAt(k).elementAt(m));
							HashMap<Edge, Vector<Path>> map = new HashMap<Edge, Vector<Path>>(1,1f);
							map.put(null,newPartition);
							nodes[i].getFromNodeMap().put((Node) fNodes[j],map);
						}
						
					} // if bubble
				} // nodesWithOutPSS
				
			}

			
			Object[] o= outEdges.toArray();
			// add possible new source
			if (outEdges.size()>= n) {
					// generate new pathes
				for (int j = 0; j < o.length; j++) {
					Edge e= (Edge) o[j];
					Path newP= new Path();
					newP.setSourceEdge(e);
					newP.setSinkEdge(e);
					newP.setTranscripts(e.getTranscripts());					
					HashMap<Edge, Vector<Path>> map= e.getHead().getFromNodeMap().remove(nodes[i]);
					if (map== null)
						map= new HashMap<Edge, Vector<Path>>(2,1f);
					Vector<Path> v= map.remove(e);
					if (v== null)
						v= new Vector<Path>(2);
					v.add(newP);
					map.put(e,v);
					e.getHead().getFromNodeMap().put(nodes[i], map);
				}
			}
							
			Iterator<Node> iterSrc= nodes[i].getFromNodes().iterator();
			while (iterSrc.hasNext()) {
				Node src= iterSrc.next();
				Iterator<Path> iterPath= nodes[i].getPathesFrom(src).iterator();	// not only one per source per inedge!
				// extend all pathes
				while (iterPath.hasNext()) {
					Path p= iterPath.next();
//					Vector<Path> split= new Vector<Path>(o.length);
					for (int j = 0; j < o.length; j++) {
						Edge e= (Edge) o[j];
						long[] newTsup= intersect(p.getTranscripts(), e.getTranscripts());
						if (isNull(newTsup))
							continue;
						Path newP= new Path();
						newP.setSourceEdge(p.getSourceEdge());
						newP.setTranscripts(newTsup);
						newP.setSinkEdge(e);
						HashMap<Edge, Vector<Path>> epMap= e.getHead().getFromNodeMap().get(p.getSourceNode());
						if (epMap== null) {
							epMap= new HashMap<Edge,Vector<Path>>(e.getHead().getInEdges().size(),1f);
							e.getHead().getFromNodeMap().put(p.getSourceNode(),epMap);
						}
						Vector<Path> v= epMap.get(e);
						if (v== null) {
							v= new Vector<Path>(1);
							epMap.put(e,v);
						}
						v.add(newP);
//						split.add(newP);
					}
					
						// split partitions if path separates
//					if (split.size()> 1) {
//						HashMap<Path,Path> partition= p.getSourceNode().getOutPartitionSet().getPartition(p.getSourceEdge());
//						partition.remove(p);	// need hashmap
//						for (int j = 0; j < split.size(); j++) 
//							partition.put(split.elementAt(j),split.elementAt(j));
//					}
				}
			}
			

		}
	}
	
	
	public void removeEdge(Edge e) {
		edgeHash.remove(new Integer(e.getTail().getSite().getPos()+ e.getHead().getSite().getPos()));
		e.getTail().removeOutEdge(e);
		e.getHead().removeInEdge(e);
	}
	
	public void removeNode(Node n) {
		nodeHash.remove(n.getSite());
	}
	
	public Node[] getNodesInGenomicOrder() {
		if (nodesInGenomicOrder == null) {
			Iterator<Node> iter= nodeHash.values().iterator();
			int cnt= 0;
			nodesInGenomicOrder= new Node[nodeHash.size()];
			while (iter.hasNext()) 
				nodesInGenomicOrder[cnt++]= iter.next();
			Arrays.sort(nodesInGenomicOrder,defaultNodeByPositionTypeComparator);
		}

		return nodesInGenomicOrder;
	}
	
	public void extractASevents(int n, HashMap<String, Integer> map) {
		
		Comparator vectorPathSourceCoordinateCompi= new Comparator() {
			@Override
			public int compare(Object o1, Object o2) {
				Vector<Path> v1= (Vector<Path>) o1;
				Vector<Path> v2= (Vector<Path>) o2;
				return (defaultNodeByPositionTypeComparator.compare(
						v1.elementAt(0).getSource(), 
						v2.elementAt(0).getSource()));
			}
		};
		
		HashMap<Node, PartitionSet> openPartitions= new HashMap<Node, PartitionSet>();
		
		Node[] nodes= getNodesInGenomicOrder();
		for (int i = 0; i < nodes.length; i++) {
			Set<Edge> inEdges= nodes[i].getInEdges();
			if (inEdges.size()>= n) {
				Object[] o= nodes[i].getBuckets().toArray();
				Arrays.sort(o, vectorPathSourceCoordinateCompi);
				for (int k = o.length- 1; k >= 0; --k) {	// close inner ones first
					Vector<Path> bucket= (Vector<Path>) o[k];
					if (bucket.size()< n)
						continue;
					Vector<Vector<Path>> partBucket= bucket.elementAt(0).getPartition().splitBucket(bucket); 
					if (bucket.size()>= n) {

						Vector<Path[]> tuples= generateTuples(n,partBucket);
						if (tuples.size()!= 0) { 
							Bubble bubble= new Bubble(this, bucket, tuples);
							bubble.extractEvents(map);
						}
						
						// does not make sense to extend
						if (bucket.size()== bucket.elementAt(0).getSource().getOutEdges().size()) { 
							if (!nodes[i].removeBucket(bucket))
								System.err.println("bucket not found.");
						} // else: all extended pathes are marked with bubble by clone()
						
						// kill partitions in super pathes
						for (int j = k-1; j>= 0; --j) {
							Vector<Path> buck2= (Vector<Path>) o[j];
							if (buck2.size()> 1) {
								PartitionSet
								if ().mergeBuckets(buck2)) {
									// remove bucket with all pathes
									openPartitions.remove(buck2.elementAt(0).getSource());
									
								}
							}
							System.currentTimeMillis();
						}
						
					}
				}
			}
			
			
			Set<Edge> outEdges= nodes[i].getOutEdges();
			Iterator<Edge> iterOE= outEdges.iterator();
			PartitionSet part= null;
			if (outEdges.size()>= n) { 	// partition for here
				part= new PartitionSet(outEdges.size());
				openPartitions.put(nodes[i], part);
			}
			HashMap<Path, Vector<Path>> extendedPathMap= new HashMap<Path, Vector<Path>>();
			while (iterOE.hasNext())  {
				Edge e= iterOE.next();
				Collection<Vector<Path>> buckets= nodes[i].getBuckets();
				Object[] o= buckets.toArray();
				Arrays.sort(o,vectorPathSourceCoordinateCompi);
				HashMap<Path, Path> pathExtendMap= new HashMap<Path, Path>();
				Vector<Path> succExtendedPathes= new Vector<Path>();	// per outEdge
				for (int j = 0; j < o.length; j++) {
					Vector<Path> pathes= (Vector<Path>) o[j];
					Iterator<Path> iterP= pathes.iterator();
					while (iterP.hasNext()) {
						Path p= iterP.next();
						long[] newTset= intersect(p.getTranscripts(), e.getHead().getTranscripts());
						if (!isNull(newTset)) {
							Path newP= (Path) p.clonePath();
							if (e.getPath()!= null)
								newP.addEdge(e);
							newP.addNode(e.getHead());
							newP.setTranscripts(newTset);
							
							pathExtendMap.put(p, newP);
							succExtendedPathes.add(newP);
							Vector<Path> vp= extendedPathMap.get(p);
							if (vp== null)
								extendedPathMap.put(p, new Vector<Path>(outEdges.size()));
							extendedPathMap.get(p).add(newP);
							
							e.getHead().addPath(newP);
						} 
					}
				}
				// new pathes
				if (outEdges.size()>= n) {
					Path newP= new Path(nodes[i]);
					if (e.getPath()!= null)
						newP.addEdge(e);
					newP.addNode(e.getHead());
					newP.setTranscripts(intersect(newP.getTranscripts(), e.getHead().getTranscripts()));
					part.addPath(newP);
					e.getHead().addToBucket(nodes[i], newP);
				}
			}
			// partition update
			Iterator<Path> iter= extendedPathMap.keySet().iterator();
			while (iter.hasNext()) {
				Path p= iter.next();				
				openPartitions.get(p.getSource()).splitPath(p, extendedPathMap.get(p));	// can be null?
			}
			// save mem
			nodes[i].bucketHash= null;
		}
		
		
	}

	public void extractASevents_old(int n, HashMap<String, Integer> map) {
		
		Comparator vectorPathSourceCoordinateCompi= new Comparator() {
			@Override
			public int compare(Object o1, Object o2) {
				Vector<Path> v1= (Vector<Path>) o1;
				Vector<Path> v2= (Vector<Path>) o2;
				return (defaultNodeByPositionTypeComparator.compare(
						v1.elementAt(0).getSource(), 
						v2.elementAt(0).getSource()));
			}
		};
		
		HashMap<Node, HashMap<Path,Path>> openBuckets= new HashMap<Node, HashMap<Path,Path>>();
		
		Node[] nodes= getNodesInGenomicOrder();
		for (int i = 0; i < nodes.length; i++) {
			Set<Edge> inEdges= nodes[i].getInEdges();
			if (inEdges.size()>= n) {
				Object[] o= nodes[i].getBuckets().toArray();
				Arrays.sort(o, vectorPathSourceCoordinateCompi);
				for (int k = o.length- 1; k >= 0; --k) {
					Vector<Path> bucket= (Vector<Path>) o[k];
					if (bucket.size()>= n) {
						boolean event= false;
						for (int j = 0; j < bucket.size(); j++) {
							if (bucket.elementAt(j).hasSpliceSite()) 
								event= true;
						}
						
						if (event) {
							Vector<Path[]> tuples= generateTuples(0,0,new int[n],bucket,new Vector<Path[]>());
							if (tuples.size()!= 0) {
								Bubble bubble= new Bubble(this, bucket, tuples);
								bubble.extractEvents(map);
							}
						}								
						// does not make sense to extend
						if (bucket.size()== bucket.elementAt(0).getSource().getOutEdges().size()) { 
							if (!nodes[i].removeBucket(bucket))
								System.err.println("bucket not found.");
						} // else: all extended pathes are marked with bubble by clone()
					}
				}
			}
			
			
			Set<Edge> outEdges= nodes[i].getOutEdges();
			Iterator<Edge> iterOE= outEdges.iterator();
			while (iterOE.hasNext())  {
				Edge e= iterOE.next();
				Collection<Vector<Path>> buckets= nodes[i].getBuckets();
				Object[] o= buckets.toArray();
				Arrays.sort(o,vectorPathSourceCoordinateCompi);
				HashMap<Path, Path> pathExtendMap= new HashMap<Path, Path>();
				Vector<Path> succExtendedPathes= new Vector<Path>();
				for (int j = 0; j < o.length; j++) {
					Vector<Path> pathes= (Vector<Path>) o[j];
					Iterator<Path> iterP= pathes.iterator();
					while (iterP.hasNext()) {
						Path p= iterP.next();
						long[] newTset= intersect(p.getTranscripts(), e.getHead().getTranscripts());
						if (!isNull(newTset)) {
							Path newP= (Path) p.clonePath();
							if (e.getPath()!= null)
								newP.addEdge(e);
							newP.addNode(e.getHead());
							newP.setTranscripts(newTset);
							
								// check for sub-/super-pathes that have successfully been extended over the same edge
								// iteration over outEdges from bucket with longest distance first
							Iterator<Path> iterSP= p.getSuperPathes().iterator();
							while (iterSP.hasNext()) {
								Path newSuper= pathExtendMap.get(iterSP.next());
								if (newSuper!= null) 
									newP.addSuperPath(newSuper);								
							}
							pathExtendMap.put(p, newP);
							succExtendedPathes.add(newP);
							
							HashMap<Path, Path> bucket= openBuckets.get(p.getSource());
							bucket.remove(p);
							bucket.put(newP, newP);
							e.getHead().addPath(newP);
						} 
					}
				}
				// new pathes
				if (outEdges.size()>= n) {
					Path newP= new Path(nodes[i]);
					if (e.getPath()!= null)
						newP.addEdge(e);
					newP.addNode(e.getHead());
					newP.setTranscripts(intersect(newP.getTranscripts(), e.getHead().getTranscripts()));
					e.getHead().addToBucket(nodes[i], newP);
					
					// is subpath of all sucessfully extended pathes
					for (int j = 0; j < succExtendedPathes.size(); j++) 
						newP.addSuperPath(succExtendedPathes.elementAt(j));
					
					HashMap<Path, Path> buck= openBuckets.get(newP.getSource());
					if (buck== null)
						openBuckets.put(newP.getSource(), new HashMap<Path, Path>());
					openBuckets.get(newP.getSource()).put(newP,newP);
				}
			}
		}
		
	}

	public Transcript[] decodeTset(long[] c) {
		int count= 0;
		Vector<Transcript> v= new Vector<Transcript>();
		for (int i = 0; i < c.length; i++) {
			int base= 0;
			long val= 1;
			int sum= 0;
			while(base< 64) {
				
				if ((c[i]&val)!= 0l) {
					v.add(trpts[count+ base]);
					sum+= val;
					if (sum== c[i])
						break;
				}
				++base;
				val*= 2;
			}
			count+= 64;
		}
		Transcript[] t= new Transcript[v.size()];
		for (int i = 0; i < t.length; i++) 
			t[i]= v.elementAt(i);
		return t;
	}
}
