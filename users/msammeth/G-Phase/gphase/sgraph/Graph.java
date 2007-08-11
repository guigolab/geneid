package gphase.sgraph;

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

import sun.reflect.ReflectionFactory.GetReflectionFactoryAction;

import gphase.io.gtf.GTFChrReader;
import gphase.model.ASEvent;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.DoubleVector;

public class Graph {
	public static void main(String[] args) {
		_070808_test();
	}
	
	static void _070808_test() {
		// 
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC0703.gtf
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_CDS.gtf		
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
		Gene[] g= reader.getGenes();
		int cnt= 0;
		HashMap<String, Integer> evMap= new HashMap<String, Integer>();
		long cumulGC= 0l, cumulEV= 0l;
		while (g!= null) {
			
			for (int i = 0; i < g.length; i++) {
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
				gr.extractASevents(2, evMap);
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
			
			try {
				reader.read();
			} catch (Exception e1) {
				e1.printStackTrace();
			}
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
			while(p> 64) {
				p-= 64;
				++cnt;
			}
			taVector[cnt]|= gphase.tools.Math.pow(2,p);
		}
		return taVector;
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
	
	public Edge createEdge(Node v, Node w) {
		Edge e= getEdge(v,w);
		if (e== null) {
			e= new Edge(v,w);
			edgeHash.put(new Integer(e.getTail().getSite().getPos()+e.getHead().getSite().getPos()),e);
		}
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
				createEdge(v,w);
			}
		}
		
		SpliceSite ss= new SpliceSite(Integer.MIN_VALUE, SpliceSite.TYPE_NOT_INITED);
		ss.setTranscripts(trpts);
		root= createNode(ss);
		Object[] o= rootHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) 
			createEdge(root, (Node) o[i]);
		
		ss= new SpliceSite(Integer.MAX_VALUE, SpliceSite.TYPE_NOT_INITED);
		ss.setTranscripts(trpts);
		leaf= createNode(ss);
		o= leafHash.keySet().toArray();
		for (int i = 0; i < o.length; i++) 
			createEdge((Node) o[i], leaf);
	}
	
	//DFS
	void contractGraph(Node n) {
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
			for (int i = 1; i < ns.size()- 1; i++) 
				removeNode((Node) ns.elementAt(i));
			
			for (int i = 0; i < edgeV.size(); i++) 
				removeEdge(edgeV.elementAt(i));
			
			Edge e= createEdge(p.getSource(), p.getSink());
			p.removeSource();
			p.removeSink();
			e.setPath(p);
		}
		
		if (outV.size()> 1) {
			Iterator<Edge> iter= outV.iterator();
			while (iter.hasNext())
				contractGraph(iter.next().getHead());
		} else 
			contractGraph(n);
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
		Iterator<Node> iter= nodeHash.values().iterator();
		int cnt= 0;
		Node[] nodes= new Node[nodeHash.size()];
		while (iter.hasNext()) 
			nodes[cnt++]= iter.next();
		Arrays.sort(nodes,defaultNodeByPositionTypeComparator);
		return nodes;
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
		
		HashMap<Node, HashMap<Path,Path>> openBuckets= new HashMap<Node, HashMap<Path,Path>>();
		
		Node[] nodes= getNodesInGenomicOrder();
		for (int i = 0; i < nodes.length; i++) {
			Set<Edge> inEdges= nodes[i].getInEdges();
			if (inEdges.size()>= n) {
				Iterator<Vector<Path>> iterP= nodes[i].getBuckets().iterator();
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
}
