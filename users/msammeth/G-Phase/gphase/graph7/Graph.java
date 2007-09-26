package gphase.graph7;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.AbstractList;
import java.util.AbstractQueue;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ConcurrentLinkedQueue;

import qalign.tools.CedricConqueror;

import com.sun.org.apache.bcel.internal.Constants;
import com.sun.org.apache.xml.internal.utils.NodeVector;

import sun.reflect.ReflectionFactory.GetReflectionFactoryAction;

import gphase.Toolbox;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASEvent;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.sgraph.Bubble;
import gphase.tools.DoubleVector;
import gphase.tools.File;
import gphase.tools.Time;

public class Graph {
	public static int counter= 0;
	
	static WriterThread writerThread= new WriterThread();
	static {
		writerThread.start();
	}
	
	protected static class WriterThread extends Thread {
		Queue<String> queue= new ConcurrentLinkedQueue<String>();
		Thread writingThread;
		boolean kill= false;
		
		public void addEvent(String evString) {
			queue.add(evString);
			writingThread.interrupt();
		}
		
		public void setKill(boolean newKill) {
			kill= newKill;
		}
		
		public void run() {
			writingThread= Thread.currentThread();
			while (true) {
				try {
					try {
						BufferedWriter buffy= new BufferedWriter(new FileWriter(outputFname, true));
						while (!queue.isEmpty()) {
							buffy.write(queue.poll()+"\n");
						}
						buffy.flush(); buffy.close();
						
					} catch (Exception e) {
						e.printStackTrace();
					}
					
					if (kill)
						return;
					
					try {
						writingThread.sleep(0);
					} catch (InterruptedException e) {
						; // :)
					}
				} catch (OutOfMemoryError e) {
					; // :)
				}
			} 
		}
	}
	
	
	public static String outputFname= "events.asta";
	
	public static void main(String[] args) { 
		_070808_test();
	}
	
	static void _070808_test() {
		boolean rusc= true;
		if (rusc)
			gphase.Constants.DATA_DIR= "/home/ugrusc/msammeth";
		
		// 
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf		
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
		gphase.tools.File file= new gphase.tools.File("/home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf");
		outputFname= "delme.asta";
		if (rusc)
			outputFname= "delme_rusc.asta";
		int n= 2;
		if (new File(outputFname).exists()) {
			new File(outputFname).delete();
			//System.out.println("File "+outputFname+" exists, check - I give up now.");
			//System.exit(0);
		}
		Species species= new Species("human");
		species.setGenomeVersion("hg18");
		boolean output= false, output2= true; 
		
		if (output2) {
			Date ti= new Date(System.currentTimeMillis());
			System.out.println("["+ti+"]  started. ");
		}
		long t0= System.currentTimeMillis();
		GTFChrReader reader= new GTFChrReader(file.getAbsolutePath()); 
		try {
			reader.read();
		} catch (Exception e1) { 
			e1.printStackTrace();
		}
		Gene[] g= reader.getGenes();
		if (output2) {
			Date ti= new Date(System.currentTimeMillis());
			System.out.println("["+ti+"] read: "+g[0].getChromosome()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
		}
		
		int cnt= 0;
		String chromo= null;
		while (g!= null) {
			long cumulGC= 0l, cumulEV= 0l;
			if (chromo== null)
				chromo= g[0].getChromosome();
			
			for (int i = 0; i < g.length; i++) {
				
				if (g[i].getTranscriptCount()== 1) {
					g[i]= null;
					continue;
				}
				g[i].setSpecies(species);

				String tmp= null;
				try {
					long t1= System.currentTimeMillis();
					Graph gr= new Graph(g[i]);
					if (output) {
						System.out.print(">  t= "+g[i].getTranscriptCount());
					}
					gr.init(n);
					long t2= System.currentTimeMillis();
					long dA= (t2-t1);
					cumulGC+= dA;
					if (output) {
						tmp= ", n= "+gr.nodeHash.size()+", e="+gr.edgeHash.size()+". aufbau "+(dA/1000)+" sec, ";
						System.out.print(tmp);
						System.out.flush();
					}
					//int oldSize= evMap.size();
					gr.getEventsByPathes(n);
					//cnt+= evMap.size()- oldSize;
					long dB= (System.currentTimeMillis()- t2);
					cumulEV+= dB;
					if (output) 
						//System.out.println("events "+(evMap.size()-oldSize));
						System.out.println(", extract "+(dB/1000)+" sec.");
				} catch (OutOfMemoryError e) {
					Time ti= new Time(System.currentTimeMillis());
					System.out.println("["+ti+"]");
					e.printStackTrace();
					System.out.println(g[i].getTranscripts()[0].getTranscriptID()+" "+g[i].getTranscriptCount()+" "+g[i].getSpliceSites().length);
					g[i]= null;
					if (tmp!= null)
						System.out.println(tmp);
				}
				
				//g[i]= null;
//				System.gc();
//				Graph g2= new Graph(g[i]);
//				g2.constructGraph();
//				Vector ev2= g2.extractASevents(2);			
//				System.currentTimeMillis();
			}
			if (output2) {
				Date ti= new Date(System.currentTimeMillis());
				System.out.println("["+ti+"] graph construct "+(cumulGC/1000)+" sec, extract events "+(cumulEV/ 1000)+" sec.");
			}
			System.gc();
			long tx= System.currentTimeMillis();
			writerThread.interrupt();
			try {
				reader.read();
			} catch (Exception e1) {
				e1.printStackTrace();
			}
			if (output2) {
				Date ti= new Date(System.currentTimeMillis());
				System.out.println("["+ti+"] read: "+chromo+" "+((System.currentTimeMillis()- tx)/ 1000)+" sec.");
			}
			chromo= null;
			g= reader.getGenes();
			System.gc();
			//Thread.currentThread().yield();
		}
		
		System.out.println("found "+cnt+" events.");
		System.out.println("took "+((System.currentTimeMillis()- t0)/1000)+" sec.");
		try {
			writerThread.setKill(true);
			writerThread.interrupt();
			writerThread.join();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
		}
		System.out.println("found "+counter+" events.");
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
							idx[i]= ++c;
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

	static TranscriptByNameComparator defaultTranscriptByNameComparator= new TranscriptByNameComparator();
	static Node.PositionTypeComparator defaultNodeByPositionTypeComparator= new Node.PositionTypeComparator();
	
	Gene gene;
	Transcript[] trpts;
	int taSize;
	HashMap<SpliceSite, Node> nodeHash= new HashMap<SpliceSite, Node>();
	HashMap<String, Edge> edgeHash= new HashMap<String, Edge>();
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
		return edgeHash.get(v.getSite().toString()+w.getSite().toString());
	}
	
	public Node createNode(SpliceSite ss) {
		Node n= getNode(ss);
		if (n== null) {
			if (ss.isSpliceSite()&& (!ss.isCanonical()))
				return null;
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
			edgeHash.put(v.getSite().toString()+w.getSite().toString(),e);
		} else 
			e.setTranscripts(unite(e.getTranscripts(),newTset));
		return e;
	}

	public Edge addEdge(Node v, Node w, long[] newTset) {
		Edge e= new Edge(v,w);
		e.setTranscripts(newTset);
		edgeHash.put(v.getSite().toString()+w.getSite().toString(),e);
		return e;
	}
	
	public void init(int n) {
		constructGraph();
		//cleanGraphByNodes();
		contractGraph(n);
	}
	
	void constructGraph() {
		Transcript[] t= gene.getTranscripts();
		HashMap<Node, Node> rootHash= new HashMap<Node, Node>(), leafHash= new HashMap<Node, Node>();
		for (int i = 0; i < t.length; i++) {
			SpliceSite[] ss= t[i].getSpliceSitesAll();
			for (int j = 1; j < ss.length; j++) {
				Node v= createNode(ss[j-1]);
				if (v== null)
					continue;
				if (j== 1&& rootHash.get(v)== null)
					rootHash.put(v,v);
				Node w= createNode(ss[j]);
				if (w== null)
					continue;
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
	
	void contractGraph(int k) {
		Vector<Edge> outV= root.getOutEdges();
		for (int i = 0; i < outV.size(); i++) {
			contractGraph(k, outV.elementAt(i));
		}
		
			// purge hashes
		Object[] keys= edgeHash.keySet().toArray();
		Object key;
		for (int i = 0; i < keys.length; i++) {
			Edge e= edgeHash.get(keys[i]);
			if (e.isContracted())
				continue;
			edgeHash.remove(keys[i]);
			e.getTail().removeOutEdge(e);
			e.getHead().removeInEdge(e);
		}
		
		keys= nodeHash.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			Node v= nodeHash.get(keys[i]);
			if (v.getInEdges().size()> 0|| v.getOutEdges().size()> 0)
				continue;
			nodeHash.remove(keys[i]);
		}
		
	}
	
	//DFS
	void contractGraph(int k, Edge e) {
		Node n= e.getHead();
		if (e.isProcessed())
			return;
		e.setProcessed(true);
		
		Vector<Edge> outV= n.getOutEdges();
		Vector<Edge> inV= n.getInEdges();
		Edge srcEdge= e, endEdge= null;
		long[] newTset= null;
		while (outV.size()< k&& inV.size()< k&& outV.size()> 0) {
			
			endEdge= outV.iterator().next();
			n= endEdge.getHead();
			outV= n.getOutEdges();
			inV= n.getInEdges();			
		}
		
			// contract
		if (endEdge!= null) {

			Edge f= new Edge(srcEdge.getTail(),endEdge.getHead());
			f.setTranscripts(srcEdge.getTranscripts());
			f.setContracted(true);
			edgeHash.put(f.getTail().getSite().toString()+f.getHead().getSite().toString(),f);	// cannot be in there
			
		} else {
			e.setContracted(true);
		}

		Object[] o= outV.toArray();
		for (int i = 0; i < o.length; i++) {
			contractGraph(k, (Edge) o[i]);
		}
	}

	//DFS
	void contractGraph_old_bugs(Node n) {
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
		}

		
		Object[] o= outV.toArray();
		for (int i = 0; i < o.length; i++) {
			contractGraph_old_bugs(((Edge) o[i]).getHead());
		}
	}
	
	public void getEventsByPathes(int n) {		
			Node[] nodes= getNodesInGenomicOrder();
			for (int i = 0; i < nodes.length; i++) {
				Vector<Edge> inEdges= nodes[i].getInEdges();
				Vector<Edge> outEdges= nodes[i].getOutEdges();
				if (inEdges.size()>= n) {
					
					Object[] fNodes= nodes[i].getFromNodes().toArray();
					Arrays.sort(fNodes,defaultNodeByPositionTypeComparator);
					for (int j = fNodes.length- 1; j >= 0; --j) {
						Node fromNode= (Node) fNodes[j];
						Vector<Path> pathes= nodes[i].getFromNodeMap().get(fNodes[j]);
						Vector<Vector<Path>> splitPathes= ((Node) fNodes[j]).splitPathes(pathes);
						// bubble
						if (splitPathes.size()>= n) {
	
								// retrieve event
							Vector<Path[]> tuples=  new Vector<Path[]>();	// !!!! killed	generateTuples(n, splitPathes); 
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
									outputEvent(ev);
	//								Integer nr= evMap.get(ev.toString());
	//								if (nr== null)
	//									nr= new Integer(0);
	//								evMap.put(ev.toString(),new Integer(nr.intValue()+1));
	//								try {
	//									FileWriter writer= new FileWriter("new.asta", true);
	//									writer.write(ev.toStringASTA()+"\n");
	//									writer.flush(); writer.close();
	//								} catch (Exception e) {
	//									e.printStackTrace();
	//								}
								}
							}	// for all tuples
							
							// merge partitions
							//nodes[i].getFromNodeMap().remove(fNodes[j]); !!! Nooooo
							((Node) fNodes[j]).mergePartitions(splitPathes);
							
						} // if bubble
					} // nodesWithOutPSS
					
				}
	
				
				Object[] o= outEdges.toArray();
				// add possible new source
				if (outEdges.size()>= n) {
						// generate new pathes
					Vector<Vector<Path>> partitions= new Vector<Vector<Path>>(outEdges.size());
					for (int j = 0; j < o.length; j++) {
						Edge e= (Edge) o[j];
						Path newP= new Path();
						newP.setSourceEdge(e);
						newP.setSinkEdge(e);
						newP.setTranscripts(e.getTranscripts());
						
						Vector<Path> v= e.getHead().getFromNodeMap().get(nodes[i]);
						if (v== null) {
							v= new Vector<Path>(2);
							e.getHead().getFromNodeMap().put(nodes[i], v);
						}
						v.add(newP);
					}
				}
								
				Iterator<Node> iterSrc= nodes[i].getFromNodes().iterator();
				while (iterSrc.hasNext()) {
					Node src= iterSrc.next();
					Vector<Path> pathes= nodes[i].getFromNodeMap().get(src);
					// extend all pathes
					for (int k = 0; k < pathes.size(); k++) {
						Path p= pathes.elementAt(k);
						if (p.getSourceNode().getOutPartitionSize()< n)
							continue;	// do not extend, cannot generate event
						for (int j = 0; j < o.length; j++) {
							Edge e= (Edge) o[j];
							long[] newTsup= intersect(p.getTranscripts(), e.getTranscripts());
							if (isNull(newTsup))
								continue;
							Path newP= new Path();
							newP.setSourceEdge(p.getSourceEdge());
							newP.setTranscripts(newTsup);
							newP.setSinkEdge(e);
							
							Vector<Path> v= e.getHead().getFromNodeMap().get(src);
							if (v== null) {
								v= new Vector<Path>(2);
								e.getHead().getFromNodeMap().put(src, v);
							}
							v.add(newP);
						}							
					}
				}	// for all src
				
				nodes[i].fromNodeMap= null;
				
			} // for all nodes
		}

	public void getEventsByPartitions(int n) {		
			Node[] nodes= getNodesInGenomicOrder();
			Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
			Vector<long[]> currentPartitionSet= new Vector<long[]>();
			for (int i = 0; i < nodes.length; i++) {
				
					// check bubble
				Vector<Edge> inEdges= nodes[i].getInEdges();
				if (inEdges.size()>= n) {
					Iterator<Edge> iterEdge= inEdges.iterator();
					while (iterEdge.hasNext()) {
						long[] unity= null;
						Edge e= iterEdge.next();
						
					}
				}				
				
				Vector<Edge> outEdges= nodes[i].getOutEdges();
				if (inEdges.size()>= n) {
					
					Object[] fNodes= nodes[i].getFromNodes().toArray();
					Arrays.sort(fNodes,defaultNodeByPositionTypeComparator);
					for (int j = fNodes.length- 1; j >= 0; --j) {
						Node fromNode= (Node) fNodes[j];
						Vector<Path> pathes= nodes[i].getFromNodeMap().get(fNodes[j]);
						Vector<Vector<Path>> splitPathes= ((Node) fNodes[j]).splitPathes(pathes);
						// bubble
						if (splitPathes.size()>= n) {
	
								// retrieve event
							Vector<Path[]> tuples=  new Vector<Path[]>();	// !!!! killed	generateTuples(n, splitPathes); 
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
									outputEvent(ev);
	//								Integer nr= evMap.get(ev.toString());
	//								if (nr== null)
	//									nr= new Integer(0);
	//								evMap.put(ev.toString(),new Integer(nr.intValue()+1));
	//								try {
	//									FileWriter writer= new FileWriter("new.asta", true);
	//									writer.write(ev.toStringASTA()+"\n");
	//									writer.flush(); writer.close();
	//								} catch (Exception e) {
	//									e.printStackTrace();
	//								}
								}
							}	// for all tuples
							
							// merge partitions
							//nodes[i].getFromNodeMap().remove(fNodes[j]); !!! Nooooo
							((Node) fNodes[j]).mergePartitions(splitPathes);
							
						} // if bubble
					} // openSrcVec
					
				}
	
				
				Object[] o= outEdges.toArray();
				// add possible new source
				if (outEdges.size()>= n) {
						// generate new pathes
					Vector<Vector<Path>> partitions= new Vector<Vector<Path>>(outEdges.size());
					for (int j = 0; j < o.length; j++) {
						Edge e= (Edge) o[j];
						Path newP= new Path();
						newP.setSourceEdge(e);
						newP.setSinkEdge(e);
						newP.setTranscripts(e.getTranscripts());
						
						Vector<Path> v= e.getHead().getFromNodeMap().get(nodes[i]);
						if (v== null) {
							v= new Vector<Path>(2);
							e.getHead().getFromNodeMap().put(nodes[i], v);
						}
						v.add(newP);
					}
				}
								
				Iterator<Node> iterSrc= nodes[i].getFromNodes().iterator();
				while (iterSrc.hasNext()) {
					Node src= iterSrc.next();
					Vector<Path> pathes= nodes[i].getFromNodeMap().get(src);
					// extend all pathes
					for (int k = 0; k < pathes.size(); k++) {
						Path p= pathes.elementAt(k);
						if (p.getSourceNode().getOutPartitionSize()< n)
							continue;	// do not extend, cannot generate event
						for (int j = 0; j < o.length; j++) {
							Edge e= (Edge) o[j];
							long[] newTsup= intersect(p.getTranscripts(), e.getTranscripts());
							if (isNull(newTsup))
								continue;
							Path newP= new Path();
							newP.setSourceEdge(p.getSourceEdge());
							newP.setTranscripts(newTsup);
							newP.setSinkEdge(e);
							
							Vector<Path> v= e.getHead().getFromNodeMap().get(src);
							if (v== null) {
								v= new Vector<Path>(2);
								e.getHead().getFromNodeMap().put(src, v);
							}
							v.add(newP);
						}							
					}
				}	// for all src
				
				nodes[i].fromNodeMap= null;
				
			} // for all nodes
		}

	public void getEventsByPathes_new(int n, HashMap<String, Integer> evMap) {		
		Node[] nodes= getNodesInGenomicOrder();
		Vector<Node> nodesWithOutPSS= new Vector<Node>(nodes.length/ 2);
		for (int i = 0; i < nodes.length; i++) {
			Vector<Edge> inEdges= nodes[i].getInEdges();
			Vector<Edge> outEdges= nodes[i].getOutEdges();
			Gene g= this.gene;
			if (inEdges.size()>= n) {
				
				Object[] fNodes= nodes[i].getFromNodes().toArray();
				Arrays.sort(fNodes,defaultNodeByPositionTypeComparator);
				for (int j = fNodes.length- 1; j >= 0; --j) {
					Node fromNode= (Node) fNodes[j];
					Vector<Path> pathes= nodes[i].getFromNodeMap().get(fNodes[j]);
					Vector<Vector<Path>> splitPathes= ((Node) fNodes[j]).splitPathes(pathes);
					// bubble
					if (splitPathes.size()>= n) {

							// retrieve event
						Vector<Path[]> tuples= generateTuples(n, splitPathes);
						for (int k = 0; k < tuples.size(); k++) {
							if (!((Node) fNodes[j]).checkTuple(tuples.elementAt(k)))
								continue;
							 
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
								outputEvent(ev);
							}
						}	// for all tuples
						
						// merge partitions
						//nodes[i].getFromNodeMap().remove(fNodes[j]); !!! Nooooo
						((Node) fNodes[j]).mergePartitions(splitPathes);
						
					} // if bubble
				} // nodesWithOutPSS
				
			}

			
			Object[] o= outEdges.toArray();
			// add possible new source
			if (outEdges.size()>= n) {
					// generate new pathes
				Vector<Vector<Path>> partitions= new Vector<Vector<Path>>(outEdges.size());
				for (int j = 0; j < o.length; j++) {
					Edge e= (Edge) o[j];
					Path newP= new Path();
					newP.setSourceEdge(e);
					newP.setSinkEdge(e);
					newP.setTranscripts(e.getTranscripts());
					
					Vector<Path> v= e.getHead().getFromNodeMap().get(nodes[i]);
					if (v== null) {
						v= new Vector<Path>(2);
						e.getHead().getFromNodeMap().put(nodes[i], v);
					}
					v.add(newP);
				}
			}
							
			Iterator<Node> iterSrc= nodes[i].getFromNodes().iterator();
			while (iterSrc.hasNext()) {
				Node src= iterSrc.next();
				Vector<Path> pathes= nodes[i].getFromNodeMap().get(src);
				// extend all pathes
				for (int k = 0; k < pathes.size(); k++) {
					Path p= pathes.elementAt(k);
					if (p.getSourceNode().getOutPartitionSize()< n)
						continue;	// do not extend, cannot generate event
					for (int j = 0; j < o.length; j++) {
						Edge e= (Edge) o[j];
						long[] newTsup= intersect(p.getTranscripts(), e.getTranscripts());
						if (isNull(newTsup))
							continue;
						Path newP= new Path();
						newP.setSourceEdge(p.getSourceEdge());
						newP.setTranscripts(newTsup);
						newP.setSinkEdge(e);
						
						Vector<Path> v= e.getHead().getFromNodeMap().get(src);
						if (v== null) {
							v= new Vector<Path>(2);
							e.getHead().getFromNodeMap().put(src, v);
						}
						v.add(newP);
					}							
				}
			}	// for all src
			
			nodes[i].setFromNodeMap(null);
			
		} // for all nodes
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
	
	private void outputEvent(final ASEvent event) {
		writerThread.addEvent(event.toStringASTA());
		++counter;
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

	void cleanGraphByNodes() {
		Object[] nodeO= nodeHash.keySet().toArray();
		for (int j = 0; j < nodeO.length; ++j) {
			Node v= (Node) nodeHash.get(nodeO[j]);
			if ((!v.getSite().isSpliceSite())|| v.getSite().isCanonical())
				continue;
				// delete
			Object[] edgeO= v.getInEdges().toArray();
			for (int i = 0; i < edgeO.length; i++) 
				removeEdge((Edge) edgeO[i]);
			edgeO= v.getOutEdges().toArray();
			for (int i = 0; i < edgeO.length; i++) 
				removeEdge((Edge) edgeO[i]);
			removeNode(v);
		}
	}
}
