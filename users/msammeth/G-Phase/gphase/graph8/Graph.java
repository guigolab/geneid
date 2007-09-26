package gphase.graph8;

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

/**
 * global partition, several versions.
 * last one has only exception in EST set, but skipping - time is not good
 * continuing with graph9
 * @author msammeth
 *
 */

public class Graph {
	public static int counter= 0;
	
	static WriterThread writerThread= new WriterThread();
	static {
		writerThread.start();
	}
	
		public static class EventExtractorThread extends Thread {
				
				Gene[] g;
				ChromosomeReaderThread upstreamThread;
				int n= 2;
				boolean output= false, output2= true; 
				Species species= new Species("human");
				
				public EventExtractorThread(ChromosomeReaderThread newUpstreamThread) {
					super();
					setName("event_extraction_thread");
					upstreamThread= newUpstreamThread;
					species.setGenomeVersion("hg18");
				}
				
				@Override
				public void run() {

					String chromo= null;
					long cumulGC= 0l, cumulGF= 0l, cumulGT= 0l, cumulEV= 0l;
					for (int i = 0; i < g.length; i++) {
						
						if (chromo== null)
							chromo= g[0].getChromosome();
												
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
							//gr.init(n);
							gr.constructGraph();
							long t2= System.currentTimeMillis();
							long dA= (t2-t1);
							cumulGC+= dA;
							if (output) {
								tmp= ", n= "+gr.nodeHash.size()+", e="+gr.edgeHash.size()+". aufbau "+(dA/1000)+" sec, ";
								System.out.print(tmp);
								System.out.flush();
							}
							
							
							gr.collapseFuzzyFlanks();
							long t22= System.currentTimeMillis();
							long dF= (t22-t2);
							cumulGF+= dF;
							if (output) {
								tmp+= ",\t n= "+gr.nodeHash.size()+", e="+gr.edgeHash.size()+". flanks "+(dF/1000)+" sec, ";
								System.out.print(tmp);
								System.out.flush();
							}
							
							gr.contractGraph(n);
							long t3= System.currentTimeMillis();
							long dT= (t3-t22);
							cumulGT+= dT;
							if (output) {
								tmp+= ",\t n= "+gr.nodeHash.size()+", e="+gr.edgeHash.size()+". extract "+(dT/1000)+" sec, ";
								System.out.print(tmp);
								System.out.flush();
							}
							//int oldSize= evMap.size();
							gr.getEventsByPartitions(n);
							//cnt+= evMap.size()- oldSize;
							long dB= (System.currentTimeMillis()- t3);
							cumulEV+= dB;
							if (output) 
								//System.out.println("events "+(evMap.size()-oldSize));
								System.out.println(", extract "+(dB/1000)+" sec.");
						} catch (OutOfMemoryError e) {
							Time ti= new Time(System.currentTimeMillis());
							System.out.println("["+ti+"]");
							e.printStackTrace();
							System.out.println(g[i].getTranscripts()[0].getTranscriptID()+" "+g[i].getTranscriptCount()+" "+g[i].getSpliceSites().length);
							if (tmp!= null)
								System.out.println(tmp);
						} finally {
							g[i]= null;
							writerThread.interrupt();
						}
						
						
					}	// for
					
					if (output2) {
						Date ti= new Date(System.currentTimeMillis());
						System.out.println("["+ti+"] "+chromo+
								" graph construct "+(cumulGC/1000)+
								" sec, fuzzy flanks "+(cumulGF/1000) +
								" sec, contraction "+(cumulGT/1000) +
								" sec, extract events "+(cumulEV/ 1000)+" sec.");
					}
					System.gc();
					writerThread.interrupt();

				}

				public boolean isOutput() {
					return output;
				}

				public void setOutput(boolean output) {
					this.output = output;
				}

				public boolean isOutput2() {
					return output2;
				}

				public void setOutput2(boolean output2) {
					this.output2 = output2;
				}

				public Gene[] getG() {
					return g;
				}

				public void setG(Gene[] g) {
					this.g = g;
				}
		}
		
		public static class ChromosomeReaderThread extends Thread {
		
		GTFChrReader reader;
		Gene[] g;
		EventExtractorThread downstreamThread;
		boolean output= false, output2= true; 
		
		public ChromosomeReaderThread(GTFChrReader newReader) {
			super();
			this.reader= newReader;
			setName("chr_reader_thread");
		}
		
		@Override
		public void run() {
			
			downstreamThread= new EventExtractorThread(this);
			downstreamThread.setOutput(output);
			downstreamThread.setOutput2(output2);
			while (true) {
				long t0= System.currentTimeMillis();
				try {
					reader.read();
				} catch (Exception e1) { 
					e1.printStackTrace();
				}
				if (reader.getUnclusteredGeneNb()== 0)
					break;
				if (output2) {
					Date ti= new Date(System.currentTimeMillis());
					System.out.println("["+ti+"] read: "+reader.getLastReadChr()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
				}
				t0= System.currentTimeMillis();
				g= reader.getGenes();
				if (g== null) {
					System.out.println(" => no genes, stopping.");	// just in case
				}
					
				if (output2) {
					Date ti= new Date(System.currentTimeMillis());
					System.out.println("["+ti+"] clustered: "+g[0].getChromosome()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
				}
				
				if (downstreamThread.isAlive())
					try {
						downstreamThread.join();
					} catch (InterruptedException e1) {
						; //:)
					}
					
				downstreamThread= new EventExtractorThread(this);
				downstreamThread.setOutput(output);
				downstreamThread.setOutput2(output2);
				downstreamThread.setG(g);
				downstreamThread.start();
				g= null;
			}
			
		}

		public EventExtractorThread getDownstreamThread() {
			return downstreamThread;
		}

		public void setDownstreamThread(EventExtractorThread downstreamThread) {
			this.downstreamThread = downstreamThread;
		}

		public Gene[] getG() {
			return g;
		}

		public void setG(Gene[] g) {
			this.g = g;
		}

		public GTFChrReader getReader() {
			return reader;
		}

		public void setReader(GTFChrReader reader) {
			this.reader = reader;
		}

		public boolean isOutput() {
			return output;
		}

		public void setOutput(boolean output) {
			this.output = output;
		}

		public boolean isOutput2() {
			return output2;
		}

		public void setOutput2(boolean output2) {
			this.output2 = output2;
		}
	}
	
	protected static class WriterThread extends Thread {
		Queue<String> queue= new ConcurrentLinkedQueue<String>();
		Thread writingThread;
		boolean kill= false;
		
		public WriterThread() {
			super();
			setName("event_writing_thread");
		}
		
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
							String s= queue.poll();
							buffy.write(s+"\n");
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
		//_070808_test();
		_240808_test_multithread(args);
	}
	
	static void _070808_test() {
			boolean rusc= false;
			if (rusc)
				gphase.Constants.DATA_DIR= "/home/ugrusc/msammeth";
			
			// 
			// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
			// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf		
			// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
			// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
			gphase.tools.File file= new gphase.tools.File("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf");
			outputFname= "delme.asta"; 
			if (rusc)
				outputFname= "delme_rusc_chr11.asta";
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
			//reader.sweepToChromosome("chr9");
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
			HashMap<String, Integer> evMap= null;//new HashMap<String, Integer>();
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
						gr.getEventsByPartitions(n);
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
				g= reader.getGenes();
				System.gc();
				if (g!= null&& g.length> 0)
					chromo= g[0].getChromosome();
				else
					chromo= null;
				if (output2&& chromo!= null) {
					Date ti= new Date(System.currentTimeMillis());
					System.out.println("["+ti+"] read: "+chromo+" "+((System.currentTimeMillis()- tx)/ 1000)+" sec.");
				}
				//Thread.currentThread().yield();
			}
			
			System.out.println("found "+cnt+" events.");
			if (evMap!= null) {
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
			}
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

	static void _240808_test_multithread(String[] args) {
		boolean rusc= false;
		if (args.length> 1&& args[1].equalsIgnoreCase("rusc"))
			rusc= true;
		if (rusc)
			gphase.Constants.DATA_DIR= "/home/ugrusc/msammeth";
		
		// 
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
		// /home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf		
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
		gphase.tools.File file= new gphase.tools.File("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf");
		outputFname= args[0];	// "/home/msammeth/graph8_EST_start_chr6.asta";
		boolean output= false, output2= true; 
//		if (rusc)
//			outputFname= "delme.asta";
		if (new File(outputFname).exists()) {
			new File(outputFname).delete();
			//System.out.println("File "+outputFname+" exists, check - I give up now.");
			//System.exit(0);
		}
		
			// init and start threads
		long t0= System.currentTimeMillis();
		if (output2) {
			Date ti= new Date(t0);
			System.out.println("["+ti+"]  started. ");
		}
		GTFChrReader reader= new GTFChrReader(file.getAbsolutePath()); 
		reader.sweepToChromosome("chr14");
		ChromosomeReaderThread readerThread= new ChromosomeReaderThread(reader);
		readerThread.setOutput(output);
		readerThread.setOutput2(output2);
		readerThread.start();
		try {
			readerThread.join();
			readerThread.getDownstreamThread().join();
		} catch (InterruptedException e1) {
			;	// :)
		}

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
	
	public static boolean intersects(long[] a, long[] b) {
		for (int i = 0; i < b.length; i++) {
			if ((a[i]& b[i])!= 0l)
				return true;
		}
		return false;
	}
	
	public static long[] unite(long[] sig1, long[] sig2) {
		assert(sig1.length== sig2.length);
		long[] reSig= new long[sig1.length];
		for (int i = 0; i < sig2.length; i++) 
			reSig[i]= sig1[i]| sig2[i];
		return reSig;
	}
	
	public static long[] xor(long[] a, long[] b) {
		long[] c= new long[a.length];
		for (int i = 0; i < c.length; i++) 
			c[i]= a[i]^ b[i];
		return c;
	}

	public static long[] without(long[] a, long[] b) {
		long[] inter= intersect(a, b);
		long[] c= new long[a.length];
		for (int i = 0; i < c.length; i++) {
			c[i]= a[i]^ inter[i];
		}
		return c;
	}
	
	
	public static long[] unite(Vector<long[]> a) {
		if (a.size()== 0)
			return null;
		if (a.size()== 1)
			return a.elementAt(0);
		long[] inter= unite(a.elementAt(0), a.elementAt(1));
		for (int i = 2; i < a.size(); i++) 
			inter= unite(inter, a.elementAt(i));
		return inter;
	}
	
	public static long[] createNullArray(int size) {
		long[] a= new long[size];
		for (int i = 0; i < a.length; i++) 
			a[i]= 0l;
		return a;
	}

	public long[] createAllArray() {
		long[] a= encodeTset(trpts);
		return a;
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

	void generateTuples(int k, Node snk, Vector<Vector<long[]>> buckets, Vector<Node> openNodes, HashMap<Transcript, Node> validSinceMap) {
	
		int[] idx= new int[k];
		for (int i = 0; i < idx.length; i++) 
			idx[i]= i;
		
		int[] combi= new int[k];
		while (idx[0]< (buckets.size()-k+ 1)) {
			
				// now all combinations
			for (int j = 0; j < combi.length; j++) 
				combi[j]= 0;
			while (true) {
				
				long[] combination= null;
				for (int j = 0; j < combi.length; j++) 
					if (combination== null)
						combination= buckets.elementAt(idx[j]).elementAt(combi[j]);
					else
						combination= unite(combination, buckets.elementAt(idx[j]).elementAt(combi[j]));
				
				Node stopNode= null;
				Transcript[] t= decodeTset(combination);
				for (int j = 0; j < t.length; j++) {
					if (stopNode== null) 
						stopNode= validSinceMap.get(t[j]);
					else {
						Node n= validSinceMap.get(t[j]);
						if (n!= null&& n.getSite().getPos()> stopNode.getSite().getPos())
							stopNode= n;
					}
				}
				
					// look for src node
				int x = openNodes.size()- 1;
				Node src= null;
				for (; x >= 0; --x) {
					if (stopNode!= null&& openNodes.elementAt(x).getSite().getPos()< stopNode.getSite().getPos())
						break;
					long[] inter= intersect(combination, openNodes.elementAt(x).getTranscripts());
					if (equalSet(inter, openNodes.elementAt(x).getTranscripts()))
						src= openNodes.remove(x);	// could additionally merge partitions and check vs < k
					else
						src= openNodes.elementAt(x);
					if (equalSet(inter, combination))
						break;
				}
				
					// extract event
				if (src!= null) {	// can happen for invalid introns kill a path from the common src
					SpliceSite[][] ss= new SpliceSite[k][];
					Transcript[][] tt= new Transcript[k][];
					for (int h = 0; h < combi.length; h++) {
						tt[h]= decodeTset(buckets.elementAt(idx[h]).elementAt(combi[h]));
						ss[h]= tt[h][0].getSpliceSitesBetween(src.getSite(), snk.getSite());
					}
					ASEvent ev= new ASEvent(tt,ss);
					boolean event= true;
					if (isRoot(src)|| isLeaf(snk))
						event= ev.isASevent();
					if (event) 
						outputEvent(ev);
				}
			
					// next combination on same partitions
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
			
				// next combination between partitions
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
	
	}

	void generateTuples_without_validity(int k, Node snk, Vector<Vector<long[]>> buckets, Vector<Node> openNodes) {
	
			int[] idx= new int[k];
			for (int i = 0; i < idx.length; i++) 
				idx[i]= i;
			
			int[] combi= new int[k];
			while (idx[0]< (buckets.size()-k+ 1)) {
				
					// now all combinations
				for (int j = 0; j < combi.length; j++) 
					combi[j]= 0;
				while (true) {
					
					long[] combination= null;
					for (int j = 0; j < combi.length; j++) 
						if (combination== null)
							combination= buckets.elementAt(idx[j]).elementAt(combi[j]);
						else
							combination= unite(combination, buckets.elementAt(idx[j]).elementAt(combi[j]));
					
						// look for src node
					int x = openNodes.size()- 1;
					Node src= null;
					for (; x >= 0; --x) {
						long[] inter= intersect(combination, openNodes.elementAt(x).getTranscripts());
						if (equalSet(inter, openNodes.elementAt(x).getTranscripts()))
							src= openNodes.remove(x);	// could additionally merge partitions and check vs < k
						else
							src= openNodes.elementAt(x);
						if (equalSet(inter, combination))
							break;
					}
					
						// extract event
					if (src!= null) {	// can happen for invalid introns kill a path from the common src
						SpliceSite[][] ss= new SpliceSite[k][];
						Transcript[][] tt= new Transcript[k][];
						for (int h = 0; h < combi.length; h++) {
							tt[h]= decodeTset(buckets.elementAt(idx[h]).elementAt(combi[h]));
							ss[h]= tt[h][0].getSpliceSitesBetween(src.getSite(), snk.getSite());
						}
						ASEvent ev= new ASEvent(tt,ss);
						boolean event= true;
						if (isRoot(src)|| isLeaf(snk))
							event= ev.isASevent();
						if (event) 
							outputEvent(ev);
					}
				
						// next combination on same partitions
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
				
					// next combination between partitions
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
	boolean canonicalSS= false;
	
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
			if (canonicalSS&& ss.isSpliceSite()&& (!ss.isCanonical()))
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
		collapseFuzzyFlanks();
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
		
		SpliceSite ss= new SpliceSite(Integer.MIN_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
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
		
		ss= new SpliceSite(Integer.MAX_VALUE, SpliceSite.TYPE_NOT_INITED, gene);
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
			long[] temp= encodeTset(tt);
			createEdge((Node) o[i], leaf, encodeTset(tt));
		}
	}
	
	void collapseFuzzyFlanks(boolean forRoot) {
		Iterator<Edge> iterEdge= null, iterEdge2;
		if (forRoot)
			iterEdge= root.getOutEdges().iterator();
		else
			iterEdge= leaf.getInEdges().iterator();
		Edge e,f;
		HashMap<Node, Vector<Edge>> map= new HashMap<Node, Vector<Edge>>();
		Vector<Edge> v;
		Node n, u;
		while (iterEdge.hasNext()) {
			e= iterEdge.next();
			if (forRoot)
				n= e.getHead();
			else 
				n= e.getTail();
			if (n.getInEdges().size()== 1&& n.getOutEdges().size()== 1) {	// TODO <k ??
				if (forRoot)
					iterEdge2= n.getOutEdges().iterator();
				else
					iterEdge2= n.getInEdges().iterator();
				while (iterEdge2.hasNext()) {
					f= iterEdge2.next();
					if (forRoot)
						u= f.getHead();
					else
						u= f.getTail();
					v= map.get(u);
					if (v== null) {
						v= new Vector<Edge>();
						map.put(u, v);
					}
					// only remove second one first, see later whether first edge still needed
					v.add(f);
				}
			}
		}
		
			// remove second edges
		Iterator<Node> iterNode= map.keySet().iterator();
		while (iterNode.hasNext()) {
			n= iterNode.next();
			v= map.get(n);
			long[] unity= null;
			if (v.size()> 1) {	// more than one chain of 2 edges
				for (int i = 0; i < v.size(); i++) {
					e= v.elementAt(i);
					if ((forRoot&& (!isRoot(e.getTail())))|| ((!forRoot)&& (!isLeaf(e.getHead())))) {		// .. because of deleted 2nd edges
						if (unity== null)
							unity= e.getTranscripts();
						else
							unity= unite(unity, e.getTranscripts());
					}
					removeEdge(e);	// nodes will be purged of the end of collapse
				}
				if (forRoot)
					createEdge(root, n, unity);
				else
					createEdge(n, leaf, unity);
			}
			
		}
		
			// check first edges
		Object[] o= null;
		if (forRoot)
			o= root.getOutEdges().toArray();
		else
			o= leaf.getInEdges().toArray();
		for (int i = 0; i < o.length; i++) {
			e= (Edge) o[i];
			if ((forRoot&& e.getHead().getOutEdges().size()== 0)|| 
					((!forRoot)&& e.getTail().getInEdges().size()== 0))	// no longer needed
				removeEdge(e);
		}		
	}
	
	void collapseFuzzyFlanks() {
		// kill divergences only by start/end
		collapseFuzzyFlanks(true);
		collapseFuzzyFlanks(false);
	}
	
	void contractGraph(int k) {
		
			// contract the rest of the nodes
		Node[] nodes= getNodesInGenomicOrder();	// have to do all for "full breaks" by edge filtering (chr2:179,380,178-179,380,190)
		for (int i = 0; i < nodes.length; i++) {
			if (nodes[i].getInEdges().size()> 0)
				continue;
			Vector<Edge> outV= nodes[i].getOutEdges();
			for (int j = 0; j < outV.size(); ++j) {
				if (!outV.elementAt(j).isContracted())
					contractGraph(k, outV.elementAt(j));
			}
		}
		nodesInGenomicOrder= null;
		
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
		Edge srcEdge= e, endEdge= null, f= null;
		long[] newTset= e.getTranscripts();	// has to be since edges are deleted
		while (outV.size()< k&& inV.size()< k&& outV.size()> 0) {
			f= outV.iterator().next();
			newTset= intersect(newTset, f.getTranscripts());
			if (isNull(newTset))
				break;	// because of deleted edges
			endEdge= f;
			n= endEdge.getHead();
			outV= n.getOutEdges();
			inV= n.getInEdges();			
		}
		
			// contract
		if (endEdge!= null) {

			f= new Edge(srcEdge.getTail(),endEdge.getHead());
			f.setTranscripts(newTset);	// srcEdge.getTranscripts() not, since edges are deleted; really TODO
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

	// search for the partitions that are 
	Partition splitPartitions(Vector<Edge> inEdges, Vector<Partition> currentPartitionSet, Vector<Vector<long[]>> splitPartitionSet) {
		Iterator<Edge> iterEdge= inEdges.iterator();
		Vector<long[]> unitedPartition= new Vector<long[]>();
		long[] unity= null;
		while (iterEdge.hasNext()) {
			Edge e= iterEdge.next();
			for (int i = 0; i < currentPartitionSet.size(); i++) {
				if (intersects(e.getTranscripts(), currentPartitionSet.elementAt(i).unity)) {
					splitPartitionSet.add(currentPartitionSet.elementAt(i).pathes);
					for (int j = 0; j < currentPartitionSet.elementAt(i).pathes.size(); j++) 
						unitedPartition.add(currentPartitionSet.elementAt(i).pathes.elementAt(j));
					if (unity== null)
						unity= currentPartitionSet.elementAt(i).unity;
					else
						unity= unite(unity, currentPartitionSet.elementAt(i).unity);
					currentPartitionSet.remove(i--);	// lemma
					break;
				}
			}
		}
		
		Partition theNewPartition= new Partition(unitedPartition, unity);
		currentPartitionSet.add(theNewPartition);
		
		return theNewPartition;
	}

	// search for the partitions that are 
	long[] splitPartitionsLast(Vector<Edge> inEdges, Vector<Vector<long[]>> currentPartitionSet, Vector<Vector<long[]>> splitPartitionSet) {
		Iterator<Edge> iterEdge= inEdges.iterator();
		Vector<long[]> unitedPartition= new Vector<long[]>();
		while (iterEdge.hasNext()) {
			Edge e= iterEdge.next();
			for (int i = 0; i < currentPartitionSet.size(); i++) {
				Vector<long[]> v= new Vector<long[]>(4,2);
				for (int j = 0; j < currentPartitionSet.elementAt(i).size(); j++) {
					long[] inter= intersect(e.getTranscripts(), currentPartitionSet.elementAt(i).elementAt(j));
					if (!isNull(inter)) {
						if (!equalSet(inter, currentPartitionSet.elementAt(i).elementAt(j)))
							System.err.println("partition error.");
						v.add(inter);
						unitedPartition.add(inter);
						currentPartitionSet.elementAt(i).remove(j--);	// == inter
						if (equalSet(inter, e.getTranscripts()))
							break;	// edge overlap multiple "pathes" of a partition
					}
				}
				if (v.size()== 0)
					continue;
				splitPartitionSet.add(v);
				if (currentPartitionSet.elementAt(i).size()== 0)
					currentPartitionSet.remove(i);
				break;	// edge can only overlap one partition
			}
		}
		currentPartitionSet.add(unitedPartition);
		
		return null;
	} 

	long[] splitPartitions_notSoOld(Vector<Edge> inEdges, Vector<long[]> currentPartitionSet, Vector<Vector<long[]>> splitPartitionSet) {
		Iterator<Edge> iterEdge= inEdges.iterator();
		long[] trinity= null;
		while (iterEdge.hasNext()) {
			long[] unity= null;
			Edge e= iterEdge.next();
			Vector<long[]> splitPartitions= new Vector<long[]>(4,2);
			for (int j = 0; j < currentPartitionSet.size(); j++) {
				long[] inter= intersect(e.getTranscripts(), currentPartitionSet.elementAt(j));
				if (isNull(inter))
					continue;
				if (unity== null)
					unity= inter;
				else
					unity= unite(unity,inter);
				splitPartitions.add(inter);
				if (equalSet(inter, currentPartitionSet.elementAt(j)))	// replace current partition
					currentPartitionSet.remove(j--);
				else
					currentPartitionSet.setElementAt(xor(currentPartitionSet.elementAt(j), inter), j);
				if (equalSet(unity, e.getTranscripts()))
					break;
			}
			splitPartitionSet.add(splitPartitions);
			if (trinity== null)
				trinity= unity;
			else
				trinity= unite(trinity, unity);
		}
		
		return trinity;
	}

	long[] splitPartitions_old(Vector<Edge> inEdges, Vector<long[]> currentPartitionSet, Vector<Vector<long[]>> splitPartitionSet) {
		Iterator<Edge> iterEdge= inEdges.iterator();
		long[] trinity= null;
		while (iterEdge.hasNext()) {
			long[] unity= null;
			Edge e= iterEdge.next();
			Vector<long[]> splitPartitions= new Vector<long[]>(4,2);
			for (int j = 0; j < currentPartitionSet.size(); j++) {
				long[] inter= intersect(e.getTranscripts(), currentPartitionSet.elementAt(j));
				if (isNull(inter))
					continue;
				if (unity== null)
					unity= inter;
				else
					unity= unite(unity,inter);
				splitPartitions.add(inter);
				if (equalSet(inter, currentPartitionSet.elementAt(j)))	// replace current partition
					currentPartitionSet.remove(j--);
				else
					currentPartitionSet.setElementAt(xor(currentPartitionSet.elementAt(j), inter), j);
				if (equalSet(unity, e.getTranscripts()))
					break;
			}
			splitPartitionSet.add(splitPartitions);
			if (trinity== null)
				trinity= unity;
			else
				trinity= unite(trinity, unity);
		}
		
		return trinity;
	}
	
	public void getEventsByPartitions(int n) {		
				if (nodeHash.size()== 0)
					return;
				Node[] nodes= getNodesInGenomicOrder();
				Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
				
				Vector<Partition> currentPartitionSet= new Vector<Partition>(trpts.length/ 2);
				Vector<long[]> initPartition= new Vector<long[]>();
				long[] all= createAllArray();
				initPartition.add(all);	// nodes[0].getTranscripts()
				currentPartitionSet.add(new Partition(initPartition, all));
				
				HashMap<Transcript, Node> validSinceMap= new HashMap<Transcript, Node>(trpts.length,1f);
				
				for (int i = 0; i < nodes.length; i++) {
	
						// check bubble
					Vector<Edge> inEdges= nodes[i].getInEdges();
					Partition unitedPartition= null;
					if (inEdges.size()>= n) {
						
							// split inPartitions in pathes
						Vector<Vector<long[]>> splitPartitionSet= new Vector<Vector<long[]>>(inEdges.size()); 
						unitedPartition= splitPartitions(inEdges, currentPartitionSet, splitPartitionSet);	// unity== nodes[i].getTranscripts(), .. deleted edges !!!
						//currentPartitionSet.add(unity);	// nodes[i].getTranscripts()																	// add new partition
						
							//generate combinations and look for lca
						generateTuples(n, nodes[i], splitPartitionSet, openSrcVec, validSinceMap);
					
						// re-surrect deleted partitions TODO in else part of indeg>= k
					} 
					
						// update valid pathes
					if (canonicalSS|| i== 0) {
						long[] unityIn= unite(nodes[i].getInPartitions());
						if (unityIn== null)
							unityIn= createNullArray(taSize);
						long[] unityOut= unite(nodes[i].getOutPartitions());
						if (unityOut== null)
							unityOut= createNullArray(taSize);
						long[] die= without(unityIn, unityOut);
						if (!isNull(die)) {
							Transcript[]  t= decodeTset(die);
							for (int j = 0; j < t.length; j++) 
								validSinceMap.remove(t[j]);
						}
						long[] born= without(unityOut, unityIn);
						if (!isNull(born)) {
							Transcript[]  t= decodeTset(born);
							for (int j = 0; j < t.length; j++) 
								validSinceMap.put(t[j],nodes[i]);
						}
					}
					
					
						// close open nodes, requires to split partitions at out edges
	//				for (int j = openSrcVec.size()- 1; j >= 0; --j) {
	//					int k = 0; 
	//					for (;k < currentPartitionSet.size(); k++) {
	//						long[] inter= intersect(openSrcVec.elementAt(j).getTranscripts(), currentPartitionSet.elementAt(k));
	//						if (equalSet(inter, openSrcVec.elementAt(j).getTranscripts()))
	//							break;
	//					}
	//					if (k< currentPartitionSet.size())
	//						openSrcVec.remove(j);
	//				}
					
					
						// splitting of edges will be handled by inedges?! deleted edges??! TODO
					Vector<Edge> outEdges= nodes[i].getOutEdges();
					if (outEdges.size()>= n) {
						openSrcVec.add(nodes[i]);
					}
					
							// split partitions
					long[] inter= null, rest= null;
					if (outEdges.size()> 1) {
						Vector<Vector<long[]>> v= new Vector<Vector<long[]>>(outEdges.size(),2);
						Vector<long[]> trinity= new Vector<long[]>(outEdges.size(),2);
							// look for common Partition
						int k= currentPartitionSet.size()- 1;
						if (unitedPartition== null) {
							for (int j = 0; j < currentPartitionSet.size(); j++) {
								if (intersects(outEdges.elementAt(0).getTranscripts(), currentPartitionSet.elementAt(j).unity)) {
									unitedPartition= currentPartitionSet.elementAt(j);
									k= j;
									break;
								}
							}
						}
						
						
						for (int j = 0; j < outEdges.size(); j++) {
							long[] unity= null; 	// all partitions combined by this edge
//							for (int k= 0; k < currentPartitionSet.size(); ++k) {	// extract from
//								if (!intersects(outEdges.elementAt(j).getTranscripts(),currentPartitionSet.elementAt(k).unity))
//										continue;
								Vector<long[]> newPartition= new Vector<long[]>();
								for (int m = 0; m < unitedPartition.pathes.size(); m++) {
									inter= intersect(unitedPartition.pathes.elementAt(m), outEdges.elementAt(j).getTranscripts());
									if (isNull(inter))
										continue;
									newPartition.add(inter);	// will not be iterated for this edge
									if (unity== null)
										unity= inter;
									else 
										unity= unite(unity, inter);
									rest= without(currentPartitionSet.elementAt(k).pathes.elementAt(m), outEdges.elementAt(j).getTranscripts());
									if (isNull(rest))
										unitedPartition.pathes.remove(m--);
									else
										unitedPartition.pathes.set(m, rest);
									if (equalSet(unity, outEdges.elementAt(j).getTranscripts()))
										break;
								}
								v.add(newPartition);
								trinity.add(unity);
								if (unitedPartition.pathes.size()== 0)
									currentPartitionSet.remove(k);
//							}
						}
						for (int j = 0; j < v.size(); j++) 
							currentPartitionSet.add(new Partition(v.elementAt(j), trinity.elementAt(j)));	// add not before, endless loops..
						
					}
					
				} // for all nodes
				
			}

	public void getEventsByPartitionsLast(int n) {		
				if (nodeHash.size()== 0)
					return;
				Node[] nodes= getNodesInGenomicOrder();
				Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
				
				Vector<Vector<long[]>> currentPartitionSet= new Vector<Vector<long[]>>(trpts.length/ 2);
				Vector<long[]> initPartition= new Vector<long[]>();
				initPartition.add(createAllArray());	// nodes[0].getTranscripts()
				currentPartitionSet.add(initPartition);
				
				HashMap<Transcript, Node> validSinceMap= new HashMap<Transcript, Node>(trpts.length,1f);
				
				for (int i = 0; i < nodes.length; i++) {

						// check bubble
					Vector<Edge> inEdges= nodes[i].getInEdges();
					if (inEdges.size()>= n) {
						
							// split inPartitions in pathes
						Vector<Vector<long[]>> splitPartitionSet= new Vector<Vector<long[]>>(inEdges.size()); 
						long[] unity= splitPartitions(inEdges, currentPartitionSet, splitPartitionSet);	// unity== nodes[i].getTranscripts(), .. deleted edges !!!
						//currentPartitionSet.add(unity);	// nodes[i].getTranscripts()																	// add new partition
						
							//generate combinations and look for lca
						generateTuples(n, nodes[i], splitPartitionSet, openSrcVec, validSinceMap);
					
						// re-surrect deleted partitions TODO in else part of indeg>= k
					} 
					
						// update valid pathes
					if (canonicalSS|| i== 0) {
						long[] unityIn= unite(nodes[i].getInPartitions());
						if (unityIn== null)
							unityIn= createNullArray(taSize);
						long[] unityOut= unite(nodes[i].getOutPartitions());
						if (unityOut== null)
							unityOut= createNullArray(taSize);
						long[] die= without(unityIn, unityOut);
						if (!isNull(die)) {
							Transcript[]  t= decodeTset(die);
							for (int j = 0; j < t.length; j++) 
								validSinceMap.remove(t[j]);
						}
						long[] born= without(unityOut, unityIn);
						if (!isNull(born)) {
							Transcript[]  t= decodeTset(born);
							for (int j = 0; j < t.length; j++) 
								validSinceMap.put(t[j],nodes[i]);
						}
					}
					
					
						// close open nodes, requires to split partitions at out edges
	//				for (int j = openSrcVec.size()- 1; j >= 0; --j) {
	//					int k = 0; 
	//					for (;k < currentPartitionSet.size(); k++) {
	//						long[] inter= intersect(openSrcVec.elementAt(j).getTranscripts(), currentPartitionSet.elementAt(k));
	//						if (equalSet(inter, openSrcVec.elementAt(j).getTranscripts()))
	//							break;
	//					}
	//					if (k< currentPartitionSet.size())
	//						openSrcVec.remove(j);
	//				}
					
					
						// splitting of edges will be handled by inedges?! deleted edges??! TODO
					Vector<Edge> outEdges= nodes[i].getOutEdges();
					if (outEdges.size()>= n) {
						openSrcVec.add(nodes[i]);
					}
					
							// split partitions
					long[] inter= null, rest= null;
					if (outEdges.size()> 1) {
						Vector<Vector<long[]>> v= new Vector<Vector<long[]>>(outEdges.size(),2);
						for (int j = 0; j < outEdges.size(); j++) {
							long[] unity= null; 	// all partitions combined by this edge
							for (int k= 0; k < currentPartitionSet.size(); ++k) {	// extract from
								Vector<long[]> newPartition= new Vector<long[]>();
								for (int m = 0; m < currentPartitionSet.elementAt(k).size(); m++) {
									inter= intersect(currentPartitionSet.elementAt(k).elementAt(m), outEdges.elementAt(j).getTranscripts());
									if (isNull(inter))
										continue;
									newPartition.add(inter);	// will not be iterated for this edge
									if (unity== null)
										unity= inter;
									else 
										unity= unite(unity, inter);
									rest= without(currentPartitionSet.elementAt(k).elementAt(m), outEdges.elementAt(j).getTranscripts());
									if (isNull(rest))
										currentPartitionSet.elementAt(k).remove(m--);
									else
										currentPartitionSet.elementAt(k).set(m, rest);
									if (equalSet(unity, outEdges.elementAt(j).getTranscripts()))
										break;
								}
								if (newPartition.size()> 0)
									v.add(newPartition);
								if (currentPartitionSet.elementAt(k).size()== 0)
									currentPartitionSet.remove(k--);
							}
						}
						for (int j = 0; j < v.size(); j++) 
							currentPartitionSet.add(v.elementAt(j));	// add not before, endless loops..
						
					}
					
				} // for all nodes
				
			}

	public void getEventsByPartitions_bug_in_valid_pathes(int n) {		
			if (nodeHash.size()== 0)
				return;
			Node[] nodes= getNodesInGenomicOrder();
			Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
			Vector<long[]> currentPartitionSet= new Vector<long[]>();
			currentPartitionSet.add(nodes[0].getTranscripts());
			for (int i = 0; i < nodes.length; i++) {
				
				
					// check bubble
				Vector<Edge> inEdges= nodes[i].getInEdges();
				if (inEdges.size()>= n) {
					
						// split inPartitions in pathes
					Vector<Vector<long[]>> splitPartitionSet= new Vector<Vector<long[]>>(inEdges.size()); 
					long[] unity= splitPartitions_old(inEdges, currentPartitionSet, splitPartitionSet);	// unity== nodes[i].getTranscripts(), .. deleted edges !!!
					currentPartitionSet.add(unity);	// nodes[i].getTranscripts()																	// add new partition
					
						//generate combinations and look for lca
					generateTuples_without_validity(n, nodes[i], splitPartitionSet, openSrcVec);
				
					// re-surrect deleted partitions TODO in else part of indeg>= k
				} 
				
				else // ? else
				
				{
					int j = 0;
					for (; j < currentPartitionSet.size(); j++) {
						long[] inter=  intersect(nodes[i].getTranscripts(), currentPartitionSet.elementAt(j));
						if (equalSet(inter, currentPartitionSet.elementAt(j))&& (!equalSet(inter, nodes[i].getTranscripts())))
							break;	// widening partition
					}
					if (j< currentPartitionSet.size()) 
						for (int k = 0; k < currentPartitionSet.size(); k++) {
							if (k== j)
								currentPartitionSet.setElementAt(nodes[i].getTranscripts(), k);
							else {
								long[] rest= without(currentPartitionSet.elementAt(k), nodes[i].getTranscripts());
								if (isNull(rest)) {
									currentPartitionSet.remove(k--);
									--j;
								} else
									currentPartitionSet.setElementAt(rest, k);
							}
						}
				}
				
					// close open nodes, requires to split partitions at out edges
//				for (int j = openSrcVec.size()- 1; j >= 0; --j) {
//					int k = 0; 
//					for (;k < currentPartitionSet.size(); k++) {
//						long[] inter= intersect(openSrcVec.elementAt(j).getTranscripts(), currentPartitionSet.elementAt(k));
//						if (equalSet(inter, openSrcVec.elementAt(j).getTranscripts()))
//							break;
//					}
//					if (k< currentPartitionSet.size())
//						openSrcVec.remove(j);
//				}
				
				
					// splitting of edges will be handled by inedges?! deleted edges??! TODO
				Vector<Edge> outEdges= nodes[i].getOutEdges();
				if (outEdges.size()>= n) {
					openSrcVec.add(nodes[i]);
						// split partitions, for killing open nodes
//					for (int j = 0; j < outEdges.size(); j++) {
//						for (int k = 0; k < currentPartitionSet.size(); ++k) {
//							long[] inter= intersect(currentPartitionSet.elementAt(k), outEdges.elementAt(j).getTranscripts());
//							if (isNull(inter))
//								continue;
//							long[] rest= without(currentPartitionSet.elementAt(k), outEdges.elementAt(j).getTranscripts());
//							currentPartitionSet.setElementAt(rest,k);
//							currentPartitionSet.add(inter);
//							break;	// can only intersect with one open partition
//						}
//					}
				}
				
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
