package gphase.graph9;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
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
import java.util.zip.GZIPOutputStream;

import qalign.tools.CedricConqueror;

import com.sun.org.apache.bcel.internal.Constants;
import com.sun.org.apache.xml.internal.utils.NodeVector;

import sun.reflect.ReflectionFactory.GetReflectionFactoryAction;

import gphase.Toolbox;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASEvent;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
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
 * problems after removing border edges on ESTs:
 * event.toString().equals("1^2-,1^2-3^4-") occur
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
				static int n= 2;
				boolean output= false, output2= true; 
				static Species species= new Species("human");
				static {
					species.setGenomeVersion("hg18");
				}
				
				public static void setSpecies(Species newSpecies) {
					species= newSpecies;
				}
				
				public EventExtractorThread(ChromosomeReaderThread newUpstreamThread) {
					super();
					setName("event_extraction_thread");
					upstreamThread= newUpstreamThread;					
				}
				
				@Override
				public void run() {

					String chromo= null;
					int evBefore= Graph.counter;
					long cumulGC= 0l, cumulGF= 0l, cumulGT= 0l, cumulEV= 0l;
					for (int i = 0; i < g.length; i++) {
						
						if (chromo== null)
							chromo= g[0].getChromosome();
												
						if (g[i].getTranscriptCount()== 1) {
							g[i]= null;
							continue;
						}
//						if (g[i].getTranscriptCount()> 5000)  {
//							if (output2) {
//								Date ti= new Date(System.currentTimeMillis());
//								System.out.println("["+ti+"] "+chromo+" skipped locus "+g[i].getTranscripts()[0].getTranscriptID());
//							}
//							continue;
//						}
						
						g[i].setSpecies(species);

						String tmp= null;
						try {
							long t1= System.currentTimeMillis();
							Graph gr= new Graph(g[i]);
							if (output) {
								System.out.println(gr.trpts[0].getTranscriptID()+"  transcripts "+g[i].getTranscriptCount());
							}
							//gr.init(n);
							gr.constructGraph();
							long t2= System.currentTimeMillis();
							long dA= (t2-t1);
							cumulGC+= dA;
							if (output) {
								tmp= gr.trpts[0].getTranscriptID()+" construct "+(dA/1000)+" n "+gr.nodeHash.size()+" e "+gr.edgeHash.size();
								System.out.println(tmp);
								System.out.flush();
							}
							
							
							//gr.collapseFuzzyFlanks();
							//gr.cleanESTborderEdges();
							long t22= System.currentTimeMillis();
							long dF= (t22-t2);
							cumulGF+= dF;
							if (output) {
								tmp= gr.trpts[0].getTranscriptID()+" flanks "+(dF/1000)+" n "+gr.nodeHash.size()+" e "+gr.edgeHash.size();
								System.out.println(tmp);
								System.out.flush();
							}
							
							gr.contractGraph(2);	// not n !!!
							long t3= System.currentTimeMillis();
							long dT= (t3-t22);
							cumulGT+= dT;
							if (output) {
								tmp= gr.trpts[0].getTranscriptID()+" contract "+(dT/1000)+" n "+gr.nodeHash.size()+" e "+gr.edgeHash.size();
								System.out.println(tmp);
								System.out.flush();
							}
							//int oldSize= evMap.size();
							
							gr.getEventsByPartitions(n);
							
							//cnt+= evMap.size()- oldSize;
							long dB= (System.currentTimeMillis()- t3);
							cumulEV+= dB;
//							if (output) 
//								//System.out.println("events "+(evMap.size()-oldSize));
//								System.out.println(" extract "+(dB/1000));
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
						int div= (int) (cumulGC+cumulGF+cumulEV)/1000;
						int frac= 0;
						if (div> 0)
							frac= (Graph.counter- evBefore)/ div;
						else
							frac= (Graph.counter- evBefore);
								
						System.out.println("["+ti+"] "+chromo+
								" graph construct "+(cumulGC/1000)+
								//" sec, fuzzy flanks "+(cumulGF/1000) +
								" sec, contraction "+(cumulGT/1000) +
								" sec, extract events "+(cumulEV/ 1000)+
								" sec, found "+(Graph.counter- evBefore)+
								" events, "+frac+" ev/sec.");
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
		boolean output= false, output2= true, checkIntrons= true; 
		
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
				
					// quickly init splice sites, no vale la pena
//				t0= System.currentTimeMillis();
//				for (int i = 0; i < g.length; i++) {
//					if (g[i].getTranscriptCount()== 1)
//						continue;
//					for (int j = 0; j < g[i].getTranscripts().length; j++) { 
//						SpliceSite[] ss= g[i].getTranscripts()[j].getSpliceSitesAll();
//						for (int x = 0; ss!= null&& x < ss.length; x++) {
//							ss[x].getTranscripts();
//						}
//					}
//				}
//				if (output2) {
//					Date ti= new Date(System.currentTimeMillis());
//					System.out.println("["+ti+"] ss init: "+g[0].getChromosome()+" "+((System.currentTimeMillis()- t0)/ 1000)+" sec.");
//				}
				
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
//				break;
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
						//BufferedWriter buffy= new BufferedWriter(new OutputStreamWriter(
						GZIPOutputStream zipperStream= new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(new File(outputFname), true)));
						while (!queue.isEmpty()) {
							String s= queue.poll();
							for (int i = 0; i < s.length(); i++) 
								zipperStream.write((byte) s.charAt(i));
							zipperStream.write('\n');
						}
						zipperStream.flush(); zipperStream.close();
						//buffy.flush(); buffy.close();
						
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

		public void run_writer() {
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

	static gphase.tools.File parseArguments(String[] args) {
		
		gphase.tools.File file= null;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("rusc")) {
				gphase.Constants.DATA_DIR= "/home/ug/msammeth";	// ugrusc, but for the nodes ug
				continue;
			}
			if (args[i].equals("-i")|| args[i].equals("--input")) {
				file= new gphase.tools.File(args[++i]);
				continue;
			}
			if (args[i].equals("-o")|| args[i].equals("--output")) {
				outputFname= args[++i]+".gz";
				continue;
			}
			if (args[i].equals("-g")|| args[i].equals("--genome")) {
				String[] s= args[++i].split("_");
				Species spe= new Species(s[0]);
				spe.setGenomeVersion(s[1]);
				EventExtractorThread.setSpecies(spe);
				continue;
			}
			if (args[i].equals("-k")|| args[i].equals("--dimension")) {
				try {
					int x= Integer.parseInt(args[++i]);
					EventExtractorThread.n= x;
				} catch (NumberFormatException e) {
					; // :)
				}
				continue;
			}
		}
		
		return file;
	}
	static void _240808_test_multithread(String[] args) {

		gphase.tools.File file= parseArguments(args);
		// 
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716.gtf		
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr11.gtf
		// /home/ug/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_splicedESTs_from_UCSC070716_chr6.gtf
		//
		// /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807.gtf
		// /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919.gtf
		// /home/ug/msammeth/annotations/mm8_0602_RefSeq_fromUCSC_070807_mRNAs_fromUCSC070919_splicedESTs_fromUCSC070919.gtf
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
			System.out.println("["+ti+"]  started, k= "+EventExtractorThread.n+"species"+EventExtractorThread.species+", input file "+file.getAbsolutePath()+", output file= "+outputFname);
		}
		GTFChrReader reader= new GTFChrReader(file.getAbsolutePath()); 
		//reader.sweepToChromosome("chr14");
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
	
	
	void generateTuples_old(int k, Node src, Node snk, Vector<Vector<long[]>> buckets) {
	
		int[] idx= new int[k];
		for (int i = 0; i < idx.length; i++) 
			idx[i]= i;
		
		int[] combi= new int[k];
		snk.partitions= new Vector(buckets.size());
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

	boolean checkValid(int[] idx, int[] combi,  Vector<Vector<Partition>> buckets) {
		// ensure there is no partition containing all elements
		// => intersection of parents is empty
		HashMap<PartitionSet, PartitionSet> partitions= (HashMap<PartitionSet, PartitionSet>) 
											buckets.elementAt(idx[0]).elementAt(combi[0]).parents.clone();
		for (int i = 1; i < idx.length; i++) {
			Object[] o= partitions.values().toArray();
			for (int j = 0; j < o.length; j++) {
				if (buckets.elementAt(idx[i]).elementAt(combi[i]).parents.get(o[j])== null)
					partitions.remove(o[j]);
				if (partitions.size()== 0)	// as soon as one is from a partition the others are not from
					return true;
			}
		}
		return false;
	}

	boolean checkValid(int[] idx, int[] idx2,  Partition[] playground) {
		// ensure there is no partition containing all elements
		// => intersection of parents is empty
		HashMap<PartitionSet, PartitionSet> partitions= (HashMap<PartitionSet, PartitionSet>) 
											playground[idx[0]].parents.clone();
		for (int i = 1; i < idx.length; i++) {
			Object[] o= partitions.values().toArray();
			for (int j = 0; j < o.length; j++) {
				if (playground[idx[i]].parents.get(o[j])== null)
					partitions.remove(o[j]);
				if (partitions.size()== 0)
					return true;
			}
		}
		
		for (int i = 0; i < idx2.length; i++) {
			Object[] o= partitions.values().toArray();
			for (int j = 0; j < o.length; j++) {
				if (playground[idx2[i]].parents.get(o[j])== null)
					partitions.remove(o[j]);
				if (partitions.size()== 0)
					return true;
			}
		}
	
		return false;
	}

	boolean checkValid_new(int pivot, int[] idx,  Partition[] playground) {
		// ensure there is no partition containing all elements
		// => intersection of parents is empty
		HashMap<PartitionSet, PartitionSet> partitions= (HashMap<PartitionSet, PartitionSet>) 
											playground[pivot].parents.clone();
		Object[] o;
		for (int i = 0; i < idx.length; i++) {
			o= partitions.values().toArray();
			for (int j = 0; j < o.length; j++) {
				if (playground[idx[i]].parents.get(o[j])== null)
					partitions.remove(o[j]);
				if (partitions.size()== 0)
					return true;
			}
		}
		
		return false;
	}
	
	void generateTuples_new(int k, Node src, Node snk, Vector<Vector<Partition>> buckets) {
	
		int s= 0;
		for (int i = 0; i < buckets.size(); i++) 
			s+= buckets.elementAt(i).size();
		Partition[] playground= new Partition[s];
		int[] borders= new int[buckets.size()];
		s= 0;		
		int t= 1;
		for (int i = 0; i < buckets.size(); i++) {
			for (int j = 0; j < buckets.elementAt(i).size(); j++) 
				playground[s++]= buckets.elementAt(i).elementAt(j);
			if (i< buckets.size()- 1)
				borders[t++]= s;
		}
		
		
		int[] idx= new int[k- 1];
		// for every pivot region
		for (int i = 1; i < borders.length; i++) {	// min 1 other 'region'
			
			int left= borders[i-1];
			
			if (left> playground.length- k)
				break;
			
			
			for (int j = left; j < borders[i]; j++) {	// pivot position
				
				if (j> playground.length- k)
					break;
				
				int x= j+1;
				for (int m = 0; m < idx.length; ++m) 
					idx[m]= x++;
				
				while (true) {

					if (checkValid(j, idx, playground)) {
						// extract event
						SpliceSite[][] ss= new SpliceSite[k][];
						Transcript[][] tt= new Transcript[k][];
						for (int h = 0; h < idx.length; h++) {
							tt[h]= decodeTset(playground[idx[h]].transcripts);
							ss[h]= tt[h][0].getSpliceSitesBetween(src.getSite(), snk.getSite());
						}
						tt[tt.length-1]= decodeTset(playground[j].transcripts);
						ss[ss.length-1]= tt[tt.length-1][0].getSpliceSitesBetween(src.getSite(), snk.getSite());
						ASEvent ev= new ASEvent(tt,ss);
						ev.setAnchors(src.getSite(), snk.getSite());
						boolean event= true;
						if (isRoot(src)|| isLeaf(snk))
							event= ev.isASevent();
						if (event) 
							outputEvent(ev);
					}

						// next
					int negPtr= idx.length-1;
					while (negPtr>= 0) {
						if (idx[negPtr]< playground.length- (k- 1- negPtr)) {
							int c= ++idx[negPtr];
							for (int m = negPtr+1; m < idx.length; m++) 
								idx[m]= ++c;
							break;
						} 
							// else
						--negPtr;
					}
					if (negPtr< 0)
						break;
			}
		}
					
	}
}

	void generateTuples(int k, Node src, Node snk, Vector<Vector<Partition>> buckets) {
	
		int s= 0;
		for (int i = 0; i < buckets.size(); i++) 
			s+= buckets.elementAt(i).size();
		if (k<= 1)
			k= s;
		Partition[] playground= new Partition[s];
		int[] borders= new int[buckets.size()];		
		s= 0;		
		int t= 1;
		for (int i = 0; i < buckets.size(); i++) {
			for (int j = 0; j < buckets.elementAt(i).size(); j++) 
				playground[s++]= buckets.elementAt(i).elementAt(j);
			if (i< buckets.size()- 1)
				borders[t++]= s;
		}
		
		int[] idx, idx2;
		// for every pivot region
		for (int i = 1; i < borders.length; i++) {	// min 1 other 'region'
			int left= borders[i-1];
			if(playground.length- left< k)
				break;
			
			// now combine all sep combinations in the pivot partition with all (k-sep) combination from the rest
			for (int sep = 1; sep < k; sep++) {	// min 1 has to be separated
				
				if (sep> (borders[i]- left)|| (k-sep)> (playground.length- borders[i]))	// no single combi possible in one of the parted sets
					continue;	// no break, no?!
				
				idx= new int[sep];	// init pivot partition combination
				int x= left;
				for (int j = 0; j < idx.length; ++j) 
					idx[j]= x++;
				
				
				idx2= new int[k-sep];
				//for (int j = left; j < borders[i]- sep; j++) {
				while (true) {	// iterate combis on pivot partition, actually idx[0]< borders[i]- sep, but performed by negPtr	
					
					x= borders[i];
					for (int j = 0; j < idx2.length; j++)	// init rest partition 
						idx2[j]= x++;
					
					while (true) {	// iterate combis on rest partition
						
						if (checkValid(idx, idx2, playground)) {
							
								// get transcripts
							Transcript[][] tt= new Transcript[k][];
							for (int h = 0; h < idx.length; h++) 
								tt[h]= decodeTset(playground[idx[h]].transcripts);
							for (int h = 0; h < idx2.length; h++) 
								tt[idx.length+h]= decodeTset(playground[idx2[h]].transcripts);
							
							// omit boundaries of ESTs
							SpliceSite srcSite= src.getSite(), snkSite= snk.getSite();
							boolean correctEST= false, startOver= false;
//							if (isRoot(src)|| isLeaf(snk)) {
//								for (int j = 0; j < tt.length; j++) {
//									int y= 0;
//									for (; y < tt[j].length; ++y) {
//										if (tt[j][y].getSource().contains("Est")) {
//											correctEST= true;
//											break;
//										}
//									}
//									if (y< tt[j].length) {
//										if (isRoot(src)&& tt[j][0].getSpliceSitesAll()[0].getPos()> srcSite.getPos())
//											srcSite= tt[j][0].getSpliceSitesAll()[0];
//										if (isLeaf(snk)&& tt[j][0].getSpliceSitesAll()[tt[j][0].getSpliceSitesAll().length- 1].getPos()< snkSite.getPos())
//											snkSite= tt[j][0].getSpliceSitesAll()[tt[j][0].getSpliceSitesAll().length- 1];
//									}
//								}
//								if (srcSite.getPos()> snkSite.getPos())
//									startOver= true;
//							}
							
							// extract event
							SpliceSite[][] ss= new SpliceSite[k][];
							if (!startOver) {
								for (int h = 0; h < idx.length; h++) 
									ss[h]= tt[h][0].getSpliceSitesBetween(srcSite, snkSite);
								for (int h = 0; h < idx2.length; h++) 
									ss[idx.length+ h]= tt[idx.length+ h][0].getSpliceSitesBetween(srcSite, snkSite);
//								if (correctEST) {	// check for 2 or more identical splice chains now
//									int j = 0;
//									for (; j < ss.length; j++) {
//										int h = j+1;
//										for (; h < ss.length; h++) {
//											if (SpliceSite.getDefaultSpliceChainComparator().compare(ss[j], ss[h])== 0)
//												break;
//										}
//										if (h< ss.length)
//											break;
//									}
//									if (j< ss.length)
//										startOver= true;
//								}
							}
							
							if (!startOver) {
								ASEvent ev= new ASEvent(tt,ss);
								ev.setAnchors(srcSite, snkSite);
//								if (ev.toString().equals("1*2^,3*4^"))
//									System.currentTimeMillis();
								boolean event= true;
								if (isRoot(src)|| isLeaf(snk))
									event= ev.isASevent();
								if (event) 
									outputEvent(ev, playground.length);	// not bucket-size
								
//								else
//									System.currentTimeMillis();
							}
						}
					
						// next combi in rest partition
						int negPtr= idx2.length-1;
						while (negPtr>= 0) {
							if (idx2[negPtr]< playground.length- ((k-sep)- negPtr)) {
								int c= ++idx2[negPtr];
								for (int m = negPtr+1; m < idx2.length; m++) 
									idx2[m]= ++c;
								break;
							} 
								// else
							--negPtr;
						}
						if (negPtr< 0)
							break;
					}
					
						// next pivot partition combination
					int negPtr= idx.length-1;
					while (negPtr>= 0) {
						if (idx[negPtr]< borders[i]- (sep- negPtr)) {
							int c= ++idx[negPtr];
							for (int m = negPtr+1; m < idx.length; m++) 
								idx[m]= ++c;
							break;
						} 
							// else
						--negPtr;
					}
					if (negPtr< 0)
						break;
				}
			}
		}
		
	}

	void generateTuples_working_pw(int k, Node src, Node snk, Vector<Vector<Partition>> buckets) {
	
		int[] idx= new int[k];
		for (int i = 0; i < idx.length; i++) 
			idx[i]= i;
		
		int[] combi= new int[k];
		snk.partitions= new Vector(buckets.size());
		while (idx[0]< (buckets.size()-k+ 1)) {

			
				// now all combinations	TODO multiple from same partition!!
			for (int j = 0; j < combi.length; j++) 
				combi[j]= 0;
			while (true) {
				
					// not necessary, no ?!
//				long[] combination= null;
//				for (int j = 0; j < combi.length; j++) 
//					if (combination== null)
//						combination= buckets.elementAt(idx[j]).elementAt(combi[j]).transcripts;
//					else
//						combination= unite(combination, buckets.elementAt(idx[j]).elementAt(combi[j]).transcripts);
					
				if (checkValid(idx, combi, buckets)) {
						// extract event
						SpliceSite[][] ss= new SpliceSite[k][];
						Transcript[][] tt= new Transcript[k][];
						for (int h = 0; h < combi.length; h++) {
							tt[h]= decodeTset(buckets.elementAt(idx[h]).elementAt(combi[h]).transcripts);
							ss[h]= tt[h][0].getSpliceSitesBetween(src.getSite(), snk.getSite());
						}
						ASEvent ev= new ASEvent(tt,ss);
						ev.setAnchors(src.getSite(), snk.getSite());
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
	boolean acceptableIntrons= true;
	
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
	
	public Edge createEdge(Node v, Node w, long[] newTset, byte type) {
		Edge e= getEdge(v,w);
		if (e== null) {
			e= new Edge(v,w);
			e.type= type;
			e.setTranscripts(newTset);
			edgeHash.put(v.getSite().toString()+w.getSite().toString(),e);
			if (acceptableIntrons) {
				if (v.getSite().isDonor()) {
					if (w.getSite().getPos()- v.getSite().getPos()< Exon.MIN_INTRON_LENGTH_HUMAN)
						e.valid= false;
					if (e.valid) {
						if (!Exon.checkAcceptableIntron(v.getSite(), w.getSite()))
							e.valid= false;
					}
				}
			}
		} else {
			e.setTranscripts(unite(e.getTranscripts(),newTset));
			if (type< e.type)
				e.type= type;
		}
		return e;
	}

	public Edge addEdge(Node v, Node w, long[] newTset) {
		Edge e= new Edge(v,w);
		e.setTranscripts(newTset);
		Edge chk= edgeHash.get(v.getSite().toString()+w.getSite().toString());
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
					createEdge(v,w,encodeTset(new Transcript[] {t[i]}),t[i].getSourceType());
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
				byte type= Byte.MAX_VALUE;
				for (int j = 0; j < tt.length; j++) { 
					tt[j]= v.elementAt(j);
	//				if (tt[j].getTranscriptID().equals("DA126574"))
	//					System.currentTimeMillis();
					if (tt[j].getSourceType()< type)
						type= tt[j].getSourceType();
				}
				
				if (type< Transcript.ID_SRC_MRNA)
					createEdge(root, (Node) o[i], encodeTset(tt), type);
				else if (((Node) o[i]).getInEdges().size()== 0)
					removeEstChain((Node) o[i], true);
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
				byte type= Byte.MAX_VALUE;
				for (int j = 0; j < tt.length; j++) { 
					tt[j]= v.elementAt(j);
					if (tt[j].getSourceType()< type)
						type= tt[j].getSourceType();
				}
				long[] temp= encodeTset(tt);
				
				if (type< Transcript.ID_SRC_MRNA)
					createEdge((Node) o[i], leaf, encodeTset(tt), type);
				else if (((Node) o[i]).getOutEdges().size()== 0)
					removeEstChain((Node) o[i], false);
			}
		}

	void constructGraph_graph9() {
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
				createEdge(v,w,encodeTset(new Transcript[] {t[i]}),t[i].getSourceType());
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
			byte type= Byte.MAX_VALUE;
			for (int j = 0; j < tt.length; j++) { 
				tt[j]= v.elementAt(j);
//				if (tt[j].getTranscriptID().equals("DA126574"))
//					System.currentTimeMillis();
				if (tt[j].getSourceType()< type)
					type= tt[j].getSourceType();
			}
			
			if (type< Transcript.ID_SRC_EST)
				createEdge(root, (Node) o[i], encodeTset(tt), type);
			else 
				removeEstChain((Node) o[i], true);
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
			byte type= Byte.MAX_VALUE;
			for (int j = 0; j < tt.length; j++) { 
				tt[j]= v.elementAt(j);
				if (tt[j].getSourceType()< type)
					type= tt[j].getSourceType();
			}
			long[] temp= encodeTset(tt);
			
			if (type< Transcript.ID_SRC_EST)
				createEdge((Node) o[i], leaf, encodeTset(tt), type);
			else 
				removeEstChain((Node) o[i], false);
		}
	}
	
	void removeEstChain(Node n, boolean fromFront) {
			
				// stop
			int min= 2;
			Transcript[] t= decodeTset(n.getTranscripts());
			if (t.length>= min)
				return;
			
				// do
			Object[] edges;
			if (fromFront)
				edges= n.getOutEdges().toArray();
			else
				edges= n.getInEdges().toArray();
			boolean doit= true;
			for (int i = 0; i < edges.length; i++) {
				Edge e= (Edge) edges[i];
				if (e.type>  Transcript.iD_SRC_REFSEQ) { // was: == ID_SRC_MRNA
					Node nextN= null;
					if (fromFront)
						nextN= e.getHead();
					else
						nextN= e.getTail();
					
					int x;
					if (fromFront)
						x= nextN.getInEdges().size();
					else
						x= nextN.getOutEdges().size();
					
					removeEdge(e);
					if (n.getInEdges().size()== 0&& n.getOutEdges().size()== 0)
						removeNode(n);
					
					if (x== 1)
						removeEstChain(nextN, fromFront);
				}
			}
		}

	void removeEstChain_rigorous(Node n, boolean fromFront) {
			
				// stop
	//		Collection c;
	//		if (fromFront)
	//			c= n.getInEdges();
	//		else
	//			c= n.getOutEdges();
	//		if (c.size()> 1)
	//			return;
			
	//		Iterator<Edge> iter;
	//		if (fromFront)
	//			iter= n.getInEdges().iterator();
	//		else
	//			iter= n.getOutEdges().iterator();
	//		while (iter.hasNext()) {
	//			if (iter.next().type< Transcript.ID_SRC_EST)
	//				return;
	//		}
				
				// do
			Object[] edges;
			if (fromFront)
				edges= n.getOutEdges().toArray();
			else
				edges= n.getInEdges().toArray();
			for (int i = 0; i < edges.length; i++) {
				Edge e= (Edge) edges[i];
				if (e.type>  Transcript.iD_SRC_REFSEQ) { // was: == ID_SRC_MRNA
					Node nextN= null;
					if (fromFront)
						nextN= e.getHead();
					else
						nextN= e.getTail();
					
					int x;
					if (fromFront)
						x= nextN.getInEdges().size();
					else
						x= nextN.getOutEdges().size();
					
					removeEdge(e);
					if (n.getInEdges().size()== 0&& n.getOutEdges().size()== 0)
						removeNode(n);
					
					if (x== 1)
						removeEstChain(nextN, fromFront);
				}
			}
		}

	void removeEstChain_graph9(Node n, boolean fromFront) {
		
			// stop
//		Collection c;
//		if (fromFront)
//			c= n.getInEdges();
//		else
//			c= n.getOutEdges();
//		if (c.size()> 1)
//			return;
		
//		Iterator<Edge> iter;
//		if (fromFront)
//			iter= n.getInEdges().iterator();
//		else
//			iter= n.getOutEdges().iterator();
//		while (iter.hasNext()) {
//			if (iter.next().type< Transcript.ID_SRC_EST)
//				return;
//		}
			
			// do
		Object[] edges;
		if (fromFront)
			edges= n.getOutEdges().toArray();
		else
			edges= n.getInEdges().toArray();
		for (int i = 0; i < edges.length; i++) {
			Edge e= (Edge) edges[i];
			if (e.type==  Transcript.ID_SRC_EST) {
				Node nextN= null;
				if (fromFront)
					nextN= e.getHead();
				else
					nextN= e.getTail();
				
				int x;
				if (fromFront)
					x= nextN.getInEdges().size();
				else
					x= nextN.getOutEdges().size();
				
				removeEdge(e);
				if (n.getInEdges().size()== 0&& n.getOutEdges().size()== 0)
					removeNode(n);
				
				if (x== 1)
					removeEstChain(nextN, fromFront);
			}
		}
	}
	
	// better on graph-level first/last two edges
	// (alignment errors)
	void collapseFuzzyFlanks(Transcript[] trpts) {
			// filter first exon
		for (int i = 0; i < trpts.length; i++) {
//			if (!trpts[i].getSource().contains("EST"))
//				continue;
			for (int j = (i+1); j < trpts.length; j++) {
				if (trpts[i].getExons().length== trpts[j].getExons().length)
					continue;
				int k = 0;
				for (; k < trpts[i].getExons().length; k++) {
					if (i> 0&& trpts[i].getExons()[k].get5PrimeEdge()!= trpts[j].getExons()[k].get3PrimeEdge())
						break;
					//if (i< trpts[i].getExons().length- 1&& trpts[i].getExons()[k].get5PrimeEdge()!= trpts[j].getExons()[k].get3PrimeEdge())
				}
			}
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
			Transcript[] t= decodeTset(e.getTranscripts());
			int x= 0;
			for (; x < t.length; x++) 
				if (!t[x].getSource().contains("Est"))
					break;
			if (x< t.length)
				continue;
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
					addEdge(root, n, unity);
				else
					addEdge(n, leaf, unity);
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
	
	void cleanESTborderEdges() {
		cleanESTborderEdges(root.getOutEdges());
		cleanESTborderEdges(leaf.getInEdges());
	}
	
	void cleanESTborderEdges(Vector<Edge> edges) {
		Iterator<Edge> iter= edges.iterator();
		while (iter.hasNext()) {
			Edge e= iter.next();
			Transcript[] t= decodeTset(e.getTranscripts());
			int i = 0;
			for (;  i < t.length; i++) 
				if (!t[i].getSource().contains("EST"))
					break;
			if (i== t.length)
				removeEdge(e);
		}
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
		boolean validity= srcEdge.valid;
		while (outV.size()< k&& inV.size()< k&& outV.size()> 0) {
			f= outV.iterator().next();
			newTset= intersect(newTset, f.getTranscripts());
			if (isNull(newTset))
				break;	// because of deleted edges
			endEdge= f;
			validity&= endEdge.valid;
			n= endEdge.getHead();
			outV= n.getOutEdges();
			inV= n.getInEdges();			
		}
		
			// contract
		if (endEdge!= null) {

			f= new Edge(srcEdge.getTail(),endEdge.getHead());
			f.setTranscripts(newTset);	// srcEdge.getTranscripts() not, since edges are deleted; really TODO
			f.setContracted(true);
			f.setProcessed(true);
			f.valid= validity;
			f.type= srcEdge.type;
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
	
		public void getEventsByPartitions_1st_Sep(int n) {		
					if (nodeHash.size()== 0)
						return;
					Node[] nodes= getNodesInGenomicOrder();
					Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
					
					HashMap<Transcript, Node> validSinceMap= new HashMap<Transcript, Node>(trpts.length,1f);
					
					long[] inter, unity, without;
					for (int i = 0; i < nodes.length; i++) {
		
							// check bubble
						Vector<Edge> inEdges= nodes[i].getInEdges();
						if (inEdges.size()>= n) {
							
							Vector<Partition> partitions= new Vector<Partition>(inEdges.size());	// TODO hashmaps
							Vector<PartitionSet> partitionSets= new Vector<PartitionSet>(inEdges.size());
							for (int j = 0; j < inEdges.size(); j++) { 
									Partition p= new Partition();
									p.transcripts= inEdges.elementAt(j).getTranscripts();
									PartitionSet s= new PartitionSet();
									p.addParent(s);
									partitions.add(p);
									partitionSets.add(s);
							}
							
								
							//for (int j = openSrcVec.size()-1; j >= 0; --j) {	// openSrcVec.elementAt(j)
							for (int j = i-1; j >= 0;  --j) {
								
								if (!intersects(nodes[j].getTranscripts(), nodes[i].getTranscripts()))
									continue;
								
									// check for invalid introns
								if (acceptableIntrons) {
									for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {
										if (!nodes[j].getOutEdges().elementAt(m).valid) {
											for (int k = 0; k < partitions.size(); k++) {
												inter= intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
												if (isNull(inter))
													continue;
												without= without(partitions.elementAt(k).transcripts, inter);
												if (isNull(without)) {
													Iterator<PartitionSet> iter= partitions.elementAt(k).parents.keySet().iterator();
													while (iter.hasNext()) {
														PartitionSet ps= iter.next();
														ps.partitions.remove(partitions.elementAt(k));
														if (ps.partitions.size()== 0)
															partitionSets.remove(ps);
													}
													partitions.remove(k--);
													if (partitions.size()== 0)
														break;
												} else 
													partitions.elementAt(k).transcripts= without;
												
											}
											if (partitions.size()== 0)
												break;
										}
										if (partitions.size()== 0)
											break;
									}
								}
								if (partitions.size()== 0)
									break;
								
								if (nodes[j].getOutEdges().size()<= 1)
									continue;
								
								Vector<Vector<Partition>> splitPathes= new Vector<Vector<Partition>>(partitions.size());
								for (int k = 0; k < partitions.size(); k++) {
									Vector<Partition> newPartitions= new Vector<Partition>();	// TODO size
									for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {
										if (!nodes[j].getOutEdges().elementAt(m).valid)
											continue;
										
										// TODO check for equalset ?
										inter= intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
										if (isNull(inter))
											continue;
										
										without= without(partitions.elementAt(k).transcripts, inter);
										if (isNull(without)) {
											newPartitions.add(partitions.remove(k--));	// just temporary remove, parent cannot disappear
											break;
										} else {
											Partition newPartition= (Partition) partitions.elementAt(k).clonePartition();
											newPartition.transcripts= inter;
											newPartitions.add(newPartition);	// new partition
											partitions.elementAt(k).transcripts= without;
										}
									}
									if (newPartitions.size()> 0)
										splitPathes.add(newPartitions);
								}
								
									// now add
								for (int k = 0; k < splitPathes.size(); k++) {
									for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
										partitions.add(splitPathes.elementAt(k).elementAt(h));
									}
								}
								
									// combinations
								if (splitPathes.size()>= n) {
									generateTuples(n, nodes[j], nodes[i], splitPathes);
								}
								
								// create a new partition set
								if (splitPathes.size()> 1) {
									PartitionSet newSet= new PartitionSet();
									partitionSets.add(newSet);
									for (int k = 0; k < splitPathes.size(); k++) {
										for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
											splitPathes.elementAt(k).elementAt(h).addParent(newSet);
										}
									}
								}
								
									// stop
								inter= intersect(nodes[j].getTranscripts(),nodes[i].getTranscripts());
								if (equalSet(inter, nodes[i].getTranscripts()))
									break;
							}
						} 
			
						
						
							// update valid pathes
	//					if (false&& canonicalSS|| i== 0) {
	//						long[] unityIn= unite(nodes[i].getInPartitions());
	//						if (unityIn== null)
	//							unityIn= createNullArray(taSize);
	//						long[] unityOut= unite(nodes[i].getOutPartitions());
	//						if (unityOut== null)
	//							unityOut= createNullArray(taSize);
	//						long[] die= without(unityIn, unityOut);
	//						if (!isNull(die)) {
	//							Transcript[]  t= decodeTset(die);
	//							for (int j = 0; j < t.length; j++) 
	//								validSinceMap.remove(t[j]);
	//						}
	//						long[] born= without(unityOut, unityIn);
	//						if (!isNull(born)) {
	//							Transcript[]  t= decodeTset(born);
	//							for (int j = 0; j < t.length; j++) 
	//								validSinceMap.put(t[j],nodes[i]);
	//						}
	//					}
						
						
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
						
						if (nodes[i].getOutEdges().size()> 1)
							openSrcVec.add(nodes[i]);
						
					} // for all nodes
					
				}

	public void getEventsByPartitions(int n) {		
		
				boolean outputCombinations= true;
				if (n<= 1) {
					outputCombinations= false;
					n= 2;
				}
		
				if (nodeHash.size()== 0)
					return;
				Node[] nodes= getNodesInGenomicOrder();
				//Vector<Node> openSrcVec= new Vector<Node>(nodes.length/ 2);
				
				HashMap<Transcript, Node> validSinceMap= new HashMap<Transcript, Node>(trpts.length,1f);
				
				long[] inter, unity, without;
				HashMap<PartitionSet, Integer> map;
				for (int i = 0; i < nodes.length; i++) {
	
						// check bubble
					Vector<Edge> inEdges= nodes[i].getInEdges();
					if (inEdges.size()>= n) {
						
						Vector<Partition> partitions= new Vector<Partition>(inEdges.size());	// TODO hashmaps
						Vector<PartitionSet> partitionSets= new Vector<PartitionSet>(inEdges.size());
						for (int j = 0; j < inEdges.size(); j++) { 
								Partition p= new Partition();
								p.transcripts= inEdges.elementAt(j).getTranscripts();
								PartitionSet s= new PartitionSet();
								p.addParent(s);
								partitions.add(p);
								partitionSets.add(s);
						}
						
							
						//for (int j = openSrcVec.size()-1; j >= 0; --j) {	// openSrcVec.elementAt(j)
						for (int j = i-1; j >= 0;  --j) {
							
							if (!intersects(nodes[j].getTranscripts(), nodes[i].getTranscripts()))
								continue;
							
								// check for invalid introns TODO outside the (indegree>= n) condition
							if (acceptableIntrons) {
								for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {
									if (!nodes[j].getOutEdges().elementAt(m).valid) {
										for (int k = 0; k < partitions.size(); k++) {
											inter= intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
											if (isNull(inter))
												continue;
											without= without(partitions.elementAt(k).transcripts, inter);
											if (isNull(without)) {
												Iterator<PartitionSet> iter= partitions.elementAt(k).parents.keySet().iterator();
												while (iter.hasNext()) {
													PartitionSet ps= iter.next();
													ps.partitions.remove(partitions.elementAt(k));
													if (ps.partitions.size()== 0)
														partitionSets.remove(ps);
												}
												partitions.remove(k--);
												if (partitions.size()== 0)
													break;
											} else 
												partitions.elementAt(k).transcripts= without;
											
										}
										if (partitions.size()== 0)
											break;
									}
									if (partitions.size()== 0)
										break;
								}
							}
							if (partitions.size()== 0)
								break;
							
							if (nodes[j].getOutEdges().size()<= 1)
								continue;
							
								// split partitions
							Vector<Vector<Partition>> splitPathes= new Vector<Vector<Partition>>(partitions.size());
							for (int k = 0; k < partitions.size(); k++) {
								Vector<Partition> newPartitions= new Vector<Partition>();	// TODO size
								for (int m = 0; m < nodes[j].getOutEdges().size(); m++) {
									if (!nodes[j].getOutEdges().elementAt(m).valid)
										continue;
									
									// TODO check for equalset ?
									inter= intersect(partitions.elementAt(k).transcripts, nodes[j].getOutEdges().elementAt(m).getTranscripts());
									if (isNull(inter))
										continue;
									
									without= without(partitions.elementAt(k).transcripts, inter);
									if (isNull(without)) {
										newPartitions.add(partitions.remove(k--));	// just temporary remove, parent cannot disappear
										break;
									} else {
										Partition newPartition= (Partition) partitions.elementAt(k).clonePartition();
										newPartition.transcripts= inter;
										newPartitions.add(newPartition);	// new partition
										partitions.elementAt(k).transcripts= without;
									}
								}
								if (newPartitions.size()> 0)
									splitPathes.add(newPartitions);
							}
							
								// now add new partitions
							for (int k = 0; k < splitPathes.size(); k++) {
								for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
									partitions.add(splitPathes.elementAt(k).elementAt(h));
								}
							}
							
								// combinations
							if (splitPathes.size()>= n) {
								if (outputCombinations)
									generateTuples(n, nodes[j], nodes[i], splitPathes);
								else {
									generateTuples(-1, nodes[j], nodes[i], splitPathes);
								}
							}
							
							// create a new partition set
							if (splitPathes.size()> 1) {
								PartitionSet newSet= new PartitionSet();
								partitionSets.add(newSet);
								map= new HashMap<PartitionSet, Integer>();
								Iterator<PartitionSet> iter;
								PartitionSet pset; 
								for (int k = 0; k < splitPathes.size(); k++) {
									for (int h = 0; h < splitPathes.elementAt(k).size(); h++) {
										iter= splitPathes.elementAt(k).elementAt(h).parents.keySet().iterator();
										while (iter.hasNext()) {
											pset= iter.next();
											if (pset.partitions.get(splitPathes.elementAt(k).elementAt(h))== null)
												continue;
											if (map.get(pset)== null)
												map.put(pset, new Integer(1));
											else
												map.put(pset, new Integer(map.get(pset).intValue()+ 1));
										}
										splitPathes.elementAt(k).elementAt(h).addParent(newSet);
									}
								}
								
								// remove un-needed partition-sets
								iter= map.keySet().iterator();
								while (iter.hasNext()) {
									pset= iter.next();
									if (pset.partitions.size()== map.get(pset).intValue()) {
										Object[] o= pset.partitions.keySet().toArray();
										for (int k = 0; k < o.length; k++) 
											((Partition) o[k]).parents.remove(pset);
										pset.partitions= null;
										partitionSets.remove(pset);
									}
								}
							}
							
							
							
								// stop
							inter= intersect(nodes[j].getTranscripts(),nodes[i].getTranscripts());
							if (equalSet(inter, nodes[i].getTranscripts()))
								break;
						}
					} 
		
					
					
						// update valid pathes
//					if (false&& canonicalSS|| i== 0) {
//						long[] unityIn= unite(nodes[i].getInPartitions());
//						if (unityIn== null)
//							unityIn= createNullArray(taSize);
//						long[] unityOut= unite(nodes[i].getOutPartitions());
//						if (unityOut== null)
//							unityOut= createNullArray(taSize);
//						long[] die= without(unityIn, unityOut);
//						if (!isNull(die)) {
//							Transcript[]  t= decodeTset(die);
//							for (int j = 0; j < t.length; j++) 
//								validSinceMap.remove(t[j]);
//						}
//						long[] born= without(unityOut, unityIn);
//						if (!isNull(born)) {
//							Transcript[]  t= decodeTset(born);
//							for (int j = 0; j < t.length; j++) 
//								validSinceMap.put(t[j],nodes[i]);
//						}
//					}
					
					
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
					
//					if (nodes[i].getOutEdges().size()> 1)
//						openSrcVec.add(nodes[i]);
					
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
//		if (e.getTail().getSite().getPos()== -224552058&& e.getHead().getSite().getPos()== -224552053)
//			System.currentTimeMillis();
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
	
	private void outputEvent(final ASEvent event, int dim) {
		writerThread.addEvent(dim+"\t"+event.toStringASTA());
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
