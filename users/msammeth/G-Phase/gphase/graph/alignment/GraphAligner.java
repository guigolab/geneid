package gphase.graph.alignment;

import java.util.Arrays;
import java.util.Comparator;
import java.util.PriorityQueue;
import java.util.Vector;

// import prefuse.data.Graph;

import gphase.algo.ASAnalyzer;
import gphase.db.EnsemblDBAdaptor;
import gphase.graph.SpliceGraph;
import gphase.graph.SpliceNode;
import gphase.graph.gui.GraphView;
import gphase.graph.gui.PFGraph;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.GraphHandler;
import gphase.model.Species;
import gphase.model.Transcript;

public class GraphAligner {

	public static void alignENSEMBLHomologGenes() {
		EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
		Graph g= adaptor.getGraphAllHomologs(EnsemblDBAdaptor.SPECIES_ISMB);
//		g.filterNonsense();
//		GraphHandler.writeOut(g, "graph_filt-nons.oos");
//		g.filterNonCodingTranscripts();
//		GraphHandler.writeOut(g, "graph_filt-nons_filt-nc.oos");
//
		//EnsemblDBAdaptor.removeNotAllHomologGenes(g);
		// "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish", "fruitfly", "mosquito"
		// 5e "human", "mouse", "rat", "cow", "dog", "chicken", "frog", "zebrafish"
		String[] speNames= new String[] {"human", "mouse", "dog"};
		EnsemblDBAdaptor.removeNotAllHomologGenes(g, speNames);
		//GraphHandler.writeOut(g, "graph_filt-nons_filt-nc_remAll.oos");

		int ctr= 0;
		for (int i = 0; i < speNames.length; i++) {
			Species sp= g.getSpeciesByName(Species.getBinomialForCommonName(speNames[i]));
			Gene[] ge= sp.getGenes();
			for (int j = 0; j < ge.length; j++) {
				//System.out.println(ge[j].getGeneID());
				++ctr;
			}
		}
		System.out.println("==\nHomolog gene nb: "+ctr+"/"+speNames.length+"="+(ctr/speNames.length));
		
		Species sp= g.getSpeciesByName(Species.getBinomialForCommonName(speNames[0]));
		Gene[] ge= sp.getGenes();
		Vector v= new Vector(ctr/speNames.length);
		for (int j = 0; j < ge.length; j++) 
			v.add(ge[j]);
		
		alignHomologs((Gene[]) gphase.tools.Arrays.toField(v));
		
	}
	
	static void alignHomologs(Gene[] hGenes) {
		for (int i = 0; i < hGenes.length; i++) {
			for (int j = i+1; j < hGenes.length; j++) {
				SpliceGraph g1= new SpliceGraph(hGenes[i].getTranscripts());
				g1.init();
				SpliceGraph g2= new SpliceGraph(hGenes[j].getTranscripts());
				g2.init();
				Mapping[] maps= align(g1, g2, -1);
				if (maps!= null)
					System.out.println(maps[0].toString());
			}
		}
	}
	
	public static void main(String[] args) {
		alignENSEMBLHomologGenes();
	}
	
	public static void test() {
		Species spec= new Species("human");
		Gene ge1= new Gene(spec, "ge1");
		ge1.setChromosome("1");
		ge1.setStrand(1);
		Transcript t11= new Transcript(ge1, "t11");
		t11.setStrand(1);
		t11.addExon(new Exon(t11, "e11_1", 1, 2));
		t11.addExon(new Exon(t11, "e11_2", 3, 4));
		t11.addExon(new Exon(t11, "e11_3", 500, 600));
		t11.addExon(new Exon(t11, "e11_4", 800, 900));
		ge1.addTranscript(t11);
		Transcript t12= new Transcript(ge1, "t12");
		t12.setStrand(1);
		t12.addExon(new Exon(t12, "e12_1", 1, 2));
		t12.addExon(new Exon(t12, "e12_2", 3, 4));
		t12.addExon(new Exon(t12, "e12_3", 500, 700));
		t12.addExon(new Exon(t12, "e12_4", 800, 900));
		ge1.addTranscript(t12);
		SpliceGraph g1= new SpliceGraph(new Transcript[] {t11,t12});
		g1.init();
		g1.getBubbles();
		
		Gene ge2= new Gene(spec, "ge2");
		ge2.setChromosome("2");
		ge2.setStrand(1);
		Transcript t21= new Transcript(ge2, "t21");
		t21.setStrand(1);
		t21.addExon(new Exon(t21, "e21_1", 1, 2));
		t21.addExon(new Exon(t21, "e21_2", 300, 400));
		t21.addExon(new Exon(t21, "e21_3", 600, 700));
		ge2.addTranscript(t21);
		Transcript t22= new Transcript(ge2, "t22");
		t22.setStrand(1);
		t22.addExon(new Exon(t22, "e22_1", 1, 2));
		t22.addExon(new Exon(t22, "e22_2", 300, 500));
		t22.addExon(new Exon(t22, "e22_3", 600, 700));
		ge2.addTranscript(t22);
		SpliceGraph g2= new SpliceGraph(new Transcript[] {t21,t22});
		g2.init();
		g2.getBubbles();
		
		Mapping[] maps= align(g1, g2);
		for (int i = 0; i < maps.length; i++) {
			System.out.println(maps[i]);
		}
		//GraphView.demo(ge1.getGeneID(), new PFGraph(g1), ge2.getGeneID(), new PFGraph(g2));
	}
	
	SpliceGraph g1, g2;
	
	public GraphAligner(SpliceGraph newG1, SpliceGraph newG2) {
		
	}
	
	public static Mapping[] align(SpliceGraph g1, SpliceGraph g2, int maxAli) {
		
		PriorityQueue q= new PriorityQueue(11, new Mapping.PriorityComparator());
		Comparator compi= new SpliceNode.PositionTypeComparator();
		SpliceNode[] listI= g1.getNodeList();
		Arrays.sort(listI, compi);
		SpliceNode[] listJ= g2.getNodeList();
		Arrays.sort(listJ, compi);
		
		if (listI.length> 0) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(listI[0], null);
			q.add(m);
		}
		if (listJ.length> 0) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(null, listJ[0]);
			q.add(m); 
		}
		if ((listI.length> 0&& listJ.length> 0)&& isAligneable(listI[0], listJ[0])) {
			Mapping m= new Mapping(g1, g2);
			m.addMapping(listI[0], listJ[0]);
			q.add(m);
		}
		Mapping map= (Mapping) q.poll();
		double ulCost= Double.MAX_VALUE;
		Vector optMaps= new Vector();
		//int maxAli= 5;
		while (map!= null&& map.getCost()<= ulCost) {
			
			if (maxAli> 0&& optMaps.size()== maxAli)
				break;
			
				// generate possibilities for next i, next j
			int nextI= -1, nextJ= -1;
			if (map.getMaxI()!= null)
				nextI= Arrays.binarySearch(listI, map.getMaxI(), compi);
			if (map.getMaxJ()!= null)
				nextJ= Arrays.binarySearch(listJ, map.getMaxJ(), compi);

			if (nextI+ 1< listI.length&& map.getMaxRel()!= 2) {
				Mapping m= null;
				try {
					m = (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(listI[nextI+ 1], null);
				q.add(m);
			}
			if (nextJ+ 1< listJ.length&& map.getMaxRel()!= 1) {
				Mapping m= null;
				try {
					m= (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(null, listJ[nextJ+ 1]);
				q.add(m);
			}
			if (((nextI+ 1< listI.length)&& (nextJ+ 1< listJ.length))&& 
					isAligneable(listI[nextI+ 1], listJ[nextJ+ 1])) {
				Mapping m= null;
				try {
					m= (Mapping) map.clone();
				} catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
				m.addMapping(listI[nextI+ 1], listJ[nextJ+ 1]);
				q.add(m);
			}
			
				// update cheapest path
			if ((nextI+ 1>= listI.length)&& (nextJ+ 1>= listJ.length)) {
				if (map.getCost()<= ulCost) {
					if (maxAli< 0&& map.getCost()< ulCost) {
						if (ulCost!= Double.MAX_VALUE)
							System.err.println("assertion failed: ulcost gets cheaper!");
						ulCost= map.getCost();
					}
					optMaps.add(map);
				}
			}
			
			map= (Mapping) q.poll();	// next
		}
		
		return (Mapping[]) gphase.tools.Arrays.toField(optMaps);
	}
	
	static boolean isAligneable(SpliceNode s1, SpliceNode s2) {
		
		if (s1== null|| s2== null)
			return true;		
		
		if ((s1.isDonor()== true&& s2.isDonor()== true)||
				(s1.isAcceptor()== true&&  s2.isAcceptor()== true)||
				(s1.isTSS()== true&&  s2.isTSS()== true)||
				(s1.isTES()== true&&  s2.isTES()== true))
			return true;
		return false;
	}
}
