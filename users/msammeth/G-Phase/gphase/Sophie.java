package gphase;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import gphase.io.gtf.GTFChrReader;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.tools.Arrays;

public class Sophie {

	public static void main(String[] args) {
		_070716_selectedExons();
	}
	
	static void _070716_selectedExons() {
			// the exons
		HashMap<DirectedRegion, Gene> exonMap= new HashMap<DirectedRegion, Gene>();
		Vector<DirectedRegion> v= new Vector<DirectedRegion>();
		v.add(new DirectedRegion(33235997,33236077,1,"chr10"));
		v.add(new DirectedRegion(33235997,33236077,-1,"chr10"));
		v.add(new DirectedRegion(143666573,143666649,1,"chr1"));
		v.add(new DirectedRegion(143666573,143666649,-1,"chr1"));
		v.add(new DirectedRegion(143663558,143663733,1,"chr1"));
		v.add(new DirectedRegion(143663558,143663733,-1,"chr1"));
		v.add(new DirectedRegion(143692139,143692331,1,"chr1"));		
		v.add(new DirectedRegion(143692139,143692331,-1,"chr1"));		
		v.add(new DirectedRegion(199596763,199596917,1,"chr1"));
		v.add(new DirectedRegion(199596763,199596917,-1,"chr1"));
		v.add(new DirectedRegion(199596204,199596395,1,"chr1"));
		v.add(new DirectedRegion(199596204,199596395,-1,"chr1"));
		v.add(new DirectedRegion(61143316,61143394,1,"chr15"));
		v.add(new DirectedRegion(61143316,61143394,-1,"chr15"));
		v.add(new DirectedRegion(35672923,35673238,1,"chr9"));
		v.add(new DirectedRegion(35672923,35673238,-1,"chr9"));
		v.add(new DirectedRegion(36662152,36662237,1,"chr2"));
		v.add(new DirectedRegion(36662152,36662237,-1,"chr2"));
		v.add(new DirectedRegion(215953779,215954048,1,"chr2"));
		v.add(new DirectedRegion(215953779,215954048,-1,"chr2"));
		v.add(new DirectedRegion(234964356,234964441,1,"chr1"));
		v.add(new DirectedRegion(234964356,234964441,-1,"chr1"));
		v.add(new DirectedRegion(61149122,61149889,1,"chr15"));
		v.add(new DirectedRegion(61149122,61149889,-1,"chr15"));
		v.add(new DirectedRegion(35672794,35672901,1,"chr9"));
		v.add(new DirectedRegion(35672794,35672901,-1,"chr9"));
		v.add(new DirectedRegion(1914020,1914060,1,"chr11"));
		v.add(new DirectedRegion(1914020,1914060,-1,"chr11"));
		v.add(new DirectedRegion(1914020,1914060,1,"chr11"));
		v.add(new DirectedRegion(88116912,88117104,-1,"chr5"));
		v.add(new DirectedRegion(100584653,100584706,1,"chr12"));
		v.add(new DirectedRegion(100584653,100584706,-1,"chr12"));
		v.add(new DirectedRegion(179362950,179363087,1,"chr2"));
		v.add(new DirectedRegion(179362950,179363087,-1,"chr2"));
		v.add(new DirectedRegion(182175795,182175886,1,"chr3"));
		v.add(new DirectedRegion(182175795,182175886,-1,"chr3"));
		v.add(new DirectedRegion(182171557,182171637,1,"chr3"));
		v.add(new DirectedRegion(182171557,182171637,-1,"chr3"));
		v.add(new DirectedRegion(29875817,29876019,1,"chr20"));
		v.add(new DirectedRegion(29875817,29876019,-1,"chr20"));
		v.add(new DirectedRegion(22721563,22721698,1,"chr15"));
		v.add(new DirectedRegion(22721563,22721698,-1,"chr15"));
		
		
			// hg18
		GTFChrReader reader= new GTFChrReader("/home/msammeth/annotations/human_hg18_RefSeqGenes_fromUCSC070716_mRNAs_fromUCSC070716_CDS.gtf");
		reader.setReadGene(true);
		reader.setChromosomeWise(true);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] g= reader.getGenes();
		DirectedRegion[] regs= (DirectedRegion[]) Arrays.toField(v);
		while (g!= null&& regs!= null) {
			
			for (int i = 0; i < g.length; i++) {
				for (int j = 0; j < regs.length; j++) {
					if (g[i].overlaps(regs[j])) {
						if (exonMap.get(regs[j])!= null)
							System.err.println("Multiple genes found.");
						exonMap.put(regs[j], g[i]);
					}
				}
			}
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			g= reader.getGenes();
		}
		
			// output
		PrintStream p= null;
		try {
			p= new PrintStream("selectedExons.html");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		p.println("<html><head></head><body>");
		p.println("<h1>Selected Exons to Micha</h1>");
		p.println("<ul>");
		for (int i = 0; i < regs.length; i++) {
			String ucscBase= AStaLaVista.UCSC_GB_CGI + AStaLaVista.UCSC_LOAD_DEFAULTS;
			ucscBase+= ";"+"org=human;db=hg18;";
			if (exonMap.get(regs[i])!= null)
				ucscBase+= AStaLaVista.UCSC_LOAD_CT+ "http://genome.imim.es/~msammeth/customTracks/hg18/domains_Pfam"+ "/" 
					+ regs[i].getChromosome()+"_"+exonMap.get(regs[i]).getNameTranscript().getTranscriptID()+ ".bed;";
			ucscBase+= "position="+regs[i].getChromosome()+":"+Math.abs(regs[i].getStart())+"-"+Math.abs(regs[i].getEnd());
			
			p.println("<li><a href=\""+ucscBase+"\">"+regs[i].toUCSCString()+(regs[i].isForward()?" forward":" reverse")+"</a></li>");
		}
		p.flush(); p.close();
		p.println("</ul>");
		p.println("</body></html>");
	}
}
