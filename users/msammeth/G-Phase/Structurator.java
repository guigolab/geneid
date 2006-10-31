

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Properties;

import javax.swing.JFrame;

import org.freehep.util.export.ExportFileType;

import gphase.algo.ASAnalyzer;
import gphase.gui.CopyOfSpliceOSigner;
import gphase.gui.pie.Pie;
import gphase.io.gtf.EncodeWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Graph;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class Structurator {
	
	final static String TABLE_EVEN_COLOR= "#FFFFFF";
	final static String TABLE_ODD_COLOR= "#BFDFDF";
	final static String HEADER_FILE= "header.ins";
	final static String STYLE_FILE= "style.ins";
	final static Color EXSKIP_COL= new Color(58, 83, 164);
	final static Color DBLSKIP_COL= new Color(100, 33, 101);
	final static Color AACC_COL= new Color(237, 34, 36);
	final static Color ADON_COL= new Color(16, 129, 64);
	final static Color IR_COL= new Color(246, 235, 22);
	final static Color OTHERS_COL= new Color(192, 192, 192);
	
	static void include(PrintStream p, String fName) {
		
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			while (buffy.ready()) 
				p.println(buffy.readLine());
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static ASVariation[][] filter(ASVariation[][] vars, int codingC) {
		ASVariation[][] filtClasses= null;
		try {
			String mName= "";
			if (codingC== ASVariation.TYPE_ALL)
				mName= "isTrue";
			else if (codingC== ASVariation.TYPE_CDS)
				mName= "isProteinCoding";
			else if (codingC== ASVariation.TYPE_UTR)
				mName= "isNotAtAllCoding";
			else if (codingC== ASVariation.TYPE_5UTR)
				mName= "isCompletelyIn5UTR";
			else if (codingC== ASVariation.TYPE_3UTR)
				mName= "isCompletelyIn3UTR";
			
			Method m = vars[0][0].getClass().getMethod(mName, null);
			filtClasses= (ASVariation[][]) Arrays.filter(vars, m);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return filtClasses;
	}
	
	static void outputPic(ASVariation as) {
	      
		  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
	      Properties props= null;
	      String creator= null;
	      Component component= new CopyOfSpliceOSigner(as);
	      component.setSize(component.getPreferredSize());
	      File f= new File("test.gif");
	      try {
			t.exportToFile(f,component,component.getParent(),props,creator);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
//	      props.put(SAVE_AS_FILE,file.getText());
//	      props.put(SAVE_AS_TYPE,currentType().getFileFilter().getDescription());
	}
	
	public static void main(String[] args) {

		//
		long t0= System.currentTimeMillis();
		boolean codingTranscripts= false;
		boolean noncodTranscripts= false;
		int filterCode= ASMultiVariation.FILTER_HIERARCHICALLY;
		int codingCode= ASVariation.TYPE_ALL;
		String fName= null;
		boolean html= false;
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-codingTranscripts"))
				codingTranscripts= true;
			if (args[i].equalsIgnoreCase("-noncodTranscripts"))
				noncodTranscripts= true;

			if (args[i].equalsIgnoreCase("-noFilter"))
				filterCode= ASMultiVariation.FILTER_NONE;
			if (args[i].equalsIgnoreCase("-structFilter"))
				filterCode= ASMultiVariation.FILTER_STRUCTURALLY;
			
			if (args[i].equalsIgnoreCase("-CDS"))
				codingCode= ASVariation.TYPE_CDS;
			if (args[i].equalsIgnoreCase("-UTR"))
				codingCode= ASVariation.TYPE_UTR;
			if (args[i].equalsIgnoreCase("-5UTR"))
				codingCode= ASVariation.TYPE_5UTR;
			if (args[i].equalsIgnoreCase("-3UTR"))
				codingCode= ASVariation.TYPE_3UTR;

			if (args[i].equalsIgnoreCase("-html"))
				html= true;
			
			if (!args[i].startsWith("-")|| args[i].contains(File.separator))
				fName= args[i];
		}
		
		if (html&& fName== null) {
			System.err.println("Cannot guess output directory for html files");
			System.exit(-1);
		}
		
			// read
		EncodeWrapper enc;
		if (fName== null) {
			enc= new EncodeWrapper(System.in);
		} else {
			enc= new EncodeWrapper(fName);
		}
		Graph g= enc.getGraph(false);
		if (g== null) {
			System.exit(-1);
		}
		if (codingTranscripts)
			g.filterNonCodingTranscripts();
		if (noncodTranscripts)
			g.filterCodingTranscripts();
		
			// get vars
		ASVariation[][] vars= g.getASVariations(filterCode);
		if (vars!= null&& vars.length!= 0) {
			vars= filter(vars, codingCode);
			vars= (ASVariation[][]) Arrays.sort2DFieldRev(vars);
		}
		if (html)
			writeHTML(vars, fName);
		else
			ASAnalyzer.outputVariations(vars, false, false, System.out);
		
		long t1= System.currentTimeMillis();
		System.out.println(t1-t0);
	}
	
	static void writeHTML(ASVariation[][] vars, String fName) {
		
		fName= fName.substring(0, fName.lastIndexOf(File.separator));
		
		PrintStream p;
		try {
			writePiePicture(vars, fName);
			writePictures(vars, fName);
			
			p= new PrintStream(fName+ File.separator+ "statistics.html");
			writeHTMLStats(vars, p);
			p.flush(); p.close();
			
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				String varStr= vars[i][0].toString();
				String varFStr= varStr.replace('/', '_')+ ".html";
				
				p= new PrintStream(fName+ File.separator+ varFStr);
				writeHTMLNumbers(vars[i], p);
				p.flush(); p.close();
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	

	static void writePiePicture(ASVariation[][] vars, String path) {

			// init pie
		Pie pie= new Pie();
		pie.setSize(new Dimension(400,400));
		pie.init();
		
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			if (vars[i][0].toString().equals("( // 1=2^)"))
				pie.AddPieSlice(vars[i].length, "exon skipping", EXSKIP_COL);
			else if (vars[i][0].toString().equals("( // 1=2^3=4^)"))
				pie.AddPieSlice(vars[i].length, "double skipping", DBLSKIP_COL);
			else if (vars[i][0].toString().equals("(1= // 2=)"))
				pie.AddPieSlice(vars[i].length, "alt acceptor", AACC_COL);
			else if (vars[i][0].toString().equals("(1^ // 2^)"))
				pie.AddPieSlice(vars[i].length, "alt donor", ADON_COL);
			else if (vars[i][0].toString().equals("( // 1^2=)"))
				pie.AddPieSlice(vars[i].length, "intron retention", IR_COL);
			else
				sum+= vars[i].length;
		}
		if (vars!= null)
			pie.AddPieSlice(sum, "other events", OTHERS_COL);
		
		
		pie.SetTitle("Distribution");
		pie.setShow_values_on_slices(1);
		
		pie.setSize(pie.getPreferredSize());
//		JFrame frame= new JFrame();
//		frame.getContentPane().add(pie);
//		frame.setVisible(true);
		try {
		      File f= new File(path+File.separator+"distribution.png");
			  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("png").get(0);
		      t.exportToFile(f,pie,pie.getParent(),null,null);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	static void writePictures(ASVariation[][] vars, String path) {
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			try {
			      String varStr= vars[i][0].toString();
			      String varFStr= varStr.replace('/', '_')+ ".gif";
			      File f= new File(path+File.separator+varFStr);
				  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			      Component component= new CopyOfSpliceOSigner(vars[i][0]);
			      component.setSize(component.getPreferredSize());
			      t.exportToFile(f,component,component.getParent(),null,null);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	static void writeHTMLStats(ASVariation[][] vars, PrintStream p) {
		int spacer= 80;
		p.println("<HTML>");
		p.println("<HEAD>");
		p.println("<TITLE>Statistics</TITLE>");
		include(p, STYLE_FILE);
		p.println("</HEAD>");
		p.println("<BODY>");

		include(p, HEADER_FILE);
		if (vars!= null)
			p.println("<div class=\"title\"><h1>STATISTICS</h1></div><br />");
		else 
			p.println("<div class=\"title\"><h1>NO AS EVENTS FOUND</h1></div><br />");
			
		p.println("<img src=\"distribution.png\">");
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			sum+= vars[i].length;
		if (vars!= null) {
			p.println("<TABLE>");
			p.println("<TR>");
			p.println("<TH align=\"right\" valign=\"middle\" >Rank</TH>");
			p.println("<TH align=\"right\" valign=\"middle\" >Count</TH>");
			p.println("<TH align=\"right\" valign=\"middle\" >Fraction</TH>");
			p.println("<TH align=\"left\" valign=\"middle\" >Structure</TH>");
			p.println("</TR>");
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				String colStr= "";
				if (i%2 == 0)
					colStr= TABLE_EVEN_COLOR;
				else
					colStr= TABLE_ODD_COLOR;
				p.println("<TR bgcolor=\""+colStr+"\">");
				p.println("\t<TD align=\"right\" valign=\"middle\" width=\"120\">#"+ (i+1)+ "</TD>");
				p.println("\t<TD align=\"right\" valign=\"middle\" width=\"120\">"+ vars[i].length+ "</TD>");
				float perc= 100f* vars[i].length/ sum; 
				String percStr= Float.toString(perc);
				int cutoff= Math.min(percStr.indexOf('.')+ 3, percStr.length());
				percStr= percStr.substring(0, cutoff)+ "%";
				p.println("\t<TD align=\"right\" valign=\"center\" width=\"120\">"+ percStr+ "</TD>");
				String varStr= vars[i][0].toString();
				String varFStr= varStr.replace('/', '_');
				p.println("\t<TD align=\"left\" valign=\"center\">"
						+ "<img src=\""+ varFStr+".gif\"><br>"
						+ "<a href=\""+ varFStr+".html\">"
						+ varStr+ "</a></TD>");
				p.println("</TR>");
			}
			p.println("</TABLE>");
		}
		p.println("</div><!-- closing main -->");
		p.println("</BODY>");
		p.println("</HTML>");
	}

	static void writeHTMLNumbers(ASVariation[] vars, PrintStream p) {
		p.println("<HTML>");
		p.println("<HEAD>");
		p.println("<TITLE>Events</TITLE>");
		include(p, STYLE_FILE);
		p.println("</HEAD>");
		p.println("<BODY>");
		
		include(p, HEADER_FILE);
		p.println("<div class=\"title\"><h1>EVENTS</h1></div><br />");
		p.println("<a href=\"statistics.html\">Back</a>");
		p.println("<TABLE border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
		for (int i = 0; i < vars.length; i++) {
			String colStr= "";
			if (i%2 == 0)
				colStr= TABLE_EVEN_COLOR;
			else
				colStr= TABLE_ODD_COLOR;
			p.println("<TR>");
			
				// event
			SpliceSite[] su= vars[i].getSpliceUniverse();
			String ucscEvent= "http://genome.ucsc.edu/cgi-bin/hgTracks?org=Homo_sapiens&db=hg17&position=chr"
				+ vars[i].getGene().getChromosome()+ ":"+Math.abs(su[0].getPos())+"-"+ Math.abs(su[su.length- 1].getPos());
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"middle\" width=\"120\"><a href=\""
					+ ucscEvent+"\">"+ (i+1)+ "</a></TD>");
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"top\" width=\"120\">"
					+ vars[i].getTranscript1().getTranscriptID()+ "<br>"
					+ vars[i].getTranscript2().getTranscriptID()+ "</TD>");

				// splice chains
			String pos1Str= "";
			for (int j = 0; j < vars[i].getSpliceChain1().length; j++)
				if (vars[i].getSpliceChain1()[j].getPos()> 0)
					pos1Str+= vars[i].getSpliceChain1()[j]+ ", ";
				else
					pos1Str+= Math.abs(vars[i].getSpliceChain1()[j].getPos())+ ", ";
			if (pos1Str.length()> 2)
				pos1Str= pos1Str.substring(0, pos1Str.length()- 2);
			String pos2Str= "";
			for (int j = 0; j < vars[i].getSpliceChain2().length; j++) 
				if (vars[i].getSpliceChain2()[j].getPos()> 0)
					pos2Str+= vars[i].getSpliceChain2()[j]+ ", ";
				else
					pos2Str+= Math.abs(vars[i].getSpliceChain2()[j].getPos())+ ", ";
			if (pos2Str.length()> 2)
				pos2Str= pos2Str.substring(0, pos2Str.length()- 2);
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"top\" width=\"120\">"
					+ vars[i].getGene().getChromosome()+ ":"+ pos1Str+ "<br>"
					+ vars[i].getGene().getChromosome()+ ":"+ pos2Str+ "</TD>");

			String varStr= vars[i].toString();
			String varFStr= varStr.replace('/', '_');
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"top\" width=\"120\">"
					+ "<img src=\""+ varFStr+".gif\"><br>"
					+ varStr+ "</TD>");
			
			p.println("</TR>");
		}
		p.println("</TABLE>");
		p.println("</div><!-- closing main -->");
		p.println("</BODY>");
		p.println("</HTML>");
	}
}
