package gphase;



import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PipedInputStream;
import java.io.PipedOutputStream;
import java.io.PrintStream;
import java.lang.reflect.Method;
import java.util.Properties;
import java.util.Vector;

import javax.swing.JFrame;

import org.freehep.util.export.ExportFileType;

import gphase.algo.ASAnalyzer;
import gphase.ext.DevNullReaderThread;
import gphase.gui.CopyOfSpliceOSigner;
import gphase.gui.pie.Pie;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class Structurator {
	final static String ANNOTATION_DIR= "annotation";
	final static String UCSC_URL= "http://genome.ucsc.edu/cgi-bin/hgTracks?";
	public static final String[] SP_UCSC_CGI_STRINGS= new String[] {
		"knownGene=dense;encodeRegions=dense;encodeGencodeGeneOct05=pack;"
		+ "gap=hide;stsMap=hide;mgcGenes=hide;exoniphy=hide;exonWalk=hide;multiz17way=hide;snp125=hide;",	// human
		"",	// chimp
		"knownGene=dense;stsMap=hide;mgcGenes=hide;snp126=hide;multiz17way=hide;uniGene_3=hide;",	// Mouse
		"stsMapRat=hide;mgcGenes=hide;multiz9way=hide;netXenTro2=hide;netHg18=hide;snp125=hide;", // Rat
		"multiz4way=hide;netHg18=hide;blastHg18KG=hide;",		// &db=canFam2
		"all_mrna=pack;",		//Cow
		"multiz7way=hide;",	// Opossum
		
		"",	//Chicken
		"",	//X.+tropicalis 
		"",	//Zebrafish
		"blastHg17KG=hide;cpgIsland=hide;blatHg16=hide;",			//Fugu
		"gaze=pack;netSelf=hide;",	//Tetraodon
		
		"flyBaseGene=dense;flyBaseNoncoding=pack;"
		+ "multiz17way=hide;blastHg17KG=hide;",	// Drosophila
		"",	//A.+gambiae
		"modelRefGene=dense;brhInparalog=hide;blastDm2FB=hide;netDm2=hide;",  //A.+mellifera
		
		"sangerGene=pack;sangerGenefinder=hide;c_briggsae_pwMaf=hide;wabaCbr=hide;axtNetCb1=hide;",	//C.+elegans
		"sgdGene=pack;sgdOther=hide;transRegCode=hide;esRegGeneToMotif=hide;multizYeast=hide;"	//S.+cerevisiae
};
	final static String UCSC_STANDARD_PAR=
			// activate
		"ruler=dense;refGene=pack;mrna=pack;"
			// deactivate
//		+"gap=hide;nscanGene=hide;intronEst=hide;rmsk=hide;"
			// all species projection on
		+"knownGene=pack;flyBaseGene=pack;modelRefGene=pack;sangerGene=pack;sgdGene=pack;" +
				"gaze=pack;flyBaseNoncoding=pack;all_mrna=pack;" +
				"encodeRegions=dense;encodeGencodeGeneOct05=pack;"
			// all species projection off
//		+"gap=hide;cpgIsland=hide;" +
//				"exoniphy=hide;exonWalk=hide;" +
//				"stsMap=hide;stsMapRat=hide;sgdOther=hide;transRegCode=hide;esRegGeneToMotif=hide;" +
//				"mgcGenes=hide;uniGene_3=hide;sangerGenefinder=hide;" +
//				"snp125=hide;snp126=hide;" +
//				"multiz17way=hide;multiz9way=hide;multiz4way=hide;multiz7way=hide;multizYeast=hide;" +
//				"netXenTro2=hide;netHg18=hide;netSelf=hide;brhInparalog=hide;blastDm2FB=hide;netDm2=hide;c_briggsae_pwMaf=hide;wabaCbr=hide;axtNetCb1=hide;" +
//				"blastHg18KG=hide;blastHg17KG=hide;blatHg16=hide;"
		+"";
	
	final static String TABLE_EVEN_COLOR= "#FFFFFF";
	final static String TABLE_ODD_COLOR= "#BFDFDF";
	final static String TABLE_HEADER_COLOR= "#00A0A0"; //"#008080";
	final static String HEADER_FILE= "header.ins";
	final static String TRAILER_FILE= "trailer.ins";
	final static String STYLE_FILE= "style.ins";
	final static Color EXSKIP_COL= new Color(58, 83, 164);
	final static Color DBLSKIP_COL= new Color(100, 33, 101);
	final static Color AACC_COL= new Color(237, 34, 36);
	final static Color ADON_COL= new Color(16, 129, 64);
	final static Color IR_COL= new Color(246, 235, 22);
	final static Color OTHERS_COL= new Color(192, 192, 192);
	final static int UCSC_FLANK= 50;
	static String speStr= null;	//Species.SP_UCSC_CGI_STRINGS[0];
	static String annStr= Species.SP_UCSC_CGI_STRINGS[0];
	static boolean outputASTA= false;
	static boolean outputGTF= false;
	
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
		boolean nmd= false;
		boolean filtGTAG= false;
		for (int i = 0; i < args.length; i++) {

				// let this first check else all double arg flags have to continue
			if (!args[i].startsWith("-")|| args[i].contains(File.separator))
				fName= args[i];

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

			if (args[i].equalsIgnoreCase("-html")) {
				html= true;
			}
			if (args[i].equalsIgnoreCase("-output")) {
				if (i+1>= args.length) {
					System.err.println("provide output format");
					break;
				}
				String oFormat= args[i+1].toUpperCase();
				if (args[i+1].equalsIgnoreCase("asta"))  
					outputASTA= true;
				else if (args[i+1].equalsIgnoreCase("gtf"))
					outputGTF= true;
				else
					System.err.println("Unknown output format "+args[i+1]);
				++i;
			}
			
			if (args[i].equalsIgnoreCase("-filtGTAG"))
				filtGTAG= true;
			if (args[i].equalsIgnoreCase("-nonmd"))
				nmd= false;
			if (args[i].equalsIgnoreCase("-nmd"))
				nmd= true;
			
			if (args[i].equalsIgnoreCase("-species")) {
				if (i+1>= args.length) {
					System.err.println("provide species name");
					break;
				}
				String[] spe= Species.SP_NAMES_COMMON;
				speStr= null;
				for (int j = 0; j < spe.length; j++) 
					if (spe[j].equalsIgnoreCase(args[i+1])) {
						speStr= Species.SP_UCSC_CGI_STRINGS[j];
						++i;
						break;
					}
				if (speStr== null) {
					System.err.println("Unknown species: "+args[i+1]);
					System.err.print("Enter one of the following: ");
					for (int j = 0; j < spe.length; j++) 
						System.err.print(spe[j]+" ");
				}
				
			}
			
			if (args[i].equalsIgnoreCase("-annotation")) {
				if (i+1>= args.length) {
					System.err.println("provide annotation name");
					break;
				}
				
				
			}

			
			if (args[i].equalsIgnoreCase("-genome")) {
				if (i+1>= args.length) {
					System.err.println("provide genome ID");
					break;
				}
				
				int p= args[i+1].indexOf("_");	// i would say mouse_mm8 ok?
				if (p>= 0) {
					speStr= args[i+1].substring(0, p);
					annStr= args[i+1].substring(p+1, args[i+1].length());
				} else 
					annStr= args[i+1];	// not tested, check for mapping to UCSC string
				++i;
			}
		}
		
		if (html) {
			if (fName== null) {
				System.err.println("Cannot guess output directory for html files");
				System.exit(-1);
			} else {
				PrintStream sysOut= System.out;
				PipedInputStream pin= new PipedInputStream();
				DevNullReaderThread nil= new DevNullReaderThread(pin);
				nil.start();
				try {
					PipedOutputStream pout= new PipedOutputStream(pin);
					PrintStream p= new PrintStream(pout);
					System.setOut(p);
					System.setErr(p);	// kill stderr, sylvain stops otherwise
				} catch (IOException e) {	
					e.printStackTrace();
				}
			}
		}
		
		// force nonmd
//		if (nmd) {
//			nmd= false;
//		}
			
			// read
		EncodeWrapper enc;
		if (fName== null) {
			enc= new EncodeWrapper(System.in);
		} else {
			enc= new EncodeWrapper(fName);
		}
		//enc.setSilent(true);
		
		if (filtGTAG) {
			enc.setSpeName(speStr);
			enc.setGenomeVer(annStr);
		}
		boolean encode= false;
//		if (fName.toUpperCase().contains("ENCODE"))
//			encode= true;
		Graph g= enc.getGraph(encode);
		enc= null;
		System.gc();
		Thread.currentThread().yield();
		
		if (g== null) {
			System.exit(-1);
		}
		if (codingTranscripts)
			g.filterNonCodingTranscripts();
		if (noncodTranscripts)
			g.filterCodingTranscripts();
		if (nmd)
			g.filterNMDTranscripts();
		
			// get vars 
		ASVariation[][] vars= g.getASVariations(filterCode);
		if (filtGTAG) {
			Vector v= new Vector();
			int cntFiltEv= 0;
			int cntAllEv= 0;
			for (int i = 0; vars!= null&& i < vars.length; i++) { 
				ASVariation[] filtClas= ASMultiVariation.filterNonGTAG(vars[i]);
				if (filtClas== null)
					cntFiltEv+= vars[i].length;
				else
					cntFiltEv+= vars[i].length- filtClas.length;
				cntAllEv+= vars[i].length;
				if (filtClas!= null&& filtClas.length> 0)
					v.add(filtClas);
			}
			vars= (ASVariation[][]) Arrays.toField(v);
			String percStr= Float.toString((cntFiltEv* 100f)/ cntAllEv);
			percStr= percStr.substring(0, percStr.indexOf('.')+ 2);
			System.out.println("filtered out "+cntFiltEv+" ("+percStr+"%) involving nonGT/AG introns.");
		}
		if (vars!= null&& vars.length!= 0) {
			vars= filter(vars, codingCode);
			vars= (ASVariation[][]) Arrays.sort2DFieldRev(vars);
		}

		int pos= fName.lastIndexOf(File.separator);
		String baseDir= "";
		if (pos>= 0)
			baseDir= fName.substring(0, pos+ 1);
		
		if (outputASTA) {
			PrintStream p;
			try {
				p= new PrintStream(baseDir+ "landscape.asta");
				writeASTA(vars, p);
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		if (outputGTF) {
			PrintStream p;
			writeGTF(vars, baseDir);
		}
		
		if (html)
			writeHTML(vars, baseDir);
		else
			ASAnalyzer.outputVariations(vars, false, false, System.out);
		
		long t1= System.currentTimeMillis();
		System.out.println("time: "+ (t1-t0)+"[msec]");
	}
	
	public static void writeASTA(ASVariation[][] vars, PrintStream p) {		
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			for (int j = 0; j < vars[i].length; j++) 
				p.println(vars[i][j].toStringASTA());
	}

	public static void writeGTF(ASVariation[][] vars, String baseDir) {
		Vector v= new Vector();
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			for (int j = 0; j < vars[i].length; j++) {
				GTFObject o= new GTFObject();
				o.setSeqname(vars[i][j].getGene().getChromosome());
				o.setSource("astalavista");
				o.setFeature("as_event");
				SpliceSite[] su= vars[i][j].getSpliceUniverse();
				int start= Math.abs(su[0].getPos());
				int end= Math.abs(su[su.length- 1].getPos());
				if (start> end) {
					int h= start;
					start= end;
					end= h;
				}
				o.setStart(start);
				o.setEnd(end);
				o.setScore(".");
				o.setStrand(vars[i][j].getGene().getStrand());
				// frame
				o.addAttribute("as_code", vars[i][j].toString());
				o.addAttribute("transcript1_id", vars[i][j].getTranscript1().toString());
				String sc= "";
				for (int k = 0; k < vars[i][j].getSpliceChain1().length; k++) 
					sc+= Math.abs(vars[i][j].getSpliceChain1()[k].getPos())+",";
				if (sc.length()> 0)
					sc= sc.substring(0, sc.length()- 1);
				o.addAttribute("splice_chain1", sc);
				
				o.addAttribute("transcript2_id", vars[i][j].getTranscript2().toString());
				sc= "";
				for (int k = 0; k < vars[i][j].getSpliceChain2().length; k++) 
					sc+= Math.abs(vars[i][j].getSpliceChain2()[k].getPos())+",";
				if (sc.length()> 0)
					sc= sc.substring(0, sc.length()- 1);
				o.addAttribute("splice_chain2", sc);
				
				v.add(o);
			}
		
		GTFWrapper gtf= new GTFWrapper(baseDir+ File.separator+ "landscape.gtf");
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v));
		gtf.setSortAttributes(new String[] {"as_code", "transcript1_id", "splice_chain1", "transcript2_id", "splice_chain2"});
		try {
			gtf.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void writeHTML(ASVariation[][] vars, String fName) {
		
		PrintStream p;
		try {
			writePiePicture(vars, fName);
			writePictures(vars, fName);
			
			p= new PrintStream(fName+ "landscape.html");
			writeLandscapeHTML(vars, p);
			p.flush(); p.close();
			
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				String varStr= vars[i][0].toString();
				String varFStr= convertFName(varStr)+ ".html";	//varStr.replace(' ', '_')+ ;
				
				p= new PrintStream(fName+ varFStr);
				writeEventsHTML(vars[i], p);
				p.flush(); p.close();
			}
			
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}
	

	static void writePiePicture(ASVariation[][] vars, String path) {

			// init pie
		Pie pie= new Pie();
		pie.setSize(new Dimension(300,300));
		pie.init();
		
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			if (vars[i][0].toString().equals("1-2^ , 0"))
				pie.AddPieSlice(vars[i].length, "exon skipping", EXSKIP_COL);
//			else if (vars[i][0].toString().equals("1-2^3-4^ , 0"))
//				pie.AddPieSlice(vars[i].length, "double skipping", DBLSKIP_COL);
			else if (vars[i][0].toString().equals("1- , 2-"))
				pie.AddPieSlice(vars[i].length, "alt acceptor", AACC_COL);
			else if (vars[i][0].toString().equals("1^ , 2^"))
				pie.AddPieSlice(vars[i].length, "alt donor", ADON_COL);
			else if (vars[i][0].toString().equals("1^2- , 0"))
				pie.AddPieSlice(vars[i].length, "intron retention", IR_COL);
			else
				sum+= vars[i].length;
		}
		if (vars!= null)
			pie.AddPieSlice(sum, "other events", OTHERS_COL);
		
		
		pie.SetTitle("Overview");
		pie.setShow_values_on_slices(1);
		
		pie.setSize(pie.getPreferredSize());
//		JFrame frame= new JFrame();
//		frame.getContentPane().add(pie);
//		frame.setVisible(true);
		try {
		      File f= new File(path+"distribution.png");
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
			      String varFStr= convertFName(varStr)+ ".gif";		//varStr.replace(' ', '_');
			      File f= new File(path+varFStr);
				  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			      Component component= new CopyOfSpliceOSigner(vars[i][0]);
			      component.setSize(component.getPreferredSize());
			      t.exportToFile(f,component,component.getParent(),null,null);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	static void headline(PrintStream p, String headline, String linkBack) {
//		<TD class="section" border=0 cellpadding=0 cellspacing=0 width=20% ALIGN=RIGHT>
//		<a href="index.html#TOP"
//		 onMouseover="window.status='TOP:INDEX';">
//		<IMG class="pnt" SRC="/g_icons/top.gif" HEIGHT=15 WIDTH=15 BORDER=0></a>
		p.println("<TABLE border=0 cellpadding=0 cellspacing=0 width=\"100%\">");
		p.println("<TR border=0 cellpadding=0 cellspacing=0>");
		if (linkBack!= null)
			p.println("<TD class=\"section\" border=0 cellpadding=0 cellspacing=0 width=\"20%\"><a href=\""+linkBack+"\"><IMG class=\"pnt\" SRC=\"/g_icons/top.gif\" HEIGHT=15 WIDTH=15 BORDER=0></a></TD>");
		p.println("<TD class=\"section\" align=\"center\">");
		p.println("<CENTER><FONT size=6 class=\"tgen\">"+headline+"</FONT></CENTER></TD>");
		p.println("</TD></TR></TABLE>");
	}
	
	static void writeLandscapeHTML(ASVariation[][] vars, PrintStream p) {
		int spacer= 80;
//		p.println("<HTML>");
//		p.println("<HEAD>");
//		p.println("<TITLE>AS Landscape</TITLE>");
//		include(p, STYLE_FILE);
//		p.println("</HEAD>");
//		p.println("<BODY>");
	
		include(p, HEADER_FILE);		

		if (vars== null)
			headline(p, "NO AS EVENTS FOUND", null);
		
		else {
			headline(p, "AS Landscape", null);
			p.println("<DIV class=\"userspace\" align=\"center\"><CENTER><br><img src=\"distribution.png\"><br><br>");
				
			if (outputASTA)
				p.println("<IMG src= \"http://genome.imim.es/astalavista/pics/dl_arrow.jpg\">"+
						"<A HREF=\"landscape.asta\">Download output (ASTA format selected)</A><br><br></CENTER></DIV>");
			if (outputGTF)
				 p.println("<IMG src= \"http://genome.imim.es/astalavista/pics/dl_arrow.jpg\">"+
						"<A HREF=\"landscape.gtf\">Download output (GTF format selected)</A><br><br></CENTER></DIV>");
		}
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			sum+= vars[i].length;
		if (vars!= null) {
			p.println("<DIV class=\"userspace\" align=\"center\">");
			p.println("<TABLE bgcolor=\"#FFFFFF\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
			p.println("<TR>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\" valign=\"middle\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Rank</b></FONT></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\" valign=\"middle\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Proportion</b></FONT></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Count</b></FONT></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Details</b></FONT></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"left\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Intron-Exon<br>Structure</b></FONT></TH>");
			p.println("</TR>");
			
			int outCnt= 0, innCnt= 0;
			int lastEvCnt= Integer.MAX_VALUE;
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				if (vars[i].length< lastEvCnt) {
					lastEvCnt= vars[i].length;
					++outCnt; innCnt= 1;
				} else 
					++innCnt;
				
					
				String colStr= "";
				if (i%2 == 0)
					colStr= TABLE_EVEN_COLOR;
				else
					colStr= TABLE_ODD_COLOR;
				p.println("<TR bgcolor=\""+colStr+"\">");
				p.println("\t<TD align=\"center\" valign=\"middle\" width=\"60\">"  
						+ "<FONT size=\"4\">"+outCnt+"</FONT>.<FONT size=\"2\">"+innCnt+ "</FONT></TD>");
				p.println("\t<TD></TD>"); 
				float perc= 100f* vars[i].length/ sum; 
				String percStr= Float.toString(perc);
				int cutoff= Math.min(percStr.indexOf('.')+ 3, percStr.length());
				percStr= percStr.substring(0, cutoff)+ "%";
				p.println("\t<TD align=\"right\" valign=\"center\" width=\"60\">"+ percStr+ "</TD>");
				p.println("\t<TD></TD>"); 
				p.println("\t<TD align=\"right\" valign=\"middle\" width=\"60\">"+ vars[i].length+ "</TD>");
				p.println("\t<TD></TD>"); 
				String varStr= vars[i][0].toString();
				String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
				varStr= varStr.replaceAll("\\s", "");
				p.println("\t<TD align=\"center\" valign=\"center\">"
						+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""+ varFStr+".html\"></FONT>Show>></a>");
				p.println("\t<TD></TD>"); 
				p.println("\t<TD align=\"left\" valign=\"center\">"
						//+ "<a href=\""+ varFStr+".html\">"
						+ "<img src=\""+ varFStr+".gif\"><br>"
						// </a>
						+ "<FONT size=\"2\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ varStr+ "</b></FONT>"
						+"</TD>");
				p.println("</TR>");
			}
			p.println("</TABLE></DIV>");
		}
		p.println("</div><!-- closing main -->");
		
		include(p, TRAILER_FILE);
	}

	static void _style_writeHTMLStats(ASVariation[][] vars, PrintStream p) {
		int spacer= 80;
		p.println("<HTML>");
		p.println("<HEAD>");
		p.println("<TITLE>AS Landscape</TITLE>");
		include(p, STYLE_FILE);
		p.println("</HEAD>");
		p.println("<BODY>");

		include(p, HEADER_FILE);
		if (vars!= null)
			p.println("<div class=\"title\"><h1>AS Landscape</h1></div><br />");
		else 
			p.println("<div class=\"title\"><h1>NO AS EVENTS FOUND</h1></div><br />");
			
		p.println("<img src=\"distribution.png\">");
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			sum+= vars[i].length;
		if (vars!= null) {
			p.println("<TABLE>");
			p.println("<TR>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"left\"><b>Rank</b></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><b>Count</b></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\"><b>Fraction</b></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Details</b></TH>");
			p.println("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Structure</b></TH>");
			p.println("</TR>");
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				String colStr= "";
				if (i%2 == 0)
					colStr= TABLE_EVEN_COLOR;
				else
					colStr= TABLE_ODD_COLOR;
				p.println("<TR bgcolor=\""+colStr+"\">");
				p.println("\t<TD align=\"left\" valign=\"middle\" width=\"60\">#"+ (i+1)+ "</TD>");
				p.println("\t<TD align=\"right\" valign=\"middle\" width=\"60\">"+ vars[i].length+ "</TD>");
				float perc= 100f* vars[i].length/ sum; 
				String percStr= Float.toString(perc);
				int cutoff= Math.min(percStr.indexOf('.')+ 3, percStr.length());
				percStr= percStr.substring(0, cutoff)+ "%";
				p.println("\t<TD align=\"right\" valign=\"center\" width=\"60\">"+ percStr+ "</TD>");
				String varStr= vars[i][0].toString();
				String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
				p.println("\t<TD align=\"center\" valign=\"center\">"
						+ "<a href=\""+ varFStr+".html\">Show>></a>");
				p.println("\t<TD align=\"left\" valign=\"center\">"
						+ "<a href=\""+ varFStr+".html\"><img src=\""+ varFStr+".gif\"></a><br>"
						+ "code "+ varStr+ "</TD>");
				p.println("</TR>");
			}
			p.println("</TABLE>");
		}
		p.println("</div><!-- closing main -->");
		p.println("</BODY>");
		p.println("</HTML>");
	}

	private static String convertFName(String in) {
		StringBuffer sb= new StringBuffer(in);
			// kill spaces
		for (int i = 0; i < sb.length(); i++) 
			if (sb.charAt(i)== ' ')
				sb.deleteCharAt(i--);
		String out= sb.toString();
		out= out.replace('^', 'd');
		out= out.replace('-', 'a');
		out= out.replace(',', 'I');
		if (out.length()> 245)	
			out= out.substring(0, 245);	// max fname len
		
		return out;
	}
	
	static void writeEventsHTML(ASVariation[] vars, PrintStream p) {
//		p.println("<HTML>");
//		p.println("<HEAD>");
//		p.println("<TITLE>Event Details</TITLE>");
//		include(p, STYLE_FILE);
//		p.println("</HEAD>");
//		p.println("<BODY>");
		
		include(p, HEADER_FILE);
		headline(p, "Event Details", null);
		String varStr= vars[0].toString();
		String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
		varStr= varStr.replaceAll("\\s", "");
		p.println("<DIV class=\"userspace\">");
		p.println("<a href=\"landscape.html\"><IMG class=\"pnt\" src=\"/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>");
		p.println("<DIV class=\"userspace\" align=\"center\">");
		p.println("<CENTER><img src=\""+ varFStr+".gif\"><br>" +
				"<FONT size=\"2\" face=\"Arial,Lucida Sans\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ varStr+ "</b></FONT></CENTER><br><br>");
		p.println("</DIV>");
		p.println("<TABLE border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
		// org=D.+melanogaster
		// org=Homo_sapiens&db=hg17
		//&clade=vertebrate&org=Mouse&db=mm8
		int ps= speStr.indexOf(";db=");	
		if (ps< 0)
			ps= speStr.length();
		String spec= speStr.substring(0, ps);
		ps= spec.indexOf("+");	
		if (ps>= 0)
			spec= spec.substring(0, ps)+ " "+ spec.substring(ps+1, spec.length());
		ps= spec.indexOf("_");	
		if (ps>= 0)
			spec= spec.substring(0, ps)+ " "+spec.substring(ps+1, spec.length());
		p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Nr</b></FONT></TD>" +
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Coordinates</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Transcript ID</b> <i>Chromosome</i>: alt. Splice Sites</FONT></TD>" +
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Genome Browser</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">");
		if (speStr!= null)
			p.print(spec+"</FONT></TD></TR>");
		p.println("</FONT></TD></TR>");
		for (int i = 0; i < vars.length; i++) {
			String colStr= "";
			if (i%2 == 0)
				colStr= TABLE_EVEN_COLOR;
			else
				colStr= TABLE_ODD_COLOR;
			
				// event
			
			p.println("<TR bgcolor=\""+colStr+"\">");			
			p.println("\t<TD valign=\"middle\" align=\"center\">"
					+ (i+1)+ "</TD>");
			p.println("\t<TD></TD>"); 
			
				// splice chains
			String pos1Str= "";
			SpliceSite[] su= vars[i].getSpliceUniverse();
			for (int j = 0; j < vars[i].getSpliceChain1().length; j++)
				if (vars[i].getSpliceChain1()[j].getPos()> 0)
					pos1Str+= vars[i].getSpliceChain1()[j].getPos()+ ", ";
				else
					pos1Str+= Math.abs(vars[i].getSpliceChain1()[j].getPos())+ ", ";
			if (pos1Str.length()> 2)
				pos1Str= pos1Str.substring(0, pos1Str.length()- 2);
			String pos2Str= "";
			for (int j = 0; j < vars[i].getSpliceChain2().length; j++) 
				if (vars[i].getSpliceChain2()[j].getPos()> 0)
					pos2Str+= vars[i].getSpliceChain2()[j].getPos()+ ", ";
				else
					pos2Str+= Math.abs(vars[i].getSpliceChain2()[j].getPos())+ ", ";
			if (pos2Str.length()> 2)
				pos2Str= pos2Str.substring(0, pos2Str.length()- 2);
			
				// swap according to rules
			String s= null;
			if (vars[i].getSpliceChain2().length== 0|| 
					(vars[i].getSpliceChain1().length> 0&& vars[i].getSpliceChain1()[0].getPos()< vars[i].getSpliceChain2()[0].getPos()))
				s= vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><b>"
					+ vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</FONT></TD>";
			else
				s= vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"<br><b>"
					+ vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"</FONT></TD>";
			
			p.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+ s);
			p.println("\t<TD></TD>"); 

				// linkout
			SpliceSite[] flanks= vars[i].getFlankingSpliceSites();
			int begin, end;
			if (flanks[0]== null)
				begin= Math.max(vars[i].getTranscript1().get5PrimeEdge(), vars[i].getTranscript2().get5PrimeEdge());
			else
				begin= flanks[0].getPos();
			if (flanks[1]== null)
				end= Math.min(vars[i].getTranscript1().get3PrimeEdge(), vars[i].getTranscript2().get3PrimeEdge());
			else
				end= flanks[1].getPos();
			if (begin< 0) {
				int h= -begin;
				begin= -end;
				end= h;
			}
			String ucscBase= UCSC_URL +"org="+ speStr+";";
			if (annStr!= null)
				ucscBase+= "db="+ annStr+ ";";
			String ucscEventFlanks= ucscBase+ "position="
			+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
			+ UCSC_STANDARD_PAR+ "hgFind.matches="+vars[i].getTranscript1()+","+vars[i].getTranscript2()+",";
			
			begin= su[0].getPos();
			end= su[su.length- 1].getPos();
			if (begin< 0) {
				int h= -begin;
				begin= -end;
				end= h;
			}
			String ucscEventVar= ucscBase+ "position="
			+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
			+ UCSC_STANDARD_PAR+ "hgFind.matches="+vars[i].getTranscript1()+","+vars[i].getTranscript2()+",";
			
			if (annStr!= null)
				p.println("\t<TD valign=\"middle\" align=\"center\" width=\"120\">"
						+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""+ ucscEventVar+"\">Show&nbsp;Alternative&nbsp;Parts>></a><br>" +
							"<a href=\""+ ucscEventFlanks+"\">Show&nbsp;Complete&nbsp;Event>></a><br></FONT></TD>");
			else
				p.println("<TD valign=\"middle\" align=\"center\" width=\"120\"></TD>");	// empty
			p.println("\t<TD></TD>"); 
			
			p.println("</TR>");
		}
		p.println("</TABLE>");
		p.println("</div><!-- closing main -->");
		
		include(p, TRAILER_FILE);
	}

	static void _style_writeHTMLNumbers(ASVariation[] vars, PrintStream p) {
		p.println("<HTML>");
		p.println("<HEAD>");
		p.println("<TITLE>Event Details</TITLE>");
		include(p, STYLE_FILE);
		p.println("</HEAD>");
		p.println("<BODY>");
		
		include(p, HEADER_FILE);
		p.println("<a href=\"statistics.html\">Back to Landscape</a>");
		p.println("<div class=\"title\"><h1>EVENT DETAILS</h1></div><br />");
		String varStr= vars[0].toString();
		String varFStr= convertFName(varStr);	//varStr.replace(' ', '_');
		p.println("<center>");
		p.println("<b>Structure:</b><br><img src=\""+ varFStr+".gif\"><br><b>code "+ varStr+ "</b><br><br>");
		p.println("</center>");
		p.println("<TABLE border=\"0\" cellspacing=\"0\" cellpadding=\"0\" width=\"100%\">");
		// org=D.+melanogaster
		// org=Homo_sapiens&db=hg17
		//&clade=vertebrate&org=Mouse&db=mm8
		int ps= speStr.indexOf("&db=");	
		if (ps< 0)
			ps= speStr.length();
		String spec= speStr.substring(0, ps);
		ps= speStr.indexOf("+");	
		if (ps>= 0)
			spec= speStr.substring(0, ps)+ speStr.substring(ps+1, speStr.length());
		p.println("<TR><TD width=\"120\" bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Event Nr</b></TD>" +
				"<TD width=\"120\" bgcolor=\""+TABLE_HEADER_COLOR+"\"><b>Chrom. Coordinates</b> (<u>Transcript ID</u> <i>chromosome</i>: var. splice sites)</TD>" +
				"<TD width=\"120\" bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><b>Show in UCSC Browser</b><br>("+spec+")</TD></TR>");
		for (int i = 0; i < vars.length; i++) {
			String colStr= "";
			if (i%2 == 0)
				colStr= TABLE_EVEN_COLOR;
			else
				colStr= TABLE_ODD_COLOR;
			
				// event
			
			p.println("<TR>");			
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"middle\" align=\"center\" width=\"120\">Event "
					+ (i+1)+ "</TD>");
			
				// splice chains
			String pos1Str= "";
			SpliceSite[] su= vars[i].getSpliceUniverse();
			for (int j = 0; j < vars[i].getSpliceChain1().length; j++)
				if (vars[i].getSpliceChain1()[j].getPos()> 0)
					pos1Str+= vars[i].getSpliceChain1()[j].getPos()+ ", ";
				else
					pos1Str+= Math.abs(vars[i].getSpliceChain1()[j].getPos())+ ", ";
			if (pos1Str.length()> 2)
				pos1Str= pos1Str.substring(0, pos1Str.length()- 2);
			String pos2Str= "";
			for (int j = 0; j < vars[i].getSpliceChain2().length; j++) 
				if (vars[i].getSpliceChain2()[j].getPos()> 0)
					pos2Str+= vars[i].getSpliceChain2()[j].getPos()+ ", ";
				else
					pos2Str+= Math.abs(vars[i].getSpliceChain2()[j].getPos())+ ", ";
			if (pos2Str.length()> 2)
				pos2Str= pos2Str.substring(0, pos2Str.length()- 2);
			
				// swap according to rules
			String s= null;
			if (vars[i].getSpliceChain2().length== 0|| 
					(vars[i].getSpliceChain1().length> 0&& vars[i].getSpliceChain1()[0].getPos()< vars[i].getSpliceChain2()[0].getPos()))
				s= vars[i].getTranscript1().getTranscriptID()+ "</u> <i>chr"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><u>"
					+ vars[i].getTranscript2().getTranscriptID()+ "</u> <i>chr"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</TD>";
			else
				s= vars[i].getTranscript2().getTranscriptID()+ "</u> <i>chr"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"<br><u>"
					+ vars[i].getTranscript1().getTranscriptID()+ "</u> <i>chr"
					+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"</TD>";
			
			p.print("\t<TD bgcolor=\""+colStr+"\" valign=\"top\" width=\"120\"><u>"+ s);
	
				// linkout
			SpliceSite[] flanks= vars[i].getFlankingSpliceSites();
			int begin, end;
			if (flanks[0]== null)
				begin= Math.max(vars[i].getTranscript1().get5PrimeEdge(), vars[i].getTranscript2().get5PrimeEdge());
			else
				begin= flanks[0].getPos();
			if (flanks[1]== null)
				end= Math.min(vars[i].getTranscript1().get3PrimeEdge(), vars[i].getTranscript2().get3PrimeEdge());
			else
				end= flanks[1].getPos();
			if (begin< 0) {
				int h= -begin;
				begin= -end;
				end= h;
			}
			String ucscEventFlanks= "http://genome.ucsc.edu/cgi-bin/hgTracks?" +
			"org="+ speStr+ "&position=chr"
			+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK);
			begin= su[0].getPos();
			end= su[su.length- 1].getPos();
			if (begin< 0) {
				int h= -begin;
				begin= -end;
				end= h;
			}
			String ucscEventVar= "http://genome.ucsc.edu/cgi-bin/hgTracks?" +
			"org="+ speStr+ "&position=chr"
			+ vars[i].getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK);
			
			p.println("\t<TD bgcolor=\""+colStr+"\" valign=\"middle\" align=\"center\" width=\"120\">"
					+ "<a href=\""+ ucscEventVar+"\">Variation>></a><br>" +
						"<a href=\""+ ucscEventFlanks+"\">Flanks>></a><br></TD>");
			
			p.println("</TR>");
		}
		p.println("</TABLE>");
		p.println("</div><!-- closing main -->");
		p.println("</BODY>");
		p.println("</HTML>");
	}
}
