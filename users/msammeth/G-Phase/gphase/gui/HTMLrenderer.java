/*
 * Created on Jun 17, 2007
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.gui;

import gphase.model.ASVariation;
import gphase.tools.File;

public class HTMLrenderer {
	ASVariation[][][] vars= null;
	File outputDir= null;
	String genome= null;
	                
	HTMLrenderer(ASVariation[][][] newVars, File newOutDir, String newGenome) {
		this.vars= newVars;
		this.outputDir= newOutDir;
		this.genome= newGenome;
	}

	protected static void writeEventsHTML(ASVariation[] vars, PrintStream p) {
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

	public static void writeHTML() {
		
		if (writePie)
			writePiePicture();			
		writeStructurePics();
		
		writeHTMLlandscape(vars, p);
		
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			String varStr= vars[i][0].toString();
			String varFStr= convertFName(varStr)+ ".html";	//varStr.replace(' ', '_')+ ;
			
			p= new PrintStream(fName+ varFStr);
			writeEventsHTML(vars[i], p);
			p.flush(); p.close();
		}
		
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
	
	
	
}
