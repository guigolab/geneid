package gphase;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import gphase.gui.SpliceOSigner;
import gphase.io.DomainToGenomeMapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.ASVariationWithRegions;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;
import gphase.tools.File;
import gphase.tools.Formatter;

// human_hg17_RefSeqGenes_fromUCSC_addDomain_addmRNA_addSplicedESTs_sorted.gtf
// human_hg17_RefSeqGenes_fromUCSC_addDomains.gtf
//  H:\annotations\human_hg17_RefSeqGenes_mRNA_CDS_fromUCSC_add_Domains.gtf
public class StructuratorClaudia extends AStaLaVista {
	
		final static boolean LINKOUT_PROBES= true;
	final static String CLAUDIA_HTML_SUM_FNAME= "claudia.html";
	final static String[] UCSC_CT_SJ= new String[] {
		"http://www.exonhit.com/bed/ASD_Cancersv1_0",
		"http://www.exonhit.com/bed/ASD_Cancersv1_0/Claudia",	// cancer chip> 300
		"http://www.exonhit.com/bed/SR",
		"http://www.exonhit.com/bed/ASD_Add1"
	};
	
	
	File clauFile= null;
	HashMap mapClauHTML= new HashMap();
	
	protected void writeHTMLeventColSchema(ASVariation var, PrintStream p) {
		super.writeHTMLeventColSchema(var, p);

		p.print("\t<TD>&nbsp;</TD>");
		Vector v= (Vector) var.getAttribute(Claudia.ID_EVENT_TAG);
		p.print("\t<TD>");
		for (int i = 0; v!= null&& i < v.size(); i++) 
			p.print("<a href=\""+CLAUDIA_HTML_SUM_FNAME+"#"+v.elementAt(i)+"\">"+v.elementAt(i)+"</a><br>");
		p.print("\t</TD>");

	}
	public void writeHTML(ASVariation[][][] vars) {
		super.writeHTML(vars);
		
			// write summary for claudia
			// hash in writEvent
//		for (int i = 0; i < vars.length; i++) {
//			for (int j = 0; j < vars[i].length; j++) {
//				Vector v= (Vector) vars[i][j][0].getAttribute(Claudia.ID_EVENT_TAG);
//				for (int x = 0; v!= null&& x < v.size(); ++x) {
//					Vector vv= (Vector) mapClauHTML.get(v.elementAt(x));
//					if (vv== null)
//						vv= new Vector();
//					vv.add(getFileName(vars[i][j][0].toStringStructureCode())+".html#"+j);
//					mapClauHTML.put(v.elementAt(x), vv);
//				}
//
//			}
//		}
		
		PrintStream p= null;
		try {
			p= new PrintStream(outDir.getAbsolutePath()+File.separator+CLAUDIA_HTML_SUM_FNAME);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		include(p, HEADER_FILE);
		headline(p, "SJ event ID summary", null);
		p.println("<DIV>");
		p.println("<a href=\""+HTML_NAME_LANDSCAPE+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
		"<FONT face=\"Arial,Lucida Sans\" size=\"1\">landscape</FONT></a>");
		p.println("<a href=\""+FNAME_REGIONS+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">"+ovlFeature+"</FONT></a><br>");
		p.print("</DIV>");

		Object[] keys= mapClauHTML.keySet().toArray();
		java.util.Arrays.sort(keys);
		for (int i = 0; i < keys.length; i++) {
			Vector v= (Vector) mapClauHTML.get(keys[i]);
			for (int j = 0; j < v.size(); j++) 
				p.print("<a name=\""+keys[i]+"\"></a><a href=\""+v.elementAt(j)+"\">"+keys[i]+"</a><br>");
		}

		p.print("<br><br>");
		p.print("</div><!-- closing main -->\n");
		
		include(p, TRAILER_FILE);
		p.flush(); p.close();

	}
	
	
	public void overlap(ASVariation[][][] vars, String featureID) {
		if (vars== null)
			return;
		
		super.overlap(vars, featureID);
		
		GTFChrReader reader= new GTFChrReader(clauFile.getAbsolutePath());
		reader.setChromosomeWise(true);
		reader.setReadGene(false);
		reader.setReadGTF(true);
		reader.setSilent(true);
		reader.setReadFeatures(new String[]{Claudia.ID_EXON_ROLE_B});
		reader.sweepToChromosome(vars[0][0][0].getTranscript1().getChromosome());
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		GTFObject[] obs= reader.getGtfObj();
		if (obs== null)
			return;
		
		HashMap mapBprobes= new HashMap();	// claudia ev IDs x V(DirectedRegion)
		for (int i = 0; i < obs.length; i++) {
			String id= null;
			id= obs[i].getSource();
			int p1= id.indexOf('.');
			int p2= id.lastIndexOf('.');
			if (p1== p2) {
				System.out.println("WARNING: strange event ID "+id);
			}
			DirectedRegion newReg= new DirectedRegion(obs[i].getStart(), obs[i].getEnd(), obs[i].getStrand());
			newReg.setChromosome(obs[i].getChromosome());
			newReg.addAttribute(Claudia.ID_EVENT_TAG, id);

			if (obs[i].getFeature().equals(Claudia.ID_EXON_ROLE_B)) {
				Vector v= (Vector) mapBprobes.remove(id);
				if (v== null) {
					v= new Vector();
				} 
				v.add(newReg);
				mapBprobes.put(id, v);	// evID
			}
		}

			// overlap
		for (int i = 0; i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				DirectedRegion[] exRegs1= vars[i][j][0].getExonicRegions1();
				DirectedRegion[] exRegs2= vars[i][j][0].getExonicRegions2();
				Object[] vals= mapBprobes.values().toArray();
					// TODO this is very inefficient !!!
				for (int k = 0; k < vals.length; k++) {	// clau events
					Vector v= (Vector) vals[k];
					int m;
					for (m = 0; exRegs1!= null&& m < v.size(); m++) {	// all b-probes of an event
						int n;
						for (n = 0; n < exRegs1.length; n++) 
							if (exRegs1[n].contains((DirectedRegion) v.elementAt(m)))
								break;
						if (n== exRegs1.length)	// at least one probe not found
							break;
					}
					if (m< v.size()&& exRegs2!= null) {	// try in second transcript..
						for (m = 0; m < v.size(); m++) {	// all b-probes of an event
							int n;
							for (n = 0; n < exRegs2.length; n++) 
								if (exRegs2[n].contains((DirectedRegion) v.elementAt(m)))
									break;
							if (n== exRegs2.length)	// at least one probe not found
								break;
						}
					}
					if (m== v.size()) {	// all found (in trpt1 or trpt2)
						for (int x = 0; x < vars[i][j].length; x++) {
							Vector vv= (Vector) vars[i][j][x].getAttribute(Claudia.ID_EVENT_TAG);
							if (vv== null)
								vv= new Vector();
							vv.add(((DirectedRegion) ((Vector) vals[k]).elementAt(0)).getAttribute(Claudia.ID_EVENT_TAG));
							vars[i][j][x].addAttribute(Claudia.ID_EVENT_TAG, vv);
						}
					}
						
				}
			}
		}
	} 
	
	public static void main(String[] args) {
		long t0= System.currentTimeMillis();
		StructuratorClaudia astalavista= new StructuratorClaudia();
		parseArguments(args, astalavista);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-sj")) {
				astalavista.setClauFile(new File(args[i+1]));
				++i;
				continue;
			}
		}
		astalavista.run();
		
		
		long t1= System.currentTimeMillis();
		System.out.println("time: "+ (t1-t0)+"[msec]"); 
		
	}
	
	public static ASVariationWithRegions[][] getEvents(File inFile) {
		return getEvents(inFile, new File("claudia"+ File.separator+ "VARIANTS_EVENTS_EXON_subst.txt.gtf"));
	}
	public static ASVariationWithRegions[][] getEvents(File inFile, File clauFile) {
		long t0= System.currentTimeMillis();
		
		System.out.println("Reading events "+clauFile.getFileNameOnly()+".");
		GTFChrReader reader= new GTFChrReader(clauFile.getAbsolutePath());
		reader.setChromosomeWise(true);
		reader.setReadGene(false);
		reader.setReadGTF(true);
		reader.setReadFeatures(new String[]{Claudia.ID_EXON_ROLE_A3,
				Claudia.ID_EXON_ROLE_A5, Claudia.ID_EXON_ROLE_B, Claudia.ID_EVENT_FEATURE});
		
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		GTFObject[] obs= reader.getGtfObj();
		
		GTFChrReader reader2= new GTFChrReader(inFile.getAbsolutePath());
		reader2.setChromosomeWise(true);
		reader2.sweepToChromosome(obs[0].getChromosome());
		reader2.setReadFeatures(new String[] {GTFObject.EXON_FEATURE_TAG, GTFObject.CDS_FEATURE_TAG, DomainToGenomeMapper.GTF_DOMAIN_FEATURE_TAG});
		//reader2.setSilent(true);
		Gene[] genes= null;
		if (obs!= null) {
			try {
				reader2.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader2.getGenes();
		}
				
		HashMap eventMap= new HashMap();
		while (obs!= null&& genes!= null) {
			
			Claudia clau= new Claudia(genes, obs);
			ASVariationWithRegions[][][] ev= clau.getClauEvents2();
			
			for (int i = 0; ev!= null&& i < ev.length; i++) {
				Vector v= (Vector) eventMap.remove(ev[i][0][0].toString());
				if (v== null)
					v= new Vector();
				for (int j = 0; j < ev[i].length; j++) {
					v.add(ev[i][j][0]);	// only the first one, non-redundant
				}
				eventMap.put(ev[i][0][0].toString(), v);
			}
				// next
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			obs= reader.getGtfObj();
			if (obs!= null) {
//				System.out.println("Reading annotation "+inFile.getFileNameOnly()+".");
//				GTFChrReader reader2= new GTFChrReader(inFile.getAbsolutePath());
//				reader2.setChromosomeWise(true);
				reader2.sweepToChromosome(obs[0].getChromosome());
//				reader2.setReadFeatures(new String[] {GTFObject.EXON_ID_TAG, DomainToGenomeMapper.GTF_DOMAIN_FEATURE_TAG});
//				reader2.setSilent(true);
				try {
					reader2.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				genes= reader2.getGenes();
			}
			System.gc();
		}
		
			// output
		ASVariationWithRegions[][] allVars= new ASVariationWithRegions[eventMap.size()][];
		Object[] keys= eventMap.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			allVars[i]= (ASVariationWithRegions[]) Arrays.toField(eventMap.get(keys[i]));
			//java.util.Arrays.sort(allVars[i]);
		}
		allVars= (ASVariationWithRegions[][]) Arrays.sortNDField((Object) allVars);
		
		System.out.println("Read in GTF events: "+Claudia.cntGTFevents);
		System.out.println("Mapped events: "+Claudia.cntMappedEvents);
		System.out.println("Extremities: "+Claudia.cntExtrem);
		System.out.println("Invalid Flanks: "+Claudia.cntInvFlanks);
		System.out.println("Failed events: "+Claudia.cntFailedEvents);
		System.out.println("Non-Domain events: "+Claudia.cntNoDomainEvents);
		System.out.println("Domain events: "+Claudia.cntDomainEvents);
		
		for (int i = 0; i < allVars.length; i++) 
			System.out.println(allVars[i][0]+"\t"+allVars[i].length);
		return allVars;
	}

	protected static void writeEventsHTML(ASVariationWithRegions[] vars, PrintStream p) {
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
							"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Schema</b></FONT>" +
							"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
							"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Domains</b></FONT>" +
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
						
							// event nr
						p.println("<TR bgcolor=\""+colStr+"\">");			
						p.println("\t<TD valign=\"middle\" align=\"center\">"
								+ (i+1)+ "</TD>");
						p.println("\t<TD></TD>"); 
						
							// picture
						String domPicName= writeDomainStructPic(vars[i]);
						p.println("\t<TD valign=\"middle\" align=\"center\"><img src=\""+ domPicName+"\">"
								+ "</TD>");
						p.println("\t<TD></TD>"); 
						
							// domains
						DirectedRegion[] dom1= vars[i].getReg1();
						DirectedRegion[] dom2= vars[i].getReg2();
						HashMap<String,Vector<DirectedRegion>> domGrpHash= new HashMap<String,Vector<DirectedRegion>>();
						p.print("\t<TD valign=\"middle\" align=\"left\">");
						for (int j = 0; j < dom1.length; j++) {
							String id= ASVariationWithRegions.stripOrderSuffix(dom1[j].getID());
							Vector<DirectedRegion> v= domGrpHash.get(id);
							if (v== null)
								v= new Vector<DirectedRegion>();
							v.add(dom1[j]);
							domGrpHash.put(id, v);
						}
						for (int j = 0; j < dom2.length; j++) {	// spacer between transcripts
							String id= ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
							Vector<DirectedRegion> v= domGrpHash.get(id);
							if (v== null)
								continue;
							v.add(null);
							domGrpHash.put(id, v);
						}
						for (int j = 0; j < dom2.length; j++) {
							String id= ASVariationWithRegions.stripOrderSuffix(dom2[j].getID());
							Vector<DirectedRegion> v= domGrpHash.get(id);
							if (v== null) {
								v= new Vector<DirectedRegion>();
								v.add(null);	// spacer
							}
							v.add(dom2[j]);
							domGrpHash.put(id, v);
						}
						
						Object[] domVal= domGrpHash.keySet().toArray();
						for (int j = 0; j < domVal.length; j++) {
							String id= (String) domVal[j];
							Vector<DirectedRegion> v= (Vector<DirectedRegion>) domGrpHash.get(id);
							Color c= SpliceOSigner.getDomainColor(v.elementAt(v.size()- 1).getID());	// to be stripped, take original
							String dotFName= writeDotPic(c);
							p.print("\t<img src=\""+ dotFName+"\">&nbsp;<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+id+"</b> ");
							if (v.elementAt(0)!= null)
								p.print(vars[i].getTranscript1().getTranscriptID()+" ");
							for (int k = 0; k < v.size(); k++) {
								if (v.elementAt(k)== null) {
									p.print("<br>"+vars[i].getTranscript2().getTranscriptID()+" ");
									continue;
								}
								p.print(Formatter.fprint(v.elementAt(k).getScore(),2)+" ");
							}
							p.print("<br>");
						}
						p.println("\t</FONT></TD><TD></TD>"); 
						
						
						
							// coordinates
						String pos1Str= "";
						SpliceSite[] su= vars[i].getSpliceUniverse();
						for (int j = 0; j < vars[i].getSpliceChain1().length; j++)
							if (vars[i].getSpliceChain1()[j].getPos()> 0)
								pos1Str+= vars[i].getSpliceChain1()[j].getPos()+ ",";
							else
								pos1Str+= Math.abs(vars[i].getSpliceChain1()[j].getPos())+ ",";
						if (pos1Str.length()> 2)
							pos1Str= pos1Str.substring(0, pos1Str.length()- 1);
						String pos2Str= "";
						for (int j = 0; j < vars[i].getSpliceChain2().length; j++) 
							if (vars[i].getSpliceChain2()[j].getPos()> 0)
								pos2Str+= vars[i].getSpliceChain2()[j].getPos()+ ",";
							else
								pos2Str+= Math.abs(vars[i].getSpliceChain2()[j].getPos())+ ",";
						if (pos2Str.length()> 2)
							pos2Str= pos2Str.substring(0, pos2Str.length()- 1);
						
							// swap according to rules
						String s= vars[i].getAttribute(Claudia.ID_EVENT_TAG)+"</b><br><i>"
								+ vars[i].getTranscript1().getTranscriptID()+ "</i> "
								+ vars[i].getGene().getChromosome()+ ":"+ pos1Str+"<br><i>"
								+ vars[i].getTranscript2().getTranscriptID()+ "</i> "
								+ vars[i].getGene().getChromosome()+ ":"+ pos2Str+"</FONT></TD>";
	
	//					if (vars[i].getSpliceChain2().length== 0|| 
	//							(vars[i].getSpliceChain1().length> 0&& vars[i].getSpliceChain1()[0].getPos()< vars[i].getSpliceChain2()[0].getPos()))
	//						s= vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
	//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"<br><b>"
	//							+ vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
	//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"</FONT></TD>";
	//					else
	//						s= vars[i].getTranscript2().getTranscriptID()+ "</b> <i>"
	//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos2Str+"<br><b>"
	//							+ vars[i].getTranscript1().getTranscriptID()+ "</b> <i>"
	//							+ vars[i].getGene().getChromosome()+ "</i>: "+ pos1Str+"</FONT></TD>";
						
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
						String ucscBase= UCSC_GB_CGI +"org="+ speStr+";";
						String refName= vars[i].getGene().getChromosome()+"_"+vars[i].getGene().getNameTranscript();
						if (annStr!= null)
							ucscBase+= "db="+ annStr+ ";";
						//ucscBase+= "&hgt.customText=http://genome.imim.es/~msammeth/tstBED/"+refName+".bed;";	// must be hg17
						// http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg17&hgt.customText=http://genome.imim.es/g
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

	public static void writeHTML(ASVariationWithRegions[][] vars, String fName) {
		
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
				StructuratorClaudia.writeEventsHTML(vars[i], p);
				p.flush(); p.close();
			}
			
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}

	public static ASVariationWithRegions[][] getEvents2(File inFile) {
				long t0= System.currentTimeMillis();
				
				System.out.println("Retrieving AS events from "+inFile.getFileNameOnly()+".");
				GTFChrReader reader= new GTFChrReader(inFile.getAbsolutePath());
				reader.setChromosomeWise(true);
				reader.setReadFeatures(new String[]{"exon", "CDS", "domain"});
				DirectedRegion regX= new DirectedRegion(81009006, 81113783, -1);
				regX.setChromosome("chr15");
				reader.setFiltRegs(new DirectedRegion[] {regX});
				try {
					reader.read();
				} catch (Exception e) {
					e.printStackTrace();
				}
				Gene[] genes= reader.getGenes();
				
				HashMap<String,Vector<ASVariation>> eventMap= new HashMap<String,Vector<ASVariation>>();
				while (genes!= null) {
					for (int i = 0; i < genes.length; i++) {
						if (!genes[i].isProteinCoding())
							continue;
	//					if (!checkGene(genes[i]))
	//						continue;
	
						ASVariation[] varsSimple= genes[i].getASVariations(ASMultiVariation.FILTER_NONE);
						if (varsSimple== null|| varsSimple.length== 0)
							continue;
						
	
							// get genomic domain regions
						HashMap<String,Vector<DirectedRegion>> domMap= new HashMap<String,Vector<DirectedRegion>>();
						Transcript[] trpts= genes[i].getTranscripts();
						for (int j = 0; j < trpts.length; j++) {
							Vector<DirectedRegion> trptDomV= new Vector<DirectedRegion>();
							if (trpts[j].getAttributes()!= null) {
								Object[] keys= trpts[j].getAttributes().keySet().toArray();
								for (int m = 0; m < keys.length; m++) {
									DirectedRegion[] diregV= (DirectedRegion[]) trpts[j].getAttribute(keys[m]);							
									int minStart= Integer.MAX_VALUE, maxEnd= Integer.MIN_VALUE;
									String lastID= null;
									for (int k = 0; k < diregV.length; k++) {
										if (diregV[k].get5PrimeEdge()< minStart)
											minStart= diregV[k].get5PrimeEdge();
										if (diregV[k].get3PrimeEdge()> maxEnd)
											maxEnd= diregV[k].get3PrimeEdge();
									}
									DirectedRegion reg= new DirectedRegion();
									reg.setStrand(trpts[j].getStrand());
									reg.setChromosome(trpts[j].getChromosome());
									reg.setStart(minStart);		// TODO check for 0-based
									reg.setEnd(maxEnd);		// -1 +1
									reg.setScore(diregV[0].getScore());
									reg.setID(diregV[0].getID());
									trptDomV.add(reg);
								}
							}
							for (int k = 0; k < trptDomV.size(); k++) {
								for (int m = k+1; m < trptDomV.size(); m++) {
									if (((DirectedRegion) trptDomV.elementAt(m)).overlaps((DirectedRegion) trptDomV.elementAt(k)))
										System.out.println("WARNING: transcript "+trpts[j].getTranscriptID()+" has overlapping domains.");
								}
							}
							
							if (trptDomV.size()> 0) {
								domMap.put(trpts[j].getTranscriptID(),trptDomV);
							}
						}
		
						// cnt events
						Vector<ASVariationWithRegions> varsDomV= new Vector<ASVariationWithRegions>(); 
						ASVariationWithRegions[] vars= new ASVariationWithRegions[varsSimple.length];
						for (int j = 0; j < vars.length; j++) 
							vars[j]= new ASVariationWithRegions(varsSimple[j], 
									domMap.get(varsSimple[j].getTranscript1().getTranscriptID()),
									domMap.get(varsSimple[j].getTranscript2().getTranscriptID()));
						Comparator domOvlCompi= new ASVariationDomain.IdentityComparator();
						for (int j = 0; vars!= null&& j < vars.length; j++) {
							DirectedRegion reg= vars[j].getRegion();
							Vector<DirectedRegion> vRegs= (Vector<DirectedRegion>) domMap.get(vars[j].getTranscript1().getTranscriptID());
							int k;
							for (k = 0; vRegs!= null&& k < vRegs.size(); k++) {
								if (vRegs.elementAt(k).overlaps(reg))
									break;
							}
							if (vRegs!= null&& k< vRegs.size()) {
								varsDomV.add(vars[j]);
								continue;
							}
							vRegs= (Vector<DirectedRegion>) domMap.get(vars[j].getTranscript2().getTranscriptID());
							for (k = 0; vRegs!= null&& k < vRegs.size(); k++) {
								if (vRegs.elementAt(k).overlaps(reg))
									break;
							}
							//if (vRegs!= null&& k< vRegs.size())
							if (vars[j].isAffectingCDS())		// ***
								varsDomV.add(vars[j]);
						}
						if (varsDomV.size()== 0)
							continue;
						vars= (ASVariationWithRegions[]) ASVariation.removeRedundancy((ASVariation[]) Arrays.toField(varsDomV), domOvlCompi);
						for (int j = 0; vars!= null&& j < vars.length; j++) {
							Vector<ASVariation> v= eventMap.remove(vars[j].toString());
							if (v== null)
								v= new Vector<ASVariation>();
							v.add(vars[j]);
							eventMap.put(vars[j].toString(), v);
						}
					}
					
					try {
						reader.read();
					} catch (Exception e) {
						e.printStackTrace();
					}
					genes= reader.getGenes();
					System.gc();
				}
				
					// output
				ASVariationWithRegions[][] allVars= new ASVariationWithRegions[eventMap.size()][];
				Object[] keys= eventMap.keySet().toArray();
				for (int i = 0; i < keys.length; i++) {
					allVars[i]= (ASVariationWithRegions[]) Arrays.toField(eventMap.get(keys[i]));
					java.util.Arrays.sort(allVars[i]);
				}
				allVars= (ASVariationWithRegions[][]) Arrays.sort2DFieldRev(allVars);
				
				for (int i = 0; i < allVars.length; i++) 
					System.out.println(allVars[i][0]+"\t"+allVars[i].length);
				return allVars;
			}

	public File getClauFile() {
		return clauFile;
	}

	public void setClauFile(File clauFile) {
		this.clauFile = clauFile;
	}
	protected void writeHTMLeventGroupColRefPair(ASVariation[] vars, PrintStream p) {
		super.writeHTMLeventGroupColRefPair(vars, p);
		p.print("\t<TD>&nbsp&nbsp</TD>\n"); 

		Vector v= (Vector) vars[0].getAttribute(Claudia.ID_EVENT_TAG);
		p.print("\t<TD>");
		for (int i = 0; v!= null&& i < v.size(); i++) {
			p.print("<a href=\""+CLAUDIA_HTML_SUM_FNAME+"#"+v.elementAt(i)+"\">"+v.elementAt(i)+"</a><br>");
			Vector vv= (Vector) mapClauHTML.get(v.elementAt(i));
			if (vv== null)
				vv= new Vector();
			vv.add(getFileName(vars[0].toCoordinates())+".html");
			mapClauHTML.put(v.elementAt(i), vv);
		}
		p.print("</TD>\n");
	}

	protected void writeHTMLeventGroupTableHeader(PrintStream p) {
		p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Nr</b></FONT></TD>" +
						"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Transcript Pairs</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Reference Pair</b> <u>all Pairs</u></FONT></TD>"+
						"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
		"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>SJ ID</b></FONT>" +
			"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
		if (ovlFeature!= null)
			p.print("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Overlap with</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">"+ovlFeature+"s</FONT></TD>"+
						"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
		
		p.print("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Genome Browser</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">");
	
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
		if (speStr!= null)
			p.print(spec+"</FONT></TD></TR>");
		p.println("</FONT></TD></TR>");
	}

	protected void writeHTMLeventGroupColUCSC(ASVariation var, PrintStream p) {
	
				// linkout
			SpliceSite[] flanks= var.getFlankingSpliceSites();
			int begin, end;
			if (flanks[0]== null)
				begin= Math.max(var.getTranscript1().get5PrimeEdge(), var.getTranscript2().get5PrimeEdge());
			else
				begin= flanks[0].getPos();
			if (flanks[1]== null)
				end= Math.min(var.getTranscript1().get3PrimeEdge(), var.getTranscript2().get3PrimeEdge());
			else
				end= flanks[1].getPos();
			if (begin< 0) {
				int h= -begin;
				begin= -end;
				end= h;
			}
			String ucscBase= UCSC_GB_CGI + UCSC_LOAD_DEFAULTS;
			ucscBase+= ";"+"org="+ speStr+";";
			if (genomeVer!= null)
				ucscBase+= "db="+ genomeVer+ ";";
			String ucscEventFlanks= ucscBase+ "position="
			+ var.getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
			+ UCSC_STANDARD_PAR+ "hgFind.matches="+var.getTranscript1()+","+var.getTranscript2()+",";
			
	//		SpliceSite[] su= var.getSpliceUniverse();
	//		begin= su[0].getPos();
	//		end= su[su.length- 1].getPos();
	//		if (begin< 0) {
	//			int h= -begin;
	//			begin= -end;
	//			end= h;
	//		}
	//		String ucscEventVar= ucscBase+ "position="
	//		+ var.getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
	//		+ UCSC_STANDARD_PAR+ "hgFind.matches="+var.getTranscript1()+","+var.getTranscript2()+",";
			
			if (genomeVer!= null) {
				p.print("\t<TD valign=\"middle\" align=\"center\" width=\"120\">"
						+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\">" +
								"<a href=\""+ ucscEventFlanks+"\" target=\"FCwin\">Show&nbsp;Event>></a>\n");	// complete event
				if (LINKOUT_DOMAINS) {
					String ucscDomains= ucscBase
						+ UCSC_LOAD_CT+ UCSC_CT_DOMAIN_URL+ "/" 
						+ var.getGene().getChromosome()+"_"+var.getGene().getNameTranscript().getTranscriptID()+ ".bed;"
						+ "position="				
						+ var.getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
						+ "hgFind.matches="+var.getTranscript1()+","+var.getTranscript2()+",";
					p.print("<a href=\""+ ucscDomains+"\" target=\"FCwin\">Show&nbsp;Domains>></a>\n");
				}

				if (LINKOUT_PROBES) {
					Vector v= (Vector) var.getAttribute(Claudia.ID_EVENT_TAG);
					if (v!= null) {
						String evFile= (String) v.elementAt(0); //123.1.1.txt for sj
//							// cancer chip
//						String bed= UCSC_CT_SJ[0];
//						evFile= evFile.substring(0, evFile.indexOf('.'));
//						int x= Integer.parseInt(evFile);
//						if (x> 300)
//							bed= UCSC_CT_SJ[1];	
						
							// sf16
						//String bed= UCSC_CT_SJ[2];
						
						// sf_hg17
						evFile= evFile.substring(0, evFile.indexOf('.'));
						String bed= UCSC_CT_SJ[3];
												
						String ucscDomains= ucscBase
							+ UCSC_LOAD_CT+  bed+ "/"
							+ evFile+ ".txt;"
							+ "mrna=pack;"
							+ "position="				
							+ var.getGene().getChromosome()+ ":"+(begin- UCSC_FLANK)+"-"+ (end+ UCSC_FLANK)+ ";"
							+ "hgFind.matches="+var.getTranscript1()+","+var.getTranscript2()+",";
						p.print("<a href=\""+ ucscDomains+"\" target=\"FCwin\">Show&nbsp;Probes>></a>\n");
					}
				}
			} else
				p.println("<TD valign=\"middle\" align=\"center\" width=\"120\"></TD>");	// empty
			p.println("\t<TD></TD>"); 
	
		}
	protected void writeHTMLeventTableHeader(PrintStream p) {
		p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Pair<br>Nr</b></FONT></TD>" +
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Transcripts</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Reference Pair</b> <u>all Pairs</u></FONT></TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
		if (ovlFeature!= null)
			p.print("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Overlap with</b></FONT>" +
					"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">"+ovlFeature+"s</FONT></TD>"+
					"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
		
		p.print("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>SJ event ID</b></FONT>" +
		"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">"+
		"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>");
		
		p.print("<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Genome Browser</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\">");
	
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
		if (speStr!= null)
			p.print(spec+"</FONT></TD></TR>");
		p.println("</FONT></TD></TR>");
	}
	
	@Override
	void writeHTMLlandscapeLinkouts(PrintStream p) {
		super.writeHTMLlandscapeLinkouts(p);
		p.print("<a href=\""+outDir.getAbsolutePath()+File.separator+CLAUDIA_HTML_SUM_FNAME+"\">SJ event summary</a><br>\n");
	}
	void writeHTMLlandscapeLinkbacks(PrintStream p) {
		super.writeHTMLlandscapeLinkbacks(p);
		p.print("<a href=\""+CLAUDIA_HTML_SUM_FNAME+"\">" +
				"<IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">SJ event summary</FONT></a><br>\n");
	}
	
	protected void writeHTMLovlLinkbacks(PrintStream p) {
		super.writeHTMLovlLinkbacks(p);
		p.print("<a href=\""+CLAUDIA_HTML_SUM_FNAME+"\">" +
				"<IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">SJ event summary</FONT></a><br>\n");
	}


}
