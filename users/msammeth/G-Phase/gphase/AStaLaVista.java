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
import java.util.HashMap;
import java.util.Properties;
import java.util.Vector;

import javax.swing.JFrame;

import org.freehep.util.export.ExportFileType;

import com.sun.corba.se.spi.ior.WriteContents;

import gphase.algo.ASAnalyzer;
import gphase.algo.ASManager;
import gphase.algo.FilterFactory;
import gphase.gui.Circle;
import gphase.gui.CopyOfSpliceOSigner;
import gphase.gui.SpliceOSigner;
import gphase.gui.pie.Pie;
import gphase.io.TabDelimitedFormatWrapper;
import gphase.io.gtf.EncodeWrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.ASVariationWithRegions;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;


// /home/msammeth/annotations/hg17_RefSeqGenes_fromUCSC_070611_CDS_mRNAs_fromUCSC_070611_CDS_addDomains.gtf
public class AStaLaVista {
	final static boolean LINKOUT_DOMAINS= true;
	
	final static String HTML_NAME_LANDSCAPE= "landscape.html";
	
	final static String UCSC_RESET_CART= "http://genome.ucsc.edu/cgi-bin/cartReset?destination=/cgi-bin/hgTables";
	final static String UCSC_LOAD_DEFAULTS= "hgS_doLoadUrl=submit;hgS_loadUrlName=http://genome.imim.es/~msammeth/customTracks/ucsc_def.txt";
	final static String UCSC_LOAD_CT= "hgt.customText=";
	final static String UCSC_CT_DOMAIN_URL= "http://genome.imim.es/~msammeth/customTracks/hg17/domains_Pfam"; // + / + RefTrpt.bed
	
	final static String FNAME_PIE= "distribution.png";
	final static String FNAME_LANDSCAPE= "landscape.html"; 
	final static String LOC_FLORET= "../pics/floret_cyan.jpg";
	final static String FNAME_REGIONS= "region.html";
	
	final static String ANNOTATION_DIR= "annotation";
	final static String UCSC_GB_CGI= "http://genome.ucsc.edu/cgi-bin/hgTracks?";
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
	final static Color ME_COL= new Color(246, 35, 122);
	final static Color OTHERS_COL= new Color(192, 192, 192);
	final static int UCSC_FLANK= 50;
	String speStr= null, genomeVer= null, ovlFeature= null;
	boolean codingTranscripts= false, noncodTranscripts= false, nmd= false, writeASTA= false, writeGTF= false, writeLandscape= true, writeHTML= false,
		ovlOnly= false, ovlExcluded= false;
	gphase.tools.File inFile= null, outDir= null;
	String[] refTrptIDs= null;
	HashMap filterMap= new HashMap();
	
		// rubbish?
	static int filterCode= ASMultiVariation.FILTER_HIERARCHICALLY;
	static int codingCode= ASVariation.TYPE_ALL;
	static String evFilterStr= "";

	boolean filtGTAG= false;	// filter GTAG introns (on event level)
	static boolean killError= false;
	HashMap mapRegions= null;
	
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

	static String usage= "AStaLaVista [flags] <input GTF file>\n\n"+
		"[options]\n"+
		"-f , --feature, -nf, --notfeature, -of, --onlyfeature <gtf feature>\noverlap with regions annotated in the gtf with <gtf feature>, e.g. \'domain\'"+
		"-g , --genome <genome version>\ngenome version in the format \'species_build\', e.g. \'human_hg17\'"+
		"--gtag\nfilter EVENTS with non-GTAG introns"+
		"-h , --html\nwrite html files"+
		"-o , --output <output directory>\nwrite output to <output directory> (default: directory of the input file)\n\n"+
		"-r , --reference <file>\nread reference transcript IDs from <file>"+
		"-w , --write <format>\nwrite output format, allowed <format> \'asta\', \'gtf\' or \'none\' (default)\n"+
		"--filter <className> <methodName>[_<method2Name>_..]\nfilter class with given method, e.g. gphase.model.Transcript isCoding or gphase.model.ASVariation hasOnlyGTAGintrons_isASevent_isContained5UTR or isAffectingCDS...\n\n";
		
	static protected void parseArguments(String[] args, AStaLaVista asta) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-i")|| args[i].equals("--input")) {
				gphase.tools.File inFile= new gphase.tools.File(args[i+1]);
				if (!inFile.exists()) {
					System.err.println("<input file> not found "+args[i+1]);
					System.exit(-1);
				}
				asta.setInFile(inFile);
				
//				if (asta.getOutDir()== null)
//					asta.setOutDir(new gphase.tools.File(inFile.getPathOnly()));
				++i;
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-f")|| args[i].equalsIgnoreCase("--feature")||
					args[i].equalsIgnoreCase("-nf")|| args[i].equalsIgnoreCase("--notfeature")||
					args[i].equalsIgnoreCase("-of")|| args[i].equalsIgnoreCase("--onlyfeature")) {
				if (i+1== args.length) {
					System.out.println("Error: provide overlap feature from GTF.");
					System.exit(-1);
				}
				asta.setOvlFeature(args[i+1]);
				if (args[i].equalsIgnoreCase("-of")|| args[i].equalsIgnoreCase("--onlyfeature"))
					asta.setOvlOnly(true);
				else if (args[i].equalsIgnoreCase("-nf")|| args[i].equalsIgnoreCase("--notfeature"))
					asta.setOvlExcluded(true);
				++i;
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-g")|| args[i].equalsIgnoreCase("--genome")) {
				if (i+1>= args.length) {
					System.err.println("no par "+args[i]);
					continue;
				}
				
				int p= args[i+1].indexOf("_");	// i would say mouse_mm8 ok?
				if (p>= 0) {
					asta.setSpeStr(args[i+1].substring(0, p));
					asta.setGenomeVer(args[i+1].substring(p+1, args[i+1].length()));
				} else
					System.out.println("ERROR: wrong format for genome "+args[i+1]);
					//annStr= args[i+1];	// not tested, check for mapping to UCSC string
				++i;
				continue;
			}
			
			if (args[i].equalsIgnoreCase("--gtag")) {
				asta.setFiltGTAG(true);
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-h")|| args[i].equalsIgnoreCase("--html")) {
				asta.setWriteHTML(true);
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-l")|| args[i].equalsIgnoreCase("--filter")) {
				asta.addFilter(args[++i], args[++i]);
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-o")|| args[i].equalsIgnoreCase("--output")) {
				if (i+1>= args.length) {
					System.err.println("-o par not recognized");
					continue;
				}
				gphase.tools.File outDir= new gphase.tools.File(args[++i]);
				if (!outDir.exists())
					System.out.println("WARNING: output folder does not exist, ignored "+outDir.getAbsolutePath());
				asta.setOutDir(outDir);
				continue;
			}			
			
			if (args[i].equalsIgnoreCase("-r")|| args[i].equalsIgnoreCase("--reference")) {
				File f= new File(args[i+1]);
				if (!f.exists())
					System.out.println("problems with reference trpt file: "+args[i+1]);
				else {
					TabDelimitedFormatWrapper reader= new TabDelimitedFormatWrapper(f.getAbsolutePath());
					try {
						reader.read();
					} catch (Exception e) {
						e.printStackTrace();
					}
					String[] ids= reader.getColumn(0);
					asta.setRefTrptIDs(ids);
				}
				++i;
				continue;
			}
			
			if (args[i].equalsIgnoreCase("-w")|| args[i].equalsIgnoreCase("--write")) {
				if (i+1>= args.length) {
					System.err.println("provide output format");
					continue;
				}
				if (args[i+1].equalsIgnoreCase("asta"))  
					asta.setWriteASTA(true);
				else if (args[i+1].equalsIgnoreCase("gtf"))
					asta.setWriteGTF(true);
				else if (!args[i+1].equalsIgnoreCase("none"))
					System.err.println("-w parameter not recognized"+args[++i]);
				continue;
			}
			
				// let this first check else all double arg flags have to continue
	
			if (args[i].equalsIgnoreCase("-codingTranscripts"))
				asta.setCodingTranscripts(true);
			if (args[i].equalsIgnoreCase("-noncodTranscripts"))
				asta.setNoncodTranscripts(true);
			if (args[i].equalsIgnoreCase("-nonmd"))
				asta.setNmd(false);
			if (args[i].equalsIgnoreCase("-nmd"))
				asta.setNmd(true);
			
//	
//			if (args[i].equalsIgnoreCase("-noFilter"))
//				filterCode= ASMultiVariation.FILTER_NONE;
//			if (args[i].equalsIgnoreCase("-structFilter"))
//				filterCode= ASMultiVariation.FILTER_STRUCTURALLY;
//			
//			if (args[i].equalsIgnoreCase("-CDS"))
//				codingCode= ASVariation.TYPE_CDS;
//			if (args[i].equalsIgnoreCase("-UTR"))
//				codingCode= ASVariation.TYPE_UTR;
//			if (args[i].equalsIgnoreCase("-5UTR"))
//				codingCode= ASVariation.TYPE_5UTR;
//			if (args[i].equalsIgnoreCase("-3UTR"))
//				codingCode= ASVariation.TYPE_3UTR;
//	
//			
//			if (args[i].equalsIgnoreCase("-filtGTAG"))
//				filtGTAG= true;
//			
//			if (args[i].equalsIgnoreCase("-species")) {
//				if (i+1>= args.length) {
//					System.err.println("provide species name");
//					break;
//				}
//				String[] spe= Species.SP_NAMES_COMMON;
//				speStr= null;
//				for (int j = 0; j < spe.length; j++) 
//					if (spe[j].equalsIgnoreCase(args[i+1])) {
//						speStr= Species.SP_UCSC_CGI_STRINGS[j];
//						++i;
//						break;
//					}
//				if (speStr== null) {
//					System.err.println("Unknown species: "+args[i+1]);
//					System.err.print("Enter one of the following: ");
//					for (int j = 0; j < spe.length; j++) 
//						System.err.print(spe[j]+" ");
//				}
//				
//			}
//			
//			if (args[i].equalsIgnoreCase("-annotation")) {
//				if (i+1>= args.length) {
//					System.err.println("provide annotation name");
//					break;
//				}
//				
//			}
//	
//		
//		if (args[i].equalsIgnoreCase("-killError")) {
//			killError= true;
//		}
		
			
	}
	}
	
	public AStaLaVista() {
	}
	
	public int run() {
		if (inFile== null) {
			System.out.println("FATAL: no input file, giving up");
			return -1;
		}
		
		if (outDir!= null) {
			String[] s= outDir.list();
			if (s.length> 0&& writeHTML) {
				System.out.print("WARNING: output dir not empty, cleaning files with digit names..");
				int ctr= 0;
				for (int i = 0; i < s.length; i++) {
					gphase.tools.File ff= new gphase.tools.File(outDir+File.separator+s[i]);
					try {
						Integer.parseInt(ff.getFileNameWithoutExtension());
						ff.delete();
						++ctr;
					} catch (NumberFormatException e) {
						; //:)
					}
				}
				System.out.println("done, removed "+ctr+" files.");
			}
		}
		
		// force nonmd
//		if (nmd) {
//			nmd= false;
//		}
		
		GTFChrReader reader= new GTFChrReader(inFile.getAbsolutePath());
		if (!reader.isApplicable()) {
			System.err.println("no strict GTF file");
			reader.reformatFile();
			System.err.println("no strict GTF file");
			System.exit(-1);
		}
		if (ovlFeature!= null)
			reader.addReadFeature(ovlFeature);
		try {
			reader.read();
		} catch (Exception e1) {
			e1.printStackTrace();
		}
		Gene[] g= reader.getGenes();
		FilterFactory filtGene= (FilterFactory) filterMap.get(Gene.class.getName());
		FilterFactory filtTrpt= (FilterFactory) filterMap.get(Transcript.class.getName());
		FilterFactory filtEvent= (FilterFactory) filterMap.get(ASVariation.class.getName());
		
		HashMap mapVars= new HashMap();
		Species spe= new Species(speStr);
		spe.setGenomeVersion(genomeVer);
		int trptRemoveGenes= 0;
		while (g!= null) {

			for (int i = 0; i < g.length; i++) 
				g[i].setSpecies(spe);		// for reading seqs
			if (filtGene!= null)
				g= (Gene[]) filtGene.filter(g);
			
			
//			if (codingTranscripts)
//				g= ASManager.filterNonCodingTranscripts(g);
//			if (noncodTranscripts)
//				g= ASManager.filterCodingTranscripts(g);
//			if (nmd)
//				g= ASManager.filterNMDTranscripts(g);	// TODO mark them instead 
			
			// get vars
			HashMap mapTrptIDs= new HashMap();
			for (int i = 0; refTrptIDs!= null&& i < refTrptIDs.length; i++) 
				mapTrptIDs.put(refTrptIDs[i], refTrptIDs[i]);
			
				// genewise
			for (int i = 0; i < g.length; i++) {
				Transcript[] trpt= g[i].getTranscripts();
				if (filtTrpt!= null)
					trpt= (Transcript[]) filtTrpt.filter(trpt);
				if (trpt== null|| trpt.length== 0) {
					++trptRemoveGenes;
					continue;
				}
				ASVariation[][][] vars= ASManager.getASVariations(trpt, mapTrptIDs);
				if (vars== null)
					continue;
				if (filtEvent!= null) {
					HashMap redEventsMap= new HashMap();
					for (int j = 0; j < vars.length; j++) 
						for (int k = 0; k < vars[j].length; k++) 
							redEventsMap.put(vars[j][k][0], vars[j][k]);
					ASVariation[] filtVars= (ASVariation[]) Arrays.toField(filtEvent.filter(redEventsMap.keySet().toArray()));
					if (filtVars== null)
						continue;
					Vector vv= new Vector();
					for (int j = 0; j < filtVars.length; j++) {
						if (redEventsMap.get(filtVars[j])!= null)
							vv.add(redEventsMap.get(filtVars[j]));
					}
					vars= ASManager.clusterStructuralEqualEvents((ASVariation[][]) Arrays.toField(vv));
				}
				if (vars== null|| vars.length== 0) 
					continue;
								
				if (ovlFeature!= null) {
					if (writeHTML) {
						Arrays.sortNDFieldRev(vars);
						writePictures(vars, outDir.getAbsolutePath());
					}
					overlap(vars, ovlFeature);
				}
				if (vars== null|| vars.length== 0) 
					continue;
				

				for (int x= 0; x < vars.length; x++) {		// add to result
					Vector v= (Vector) mapVars.remove(vars[x][0][0].toStringStructureCode());
					if (v== null)
						v= new Vector();
					for (int j = 0; j < vars[x].length; j++) 
						v.add(vars[x][j]);
					mapVars.put(vars[x][0][0].toStringStructureCode(), v);
				}
			}
			
				// 
//			ASVariation[][][] vars= ASManager.getASVariations(g, mapTrptIDs);	// , filterCode
//			if (ovlFeature!= null) {
//				if (writeHTML) {
//					Arrays.sortNDFieldRev(vars);
//					writePictures(vars, outDir.getAbsolutePath());
//				}
//				overlap(vars, ovlFeature);
//			}
//			if (vars== null) {
//				try {
//					reader.read();
//				} catch (Exception e1) {
//					e1.printStackTrace();
//				}
//				g= reader.getGenes();
//				continue;
//			}
//			
//			
//			
//			
//			if (filtGTAG) {
//				int before= ASManager.countDifferentEvents(vars);
//				vars= ASManager.filterNonGTAGintrons(vars);
//				int after= ASManager.countDifferentEvents(vars);
//				
//				String percStr= Float.toString(((before- after)* 100f)/ before);
//				percStr= percStr.substring(0, percStr.indexOf('.')+ 2);
//				System.out.println("filtered out "+(before-after)+" ("+percStr+"%) involving nonGT/AG introns.");
//			}
//			
//			if (ovlOnly) {
//				int before= ASManager.countDifferentEvents(vars);
//				Method m= null;
//				try {
//					m= ASVariationWithRegions.class.getMethod("isASVariationWithRegions", ASVariation.class);
//				} catch (Exception e) {
//					e.printStackTrace();
//				} 
//				vars= ASManager.filter(vars, m);
//				int after= ASManager.countDifferentEvents(vars);
//				
//				String percStr= Float.toString(((before- after)* 100f)/ before);
//				percStr= percStr.substring(0, percStr.indexOf('.')+ 2);
//				System.out.println("filtered out "+(before-after)+" ("+percStr+"%) not overlapping with a "+ovlFeature+".");
//			} else if (ovlExcluded) {
//				int before= ASManager.countDifferentEvents(vars);
//				Method m= null;
//				try {
//					m= ASVariationWithRegions.class.getMethod("isNotASVariationWithRegions", ASVariation.class);
//				} catch (Exception e) {
//					e.printStackTrace();
//				} 
//				vars= ASManager.filter(vars, m);
//				int after= ASManager.countDifferentEvents(vars);
//				
//				String percStr= Float.toString(((before- after)* 100f)/ before);
//				percStr= percStr.substring(0, percStr.indexOf('.')+ 2);
//				System.out.println("filtered out "+(before-after)+" ("+percStr+"%) overlapping with a "+ovlFeature+".");
//			}
//			
//			
//			if (vars!= null&& vars.length!= 0) {
//				try {
//					vars= ASVariation.filter(vars, ASVariation.class.getMethod(evFilterStr, null));
//				} catch (Exception e) {
//					;//e.printStackTrace();
//				}
//				for (int i = 0; i < vars.length; i++) {		// add to result
//					Vector v= (Vector) mapVars.remove(vars[i][0][0].toStringStructureCode());
//					if (v== null)
//						v= new Vector();
//					for (int j = 0; j < vars[i].length; j++) 
//						v.add(vars[i][j]);
//					mapVars.put(vars[i][0][0].toStringStructureCode(), v);
//				}
//			}
			
			try {
				reader.read();
			} catch (Exception e1) {
				e1.printStackTrace();
			}
			g= reader.getGenes();
			
		}
		
		if (filtGene!= null|| filtTrpt!= null|| filtEvent!= null)
			System.out.println("Filters");
		if (filtGene!= null)
			System.out.println("Gene level: "+filtGene.toStringStats());
		if (filtTrpt!= null) {
			System.out.println("Transcript level: "+filtTrpt.toStringStats());
			System.out.println("\tremoved "+trptRemoveGenes+ " genes.");
		}
		if (filtEvent!= null)
			System.out.println("Event level: "+filtEvent.toStringStats());
		
		ASVariation[][][] vars= (ASVariation[][][]) Arrays.toField(mapVars.values());
		vars= (ASVariation[][][]) Arrays.sortNDFieldRev(vars);
		ASManager.printStats(vars, System.out);
		
		if (writeASTA) 
			writeASTA(vars);
		
		if (writeGTF) 
			writeGTF(vars);
		
		if (writeHTML)
			writeHTML(vars);
		if (writeLandscape) {
			String fName= inFile.getPathOnly();
			PrintStream p= System.out;
			if (outDir!= null) {
				fName= outDir.getAbsolutePath();
				fName+= File.separator+ inFile.getFileNameWithoutExtension();
				if (filtGene!= null) {
					fName+= filtGene.getC().getName();
					for (int i = 0; i < filtGene.getMNames().length; i++) 
						fName+= "_"+ filtGene.getMNames()[i];
				}
				if (filtTrpt!= null) {
					fName+= filtTrpt.getC().getName();
					for (int i = 0; i < filtTrpt.getMNames().length; i++) 
						fName+= "_"+ filtTrpt.getMNames()[i];
					
				}
				if (filtEvent!= null) {
					fName+= filtEvent.getC().getName();
					for (int i = 0; i < filtEvent.getMNames().length; i++) 
						fName+= "_"+ filtEvent.getMNames()[i];
				}
				fName+= ".landscape";
				try {
					p= new PrintStream(fName);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			}
			ASManager.printLandscape(vars, p, null);	// "1* , 0"
		}
		
		long t1= System.currentTimeMillis();
		return 0;
	}
	
	public static void main(String[] args) {

		//
		long t0= System.currentTimeMillis();
		AStaLaVista astalavista= new AStaLaVista();
		parseArguments(args, astalavista);
		astalavista.run();
		
		
		long t1= System.currentTimeMillis();
		System.out.println("time: "+ (t1-t0)+"[msec]");
	}
	
	public void writeASTA(ASVariation[][][] vars) {
		PrintStream p= null;
		try {
			p= new PrintStream(outDir+File.separator+ "landscape.asta");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			return;
		}
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			for (int j = 0; j < vars[i].length; j++) 
				//for (int k = 0; k < vars[i][j].length; k++) 
					p.println(vars[i][j][0].toStringASTA());
		p.flush(); p.close();
	}

	public void writeGTF(ASVariation[][][] vars) {
		Vector v= new Vector();
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			for (int j = 0; j < vars[i].length; j++) {
				for (int k = 0; k < vars[i][j].length; k++) {
					
					GTFObject o= new GTFObject();
					o.setSeqname(vars[i][j][k].getGene().getChromosome());				
					o.setSource("astalavista");
					o.setFeature("as_event");
					SpliceSite[] su= vars[i][j][k].getSpliceUniverse();
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
					o.setStrand(vars[i][j][k].getGene().getStrand());
					// frame
					o.addAttribute("as_code", vars[i][j][k].toString());
					o.addAttribute("transcript1_id", vars[i][j][k].getTranscript1().toString());
					String sc= "";
					for (int m = 0; m < vars[i][j][k].getSpliceChain1().length; m++) 
						sc+= Math.abs(vars[i][j][k].getSpliceChain1()[m].getPos())+",";
					if (sc.length()> 0)
						sc= sc.substring(0, sc.length()- 1);
					o.addAttribute("splice_chain1", sc);
					
					o.addAttribute("transcript2_id", vars[i][j][k].getTranscript2().toString());
					sc= "";
					for (int m = 0; m < vars[i][j][k].getSpliceChain2().length; m++) 
						sc+= Math.abs(vars[i][j][k].getSpliceChain2()[m].getPos())+",";
					if (sc.length()> 0)
						sc= sc.substring(0, sc.length()- 1);
					o.addAttribute("splice_chain2", sc);
					
					v.add(o);
				}
			}
		
		GTFWrapper gtf= new GTFWrapper(outDir+ File.separator+ "landscape.gtf");
		gtf.setGtfObj((GTFObject[]) Arrays.toField(v));
		gtf.setSortAttributes(new String[] {"as_code", "transcript1_id", "splice_chain1", "transcript2_id", "splice_chain2"});
		try {
			gtf.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public void writeHTML(ASVariation[][][] vars) {
		
		PrintStream p;
		try {
			if (ovlFeature!= null) 
				mapRegions= new HashMap();
			ASManager.sortLandscape(vars);
			writeHTMLlandscapePie(vars);
			
			p= new PrintStream(outDir.getAbsolutePath()+File.separator+FNAME_LANDSCAPE);
			writeHTMLlandscape(vars, p);
			p.flush(); p.close();
			
			writePicturesWithRegions(vars, outDir.getAbsolutePath());	// do before
			for (int i = 0; vars!= null&& i < vars.length; i++) {
				String varFStr= getFileName(vars[i][0][0].toStringStructureCode())+ ".html";	//varStr.replace(' ', '_')+ ;
				
				p= new PrintStream(outDir.getAbsolutePath()+File.separator + varFStr);
				writeHTMLeventGroup(vars[i], p);
				p.flush(); p.close();
			}
			if (ovlFeature!= null) {
				writeHTMLovlSummary();
				SpliceOSigner.writeOutDomainColorMap();
			}
			
		} catch (Throwable e) {
			e.printStackTrace();
		}
	}


	protected void writeHTMLovlLinkbacks(PrintStream p) {
		p.println("<a href=\""+HTML_NAME_LANDSCAPE+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
		"<FONT face=\"Arial,Lucida Sans\" size=\"1\">landscape</FONT></a>");
	}
	
	protected void writeHTMLovlSummary() {
		PrintStream p= null;
		try {
			p= new PrintStream(outDir.getAbsolutePath()+File.separator+FNAME_REGIONS);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		include(p, HEADER_FILE);
		headline(p, ovlFeature+" summary", null);
		p.println("<DIV>");
		writeHTMLovlLinkbacks(p);
		p.print("</DIV>");

//		p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>"+ovlFeature+" ID</b></FONT></TD>" +
//				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
//		"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event Count</b></FONT>" +
//				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
//		"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Events</b></FONT></TR>");

		p.print("<TABLE align=\"center\" width=\"100%\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">\n");	// width=\"100%\"
		Object[] regIDs= mapRegions.keySet().toArray();
		java.util.Arrays.sort(regIDs);
		for (int i = 0; i < regIDs.length; i++) {
			String stripID= ASVariationWithRegions.stripOrderSuffix((String) regIDs[i]);
			p.print("<TR>\n");
			Vector v= (Vector) mapRegions.get(regIDs[i]);
			p.print("\t<TR><TD><a name=\""+stripID+"\"><img src=\""+writeDotPic(SpliceOSigner.getDomainColor((String) regIDs[i]))+"\">&nbsp;"+
					stripID+"</a></TD><TD></TD>\n");
			p.print("\t<TD>"+v.size()+"</TD><TD></TD>\n");
			p.print("\t<TD>");
			for (int j = 0; j < v.size(); j++) 
				p.print(v.elementAt(j)+"&nbsp;");
			p.print("</TD>\n");
		}
		 
		p.print("</TABLE>\n");
		p.print("<br><br>");
		p.print("</div><!-- closing main -->\n");
		
		include(p, TRAILER_FILE);
		p.flush(); p.close();
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



	void writeHTMLlandscapePie(ASVariation[][][] vars) {

			// init pie
		Pie pie= new Pie();
		pie.setSize(new Dimension(300,300));
		pie.init();
		
		int sum= 0;
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			if (vars[i][0][0].toStringStructureCode().equals("1-2^ , 0"))
				pie.AddPieSlice(vars[i].length, "exon skipping", EXSKIP_COL);
//			else if (vars[i][0].toString().equals("1-2^3-4^ , 0"))
//				pie.AddPieSlice(vars[i].length, "double skipping", DBLSKIP_COL);
			else if (vars[i][0][0].toStringStructureCode().equals("1- , 2-"))
				pie.AddPieSlice(vars[i].length, "alt acceptor", AACC_COL);
			else if (vars[i][0][0].toStringStructureCode().equals("1^ , 2^"))
				pie.AddPieSlice(vars[i].length, "alt donor", ADON_COL);
			else if (vars[i][0][0].toStringStructureCode().equals("1^2- , 0"))
				pie.AddPieSlice(vars[i].length, "intron retention", IR_COL);
			else if (vars[i][0][0].toStringStructureCode().equals("1-2^ , 3-4^"))
				pie.AddPieSlice(vars[i].length, "mutually exclusive", ME_COL);
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
		      File f= new File(outDir.getAbsolutePath()+File.separator+FNAME_PIE);
			  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("png").get(0);
			  pie.setPreferredSize(pie.getPreferredSize());
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



	void writePictures(ASVariation[][][] vars, String path) {
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			try {
			      String varStr= vars[i][0][0].toStringStructureCode();
			      String varFStr= getFileName(varStr)+ ".gif";		//varStr.replace(' ', '_');
			      File f= new File(path+File.separator+varFStr);
			      if (f.exists())
			    	  continue;
				  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			      Component component= new SpliceOSigner(vars[i][0][0]);
			      component.setSize(component.getPreferredSize());
			      t.exportToFile(f,component,component.getParent(),null,null);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

 

	void writePicturesWithRegions(ASVariation[][][] vars, String path) {
		for (int i = 0; vars!= null&& i < vars.length; i++) 
			for (int j = 0; j < vars[i].length; j++) 
				for (int k = 0; k < vars[i][j].length; k++) 
			
			try {
			      String varStr= vars[i][j][k].toString();
			      String varFStr= getFileName(varStr)+ ".gif";		//varStr.replace(' ', '_');
			      File f= new File(path+File.separator+varFStr);
			      if (f.exists())
			    	  continue;
				  ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			      Component component= new SpliceOSigner(vars[i][j][k]);
			      component.setSize(component.getPreferredSize());
			      t.exportToFile(f,component,component.getParent(),null,null);
			} catch (Exception e) {
				e.printStackTrace();
			}
		
	}
	
	static void headline(PrintStream p, String headline, String linkBack) {
//		<TD class="section" border=0 cellpadding=0 cellspacing=0 width=20% ALIGN=RIGHT>
//		<a href="index.html#TOP"
//		 onMouseover="window.status='TOP:INDEX';">
//		<IMG class="pnt" SRC="http://genome.imim.es/g_icons/top.gif" HEIGHT=15 WIDTH=15 BORDER=0></a>
		p.println("<TABLE border=0 cellpadding=0 cellspacing=0 width=\"100%\">");
		p.println("<TR border=0 cellpadding=0 cellspacing=0>");
		if (linkBack!= null)
			p.println("<TD class=\"section\" border=0 cellpadding=0 cellspacing=0 width=\"20%\"><a href=\""+linkBack+"\"><IMG class=\"pnt\" SRC=\"http://genome.imim.es/g_icons/top.gif\" HEIGHT=15 WIDTH=15 BORDER=0></a></TD>");
		p.println("<TD class=\"section\" align=\"center\">");
		p.println("<CENTER><FONT size=6 class=\"tgen\">"+headline+"</FONT></CENTER></TD>");
		p.println("</TD></TR></TABLE>");
	}
	
	void writeHTMLlandscapeLinkouts(PrintStream p) {
		if (writeASTA)
			p.println("<IMG src= \"http://genome.imim.es/astalavista/pics/dl_arrow.jpg\">"+
					"<A HREF=\"landscape.asta\">Download output (ASTA format selected)</A><br><br></CENTER></DIV>");
		if (writeGTF)
			 p.println("<IMG src= \"http://genome.imim.es/astalavista/pics/dl_arrow.jpg\">"+
					"<A HREF=\"landscape.gtf\">Download output (GTF format selected)</A><br><br></CENTER></DIV>");
	}
	
	void writeHTMLlandscapeLinkbacks(PrintStream p) {
		if (ovlFeature!= null) {
			p.println("<a href=\""+FNAME_REGIONS+"\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\">" +
					"<FONT face=\"Arial,Lucida Sans\" size=\"1\">"+ovlFeature+"</FONT></a>");
		}
	}
	
	void writeHTMLlandscape(ASVariation[][][] vars, PrintStream p) {
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
			p.print("<DIV>\n");
			writeHTMLlandscapeLinkbacks(p);
			p.print("</DIV>");
			p.println("<DIV align=\"center\"><CENTER><br><img src=\"distribution.png\"><br><br>");

			writeHTMLlandscapeLinkouts(p);
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
				String varFStr= getFileName(vars[i][0][0].toStringStructureCode());	//varStr.replace(' ', '_');
				p.println("\t<TD align=\"center\" valign=\"center\">"
						+ "<FONT face=\"Arial,Lucida Sans\" size=\"2\"><a href=\""+ varFStr+".html\"></FONT>Show>></a>");
				p.println("\t<TD></TD>"); 
				p.println("\t<TD align=\"left\" valign=\"center\">"
						//+ "<a href=\""+ varFStr+".html\">"
						+ "<img src=\""+ varFStr+".gif\"><br>"
						// </a>
						+ "<FONT size=\"2\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ vars[i][0][0].toStringStructureCode()+ "</b></FONT>"
						+"</TD>");
				p.println("</TR>");
			}
			p.println("</TABLE></DIV>");
		}
		p.print("<br><br>");
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

	protected static String getFileName(String in) {
		Integer name= (Integer) fileNameMap.get(in);
		if (name== null) {
			name= new Integer(++eventTypeNr);
			fileNameMap.put(in, name);
		}
		return name.toString();
	}



	protected static String convertFName(String in) {
		StringBuffer sb= new StringBuffer(in);
			// kill spaces
		for (int i = 0; i < sb.length(); i++) 
			if (sb.charAt(i)== ' ')
				sb.deleteCharAt(i--);
		String out= sb.toString();
		out= out.replace('^', 'd');
		out= out.replace('-', 'a');
		out= out.replace(',', 'I');
		Integer name= (Integer) fileNameMap.get(out);
		if (name== null) {
			name= new Integer(++eventTypeNr);
			fileNameMap.put(out, name); 
		}
		return name.toString();
//		if (out.length()> 245)	{
//			System.out.println("WARNING: shorted filename "+out);
//			out= out.substring(0, 245);	// max fname len
//		}
//		
//		return out;
	}
	
	protected void writeHTMLeventGroupColRefPair(ASVariation[] vars, PrintStream p) {
		
		p.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+ 
				vars[0].getTranscript1()+" "+vars[0].getTranscript2()
				+"</b><br>\n");
		String fName= getFileName(vars[0].toCoordinates())+".html";
		p.print("<a href=\""+fName+"\">all "+vars.length+" transcript pairs</a></TD>\n");
		
			// event details
		try {
			PrintStream pp= new PrintStream(outDir.getAbsolutePath()+File.separator+fName);
			writeHTMLevent(vars, pp);
			pp.flush(); pp.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	protected void writeHTMLeventColTPair(ASVariation var, PrintStream p) {
		p.print("\t<TD valign=\"top\"<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+ 
				var.getTranscript1()+" "+var.getTranscript2()
				+"</b></TD>\n");
	}
	
	protected void writeHTMLevent(ASVariation[] vars, PrintStream p) {

		java.util.Arrays.sort(vars, new ASVariation.SpliceStringComparator());
		include(p, HEADER_FILE);
		headline(p, "Event Details", null);

		// splice chains
		String pos1Str= ""; 
		SpliceSite[] su= vars[0].getSpliceUniverse();
		for (int j = 0; j < vars[0].getSpliceChain1().length; j++)
			if (vars[0].getSpliceChain1()[j].getPos()> 0)
				pos1Str+= vars[0].getSpliceChain1()[j].getPos()+ " ";
			else
				pos1Str+= Math.abs(vars[0].getSpliceChain1()[j].getPos())+ " ";
		String pos2Str= "";
		for (int j = 0; j < vars[0].getSpliceChain2().length; j++) 
			if (vars[0].getSpliceChain2()[j].getPos()> 0)
				pos2Str+= vars[0].getSpliceChain2()[j].getPos()+ " "; 
			else
				pos2Str+= Math.abs(vars[0].getSpliceChain2()[j].getPos())+ " ";
		
			// position
		String s= vars[0].getGene().getChromosome()+":"+pos1Str+","+pos2Str;
		String varFStr= getFileName(vars[0].toStringStructureCode());	//varStr.replace(' ', '_');
		
		
		p.println("<DIV>");
		p.println("<a href=\""+varFStr+".html\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Structure Group\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">Structure Group</FONT></a></DIV>");
		p.println("<DIV align=\"center\">");	// class=\"userspace\"  makes it white
		p.println("<CENTER><img src=\""+ varFStr+".gif\"><br>" +
				"<FONT size=\"2\" face=\"Courier\"><b>"+ s+ "</b></FONT></CENTER><br><br>");
		p.println("</DIV>");
		
		p.print("<DIV align=\"center\">");
		p.print("<TABLE align=\"center\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">\n");
		writeHTMLeventTableHeader(p);
		HashMap mapOvlStruct= new HashMap();
		for (int i = 0; i < vars.length; i++) {
			String colStr= "";
			if (i%2 == 0)
				colStr= TABLE_EVEN_COLOR;
			else
				colStr= TABLE_ODD_COLOR;
			p.println("<TR bgcolor=\""+colStr+"\">");
			
				// counter
			p.println("\t<a name=\""+ (i+1)+"\"></a><TD valign=\"middle\" align=\"center\">"
					+ (i+1)+ "</TD>");
			p.println("\t<TD></TD>"); 
			
			writeHTMLeventColTPair(vars[i], p);
			p.println("\t<TD></TD>"); 
			if (ovlFeature!= null) {
				writeHTMLeventColSchema(vars[i],p);
				if (vars[i] instanceof ASVariationWithRegions&& ((ASVariationWithRegions) vars[i]).getRegions().length> 0&&
						(mapOvlStruct.get(((ASVariationWithRegions) vars[i]).toStringRegionColorCode())== null)) {
					mapOvlStruct.put(((ASVariationWithRegions) vars[i]).toStringRegionColorCode(),
							((ASVariationWithRegions) vars[i]).toStringRegionColorCode());
					String[] ids= ((ASVariationWithRegions) vars[i]).getRegionIDsNonRedundant();
					for (int j = 0; j < ids.length; j++) {
						String stripID= ASVariationWithRegions.stripOrderSuffix(ids[j]);
						Vector v= (Vector) mapRegions.remove(ids[j]);
						if (v== null)
							v= new Vector();
						String fName= getFileName(vars[0].toCoordinates());
						v.add("<a href=\""+fName+".html#"+(i+1)+"\">" +
								"<img src=\""+getFileName(vars[i].toString())+".gif\"></a>");
						mapRegions.put(ids[j], v);
					}
				}
				p.println("\t<TD></TD>"); 
			}
			writeHTMLeventGroupColUCSC(vars[i], p);
			p.print("</TR>\n");
		}
		p.print("</TABLE>\n");
		p.print("<br><br>");
		p.print("</DIV>");
		p.print("</div><!-- closing main -->\n");
		
		include(p, TRAILER_FILE);
		 
	}
	
	protected void writeHTMLeventGroupTableHeader(PrintStream p) {
		p.print("<TR><TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"center\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Event<br>Nr</b></FONT></TD>" +
						"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\" align=\"right\">&nbsp&nbsp</TD>"+
				"<TD bgcolor=\""+TABLE_HEADER_COLOR+"\"><FONT face=\"Arial,Lucida Sans\" size=\"4\" color=\"#FFFFFF\"><b>Transcript Pairs</b></FONT>" +
						"<br><FONT face=\"Arial,Lucida Sans\" size=\"2\" color=\"#FFFFFF\"><b>Reference Pair</b> <u>all Pairs</u></FONT></TD>"+
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
		} else
			p.println("<TD valign=\"middle\" align=\"center\" width=\"120\"></TD>");	// empty
		p.println("\t<TD></TD>"); 

	}
	
	protected void writeHTMLeventGroupColSchema(ASVariation[] vars, PrintStream p) {
		p.print("<TD align=\"left\" valign=\"center\">");
		HashMap mapReg= new HashMap();
		for (int i = 0; i < vars.length; i++) {
			ASVariation var= vars[i];
			if (var instanceof ASVariationWithRegions) {
				String[] regIDs= ((ASVariationWithRegions) var).getRegionIDsNonRedundant();
				for (int j = 0; j < regIDs.length; j++) {
					String id= regIDs[j];
					String idStrip= ASVariationWithRegions.stripOrderSuffix(id); 
					if (mapReg.get(idStrip)!= null)
						continue;
					p.print("<br><img src=\""+writeDotPic(SpliceOSigner.getDomainColor(id))+"\">" +
							"&nbsp;<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b><a href=\""+
							FNAME_REGIONS+"#"+idStrip+"\">"+idStrip+"</b></FONT>\n");
					mapReg.put(idStrip, idStrip);
				} 
			}
		}
		p.print("</TD>\n");
	}



	protected void writeHTMLeventColSchema(ASVariation var, PrintStream p) {
		String s= var.toString();
		String f= getFileName(var.toString());
		p.print("<TD align=\"left\" valign=\"center\"><img src=\""+getFileName(var.toString())+".gif\">\n");
		if (var instanceof ASVariationWithRegions) {
			DirectedRegion[] regs= ((ASVariationWithRegions) var).getRegions();
			HashMap mapReg= new HashMap();
			for (int i = 0; i < regs.length; i++) {
				String id= regs[i].getID();
				String idStrip= ASVariationWithRegions.stripOrderSuffix(id); 
				if (mapReg.get(idStrip)!= null)
					continue;
				p.print("<br><img src=\""+writeDotPic(SpliceOSigner.getDomainColor(id))+"\">" +
						"&nbsp;<FONT face=\"Arial,Lucida Sans\" size=\"2\"><b>"+
						"<a href=\""+FNAME_REGIONS+"#"+idStrip+"\">"+idStrip+"</a></b></FONT>\n");
				mapReg.put(idStrip, idStrip);
			} 
		}
		p.print("</TD>\n");
	}
	
	protected void writeHTMLeventGroup(ASVariation[][] vars, PrintStream p) {
		
		include(p, HEADER_FILE);
		headline(p, "Event Structure", null);
		String varFStr= getFileName(vars[0][0].toStringStructureCode());	//varStr.replace(' ', '_');
		String varStr= vars[0][0].toStringStructureCode().replaceAll("\\s", "");
		
		
		p.println("<DIV>");
		p.println("<a href=\"landscape.html\"><IMG class=\"pnt\" src=\"http://genome.imim.es/g_icons/top.gif\" height=\"15\" width=\"15\" border=\"0\" alt=\"Landscape\">" +
				"<FONT face=\"Arial,Lucida Sans\" size=\"1\">Landscape</FONT></a></DIV>");
		p.println("<DIV align=\"center\">");	// class=\"userspace\"  makes it white
		p.println("<CENTER><img src=\""+ varFStr+".gif\"><br>" +
				"<FONT size=\"2\" face=\"Arial,Lucida Sans\">code: </FONT><FONT size=\"2\" face=\"Courier\"><b>"+ varStr+ "</b></FONT></CENTER><br><br>");
		p.println("</DIV>");
		
		
		p.print("<DIV align=\"center\">");
		p.print("<TABLE align=\"center\" border=\"0\" cellspacing=\"0\" cellpadding=\"0\">\n");	// width=\"100%\"
		// org=D.+melanogaster
		// org=Homo_sapiens&db=hg17
		//&clade=vertebrate&org=Mouse&db=mm8
		
		writeHTMLeventGroupTableHeader(p);
		
		for (int i = 0; i < vars.length; i++) {
			String colStr= "";
			if (i%2 == 0)
				colStr= TABLE_EVEN_COLOR;
			else
				colStr= TABLE_ODD_COLOR;
			p.println("<TR bgcolor=\""+colStr+"\">");
			
				// counter
			p.println("\t<TD valign=\"middle\" align=\"center\">"
					+ (i+1)+ "</TD>");
			p.println("\t<TD>&nbsp&nbsp</TD>"); 
			
			writeHTMLeventGroupColRefPair(vars[i],p);
			p.println("\t<TD>&nbsp&nbsp</TD>"); 
			if (ovlFeature!= null) {
				writeHTMLeventGroupColSchema(vars[i],p);
				p.println("\t<TD>&nbsp&nbsp</TD>"); 
			}
			writeHTMLeventGroupColUCSC(vars[i][0],p);
			p.print("</TR>\n");
		}
		p.print("</TABLE>\n");
		p.print("<br><br>");
		p.print("</DIV>");
		p.print("</div><!-- closing main -->\n");
		
		include(p, TRAILER_FILE);
	}

	static int eventTypeNr= 0;
	static HashMap fileNameMap= new HashMap();

	public String writeDotPic(Color c) {
		try {
			String imgID= getFileName(Integer.toHexString(c.getRGB()))+".gif";
			File f= new File(outDir.getAbsolutePath()+File.separator+imgID);
			if (f.exists())
				return imgID;
			ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
			Component component= new Circle(c,8);
		    component.setSize(component.getPreferredSize());
		    t.exportToFile(f,component,component.getParent(),null,null);
		    return imgID;
		} catch (Exception e) {
			// TODO: handle exception
		}
		return null;
	}



	public void overlap(ASVariation[][][] vars, String featureID) {
		for (int i = 0; vars!= null&& i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				for (int k = 0; k < vars[i][j].length; k++) {
					Object o1= null, o2= null;
					if (featureID.equals(GTFObject.CDS_FEATURE_TAG)) {
						o1= new DirectedRegion[] {vars[i][j][k].getTranscript1().getTranslations()[0]};
						o2= new DirectedRegion[] {vars[i][j][k].getTranscript2().getTranslations()[0]};
					} else {
						o1= vars[i][j][k].getTranscript1().getAttribute(featureID);
						o2= vars[i][j][k].getTranscript2().getAttribute(featureID);
					}
					if (o1== null&& o2== null)
						continue;
					DirectedRegion[] v1, v2;
					try {
						v1= (DirectedRegion[]) o1;
						v2= (DirectedRegion[]) o2;
					} catch (ClassCastException e) {
						continue;
					}
					ASVariationWithRegions var= new ASVariationWithRegions(
							vars[i][j][k], v1, v2); 
					if (var.getRegions()!= null&& var.getRegions().length> 0)
						vars[i][j][k]= var;
				}
			}
		}
	}



	public boolean isCodingTranscripts() {
		return codingTranscripts;
	}



	public void setCodingTranscripts(boolean codingTranscripts) {
		this.codingTranscripts = codingTranscripts;
	}



	public String getGenomeVer() {
		return genomeVer;
	}



	public void setGenomeVer(String genomeVer) {
		this.genomeVer = genomeVer;
	}



	public gphase.tools.File getInFile() {
		return inFile;
	}



	public void setInFile(gphase.tools.File inFile) {
		this.inFile = inFile;
	}



	public boolean isNoncodTranscripts() {
		return noncodTranscripts;
	}



	public void setNoncodTranscripts(boolean noncodTranscripts) {
		this.noncodTranscripts = noncodTranscripts;
	}



	public gphase.tools.File getOutDir() {
		return outDir;
	}



	public void setOutDir(gphase.tools.File outDir) {
		this.outDir = outDir;
	}



	public String getSpeStr() {
		return speStr;
	}



	public void setSpeStr(String speStr) {
		this.speStr = speStr;
	}



	public String getOvlFeature() {
		return ovlFeature;
	}



	public void setOvlFeature(String ovlFeature) {
		this.ovlFeature = ovlFeature;
	}



	public boolean isWriteHTML() {
		return writeHTML;
	}



	public void setWriteHTML(boolean writeHTML) {
		this.writeHTML = writeHTML;
	}



	public boolean isWriteGTF() {
		return writeGTF;
	}



	public void setWriteGTF(boolean writeGTF) {
		this.writeGTF = writeGTF;
	}



	public boolean isWriteASTA() {
		return writeASTA;
	}



	public void setWriteASTA(boolean writeASTA) {
		this.writeASTA = writeASTA;
	}



	public boolean isNmd() {
		return nmd;
	}



	public void setNmd(boolean nmd) {
		this.nmd = nmd;
	}



	public String[] getRefTrptIDs() {
		return refTrptIDs;
	}



	public void setRefTrptIDs(String[] refTrptIDs) {
		this.refTrptIDs = refTrptIDs;
	}



	public boolean isFiltGTAG() {
		return filtGTAG;
	}



	public void setFiltGTAG(boolean filtGTAG) {
		this.filtGTAG = filtGTAG;
	}



	public boolean isOvlExcluded() {
		return ovlExcluded;
	}



	public void setOvlExcluded(boolean ovlExcluded) {
		this.ovlExcluded = ovlExcluded;
	}



	public boolean isOvlOnly() {
		return ovlOnly;
	}



	public void setOvlOnly(boolean ovlOnly) {
		this.ovlOnly = ovlOnly;
	}
	
	public void addFilter(String className, String filterString) {
		if (filterMap== null)
			filterMap= new HashMap();
		
		FilterFactory fac= (FilterFactory) filterMap.get(className);
		if (fac== null) {
			fac= new FilterFactory(className, filterString);
			filterMap.put(className, fac);
		} else
			fac.addMethodString(filterString);
	}
}
