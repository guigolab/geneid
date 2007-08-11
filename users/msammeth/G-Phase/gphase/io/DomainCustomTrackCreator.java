package gphase.io;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Vector;

import gphase.gui.SpliceOSigner;
import gphase.io.bed.BEDobject;
import gphase.io.bed.BEDwrapper;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Transcript;
import gphase.tools.Array;
import gphase.tools.Arrays;
import gphase.tools.File;
import gphase.tools.Formatter;

public class DomainCustomTrackCreator {
	static final String errorMsg= "usage: DomainToGenomeMapper [options] <inputFile>\n\n"+
	"where\n"+
	"<inputFile>\ta GTF file with annotated domains\n"+
	"and [options] may be:\n"+
	"-out <outputFile>\tthe name of the output file (default name is inputFile.bed)"+
	"\n\nmicha, may 07";
	public static final int VISIBILITY_HIDE= 0, VISIBILITY_DENSE= 1, VISIBILITY_FULL= 2, VISIBILITY_PACK= 3, VISIBILITY_SQUISH= 4;
	public static final String BED_FILE_EXTENSION= "bed";
	
	public static void main(String[] args) {
		if (args.length< 1) {
			System.err.println(errorMsg);
			System.exit(-1);
		}
		File inFile= new File(args[args.length- 1]);
		if (!inFile.exists())
			System.err.println("Input file not valid.");
		String outFname= inFile.getPathOnly()+File.separator+
			inFile.getFileNameOnly()+ "."+ BED_FILE_EXTENSION;
		for (int i = 0; i < args.length- 2; i++) {
			if (args[i].equalsIgnoreCase("-out")&& i+1< args[i].length()) {
				outFname= args[i+1];
				++i;
				continue;
			}
			
		}
		File outFile= new File(outFname);
		File.checkForOverwrite(System.out, outFile);
		
		writeBEDfiles(inFile, outFile);
		SpliceOSigner.writeOutDomainColorMap();
	}
	
	
	File outFile= null, inFile= null;
	Transcript[] trpts= null;
	HashMap<String,DirectedRegion[]> domMap= new HashMap<String,DirectedRegion[]>();
	
	public static String getBrowserHeader(DirectedRegion position,
			String[] hide, String[] dense, String[] pack, String[] squish, String[] full) {
		String s= "";
		if (position!= null) 
			s+= "browser "+position.getChromosome()+":"+Math.abs(position.getStart())+"-"+Math.abs(position.getEnd())+"\n";
		if (hide!= null) {
			s+= "browser ";
			for (int i = 0; i < hide.length; i++) {
				s+= hide[i];
				if (i< hide.length)
					s+= " ";
			}
			s+= "\n";
		}
		if (dense!= null) {
			s+= "browser ";
			for (int i = 0; i < dense.length; i++) {
				s+= dense[i];
				if (i< dense.length)
					s+= " ";
			}
			s+= "\n";
		}
		if (pack!= null) {
			s+= "browser ";
			for (int i = 0; i < pack.length; i++) {
				s+= pack[i];
				if (i< pack.length)
					s+= " ";
			}
			s+= "\n";
		}
		if (squish!= null) {
			s+= "browser ";
			for (int i = 0; i < squish.length; i++) {
				s+= squish[i];
				if (i< squish.length)
					s+= " ";
			}
			s+= "\n";
		}
		if (full!= null) {
			s+= "browser ";
			for (int i = 0; i < full.length; i++) {
				s+= full[i];
				if (i< full.length)
					s+= " ";
			}
			s+= "\n";
		}
		return s;
	}
	public static String getTrackHeader(
			String trName, String trDescription, int trDispMode, Color trCol, boolean trItemRGB, boolean trUseScore,
			String trGroupName, int trPriority, int trOffset, String trURL, String trHtmlURL) {
		String s= "";
			s+= "track ";
				// name=<track_label> - Defines the track label that will be displayed to the left of the track in the Genome Browser window, and also the label of the track control at the bottom of the screen. The name can consist of up to 15 characters, and must be enclosed in quotes if the text contains spaces. We recommend that the track_label be restricted to alpha-numeric characters and spaces to avoid potential parsing problems. The default value is "User Track".
			if (trName!= null)
				s+= "name=\""+trName+"\" ";
		    	//description=<center_label> - Defines the center label of the track in the Genome Browser window. The description can consist of up to 60 characters, and must be enclosed in quotes if the text contains spaces. The default value is "User Supplied Track".
			if (trDescription!= null)
				s+= "description=\""+trName+"\" ";
				// visibility=<display_mode> - Defines the initial display mode of the annotation track. Values for display_mode include: 0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish. The numerical values or the words can be used, i.e. full mode may be specified by "2" or "full". The default is "1".
			if (trDispMode>= 0)
				s+= "visibility=\""+trDispMode+"\" ";
				// color=<RRR,GGG,BBB> - Defines the main color for the annotation track. The track color consists of three comma-separated RGB values from 0-255. The default value is 0,0,0 (black).
			if (trCol!= null)
				s+= "color=\""+trCol.getRed()+","+trCol.getGreen()+","+trCol.getBlue()+"\" ";
				// itemRgb=On - If this attribute is present and is set to "On", the Genome Browser will use the RGB value shown in the itemRgb field in each data line of the associated BED track to determine the display color of the data on that line.
			if (trItemRGB)
				s+= "itemRgb=\"On\" ";
				// useScore=<use_score> - If this attribute is present and is set to 1, the score field in each of the track's data lines will be used to determine the level of shading in which the data is displayed. The track will display in shades of gray unless the color attribute is set to 100,50,0 (shades of brown) or 0,60,120 (shades of blue). The default setting for useScore is "0".
			if (trUseScore)
				s+= "useScore=\"1\" ";
			else
				s+= "useScore=\"0\" ";
				// group=<group> - Defines the annotation track group in which the custom track will display in the Genome Browser window. By default, group is set to "user", which causes custom tracks to display at the top of the window.
			if (trGroupName!= null)
				s+= "group=\""+trGroupName+"\" ";
				// priority=<priority> - When the group attribute is set, defines the display position of the track relative to other tracks within the same group in the Genome Browser window. If group is not set, the priority attribute defines the track's order relative to other custom tracks displayed in the default group, "user".
			if (trPriority>= 0)
				s+= "priority=\""+trPriority+"\" ";
				// offset=<offset> - Defines a number to be added to all coordinates in the annotation track. The default is "0".
			if (trOffset> 0)
				s+= "offset=\""+trOffset+"\" ";
				// url=<external_url> - Defines a URL for an external link associated with this track. This URL will be used in the details page for the track. Any '$$' in this string this will be substituted with the item name. There is no default for this attribute.
			if (trURL!= null)
				s+= "url=\""+trURL+"\" ";
				// htmlUrl=<external_url> - Defines a URL for an HTML description page to be displayed with this track. There is no default for this attribute. A template for a standard format HTML track description is here.
			if (trHtmlURL!= null)
				s+= "url=\""+trHtmlURL+"\" ";
			
			return s;
	   }
	public static strictfp void writeBEDfiles(File inFile, File outDir) {
		if (!outDir.exists())
			outDir.mkdir();
		else {
			String[] list= outDir.list();
			for (int i = 0; i < list.length; i++) 
				new File(outDir.getAbsolutePath()+ File.separator+ list[i]).delete();
		}
		GTFChrReader reader= new GTFChrReader(inFile.getAbsolutePath());
		reader.setReadAllLines();
		reader.setReadFeatures(new String[] {"exon", "CDS", "domain"});

		reader.setChromosomeWise(true);
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		
		while (genes!= null) {
			for (int i = 0; i < genes.length; i++) {
				String geneBedName= outDir.getAbsolutePath()+ File.separator 
				+genes[i].getChromosome()+"_"+genes[i].getNameTranscript().getTranscriptID()+ ".bed";
					// browser header for whole gene
				try {
					BufferedWriter buffy= new BufferedWriter(new FileWriter(geneBedName, true));
					buffy.write(getBrowserHeader(genes[i], new String[] {"all"}, null, new String[] {"RefSeq"}, null,null));
					buffy.flush(); buffy.close();
				} catch (Exception e) {
					e.printStackTrace();
				}
				Transcript[] trpts= genes[i].getTranscripts();
				Vector<BEDobject> v= new Vector<BEDobject>();
				for (int j = 0; j < trpts.length; j++) {
					if (trpts[j].getAttributes()== null|| trpts[j].getAttributes().size()== 0)
						continue;
					
					Vector domV= new Vector();
					DirectedRegion[] doms= (DirectedRegion[]) trpts[j].getAttribute(DomainToGenomeMapper.GTF_DOMAIN_FEATURE_TAG);	// one genomic reg for each domain
					HashMap mapBEDs= new HashMap();

						// construct boundaries (start!!) of bed
					Vector exonicDomsV= new Vector();
					for (int k = 0; k < doms.length; k++) {	
						int len= doms[k].getLength();	// better: satisfies UCSC condition start+ len= end 
						int absStart= Math.abs(doms[k].getStart());
						int score= 1000;								// adapt score
						if (doms[k].getScore()> 0d)
							score= (int) Math.abs(Math.log(doms[k].getScore())* 10);
						if (score> 1000)
							score= 1000;
						if (mapBEDs.get(doms[k].getID())== null) {
							BEDobject obj= new BEDobject(trpts[j].getChromosome(),trpts[j].getStrand());
							obj.setStart(absStart);
							obj.setEnd(absStart+ len);
							obj.setScore(score);
							obj.setName(doms[k].getID()+"_"+Formatter.fprint(doms[k].getScore(),2));
							obj.setCol(SpliceOSigner.getDomainColor(doms[k].getID()));
							mapBEDs.put(doms[k].getID(), obj);
						} else {
							BEDobject obj= (BEDobject) mapBEDs.get(doms[k].getID());
							if (absStart< obj.getStart())
								obj.setStart(absStart);
							if (absStart+ len> obj.getEnd())
								obj.setEnd(absStart+ len);
							if (obj.getScore()!= score)
								System.out.println("WARNING: score not matching "+trpts[j]+" "+
									doms[k].getID()+" "+doms[k].getScore()+" != "+obj.getScore());
						}
						
							// tokenize to exonic Regions only
						DirectedRegion[] exRegs= DirectedRegion.intersect(new DirectedRegion[] {doms[k]},
								trpts[j].getExons());
						for (int m = 0; m < exRegs.length; m++) {
							exRegs[m].setID(doms[k].getID());
							exRegs[m].setChromosome(trpts[j].getChromosome());
							exRegs[m].setScore(doms[k].getScore());
							exonicDomsV.add(exRegs[m]);
						}
					}
					
					
						// now add blocks
					doms= (DirectedRegion[]) Arrays.toField(exonicDomsV);
					for (int k = 0; k < doms.length; k++) {	// each domain a line in the track of the transcript
						BEDobject obj= (BEDobject) mapBEDs.get(doms[k].getID());
						obj.addBlockStart(Math.abs(doms[k].getStart())- obj.getStart(), doms[k].getLength());							
						obj.setBlockCount(obj.getBlockCount()+ 1);
					}

					Object[] vals= mapBEDs.values().toArray();
					for (int k = 0; k < vals.length; k++) {
						BEDobject obj= (BEDobject) vals[k];
						if (trpts[j].getTranslations()!= null) {
							int cdsStart= Math.abs(trpts[j].getTranslations()[0].getStart());
							int cdsEnd= Math.abs(trpts[j].getTranslations()[0].getEnd());
							int thickStart= cdsStart, thickEnd= cdsEnd;
							if (cdsStart< obj.getStart())
								thickStart= obj.getStart();
							if (cdsEnd> obj.getEnd())
								thickEnd= obj.getEnd();
							obj.setThickStart(thickStart);
							obj.setThickEnd(thickEnd);
						}
						domV.add(obj);
					}
					// end of domains
					
					
					// transcript based tracks
					BEDobject[] beds=(BEDobject[]) Arrays.toField(domV);
					try {
						BufferedWriter buffy= new BufferedWriter(new FileWriter(geneBedName, true));
						buffy.write("\n\n");
						buffy.write(getTrackHeader(trpts[j].getTranscriptID()+"_val", null, VISIBILITY_PACK, null, false, true,
								null, -1, -1, null, null)+"\n");
						buffy.flush(); buffy.close();
						BEDwrapper writer= new BEDwrapper(geneBedName);
							writer.setBeds(beds);
						writer.write(true);
						
						buffy= new BufferedWriter(new FileWriter(geneBedName, true));
						buffy.write("\n\n");
						buffy.write(getTrackHeader(trpts[j].getTranscriptID()+"_col", null, VISIBILITY_PACK, null, true, false,
								null, -1, -1, null, null)+"\n");
						buffy.flush(); buffy.close();
						writer= new BEDwrapper(geneBedName);
							writer.setBeds(beds);
						writer.write(true);
					} catch (Exception e) {
						e.printStackTrace();
					}
					
					for (int k = 0; k < beds.length; k++) {
						beds[k].setName(trpts[j].getTranscriptID()+ "_"+ beds[k].getName());
						v.add(beds[k]);
					}
				}	// end of transcripts

				if (v.size()> 0) {
					BEDobject[] beds= (BEDobject[]) Arrays.toField(v);
					try {
						BufferedWriter buffy= new BufferedWriter(new FileWriter(geneBedName, true));
						buffy.write("\n\n");
						buffy.write(getTrackHeader("all_val", null, VISIBILITY_SQUISH, null, false, true,
								null, -1, -1, null, null)+"\n");
						buffy.flush(); buffy.close();
						BEDwrapper writer= new BEDwrapper(geneBedName);
							writer.setBeds(beds);
						writer.write(true);
	
						buffy= new BufferedWriter(new FileWriter(geneBedName, true));
						buffy.write("\n\n");
						buffy.write(getTrackHeader("all_col", null, VISIBILITY_SQUISH, null, true, false,
								null, -1, -1, null, null)+"\n");
						buffy.flush(); buffy.close();
						writer= new BEDwrapper(geneBedName);
							writer.setBeds(beds);
						writer.write(true);
					} catch (Exception e) {
						e.printStackTrace();
					}
				}
			}	// genes
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
			System.gc();
			Thread.yield();
		}

	}

	public static void writeBEDfile(File inFile, File outFile) {
		GTFChrReader reader= new GTFChrReader(inFile.getAbsolutePath());
		reader.setReadFeatures(new String[] {"exon", "CDS", "domain"});
		BEDwrapper writer= new BEDwrapper(outFile.getAbsolutePath());
	
		reader.setChromosomeWise(true);
		try {
			reader.read();
			BufferedWriter buffy= new BufferedWriter(new FileWriter(outFile.getAbsolutePath()));
			buffy.write(getTrackHeader("Domains", null, VISIBILITY_PACK, null, false, true,
						null, -1, -1, null, null));
			buffy.flush(); buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		Gene[] genes= reader.getGenes();
		
		while (genes!= null) {
			for (int i = 0; i < genes.length; i++) {
				Transcript[] trpts= genes[i].getTranscripts();
				for (int j = 0; j < trpts.length; j++) {
					if (trpts[j].getAttributes()== null|| trpts[j].getAttributes().size()== 0)
						continue;
					Object[] keys= trpts[j].getAttributes().keySet().toArray();
					BEDobject[] beds= new BEDobject[keys.length];
					for (int k = 0; k < keys.length; k++) {
						DirectedRegion[] doms= (DirectedRegion[]) trpts[j].getAttribute(keys[k]);
						int minStart= Integer.MAX_VALUE, maxEnd= -1; 
						for (int m = 0; m < doms.length; m++) {	// each domain a line
							int len= doms[m].getLength();
							int absStart= Math.abs(doms[m].getStart());
							if (absStart< minStart)
								minStart= absStart;
							if (absStart+ len> maxEnd)
								maxEnd= absStart+ len;
						}
						BEDobject obj= new BEDobject(
								trpts[j].getChromosome(), 
								trpts[j].getStrand(),
								minStart, maxEnd);
						obj.setName(trpts[j].getTranscriptID()+":"+doms[0].getID());
						int score= 1000;
						if (doms[0].getScore()> 0d)
							score= (int) Math.abs(Math.log(doms[0].getScore())* 10);
						if (score> 1000)
							score= 1000;
						obj.setScore(score);
						int cdsStart= Math.abs(trpts[j].getTranslations()[0].getStart());
						int cdsEnd= Math.abs(trpts[j].getTranslations()[0].getEnd());
						int thickStart= cdsStart, thickEnd= cdsEnd;
						if (cdsStart< obj.getStart())
							thickStart= obj.getStart();
						if (cdsEnd> obj.getEnd())
							thickEnd= obj.getEnd();
						obj.setThickStart(thickStart);
						obj.setThickEnd(thickEnd);
	
						obj.setBlockCount(doms.length);
						for (int m = 0; m < doms.length; m++) 	// each domain a line
							obj.addBlockStart(Math.abs(doms[m].getStart())- obj.getStart(), doms[m].getLength());
						
						beds[k]= obj;
					}
					
					writer.setBeds(beds);
					writer.write(true);
				}
			}
			
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			genes= reader.getGenes();
			System.gc();
			Thread.yield();
		}
	
	}
}

