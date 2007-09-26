package gphase.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;
import java.util.Vector;

import com.sun.org.apache.xml.internal.utils.StringToIntTable;

import gphase.model.DirectedRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Region;
import gphase.model.SpliceSite;

public class GPhase {
	static DefaultRegion[] encRegions= null;
	
	public static String mapFileName= "encode/encode_regions_coords_NCBI35.txt";
	public static DirectedRegion convertToEncodeCoord(DirectedRegion region) {
		DefaultRegion[] reg= getEncodeRegions();
			// normalize coordinates
		DefaultRegion regionN= new DefaultRegion(Math.abs(region.getStart()), Math.abs(region.getEnd()));
		regionN.setChromosome(region.getChromosome());
		for (int i = 0; i < reg.length; i++) 
			if (reg[i].contains(regionN)) {
				region.setStart(regionN.getStart()- reg[i].getStart());
				region.setEnd(regionN.getEnd()- reg[i].getStart());
				region.setID(reg[i].getID());
				return region;
			}
		return null;
	}

	public static DefaultRegion getEncodeRegion(int pos) {
		DefaultRegion[] reg= getEncodeRegions();
		for (int i = 0; i < reg.length; i++) {
			if (reg[i].getStart()<= pos&& reg[i].getEnd()>= pos)
				return reg[i];
		}
		return null;
	}
	
	public static DefaultRegion convertToEncodeCoord(SpliceSite ss) {
		DefaultRegion[] reg= getEncodeRegions();
			// normalize coordinates
		int eins= Math.abs(ss.getPos());
		int start= Math.abs(ss.getPos());
		if (ss.isDonor()) 
			start= Math.abs(ss.getPos())-1;
		int end= start+1;
		DefaultRegion regionN= new DefaultRegion(start, end);
		regionN.setChromosome(ss.getGene().getChromosome());
		for (int i = 0; i < reg.length; i++) 
			if (reg[i].contains(regionN)) {
				regionN.setStart(regionN.getStart()- reg[i].getStart());
				regionN.setEnd(regionN.getEnd()- reg[i].getStart());
				regionN.setID(reg[i].getID());
				return regionN;
			}
		return null;
	}
	
	public static DefaultRegion[] getEncodeRegions() {
		if (encRegions== null) {
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(new File(mapFileName).getAbsolutePath()));
				String line= null;
				Vector regV= new Vector(44);
				while (buffy.ready()) {
					line= buffy.readLine().trim();
					if (line.length()< 1|| line.startsWith("#"))
						continue;
					// #name	chrom	chromStart	chromEnd
					StringTokenizer toki= new StringTokenizer(line);
					if (toki.countTokens()!= 4)
						System.err.println("Error parsing line "+line);
					String encName= toki.nextToken();
					String chrom= toki.nextToken();
					chrom= chrom.substring(chrom.indexOf("chr")+ 3);
					int start= Integer.parseInt(toki.nextToken());
					int end= Integer.parseInt(toki.nextToken());
					DefaultRegion reg= new DefaultRegion(start, end);
					reg.setChromosome(chrom);
					reg.setID(encName);
					regV.add(reg);
				}
				encRegions= (DefaultRegion[]) Arrays.toField(regV);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return encRegions;
	}
}
