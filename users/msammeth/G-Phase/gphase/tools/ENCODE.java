package gphase.tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.Vector;

import com.sun.org.apache.xml.internal.utils.StringToIntTable;

import gphase.io.gtf.GTFObject;
import gphase.model.DirectedRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Region;
import gphase.model.SpliceSite;

public class ENCODE {
	static DefaultRegion[] encRegions= null;
	
	public static String mapFileName= "encode/encode_regions_coords_NCBI35.txt";
	HashMap hm= null;
	
	
	public static String[] getEncodeRegionNames() {
		String[] regNames= new String[getEncodeRegions().length];
		for (int i = 0; i < regNames.length; i++) 
			regNames[i]= getEncodeRegions()[i].getID();
		return regNames;
	}
	
	public ENCODE() {
		DefaultRegion[] reg= getEncodeRegions();
		hm= new HashMap();
		for (int i = 0; i < reg.length; i++) {
			Vector v= (Vector) hm.get(reg[i].getChromosome());
			if (v== null)
				v= new Vector();
			v.add(reg[i]);
			hm.put(reg[i].getChromosome(), v);
		}
	}
	
	public boolean isInEncode(Region reg) {
		Vector v= (Vector) hm.get(reg.getChromosome());
		for (int i = 0; v!= null&& i < v.size(); i++) {
			Region r= (Region) v.elementAt(i);
			if (r.contains(reg))
				return true;
		}
		return false;
	}
	
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

	public static DirectedRegion convertToChromosomalCoord(DirectedRegion region) {
		DefaultRegion[] reg= getEncodeRegions();
			// normalize coordinates
		Region regionN= region.getAbsoluteRegion();
		for (int i = 0; i < reg.length; i++) 
			if (reg[i].getID().equals(region.getChromosome())) {
				region.setStart(regionN.getStart()+ reg[i].getStart());
				region.setEnd(regionN.getEnd()+ reg[i].getStart());
				region.setID(reg[i].getChromosome());
				return region;
			}
		return null;
	}

	public static GTFObject convertToEncodeCoord(GTFObject o) {
		
		DefaultRegion[] reg= getEncodeRegions();
		
			// normalize coordinates
		DefaultRegion regionN= new DefaultRegion(o.getStart(), o.getEnd());
		regionN.setChromosome(o.getSeqname());
		for (int i = 0; i < reg.length; i++) 
			if (reg[i].contains(regionN)) {
				o.setStart(regionN.getStart()- reg[i].getStart());
				o.setEnd(regionN.getEnd()- reg[i].getStart());
				o.setSeqname(reg[i].getID());
				return o;
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
