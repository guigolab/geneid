package gphase.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import gphase.Toolbox;
import gphase.io.gtf.GTFChrReader;
import gphase.io.gtf.GTFObject;
import gphase.model.DirectedRegion;
import gphase.tools.File;

public class OverlapComparator {
	
	File file1= null, file2= null;
	
	
	public static void main(String[] args) {
		OverlapComparator myComparator= new OverlapComparator();
		myComparator.setFile1(new File(args[0]));
		myComparator.setFile2(new File(args[1]));
		myComparator.run();
	}
	
	public OverlapComparator() {
		
	}
	
	public void run() {
		
		GTFChrReader reader= new GTFChrReader(file1.getAbsolutePath());
		reader.setReadGene(false);
		reader.setReadGTF(true);
		reader.setReadFeatures(new String[] {"mpss_tag_1", "mpss_tag_2", "sage_tag_1", "sage_tag_2", "sage_tag_3", "sage_tag_4", "mpss_tag_3", "mpss_tag_4"});
		try {
			reader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		GTFObject[] obs= reader.getGtfObj();
		Comparator<GTFObject> gtfStartEndStrandCompi= GTFObject.getDefaultPositionComparator();
		long astaBytesRead= 0l, maxAstaBytesRead= 0l;
		HashMap<String, Long> astaChrRead= new HashMap<String, Long>();
		String fName= Toolbox.checkFileExists(file1.getPathOnly()+File.separator+file1.getFileNameOnly()+"_overlapWith_"+file2.getFileNameOnly()+"."+file1.getExtension());
		File outFile1= new File(fName); 
		fName= Toolbox.checkFileExists(file2.getPathOnly()+File.separator+file2.getFileNameOnly()+"_overlapWith_"+file1.getFileNameOnly()+"."+file2.getExtension());
		File outFile2= new File(fName); 
		BufferedReader buffy= null;
		InputStream iStream= null;
		
		while (obs!= null) {
			
			Arrays.sort(obs, gtfStartEndStrandCompi);
			String chr1= obs[0].getChromosome();
			try {
				iStream= new FileInputStream(file2);
			} catch (FileNotFoundException e3) {
				e3.printStackTrace();
			}
			
			if (astaChrRead.get(chr1)!= null) 
				astaBytesRead= astaChrRead.get(chr1).longValue();
			else
				astaBytesRead= maxAstaBytesRead;
			
			try {
				iStream.skip(astaBytesRead);
			} catch (IOException e2) {
				e2.printStackTrace();
			}
			buffy= new BufferedReader(new InputStreamReader(iStream));
			
				// skip
			String line= null;
			try {
				line= buffy.readLine();
			} catch (IOException e1) {
				e1.printStackTrace();
			}
			astaBytesRead+= line.length()+ 1;
			GTFObject[] o= ASTA3reader.getGTFObject_exclusiveExonicAreas(line);
			String lastChr= null;
			while (!chr1.equals(o[0].getChromosome())) {
				if (!lastChr.equals(o[0].getChromosome())) {
					if (astaChrRead.get(o[0].getChromosome())== null) 
						astaChrRead.put(o[0].getChromosome(), new Long(astaBytesRead- line.length()- 1));
					lastChr= o[0].getChromosome();
				}
				
				try {
					line= buffy.readLine();
				} catch (IOException e) {
					e.printStackTrace();
				}
				astaBytesRead+= line.length()+ 1;
				o= ASTA3reader.getGTFObject_exclusiveExonicAreas(line);
			}
			
				// read
			HashMap<GTFObject, GTFObject> v1= new HashMap<GTFObject,GTFObject>(obs.length/ 2),
				v2= new HashMap<GTFObject,GTFObject>(obs.length/ 2);
			while (chr1.equals(o[0].getChromosome())&& line!= null) {
				GTFObject baseO= ASTA3reader.getGTFObject(line);
				HashMap<GTFObject, GTFObject> map= new HashMap<GTFObject, GTFObject>();
				for (int x = 0; x < o.length; x++) {
					int idx= Arrays.binarySearch(obs, o[x], gtfStartEndStrandCompi);
					if (idx< 0)
						idx= -(idx+1);
					if (idx== obs.length)
						--idx;
					for (int i = idx; i >= 0; --i) {
						if (map.get(obs[i])!= null)
							continue;
						if (obs[i].getEnd()< o[x].getStart()|| obs[i].getStrand()!= o[x].getStrand())
							break;
						boolean overlaps= checkOverlap(obs[i], o[x]);
						if (overlaps) {
							if (v1.get(obs[i])== null)
								v1.put(obs[i], obs[i]);
							if (v2.get(baseO)== null)
								v2.put(baseO, baseO);
							map.put(obs[i], obs[i]);
						}
					}
					for (int i = idx; i < obs.length; ++i) {
						if (map.get(obs[i])!= null)
							continue;
						if (obs[i].getStart()> o[x].getEnd()|| obs[i].getStrand()!= o[x].getStrand())
							break;
						boolean overlaps= checkOverlap(obs[i], o[x]);
						if (overlaps) {
							if (v1.get(obs[i])== null)
								v1.put(obs[i], obs[i]);
							if (v2.get(baseO)== null)
								v2.put(baseO, baseO);
							map.put(obs[i], obs[i]);
						}
					}
				}

				
				try {
					line= buffy.readLine();
				} catch (IOException e) {
					e.printStackTrace();
				}
				if (line!= null)  {
					astaBytesRead+= line.length()+ 1;
					o= ASTA3reader.getGTFObject_exclusiveExonicAreas(line);
				}
			}
			if (line!= null)
				astaBytesRead-= line.length()+ 1;
			maxAstaBytesRead= Math.max(astaBytesRead, maxAstaBytesRead);
			
				// output
			try {
				GTFChrReader writer= new GTFChrReader(outFile1.getAbsolutePath());
				GTFObject[] outObs= new GTFObject[v1.size()];
				Object[] oo= v1.values().toArray();
				for (int i = 0; i < oo.length; i++) 
					outObs[i]= (GTFObject) oo[i];
				writer.setGtfObj(outObs);
				writer.write(true);
				
				writer= new GTFChrReader(outFile2.getAbsolutePath());
				outObs= new GTFObject[v2.size()];
				oo= v2.values().toArray();
				for (int i = 0; i < oo.length; i++) 
					outObs[i]= (GTFObject) oo[i];
				writer.setGtfObj(outObs);
				writer.write(true);
			} catch (Exception e) {
				e.printStackTrace();
			}
				
				// read next
			try {
				reader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			obs= reader.getGtfObj();
		}
	}

	private boolean checkOverlap(GTFObject object, GTFObject o) {
		final String GTF_ASEVENT_ID= "asEvent";
		final String GTF_TAG_ID= "tags";
		if (o.overlaps(object)) {
			String s= object.removeAttribute(GTF_ASEVENT_ID);
			if (s== null)
				s= "1";
			else
				s= Integer.toString(Integer.parseInt(s)+ 1);
			//s+= o.getAttribute(GTFObject.ATTRIBUTE_EVENT_STRUCTURE)+ " "+o.getChromosome()+":"+o.getStart()+"-"+o.getEnd()+ ", ";
			object.addAttribute(GTF_ASEVENT_ID, s);
			
			s= o.removeAttribute(GTF_TAG_ID);
			if (s== null)
				s= "1";
			else
				s= Integer.toString(Integer.parseInt(s)+ 1);
			//s+= object.getFeature()+" "+object.getChromosome()+ ":"+ object.getStart()+ "-"+ object.getEnd()+", ";
			o.addAttribute(GTF_TAG_ID, s);
			
			return true;
		}
		
		return false;
	}

	public File getFile1() {
		return file1;
	}

	public void setFile1(File file1) {
		this.file1 = file1;
	}

	public File getFile2() {
		return file2;
	}

	public void setFile2(File file2) {
		this.file2 = file2;
	}
}
