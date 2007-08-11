package gphase.db;

import gphase.Constants;
import gphase.io.TabDelimitedFormatWrapper;
import gphase.io.gtf.GTFObject;
import gphase.io.gtf.GTFWrapper;
import gphase.model.Species;
import gphase.model.Translation;
import gphase.tools.Arrays;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.Collator;
import java.util.HashMap;
import java.util.StringTokenizer;
import java.util.Vector;

public class MapTable {

	public static String[] notFoundIDs= null;
	static String[][] PRIMARY_GENE_ID= new String[][] {
		new String[]{"knownGene.name", "refGene.name", "ensGene.name", "mgcGenes.name"}, // human
		new String[]{""}, // chimp
		new String[]{"knownGene.name", "refGene.name", "ensGene.name", "mgcGenes.name"},	// mouse
		new String[]{""}, // rat
		new String[]{""}, // dog
		new String[]{""}, // cow
		new String[]{""}, // opossum
		
		new String[]{""}, // chicken
		new String[]{""}, // frog
		new String[]{""}, // zebrafish
		new String[]{""}, // fugu
		new String[]{""}, // tetraodon
		new String[]{"refGene.name", "Xref.name"}, // FlyBase
		new String[]{""}, // mosquito
		new String[]{""}, // honeybee
		new String[]{"refGene.name", "sangerGene.name"}, // worm
		new String[]{""}, // yeast
		new String[]{""}, // seasquirt
		new String[]{""}, // seasquirt2
		new String[]{""}, // platypus
		new String[]{""} // cress
	};
	
	public static void main(String[] args) {
		//addXRefsUCSC();
		addXRefsHPRD("human", "RefSeq");
	}
	
	public static void addXRefsUCSC() {
//		addXRefsToGTF("worm", "RefSeq");
//		addXRefsToGTF("worm", "WormBase");
//		addXRefsToGTF("fruitfly", "RefSeq");
//		addXRefsToGTF("fruitfly", "FlyBase");
//		addXRefsToGTF("mouse", "RefSeq");
//		addXRefsToGTF("mouse", "EnsEmbl");
		addXRefsToGTF("human", "RefSeq");
//		addXRefsToGTF("human", "EnsEmbl");
//		addXRefsToGTF("human", "UCSCGenes");
	}
	
	public static void addXRefsToGTF(String commonSpeName, String annoName) {
			String XREF_SFX= "_xref";
			
			// get maptable
			String[] list= new File(Constants.MAPTABLES_SUBDIR).list();
			String fName= null;
			for (int i = 0; i < list.length; i++) {
				if (list[i].startsWith(commonSpeName)&& list[i].contains(annoName)&&
						list[i].endsWith(".map")) {
					if (fName!= null)
						System.err.println("\tMore than one maptable found for ("+commonSpeName
								+ annoName+"): "+ fName+","+list[i]);
					fName= list[i];
				}
			}
			if (fName== null) {
				System.err.println("\tNo maptable found for "+commonSpeName+ " annotated in "+ annoName);
				return;
			}
			System.out.println("Reading maptable "+fName);
			TabDelimitedFormatWrapper tabReader= new TabDelimitedFormatWrapper(
					Constants.MAPTABLES_SUBDIR+ File.separator+ fName);
			try {
				tabReader.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			String[][] table= tabReader.getTable();
			String[] invalid= new String[] {"NULL", "N/A", "0", "MISS"};
			HashMap mapTrpt= new HashMap(table.length);
			for (int i = 1; i < table.length; i++) {		// 0 line contains defs
				HashMap mapXRefs= new HashMap(table[0].length);
				mapTrpt.put(table[i][0].trim(), mapXRefs);
				for (int j = 1; j < table[i].length; j++) {		// 0 col contains ID
					String[] tokens= null;
					if (table[i][j].indexOf(',')>= 0)
						tokens= table[i][j].split(",");
					else
						tokens= new String[] {table[i][j]};
					
					int cnt= 0;
					for (int k = 0; k < tokens.length; k++) {
						
							// throw out descriptions
						if (tokens[k].indexOf(' ')>= 0)
						//if (tokens[k].length()> 25)
							continue;
						
							// in refseq files only refseq IDs, for sylvain
						if (annoName.equalsIgnoreCase("RefSeq")&&
								(Translation.getProteinID(tokens[k])!= Translation.REFSEQ_ID))
							continue;
						
						
							// throw out numbers
						try {
							Integer.parseInt(tokens[k]);
							continue;
						} catch (NumberFormatException e) {
							; // :)
						}
						
						String tst= tokens[k].trim().toUpperCase();
						int c;	// check invalid
						for (c = 0; c < invalid.length; c++) 
							if (invalid[c].equals(tst))
								break;
						if (c< invalid.length)
							continue;
							
							// check against already added
						Object[] o= mapXRefs.values().toArray();
						for (c = 0; c < o.length; c++) {
							String former= ((String) o[c]).toUpperCase();
							if (former.contains(tst))
								break;
						}
						if (c< o.length)
							continue;
						
							// add
						int x= j;
						if (x>= table[0].length)
							x= table[0].length- 1;
						if (cnt> 0)
							mapXRefs.put(table[0][x]+cnt++, tokens[k]);
						else
							mapXRefs.put(table[0][x], tokens[k]);
					}
				}
			}
			
			// get gtf
			list= new File(Constants.SUBDIR_ANNOTATION_UCSC).list();
			fName= null;
			for (int i = 0; i < list.length; i++) {
				if (list[i].startsWith(commonSpeName)&& list[i].contains(annoName)
						&& list[i].endsWith(".gtf")) {
					if (fName!= null)
						System.err.println("\tMore than one annotation found for ("+commonSpeName
								+ annoName+"): "+ fName+","+list[i]);
					fName= list[i];
				}
			}
			if (fName== null) {
				System.err.println("\tNo annotation found for "+commonSpeName+ " annotated in "+ annoName);
				return;
			}
			
			System.out.println("Reading gtf "+fName);
			GTFWrapper gtfWrapper= new GTFWrapper(
					Constants.SUBDIR_ANNOTATION_UCSC+ File.separator+ fName);
			try {
				gtfWrapper.read();
			} catch (Exception e) {
				e.printStackTrace();
			}
			GTFObject[] gtfObjs= gtfWrapper.getGtfObj();
			for (int i = 0; i < gtfObjs.length; i++) {
				HashMap mapXRefs= (HashMap) mapTrpt.get(gtfObjs[i].getTranscriptID());
				if (mapXRefs== null)
					continue;
				Object[] vals= mapXRefs.values().toArray();
				String key= gtfObjs[i].getTranscriptID();
				String valStr= "";
				for (int j = 0; j < vals.length; j++) {
					if (key.equals(vals[j]))
						continue;
					valStr+= vals[j]+" ";
				}
				if (valStr.length()== 0)
					continue;
				gtfObjs[i].addAttribute("other_id", valStr.substring(0, valStr.length()- 1));
				
	//			Object[] keys= mapXRefs.keySet().toArray();
	//			for (int j = 0; j < keys.length; j++) {
	//				String key= ((String) keys[j]).replace('.', '_');
	//				gtfObjs[i].addAttribute(key, (String) mapXRefs.get(keys[j]));
	//			}
			}
			
				// write gtf
			gtfWrapper.setGtfObj(gtfObjs);
			gtfWrapper.setFileName(gtfWrapper.getFileName()+ XREF_SFX);
			System.out.println("Writing maptable "+gtfWrapper.getFileName());
			gtfWrapper.setSortAttributes(new String[] {GTFObject.GENE_ID_TAG, GTFObject.TRANSCRIPT_ID_TAG, GTFObject.EXON_ID_TAG});
			try {
				gtfWrapper.write();
			} catch (Exception e) {
				e.printStackTrace();
			}		
		}

	static void addXRefsHPRD(String commonSpeName, String annoName) {
		String SFX= "_HPRD";
		
		// get maptable
		String[] list= new File(Constants.MAPTABLES_SUBDIR).list();
		String fName= null;
		for (int i = 0; i < list.length; i++) {
			if (list[i].startsWith(commonSpeName)&& list[i].contains(annoName)&& 
					list[i].contains("HPRD")&& list[i].endsWith(".map")) {
				if (fName!= null)
					System.err.println("\tMore than one maptable found for ("+commonSpeName
							+ annoName+"): "+ fName+","+list[i]);
				fName= list[i];
			}
		}
		if (fName== null) {
			System.err.println("\tNo maptable found for "+commonSpeName+ " annotated in "+ annoName);
			return;
		}
		System.out.println("Reading maptable "+fName);
		TabDelimitedFormatWrapper tabReader= new TabDelimitedFormatWrapper(
				Constants.MAPTABLES_SUBDIR+ File.separator+ fName);
		try {
			tabReader.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		String[][] table= tabReader.getTable();
		String[] invalid= new String[] {"NULL", "N/A", "0", "MISS"};
		HashMap mapTrpt= new HashMap(table.length);
		for (int i = 1; i < table.length; i++) {		// 0 line contains defs
			HashMap mapXRefs= new HashMap(table[0].length);
			mapTrpt.put(table[i][0].trim(), mapXRefs);
			for (int j = 1; j < table[i].length; j++) {		// 0 col contains ID
				String[] tokens= null;
				if (table[i][j].indexOf(',')>= 0)
					tokens= table[i][j].split(",");
				else
					tokens= new String[] {table[i][j]};
				
				int cnt= 0;
				for (int k = 0; k < tokens.length; k++) {

					String tst= tokens[k].trim().toUpperCase();
					int c;	// check invalid
					for (c = 0; c < invalid.length; c++) 
						if (invalid[c].equals(tst))
							break;
					if (c< invalid.length)
						continue;
						
						// check against already added
					Object[] o= mapXRefs.values().toArray();
					for (c = 0; c < o.length; c++) {
						String former= ((String) o[c]).toUpperCase();
						if (former.contains(tst))
							break;
					}
					if (c< o.length)
						continue;
					
						// add
					int x= j;
					if (x>= table[0].length)
						x= table[0].length- 1;
					if (cnt> 0)
						mapXRefs.put(table[0][x]+cnt++, tokens[k]);
					else
						mapXRefs.put(table[0][x], tokens[k]);
				}
			}
		}
		
		// get gtf
		list= new File(Constants.SUBDIR_ANNOTATION_UCSC).list();
		fName= null;
		for (int i = 0; i < list.length; i++) {
			if (list[i].startsWith(commonSpeName)&& list[i].contains(annoName)
					&& list[i].endsWith(".gtf_xref")) {
				if (fName!= null)
					System.err.println("\tMore than one annotation found for ("+commonSpeName
							+ annoName+"): "+ fName+","+list[i]);
				fName= list[i];
			}
		}
		if (fName== null) {
			System.err.println("\tNo annotation found for "+commonSpeName+ " annotated in "+ annoName);
			return;
		}
		
		System.out.println("Reading gtf "+fName);
		GTFWrapper gtfWrapper= new GTFWrapper(
				Constants.SUBDIR_ANNOTATION_UCSC+ File.separator+ fName);
		try {
			gtfWrapper.read();
		} catch (Exception e) {
			e.printStackTrace();
		}
		GTFObject[] gtfObjs= gtfWrapper.getGtfObj();
		for (int i = 0; i < gtfObjs.length; i++) {
			HashMap mapXRefs= (HashMap) mapTrpt.get(gtfObjs[i].getTranscriptID());
			if (mapXRefs== null)
				continue;
			Object[] vals= mapXRefs.values().toArray();
			String key= gtfObjs[i].getTranscriptID();
			String valStr= gtfObjs[i].getAttribute("other_id");
			if (valStr== null)
				valStr= "";
			else
				valStr+= " ";
			for (int j = 0; j < vals.length; j++) {
				if (key.equals(vals[j]))
					continue;
				valStr+= "HPRD:"+vals[j]+" ";
			}
			if (valStr.length()== 0)
				continue;
			gtfObjs[i].addAttribute("other_id", valStr.substring(0, valStr.length()- 1));
			
//			Object[] keys= mapXRefs.keySet().toArray();
//			for (int j = 0; j < keys.length; j++) {
//				String key= ((String) keys[j]).replace('.', '_');
//				gtfObjs[i].addAttribute(key, (String) mapXRefs.get(keys[j]));
//			}
		}
		
			// write gtf
		gtfWrapper.setGtfObj(gtfObjs);
		gtfWrapper.setFileName(gtfWrapper.getFileName()+ SFX);
		System.out.println("Writing maptable "+gtfWrapper.getFileName());
		gtfWrapper.setSortAttributes(new String[] {GTFObject.GENE_ID_TAG, GTFObject.TRANSCRIPT_ID_TAG, GTFObject.EXON_ID_TAG});
		try {
			gtfWrapper.write();
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	
	public static void addPDBInfo() {
	
		String commonSpeName= "human";
		String[] ids= null, pdb= null;
			// read in additional maps
		try {
			TabDelimitedFormatWrapper reader= new TabDelimitedFormatWrapper(
					new File("maptables"+File.separator+"elza.maptable").getAbsolutePath());
			reader.read();
			pdb= reader.getColumn(0);
			ids= reader.getColumn(1);
			for (int i = 0; i < ids.length; i++) 
				ids[i]= ids[i].substring(ids[i].indexOf("via")+ 4, ids[i].length());
			System.out.println("read "+ids.length+" additional ids.");
		} catch (Exception e) {
			e.printStackTrace();
		}
		
			// add to maptables
		String[] specAnnos= Species.SPECIFIC_ANNOTATIONS[Species.getSpeciesNumber(commonSpeName)];
		for (int i = 0; i < specAnnos.length; i++) {
			System.out.print("\tadding to "+specAnnos[i]+" ");
			byte[] found= new byte[ids.length];
			for (int j = 0; j < found.length; j++) 
				found[j]= 0;
			
			String[] list= new File(Constants.MAPTABLES_SUBDIR).list();
			String fName= null;
			for (int j = 0; j < list.length; j++) {
				if (list[j].startsWith(commonSpeName)&& list[j].contains(specAnnos[i])) {
					if (fName!= null)
						System.err.println("More than one maptable found for ("+commonSpeName
								+ specAnnos[i]+"): "+ fName+","+list[j]);
					fName= list[j];
				}
			}
			if (fName== null) {
				System.err.println("No maptable found for "+commonSpeName+", "+specAnnos[i]);
				continue;
			}
			
			try {
				String inName= Constants.MAPTABLES_SUBDIR+ File.separator+fName;
				File inFile= new File(inName);
				File modInFile= new File(inName+"_original");
				inFile.renameTo(modInFile);
				File outFile= new File(inName);
				BufferedReader buffy= new BufferedReader(new FileReader(modInFile));
				BufferedWriter writer= new BufferedWriter(new FileWriter(outFile));
				
				long bytesRead= 0l;
				long size= modInFile.length();
				int lastPerc= 0;
				while (buffy.ready()) {
					String line= buffy.readLine();
					bytesRead+= line.length()+ 1;
					int perc= (int) ((bytesRead* 10d)/ size);
					if (perc> lastPerc) {
						++lastPerc;
						System.out.print("*");
						System.out.flush();
					}
					
					String[] tokens= line.split("\t");
					for (int j = 0; j < ids.length; j++) 
						if (found[j]== 0) { // assuming pdb:otherID = 1:1
							int k;
							for (k = 0; k < tokens.length; k++) {
								if (tokens[k].contains(ids[j])) {
									line+= "\t"+ pdb[j];
									++found[j];
								}
								break;	// assuming only one pdb id maps per line..
							}
							if (k< tokens.length)
								break;	
						}
					writer.write(line+"\n");
				}
				writer.flush(); writer.close();
				buffy.close();
				
			} catch (Exception e) {
				System.out.println();
				e.printStackTrace();
			}
			System.out.println();
			int cntFound= 0;
			for (int j = 0; j < found.length; j++) 
				if (found[j]> 0)
					++cntFound;
			System.out.println("\tsuccessfully mapped "+cntFound+" of "+found.length+" ids.");
		}
		
	}
	
	public static String[] getPrimaryGeneID(String commonSpeName, String annoName, String[] someID) {
		
		String priGeneID= PRIMARY_GENE_ID[Species.getSpeciesNumber(commonSpeName)][Species.getAnnotationNumber(commonSpeName, annoName)];		
		Vector v= null;
		int colNr= 0;
		int[] hitCount= new int[someID.length];
		for (int i = 0; i < hitCount.length; i++) 
			hitCount[i]= 0;
		int speNr= Species.getSpeciesNumber(commonSpeName);
		System.out.print("\tconverting to gene ids "+annoName+" ");
		System.out.flush();
		try {
			String[] list= new File(Constants.MAPTABLES_SUBDIR).list();
			String fName= null;
			for (int i = 0; i < list.length; i++) {
				if (list[i].startsWith(commonSpeName)&& list[i].contains(annoName)) {
					if (fName!= null)
						System.err.println("\tMore than one maptable found for ("+commonSpeName
								+ annoName+"): "+ fName+","+list[i]);
					fName= list[i];
				}
			}
			if (fName== null) {
				System.err.println("\tNo maptable found for "+commonSpeName);
				return null;
			}
			File f= new File(Constants.MAPTABLES_SUBDIR+ File.separator+fName);
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			String line= buffy.readLine();
			String[] tokens= line.split("\t");
			for (int i = 0; i < tokens.length; i++) {
				if (tokens[i].contains(priGeneID))
					break;
				++colNr;
			}
			if ((colNr+1)> tokens.length)
				System.err.println("\tPrimary key column "+priGeneID+" not found.");
			v= new Vector();
			int lineCtr= 0;
			long bytesRead= 0l;
			int lastPerc= 0;
			long fSize= f.length();
			while (buffy.ready()) {
				line= buffy.readLine();
				if (line== null)
					break;
				bytesRead= line.length()+ 1;
				++lineCtr;
				if ((bytesRead* 10l/ fSize)> lastPerc) {	// (lineCtr%1000== 0) {
					System.out.print("*");
					System.out.flush();
					++lastPerc;
				}
				
					// compare for geneID
				for (int i = 0; i < someID.length; i++) {	// check only not found ones? hitCount[i]< 1 ???
					StringTokenizer toki= new StringTokenizer(line, "\t,");
					while (toki.hasMoreTokens()) {
						String token= toki.nextToken();
						if (token.equals(someID[i])) {
							v.add(line);	
							++hitCount[i];
							break;	// find an ID only once per line
						}
					}
				}
			}
			buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		System.out.println();
		
		String s= "";
		String sDbl= "";
		int cnt= 0;
		int dblCnt= 0;
		Vector notFoundV= new Vector();
		for (int i = 0; i < hitCount.length; i++) {
			if (hitCount[i]== 0) {
				s+= someID[i]+", ";
				++cnt;
				notFoundV.add(someID[i]);
			}
			else if (hitCount[i]> 1) {
				sDbl+= someID[i]+" ["+hitCount[i]+"], ";
				++dblCnt;
			}
				
		}
		notFoundIDs= (String[]) Arrays.toField(notFoundV);
		if (s.length()> 0)
			System.out.println("\tNo map found for "+cnt+" IDs");	// +s
		if (sDbl.length()> 0)
			System.out.println("\tMultiple map found for "+dblCnt+" IDs: "+sDbl.substring(0, sDbl.length()- 2));
		
		Vector geneIDvec= new Vector(v.size());
		for (int i = 0; i < v.size(); i++) {
			String[] cols= ((String) v.elementAt(i)).split("\t");
			String geneID= null;
			try {
				geneID= cols[colNr];
				if (commonSpeName.equals("fruitfly")&& annoName.equals("FlyBase"))
					geneID= geneID.substring(0, geneID.indexOf("-"));	// CG10033-RB
				Arrays.addUnique(geneIDvec, geneID, Collator.getInstance());
			} catch (ArrayIndexOutOfBoundsException e) {
				System.err.println("\tNot found ID "+priGeneID+" (col "+colNr+") in:\n"+v.elementAt(i));
			}
		}
		
		return (String[]) Arrays.toField(geneIDvec);
	}

	public static String[] getNotFoundIDs() {
		return notFoundIDs;
	}
}
