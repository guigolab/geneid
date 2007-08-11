package gphase.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.StringTokenizer;

import gphase.Constants;
import gphase.model.Gene;
import gphase.model.GeneHomology;
import gphase.model.Species;

public class OrthologousGeneWrapper {
	
	String fName= null;
	int lineCtr= 0;
	Species spe1= null, spe2= null;
	boolean rev= false;
	
	public OrthologousGeneWrapper(Species spe1, Species spe2) {

		this.spe1= spe1;
		this.spe2= spe2;
		String fName1= spe1.getAbbreviatedName()+ "-"+ spe2.getAbbreviatedName()+".UBRH";
		String fName2= spe2.getAbbreviatedName()+ "-"+ spe1.getAbbreviatedName()+".UBRH";
		String fPath= Constants.DATA_DIR+ File.separator+ Constants.SEQUENCES_SUBDIR+  File.separator+ 
			"caipirinha"+ File.separator+ "UBRH";
		String[] list= new File(fPath).list();
		for (int i = 0; i < list.length; i++) {
			if (list[i].equalsIgnoreCase(fName1)) {
				this.fName= fName1;
				break;
			}
			if (list[i].equalsIgnoreCase(fName2)) {
				this.fName= fName2;
				rev= true;
				break;
			}
		}
			
		assert(fName!= null);
		this.fName= fPath+File.separator+fName;	// make absolute path
	}
	
	public GeneHomology getNextGeneHomology() {
		String line= null;
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(fName));
			for (int i = 0; i < lineCtr; i++) 
				buffy.readLine();
			line= buffy.readLine();
			++lineCtr;
			buffy.close();
			
			StringTokenizer toki= new StringTokenizer(line);
			assert(toki.countTokens()== 2);
			Gene ge1= new Gene(spe1, toki.nextToken());
			Gene ge2= new Gene(spe2, toki.nextToken());
			if (rev) {
				Gene h= ge1;
				ge1= ge2;
				ge2= h;
			}
			GeneHomology homology= new GeneHomology(ge1, ge2);
			return homology;
		
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}

	}
}
