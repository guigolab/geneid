package gphase.io;

import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gphase.io.gtf.GTFObject;

public class ASTA3reader {
	
	public static GTFObject getGTFObject(String astaString) {
		String[] tokens= astaString.split("\t");
		GTFObject o= new GTFObject();
		o.setFeature(GTFObject.FEATURE_ASEVENT);
		String[] h1= tokens[1].split(":");
		o.setSeqname(h1[0]);
		String[] h2= h1[1].split("-");
		o.setStart(Integer.parseInt(h2[0]));
		o.setEnd(Integer.parseInt(h2[1]));
		
			// get strand
		if (tokens[2].equals("-"))
			o.setStrand(-1);
		else if (tokens[2].equals("+"))
			o.setStrand(1);
		else
			System.err.println("strand uninterpretable "+tokens[2]);
		
		o.addAttribute(GTFObject.ATTRIBUTE_EVENT_STRUCTURE, tokens[0]);
		
		return o;
	}

	public static GTFObject getGTFObject_variableArea(String astaString) {
		String[] tokens= astaString.split("\t");
		GTFObject o= new GTFObject();
		o.setFeature(GTFObject.FEATURE_ASEVENT);
		String[] h1= tokens[1].split(":");
		o.setSeqname(h1[0]);
		String[] h2= h1[1].split("-");
		
			// splice Chains
		int min= Integer.MAX_VALUE, max= Integer.MIN_VALUE;
		for (int idx= 4; idx< tokens.length; idx+= 2) {
			if (tokens[idx].length()== 0) 
				continue;			
			String[] ss= tokens[idx].split(",");
			int h= Integer.parseInt(ss[0]);
			if (h< min)
				min= h;
			h= Integer.parseInt(ss[ss.length-1]);
			if (h> max)
				max= h;
			idx+= 2;
		}
		o.setStart(min);
		o.setEnd(max);
		
			// get strand
		if (tokens[2].equals("-"))
			o.setStrand(-1);
		else if (tokens[2].equals("+"))
			o.setStrand(1);
		else
			System.err.println("strand uninterpretable "+tokens[2]);
		
		o.addAttribute(GTFObject.ATTRIBUTE_EVENT_STRUCTURE, tokens[0]);
		
		return o;
	}

	// works probably only pw
	public static GTFObject[] getGTFObject_exclusiveExonicAreas(String astaString) {
		String[] tokens= astaString.split("\t");
		GTFObject o= new GTFObject();
		o.setFeature(GTFObject.FEATURE_ASEVENT);
		String[] h1= tokens[1].split(":");
		o.setSeqname(h1[0]);
		String[] h2= h1[1].split("-");

		// get strand
		if (tokens[2].equals("-"))
			o.setStrand(-1);
		else if (tokens[2].equals("+"))
			o.setStrand(1);
		else
			System.err.println("strand uninterpretable "+tokens[2]);

			// parse structure
		Vector<GTFObject> v= new Vector<GTFObject>(2);
		String[] structChains= tokens[0].split(",");
		String[][] struct= new String[structChains.length][];
		for (int i = 0; i < struct.length; i++) {
			StringTokenizer toki= new StringTokenizer(structChains[i], "*^-;",true);
			if (toki.countTokens()== 1)	// 0
				struct[i]= new String[0];
			else
				struct[i]= new String[toki.countTokens()];
			for (int j = 0; j < struct[i].length; j++) 
				struct[i][j]= toki.nextToken();
		}
		
		int[][] spliceChains= new int[structChains.length][];
		for (int i = 0; i < struct.length; i++) {
			String[] t= tokens[4+(i*2)].split(",");
			if (t.length> 0&& t[0].equals(""))
				spliceChains[i]= new int[0];
			else
				spliceChains[i]= new int[t.length];
			for (int j = 0; j < spliceChains[i].length; j++) 
				spliceChains[i][j]= Integer.parseInt(t[j]);			
		}
		
		
		
		int[] idx= new int[struct.length];
		for (int i = 0; i < idx.length; i++) 
			idx[i]= 0;
		int lastIdx= 0, currIdx= 1;
		String lastType= "";
		while (true) {
			for (int i = 0; i < struct.length; i++) {
				if (idx[i]== struct[i].length)
					continue;
				if (Integer.parseInt(struct[i][idx[i]])== currIdx) {
					if (currIdx> 1) {
						if ((((lastType.equals("*")|| lastType.equals("-"))&& (struct[i][idx[i]+1].equals("*")|| struct[i][idx[i]+1].equals("-")))||
								((lastType.equals("^")|| lastType.equals(";"))&& (struct[i][idx[i]+1].equals("^")|| struct[i][idx[i]+1].equals(";"))))|| // same type, must be from diff splice chains; 
								(lastIdx== i)) {	//else: if its from different splice chains we have ovl exons
							
							boolean skip= false;
							if (lastType.equals("^")|| lastType.equals(";")) {	// introns
								int j = 0;
								for (; j < idx.length; j++) {
									if (j== i)
										continue;
									if (struct[j].length== 0|| ((!struct[j][1].equals("*"))&& (!struct[j][struct[j].length-1].equals(";")))||
										(struct[j][1].equals("*")&& idx[j]> 0)|| (struct[j][struct[j].length-1].equals(";")&& idx[j]< struct[j].length)) {
										break;
									}
								}
								if (j== idx.length)
									skip= true;
							}

							if (!skip) {
								GTFObject obj= (GTFObject) o.clone();
								obj.setStart(spliceChains[lastIdx][(idx[lastIdx]/2) -1]);
								obj.setEnd(spliceChains[i][idx[i]/2]);
								v.add(obj);
							}
						}
					}
					lastIdx= i;
					lastType= struct[i][idx[i]+1];
					for (int j = 0; j < idx.length; j++) {	// for identical positions
						if (i== j)
							continue;
						if (idx[j]< struct[j].length&& Integer.parseInt(struct[j][idx[j]])== currIdx)
							idx[j]+= 2;
					}
					++currIdx;
					idx[i]+= 2;
				}
			}
			
			int c= 0;
			for (int i = 0; i < idx.length; i++) 
				if (idx[i]== struct[i].length)
					++c;
			if (c== struct.length)
				break;
		}
		
		GTFObject[] obs= new GTFObject[v.size()];
		for (int i = 0; i < obs.length; i++) 
			obs[i]= v.elementAt(i);
		
		return obs;
	}
}
