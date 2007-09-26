package gphase.io;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import gphase.io.gtf.GTFObject;

public class ASTA2reader {
	
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
			// look for empty schain
		boolean strandFound= false;
		for (int i = 0; i < tokens[3].length(); i++) {
			if (Character.isLetter(tokens[3].charAt(i))) {
				o.setStrand(-1);
				strandFound= true;
				break;
			}
		}
		if (!strandFound) {
			for (int i = 0; i < tokens[5].length(); i++) {
				if (Character.isLetter(tokens[3].charAt(i))) {
					o.setStrand(1);
					strandFound= true;
					break;
				}
			}
		}
		
			// look for first smaller splice site or shorter chain
		String[] tokens3= tokens[3].split(","); 
		String[] tokens5= tokens[5].split(",");
		for (int i = 0; i < tokens5.length; i++) {
			if (tokens3.length<= i) {
				o.setStrand(1);
				strandFound= true;
			}
			if (tokens5.length<= i) {
				o.setStrand(-1);
				strandFound= true;
			}
			if (strandFound)
				break;
			int v3= Integer.parseInt(tokens3[i]);
			int v5= Integer.parseInt(tokens5[i]);
			if (v3< v5) {
				o.setStrand(1);
				strandFound= true;
			} else if (v5> v3) {
				o.setStrand(-1);
				strandFound= true;
			}
			if (strandFound)
				break;
		}
		
		return o;
	}
}
