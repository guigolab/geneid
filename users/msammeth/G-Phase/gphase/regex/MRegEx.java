package gphase.regex;

import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * multi commands in brackets []
 * / ascending
 * \ falling
 * 3+ number of transcripts
 * e.g., (1-2^)[/3+]3-[\3+]
 * @author micha
 */
public class MRegEx {

	String regEx= null;
	
	public MRegEx(String newRegEx) {
		this.regEx= newRegEx;
	}
	
	String[] parse() {
		if (regEx== null)
			System.err.println("no reg expression");
		
		int last= -1;	
		Vector exprV= new Vector();
		Vector mQuantV= new Vector();
		
		for (int pos = 0; pos < regEx.length(); pos++) {	// tokenize multi expressions
			if (regEx.charAt(pos)== '[') { 
				for (int i = pos; i < regEx.length(); i++) {
					if (regEx.charAt(i)== ']') {
						exprV.add(regEx.substring(last+1, pos));
						mQuantV.add(regEx.substring(pos+1, i));
						pos= last= i;
						break;
					}
				}
			}
		}
		
		String[] result= null;
		Pattern patty= Pattern.compile("\\D*(\\d+)\\D?");
		for (int i = 0; i < exprV.size(); i++) {
			String quant= (String) mQuantV.elementAt(i);	// get quantifier
			Matcher matty= patty.matcher(quant);
			if (!matty.matches())
				System.err.println("invalid mExpr: "+ quant);
			int q= Integer.parseInt(matty.group(1));
			if (result== null) {
				result= new String[q];
				for (int x = 0; x < result.length; x++) 
					result[x]= "";
			}
		
			String mpatt= (String) exprV.elementAt(i);
			matty= patty.matcher(mpatt); 
				// assuming ascending /
			for (int j = 0; j < q; j++) {
				if (j== 0) {
					result[j]= mpatt;
					continue;
				}
				int pos= 0;
				while (matty.find(pos)) {
					int base= Integer.parseInt(matty.group(1));
					result[j]+= mpatt.substring(pos, matty.start(1))+
								(j+1)* base;
					pos= matty.end(1);
				}
				result[j]+= mpatt.substring(pos, mpatt.length());
			}
		}
		
		return result;
	}
	
	public static void main(String[] args) {
		MRegEx myMReg= new MRegEx("(1-+,2'+)[/3+]");
		myMReg.getAutomatons();
		System.currentTimeMillis();
	}
	
	public Automaton[] getAutomatons() {
		
		String[] regExs= parse();
		Automaton[] autos= new Automaton[regExs.length];
		for (int i = 0; i < regExs.length; i++) {
			RegExp regEx= new RegExp(regExs[i]);
			autos[i]= regEx.getAutomaton();
		}
		
		return autos;
	}
}
