package repeat.algorithm;
import qalign.algo.CostTable;
import qalign.model.MultipleAlignmentModel;

/*
 * Created on Nov 26, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

/**
 * 
 * 
 * @author micha
 */
public class FragmentScoreAssigner {

	/**
	 * 
	 */
	public FragmentScoreAssigner() {
		super();
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) {
	
		System.out.println(score("ATTC", "AGGC", CostTable.BLOSUM62));	
	}
	
	public static int score(String str1, String str2, int type) {
		
		if (str1.length()!= str2.length())
			throw new IllegalArgumentException();
			
		int sum= 0;
		for (int i= 0; i < str1.length(); i++) 
			sum+= score(str1.charAt(i), str2.charAt(i), type);
		return sum;
	}


	public static byte score(char c1, char c2, int type) {
		
		return score(c1,c2,CostTable.getCostTable(type));
	}

	public static byte score(char c1, char c2, byte[][] table) {
		
		if (MultipleAlignmentModel.isGapChar(c1))
			c1= '-';
		if (MultipleAlignmentModel.isGapChar(c2))
			c2= '-';
		if (c1== c2)
			return 0;
		
		int p1,p2;			// lookup
		for (p1= 1; p1 < table.length; p1++) 
			if ((table[p1]!= null)&& (table[p1][0]== (char) c1))
				break;
					
		for (p2= 1; p2 < table.length; p2++) 
			if ((table[p2]!= null)&& (table[p2][0]== (char) c2))
				break;
				
		if (p2< p1) {		// swap if necessary (for just one triangle)
			int pos= p1;
			p1= p2;
			p2= pos;
			char cc= c1;
			c1= c2;
			c2= cc;
		}
		
		return (table[p2][p1]);
	}
}
