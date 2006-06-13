/*
 * Created on Feb 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import com.sun.org.apache.xalan.internal.xsltc.runtime.Hashtable;
import com.sun.org.apache.xerces.internal.impl.xs.SubstitutionGroupHandler;

/**
 * 
 * 
 * @author msammeth
 */
public class ASVariation implements Serializable {

	static final long serialVersionUID = 2433674838499118768L;
	
	Transcript trans1= null;
	Transcript trans2= null;
	SpliceSite[] spliceChain1= null;	// sorted!
	SpliceSite[] spliceChain2= null;	// sorted!
	int degree= -1;
	ASEvent[] asEvents= null; // events
	
	public static class VarLengthFieldComparator implements Comparator {
		public int compare(Object arg0, Object arg1) {
			ASVariation[] v1= (ASVariation[]) arg0;
			ASVariation[] v2= (ASVariation[]) arg1;
			if (v1[0].toString().length()> v2[0].toString().length())
				return 1;
			if (v1[0].toString().length()> v2[0].toString().length())
				return -1;
			return 0;
		}
	}

	public int getLengthDiff(boolean exon) {
		int[] a= getLength(exon);
		int diffA= Math.abs(a[0]- a[1]);
		return diffA;
	}
	public int[] getLength(boolean exon) {
		
		int value= 0;
		SpliceSite[] borders= getBorderingSpliceSites();
		SpliceSite[][] schains2= new SpliceSite[][] {spliceChain1, spliceChain2};
		int[] result= new int[schains2.length];
		for (int i = 0; i < schains2.length; i++) {
			for (int j = 0; j < schains2[i].length; j++) {
				
				if ((schains2[i][j].isDonor()&& !exon)|| 	// check whether intron/exon relevant
						(schains2[i][j].isAcceptor()&& exon))
					continue;
				
				if (j== 0) {
					if (borders[0].isDonor()!= schains2[i][j].isDonor())
						System.err.println("unequal bordering ss's!");
					result[i]+= Math.abs(schains2[i][j].getPos()- borders[0].getPos());	// not +1 here, both are donors!
				} else {
					if (schains2[i][j-1].isDonor()== schains2[i][j].isDonor())
						System.err.println("other ss type requested (inconsistent splice chain)");
					result[i]+= schains2[i][j].getPos()- schains2[i][j-1].getPos()+ ((exon)?1:-1);	// intron -1, exon +1, draw
				}
			}
				// last case
			if ((borders[1].isAcceptor()&& !exon)|| (borders[1].isDonor()&& exon))
				continue;
			if (schains2[i]!= null&& schains2[i].length> 0)
				result[i]+= Math.abs(borders[1].getPos()- schains2[i][schains2[i].length- 1].getPos());
			else
				result[i]+= Math.abs(borders[1].getPos()- borders[0].getPos())+ ((exon)?-1:1);
		}
		
		return result;
	}
	
	/**
	 * returns the first and the last splice site in the splice universe
	 * @return
	 */
	public SpliceSite[] getBorderingSpliceSites() {
		SpliceSite[] su= getSpliceUniverse();
		return new SpliceSite[] {su[0], su[su.length- 1]};
	}

	
	public static ASVariation[][] sort2DForLengthOfVariation(ASVariation[][] vars) {
		Arrays.sort(vars, new VarLengthFieldComparator());
		return vars;
	}
	                            
	public static class SpliceStringComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			if (((ASVariation) arg0).toString().equals(((ASVariation) arg1).toString()))
				return 0;
			return -1;
		}
	
		/**
		 * Proceeds as follows:
		 * checks for the same start in each splice chain and then for the same 
		 * qualitative order of the splicesites in 2* both chains.
		 * 
		 * @return <code>0</code> if structure is equal, otherwise <code>-1</code>  
		 */
		public int compare_old(Object arg0, Object arg1) {
			
			SpliceSite[] sc11= ((ASVariation) arg0).getSpliceChain1();
			SpliceSite[] sc12= ((ASVariation) arg0).getSpliceChain2();
			SpliceSite[] sc21= ((ASVariation) arg1).getSpliceChain1();
			SpliceSite[] sc22= ((ASVariation) arg1).getSpliceChain2();
			
				// init corresponding splice chains
			SpliceSite[] first1= null;
			if (sc11!= null&& sc11.length> 0) {
				if (sc12!= null&& sc12.length> 0) {
					if (sc11[0].getPos()< sc12[0].getPos())
						first1= sc11;
					else
						first1= sc12;
				} else 
					first1= sc11;								
			} else
				first1= sc12;
			SpliceSite[] first2= null;
			if (sc21!= null&& sc21.length> 0) {
				if (sc22!= null&& sc22.length> 0) {
					if (sc21[0].getPos()< sc22[0].getPos())
						first2= sc21;
					else
						first2= sc22;
				} else 
					first2= sc21;								
			} else
				first2= sc22;
			HashMap correspond= new HashMap(2);
			correspond.put(first1, first2);
			SpliceSite[] scnd1= (first1==sc11)?sc12:sc11;
			SpliceSite[] scnd2= (first2==sc21)?sc22:sc21;
			correspond.put(scnd1, scnd2);
			
				// project splice chains into common array
			if (((first1== null)!= (first2== null))||((scnd1== null)!= (scnd2== null))||  
					(first1!= null&& first2!= null&& first1.length!= first2.length)|| 
					(first1!= null&& first2!= null&& first1.length> 0&& first2.length> 0&& first1[0].isDonor()!= first2[0].isDonor()) || 
					(scnd1!= null&& scnd2!= null&& scnd1.length!= scnd2.length)|| 
					(scnd1!= null&& scnd2!= null&& scnd1.length> 0&& scnd2.length> 0&& scnd1[0].isDonor()!= scnd2[0].isDonor())) 
				return (-1);
			Vector sc1Vec= new Vector(first1.length+ ((scnd1== null)?0:scnd1.length));
			int pos1= 0, pos2= 0;
			while((sc11!= null&& pos1< sc11.length)|| (sc12!= null&& pos2< sc12.length)) {
				if (sc11== null|| pos1>= sc11.length) {
					sc1Vec.add(sc12);
					++pos2;
				} else if (sc12== null|| pos2>= sc12.length) {
					sc1Vec.add(sc11);
					++pos1;
				} else if (sc11[pos1].getPos()< sc12[pos2].getPos()) {
					sc1Vec.add(sc11);
					++pos1;
				} else {
					sc1Vec.add(sc12);
					++pos2;
				}
			}
			Vector sc2Vec= new Vector(first2.length+ ((scnd2== null)?0:scnd2.length));
			pos1= 0; pos2= 0;
			while((sc21!= null&& pos1< sc21.length)|| (sc22!= null&& pos2< sc22.length)) { 
				if (sc21== null|| pos1>= sc21.length) {
					sc2Vec.add(sc22);
					++pos2;
				} else if (sc22== null|| pos2>= sc22.length) {
					sc2Vec.add(sc21);
					++pos1;
				} else if (sc21[pos1].getPos()< sc22[pos2].getPos()) {
					sc2Vec.add(sc21);
					++pos1;
				} else {
					sc2Vec.add(sc22);
					++pos2;
				} 
			}
			
			
				// main comparison loop
			while(pos1< sc1Vec.size()&& pos2< sc2Vec.size()) 
				if(correspond.get(sc1Vec.elementAt(pos1++))!= sc2Vec.elementAt(pos2++))
					return (-1);
			
			return 0;
		}
		
	}

	public static class SpliceChainComparator extends StructureComparator {
		
		public int compare(Object arg0, Object arg1) {
			
			SpliceSite[] s1= (SpliceSite[]) arg0;
			SpliceSite[] s2= (SpliceSite[]) arg1;
			
			Comparator compi= new SpliceSite.PositionComparator();
			if (s1== null|| s1.length< 1) {
				if (s2== null|| s2.length< 1) 
					return 0;
				else 
					return -1;
			} else {
				if (s2== null|| s2.length< 1)
					return 1;
				else {
					if (s1.length!= s2.length)
						return -1;
					for (int i = 0; i < s2.length; i++) {
						if (compi.compare(s1[i], s2[i])!= 0)
							return compi.compare(s1[i], s2[i]);
					}
					return 0;
				}
			}
		}
	}
	public static class UniqueCodingComparator extends StructureComparator {
		
		public int compare(Object arg0, Object arg1) {
			
			int eq= super.compare(arg0, arg1);
			if (eq!= 0)
				return -1;	// structurally different
			
			ASVariation as1= (ASVariation) arg0;
			ASVariation as2= (ASVariation) arg1;
//			if ((as1.isProteinCoding()!= as2.isProteinCoding())
//				|| (as1.isPartiallyCoding()!= as2.isPartiallyCoding())
//				|| (as1.isNotAtAllCoding()!= as2.isNotAtAllCoding()))
//				return -1;
			
				// priority list
			if (as1.isProteinCoding()^ as2.isProteinCoding()) 
				return 1;	// remove
			if (as1.isNotAtAllCoding()^ as2.isNotAtAllCoding()) 
				return -1;	// remove
			return 0;	// both protein coding, partially or not at all coding, keep one

		}
	}
	
	/**^
	 * @deprecated
	 *
	 */
	private static void trimAlternative() {
//		
//		int min= Integer.MAX_VALUE;
//		int max= Integer.MIN_VALUE;
//		for (int k = 0; k < ss2.length; k++) 
//			if (ss2[k].length> 0) {
//				min= Math.min(min, ss2[k][0].getPos());
//				max= Math.max(max, ss2[k][ss2[k].length- 1].getPos());
//			}
//		boolean omitStart= false, omitEnd= false;
//		// splice events only between conserved splice sites
//		if (transcripts[i].getPredSpliceSite(min)== null|| transcripts[j].getPredSpliceSite(min)== null)
//			omitStart= true;
//		if (transcripts[i].getSuccSpliceSite(max)== null|| transcripts[j].getSuccSpliceSite(max)== null)
//			omitEnd= true;
		
		// splice events only with splice sites covered by both transcripts
//		if (cluster[i].length> 0) {
//			if (cluster[i][0].getPos()< transcripts[j].getStart())
//				omitStart= true;
//			if (cluster[i][cluster[i].length- 1].getPos()> transcripts[j].getEnd())
//				omitEnd= true;
//		} else  {	// empty splice chain
//			if (cluster[j].length== 0|| 		// both 0-lenth
//					(cluster[j][cluster[j].length- 1].getPos()< transcripts[i].getStart()||	// check for out of range
//					cluster[j][0].getPos()> transcripts[i].getEnd()))
//					continue;
//		}
//		if (cluster[j].length> 0) {
//			if (cluster[j][0].getPos()< transcripts[i].getStart())
//				omitStart= true;
//			if (cluster[j][cluster[j].length- 1].getPos()> transcripts[i].getEnd())
//				omitEnd= true;
//		} else  {	// empty splice chain
//			if (cluster[i].length== 0|| 		// both 0-lenth
//					(cluster[i][cluster[i].length- 1].getPos()< transcripts[j].getStart()||	// check for out of range
//					cluster[i][0].getPos()> transcripts[j].getEnd()))
//					continue;
//		}
			
	}
	/**
	 * trims splice chains from both sides to the first splice site
	 * that is covered by exonic positions in all other splice chains.
	 *  
	 * @param vars
	 * @return
	 */
	public static SpliceSite[][] trim(SpliceSite[][] vars, Transcript[] trans) {
		
			// check for trimming left, right edge 
		boolean end5= ((vars[0].length> 0&& trans[0].getPredSpliceSite(vars[0][0])== null) ||
				(vars[1].length> 0&& trans[1].getPredSpliceSite(vars[1][0])== null))? true: false;
		boolean end3= ((vars[0].length> 0&& trans[0].getSuccSpliceSite(vars[0][vars[0].length- 1])== null) ||
				(vars[1].length> 0&& trans[1].getSuccSpliceSite(vars[1][vars[1].length- 1])== null))? true: false;				
		if (!end5&& !end3)	// nothing to do
			return vars;
		
			// create splice universe
		Vector suVec= new Vector();
		for (int i = 0; i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				int k;
				for (k = 0; k < suVec.size(); k++) 
					if (((SpliceSite) suVec.elementAt(k)).getPos()== vars[i][j].getPos())
						break;
				if (k== suVec.size())
					suVec.add(vars[i][j]);
			}
		}
		Vector tsVec= new Vector();	// add tss, tes
		if (end5) {
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].get5PrimeEdge())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].get5PrimeEdge()));
			}
		}
		if (end3) {
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].get3PrimeEdge())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].get3PrimeEdge()));
			}
		}
		int[] su= new int[suVec.size()+ tsVec.size()];	// convert to int[]
		for (int i = 0; i < suVec.size(); i++) 
			su[i]= ((SpliceSite) suVec.elementAt(i)).getPos();
		for (int i = 0; i < tsVec.size(); i++) 
			su[i+suVec.size()]= ((Integer) tsVec.elementAt(i)).intValue();
		
		if (su.length<= 0)
			return vars;	// empty splice chains
		
		Comparator compi= new AbstractSite.PositionComparator();
		Arrays.sort(su);
		int min= su[0];
		int max= su[su.length- 1];
		
			// trim end5
		while (end5) {
			// check for exonic area
			int i;
			for (i = 0; i < vars.length; i++) {
				if (Transcript.getSpliceSiteByPos(vars[i], min)== null) {	
//					if (vars[i]== null)	// if not containing edge site ----> bullshit, intron retention in 3'UTR
//						break;
					SpliceSite pred= Transcript.getPredSpliceSite(vars[i], min);
					SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], min);
					if (pred== null) {
						if (succ== null) {
							if (!trans[i].contains(min))
								break;
						} else {
							if (trans[i].isUpstream(min)|| !succ.isDonor())
								break;
						}
					} else {
						if (succ== null) {
							if (!pred.isAcceptor()|| trans[i].isDownstream(min))
								break;
						} else {
							if (!pred.isAcceptor()|| !succ.isDonor())
								break;
						}
					}
				}
			}
			if (i== vars.length)	// finish: all exonic
				break;
					
			int s= Transcript.getSuccPos(su, min);		// containing the splice site, accepted, next
			if (s== 0) {
				min= 0;
				break;
			}
			min= s;
		}
				
		// trim end3
		while (end3) {
			// check for exonic area
			int i;
			for (i = 0; i < vars.length; i++) {
				if (Transcript.getSpliceSiteByPos(vars[i], max)== null) {
//					if (vars[i].length== 0)	// cannot succeed in containing ----> bullshit, intron retention in 3'UTR
//						break;
					SpliceSite pred= Transcript.getPredSpliceSite(vars[i], max);
					SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], max);
					if (pred== null) {
						if (succ== null) {
							if (!trans[i].contains(max))
								break;
						} else {
							if (trans[i].isUpstream(max)|| !succ.isDonor())
								break;
						}
					} else {
						if (succ== null) {
							if (!pred.isAcceptor()|| trans[i].isDownstream(max))
								break;
						} else {
							if (!pred.isAcceptor()|| !succ.isDonor())
								break;
						}
					}
				}
			}
			if (i== vars.length)
				break;
					
			int s= Transcript.getPredPos(su, max);		// containing the splice site, accepted, next
			if (s== 0) {
				max= 0;
				break;
			}					
			max= s;
		}
		
		// cut
		if (min== 0|| max== 0) {
			return new SpliceSite[vars.length][0];
		}
		boolean tss= false, tes= false;	// check if trimmed flank contains a tss/tes
		for (int i = 0; i < vars.length; i++) {
			if (trans[i].get5PrimeEdge()== min)
				tss= true;
			if (trans[i].get3PrimeEdge()== max)
				tes= true;
		}
		SpliceSite[][] trimmed= new SpliceSite[vars.length][];
		for (int i = 0; i < vars.length; i++) {
			Vector splicVec= new Vector();
			int j= 0;
			if (end5) {
				for (j = 0; j < vars[i].length; j++) {
					if ((vars[i][j].getPos()> min)
					|| (tss&& vars[i][j].getPos()== min&& vars[i][j].isDonor()))	// conserve end5 end coincidence w tss
						break;
				}
			} else {
				for (j = 0; j < vars[i].length; j++) 
					if (vars[i][j].getPos()>= min)
						break;
			}
			
			if (end3)
				for (; j< vars[i].length&& vars[i][j].getPos()< max; j++) 
					splicVec.add(vars[i][j]);
			else 
				for (; j< vars[i].length&& vars[i][j].getPos()<= max; j++) 
					splicVec.add(vars[i][j]);
			if (end3&& tes&& j< vars[i].length&& vars[i][j].isAcceptor()&& vars[i][j].getPos()== max)
				splicVec.add(vars[i][j]);	// conserve end3 edge coincidence with tes
			trimmed[i]= (SpliceSite[]) gphase.tools.Arrays.toField(splicVec);
			if (trimmed[i]== null)
				trimmed[i]= new SpliceSite[0];
		}
		
		return trimmed;
	}
	/**
	 * trims splice chains from both sides to the first splice site
	 * that is covered by exonic positions in all other splice chains.
	 *  
	 * @param vars
	 * @return
	 */
	public static SpliceSite[][] trimOld(SpliceSite[][] vars, Transcript[] trans) {
		
			// check for trimming left, right edge 
		boolean left= ((vars[0].length> 0&& trans[0].getPredSpliceSite(vars[0][0])== null) ||
				(vars[1].length> 0&& trans[1].getPredSpliceSite(vars[1][0])== null))? true: false;
		boolean right= ((vars[0].length> 0&& trans[0].getSuccSpliceSite(vars[0][vars[0].length- 1])== null) ||
				(vars[1].length> 0&& trans[1].getSuccSpliceSite(vars[1][vars[1].length- 1])== null))? true: false;				
		if (!left&& !right)	// nothing to do
			return vars;
		
			// create splice universe
		Vector suVec= new Vector();
		for (int i = 0; i < vars.length; i++) {
			for (int j = 0; j < vars[i].length; j++) {
				int k;
				for (k = 0; k < suVec.size(); k++) 
					if (((SpliceSite) suVec.elementAt(k)).getPos()== vars[i][j].getPos())
						break;
				if (k== suVec.size())
					suVec.add(vars[i][j]);
			}
		}
		Vector tsVec= new Vector();	// add tss, tes
		if (left) {
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].getStart())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].getStart()));
			}
		}
		if (right) {
			for (int i = 0; i < trans.length; i++) {
				int k;
				for (k = 0; k < tsVec.size(); k++) 
					if (((Integer) tsVec.elementAt(k)).intValue()== trans[i].getEnd())
						break;
				if (k== tsVec.size())
					tsVec.add(new Integer(trans[i].getEnd()));
			}
		}
		int[] su= new int[suVec.size()+ tsVec.size()];
		for (int i = 0; i < suVec.size(); i++) 
			su[i]= ((SpliceSite) suVec.elementAt(i)).getPos();
		for (int i = 0; i < tsVec.size(); i++) 
			su[i+suVec.size()]= ((Integer) tsVec.elementAt(i)).intValue();
		
		if (su.length<= 0)
			return vars;	// empty splice chains
		
		Comparator compi= new AbstractSite.PositionComparator();
		Arrays.sort(su);
		int min= su[0];
		int max= su[su.length- 1];
		
			// trim left
		while (left) {
			// check for exonic area
			int i;
			for (i = 0; i < vars.length; i++) {
				if (Transcript.getSpliceSiteByPos(vars[i], min)== null) {	
					if (vars[i]== null)	// if not containing edge site
						break;
					SpliceSite pred= Transcript.getPredSpliceSite(vars[i], min);
					SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], min);
					if (pred== null) {
						if (succ== null) {
							if (min< trans[i].getStart()|| min> trans[i].getEnd())
								break;
						} else {
							if (min< trans[i].getStart()|| !succ.isDonor())
								break;
						}
					} else {
						if (succ== null) {
							if (!pred.isAcceptor()|| min> trans[i].getEnd())
								break;
						} else {
							if (!pred.isAcceptor()|| !succ.isDonor())
								break;
						}
					}
				}
			}
			if (i== vars.length)	// finish: all exonic
				break;
					
			int s= Transcript.getSuccPos(su, min);		// containing the splice site, accepted, next
			if (s== 0) {
				min= 0;
				break;
			}
			min= s;
		}
				
		// trim right
		while (right) {
			// check for exonic area
			int i;
			for (i = 0; i < vars.length; i++) {
				if (Transcript.getSpliceSiteByPos(vars[i], max)== null) {
					if (vars[i].length== 0)	// cannot succeed in containing
						break;
					SpliceSite pred= Transcript.getPredSpliceSite(vars[i], max);
					SpliceSite succ= Transcript.getSuccSpliceSite(vars[i], max);
					if (pred== null) {
						if (succ== null) {
							if (max< trans[i].getStart()|| max> trans[i].getEnd())
								break;
						} else {
							if (max< trans[i].getStart()|| !succ.isDonor())
								break;
						}
					} else {
						if (succ== null) {
							if (!pred.isAcceptor()|| max> trans[i].getEnd())
								break;
						} else {
							if (!pred.isAcceptor()|| !succ.isDonor())
								break;
						}
					}
				}
			}
			if (i== vars.length)
				break;
					
			int s= Transcript.getPredPos(su, max);		// containing the splice site, accepted, next
			if (s== 0) {
				max= 0;
				break;
			}					
			max= s;
		}
		
		// cut
		if (min== 0|| max== 0) {
			return new SpliceSite[vars.length][0];
		}
		boolean tss= false, tes= false;	// check if trimmed flank contains a tss/tes
		for (int i = 0; i < vars.length; i++) {
			if (trans[i].getStart()== min)
				tss= true;
			if (trans[i].getEnd()== max)
				tes= true;
		}
		SpliceSite[][] trimmed= new SpliceSite[vars.length][];
		for (int i = 0; i < vars.length; i++) {
			Vector splicVec= new Vector();
			int j= 0;
			if (left) {
				for (j = 0; j < vars[i].length; j++) {
					if ((vars[i][j].getPos()> min)
					|| (tss&& vars[i][j].getPos()== min&& vars[i][j].isDonor()))	// conserve left end coincidence w tss
						break;
				}
			} else {
				for (j = 0; j < vars[i].length; j++) 
					if (vars[i][j].getPos()>= min)
						break;
			}
			
			if (right)
				for (; j< vars[i].length&& vars[i][j].getPos()< max; j++) 
					splicVec.add(vars[i][j]);
			else 
				for (; j< vars[i].length&& vars[i][j].getPos()<= max; j++) 
					splicVec.add(vars[i][j]);
			if (right&& tes&& j< vars[i].length&& vars[i][j].isAcceptor()&& vars[i][j].getPos()== max)
				splicVec.add(vars[i][j]);	// conserve right edge coincidence with tes
			trimmed[i]= (SpliceSite[]) gphase.tools.Arrays.toField(splicVec);
			if (trimmed[i]== null)
				trimmed[i]= new SpliceSite[0];
		}
		
		return trimmed;
	}
	public boolean checkFramePreservation() {
		
		SpliceSite[] flanks= getFlankingSpliceSites();

		SpliceSite[] su1= new SpliceSite[spliceChain1.length+ 2];	// position arrays containing flanks
		for (int i = 1; i < su1.length- 1; i++) 
			su1[i]= spliceChain1[i-1];
		su1[0]= flanks[0];
		su1[su1.length- 1]= flanks[1];
		int exLen1= 0;
		for (int i = 1; i < su1.length; i++) 
			if (su1[i].isDonor()) {
				int x;
				if (i== 1&& su1[0]== null)
					x= trans1.getStart();
				else
					x= su1[i].getPos();
				exLen1+= su1[i].getPos()- x+ 1;
			}
		
		SpliceSite[] su2= new SpliceSite[spliceChain2.length+ 2];	// position arrays containing flanks
		for (int i = 1; i < su2.length- 1; i++) 
			su2[i]= spliceChain2[i-1];
		su2[0]= flanks[0];
		su2[su2.length- 1]= flanks[1];
		int exLen2= 0;
		for (int i = 1; i < su2.length; i++) 
			if (su2[i].isDonor())
				exLen2+= su2[i].getPos()- su2[i].getPos()+ 1;
		
		if ((exLen1% 3)!= (exLen2% 3))
			return false;
		return true;		
	}

	public static class StructureComparator implements Comparator {
		
		/**
		 * @return <code>0</code> if both objects are equal, <code>1</code>
		 * otherwise
		 */
		public int compare_old(Object arg0, Object arg1) {
			
			ASVariation as1= (ASVariation) arg0;
			ASVariation as2= (ASVariation) arg1;
			
				// find analogs
			SpliceSite[] sc1= null;
			if (as1.spliceChain1!= null && as1.spliceChain1.length> 0) 
				sc1= as1.spliceChain1;
			else
				sc1= as1.spliceChain2;	// both cannot be empty
			SpliceSite[] sc2= as1.getOtherSpliceChain(sc1);
			
			SpliceSite[] sc1analog;
			if (as2.spliceChain1!= null&& as2.spliceChain1.length> 0
					&& as2.spliceChain1[0]== sc1[0])	// only one sc of as2 can start with the same splice site than sc1 !!
				sc1analog= as2.spliceChain1;
			else
				sc1analog= as2.spliceChain2;
			SpliceSite[] sc2analog= as2.getOtherSpliceChain(sc1analog);

				// compare analogs
			Comparator compi= new AbstractSite.PositionComparator();
			if (sc1.length!= sc1analog.length)
				return (-1);
			for (int i = 0; i < sc1.length; i++)	// first pair cannot be empty 
				for (int j = 0; j < sc1analog.length; j++) 
					if (compi.compare(sc1[i], sc1analog[i])!= 0)
						return (-1);

			if ((sc2== null^ sc2analog== null))		// catch nullpointers first..
				return (-1);
			if (sc2== null&& sc2analog== null)
				return 0;
			if (sc2.length!= sc2analog.length)
				return (-1);
			for (int i = 0; i < sc2.length; i++)	 
				for (int j = 0; j < sc2analog.length; j++) 
					if (compi.compare(sc2[i], sc2analog[i])!= 0)
						return (-1);
			return 0;
		}

		/**
		 * @return <code>0</code> if both objects are equal, <code>-1</code>
		 * otherwise
		 */
		public int compare(Object arg0, Object arg1) {
			
			ASVariation as1= (ASVariation) arg0;
			ASVariation as2= (ASVariation) arg1;
			
			Comparator compi= new SpliceChainComparator();
			SpliceSite[][] s1= new SpliceSite[][] {as1.spliceChain1, as1.spliceChain2};
			SpliceSite[][] s2= new SpliceSite[][] {as2.spliceChain1, as2.spliceChain2};
			SpliceSite[] su1= as1.getSpliceUniverse();
			SpliceSite[] su2= as2.getSpliceUniverse();
			if (((compi.compare(s1[0], s2[0])== 0)&& (compi.compare(s1[1], s2[1])== 0))||
				((compi.compare(s1[1], s2[0])== 0)&& (compi.compare(s1[0], s2[1])== 0)))
				return 0;
			return -1;
			
			// compare splice universees
			// night of 6.6.06: NO! dbl exon skip/mut exclusive			
//			SpliceSite[] su1= as1.getSpliceUniverse();
//			SpliceSite[] su2= as2.getSpliceUniverse();
//			if (su1.length!= su2.length)
//				return -1;
			
				// check also flanking sites (not tss, tes)
				// night of 6.6.06: no longer check flanking sites
//			SpliceSite[] flank= as1.getFlankingSpliceSites();
//			if (flank[0]!= null)
//				su1= (SpliceSite[]) gphase.tools.Arrays.add(su1, flank[0]);
//			if (flank[1]!= null)
//				su1= (SpliceSite[]) gphase.tools.Arrays.add(su1, flank[1]);
//			Comparator compi= new AbstractSite.PositionComparator();
//			Arrays.sort(su1, compi);
//			flank= as2.getFlankingSpliceSites();
//			if (flank[0]!= null)
//				su2= (SpliceSite[]) gphase.tools.Arrays.add(su2, flank[0]);
//			if (flank[1]!= null)
//				su2= (SpliceSite[]) gphase.tools.Arrays.add(su2, flank[1]);
//			Arrays.sort(su2, compi);
//			if (su1.length!= su2.length)
//				return -1;
			
//			for (int i = 0; i < su2.length; i++) 
//				if (su1[i].getPos()!= su2[i].getPos())
//					return -1;
//			return 0;
		}
	}
	public static class CodingHierarchyFilter extends StructureComparator {
			
			public int compare(Object arg0, Object arg1) {
				
				int eq= super.compare(arg0, arg1);
				if (eq!= 0)
					return -1;	// structurally different
				
				ASVariation as1= (ASVariation) arg0;
				ASVariation as2= (ASVariation) arg1;
	//			if ((as1.isProteinCoding()!= as2.isProteinCoding())
	//				|| (as1.isPartiallyCoding()!= as2.isPartiallyCoding())
	//				|| (as1.isNotAtAllCoding()!= as2.isNotAtAllCoding()))
	//				return -1;
				
					// priority list
				String t11= as1.getTranscript1().getTranscriptID();
				String t12= as1.getTranscript2().getTranscriptID();
				String t21= as2.getTranscript1().getTranscriptID();
				String t22= as2.getTranscript2().getTranscriptID();
				if (as1.isProteinCoding()) 
					return (as2.isProteinCoding())?0:1;	// keep only first
				if (as2.isProteinCoding())
					return 2;		// keep second
				if (as1.isPartiallyCoding())
					return as2.isPartiallyCoding()?0:1;	// keep first
				if (as2.isPartiallyCoding())
					return 2;
				return 0;	// keep one (whichever)
			}
		}

	public static class CodingComparator extends StructureComparator {
			
			public int compare(Object arg0, Object arg1) {
				
				int eq= super.compare(arg0, arg1);
				if (eq!= 0)
					return -1;	// structurally different
				
				ASVariation as1= (ASVariation) arg0;
				ASVariation as2= (ASVariation) arg1;
	//			if ((as1.isProteinCoding()!= as2.isProteinCoding())
	//				|| (as1.isPartiallyCoding()!= as2.isPartiallyCoding())
	//				|| (as1.isNotAtAllCoding()!= as2.isNotAtAllCoding()))
	//				return -1;
				
					// priority list
				if ((as1.isProteinCoding()&& as2.isProteinCoding())||
					(as1.isPartiallyCoding()&& as2.isPartiallyCoding())||
					(as1.isNotAtAllCoding()&& as2.isNotAtAllCoding()))
					return 0;
				return 1;	
			}
		}

	public Gene getGene() {
		Gene g1= trans1.getGene();
		Gene g2= trans2.getGene();
		if (g1!= g2)
			System.err.println("Gene mismatch between "+trans1+ " and "+trans2+"!");
		return g1;
	}
	/**
	 * Checks for OID.. 
	 * @param schain
	 * @return
	 */
	public SpliceSite[] getOtherSpliceChain(SpliceSite[] schain) {
		if (schain== spliceChain1)
			return spliceChain2;
		if (schain!= spliceChain2)
			System.err.println("SpliceChain "+schain+ " does not match ASVariation!");
		return spliceChain1;
	}
	public ASVariation(Transcript newTID1, Transcript newTID2, SpliceSite[] newSChain1, SpliceSite[] newSChain2) {		
		trans1= newTID1;
		trans2= newTID2;		
		Arrays.sort(newSChain1, new SpliceSite.PositionComparator());
		Arrays.sort(newSChain2, new SpliceSite.PositionComparator());
		spliceChain1= newSChain1;
		spliceChain2= newSChain2;
		markAS(spliceChain1);
		markAS(spliceChain2);
	}
	
	void markAS(SpliceSite[] schain) {
		for (int i = 0; i < schain.length; i++) {
			schain[i].setConstitutive(false);
			//schain[i].addASVar(this);
		}
	}
	
	public SpliceSite getSpliceSite(int pos) throws IllegalArgumentException {
		
		SpliceSite ss1= null;
		if (trans1!= null)
			ss1= trans1.getSpliceSite(pos);
		SpliceSite ss2= trans2.getSpliceSite(pos);
		if (ss1!= null&& ss2!= null&& ss1!= ss2) {
			if (ss1.isDonor()!= ss2.isDonor()) {
				System.err.println("Donor/acceptor mix!");
				throw new IllegalArgumentException();
			} else 
				System.err.println("Inconsitency in SpliceSite objects: multiple instances of same object!");
		}
		return ((ss1== null)?ss2:ss1);
	}
	
	public boolean contains(int[][] pattern){
		
		if (pattern== null|| pattern.length!= 2)		
			return false;
		
		
	}
	
	public int[][] getSSRelPos() {
		
//		int[][] pos= new int[2][];
//		pos[0]= (spliceChain1== 0)? new int[0]:new int[spliceChain1.length];
//		pos[1]= (spliceChain2== 0)? new int[0]:new int[spliceChain2.length];
//			
//		int ltt= 1;
//		int c1= 0, c2= 0;
//		int p1= 0, p2= 0;
//		while ((spliceChain1!= null&& p1< spliceChain1.length)||
//				(spliceChain2!= null&& p2< spliceChain2.length)) {
//			if (spliceChain1!= null&& p1< spliceChain1.length)
//				if (spliceChain2!= null&& p2< spliceChain2.length)
//					if (spliceChain1[p1].getPos()< spliceChain2[p2].getPos())
//						pos[c1++]= ltt;
//					else
//						c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
//				else
//					c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
//			else
//				c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
//				
//			ltt++;
//		}
		return null;
	}
	
			/**
			 * A splice variation is defined to be protein coding if the 
			 * constitutive flanking splice sites lie within any translation
			 * of the respective transcript.
			 */
			public boolean isProteinCoding() {
			
				SpliceSite[] su= getSpliceUniverse();
				if (trans1.getTranslation(su[0].getPos(), su[su.length- 1].getPos())== null)
						return false;
				if (trans2.getTranslation(su[0].getPos(), su[su.length- 1].getPos())== null)
						return false;
		
				return true;
			}
			
			/**
			 * A variation that is covered fully by a CDS in exactly one of the transcripts
			 */
			public boolean isSemiProteinCoding() {
				SpliceSite[] su= getSpliceUniverse();			
				if ((trans1.getTranslation(su[0].getPos(), su[su.length- 1].getPos())== null)^
					(trans2.getTranslation(su[0].getPos(), su[su.length- 1].getPos())== null))
						return true;
		
				return false;
			}
			
			public boolean checkUniqueCoding() {
				if (isProteinCoding()^ isNotAtAllCoding()^ isPartiallyCoding()) {
					return true;
				}
				System.err.println("inconsistent coding type "+isProteinCoding()+","+isPartiallyCoding()+","+isNotAtAllCoding());
				outputDetail(System.err);
				if (trans1.translations!= null)
					System.err.println(trans1.translations[0]);
				else
					System.err.println("null");
				if (trans2.translations!= null)
					System.err.println(trans2.translations[0]);
				else
					System.err.println("null");
				return false;
			}
	
			/**
			 * A splice variation is defined to be protein coding if the 
			 * constitutive flanking splice sites lie within any translation
			 * of the respective transcript.
			 */
			public boolean isProteinCoding_new() {
			
				// dropped, new criterion
//				int[] flanks= getFlankingSpliceSites();
//				int x= flanks[0];
//				int y= flanks[1];
//				Translation[] t1= trans1.getTranslation(x, y);
//				Translation[] t2= trans2.getTranslation(x, y);
//				if (t1== null^ t2== null)
//					;// System.out.println("just one protein coding in AS");
//				if (t1== null|| t2== null)
//					return false;
//				return true;
				
				System.currentTimeMillis();
				int x,y;
				if (spliceChain1!= null && spliceChain1.length> 0) {
					if (trans1.isForward()) { 
						x= spliceChain1[0].isDonor()?spliceChain1[0].getPos()-1:spliceChain1[0].getPos()+2;
						y= spliceChain1[spliceChain1.length- 1].isDonor()?spliceChain1[spliceChain1.length- 1].getPos()-1:
							spliceChain1[spliceChain1.length- 1].getPos()+2;
					} else {
						x= spliceChain1[0].isDonor()?spliceChain1[0].getPos()+2:spliceChain1[0].getPos()-1;
						y= spliceChain1[spliceChain1.length- 1].isDonor()?spliceChain1[spliceChain1.length- 1].getPos()+2:
							spliceChain1[spliceChain1.length- 1].getPos()-1;
					}
					if (trans1.getTranslation(x,y)== null)
						return false;
				}

				if (spliceChain2!= null && spliceChain2.length> 0) { 
						
					if (trans2.isForward()) { 
						x= spliceChain2[0].isDonor()?spliceChain2[0].getPos()-1:spliceChain2[0].getPos()+2;
						y= spliceChain2[spliceChain2.length- 1].isDonor()?spliceChain2[spliceChain2.length- 1].getPos()-1:
							spliceChain2[spliceChain2.length- 1].getPos()+2;
					} else {
						x= spliceChain2[0].isDonor()?spliceChain2[0].getPos()+2:spliceChain2[0].getPos()-1;
						y= spliceChain2[spliceChain2.length- 1].isDonor()?spliceChain2[spliceChain2.length- 1].getPos()+2:
							spliceChain2[spliceChain2.length- 1].getPos()-1;
					}
					if (trans2.getTranslation(x,y)== null)
						return false;
				}

				return true;
			}

	public SpliceSite[] getFlankingSpliceSites()  {
		SpliceSite pred= null;
		SpliceSite succ= null;
		if (spliceChain1!= null&& spliceChain1.length> 0) {	// one of both transcripts has to provide a non-empty splicechain
			pred= trans1.getPredSpliceSite(spliceChain1[0]);			
			succ= trans1.getSuccSpliceSite(spliceChain1[spliceChain1.length- 1]);
		} 
		if (spliceChain2!= null&& spliceChain2.length> 0) {
			SpliceSite pred2= trans2.getPredSpliceSite(spliceChain2[0]);
			SpliceSite succ2= trans2.getSuccSpliceSite(spliceChain2[spliceChain2.length- 1]);
			if (pred== null)
				pred= trans2.getPredSpliceSite(spliceChain2[0]);
			else
				if (pred2!= null&& pred.getPos()!= pred2.getPos())	// can happen now, due to new trimming
					;// System.err.println("Not same flank "+trans2.getGene().getStableID());
			if (succ== null)
				succ= trans2.getSuccSpliceSite(spliceChain2[spliceChain2.length- 1]);
			else
				if (succ2!= null&& succ.getPos()!= succ2.getPos())
					;// System.err.println("Not same flank "+trans1.getGene().getStableID());
				
		}
		
		SpliceSite x= null, y= null;
		SpliceSite[] su= getSpliceUniverse();
		if (pred== null)
			;//System.err.println("no predecessor!");	//x= su[0].getPos()- 1;
		else
			x= pred;
		if (succ== null)
			;//System.err.println("no successor!");	//y= su[su.length- 1].getPos()- 1;
		else
			y= succ;

		return new SpliceSite[] {x,y};
	}
	/**
	 * A splice variation is defined to be strictly not protein coding if the 
	 * constitutive flanking splice sites lie both outside of any translation.
	 */
	public boolean isNotAtAllCoding() {

//		int[] flanks= getFlankingSpliceSites();
//		int x= flanks[0];
//		int y= flanks[1];
//		
//		Translation[] t1= trans1.getTranslation(x, y);
//		Translation[] t2= trans2.getTranslation(x, y);
//		if (t1== null&& t2== null) {
//			return true;
//		}
//		return false;
		
		SpliceSite[] su= getSpliceUniverse();
		for (int i = su[0].getPos(); i <= su[su.length- 1].getPos(); i++) {
			if (trans1.isCDS(i)|| trans2.isCDS(i))
				return false;
		}

		return true;

	}
	
	public boolean isTrue() {
		return true;
	}

	/**
	 * A splice variation is defined to be PARTIALLY PROTEIN CODING
	 * if at least one transcript provides a CDS area and 
	 * at least one transcript provides a non-CDS area.
	 */
	public boolean isPartiallyCoding() {

//		int[] flanks= getFlankingSpliceSites();
//		int x= flanks[0];
//		int y= flanks[1];
		
		SpliceSite[] su= getSpliceUniverse();
		int x= su[0].getPos();
		int y= su[su.length- 1].getPos();
		
			// look for coding part
		int i;
		for (i = x; i <= y; i++) 
			if (trans1.isCDS(i)|| trans2.isCDS(i))
				break;
		boolean cds= (i<= y)? true: false;
		
		// look for non coding part
		for (i = x; i <= y; i++) 
			if (!trans1.isCDS(i)|| !trans2.isCDS(i))
				break;
		boolean not= (i<= y)? true: false;
		

		if (cds&& not) {
			return true;
		}
		return false;
	}
	
	
	/**
	 * @return <code>true</code> if it includes the first splice site in both 
	 * transcripts
	 */
	public boolean includesFirstSpliceSite() {
		Comparator compi= new SpliceSite.PositionComparator();
		if (spliceChain1!= null&& spliceChain1.length> 0&& compi.compare(spliceChain1[0], trans1.getSpliceChain()[0])!= 0)
			return false;
		if (spliceChain2!= null&& spliceChain2.length> 0&& compi.compare(spliceChain2[0], trans2.getSpliceChain()[0])!= 0)
			return false;
		return true;
	}
	
	/**
	 * @return <code>true</code> if it includes the last splice site in both 
	 * transcripts
	 */
	public boolean includesLastSpliceSite() {
		Comparator compi= new SpliceSite.PositionComparator();
		if (spliceChain1!= null&& spliceChain1.length> 0&& 
				compi.compare(spliceChain1[0], trans1.getSpliceChain()[trans1.getSpliceChain().length])!= 0)
			return false;
		if (spliceChain2!= null&& spliceChain2.length> 0&& 
				compi.compare(spliceChain2[0], trans2.getSpliceChain()[trans2.getSpliceChain().length])!= 0)
			return false;
		return true;
	}	
	
	public boolean is5UTR() {

		SpliceSite pred= null;
		SpliceSite succ= null;
		if (spliceChain1!= null&& spliceChain1.length> 0) {	// one of both transcripts has to provide a non-empty splicechain
			pred= trans1.getPredSpliceSite(spliceChain1[0]);			
			succ= trans1.getSuccSpliceSite(spliceChain1[spliceChain1.length- 1]);
		} 
		if (spliceChain2!= null&& spliceChain2.length> 0) {
			if (pred== null)
				pred= trans2.getPredSpliceSite(spliceChain2[0]);
			else
				if (pred.getPos()!= trans2.getPredSpliceSite(spliceChain2[0]).getPos())
					System.err.println("Not same flank "+trans2.getGene().getStableID());
			if (succ== null)
				succ= trans2.getSuccSpliceSite(spliceChain2[spliceChain2.length- 1]);
			else
				if (succ.getPos()!= trans2.getSuccSpliceSite(spliceChain2[spliceChain2.length- 1]).getPos())
					System.err.println("Not same flank "+trans1.getGene().getStableID());
				
		}
		
		int x, y;
		SpliceSite[] su= getSpliceUniverse();
		if (pred== null)
			x= su[0].getPos()- 1;
		else
			x= pred.getPos();
		if (succ== null)
			y= su[su.length- 1].getPos()- 1;
		else
			y= succ.getPos();
		
		
		if (trans1.is5UTR(x)&& trans2.is5UTR(x)&& trans1.is5UTR(y)&& trans2.is5UTR(y))
			return true;
		return false;
	}	
	
	/**
	 * @deprecated check
	 * @return
	 */
	// muy feo :(
	// may the force be with you!
	public ASEvent[] getASEvents() {

		if (asEvents == null) {
			
			Comparator compi= new SpliceSite.PositionComparator();
			SpliceSite[] t1= trans1.getSpliceChain();	// complete chains of splice sites
			SpliceSite[] t2= trans2.getSpliceChain();
			
			
				// get exon areas
			Vector exVec1= new Vector(), exVec2= new Vector();
			if (spliceChain1!= null&& spliceChain1.length> 0&& spliceChain1[0].isDonor()) 	// first exons
				exVec1.add(new DummyRegion(trans1.getOtherSideOfExon(
						spliceChain1[0].getPos()), 
						spliceChain1[0].getPos()));
			if (spliceChain2!= null&& spliceChain2.length> 0&& spliceChain2[0].isDonor()) 	
				exVec2.add(new DummyRegion(trans2.getOtherSideOfExon(
						spliceChain2[0].getPos()), 
						spliceChain2[0].getPos()));
			
			for (int i = 0; i < spliceChain1.length; i++) {		// following exons
				if (spliceChain1[i].isDonor())
					continue;
				
				int start= spliceChain1[i].getPos();
				int end= -1;
				if (i+1< spliceChain1.length)
					end= spliceChain1[i+1].getPos();
				else
					end= trans1.getOtherSideOfExon(spliceChain1[i].getPos());	// last exon
				exVec1.add(new DummyRegion(start, end)); 
			}
			for (int i = 0; i < spliceChain2.length; i++) {
				if (spliceChain2[i].isDonor())
					continue;
				
				int start= spliceChain2[i].getPos();
				int end= -1;
				if (i+1< spliceChain2.length)
					end= spliceChain2[i+1].getPos();
				else
					end= trans2.getOtherSideOfExon(spliceChain2[i].getPos());	// last exon
				exVec2.add(new DummyRegion(start, end)); 
			}
			
			
			
				// generate events
			Vector aseVec= new Vector();
			for (int i = 0; i < exVec1.size(); i++) {
				
				Vector overlap= new Vector();
				for (int j = 0; j < exVec2.size(); j++)			// collect all overlapping 
					if (((DummyRegion) exVec1.elementAt(i)).overlaps((DummyRegion) exVec2.elementAt(j)))
						overlap.add(exVec2.elementAt(j));
		        
				if (overlap.size()< 1) { 
					if (!(trans1.getTSSPos()== ((Region) exVec1.elementAt(i)).getStart())	// exclude TSS / TSE
						&& !(trans1.getTESPos()== ((Region) exVec1.elementAt(i)).getEnd())
						&& !(((Region) exVec1.elementAt(i)).getStart()< trans2.getTSSPos())	// exclude alternate TSS /TSE
						&& !(((Region) exVec1.elementAt(i)).getEnd()> trans2.getTESPos())) {			// exon skipping or mutually exclusive
							try {
								aseVec.add(new ASEvent(this, 1, 
										getSpliceSite(((Region) exVec1.elementAt(i)).getStart()), 
										getSpliceSite(((Region) exVec1.elementAt(i)).getEnd())));
								if (aseVec.size()> 1&& 
										((ASEvent) aseVec.elementAt(aseVec.size()- 2)).isExonSkipping()&&
										((ASEvent) aseVec.elementAt(aseVec.size()- 2)).getTranscriptNr()== 2) {	// mutually exclusive
									((ASEvent) aseVec.elementAt(aseVec.size()- 2)).setType(ASEvent.TYPE_MUTUALLY_EXCLUSIVE);
									((ASEvent) aseVec.elementAt(aseVec.size()- 1)).setType(ASEvent.TYPE_MUTUALLY_EXCLUSIVE);
								}
							} catch (IllegalArgumentException e) {;}
					}						
					continue;
				} 
				
				if (overlap.size()> 1) {			// intron retention (1+) and maybe alternate donor/acceptor (at start/end of chain)
					for (int j = 0; j < overlap.size()- 1; j++) 	
						try {
							aseVec.add(new ASEvent(this, 2, 			// 1+ intron retention
									getSpliceSite(((Region) overlap.elementAt(j)).getEnd()), 
									getSpliceSite(((Region) overlap.elementAt(j+1)).getStart())));
						} catch (IllegalArgumentException e) {;}
					if (((Region) exVec1.elementAt(i)).getStart()!= ((Region) overlap.elementAt(0)).getStart()
							&& !(trans1.getTSSPos()== ((Region) exVec1.elementAt(i)).getStart())
							&& !(trans2.getTSSPos()== ((Region) overlap.elementAt(0)).getStart()))
						try {
							aseVec.add(new ASEvent(this, 0, 									// alternate acceptor
									getSpliceSite(((Region) exVec1.elementAt(i)).getStart()), 
									getSpliceSite(((Region) overlap.elementAt(0)).getStart())));
						} catch (IllegalArgumentException e) {;}
					if (((Region) exVec1.elementAt(i)).getEnd()!= ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()
							&& !(trans1.getTESPos()== ((Region) exVec1.elementAt(i)).getEnd())
							&& !(trans2.getTESPos()== ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()))
						try {
							aseVec.add(new ASEvent(this, 0, 									// alternate donor
									getSpliceSite(((Region) exVec1.elementAt(i)).getEnd()), 
									getSpliceSite(((Region) overlap.elementAt(overlap.size()- 1)).getEnd())));
						} catch (IllegalArgumentException e) {;}
					continue;
				}
				
					// (overlap.size()== 1): alternate donor/acceptor
				Vector overlap2= new Vector();	// check for intron retention in other sequence
				for (int k = 0; k < exVec1.size(); k++) 
					if (((DummyRegion) overlap.elementAt(0)).overlaps((DummyRegion) exVec1.elementAt(k)))
						overlap2.add(exVec1.elementAt(k));
				if (overlap2.size()> 1)
					continue;
				
				if (((Region) exVec1.elementAt(i)).getStart()!= ((Region) overlap.elementAt(0)).getStart()
					&& !(trans1.getTSSPos()== ((Region) exVec1.elementAt(i)).getStart())
					&& !(trans2.getTSSPos()== ((Region) overlap.elementAt(0)).getStart()))
					try {
						aseVec.add(new ASEvent(this, 0, 			// alternate acceptor
								getSpliceSite(((Region) exVec1.elementAt(i)).getStart()), 
								getSpliceSite(((Region) overlap.elementAt(0)).getStart())));
					} catch (IllegalArgumentException e) {;}
					
						
				if (((Region) exVec1.elementAt(i)).getEnd()!= ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()
					&& !(trans1.getTESPos()== ((Region) exVec1.elementAt(i)).getEnd())
					&& !(trans2.getTESPos()== ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()))
					try {
						aseVec.add(new ASEvent(this, 0, 			// alternate donor
								getSpliceSite(((Region) exVec1.elementAt(i)).getEnd()), 
								getSpliceSite(((Region) overlap.elementAt(overlap.size()- 1)).getEnd())));
					} catch (IllegalArgumentException e) {;}
			}

				// check exon skippings in other sequence
			for (int i = 0; i < exVec2.size(); i++) {
				
				Vector overlap= new Vector();
				for (int j = 0; j < exVec1.size(); j++)			// collect all overlapping 
					if (((DummyRegion) exVec2.elementAt(i)).overlaps((DummyRegion) exVec1.elementAt(j)))
						overlap.add(exVec1.elementAt(j));
		        
				if ((overlap.size()< 1) 
					&& !(trans2.getTSSPos()== ((Region) exVec2.elementAt(i)).getStart())
					&& !(trans2.getTESPos()== ((Region) exVec2.elementAt(i)).getEnd())
					&& !(((Region) exVec2.elementAt(i)).getStart()< trans1.getTSSPos())	// exclude alternate TSS /TSE
					&& !(((Region) exVec2.elementAt(i)).getEnd()> trans1.getTESPos())) {			// exon skipping or mutually exclusive
					
					try {
						ASEvent ev= new ASEvent(this, 2, 
								getSpliceSite(((Region) exVec2.elementAt(i)).getStart()), 
								getSpliceSite(((Region) exVec2.elementAt(i)).getEnd()));
						compi= new ASEvent.PositionComparator();
						int b= Arrays.binarySearch(overlap.toArray(), ev, compi);
						if (b> 0) {
							System.err.println("exon skip event already contained");
							continue;
						}
						b= (b+ 1)* (-1);
						aseVec.insertElementAt(ev, b);
						if (b> 0&& 
								((ASEvent) aseVec.elementAt(b- 1)).isExonSkipping()&&
								((ASEvent) aseVec.elementAt(b- 1)).getTranscriptNr()== 1) {	// mutually exclusive
							((ASEvent) aseVec.elementAt(b- 1)).setType(ASEvent.TYPE_MUTUALLY_EXCLUSIVE);
							((ASEvent) aseVec.elementAt(b)).setType(ASEvent.TYPE_MUTUALLY_EXCLUSIVE);
						}
					} catch (IllegalArgumentException e) {;}

				}
				
				if (overlap.size()> 1) {			// intron retention (1+) and maybe alternate donor/acceptor (at start/end of chain)
					for (int j = 0; j < overlap.size()- 1; j++)
						try {
							aseVec.add(new ASEvent(this, 1, 			// 1+ intron retention
									getSpliceSite(((Region) overlap.elementAt(j)).getEnd()), 
									getSpliceSite(((Region) overlap.elementAt(j+1)).getStart())));
						} catch (IllegalArgumentException e) {;}
					if (((Region) exVec2.elementAt(i)).getStart()!= ((Region) overlap.elementAt(0)).getStart()
							&& !(trans2.getTSSPos()== ((Region) exVec2.elementAt(i)).getStart())
							&& !(trans1.getTSSPos()== ((Region) overlap.elementAt(0)).getStart()))
						try {
							aseVec.add(new ASEvent(this, 0, 									// alternate acceptor
									getSpliceSite(((Region) exVec2.elementAt(i)).getStart()), 
									getSpliceSite(((Region) overlap.elementAt(0)).getStart())));
						} catch (IllegalArgumentException e) {;}
					if (((Region) exVec2.elementAt(i)).getEnd()!= ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()
							&& !(trans2.getTESPos()== ((Region) exVec2.elementAt(i)).getEnd())
							&& !(trans1.getTESPos()== ((Region) overlap.elementAt(overlap.size()- 1)).getEnd()))
						try {
							aseVec.add(new ASEvent(this, 0, 									// alternate donor
									getSpliceSite(((Region) exVec2.elementAt(i)).getEnd()), 
									getSpliceSite(((Region) overlap.elementAt(overlap.size()- 1)).getEnd())));
						} catch (IllegalArgumentException e) {;}
						
					continue;
				}
				
			}				
			
			asEvents = ASEvent.toArray(aseVec); 
		}

		return asEvents;
	}
	
	public String toASEventString() {
		String result= "";
		for (int i = 0; i < getASEvents().length; i++) 
			result+= getASEvents()[i];
		
		return result;
	}
	
	
	public boolean containsTranscripts(Transcript tID1, Transcript tID2) {
		if ((tID1== trans1&& tID2== trans2)||
				(tID2== trans1&& tID1== trans2))
			return true;
		return false;
	}
	
	public void outputDetail(PrintStream stream) {
		
		String s= null;
		if (isProteinCoding())
			s= "prot";
		else if (isNotAtAllCoding())
			s= "non";
		else if (isPartiallyCoding())
			s= "part";
		stream.println(toString()+"\t"+trans1+"\t"+trans2+"\t"+s);
		
		stream.print(getGene().getChromosome()+"\t");
		for (int j = 0; j < spliceChain1.length; j++) 
			stream.print(spliceChain1[j]+" ");
		stream.println();
		
		stream.print(getGene().getChromosome()+"\t");
		for (int j = 0; j < spliceChain2.length; j++) 
			stream.print(spliceChain2[j]+" ");		
		stream.println();
	}
	
	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toString() {
		String c1= "";
		String c2= "";
		
		int ltt= 1;
		int p1= 0, p2= 0;
		while ((spliceChain1!= null&& p1< spliceChain1.length)||
				(spliceChain2!= null&& p2< spliceChain2.length)) {
			if (spliceChain1!= null&& p1< spliceChain1.length)
				if (spliceChain2!= null&& p2< spliceChain2.length)
					if (spliceChain1[p1].getPos()< spliceChain2[p2].getPos())
						c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
					else
						c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
				else
					c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
			else
				c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "=";
				
			ltt++;
		}

		if (c2.length()< c1.length()) {	// print shortest first
			String h= c1;
			c1= c2;
			c2= h;
		} else if (c1.length()== c2.length()) {	// sort according to first number
			int x= 0;
			while (Character.isDigit(c1.charAt(x++)));
			int n1= Integer.parseInt(c1.substring(0, x-1));
			x= 0;
			while (Character.isDigit(c2.charAt(x++)));
			int n2= Integer.parseInt(c2.substring(0, x-1));
			if (n2< n1)  {
				String h= c1;
				c1= c2;
				c2= h;
			}
		}
		return "("+ c1+ " // "+ c2+ ")";
	}
	
	/**
	 * @return a chain of symbols representing the single AS events contained in 
	 * the variation.
	 */
	public String toStringASEvents() {
		
		String result= "";
		for (int i = 0; i < getASEvents().length; i++) 
			result+= asEvents[i]+ "-";
		if (result.length()> 1)
			result= result.substring(0, result.length()- 1);
		
		return result;
	}
	
	public int getDegree() {
		int x= 0; 
		if (spliceChain1!= null)
			x+= spliceChain1.length;
		if (spliceChain2!= null)
			x+= spliceChain2.length;
		return x;
	}
	
	public String toCodingInfoString() {
		if (isProteinCoding())
			return "prot";
		if (isPartiallyCoding())
			return "part";
		if (isNotAtAllCoding())
			return "not";
		return "";
	}
	
	public String toBitString() {
		int[][] bMatrix= toBitMatrix();
		String s= "";
		for (int i = 0; i < bMatrix.length; i++) {
			if (bMatrix[i][0]== 0) {
				if (bMatrix[i][1]== 0)
					s+= "A";
				else
					s+= "B";
			}  else {
				if (bMatrix[i][1]== 0)
					s+= "C";
				else
					s+= "D";
			}
		}
		
		return s;
	}
	
	public SpliceSite[][] getSpliceSites(int bitStart, int bitEnd) {
		int[][] bMatrix= toBitMatrix();
		int start1SS= 0, start2SS= 0, end1SS= 0, end2SS= 0;
		for (int i = 1; i <= bitStart; i++) {	// count skipped ss until first ss
			if (bMatrix[i][0]!= bMatrix[i-1][0])
				++start1SS;
			if (bMatrix[i][1]!= bMatrix[i-1][1])
				++start2SS;
		}
		for (int i = (bitStart+1); i <= bitEnd; i++) {	// SSs taken 
			if (bMatrix[i][0]!= bMatrix[i-1][0])
				++end1SS;
			if (bMatrix[i][1]!= bMatrix[i-1][1])
				++end2SS;
		}
		
		SpliceSite[][] result= new SpliceSite[2][];
		result[0]= new SpliceSite[end1SS];
		for (int i = 0; i < result[0].length; i++) 
			result[0][i]= spliceChain1[start1SS+ i];	// +1 next, -1 0-based
		result[1]= new SpliceSite[end2SS];
		for (int i = 0; i < result[1].length; i++) 
			result[1][i]= spliceChain2[start2SS+ i];
		
		return result;
	}
	
	public int[][] toBitMatrix() {

		int[][] bitMatrix= new int[getDegree()+ 1][];
		for (int i = 0; i < bitMatrix.length; i++) {
			bitMatrix[i]= new int[2];
		}
		SpliceSite ref= null;
		if (spliceChain1!= null&& spliceChain1.length> 0) 
			ref= spliceChain1[0];
		else 
			ref= spliceChain2[0];
		
		boolean exonic1, exonic2;
		if (ref.isDonor()) {
			exonic1= true;
			exonic2= true;
		} else {
			exonic1= false;
			exonic2= false;
		}
		int bit= 0;
		int p1= 0, p2= 0;
		while ((spliceChain1!= null&& p1< spliceChain1.length)||
				(spliceChain2!= null&& p2< spliceChain2.length)) {
			
			bitMatrix[bit][0]= (exonic1)?1:0;
			bitMatrix[bit++][1]= (exonic2)?1:0;

			if (spliceChain1!= null&& p1< spliceChain1.length) {
				if (spliceChain2!= null&& p2< spliceChain2.length) {
					if (spliceChain1[p1].getPos()< spliceChain2[p2].getPos()) {
						exonic1= !exonic1;
						++p1;
					} else {
						exonic2= !exonic2;
						++p2;
					}
				} else {
					exonic1= !exonic1;
					++p1;
				}
			} else {
				exonic2= !exonic2;
				++p2;
			}
		}
		bitMatrix[bit][0]= (exonic1)?1:0;
		bitMatrix[bit++][1]= (exonic2)?1:0;
		
		return bitMatrix;
	}

	/**
	 * String representation with letters.
	 * @deprecated runs out of letters for long splicing variations
	 */
	public String toStringAlpha() {
		String c1= "";
		String c2= "";
		
		char ltt= 65;
		int p1= 0, p2= 0;
		while ((spliceChain1!= null&& p1< spliceChain1.length)||
				(spliceChain2!= null&& p2< spliceChain2.length)) {
			char lttSmall= (char) (ltt+ 32);
			if (spliceChain1!= null&& p1< spliceChain1.length)
				if (spliceChain2!= null&& p2< spliceChain2.length)
					if (spliceChain1[p1].getPos()< spliceChain2[p2].getPos())
						c1+= Character.toString(spliceChain1[p1++].isDonor()?ltt:lttSmall);
					else
						c2+= Character.toString(spliceChain2[p2++].isDonor()?ltt:lttSmall);
				else
					c1+= Character.toString(spliceChain1[p1++].isDonor()?ltt:lttSmall);
			else
				c2+= Character.toString(spliceChain2[p2++].isDonor()?ltt:lttSmall);
				
			ltt++;
		}
				
		return c1+ "\n"+ c2;
	}
	
	
	
	/**
	 * Representation with absolute (chromosomal) coordinates for each splice site. 
	 */
	public String toStringCoordinates() {

		String result= getDegree()+ ":( ";
		for (int j = 0; j < spliceChain1.length; j++) {
			result+= spliceChain1[j].getPos();
			if (spliceChain1[j].isDonor())
				result+= "_";
			else
				result+= "-";
		}
		result+= " / ";
		for (int j = 0; j < spliceChain2.length; j++) {
			result+= spliceChain2[j].getPos();
			if (spliceChain2[j].isDonor())
				result+= "_";
			else
				result+= "-";
		}
		result+= ")";
		
		return result;
	}
	
	public SpliceSite[] getSpliceChain(Transcript transcriptID) {
		if (transcriptID== trans1)
			return spliceChain1;
		if (transcriptID== trans2)
			return spliceChain2;
		return null;
	}
	
	public SpliceSite[] getSpliceUniverse() {
		SpliceSite[] suniv= new SpliceSite[spliceChain1.length+ spliceChain2.length];
		int s1= 0; int s2= 0;
		for (int i = 0; i < suniv.length; i++) {
			if (s2>= spliceChain2.length|| 
					(s1< spliceChain1.length&& spliceChain1[s1].getPos()< spliceChain2[s2].getPos()))
				suniv[i]= spliceChain1[s1++];
			else
				suniv[i]= spliceChain2[s2++];
		}
		return suniv;
	}
	public SpliceSite[] getSpliceChain1() {
		return spliceChain1;
	}
	public SpliceSite[] getSpliceChain2() {
		return spliceChain2;
	}
	public Transcript getTranscript1() {
		return trans1;
	}
	public Transcript getTranscript2() {
		return trans2;
	}
}
