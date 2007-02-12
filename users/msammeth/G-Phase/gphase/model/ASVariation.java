/*
 * Created on Feb 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.tools.ENCODE;

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

	public static final String ID_PURE_AD= "(1^ // 2^)";
	public static final String ID_PURE_AA= "(1= // 2=)";
	
	static final long serialVersionUID = 2433674838499118768L;
	
	public static final int TYPE_ALL= 0;
	public static final int TYPE_CDS= 1;
	public static final int TYPE_UTR= 2;
	public static final int TYPE_5UTR= 3;
	public static final int TYPE_3UTR= 4;
	
	public static final String CODE_INTRON_RETENTION= "1^2- , 0";
	
	byte ssRegionID3UTR= 0;
	byte ssRegionIDCDS= 0;
	byte ssRegionID5UTR= 0;
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

	/**
	 * 
	 * @param reference
	 * @param intersect
	 * @return ASVariation[2][], common vars, unique vars in reference
	 */
	public static ASVariation[][] commonVariations(ASVariation[] reference, ASVariation[] intersect) {
		
		if (reference== null|| intersect== null) {
			ASVariation[][] res= new ASVariation[2][];
			res[0]= new ASVariation[0];
			if (reference== null) 
				res[1]= new ASVariation[0];
			else 
				res[1]= reference;
			return res;	
		}
		
		Comparator compi= new StructureComparator();
		Vector comV= new Vector();
		Vector difV= new Vector();
		for (int i = 0; i < reference.length; i++) {
			int j;
			for (j = 0; j < intersect.length; j++) 
				if (compi.compare(reference[i], intersect[j])== 0)
					break;
			
			if (j< intersect.length)
				comV.add(reference[i]);
			else
				difV.add(reference[i]);
		}
		
		ASVariation[][] res= new ASVariation[2][];
		res[0]= (ASVariation[]) gphase.tools.Arrays.toField(comV);
		res[1]= (ASVariation[]) gphase.tools.Arrays.toField(difV);
		return res;
	}
	
	public int getLengthDiff(boolean exon) {
		int[] a= getLength(exon);
		int diffA= Math.abs(a[0]- a[1]);
		return diffA;
	}
	
	public boolean isIntronRetention() {
		if (toString().equals(CODE_INTRON_RETENTION))
			return true;
		return false;
	}
	
	public DirectedRegion getRegion() {
		SpliceSite[] su= getSpliceUniverse();
		int min= su[0].getPos();
		int max= su[su.length- 1].getPos();
		if (Math.abs(min)> Math.abs(max)) {
			int h= min;
			min= max;
			max= h;
		}
		
		DirectedRegion reg= new DirectedRegion(min, max, trans1.getStrand());
		reg.setChromosome(trans1.getChromosome());
		return reg;
	}
	
	public int getLengthDiff_relative12(boolean exon) {
		int[] a= getLength(exon);
		int diffA= a[0]- a[1];
		return diffA;
	}
	public int getCommonLength() {
		
		int value= 0;
		SpliceSite[] borders= getBorderingSpliceSites();
		SpliceSite[][] schains2= new SpliceSite[][] {spliceChain1, spliceChain2};
		
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite[] flanks= getFlankingSpliceSites();
		Vector ssV= new Vector();
		if (flanks[0]!= null)
			ssV.add(flanks[0]);
		for (int i = 0; i < su.length; i++) 
			ssV.add(su[i]);
		if (flanks[1]!= null)
			ssV.add(flanks[1]);
	
		int commonLength= 0;
		for (int i = 0; i < ssV.size(); i++) {
			SpliceSite ss= (SpliceSite) ssV.elementAt(i);
			if (ss.isDonor())
				continue;
			
				// isolate first cluster
			Transcript tx= null;
			if (trans1.containsSS(ss))
				tx= trans1;
			else
				tx= trans2;
			int startCluster= ss.getPos();
			int endCluster= 0;
			for (int j = i+1; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss2.isDonor()&& tx.contains(ss2.getPos())) {
					endCluster= ss2.getPos();
					break;
				}
			}
			if (endCluster== 0)
				continue;
			
				// isolate second cluster
			int start2Cluster= 0, end2Cluster= 0;
			Transcript ty= (tx== trans1?trans2:trans1);
			for (int j = 0; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss.getPos()< startCluster)
					continue;
				if (ss2.isAcceptor()&& ty.containsSS(ss2)) 
					start2Cluster= ss2.getPos();
				if (start2Cluster!= 0&& ss2.isDonor()&& ty.containsSS(ss2)) {
					end2Cluster= ss2.getPos();
					break;
				}
			}
			if (start2Cluster== 0|| end2Cluster== 0)
				continue;
			
				// check overlap
			if (start2Cluster>= startCluster&& start2Cluster<= endCluster) {
				commonLength+= Math.min(endCluster, end2Cluster)- Math.max(startCluster, start2Cluster)+ 1;
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< Math.max(endCluster, end2Cluster)&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
			} else {
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< endCluster&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
			}
		}
		
		return commonLength;
	}

	public int getDiffLength() {
		
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite[] flanks= getFlankingSites();
		Vector ssV= new Vector();
		ssV.add(flanks[0]);
		for (int i = 0; i < su.length; i++) 
			ssV.add(su[i]);
		ssV.add(flanks[1]);
	
		int length= 0;
		for (int i = 0; i < ssV.size(); i++) {
			
			SpliceSite ss= (SpliceSite) ssV.elementAt(i);
			if (ss.isDonor())
				continue;
			
				// isolate first cluster
			Transcript tx= null;
			if (ss== ssV.elementAt(0)|| ss== ssV.elementAt(ssV.size()- 1)|| trans1.containsSS(ss))
				tx= trans1;
			else
				tx= trans2; 
			int startCluster= ss.getPos();
			int endCluster= 0;
			for (int j = i+1; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss2.isDonor()&& (tx.containsSS(ss2)|| ss2== flanks[1])) {
					endCluster= ss2.getPos();
					break;
				}
			}
			if (endCluster== 0)
				continue;
			
				// isolate second cluster
			int start2Cluster= 0, end2Cluster= 0;
			Transcript ty= (tx== trans1?trans2:trans1);
			for (int j = i; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss.getPos()< startCluster)
					continue;
				if (ss2.isAcceptor()&& (ty.containsSS(ss2)|| ss2== flanks[0])) 
					start2Cluster= ss2.getPos();
				if (start2Cluster!= 0&& ss2.isDonor()&& (ty.containsSS(ss2)|| ss2== flanks[1])) {
					end2Cluster= ss2.getPos();
					break;
				}
			}
			if ((start2Cluster== 0|| end2Cluster== 0)||  
					start2Cluster< startCluster|| start2Cluster> endCluster){	// no 2nd cluster or no overlap
				length+= endCluster- startCluster;
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< endCluster&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
					
			} else {
				int common= Math.min(endCluster, end2Cluster)- Math.max(startCluster, start2Cluster)+ 1;
				int diff= Math.max(endCluster, end2Cluster)- Math.min(startCluster, start2Cluster)+ 1- common;
				length+= diff;
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< Math.max(endCluster, end2Cluster)&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
			}
		}
		
		return length;
	}

	public int getDiffLength_old() {
		
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite[] flanks= getFlankingSpliceSites();
		Vector ssV= new Vector();
		if (flanks[0]!= null)
			ssV.add(flanks[0]);
		for (int i = 0; i < su.length; i++) 
			ssV.add(su[i]);
		if (flanks[1]!= null)
			ssV.add(flanks[1]);

		int length= 0;
		for (int i = 0; i < ssV.size(); i++) {
			SpliceSite ss= (SpliceSite) ssV.elementAt(i);
			if (ss.isDonor())
				continue;
			
				// isolate first cluster
			Transcript tx= null;
			if (trans1.containsSS(ss))
				tx= trans1;
			else
				tx= trans2;
			int startCluster= ss.getPos();
			int endCluster= 0;
			for (int j = i+1; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss2.isDonor()&& tx.containsSS(ss2)) {
					endCluster= ss2.getPos();
					break;
				}
			}
			if (endCluster== 0)
				continue;
			
				// isolate second cluster
			int start2Cluster= 0, end2Cluster= 0;
			Transcript ty= (tx== trans1?trans2:trans1);
			for (int j = i; j < ssV.size(); j++) {
				SpliceSite ss2= (SpliceSite) ssV.elementAt(j);
				if (ss.getPos()< startCluster)
					continue;
				if (ss2.isAcceptor()&& ty.containsSS(ss2)) 
					start2Cluster= ss2.getPos();
				if (start2Cluster!= 0&& ss2.isDonor()&& ty.containsSS(ss2)) {
					end2Cluster= ss2.getPos();
					break;
				}
			}
			if ((start2Cluster== 0|| end2Cluster== 0)||  
					start2Cluster< startCluster|| start2Cluster> endCluster){	// no 2nd cluster or no overlap
				length+= endCluster- startCluster;
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< endCluster&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
					
			} else {
				int common= Math.min(endCluster, end2Cluster)- Math.max(startCluster, start2Cluster)+ 1;
				int diff= Math.max(endCluster, end2Cluster)- Math.min(startCluster, start2Cluster)+ 1- common;
				length+= diff;
				SpliceSite s= (SpliceSite) ssV.elementAt(i);
				while(s.getPos()< Math.max(endCluster, end2Cluster)&& i< ssV.size()- 1)
					s= (SpliceSite) ssV.elementAt(++i);
			}
		}
		
		return length;
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

	/*
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
	public int trim(boolean end5) {
		
		int[] pos1= new int[trans1.getSpliceChain().length+ 2];
		pos1[0]= trans1.getTSSPos();
		for (int i = 1; i < pos1.length- 1; i++) 
			pos1[i]= trans1.getSpliceChain()[i-1].getPos();
		pos1[pos1.length- 1]= trans1.getTESPos();
		int[] pos2= new int[trans2.getSpliceChain().length+ 2];
		pos2[0]= trans2.getTSSPos();
		for (int i = 1; i < pos2.length- 1; i++) 
			pos2[i]= trans2.getSpliceChain()[i-1].getPos();
		pos2[pos2.length- 1]= trans2.getTESPos();
		
		if (end5) {
			int min= Integer.MAX_VALUE;
			for (int i = 0; i < pos1.length; i+= 2) {
				int x= Arrays.binarySearch(pos2, pos1[i]);
				if (x> 0|| (++x% 2!= 0)) {	// pos: don or acc pos ok, neg: insert before acc (after a donor)	
					min= pos1[i];
					break;
				}
			}
			for (int i = 0; i < pos2.length; i+= 2) {
				int x= Arrays.binarySearch(pos1, pos2[i]);
				if (x> 0|| (++x% 2!= 0)) {	// pos: don or acc pos ok, neg: insert before acc (after a donor)	
					min= Math.min(min, pos2[i]);
					break;
				}
			}
			return min;
		} else {
			int max= Integer.MIN_VALUE;
			for (int i = pos1.length- 1; i >= 0; i-= 2) {
				int x= Arrays.binarySearch(pos2, pos1[i]);
				if (x> 0|| (++x% 2!= 0)) {	// pos: don or acc pos ok, neg: insert before acc (after a donor)	
					max= pos1[i];
					break;
				}
			}
			for (int i = pos2.length- 1; i >=0; i-= 2) {
				int x= Arrays.binarySearch(pos1, pos2[i]);
				if (x> 0|| (++x% 2!= 0)) {	// pos: don or acc pos ok, neg: insert before acc (after a donor)	
					max= Math.max(max, pos2[i]);
					break;
				}
			}
			return max;
		}
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

	public static class StructureTSSComparator implements Comparator {
			
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
					((compi.compare(s1[1], s2[0])== 0)&& (compi.compare(s1[0], s2[1])== 0))) {
					int a1= as1.trans1.get5UTR(false);
					int a2= as1.trans1.get5UTR(false);
					int b1= as1.trans2.get5UTR(false);
					int b2= as1.trans2.get5UTR(false);
					
					return 0;
				}
				return -1;
				
			}
		}

	public static class HierarchyFilter extends StructureComparator {
			
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

				if (as1.isCodingFunc()|| as1.is5UTRFunc()|| as1.is3UTRFunc()) 
					return (as2.isProteinCoding())?0:1;	// keep only first
				if (as2.isCodingFunc()|| as2.is5UTRFunc()|| as2.is3UTRFunc())
					return 2;		// keep second
				
				if (as1.isProteinCoding()|| as1.isCompletelyIn5UTR()|| as1.isCompletelyIn3UTR()) 
					return (as2.isProteinCoding())?0:1;	// keep only first
				if (as2.isProteinCoding()|| as2.isCompletelyIn5UTR()|| as2.isCompletelyIn3UTR())
					return 2;		// keep second

					// separate from spooky twilight zone
				if (as1.is_twilight_CDS()|| as1.is_twilight_5UTR()|| as1.is_twilight_3UTR()) 
					return (as2.isProteinCoding())?0:1;	// keep only first
				if (as2.is_twilight_CDS()|| as2.is_twilight_5UTR()|| as2.is_twilight_3UTR())
					return 2;		// keep second
				return 0;	// keep one (whichever)
			}
		}

	/**
	 * sorts according to 
	 * 		shortest schain first
	 * 		lowest ss pos first
	 * @author micha
	 *
	 */
	public static class SpliceSiteOrderComparator implements Comparator {
			
			public int compare(Object arg0, Object arg1) {
				
					// exclude null pointers
				if (arg0== null&& arg1== null)
					return 0;	// cannot happen
				else {
					if (arg0== null)
						return -1;
					else if (arg1== null)
						return 1;
				}
			
				SpliceSite[] ss0= (SpliceSite[]) arg0;
				SpliceSite[] ss1= (SpliceSite[]) arg1;

				if (ss0.length< ss1.length)
					return -1;
				else if (ss1.length< ss0.length)
					return 1;
				
					// 2 schains of equal length
				if (ss0[0].getPos()< ss1[0].getPos())
					return -1;
				return 1;	// else, cannot be equal
	
			}
		}

	public static class RedundancyHierarchyFilter extends StructureComparator {
			
			public int compare(Object arg0, Object arg1) {
				
				int eq= super.compare(arg0, arg1);
				if (eq!= 0)
					return -1;	// structurally different
				
				ASVariation as1= (ASVariation) arg0;
				ASVariation as2= (ASVariation) arg1;
				if (as1.is_non_classifiable()) {
					return (as2.is_non_classifiable())?0:2;	// keep only secd
				} else {
					return (as2.is_non_classifiable())?1:-1;	// keep only secd
				}
			}
		}

	public static class HierarchyMergeFlags extends StructureComparator {
		
		public int compare(Object arg0, Object arg1) {
			
			int eq= super.compare(arg0, arg1);
			if (eq!= 0)
				return -1;	// structurally different
			
			ASVariation as1= (ASVariation) arg0;
			ASVariation as2= (ASVariation) arg1;
			
			if (as2.is_affecting_5UTR())
				as1.setSsRegionID5UTR((byte) 1);
			if (as2.is_affecting_CDS())
				as1.setSsRegionID5UTR((byte) 1);
			if (as2.is_affecting_3UTR())
				as1.setSsRegionID5UTR((byte) 1);
			
			return 2;
		}
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
			//schain[i].setConstitutive(false);
			schain[i].addASVar(this);
		}
	}
	
	public void removeFromASS() {
		SpliceSite[] sc1= trans1.getSpliceChain();
		for (int i = 0; sc1!= null&& i < sc1.length; i++) 
			sc1[i].removeASVar(this);
		SpliceSite[] sc2= trans2.getSpliceChain();
		for (int i = 0; sc2!= null&& i < sc2.length; i++) 
			sc2[i].removeASVar(this);
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
	 * Returns the respective region
	 * @param reg -1= 5'UTR, 0= CDS, 1= 3'UTR
	 * @return the genomic start end end positions of the respective regions
	 * 
	 */
	public AbstractRegion[] getASRegion(int reg) {
		if ((trans1.isNonCoding()&& trans2.isNonCoding())||
				reg< -1|| reg> 1)
			return null;
		
		
		if (!trans1.isNonCoding()) {
			trans1.getTranslations();
		}
		return null;
	}
	
	public boolean isNotProteinCoding() {
		return !(isProteinCoding());
	}

	public boolean isNotProteinCoding_1cover() {
		return !(isProteinCoding_1cover());
	}
	
	
	// no coding transcr or contradicting locations of 2 CDSs 
	public boolean isTwilightZone() {
		if (trans1.getTranslations()== null&& trans2.getTranslations()== null)
			return true;
		if (isNotProteinCoding()&& (!isCompletelyIn5UTR())&& (!isCompletelyIn3UTR()))
			return true;	// improbable to happen though ?!
		return false;
	}
	
	public boolean isUTRnotLastIntron() {
		if (!isCompletelyIn5UTR())
			return false;
			// exclude that affects last intron
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite[] lintron1= trans1.getLastUTRIntron();
		SpliceSite[] lintron2= trans2.getLastUTRIntron();
		if (lintron1== null|| lintron2== null)
			return false;
		for (int i = su.length- 1; i >=0; --i) 
			if ((su[i].getPos()== lintron1[0].getPos()|| su[i].getPos()== lintron1[1].getPos())
					|| (su[i].getPos()== lintron2[0].getPos()|| su[i].getPos()== lintron2[1].getPos()))
				return false;
		return true;
	}
	
	public boolean is_twilight_CDS() {
		if (!isTwilightZone())
			return false;
		
		SpliceSite[] su= getSpliceUniverse();
			// sharing SS with CDS
//		for (int i = 0; trans1.getTu()!= null&& i < trans1.getTu().length; i++) {
//			for (int j = 0; j < su.length; j++) {
//				if (trans1.getTu()[i].contains(su[j]))
//					return true;
//			}
//		}
//		for (int i = 0; trans2.getTu()!= null&& i < trans2.getTu().length; i++) {
//			for (int j = 0; j < su.length; j++) {
//				if (trans2.getTu()[i].contains(su[j]))
//					return true;
//			}
//		}
		
		boolean b= false; 
		for (int i = 0; trans1.getTu()!= null&& i < trans1.getTu().length; i++) {
			if (trans1.getTu()[i].contains(su[0].getPos())&& trans1.getTu()[i].contains(su[su.length- 1].getPos())) {
				if (!b)
					b= true;
			} else {
				if (b)
					return false;	// contained in one CDS but not in the other one
			}
		}		
		for (int i = 0; trans2.getTu()!= null&& i < trans2.getTu().length; i++) {
			if (trans2.getTu()[i].contains(su[0].getPos())&& trans2.getTu()[i].contains(su[su.length- 1].getPos())) {
				if (!b)
					b= true;
			} else {
				if (b)
					return false;	// contained in one CDS but not in the other one
			}
		}
		return b;
	}
	
	public boolean is_twilight_5UTR() {
		
		if (!isTwilightZone())
			return false;
		
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite rightS= su[su.length- 1];
		boolean b= false;
		for (int i = 0; trans1.getTu()!= null&& i < trans1.getTu().length; i++) {
			if (rightS.getPos()< trans1.getTu()[i].get5PrimeEdge()) {
				if (!b) b= true;
			} else {
				if (b) return false;
			}
		}		
		for (int i = 0; trans2.getTu()!= null&& i < trans2.getTu().length; i++) {
			if (rightS.getPos()< trans2.getTu()[i].get5PrimeEdge()) {
					if (!b) b= true;
			} else {
				if (b) return false;
			}
		}		
		return b;
	}
	
	public boolean is_twilight_3UTR() {
		if (!isTwilightZone())
			return false;
		
		SpliceSite[] su= getSpliceUniverse();
		SpliceSite leftS= su[0];
		boolean b= false;
		for (int i = 0; trans1.getTu()!= null&& i < trans1.getTu().length; i++) {
			if (leftS.getPos()> trans1.getTu()[i].get3PrimeEdge()) {
				if (!b) b= true;
			} else {
				if (b) return false;
			}
		}		
		for (int i = 0; trans2.getTu()!= null&& i < trans2.getTu().length; i++) {
			if (leftS.getPos()> trans2.getTu()[i].get3PrimeEdge()) {
					if (!b) b= true;
			} else {
				if (b) return false;
			}
		}		
		return b;
	}
	
	public boolean isTwilightSpooky() {
		if (isTwilightZone()&& !is_twilight_CDS()&& !is_twilight_5UTR()&& !is_twilight_3UTR())
			return true;
		return false;
	}
	
	/**
	 * A splice variation is defined to be protein coding if the 
	 * splicing influences the CDS in one transcript, ie, at least one
	 * ss is contained in a CDS
	 */
	public boolean isProteinCoding() {
	
		SpliceSite[] su1= getSpliceChain1();
		SpliceSite[] su2= getSpliceChain2();
		
		if (su1!= null&& su1.length> 0) {
			for (int i = 0; trans1.getTranslations()!= null&& i < su1.length; i++) {
				if (trans1.getTranslations()[0].contains(su1[i].getPos()))
					return true;
			}
			for (int i = 0; trans2.getTranslations()!= null&& i < su1.length; i++) {
				if (trans2.getTranslations()[0].contains(su1[i].getPos()))
					return true;
			}
		} 
		
		if (su2!= null&& su2.length> 0) {
			for (int i = 0; i < su2.length&& trans2.getTranslations()!= null; i++) {
				if (trans2.getTranslations()[0].contains(su2[i].getPos()))
					return true;
			}
			for (int i = 0; i < su2.length&& trans1.getTranslations()!= null; i++) {
				if (trans1.getTranslations()[0].contains(su2[i].getPos()))
					return true;
			}
		} 
		
		return false;
	}

	/**
	 * A splice variation is defined to be protein coding if the 
	 * splicing influences the CDS in one transcript, ie, at least one
	 * ss is contained in a CDS
	 */
	public boolean isProteinCoding_1cover() {
	
		SpliceSite[] su1= getSpliceChain1();
		SpliceSite[] su2= getSpliceChain2();
		
		if (su1!= null&& su1.length> 0) {
			for (int i = 0; trans1.getTranslations()!= null&& i < su1.length; i++) {
				if ((trans1.getTranslations()[0].contains(su1[i].getPos()))&&
						(trans1.getTranslations()[0].contains(su1[su1.length- 1].getPos())))
					return true;
			}
			for (int i = 0; trans2.getTranslations()!= null&& i < su1.length; i++) {
				if ((trans2.getTranslations()[0].contains(su1[i].getPos()))&& 
						(trans2.getTranslations()[0].contains(su1[su1.length- 1].getPos())))
					return true;
			}
		} 
		
		if (su2!= null&& su2.length> 0) {
			for (int i = 0; i < su2.length&& trans2.getTranslations()!= null; i++) {
				if ((trans2.getTranslations()[0].contains(su2[i].getPos()))&&
						(trans2.getTranslations()[0].contains(su2[su2.length- 1].getPos())))
					return true;
			}
			for (int i = 0; i < su2.length&& trans1.getTranslations()!= null; i++) {
				if ((trans1.getTranslations()[0].contains(su2[i].getPos()))&&
						(trans1.getTranslations()[0].contains(su2[su2.length- 1].getPos())))
					return true;
			}
		} 
		
		return false;
	}

	/**
			 * A splice variation is defined to be protein coding if the 
			 * constitutive flanking splice sites lie within any translation
			 * of the respective transcript.
			 */
			public boolean isProteinCoding_old_publ() {
			
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

	public int[] getFlankingPos() {
				
				SpliceSite[] borders= getBorderingSpliceSites();
				SpliceSite pred= trans1.getPredSpliceSite(borders[0].getPos());			
				SpliceSite succ= trans1.getSuccSpliceSite(borders[1].getPos());
				SpliceSite pred2= trans2.getPredSpliceSite(borders[0].getPos());
				SpliceSite succ2= trans2.getSuccSpliceSite(borders[1].getPos());
			
				int start;
				if (pred== null|| pred2== null|| pred!= pred2) {
					start= trim(true);
				} else {
					start= pred.getPos();
				}
					
				int end;
				if (succ== null|| succ2== null|| succ!= succ2) {
					end= trim(false);
				} else {
					end= succ.getPos();
				}
			
				return new int[] {start, end};
			}

	public SpliceSite[] getFlankingSites() {
		
		SpliceSite[] borders= getBorderingSpliceSites();
		SpliceSite pred= trans1.getPredSpliceSite(borders[0].getPos());			
		SpliceSite succ= trans1.getSuccSpliceSite(borders[1].getPos());
		SpliceSite pred2= trans2.getPredSpliceSite(borders[0].getPos());
		SpliceSite succ2= trans2.getSuccSpliceSite(borders[1].getPos());

		SpliceSite[] result= new SpliceSite[2];
		if (pred== null|| pred2== null|| pred!= pred2) {
			result[0]= new SpliceSite(null, trim(true), false);	// trimmed tss
		} else {
			result[0]= pred;
		}
			
		if (succ== null|| succ2== null|| succ!= succ2) {
			result[1]= new SpliceSite(null, trim(false), true);
		} else {
			result[1]= succ;
		}

		return result;
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
	
	public boolean is5UTRFunc() {
		return (isCompletelyIn5UTR()&& trans1.getTranslations()!= null&& trans2.getTranslations()!= null);
	}
	
	public boolean is5UTRMaxTranscriptSS() {
		if (ssRegionID5UTR== 0) {
			SpliceSite[] su= getSpliceUniverse();
			int i;
			for (i = 0; i < su.length; i++) 
				if (su[i].is5UTRMaxTranscript()) {
					ssRegionID5UTR= 1;
					break;
				}
			if (i== su.length)
				ssRegionID5UTR= -1;
		}
		return (ssRegionID5UTR> 0);
	}
	
	public boolean is3UTRMaxTranscriptSS() {
		if (ssRegionID3UTR== 0) {
			SpliceSite[] su= getSpliceUniverse();
			int i;
			for (i = 0; i < su.length; i++) 
				if (su[i].is3UTRMaxTranscript()) {
					ssRegionID3UTR= 1;
					break;
				}
			if (i== su.length)
				ssRegionID3UTR= -1;
		}
		return (ssRegionID3UTR> 0);
	}

	public boolean is_all() {
		return true;
	}
	
	public DirectedRegion[] getAlternativeIntrons() {
		AbstractSite[] flanks= getFlankingSites();
		Vector regV= new Vector();		
		for (int i = 0; i < getSpliceChain1().length; i++) {
			if (getSpliceChain1()[i].isDonor())
				continue;
			int start;
			if (i== 0)
				start= flanks[0].getPos();
			else
				start= getSpliceChain1()[i-1].getPos();
			
			++start;
			int end= getSpliceChain1()[i].getPos()- 1;
			DirectedRegion reg= new DirectedRegion(start, end, getTranscript1().getStrand());
			reg.setSpecies(getTranscript1().getSpecies());
			reg.setChromosome(getTranscript1().getChromosome());
			regV.add(reg);
		}
		for (int i = 0; i < getSpliceChain2().length; i++) {
			if (getSpliceChain2()[i].isDonor())
				continue;
			int start;
			if (i== 0)
				start= flanks[0].getPos();
			else
				start= getSpliceChain2()[i-1].getPos();
				
			++start;
			int end= getSpliceChain2()[i].getPos()- 1;
			DirectedRegion reg= new DirectedRegion(start, end, getTranscript2().getStrand());
			reg.setSpecies(getTranscript1().getSpecies());
			reg.setChromosome(getTranscript1().getChromosome());
			regV.add(reg);
		}
		if ((getSpliceChain1().length> 0&& getSpliceChain1()[getSpliceChain1().length-1].isDonor())
				|| (getSpliceChain2().length> 0&& getSpliceChain2()[getSpliceChain2().length-1].isDonor())) {
			int start;
			if (getSpliceChain1().length< 1)
				start= flanks[0].getPos();
			else
				start= getSpliceChain1()[getSpliceChain1().length- 1].getPos();

			++start;
			int end= flanks[1].getPos()- 1;
			DirectedRegion reg= new DirectedRegion(start, end,getTranscript1().getStrand());
			reg.setSpecies(getTranscript1().getSpecies());
			reg.setChromosome(getTranscript1().getChromosome());
			regV.add(reg);
			
			if (getSpliceChain2().length< 1)
				start= flanks[0].getPos();
			else
				start= getSpliceChain2()[getSpliceChain2().length- 1].getPos();

			++start;
			end= flanks[1].getPos()- 1;
			reg= new DirectedRegion(start, end, getTranscript2().getStrand());
			reg.setSpecies(getTranscript1().getSpecies());
			reg.setChromosome(getTranscript1().getChromosome());
			regV.add(reg);
		}
			
		return (DirectedRegion[]) gphase.tools.Arrays.toField(regV);
	}
	
	public boolean has_nonGTAG_Intron() {
		DirectedRegion[] regs= getAlternativeIntrons();
		if (regs== null)
			return false;
		for (int i = 0; i < regs.length; i++) {
			String seq= Graph.readSequence(regs[i]);
			if (seq.length()< 4)// GTAG
				return true;
			String don= seq.substring(0, 2);
			String acc= seq.substring(seq.length()- 2, seq.length());
			if ((!don.equalsIgnoreCase("GT"))|| (!acc.equalsIgnoreCase("AG")))
				return true;
		}
		return false;		
	}
	
	public boolean isAlignmentArtefact() {
		return false;
	}
	
	public boolean is_affecting_3UTR() {
		if (ssRegionID3UTR== 0) {
			//SpliceSite[] su= getSpliceUniversePlusFlanks();
			SpliceSite[] su= getSpliceUniverse();
			int i;
			for (i = 0; i < su.length; i++) 
				if ((trans1.getTranslations()!= null&& trans1.getTranslations()[0].get3PrimeEdge()< su[i].getPos())||
						(trans2.getTranslations()!= null&& trans2.getTranslations()[0].get3PrimeEdge()< su[i].getPos())) {
					ssRegionID3UTR= 1;
					break;
				}
			if (i== su.length)
				ssRegionID3UTR= -1;
		}
		return (ssRegionID3UTR> 0);
	}

	public boolean is_affecting_5UTR() {
		if (ssRegionID5UTR== 0) {
			SpliceSite[] su= getSpliceUniverse();	//getSpliceUniversePlusFlanks();
			int i;
			for (i = 0; i < su.length; i++) 
				if ((trans1.getTranslations()!= null&& trans1.getTranslations()[0].get5PrimeEdge()> su[i].getPos())||
						(trans2.getTranslations()!= null&& trans2.getTranslations()[0].get5PrimeEdge()> su[i].getPos())) {
					ssRegionID5UTR= 1;
					break;
				}
			if (i== su.length)
				ssRegionID5UTR= -1;
		}
		return (ssRegionID5UTR> 0);
	}
	
	public boolean isCodingFunc() {
		return (isProteinCoding()&& trans1.getTranslations()!= null&& trans2.getTranslations()!= null);
	}
	
	public boolean isCodingMaxTranscriptSS() {
		if (ssRegionIDCDS== 0) {
			SpliceSite[] su= getSpliceUniverse();
			int i;
			for (i = 0; i < su.length; i++) 
				if (su[i].isCDSMaxTranscript()) {
					ssRegionIDCDS= 1;
					break;
				}
			if (i== su.length)
				ssRegionIDCDS= -1;
		}
		return (ssRegionIDCDS> 0);
	}

	public boolean is_affecting_CDS() {
		if (ssRegionIDCDS== 0) {
			SpliceSite[] su= getSpliceUniverse();
			int i;
			for (i = 0; i < su.length; i++) 
				if ((trans1.getTranslations()!= null&& trans1.getTranslations()[0].contains(su[i].getPos()))||
						(trans2.getTranslations()!= null&& trans2.getTranslations()[0].contains(su[i].getPos()))) {
					ssRegionIDCDS= 1;
					break;
				}
			if (i== su.length)
				ssRegionIDCDS= -1;
		}
		return (ssRegionIDCDS> 0);
	}
	
	public boolean is_contained_in_CDS() {
		if (is_affecting_CDS()&& !is_affecting_5UTR()&& !is_affecting_3UTR())
			return true;
		return false;
	}
	
	public boolean is_contained_in_5UTR() {
		if (!is_affecting_CDS()&& is_affecting_5UTR()&& !is_affecting_3UTR())
			return true;
		return false;
	}
	
	public boolean is_contained_in_3UTR() {
		if (!is_affecting_CDS()&& !is_affecting_5UTR()&& is_affecting_3UTR())
			return true;
		return false;
	}
	
	public boolean isTwilightMaxTranscriptSS() {
		int ctr= 0;
		if (is5UTRMaxTranscriptSS())
			++ctr;
		if (isCodingMaxTranscriptSS())
			++ctr;
		if (is3UTRMaxTranscriptSS())
			++ctr;
		return (ctr> 1);
	}

	public boolean isTwilight() {
		int ctr= 0;
		if (is_affecting_5UTR())
			++ctr;
		if (is_affecting_CDS())
			++ctr;
		if (is_affecting_3UTR())
			++ctr;
		return (ctr> 1);
	}
	
	public boolean is_twilight_5UTR_CDS() {
		if (is_affecting_5UTR()&& is_affecting_CDS())
			return true;
		return false;
	}
	
	public boolean is_twilight_CDS_3UTR() {
		if (is_affecting_3UTR()&& is_affecting_CDS())
			return true;
		return false;
	}
	
	public boolean is_twilight_5UTR_3UTR() {
		if (is_affecting_3UTR()&& is_affecting_5UTR())
			return true;
		return false;
	}

	public boolean is_twilight_5UTR_CDS_3UTR() {
		if (is_affecting_3UTR()&& is_affecting_5UTR()&& is_affecting_CDS())
			return true;
		return false;
	}
	
	public boolean isNoneMaxTranscriptSS() {
		int ctr= 0;
		if (is5UTRMaxTranscriptSS())
			++ctr;
		if (isCodingMaxTranscriptSS())
			++ctr;
		if (is3UTRMaxTranscriptSS())
			++ctr;
		return (ctr== 0);
	}

	public boolean is_non_classifiable() {
		int ctr= 0;
		if (is_affecting_5UTR())
			++ctr;
		if (is_affecting_CDS())
			++ctr;
		if (is_affecting_3UTR())
			++ctr;
		return (ctr== 0);
	}
	
	public boolean isTwilightFunc() {
		return (isTwilightZone()&& trans1.getTranslations()!= null&& trans2.getTranslations()!= null);
	}
	
	public boolean is3UTRFunc() {
		return (isCompletelyIn3UTR()&& trans1.getTranslations()!= null&& trans2.getTranslations()!= null);
	}
	
	public boolean isCompletelyIn5UTR() {
	
			SpliceSite[] su= getSpliceUniverse();
			boolean t1= (trans1.is5UTR(su[0].getPos())&& trans1.is5UTR(su[su.length- 1].getPos()));
			boolean t2= (trans2.is5UTR(su[0].getPos())&& trans2.is5UTR(su[su.length- 1].getPos()));
			
			if ((t1&& t2)|| (t1&& trans2.getTranslations()== null)|| (t2&& trans1.getTranslations()== null))
				return true;
			return false;
		}

	public boolean isCompletelyIn3UTR() {
	
			SpliceSite[] su= getSpliceUniverse();
			boolean t1= trans1.is3UTR(su[0].getPos())&& trans1.is3UTR(su[su.length- 1].getPos());
			boolean t2= trans2.is3UTR(su[0].getPos())&& trans2.is3UTR(su[su.length- 1].getPos());
			
			if ((t1&& t2)|| (t1&& trans2.getTranslations()== null)|| (t2&& trans1.getTranslations()== null))
				return true;
			return false;
		}

	public boolean isCompletelyInCDS() {
	
			SpliceSite[] su= getSpliceUniverse();
			boolean t1= trans1.isCDS(su[0].getPos())&& trans1.isCDS(su[su.length- 1].getPos());
			boolean t2= trans2.isCDS(su[0].getPos())&& trans2.isCDS(su[su.length- 1].getPos());
			
			if ((t1&& t2)|| t1&& (trans2.getTranslations()== null)|| t2&& trans1.getTranslations()== null)
				return true;
			return false;
		}

	public boolean isTouching5UTR() {
	
		SpliceSite[] su= getSpliceUniverse();
		if (trans1.is5UTR(su[0].getPos())|| trans2.is5UTR(su[0].getPos())|| trans1.is5UTR(su[su.length- 1].getPos())|| trans2.is5UTR(su[su.length- 1].getPos()))
			return true;
		return false;
	}

	public boolean isTouching3UTR() {

		SpliceSite[] su= getSpliceUniverse();
		if (trans1.is3UTR(su[0].getPos())|| trans2.is3UTR(su[0].getPos())|| trans1.is3UTR(su[su.length- 1].getPos())|| trans2.is3UTR(su[su.length- 1].getPos()))
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
	
	public String toCoordinates() {
		String s= trans1.getTranscriptID()+ ":";
		for (int i = 0; spliceChain1!= null&& i < spliceChain1.length; i++) 
			s+= spliceChain1[i]+" ";
		s+= "\t"+ trans2.getTranscriptID();
		for (int i = 0; spliceChain2!= null&& i < spliceChain2.length; i++) 
			s+= spliceChain2[i]+" ";
		return s;
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
						c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
					else
						c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
				else
					c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
			else
				c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
				
			ltt++;
		}

//		if (c2.length()< c1.length()) {	// print shortest first
//			String h= c1;
//			c1= c2;
//			c2= h;
//		} else if (c1.length()== c2.length()) {	// sort according to first number
//			int x= 0;
//			while (Character.isDigit(c1.charAt(x++)));
//			int n1= Integer.parseInt(c1.substring(0, x-1));
//			x= 0;
//			while (Character.isDigit(c2.charAt(x++)));
//			int n2= Integer.parseInt(c2.substring(0, x-1));
//			if (n2< n1)  {
//				String h= c1;
//				c1= c2;
//				c2= h;
//			}
//		}
		if (c1.equals("")|| c2.equals("")) {
			if (c1.equals("")) {
				c1= "0";
				String h= c1;
				c1= c2;
				c2= h;
			} else
				c2= "0";
		} else {
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
		return c1+ " , "+ c2;
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
	
	/**
	 * 
	 * @param bitStart
	 * @param bitEnd
	 * @return null if there is no flanking site (TSS/TES)
	 */
	public DirectedRegion getBitRegion(int bitPos) {

		int[][] bMatrix= toBitMatrix();

			// start
		int start= 0;
		int start1SS= -1, start2SS= -1, end1SS= 0, end2SS= 0;
		if (bitPos< 1) {	// get last SS, max(TSS)
			SpliceSite pred= getFlankingSpliceSites()[0];
			if (pred== null)
				start= Math.max(trans1.get5PrimeEdge(), trans2.get5PrimeEdge());
			else
				start= pred.getPos();
		} else {			// find pred SS in schain
			for (int i = 1; i <= bitPos; i++) {	// count skipped ss until first ss
				if (bMatrix[i][0]!= bMatrix[i-1][0])
					++start1SS;
				if (bMatrix[i][1]!= bMatrix[i-1][1])
					++start2SS;
			}
			if (start1SS>= 0) {
				if (start2SS>= 0) 
					start= Math.max(spliceChain1[start1SS].getPos(), spliceChain2[start2SS].getPos());
				else 
					start= spliceChain1[start1SS].getPos();
			} else {
				start= spliceChain2[start2SS].getPos();
			}
		}
		
			// end
		int end= 0;
		if (bitPos+ 1>= bMatrix.length) {
			SpliceSite succ= getFlankingSpliceSites()[1];
			if (succ== null)
				end= Math.min(trans1.get3PrimeEdge(), trans2.get3PrimeEdge());
			else
				end= succ.getPos();
		} else {
			
			if (bMatrix[bitPos][0]!= bMatrix[bitPos+ 1][0]) {
				end= spliceChain1[start1SS+ 1].getPos();
				if (bMatrix[bitPos][1]!= bMatrix[bitPos+ 1][1]) 
					System.err.println("assertion failed, 2 cannot change in 1 bitmap col");
			} else {
				if (bMatrix[bitPos][1]!= bMatrix[bitPos+ 1][1]) 
					end= spliceChain2[start2SS+ 1].getPos();
				else
					System.err.println("assertion failed, 1 has to change in a bitmap col");
			}
		}
		
		if (Math.abs(start)> Math.abs(end)) {
			int h= start;
			start= end;
			end= h;
		}
		DirectedRegion dir= new DirectedRegion(start, end, getGene().getStrand());
		dir.setChromosome(getGene().getChromosome());
		return dir;
	}

	/**
	 * 
	 * @param bitStart
	 * @param bitEnd
	 * @return null if there is no flanking site (TSS/TES)
	 */
	public SpliceSite[][] getSpliceSites(int bitStart, int bitEnd) {
		int[][] bMatrix= toBitMatrix();
		int start1SS= 0, start2SS= 0, end1SS= 0, end2SS= 0;
		for (int i = 1; i <= bitStart; i++) {	// count skipped ss until first ss
			if (bMatrix[i][0]!= bMatrix[i-1][0])
				++start1SS;
			if (bMatrix[i][1]!= bMatrix[i-1][1])
				++start2SS;
		}
		boolean lowerFlank= false, upperFlank= false;
		for (int i = (bitStart+1); i <= bitEnd; i++) {	// SSs taken 
			if (bMatrix[i][0]!= bMatrix[i-1][0])
				++end1SS;
			if (bMatrix[i][1]!= bMatrix[i-1][1])
				++end2SS;
		}
		
		SpliceSite[][] result= new SpliceSite[2][];
		result[0]= new SpliceSite[end1SS];
		for (int i= 0; i < result[0].length; i++) 
			result[0][i]= spliceChain1[start1SS+ i];	// +1 next, -1 0-based
		result[1]= new SpliceSite[end2SS];
		for (int i= 0; i < result[1].length; i++) 
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
	
	
	public DirectedRegion getVariationArea() {
		SpliceSite[] su= getSpliceUniverse();
		DirectedRegion reg= new DirectedRegion(su[0].getPos(), su[su.length- 1].getPos(), trans1.getStrand());
		reg.setChromosome(trans1.getChromosome());
		return reg;
	}
	
	public DirectedRegion[] getVariableRegion() {
		return getVariableRegion("BC");
	}
	
	public DirectedRegion[] getVariableRegion1() {
		return getVariableRegion("B");
	}
	
	public DirectedRegion[] getVariableRegion2() {
		return getVariableRegion("C");
	}
	
	public DirectedRegion[] getVariableRegion(String id) {
		String bit= toBitString();
		Vector regV= new Vector();
		for (int i = 0; i < bit.length(); i++) {
			if (id.contains(Character.toString(bit.charAt(i)))) {
					// create variable area
				SpliceSite[][] ss= getSpliceSites(i-1, i+1);
				int min= 0, max= 0;
				if (ss[0]!= null&& ss[0].length> 0) {
					if (ss[1]!= null&& ss[1].length> 0) {
						min= Math.min(Math.abs(ss[0][0].getPos()), 
								Math.abs(ss[1][0].getPos()));
						max= Math.max(Math.abs(ss[1][ss[1].length- 1].getPos()), 
								Math.abs(ss[0][ss[0].length- 1].getPos()));
					} else {
						min= Math.abs(ss[0][0].getPos());
						max= Math.abs(ss[0][ss[0].length- 1].getPos());
					}
				} else {
					if (ss[1]!= null&& ss[1].length> 0) {
						min= Math.abs(ss[1][0].getPos());
						max= Math.abs(ss[1][ss[1].length- 1].getPos());
					} else {
						System.err.println("empty event");
					}
				}
				if (min> max) {
					int h= min;
					min= max;
					max= h;
				}
				
					// create regions
				DirectedRegion reg= new DirectedRegion(min, max, getGene().getStrand());
				reg.setChromosome(getGene().getChromosome());
				regV.add(reg);
			}
		}
		return (DirectedRegion[]) gphase.tools.Arrays.toField(regV);
	}

	public DirectedRegion[] getExonicOnlyRegion() {
		String bit= toBitString();
		Vector regV= new Vector();
		for (int i = 0; i < bit.length(); i++) {
			if (bit.charAt(i)== 'D') {
					// create variable area
				DirectedRegion dir= getBitRegion(i);
				regV.add(dir);
			}
		}
		return (DirectedRegion[]) gphase.tools.Arrays.toField(regV);
	}

	public int getMinVarSSPos() {
		if (spliceChain1!= null&& spliceChain1.length> 0) {
			if (spliceChain2!= null&& spliceChain2.length> 0)
				return Math.min(spliceChain1[0].getPos(), spliceChain2[0].getPos());
			else
				return spliceChain1[0].getPos();
		} else {
			if (spliceChain2!= null&& spliceChain2.length> 0)
				return spliceChain2[0].getPos();
			else
				System.err.println("empty variation");
		}
		return -1;
	}

	
		// genomic regions, pos
	DefaultRegion[] getExonicRegions(SpliceSite[] sc) {
		int min, max;
		Vector regV= new Vector();
		for (int i = 0; sc!= null&& i < sc.length; i++) {
			if (!sc[i].isDonor())
				continue;
			if (i== 0) 
				min= Math.abs(getMinVarSSPos());
			else 
				min= Math.abs(sc[i-1].getPos());
			max= sc[i].getPos();
			
			if (min== max)	// empty
				continue;
			
			if (min> max) {	// swap on neg strand
				int h= min;
				min= max;
				max= h;
			}
			DefaultRegion reg= new DefaultRegion(min, max);
			reg.setChromosome(getGene().getChromosome());
			regV.add(reg);
		}
		
		return (DefaultRegion[]) gphase.tools.Arrays.toField(regV);
	}
	
	// returns genomic coordinates (always pos)
	public AbstractRegion[][] getVariableRegions() {
		
		AbstractRegion[][] result= new AbstractRegion[2][];
		result[0]= getExonicRegions(spliceChain1);
		result[1]= getExonicRegions(spliceChain2);
		return result;
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

	/**
	 * Representation with absolute (chromosomal) coordinates for each splice site. 
	 */
	public String toStringUCSC() {
		
		String result= "chr"+ trans1.getChromosome()+ ":";
		SpliceSite[] su= getSpliceUniverse();
		int left= Math.abs(su[0].getPos());
		int right= Math.abs(su[su.length- 1].getPos());
		if (left> right) {
			int h= left;
			left= right;
			right= h;
		}
		result+= left+ "-"+ right;
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

	public SpliceSite[] getSpliceUniversePlusFlanks() {
		SpliceSite[] suniv= new SpliceSite[spliceChain1.length+ spliceChain2.length];
		SpliceSite[] flanks= getFlankingSites();	// getFlankingSpliceSites();
		suniv[0]= flanks[0];
		int s1= 0; int s2= 0;
		for (int i = 1; i < suniv.length- 1; i++) {
			if (s2>= spliceChain2.length|| 
					(s1< spliceChain1.length&& spliceChain1[s1].getPos()< spliceChain2[s2].getPos()))
				suniv[i]= spliceChain1[s1++];
			else
				suniv[i]= spliceChain2[s2++];
		}
		suniv[suniv.length- 1]= flanks[1];
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

	public void setSsRegionID3UTR(byte ssRegionID3UTR) {
		this.ssRegionID3UTR = ssRegionID3UTR;
	}

	public void setSsRegionID5UTR(byte ssRegionID5UTR) {
		this.ssRegionID5UTR = ssRegionID5UTR;
	}

	public String getSsRegionIDs() {
		String s= ssRegionID5UTR+" "+ssRegionIDCDS+" "+ssRegionID3UTR;
		return s;
	}

	public void setSsRegionIDCDS(byte ssRegionIDCDS) {
		this.ssRegionIDCDS = ssRegionIDCDS;
	}
}
