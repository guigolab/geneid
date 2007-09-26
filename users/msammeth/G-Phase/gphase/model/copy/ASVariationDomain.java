package gphase.model.copy;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

public class ASVariationDomain extends ASVariation implements Comparable {
	
	DirectedRegion[] dom1= null, dom2= null;
	HashMap<String,String> domSymbols= null;
	String domCode= null;
	int[] dom1Starts, dom1Ends, dom2Starts, dom2Ends;
	static StructureComparator compi= null;
	
	public static StructureComparator getDefaultComparator() {
		if (compi == null) {
			compi = new StructureComparator();
		}

		return compi;
	}
	
	public int getDomainDegree() {
		return dom1.length+ dom2.length;
	}
	
	public static String stripOrderSuffix(String domName) {
		if (domName== null|| domName.indexOf('-')< 0)
			return domName;
		int p= domName.lastIndexOf('-');
		try {
			Integer.parseInt(domName.substring(p+1));	// only numbers 
		} catch (Exception e) {
			return domName;
		}
		return domName.substring(0, p);
	}
	public ASVariationDomain(ASVariation var, Vector<DirectedRegion> regs1, Vector<DirectedRegion> regs2) {
		ssRegionID3UTR= var.ssRegionID3UTR;
		ssRegionIDCDS= var.ssRegionIDCDS;
		ssRegionID5UTR= var.ssRegionID5UTR;
		trans1= var.trans1;
		trans2= var.trans2;
		spliceChain1= var.spliceChain1;	// sorted!
		spliceChain2= var.spliceChain2;	// sorted!
		degree= var.degree;
		asEvents= var.asEvents; // events
		stringRep= var.stringRep;
		if (var.attributes!= null)
			attributes= (HashMap) var.attributes.clone();
		
		overlap(regs1, regs2);
		toStringDomainCode();
	}
	
	public int compareTo(Object o) {
		return getDefaultComparator().compare(this, o);
	}
	
	public static class StructureComparator extends ASVariation.StructureComparator {
		
		@Override
		public int compare(Object arg0, Object arg1) {
				// do that anyway to init ovl chains
			ASVariationDomain as1= (ASVariationDomain) arg0;
			ASVariationDomain as2= (ASVariationDomain) arg1;
			
			int[] domStarts11= as1.getDom1Starts();
			int[] domEnds11= as1.getDom1Ends();
			int[] domStarts12= as1.getDom2Starts();
			int[] domEnds12= as1.getDom2Ends();
			
			int[] domStarts21= as2.getDom1Starts();
			int[] domEnds21= as2.getDom1Ends();
			int[] domStarts22= as2.getDom2Starts();
			int[] domEnds22= as2.getDom2Ends();
	
			// schain match between ASs nuked, should be handled in Constructor
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;	// structurally different
			
				// the ones with many domains first
			if (as1.getDomainDegree()< as2.getDomainDegree())
				return -1;
			if (as1.getDomainDegree()> as2.getDomainDegree())
				return 1;
						
			
			
				// compare
			if (domStarts11.length< domStarts21.length)
				return -1;
			if (domStarts11.length> domStarts21.length)
				return 1;
			if (domStarts12.length< domStarts22.length)
				return -1;
			if (domStarts12.length> domStarts22.length)
				return 1;
			
			for (int i = 0; i < domStarts11.length; i++) {
				if (domStarts11[i]< domStarts21[i])
					return -1;
				if (domStarts11[i]> domStarts21[i])
					return 1;
			}
			for (int i = 0; i < domEnds11.length; i++) {
				if (domEnds11[i]< domEnds21[i])
					return -1;
				if (domEnds11[i]> domEnds21[i])
					return 1;
			}
			for (int i = 0; i < domStarts12.length; i++) {
				if (domStarts12[i]< domStarts22[i])
					return -1;
				if (domStarts12[i]> domStarts22[i])
					return 1;
			}
			for (int i = 0; i < domEnds12.length; i++) {
				if (domEnds12[i]< domEnds22[i])
					return -1;
				if (domEnds12[i]> domEnds22[i])
					return 1;
			}
			
			// TODO check for domain identity !!!
			
			return 0; 	// equal
		}
	}

	public static class IdentityComparator extends ASVariation.IdentityComparator {
		
		@Override
		public int compare(Object arg0, Object arg1) {
				// do that anyway to init ovl chains
			ASVariationDomain as1= (ASVariationDomain) arg0;
			ASVariationDomain as2= (ASVariationDomain) arg1;
			
			int[] domStarts11= as1.getDom1Starts();
			int[] domEnds11= as1.getDom1Ends();
			int[] domStarts12= as1.getDom2Starts();
			int[] domEnds12= as1.getDom2Ends();
			
			int[] domStarts21= as2.getDom1Starts();
			int[] domEnds21= as2.getDom1Ends();
			int[] domStarts22= as2.getDom2Starts();
			int[] domEnds22= as2.getDom2Ends();
	
			// schain match between ASs nuked, should be handled in Constructor
			
			int res= super.compare(arg0, arg1);
			if (res!= 0)
				return res;	// structurally different
			
				// the ones with many domains first
			if (as1.getDomainDegree()< as2.getDomainDegree())
				return -1;
			if (as1.getDomainDegree()> as2.getDomainDegree())
				return 1;
						
			
			
				// compare
			if (domStarts11.length< domStarts21.length)
				return -1;
			if (domStarts11.length> domStarts21.length)
				return 1;
			if (domStarts12.length< domStarts22.length)
				return -1;
			if (domStarts12.length> domStarts22.length)
				return 1;
			
			for (int i = 0; i < domStarts11.length; i++) {
				if (domStarts11[i]< domStarts21[i])
					return -1;
				if (domStarts11[i]> domStarts21[i])
					return 1;
			}
			for (int i = 0; i < domEnds11.length; i++) {
				if (domEnds11[i]< domEnds21[i])
					return -1;
				if (domEnds11[i]> domEnds21[i])
					return 1;
			}
			for (int i = 0; i < domStarts12.length; i++) {
				if (domStarts12[i]< domStarts22[i])
					return -1;
				if (domStarts12[i]> domStarts22[i])
					return 1;
			}
			for (int i = 0; i < domEnds12.length; i++) {
				if (domEnds12[i]< domEnds22[i])
					return -1;
				if (domEnds12[i]> domEnds22[i])
					return 1;
			}
			
			// TODO check for domain identity !!!
			
			return 0; 	// equal
		}
	}

	/**
	 * input regions must be sorted 
	 * @param regs
	 * @return
	 */
	protected void overlap(Vector<DirectedRegion> regs1, Vector<DirectedRegion> regs2) {
		
			//int[] sp1= SpliceSite.getPositions(spliceChain1); 
			//int[] sp2= SpliceSite.getPositions(spliceChain2); 
			int[] su= SpliceSite.getPositions(getSpliceUniverse());
				
			DirectedRegion reg= getRegion();
			Vector<Integer> domStarts11= new Vector<Integer>(), domEnds11= new Vector<Integer>(),
				domStarts12= new Vector<Integer>(), domEnds12= new Vector<Integer>();
			Vector<DirectedRegion> dom1V= new Vector<DirectedRegion>();
			Vector<DirectedRegion> dom2V= new Vector<DirectedRegion>();
			for (int i = 0; regs1!= null&& i < regs1.size(); i++) {
				if (!regs1.elementAt(i).overlaps(reg))
					continue;
				domStarts11.add(new Integer(gphase.tools.Arrays.convertInsertionPoint(
						Arrays.binarySearch(su, regs1.elementAt(i).get5PrimeEdge()))));
				domEnds11.add(new Integer(gphase.tools.Arrays.convertInsertionPoint(
						Arrays.binarySearch(su, regs1.elementAt(i).get3PrimeEdge()))));
				dom1V.add(regs1.elementAt(i));
			}
			for (int i = 0; regs2!= null&& i < regs2.size(); i++) {
				if (!regs2.elementAt(i).overlaps(reg))
					continue;				
				domStarts12.add(new Integer(gphase.tools.Arrays.convertInsertionPoint(
						Arrays.binarySearch(su, regs2.elementAt(i).get5PrimeEdge()))));
				domEnds12.add(new Integer(gphase.tools.Arrays.convertInsertionPoint(
						Arrays.binarySearch(su, regs2.elementAt(i).get3PrimeEdge()))));
				dom2V.add(regs2.elementAt(i));
			}
	
			if (domStarts11!= null)
				dom1Starts= new int[domStarts11.size()];
			else
				dom1Starts= new int[0];
			for (int i = 0; domStarts11!= null&& i < domStarts11.size(); i++) 
				dom1Starts[i]= domStarts11.elementAt(i).intValue();
			if (domEnds11!= null)
				dom1Ends= new int[domStarts11.size()];
			else
				dom1Ends= new int[0];
			for (int i = 0; domEnds11!= null&& i < domEnds11.size(); i++) 
				dom1Ends[i]= domEnds11.elementAt(i).intValue();
			if (domStarts12!= null)
				dom2Starts= new int[domStarts12.size()];
			else
				dom2Starts= new int[0];
			for (int i = 0; domStarts12!= null&& i < domStarts12.size(); i++) 
				dom2Starts[i]= domStarts12.elementAt(i).intValue();
			if (domEnds12!= null)
				dom2Ends= new int[domEnds12.size()];
			else
				dom2Ends= new int[0];
			for (int i = 0; domEnds12!= null&& i < domEnds12.size(); i++) 
				dom2Ends[i]= domEnds12.elementAt(i).intValue();
			
			dom1= (DirectedRegion[]) gphase.tools.Arrays.toField(dom1V);
			if (dom1== null)
				dom1= new DirectedRegion[0];
			dom2= (DirectedRegion[]) gphase.tools.Arrays.toField(dom2V);
			if (dom2== null)
				dom2= new DirectedRegion[0];

			Vector v= new Vector();
			v.add(dom1Ends); v.add(dom1);
			gphase.tools.Arrays.synchroneousSort(dom1Starts, v);
			orderNeighbors(dom1Starts, dom1Ends, dom1);
			v= new Vector();
			v.add(dom2Ends); v.add(dom2);
			gphase.tools.Arrays.synchroneousSort(dom2Starts, v);
			orderNeighbors(dom2Starts, dom2Ends, dom2);
	}
	
	// pfusch
	private void orderNeighbors(int[] start, int[] end, DirectedRegion[] reg) {
		for (int i = 0; i < start.length- 1; i++) {	// pfusch
			if (start[i]== start[i+1]&& 
					end[i]> end[i+1]) {
				int h= start[i];
				start[i]= start[i+1];
				start[i+1]= h;
				h= end[i];
				end[i]= end[i+1];
				end[i+1]= h;
				DirectedRegion r= reg[i];
				reg[i]= reg[i+1];
				reg[i+1]= r;
			}
		}

	}
	
	/**
	 * domain name with or with-out final numbering
	 * @param domName
	 * @return
	 */
	public String getDomainSymbol(String domName) {
		if (domSymbols == null) {	// build up hash
			domSymbols = new HashMap<String,String>();
			char c= 'A';
			for (int i = 0; i < dom1.length; i++) {
				String id= dom1[i].getID();
				String baseID= stripOrderSuffix(id);
				String symbol= domSymbols.remove(baseID);
				if (symbol== null) 
					symbol= new Character(c++).toString();
				domSymbols.put(id, symbol);
				domSymbols.put(baseID, symbol);
			}
			for (int i = 0; i < dom2.length; i++) {
				String id= dom2[i].getID();
				String baseID= stripOrderSuffix(id);
				String symbol= domSymbols.remove(baseID);
				if (symbol== null) 
					symbol= new Character(c++).toString();
				domSymbols.put(id, symbol);
				domSymbols.put(baseID, symbol);
			}
			
		}

		return domSymbols.get(domName);
	}
	
	private String toStringDomainCode_addDomain(String baseStr, DirectedRegion[] domains, int[] domStarts, int[] domEnds, int scPos) {
		
		if (domStarts== null|| domStarts.length== 0)
			return baseStr;

		int i= 0, j= 0;
		while(true) {
			String ins1= "", ins2= "";
			for (; i < domStarts.length; i++) {
				if (domStarts[i]== scPos) {
					ins1= getDomainSymbol(domains[i++].getID())+"[";
					break;	// has to close 
				}
			}
			for (; j < domEnds.length; j++) {
				if (domEnds[j]== scPos) {
					ins2= getDomainSymbol(domains[j++].getID())+"]";
					break;	// has to open another one 
				}
			}
			
			if (ins1.length()> 0&& ins2.length()> 0&& ins1.charAt(0)< ins2.charAt(0))
				baseStr+= ins1+ ins2;	
			else
				baseStr+= ins2+ ins1;
				
			
				
			if (i== domStarts.length&& j== domEnds.length)
				break;
		}
		return baseStr;
	}
	
	public String toStringDomainCode() {
		if (domCode== null) {
			String c1= "";
			String c2= "";
			
			int ltt= 1;
			int p1= 0, p2= 0;
			while ((spliceChain1!= null&& p1< spliceChain1.length)||
					(spliceChain2!= null&& p2< spliceChain2.length)) {
				if (spliceChain1!= null&& p1< spliceChain1.length) {
					if (spliceChain2!= null&& p2< spliceChain2.length) {
						if (spliceChain1[p1].getPos()< spliceChain2[p2].getPos()) {
							c1= toStringDomainCode_addDomain(c1, dom1, dom1Starts, dom1Ends, getSplicePos1(getSpliceUniverse()[p1]));
							c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
						} else {
							c2= toStringDomainCode_addDomain(c2, dom2, dom2Starts, dom2Ends, getSplicePos1(getSpliceUniverse()[p2]));
							c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
						}
					} else {
						c1= toStringDomainCode_addDomain(c1, dom1, dom1Starts, dom1Ends, getSplicePos1(getSpliceUniverse()[p1]));
						c1+= spliceChain1[p1++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
					}
				} else {
					c2= toStringDomainCode_addDomain(c2, dom2, dom2Starts, dom2Ends, getSplicePos1(getSpliceUniverse()[p2]));
					c2+= spliceChain2[p2++].isDonor()?Integer.toString(ltt)+ "^":Integer.toString(ltt)+ "-";
				}
					
				ltt++;
			}
	
			if (c2.equals("")) {
				if (dom2.length> 0) 
					c2= getDomainSymbol(dom2[0].getID())+"[0"+
					getDomainSymbol(dom2[0].getID())+"]";
				else
					c2= "0";
			} 
			
			if (dom1.length> 0) {
				if (dom1Starts[dom1Starts.length- 1]== spliceChain1.length) 
					c1+= getDomainSymbol(dom1[dom1.length- 1].getID())+ "[";
				if (dom1Ends[dom1Ends.length- 1]== spliceChain1.length)
					c1+= getDomainSymbol(dom1[dom1.length- 1].getID())+ "]";
			}
				
			if (spliceChain2.length> 0&& dom2.length> 0) {
				if (dom2Starts[dom2Starts.length- 1]== spliceChain2.length)
					c2+= getDomainSymbol(dom2[dom2.length- 1].getID())+ "[";
				if (dom2Ends[dom2Ends.length- 1]== spliceChain2.length)
					c2+= getDomainSymbol(dom2[dom2.length- 1].getID())+ "]";
			}
			
			domCode= c1+ " , "+ c2;
		}
		return domCode;
		
	}

	public DirectedRegion[] getDom1() {
		return dom1;
	}

	public DirectedRegion[] getDom2() {
		return dom2;
	}
	
	public DirectedRegion[] getDomains() {
		DirectedRegion[] reg= new DirectedRegion[dom1.length+ dom2.length];
		for (int i = 0; i < dom1.length; i++) 
			reg[i]= dom1[i];
		for (int i = 0; i < dom2.length; i++) 
			reg[i+dom1.length]= dom2[i];
		return reg;
	}
	
	public String[] getDomainIDs() {
		String[] reg= new String[dom1.length+ dom2.length];
		for (int i = 0; i < dom1.length; i++) 
			reg[i]= dom1[i].getID();
		for (int i = 0; i < dom2.length; i++) 
			reg[i+dom1.length]= dom2[i].getID();
		return reg;
	}
	public int[] getDom1Ends() {
		return dom1Ends;
	}
	public int[] getDom1Starts() {
		return dom1Starts;
	}
	public int[] getDom2Ends() {
		return dom2Ends;
	}
	public int[] getDom2Starts() {
		return dom2Starts;
	}

}
