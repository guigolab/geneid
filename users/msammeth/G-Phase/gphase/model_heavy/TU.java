package gphase.model_heavy;

import java.util.Comparator;

import gphase.tools.Arrays;

public class TU extends DirectedRegion {

	Gene gene;
	Transcript[] trans= null;
	SpliceSite[] ss= null;
	
	public TU(Gene newGene, Transcript newTranscript) {
		this.gene= newGene;
		setStrand(gene.getStrand());
		addTranscript(newTranscript);
	}
	
	public boolean contains(SpliceSite s) {
		Comparator compi= new SpliceSite.PositionComparator();
		if (ss== null)
			return false;
		int p= java.util.Arrays.binarySearch(ss, s, compi);
		return (p>= 0);
	}
	
	public void addTranscript(Transcript newTranscript) {
		
			// update boundaries of CDS
		if (newTranscript.isCoding()) {
			Translation t= newTranscript.getTranslations()[0];
			if (Math.abs(t.getStart())< Math.abs(getStart()))
				setStart(t.getStart());
			if (Math.abs(t.getEnd())> Math.abs(getEnd()))
				setEnd(t.getEnd());
		}
		
			// add SS
		if (newTranscript.isCoding()) {
			SpliceSite[] s= newTranscript.getSpliceChain();
			Comparator compi= new SpliceSite.PositionComparator();
			for (int i = 0; i < s.length; i++) {
				if (ss== null) 
					ss= new SpliceSite[] {s[i]};
				else {
					int pos= java.util.Arrays.binarySearch(ss, s[i], compi);
					if (pos< 0)
						ss= (SpliceSite[]) Arrays.insert(ss, s[i], pos);
				}
			}
		}
		
		newTranscript.addTU(this);
		trans= (Transcript[]) Arrays.add(trans, newTranscript);
	}

	public boolean removeTranscript(Transcript ntrans) {
		
		ntrans.removeTU(this);
		
			// remove 
		Transcript[] newTranscripts= new Transcript[trans.length- 1];
		int pos= 0;
		boolean flag= false;
		for (int i = 0; i < trans.length; i++) 
			if (trans[i]!= ntrans)
				newTranscripts[pos++]= trans[i];
			else
				flag= true;
		
		if (flag) 
			trans= newTranscripts;

			// remove SSs
		if (ntrans.isCoding()) {
			SpliceSite[] s= ntrans.getSpliceChain();
			for (int i = 0; i < s.length; i++) {
				int j;
				for (j = 0; j < trans.length; j++) {
					if (trans[j]== ntrans)
						continue;
					int k;
					for (k = 0; k < s[i].getTranscripts().length; k++) {
						if (trans[j]== s[i].getTranscripts()[k]&& trans[j].isCoding())
							break;
					}
					if (k< s[i].getTranscripts().length)
						break;
				}
				if (j>= trans.length)
					ss= (SpliceSite[]) Arrays.remove(ss, s[i]);
			}
		}
		
			// update region
		int minStart= Integer.MAX_VALUE, maxEnd= 0;
		for (int i = 0; i < trans.length; i++) {
			if (trans[i].isCoding()) {
				Translation t= trans[i].getTranslations()[0];
				if (Math.abs(t.getStart())< Math.abs(minStart))
					minStart= t.getStart();
				if (Math.abs(t.getEnd())> Math.abs(maxEnd))
					maxEnd= t.getEnd();
			}
		}
		setStart(minStart);
		setEnd(maxEnd);
		
		return flag;
	}

	public Transcript[] getTranscripts() {
		return trans;
	}
	
	
}
