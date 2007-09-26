/*
 * Created on Feb 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model;

import gphase.AStaLaVista;
import gphase.tools.ENCODE;
import gphase.tools.IntVector;

import java.io.PrintStream;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.text.Collator;
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
public class ASEvent implements Serializable {

	public static class SpliceChainComparator implements Comparator<SpliceSite[]> {
		public int compare(SpliceSite[] o1, SpliceSite[] o2) {
			if (o1== null)	// empty chains at the end
				return 1;
			if (o2== null)
				return -1;

			int min= Math.min(o1.length, o2.length);
			for (int i = 0; i < min; i++) {
				if (o1[i].getPos()< o2[i].getPos())
					return -1;
				else if (o2[i].getPos()< o1[i].getPos())
					return 1;
			}
			
			if (o1.length> min)
				return 1;
			if (o2.length> min)
				return -1;
				// equal length and equal ss
			return 0;
		}
	}
	
	public static SpliceChainComparator defaultSpliceChainComparator= new SpliceChainComparator();

	
	byte ssRegionID3UTR= 0;
	byte ssRegionIDCDS= 0;
	byte ssRegionID5UTR= 0;
	Transcript[][] trpts;
	SpliceSite[][] spliceChains;	// sorted!
	String stringRep= null;
	SpliceSite src, snk;
	public ASEvent(Transcript[][] newTrpts, SpliceSite[][] ss) {
		SpliceSite[] firsts= new SpliceSite[ss.length];
		HashMap<SpliceSite[], Transcript[]> map= new HashMap<SpliceSite[], Transcript[]>(ss.length, 1f);
		for (int i = 0; i < ss.length; i++) {
			map.put(ss[i], newTrpts[i]);	// null still can only be one, we hope
		}
		Arrays.sort(ss, defaultSpliceChainComparator);
		
		trpts= new Transcript[newTrpts.length][];
		for (int i = 0; i < ss.length; i++) {
			trpts[i]= map.get(ss[i]);
		}
		spliceChains= ss;
	}
	
	/**
	 * 0= not included
	 * 1= overlapping
	 * 2= included
	 * @return
	 */
	public int inCDS() {
		int min= Integer.MAX_VALUE;
		int max= Integer.MIN_VALUE;
		for (int i = 0; i < spliceChains.length; i++) {
			if (spliceChains[i].length== 0)
				continue;
			if (spliceChains[i][0].getPos()< min)
				min= spliceChains[i][0].getPos();
			if (spliceChains[i][spliceChains[i].length-1].getPos()> max)
				max= spliceChains[i][spliceChains[i].length-1].getPos();
		}
		
		DirectedRegion[] cdss= trpts[0][0].getGene().getCDSS();	// assume non-overlapping
		for (int j = 0; j < cdss.length; j++) {
			boolean minContained= cdss[j].contains(min);
			boolean maxContained= cdss[j].contains(max);
			if (minContained&& maxContained)
				return 2;
			if (minContained|| maxContained || (min< cdss[j].get5PrimeEdge()&& max> cdss[j].get3PrimeEdge()))
				return 1;
			
		}
		return 0;
	}
	
	public void setAnchors(SpliceSite src, SpliceSite snk) {
		this.src= src;
		this.snk= snk;
	}
	
	public boolean isASevent() {
		
			// common boundaries
		int min= Integer.MIN_VALUE;
		int max= Integer.MAX_VALUE;
		for (int i = 0; i < trpts.length; i++) {
			// TODO chk? another loop over all transcripts needed? guess not
			min= Math.max(min, trpts[i][0].get5PrimeEdge());
			max= Math.min(max, trpts[i][0].get3PrimeEdge());
		}
		
		for (int i = 0; i < spliceChains.length; i++) {
			for (int j = 0; j < spliceChains[i].length; j++) {
				if (!spliceChains[i][j].isSpliceSite())
					continue;
				if (spliceChains[i][j].getPos()>= min&& spliceChains[i][j].getPos()<= max)
					return true;
			}
		}
		
		return false;
	}

	// checks only first transcript w sSs, see mutex 1st exons
	public boolean isASevent_not_ok() {
		int x;
		for (x = 0; x < spliceChains.length; x++) 
			if (spliceChains[x].length> 0)
				break;
		for (int i = 0; i < spliceChains[x].length; i++) {
			if (!spliceChains[x][i].isSpliceSite())
				continue;
			int j;
			for (j = 0; j < trpts.length; j++) {
				if (j== x)
					continue;
				if (!trpts[j][0].contains(spliceChains[x][i].getPos()))
					break;
			}
			if (j== trpts.length)	// one ss common to all found
				return true;
		}
		
		return false;
	}

	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toString() {
		
		if (stringRep== null) {
			int[] p= new int[spliceChains.length];
			StringBuffer[] sb= new StringBuffer[spliceChains.length];
			for (int i = 0; i < p.length; i++) { 
				p[i]= 0;
				sb[i]= new StringBuffer();
			}
			int cnt= 1;
			
			while (true) {
				IntVector nextI= new IntVector();
				int nextVal= Integer.MAX_VALUE;
				for (int i = 0; i < spliceChains.length; i++) {
					if (p[i]== spliceChains[i].length)
						continue;
					if (spliceChains[i][p[i]].getPos()< nextVal) {						
						nextI= new IntVector();
						nextI.add(i);
						nextVal= spliceChains[i][p[i]].getPos();
					} else if (spliceChains[i][p[i]].getPos()== nextVal)
						nextI.add(i);
				}
				
				for (int i = 0; i < nextI.size(); i++) {
					sb[nextI.get(i)].append(cnt);
					sb[nextI.get(i)].append(spliceChains[nextI.get(i)][p[nextI.get(i)]++].getSiteSymbol());
				}
				++cnt;	// now increment, for identical positions having same index
				
				int x= 0;
				for (; x < p.length; x++) 
					if (p[x]< spliceChains[x].length)
						break;
				if (x== p.length)
					break;
			}

			StringBuffer stringBuf= new StringBuffer();
			for (int i = 0; i < sb.length; i++) {
				if (sb[i].length()> 0)
					stringBuf.append(sb[i]);
				else 
					stringBuf.append('0');
				stringBuf.append(",");
			}
			stringBuf.deleteCharAt(stringBuf.length()- 1);
			stringRep= stringBuf.toString();
		}
		
		return stringRep;
	}

	public String[] getSources() {
		String[] sources= new String[trpts.length];
		for (int i = 0; i < sources.length; i++) {
			byte type= trpts[i][0].getSourceType();
			for (int j = 1; type!= 0&& j < trpts[i].length; j++) 
				if (trpts[i][j].getSourceType()< type)
					type= trpts[i][j].getSourceType();
			sources[i]= Transcript.TYPE_TO_ID[type];
		}
		return sources;
	}
	
	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toStringASTA() {
	
			// build final string
		StringBuffer result= new StringBuffer(Integer.toString(trpts.length));
		String evCode= toString();
		result.append("\t");
		result.append(evCode);
		result.append("\t");
		result.append("CDS=");
		result.append(inCDS());
		result.append("\t");
		result.append(trpts[0][0].getChromosome());
		result.append(":");
		int srcPos= src.getPos();
		if (srcPos== Integer.MIN_VALUE) {
			srcPos= Integer.MAX_VALUE;
			for (int i = 0; i < spliceChains.length; i++) 
				if (spliceChains[i].length> 0&& spliceChains[i][0].getPos()< srcPos)
					srcPos= spliceChains[i][0].getPos();
		}
		srcPos-= 1;
		
		int snkPos= snk.getPos();
		if (snkPos== Integer.MAX_VALUE) {
			snkPos= Integer.MIN_VALUE;
			for (int i = 0; i < spliceChains.length; i++) 
				if (spliceChains[i].length> 0&& spliceChains[i][spliceChains[i].length-1].getPos()> snkPos)
					snkPos= spliceChains[i][spliceChains[i].length-1].getPos();
		}
		snkPos+= 1;
			
		if (srcPos< 0) {
			result.append(Math.abs(snkPos));
			result.append("-");
			result.append(Math.abs(srcPos));
			result.append("\t");
			result.append("-");
		} else {
			result.append(srcPos);
			result.append("-");
			result.append(snkPos);
			result.append("\t");
			result.append("+");
		}
		
		
		for (int i = 0; i < spliceChains.length; i++) {
			result.append("\t");
			if (spliceChains[i].length> 0) {
				for (int x = 0; x < spliceChains[i].length; x++) {	// do not reverse order of SSs <-> breaks syntheny with rel. positions 
					if (spliceChains[i][x].getPos()>= 0)
						result.append(spliceChains[i][x].getPos());
					else
						result.append(Math.abs(spliceChains[i][x].getPos()));
					result.append(",");
				}
				result.deleteCharAt(result.length()- 1);
			}
		}

		for (int i = 0; i < trpts.length; i++) {
			result.append("\t");
			for (int j = 0; j < trpts[i].length; j++) { 
				result.append(trpts[i][j].getTranscriptID());
				result.append(",");
			}
			result.deleteCharAt(result.length()- 1);
		}
		
		String[] sources= getSources();
		for (int i = 0; i < sources.length; i++) {
			result.append("\t");
			result.append(sources[i]);
		}
		
		return result.toString();
	}

	public SpliceSite getSnk() {
		return snk;
	}

	public void setSnk(SpliceSite snk) {
		this.snk = snk;
	}

	public SpliceSite getSrc() {
		return src;
	}

	public void setSrc(SpliceSite src) {
		this.src = src;
	}
	
	
}
