/*
 * Created on Feb 23, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.intransparent.init;

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

	byte ssRegionID3UTR= 0;
	byte ssRegionIDCDS= 0;
	byte ssRegionID5UTR= 0;
	Transcript[][] trpts;
	SpliceSite[][] spliceChains;	// sorted!
	String stringRep= null;

	public ASEvent(Transcript[][] newTrpts, SpliceSite[][] ss) {
		SpliceSite[] firsts= new SpliceSite[ss.length];
		HashMap<SpliceSite, Transcript[]> map= new HashMap<SpliceSite, Transcript[]>();
		HashMap<SpliceSite, SpliceSite[]> map2= new HashMap<SpliceSite, SpliceSite[]>();
		for (int i = 0; i < ss.length; i++) {
			if (ss[i].length> 0) 
				firsts[i]= ss[i][0];
		  	else
				firsts[i]= new SpliceSite(Integer.MAX_VALUE,SpliceSite.TYPE_NOT_INITED);	// empty at end
			map.put(firsts[i], newTrpts[i]);
			map2.put(firsts[i], ss[i]);
		}
		Arrays.sort(firsts, SpliceSite.getDefaultPositionTypeComparator());
		trpts= new Transcript[newTrpts.length][];
		spliceChains= new SpliceSite[ss.length][];
		for (int i = 0; i < firsts.length; i++) {
			trpts[i]= map.get(firsts[i]);
			spliceChains[i]= map2.get(firsts[i]);
		}
	}
	
	public boolean isASevent() {
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
					sb[nextI.get(i)].append(cnt++);
					sb[nextI.get(i)].append(spliceChains[nextI.get(i)][p[nextI.get(i)]++].getSiteSymbol());
				}
				
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

	/**
	 * One line schematical representation of the splicing variation.
	 */
	public String toStringASTA() {
	
			// build final string
		String evCode= toString();
		StringBuffer result= new StringBuffer(evCode+ "\t"+ trpts[0][0].getChromosome());
		
		for (int i = 0; i < spliceChains.length; i++) {
			result.append("\t");
			for (int j = 0; j < trpts[i].length; j++) 
				result.append(trpts[i][j].getTranscriptID()+ ",");
			result.deleteCharAt(result.length()- 1);
			result.append("\t");
			if (spliceChains[i].length> 0) {
				if (spliceChains[i][0].getPos()>= 0)
					for (int x = 0; x < spliceChains[i].length; x++) 
						result.append(spliceChains[i][x].getPos()+",");
				else
					for (int x = spliceChains[i].length- 1; x >= 0; --x)
					//for (int x = 0; x < spliceChains[i].length; x++)	// also forward, for downward compatibility and ASTAcompator
						result.append(Math.abs(spliceChains[i][x].getPos())+",");
				result.deleteCharAt(result.length()- 1);
			}
		}
		
		return result.toString();
	}
	
	
}
