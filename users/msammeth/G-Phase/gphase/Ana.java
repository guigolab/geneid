package gphase;

import java.util.Arrays;
import java.util.Vector;

import gphase.algo.ASAnalyzer;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractSite;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.SpliceSite;

public class Ana {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		testAllSS(ASAnalyzer.getGraph(ASAnalyzer.INPUT_ENCODE));
	}
	
	public static void testAllSS(Graph g) {
		
			// get es
		Gene[] ge= g.getGenes();
		Vector ssV= new Vector();
		Vector ssFV= new Vector();
		for (int i = 0; i < ge.length; i++) {
			SpliceSite[] ss= ge[i].getSpliceSites();
			for (int j = 0; ss!= null&& j < ss.length; j++) 
				if (ss[j].isAcceptor()) {
					ssV.add(ss[j]);
					// ssFV
				}
		}
		
			// get 3'ss
		SpliceSite[] ss3= (SpliceSite[]) gphase.tools.Arrays.toField(ssV);
		SpliceSite[] ssF= (SpliceSite[]) gphase.tools.Arrays.toField(ssFV);
		
			// read seq
		int cntH1= 0, cntH2= 0;
		Vector sV= new Vector();
		for (int i = 0; i < ss3.length; i++) {
			String s= Graph.readSequence(ss3[i], 4, 0);
			if (s.equalsIgnoreCase("acagg")) {
				sV.add(s+" : " /*+ (ss3[i].getPos()- ssF[i].getPos()) + " : "*/+ ss3[i].getGene().getChromosome()+" "+ss3[i].getPos());
				++cntH1;
			} else if (s.equalsIgnoreCase("gcagg")) {
				sV.add(s+" : " /*+ (ss3[i].getPos()- ssF[i].getPos()) + " : "*/+ ss3[i].getGene().getChromosome()+" "+ss3[i].getPos());
				++cntH2;
			}
		}
		Object[] sA= sV.toArray();
		Arrays.sort(sA);
		for (int i = 0; i < sA.length; i++) 
			System.out.println(sA[i]);
		
		System.out.println(cntH1+" + "+cntH2+" / "+ss3.length);
		System.out.println(((float) cntH1/ ss3.length)+ " + "+((float) cntH2/ ss3.length));
	}

	public static void testESss(Graph g) {
		
			// get es
		ASVariation[][] vars= g.getASVariations(ASMultiVariation.FILTER_STRUCTURALLY);
		ASVariation[] esVars= null;
		for (int i = 0; i < vars.length; i++) {
			if (vars[i][0].toString().equals("1-2^ , 0")) {
				esVars= vars[i];
				break;
			} else
				;// System.out.println(vars[i][0]);
		}
		
			// get 3'ss
		SpliceSite[] ss3= new SpliceSite[esVars.length];
		int[] fs= new int[esVars.length];
		for (int i = 0; i < ss3.length; i++) {
			ss3[i]= esVars[i].getSpliceUniverse()[0];
			fs[i]= esVars[i].getFlankingPos()[0];
		}
		
			// read seq
		int cntH1= 0, cntH2= 0;
		Vector sV= new Vector();
		for (int i = 0; i < ss3.length; i++) {
			String s= Graph.readSequence(ss3[i], 4, 0);
			if (s.equalsIgnoreCase("acagg")) {
				sV.add(s+" : "+(ss3[i].getPos()-fs[i])+ " : " + ss3[i].getGene().getChromosome()+" "+ss3[i].getPos());
				++cntH1;
			} else if (s.equalsIgnoreCase("gcagg")) {
				sV.add(s+" : "+(ss3[i].getPos()-fs[i])+ " : " + ss3[i].getGene().getChromosome()+" "+ss3[i].getPos());
				++cntH2;
			}
		}
		Object[] sA= sV.toArray();
		Arrays.sort(sA);
		for (int i = 0; i < sA.length; i++) 
			System.out.println(sA[i]);
		
		System.out.println(cntH1+" + "+cntH2+" / "+ss3.length);
		System.out.println(((float) cntH1/ ss3.length)+ " + "+((float) cntH2/ ss3.length));
	}

}
