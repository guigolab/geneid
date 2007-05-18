package gphase;

import gphase.algo.AlgoHandler;
import gphase.io.gtf.GTFChrReader;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.SpliceSite;
import gphase.tools.Formatter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.lang.reflect.Method;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Vector;

import sun.management.counter.Variability;
import sun.security.action.GetLongAction;

import com.sun.net.ssl.SSLContext;
import com.sun.org.apache.bcel.internal.generic.GETSTATIC;

public class Noboru {

	public static final String SUBDIR_NOBORU= "noboru";
	public static final String BASE_ORDER= "ACGT";
	public static final String SFX_DONOR_TABLE= "_freq_don";
	public static final String SFX_ACCEPTOR_TABLE= "_freq_acc";
	
	
	String speName= null;
	String annoName= null;
	double[][] donTab= null;
	double donMinT= -1d, donMaxT= -1d;
	double[][] accTab= null;
	double accL1= -1d, accL2= -1d, accH1= -1d, accH2= -1d;
	String inFile= null;
	


	public static void testNoboru() {
		
		String speName= "fruitfly";
		String annoName= "RefSeq";
		
		HashMap localHash= new HashMap();
		try {
			String fName= Constants.getLatestUCSCAnnotation(speName, annoName, null);
			System.out.println("Getting IR splice sites from "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) 
					_03_getIRspliceSites(genes[i], null, localHash);
				reader.read();
				genes= reader.getGenes();
			}			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
	
		Noboru nobi= new Noboru(speName, annoName);
		System.out.println("high_IRF_donors\n");
		double[] scores= nobi.scoreSpliceSite((SpliceSite[]) gphase.tools.Arrays.toField(
				(Vector) localHash.get("high_IRF_donors")));
		for (int i = 0; i < scores.length; i++) 
			System.out.println(scores[i]);
		System.out.println("\n\nhigh_IRF_acceptors\n");
		scores= nobi.scoreSpliceSite((SpliceSite[]) gphase.tools.Arrays.toField(
				(Vector) localHash.get("high_IRF_acceptors")));
		for (int i = 0; i < scores.length; i++) 
			System.out.println(scores[i]);

		

	}
	public static void main(String[] args) {
		try {
			testNoboru();
			if (1== 1)
				System.exit(0);
			
			Class[] par= new Class[] {Gene.class, String.class, HashMap.class};
				//"_02_outputSpliceSites"
			Method m= Noboru.class.getMethod("_03_getIRspliceSiteSequences", par);
			String fName= Constants.getLatestUCSCAnnotation("fruitfly", "RefSeq", null);
			AlgoHandler._00_mainLoopChromosomes(fName, new Method[] {m});
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void _02_outputSpliceSites(Gene g, String fName, HashMap map) {
		final String ID_CONST_DONOR= "const_donor"; 
		final String ID_CONST_ACCEPTOR= "const_acceptor"; 
		final String ID_ALT_DONOR= "alt_donor"; 
		final String ID_ALT_ACCEPTOR= "alt_acceptor";
		final String ID_SS_POS_STRAND= "ss_pos";
		final String ID_SS_NEG_STRAND= "ss_neg";
		final String ID_SS_ALT= "ss_alt";
		final String ID_SS_CONST= "ss_const";
		
		boolean gtagOnly= true;
		boolean constSS= true;
		boolean altSS= false;
		
		Vector constDonorV= (Vector) map.get(ID_CONST_DONOR);
		if (constDonorV== null) {
			constDonorV= new Vector();
			map.put(ID_CONST_DONOR, constDonorV);
		}
		Vector constAcceptorV= (Vector) map.get(ID_CONST_ACCEPTOR);
		if (constAcceptorV== null) {
			constAcceptorV= new Vector();
			map.put(ID_CONST_ACCEPTOR, constAcceptorV);
		}
		Vector altDonorV= (Vector) map.get(ID_ALT_DONOR);
		if (altDonorV== null) {
			altDonorV= new Vector();
			map.put(ID_ALT_DONOR, altDonorV);
		}
		Vector altAcceptorV= (Vector) map.get(ID_ALT_ACCEPTOR);
		if (altAcceptorV== null) {
			altAcceptorV= new Vector();
			map.put(ID_ALT_ACCEPTOR, altAcceptorV);
		}
		Integer valPos= (Integer) map.get(ID_SS_POS_STRAND);
		if (valPos== null)
			valPos= new Integer(0);
		Integer valNeg= (Integer) map.get(ID_SS_NEG_STRAND);
		if (valNeg== null)
			valNeg= new Integer(0);
		Integer valConst= (Integer) map.get(ID_SS_CONST);
		if (valConst== null)
			valConst= new Integer(0);
		Integer valAlt= (Integer) map.get(ID_SS_ALT);
		if (valAlt== null)
			valAlt= new Integer(0);
		
		if (g!= null) {
			g.getASVariations(ASMultiVariation.FILTER_NONE);
			SpliceSite[] ss= g.getSpliceSites();
			int cntPos= valPos.intValue(), cntNeg= valNeg.intValue(), 
				cntConst= valConst.intValue(), cntAlt= valAlt.intValue();
			for (int i = 0; ss!= null&& i < ss.length; i++) {
				if (ss[i].getPos()>= 0) 
					map.put(ID_SS_POS_STRAND, new Integer(++cntPos));
				else
					map.put(ID_SS_NEG_STRAND, new Integer(++cntNeg));
	
				if (ss[i].isConstitutive()) {
					map.put(ID_SS_CONST, new Integer(++cntConst));
					if (!constSS)
						continue;
				} else {
					map.put(ID_SS_ALT, new Integer(++cntAlt));
					if (!altSS)
						continue;
				}
	
				DirectedRegion reg= ss[i].getShapiroRegion();
				String seq= Graph.readSequence(reg);
				
				int don5Extent= 3;
				int acc5Extent= 11;
				if (ss[i].isDonor()) {
					String dinuc= seq.substring(don5Extent, don5Extent+2).toLowerCase();
					if (gtagOnly&& !dinuc.equals("gt"))
						continue;
					seq=
						seq.substring(0, don5Extent).toUpperCase()+
						dinuc+
						seq.substring(don5Extent+2, seq.length()).toUpperCase();
					if (ss[i].isConstitutive())
						constDonorV.add(ss[i]);
					else
						altDonorV.add(ss[i]);
				} else {
					String dinuc= seq.substring(acc5Extent, acc5Extent+ 2).toLowerCase();
					if (gtagOnly&& !dinuc.equals("ag"))
						continue;
					seq=
						seq.substring(0, acc5Extent).toUpperCase()+
						dinuc+
						seq.substring(acc5Extent+2, seq.length()).toUpperCase();
					if (ss[i].isConstitutive())
						constAcceptorV.add(ss[i]);
					else
						altAcceptorV.add(ss[i]);
				}
	
			}
		}
		
			// output
		if (fName!= null) {
			int sum= valPos.intValue()+ valNeg.intValue();
			if (sum!= valConst.intValue()+ valAlt.intValue())
				System.err.println("Error: sums different.");
			System.out.println("Found "+ sum+ " SSs, "
					+ valPos+ " pos ("+Formatter.fprint((valPos* 100d)/ sum, 2)+ "%), "
					+ valNeg+ " neg ("+Formatter.fprint((valNeg* 100d)/ sum, 2)+ "%), "
					+ valConst+ " const ("+Formatter.fprint((valConst* 100d)/ sum, 2)+ "%), "
					+ valAlt+ " alt ("+Formatter.fprint((valAlt* 100d)/ sum, 2)+ "%), "
			);
			try {
				PrintStream p= new PrintStream(fName+"_donors_const-"+constSS+"_alt-"+altSS+"_gtag-"+gtagOnly);
				for (int i = 0; i < constDonorV.size(); i++) {
					p.print(">constDon"+i+"\n");
					p.print(constDonorV.elementAt(i)+"\n");
					p.print("\n");
				}
				for (int i = 0; i < altDonorV.size(); i++) {
					p.print(">altDon"+i+"\n");
					p.print(altDonorV.elementAt(i)+"\n");
					p.print("\n");
				}
				p.flush(); p.close();
				
				p= new PrintStream(fName+"_acceptors_const-"+constSS+"_alt-"+altSS+"_gtag-"+gtagOnly);
				for (int i = 0; i < constAcceptorV.size(); i++) {
					p.print(">constAcc"+i+"\n");
					p.print(constAcceptorV.elementAt(i)+"\n");
					p.print("\n");
				}
				for (int i = 0; i < altAcceptorV.size(); i++) {
					p.print(">altAcc"+i+"\n");
					p.print(altAcceptorV.elementAt(i)+"\n");
					p.print("\n");
				}
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	public static void _02_outputSpliceSiteSequences(Gene g, String fName, HashMap map) {
		final String ID_CONST_DONOR= "const_donor"; 
		final String ID_CONST_ACCEPTOR= "const_acceptor"; 
		final String ID_ALT_DONOR= "alt_donor"; 
		final String ID_ALT_ACCEPTOR= "alt_acceptor";
		final String ID_SS_POS_STRAND= "ss_pos";
		final String ID_SS_NEG_STRAND= "ss_neg";
		final String ID_SS_ALT= "ss_alt";
		final String ID_SS_CONST= "ss_const";
		
		boolean gtagOnly= true;
		boolean constSS= true;
		boolean altSS= false;
		
		Vector constDonorV= (Vector) map.get(ID_CONST_DONOR);
		if (constDonorV== null) {
			constDonorV= new Vector();
			map.put(ID_CONST_DONOR, constDonorV);
		}
		Vector constAcceptorV= (Vector) map.get(ID_CONST_ACCEPTOR);
		if (constAcceptorV== null) {
			constAcceptorV= new Vector();
			map.put(ID_CONST_ACCEPTOR, constAcceptorV);
		}
		Vector altDonorV= (Vector) map.get(ID_ALT_DONOR);
		if (altDonorV== null) {
			altDonorV= new Vector();
			map.put(ID_ALT_DONOR, altDonorV);
		}
		Vector altAcceptorV= (Vector) map.get(ID_ALT_ACCEPTOR);
		if (altAcceptorV== null) {
			altAcceptorV= new Vector();
			map.put(ID_ALT_ACCEPTOR, altAcceptorV);
		}
		Integer valPos= (Integer) map.get(ID_SS_POS_STRAND);
		if (valPos== null)
			valPos= new Integer(0);
		Integer valNeg= (Integer) map.get(ID_SS_NEG_STRAND);
		if (valNeg== null)
			valNeg= new Integer(0);
		Integer valConst= (Integer) map.get(ID_SS_CONST);
		if (valConst== null)
			valConst= new Integer(0);
		Integer valAlt= (Integer) map.get(ID_SS_ALT);
		if (valAlt== null)
			valAlt= new Integer(0);
		
		if (g!= null) {
			g.getASVariations(ASMultiVariation.FILTER_NONE);
			SpliceSite[] ss= g.getSpliceSites();
			int cntPos= valPos.intValue(), cntNeg= valNeg.intValue(), 
				cntConst= valConst.intValue(), cntAlt= valAlt.intValue();
			for (int i = 0; ss!= null&& i < ss.length; i++) {
				if (ss[i].getPos()>= 0) 
					map.put(ID_SS_POS_STRAND, new Integer(++cntPos));
				else
					map.put(ID_SS_NEG_STRAND, new Integer(++cntNeg));

				if (ss[i].isConstitutive()) {
					map.put(ID_SS_CONST, new Integer(++cntConst));
					if (!constSS)
						continue;
				} else {
					map.put(ID_SS_ALT, new Integer(++cntAlt));
					if (!altSS)
						continue;
				}

				DirectedRegion reg= ss[i].getShapiroRegion();
				String seq= Graph.readSequence(reg);
				
				int don5Extent= SpliceSite.NOBORU_DON5_EXTENT;
				int acc5Extent= SpliceSite.NOBORU_ACC5_EXTENT;
				if (ss[i].isDonor()) {
					String dinuc= seq.substring(don5Extent, don5Extent+2).toLowerCase();
					if (gtagOnly&& !dinuc.equals("gt"))
						continue;
					seq=
						seq.substring(0, don5Extent).toUpperCase()+
						dinuc+
						seq.substring(don5Extent+2, seq.length()).toUpperCase();
					if (ss[i].isConstitutive())
						constDonorV.add(seq);
					else
						altDonorV.add(seq);
				} else {
					String dinuc= seq.substring(acc5Extent, acc5Extent+ 2).toLowerCase();
					if (gtagOnly&& !dinuc.equals("ag"))
						continue;
					seq=
						seq.substring(0, acc5Extent).toUpperCase()+
						dinuc+
						seq.substring(acc5Extent+2, seq.length()).toUpperCase();
					if (ss[i].isConstitutive())
						constAcceptorV.add(seq);
					else
						altAcceptorV.add(seq);
				}

			}
		}
		
			// output
		if (fName!= null) {
			int sum= valPos.intValue()+ valNeg.intValue();
			if (sum!= valConst.intValue()+ valAlt.intValue())
				System.err.println("Error: sums different.");
			System.out.println("Found "+ sum+ " SSs, "
					+ valPos+ " pos ("+Formatter.fprint((valPos* 100d)/ sum, 2)+ "%), "
					+ valNeg+ " neg ("+Formatter.fprint((valNeg* 100d)/ sum, 2)+ "%), "
					+ valConst+ " const ("+Formatter.fprint((valConst* 100d)/ sum, 2)+ "%), "
					+ valAlt+ " alt ("+Formatter.fprint((valAlt* 100d)/ sum, 2)+ "%), "
			);
			try {
				PrintStream p= new PrintStream(fName+"_donors_const-"+constSS+"_alt-"+altSS+"_gtag-"+gtagOnly);
				for (int i = 0; i < constDonorV.size(); i++) {
					p.print(">constDon"+i+"\n");
					p.print(constDonorV.elementAt(i)+"\n");
					p.print("\n");
				}
				for (int i = 0; i < altDonorV.size(); i++) {
					p.print(">altDon"+i+"\n");
					p.print(altDonorV.elementAt(i)+"\n");
					p.print("\n");
				}
				p.flush(); p.close();
				
				p= new PrintStream(fName+"_acceptors_const-"+constSS+"_alt-"+altSS+"_gtag-"+gtagOnly);
				for (int i = 0; i < constAcceptorV.size(); i++) {
					p.print(">constAcc"+i+"\n");
					p.print(constAcceptorV.elementAt(i)+"\n");
					p.print("\n");
				}
				for (int i = 0; i < altAcceptorV.size(); i++) {
					p.print(">altAcc"+i+"\n");
					p.print(altAcceptorV.elementAt(i)+"\n");
					p.print("\n");
				}
				p.flush(); p.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	public static void _03_getIRspliceSiteSequences(Gene g, String fName, HashMap map) {
		final String ID_HIGH_IRF_DONORS= "high_IRF_donors";
		final String ID_HIGH_IRF_ACCEPTORS= "high_IRF_acceptors";
		final String ID_LOW_IRF_DONORS= "low_IRF_donors";
		final String ID_LOW_IRF_ACCEPTORS= "low_IRF_acceptors";
		final String ID_HIGH_IRF_INTRONS= "high_IRF_introns";
		final String ID_LOW_IRF_INTRONS= "low_IRF_introns";
		final String ID_NB_IR_EVENTS= "nb_ir_events";
		final String ID_NB_SS_EQ_IRF_FREQ= "nb_ss_eq_irf_freq";
		
		Vector highIRFdonorV= (Vector) map.get(ID_HIGH_IRF_DONORS);
		if (highIRFdonorV== null) {
			highIRFdonorV= new Vector();
			map.put(ID_HIGH_IRF_DONORS, highIRFdonorV);
		}
		Vector highIRFacceptorV= (Vector) map.get(ID_HIGH_IRF_ACCEPTORS);
		if (highIRFacceptorV== null) {
			highIRFacceptorV= new Vector();
			map.put(ID_HIGH_IRF_ACCEPTORS, highIRFacceptorV);
		}
		Vector lowIRFdonorV= (Vector) map.get(ID_LOW_IRF_DONORS);
		if (lowIRFdonorV== null) {
			lowIRFdonorV= new Vector();
			map.put(ID_LOW_IRF_DONORS, lowIRFdonorV);
		}
		Vector lowIRFacceptorV= (Vector) map.get(ID_LOW_IRF_ACCEPTORS);
		if (lowIRFacceptorV== null) {
			lowIRFacceptorV= new Vector();
			map.put(ID_LOW_IRF_ACCEPTORS, lowIRFacceptorV);
		}
		Vector highIRFintronV= (Vector) map.get(ID_HIGH_IRF_INTRONS);
		if (highIRFintronV== null) {
			highIRFintronV= new Vector();
			map.put(ID_HIGH_IRF_INTRONS, highIRFintronV);
		}
		Vector lowIRFintronV= (Vector) map.get(ID_LOW_IRF_INTRONS);
		if (lowIRFintronV== null) {
			lowIRFintronV= new Vector();
			map.put(ID_LOW_IRF_INTRONS, lowIRFintronV);
		}
		Integer nbIRevents= (Integer) map.get(ID_NB_IR_EVENTS);
		if (nbIRevents== null)
			nbIRevents= new Integer(0);
		Integer nbEqIRFfreq= (Integer) map.get(ID_NB_SS_EQ_IRF_FREQ);
		if (nbEqIRFfreq== null)
			nbEqIRFfreq= new Integer(0);
		
		boolean gtagOnly= true;
		if (g!= null) {
			ASVariation[] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			Vector v= new Vector();
			for (int i = 0; vars!= null&& i < vars.length; i++) 
				if (vars[i].isIntronRetention())
					v.add(vars[i]);
			map.put(ID_NB_IR_EVENTS, new Integer(nbIRevents.intValue()+ v.size()));
			int cntEQfreqIR= nbEqIRFfreq.intValue();
			for (int i = 0; i < v.size(); i++) {
				SpliceSite[] ss= ((ASVariation) v.elementAt(i)).getSpliceUniverse();
				int prime5= ss[0].getPos()+ 1;
				int prime3= ss[1].getPos()- 1;
				DirectedRegion regIntron= new DirectedRegion(prime5, prime3, ss[0].getTranscripts()[0].getStrand());
				regIntron.setChromosome(ss[0].getTranscripts()[0].getChromosome());
				regIntron.setSpecies(ss[0].getTranscripts()[0].getSpecies());
	
					// determine IRF
				int cntRetained= 0, cntNotRetained= 0;
				for (int j = 0; j < g.getTranscripts().length; j++) {
					DirectedRegion[] introns= g.getTranscripts()[j].getIntrons();
					for (int k = 0; introns!= null&& k < introns.length; k++) 
						if (introns[k].equals(regIntron)) {
							++cntNotRetained;
							break;
						}
				}
				for (int j = 0; j < g.getTranscripts().length; j++) {
					DirectedRegion[] exons= g.getTranscripts()[j].getExons();
					for (int k = 0; k < exons.length; k++) 
						if (exons[k].contains(regIntron)) {
							++cntRetained;
							break;
						}
				}
				
					// get sequence				
				int don5Extent= SpliceSite.NOBORU_DON5_EXTENT;
				DirectedRegion reg= ss[0].getNoboruRegion();
				String seqDon= Graph.readSequence(reg);
				String dinuc= seqDon.substring(don5Extent, don5Extent+2).toLowerCase();
				if (gtagOnly&& !dinuc.equals("gt"))
					continue;
				seqDon=
					seqDon.substring(0, don5Extent).toUpperCase()+
					dinuc+
					seqDon.substring(don5Extent+2, seqDon.length()).toUpperCase();
				
				int acc5Extent= SpliceSite.NOBORU_ACC5_EXTENT;				
				reg= ss[1].getNoboruRegion();
				String seqAcc= Graph.readSequence(reg);
				dinuc= seqAcc.substring(acc5Extent, acc5Extent+ 2).toLowerCase();
				if (gtagOnly&& !dinuc.equals("ag"))
					continue;
				seqAcc=
					seqAcc.substring(0, acc5Extent).toUpperCase()+
					dinuc+
					seqAcc.substring(acc5Extent+2, seqAcc.length()).toUpperCase();
	
					// test
//				Noboru nob= new Noboru("fruitfly", "RefSeq");
//				double[] scores= nob.scoreSpliceSite(ss);
				
				
					// sort/add
				if (cntNotRetained> cntRetained) {
					lowIRFdonorV.add(seqDon);
					lowIRFacceptorV.add(seqAcc);
					lowIRFintronV.add(regIntron.getLength()+"\t"+
							((ASVariation) v.elementAt(i)).is_affecting_5UTR()+ "\t"+
							((ASVariation) v.elementAt(i)).is_affecting_CDS()+"\t"+
							((ASVariation) v.elementAt(i)).is_affecting_3UTR());
				} else  if (cntRetained> cntNotRetained) {
					highIRFdonorV.add(seqDon);
					highIRFacceptorV.add(seqAcc);
					highIRFintronV.add(regIntron.getLength()+"\t"+
							((ASVariation) v.elementAt(i)).is_affecting_5UTR()+ "\t"+
							((ASVariation) v.elementAt(i)).is_affecting_CDS()+"\t"+
							((ASVariation) v.elementAt(i)).is_affecting_3UTR());
				} else {	// equal
					map.put(ID_NB_SS_EQ_IRF_FREQ, new Integer(++cntEQfreqIR));
				}
			}
		}
		
		
		if (fName!= null) {
			try {
				System.out.println("Found "+nbEqIRFfreq.intValue()+" IR events with equal freqs.");
				
				BufferedWriter buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_donors_gtag-"+gtagOnly));
				for (int i = 0; i < lowIRFdonorV.size(); i++) {
					buffy.write(">Donor|don"+i+"|loIRF|gtag-"+gtagOnly+"\n");
					buffy.write(lowIRFdonorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_acceptors_gtag-"+gtagOnly));
				for (int i = 0; i < lowIRFacceptorV.size(); i++) { 
					buffy.write(">Acceptor|acc"+i+"|loIRF|gtag-"+gtagOnly+"\n");
					buffy.write(lowIRFacceptorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_introns_gtag-"+gtagOnly));
				buffy.write("intronNr\tlength\taff_5UTR\taff_CDS\taff_3UTR\n");
				for (int i = 0; i < lowIRFintronV.size(); i++) { 
					buffy.write("intron"+i+"\t"+lowIRFintronV.elementAt(i)+"\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_donors_gtag-"+gtagOnly));
				for (int i = 0; i < highIRFdonorV.size(); i++) {
					buffy.write(">Donor|don"+i+"|hiIRF|gtag-"+gtagOnly+"\n");
					buffy.write(highIRFdonorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_acceptors_gtag-"+gtagOnly));
				for (int i = 0; i < highIRFacceptorV.size(); i++) { 
					buffy.write(">Acceptor|acc"+i+"|hiIRF|gtag-"+gtagOnly+"\n");
					buffy.write(highIRFacceptorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_introns_gtag-"+gtagOnly));
				buffy.write("intronNr\tlength\taff_5UTR\taff_CDS\taff_3UTR\n");
				for (int i = 0; i < highIRFintronV.size(); i++) { 
					buffy.write("intron"+i+"\t"+highIRFintronV.elementAt(i)+"\n");
				}
				buffy.flush(); buffy.close();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	public static void _03_getIRspliceSites(Gene g, String fName, HashMap map) {
		final String ID_HIGH_IRF_DONORS= "high_IRF_donors";
		final String ID_HIGH_IRF_ACCEPTORS= "high_IRF_acceptors";
		final String ID_LOW_IRF_DONORS= "low_IRF_donors";
		final String ID_LOW_IRF_ACCEPTORS= "low_IRF_acceptors";
		final String ID_HIGH_IRF_INTRONS= "high_IRF_introns";
		final String ID_LOW_IRF_INTRONS= "low_IRF_introns";
		final String ID_NB_IR_EVENTS= "nb_ir_events";
		final String ID_NB_SS_EQ_IRF_FREQ= "nb_ss_eq_irf_freq";
		
		Vector highIRFdonorV= (Vector) map.get(ID_HIGH_IRF_DONORS);
		if (highIRFdonorV== null) {
			highIRFdonorV= new Vector();
			map.put(ID_HIGH_IRF_DONORS, highIRFdonorV);
		}
		Vector highIRFacceptorV= (Vector) map.get(ID_HIGH_IRF_ACCEPTORS);
		if (highIRFacceptorV== null) {
			highIRFacceptorV= new Vector();
			map.put(ID_HIGH_IRF_ACCEPTORS, highIRFacceptorV);
		}
		Vector lowIRFdonorV= (Vector) map.get(ID_LOW_IRF_DONORS);
		if (lowIRFdonorV== null) {
			lowIRFdonorV= new Vector();
			map.put(ID_LOW_IRF_DONORS, lowIRFdonorV);
		}
		Vector lowIRFacceptorV= (Vector) map.get(ID_LOW_IRF_ACCEPTORS);
		if (lowIRFacceptorV== null) {
			lowIRFacceptorV= new Vector();
			map.put(ID_LOW_IRF_ACCEPTORS, lowIRFacceptorV);
		}
		Vector highIRFintronV= (Vector) map.get(ID_HIGH_IRF_INTRONS);
		if (highIRFintronV== null) {
			highIRFintronV= new Vector();
			map.put(ID_HIGH_IRF_INTRONS, highIRFintronV);
		}
		Vector lowIRFintronV= (Vector) map.get(ID_LOW_IRF_INTRONS);
		if (lowIRFintronV== null) {
			lowIRFintronV= new Vector();
			map.put(ID_LOW_IRF_INTRONS, lowIRFintronV);
		}
		Integer nbIRevents= (Integer) map.get(ID_NB_IR_EVENTS);
		if (nbIRevents== null)
			nbIRevents= new Integer(0);
		Integer nbEqIRFfreq= (Integer) map.get(ID_NB_SS_EQ_IRF_FREQ);
		if (nbEqIRFfreq== null)
			nbEqIRFfreq= new Integer(0);
		
		boolean gtagOnly= true;
		if (g!= null) {
			ASVariation[] vars= g.getASVariations(ASMultiVariation.FILTER_HIERARCHICALLY);
			Vector v= new Vector();
			for (int i = 0; vars!= null&& i < vars.length; i++) 
				if (vars[i].isIntronRetention())
					v.add(vars[i]);
			map.put(ID_NB_IR_EVENTS, new Integer(nbIRevents.intValue()+ v.size()));
			int cntEQfreqIR= nbEqIRFfreq.intValue();
			for (int i = 0; i < v.size(); i++) {
				SpliceSite[] ss= ((ASVariation) v.elementAt(i)).getSpliceUniverse();
				int prime5= ss[0].getPos()+ 1;
				int prime3= ss[1].getPos()- 1;
				DirectedRegion regIntron= new DirectedRegion(prime5, prime3, ss[0].getTranscripts()[0].getStrand());
				regIntron.setChromosome(ss[0].getTranscripts()[0].getChromosome());
				regIntron.setSpecies(ss[0].getTranscripts()[0].getSpecies());

					// determine IRF
				int cntRetained= 0, cntNotRetained= 0;
				for (int j = 0; j < g.getTranscripts().length; j++) {
					DirectedRegion[] introns= g.getTranscripts()[j].getIntrons();
					for (int k = 0; introns!= null&& k < introns.length; k++) 
						if (introns[k].equals(regIntron)) {
							++cntNotRetained;
							break;
						}
				}
				for (int j = 0; j < g.getTranscripts().length; j++) {
					DirectedRegion[] exons= g.getTranscripts()[j].getExons();
					for (int k = 0; k < exons.length; k++) 
						if (exons[k].contains(regIntron)) {
							++cntRetained;
							break;
						}
				}
				
					// get sequence				
				int don5Extent= SpliceSite.NOBORU_DON5_EXTENT;
				DirectedRegion reg= ss[0].getNoboruRegion();
				String seqDon= Graph.readSequence(reg);
				String dinuc= seqDon.substring(don5Extent, don5Extent+2).toLowerCase();
				if (gtagOnly&& !dinuc.equals("gt"))
					continue;
				seqDon=
					seqDon.substring(0, don5Extent).toUpperCase()+
					dinuc+
					seqDon.substring(don5Extent+2, seqDon.length()).toUpperCase();
				
				int acc5Extent= SpliceSite.NOBORU_ACC5_EXTENT;				
				reg= ss[1].getNoboruRegion();
				String seqAcc= Graph.readSequence(reg);
				dinuc= seqAcc.substring(acc5Extent, acc5Extent+ 2).toLowerCase();
				if (gtagOnly&& !dinuc.equals("ag"))
					continue;
				seqAcc=
					seqAcc.substring(0, acc5Extent).toUpperCase()+
					dinuc+
					seqAcc.substring(acc5Extent+2, seqAcc.length()).toUpperCase();

					// sort/add
				if (cntNotRetained> cntRetained) {
					lowIRFdonorV.add(ss[0]);
					lowIRFacceptorV.add(ss[1]);
					lowIRFintronV.add(regIntron);
				} else  if (cntRetained> cntNotRetained) {
					highIRFdonorV.add(ss[0]);
					highIRFacceptorV.add(ss[1]);
					highIRFintronV.add(regIntron);
				} else {	// equal
					map.put(ID_NB_SS_EQ_IRF_FREQ, new Integer(++cntEQfreqIR));
				}
			}
		}
		
		
		if (fName!= null) {
			try {
				System.out.println("Found "+nbEqIRFfreq.intValue()+" IR events with equal freqs.");
				
				BufferedWriter buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_donors_gtag-"+gtagOnly));
				for (int i = 0; i < lowIRFdonorV.size(); i++) {
					buffy.write(">Donor|don"+i+"|loIRF|gtag-"+gtagOnly+"\n");
					buffy.write(lowIRFdonorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_acceptors_gtag-"+gtagOnly));
				for (int i = 0; i < lowIRFacceptorV.size(); i++) { 
					buffy.write(">Acceptor|acc"+i+"|loIRF|gtag-"+gtagOnly+"\n");
					buffy.write(lowIRFacceptorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_loIRF_introns_gtag-"+gtagOnly));
				buffy.write("intronNr\tlength\taff_5UTR\taff_CDS\taff_3UTR\n");
				for (int i = 0; i < lowIRFintronV.size(); i++) { 
					buffy.write("intron"+i+"\t"+lowIRFintronV.elementAt(i)+"\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_donors_gtag-"+gtagOnly));
				for (int i = 0; i < highIRFdonorV.size(); i++) {
					buffy.write(">Donor|don"+i+"|hiIRF|gtag-"+gtagOnly+"\n");
					buffy.write(highIRFdonorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_acceptors_gtag-"+gtagOnly));
				for (int i = 0; i < highIRFacceptorV.size(); i++) { 
					buffy.write(">Acceptor|acc"+i+"|hiIRF|gtag-"+gtagOnly+"\n");
					buffy.write(highIRFacceptorV.elementAt(i)+"\n\n");
				}
				buffy.flush(); buffy.close();
				
				buffy= new BufferedWriter(new FileWriter(fName+"_hiIRF_introns_gtag-"+gtagOnly));
				buffy.write("intronNr\tlength\taff_5UTR\taff_CDS\taff_3UTR\n");
				for (int i = 0; i < highIRFintronV.size(); i++) { 
					buffy.write("intron"+i+"\t"+highIRFintronV.elementAt(i)+"\n");
				}
				buffy.flush(); buffy.close();
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}
	
	
	public Noboru(String newSpeName, String newAnnoName) {
		this.speName= newSpeName;
		this.annoName= newAnnoName;
	}
	
	/** 
	 * data for 5ss statistics -3 +6
	 * mint=sum of all lowest nt frequencies along -3 +6 
	 * maxt=sum of all highest nt frequencies along -3 +6 
	 * nt frequencies at each pos from -3 to +6
	 **/
	public double[][] getDonorFreqTable() {
		if (donTab == null) {
			String fName= Constants.getAnnotationFile(speName, annoName);
			File f= new File(fName+SFX_DONOR_TABLE);
			if (!f.exists()) { 
				generateFreqTables();
				writeFreqTables();	// dkslfgjklas klsdjf
			} else {
				donTab= new double[9][];
				double[] minMax= readFreqTable(donTab, f);
				donMinT= minMax[0];
				donMaxT= minMax[1];
			}
		}
	
		return donTab;
	}

	/** 
	 * data for 5ss statistics -3 +6
	 * mint=sum of all lowest nt frequencies along -3 +6 
	 * maxt=sum of all highest nt frequencies along -3 +6 
	 * nt frequencies at each pos from -3 to +6
	 **/
	public double[][] getAcceptorFreqTable() {
		if (accTab == null) {
			String fName= Constants.getAnnotationFile(speName, annoName);
			File f= new File(fName+SFX_ACCEPTOR_TABLE);
			if (!f.exists()) {
				generateFreqTables();
				writeFreqTables();
			} else {
				accTab= new double[14][];
				double[] minMax= readFreqTable(accTab, f);
	//			($t, $l1, $l2, $h1, $h2)=split/\s+/;
				accL1= minMax[0];
				accL2= minMax[1];
				accH1= minMax[2];
				accH2= minMax[3];
			}
		}

		return accTab;
	}

	void generateFreqTables() {
		try {
			String fName= Constants.getLatestUCSCAnnotation(speName, annoName, null);
			System.out.println("Generating freq table for "+fName);
			GTFChrReader reader= new GTFChrReader(fName);
			reader.setChromosomeWise(true);
			reader.setReadGene(true);
			reader.read();
			Gene[] genes= reader.getGenes();
			HashMap localHash= new HashMap();
			while (genes!= null) {
				for (int i = 0; i < genes.length; i++) 
					_02_outputSpliceSites(genes[i], null, localHash);
				reader.read();
				genes= reader.getGenes();
			}			
			
			final String ID_CONST_DONOR= "const_donor"; 
			final String ID_CONST_ACCEPTOR= "const_acceptor";
			Vector constDonor= (Vector) localHash.get(ID_CONST_DONOR);
			Vector constAcceptor= (Vector) localHash.get(ID_CONST_ACCEPTOR);
			countFreq((SpliceSite[]) gphase.tools.Arrays.toField(constDonor));
			countFreq((SpliceSite[]) gphase.tools.Arrays.toField(constAcceptor));
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/**
	 * 		# data for 5ss statistics -3 +6
	 * # mint=sum of all lowest nt frequencies along -3 +6
	 * # maxt=sum of all highest nt frequencies along -3 +6
	 * # nt frequencies at each pos from -3 to +6
	 * 
	 */
	// TODO perl porting, make strucuture nice
	strictfp void countFreq(SpliceSite[] ss) {

		if (ss== null|| ss.length< 1)
			return;
		
			// count frequencies
		int totDonors= 0, totAcceptors= 0;
		boolean initDon= false;
		boolean initAcc= false;
		for (int i = 0; i < ss.length; i++) {
			if (ss[i].isDonor()&& donTab== null) {
				initDon= true;
				donTab= new double[ss[i].getShapiroRegion().getLength()][];
				for (int j = 0; j < donTab.length; j++) {
					donTab[j]= new double[4];
					for (int k = 0; k < donTab[i].length; k++) 
						donTab[j][k]= 0d;
				}
			} else if (ss[i].isAcceptor()&& accTab== null) {
				initAcc= true;
				accTab= new double[ss[i].getShapiroRegion().getLength()][];
				for (int j = 0; j < accTab.length; j++) {
					accTab[j]= new double[4];
					for (int k = 0; k < accTab[i].length; k++) 
						accTab[j][k]= 0d;
				}
			}
			
			String seq= Graph.readSequence(ss[i].getShapiroRegion());
			for (int j = 0; j < seq.length(); j++) {
				for (int k = 0; k < BASE_ORDER.length(); k++) 
					if (BASE_ORDER.charAt(k)== seq.charAt(j)) {
						if (ss[i].isDonor())
							++donTab[j][k];
						else
							++accTab[j][k];
						break;
					}
				
			}
			if (ss[i].isDonor())
				++totDonors;
			else
				++totAcceptors;
		}
		
			// convert to fractions
		if (initDon) {
			donMinT= 0d; donMaxT= 0d;
			for (int i = 0; i < donTab.length; i++) {
				int min= Integer.MAX_VALUE;
				int max= Integer.MIN_VALUE;
				for (int j = 0; j < donTab[i].length; j++) {
					min= Math.min(min, (int) donTab[i][j]);
					max= Math.max(max, (int) donTab[i][j]);		// concat here for consensus
						//100*$freq[$i]->{A}/$total, "\n";
					donTab[i][j]= 100d* donTab[i][j]/ totDonors;
				}
					// $mint+=100*$min/$total;
				donMinT+= 100d* min/ totDonors;
				donMaxT+= 100d* max/ totDonors;
			}
		}
		if (initAcc) {
			double[] min10= new double[accTab.length- 4];
			double[] max10= new double[accTab.length- 4];
			for (int i = 0; i < accTab.length- 4; i++) {
				int min= Integer.MAX_VALUE;
				int max= Integer.MIN_VALUE;
				for (int j = 0; j < accTab[i].length; j++) {
					min= Math.min(min, (int) accTab[i][j]);
					max= Math.max(max, (int) accTab[i][j]);		// concatenate here letters for consensus seq
						//100*$freq[$i]->{A}/$total, "\n";
					accTab[i][j]= 100d* accTab[i][j]/ totAcceptors;
				}
					// $mint+=100*$min/$total;
				min10[i]= 100d* min/ totAcceptors;
				max10[i]= 100d* max/ totAcceptors;
			}
			
			Arrays.sort(min10);
			accL1= 0d;
			for (int i = 0; i < min10.length; i++) 
				accL1+= min10[i];
	//	    @top8=reverse sort {$a <=> $b} @high;
	//	    for(my $k=0; $k<8; $k++){ $h1+=$top8[$k]; }
			Arrays.sort(max10);
			accH1= 0d;
			for (int i = max10.length- 1; i> 1; --i) 
				accH1+= max10[i];
			
			accL2= 0d; accH2= 0d;
			for (int i = accTab.length- 4; i < accTab.length; i++) {
				int min= Integer.MAX_VALUE;
				int max= Integer.MIN_VALUE;
				for (int j = 0; j < accTab[i].length; j++) {
					min= Math.min(min, (int) accTab[i][j]);
					max= Math.max(max, (int) accTab[i][j]);		// concatenate here letters for consensus seq
						//100*$freq[$i]->{A}/$total, "\n";
					accTab[i][j]= 100d* accTab[i][j]/ totAcceptors;
				}
					// $mint+=100*$min/$total;
				accL2+= 100d* min/ totAcceptors;
				accH2+= 100d* max/ totAcceptors;
			}
		}
	}
	
	
	void writeFreqTables() {
		try {
			BufferedWriter writer= new BufferedWriter(new FileWriter(getInFile()+SFX_DONOR_TABLE));
			for (int i = 0; i < donTab.length; i++) 
				for (int j = 0; j < donTab[i].length; j++) 
					writer.write("pos "+i+" freq"+BASE_ORDER.charAt(j)+" "+donTab[i][j]+"\n");
			writer.write("mint "+donMinT+" maxt "+donMaxT+"\n");
			writer.flush(); writer.close();
			
			writer= new BufferedWriter(new FileWriter(getInFile()+SFX_ACCEPTOR_TABLE));
			for (int i = 0; i < accTab.length; i++) 
				for (int j = 0; j < accTab[i].length; j++) 
					writer.write("pos "+i+" freq"+BASE_ORDER.charAt(j)+" "+accTab[i][j]+"\n");
			//			($t, $l1, $l2, $h1, $h2)=split/\s+/;
			writer.write("l1 "+accL1+" "+accL2+" "+accH1+" "+accH2+"\n");
			writer.flush(); writer.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	double[] readFreqTable(double[][] tab, File f) {
		
		double[] result= null;
		for (int i = 0; i < tab.length; i++) 
			tab[i]= new double[4];		
		try {
			BufferedReader reader= new BufferedReader(new FileReader(f));
			Vector distrV= new Vector();
			int ctr= 0;
			int cnt= 0;
			while (reader.ready()) {
				String line= reader.readLine();
				String[] tokens= line.split("\\s+");				
				//pos 0 freqA 36.1182247491953
				if (tokens[0].equals("pos")) {
					int pos= Integer.parseInt(tokens[1]);
					char lett= tokens[2].charAt(tokens[2].length()- 1);
					for (int i = 0; i < BASE_ORDER.length(); i++) 
						if (BASE_ORDER.charAt(i)== lett) {
							tab[pos][i]= Double.parseDouble(tokens[3]);
							break;
						}
				} else if (tokens[0].equals("mint")) {
					result= new double[2];
					result[0]= Double.parseDouble(tokens[1]);
					result[1]= Double.parseDouble(tokens[3]);
//					l1 67.1490944851569 11.3057541810082 413.248301234012 313.96098318118
				} else if (tokens[0].equals("l1")) {	
					result= new double[4];
					result[0]= Double.parseDouble(tokens[1]);
					result[1]= Double.parseDouble(tokens[2]);
					result[2]= Double.parseDouble(tokens[3]);
					result[3]= Double.parseDouble(tokens[4]);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return result;
	}
//	# data for 5ss statistics -3 +6
//	# mint=sum of all lowest nt frequencies along -3 +6
//	# maxt=sum of all highest nt frequencies along -3 +6
//	# nt frequencies at each pos from -3 to +6
	public strictfp double[] scoreSpliceSite(SpliceSite[] ss) {
//		print STDERR "xxxxxxxxxxxxxxxxxCAGGTAAGTxxxxxxxxxxxxxx\n";
//		print STDERR "xxxxxxxxxxTTTTTTTTTTCAGGxxxxxxxxxxxxxxxx\n";
//
//		if($site eq "5ss+2+5"){$start=22; $end=26; }
//		if($site eq "5ss"){$start=17; $end=26;}
//		if($site eq "3ss"){$start=10; $end=20; }
//
		double[] scores= new double[ss.length];
		for (int i = 0; i < ss.length; i++) {
			double[][] table= null;
			if (ss[i].isDonor()) {
				table= getDonorFreqTable();
			} else {
				table= getAcceptorFreqTable();
			}
			
			String seq= Graph.readSequence(ss[i].getShapiroRegion());
			
			
			double score= -1d;
			if (ss[i].isDonor()) {
				double sum= 0d;
				for (int j = 0; j < seq.length(); j++) {
					for (int k = 0; k < BASE_ORDER.length(); k++) 
						if (BASE_ORDER.charAt(k)== seq.charAt(j)) {
							sum+= table[j][k];
							break;
						}
				}
			    score=100*(sum- donMinT)/(donMaxT-donMinT);
			    //			  acceptors
			} else {
				double[] surface= new double[seq.length()- 4];
				for (int j = 0; j < seq.length()- 4; j++) {
					for (int k = 0; k < BASE_ORDER.length(); k++) 
						if (BASE_ORDER.charAt(k)== seq.charAt(j)) {
							surface[j]= table[j][k];
							break;
						}
				}
				Arrays.sort(surface);
				double sum= 0d;
				for (int j = surface.length- 1; j > surface.length- 9; --j) 
					sum+= surface[j];
				
				
	//		    for($i=20, $j=10; $i<24; $i++, $j++){
				double sum2= 0d;
				for (int j = seq.length()- 4; j < seq.length(); j++) {
					for (int k = 0; k < BASE_ORDER.length(); k++) 
						if (BASE_ORDER.charAt(k)== seq.charAt(j)) {
							sum2+= table[j][k];
							break;
						}
				}
//			    $score=100*( (($t1-$l1)/($h1-$l1))+(($t2-$l2)/($h2-$l2)) )/2;
			    score= 100*( ((sum- accL1)/(accH1-accL1))+
			    		((sum2- accL2)/ (accH2- accL2)) ) / 2;
			}
			
		    scores[i]= score;
		}
		
		return scores;
	}
	public String getInFile() {
		if (inFile == null) {
			inFile = Constants.getLatestUCSCAnnotation(speName, annoName, null);
		}

		return inFile;
	}
}
