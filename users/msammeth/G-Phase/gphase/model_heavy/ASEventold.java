/*
 * Created on Mar 6, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model_heavy;

import java.io.Serializable;
import java.util.Comparator;
import java.util.Vector;

/**
 * alternate donor/acceptor are overlapping, if not -> two exon skipping
 * 
 * 
 * @author msammeth
 */
public class ASEventold implements Serializable {

	public static class PositionComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			ASEventold asEv1= (ASEventold) arg0;
			ASEventold asEv2= (ASEventold) arg1;
			
			int ss1= Math.min(asEv1.ss1.getPos(),asEv1.ss2.getPos());
			int ss2= Math.min(asEv2.ss1.getPos(),asEv2.ss2.getPos());
			
			if (ss1< ss2)
				return (-1);
			if (ss1> ss2)
				return 1;
			return 0;
		}
	}
	
	public final static String TYPE_ALTERNATE_ACCEPTOR= "AA";
	public final static String TYPE_ALTERNATE_DONOR= "AD";
	public final static String TYPE_EXON_SKIPPING= "ES";
	public final static String TYPE_MUTUALLY_EXCLUSIVE= "ME";	// mutually exclusive exons: 2 events
	public final static String TYPE_INTRON_RETENTION= "IR";
	
	ASVariation variation= null;
	int transcriptNr= -1; // 0= both, 1= pattern.trans1, 2= pattern.trans2
	SpliceSite ss1= null, ss2= null; // transcriptNr= 0: ss1 in pattern.trans1, ss2 in pattern.trans2
	String type= null;
	
	public static ASEventold[] toArray(Vector v) {
		ASEventold[] result= new ASEventold[v.size()];
		for (int i = 0; i < result.length; i++) 
			result[i]= (ASEventold) v.elementAt(i);
		
		return result;
	}
	
	
	public ASEventold(ASVariation newVar, int newTransNr, SpliceSite newSS1, SpliceSite newSS2) {
		
		this.variation= newVar;
		this.transcriptNr= newTransNr;
		this.ss1= newSS1;
		this.ss2= newSS2;
		if (transcriptNr> 0&& newSS2.getPos()< newSS1.getPos()) {	// order ss in same transcript
			this.ss1= newSS2;
			this.ss2= newSS1;
		}
		
		this.type= getType();
	}
	
	String getType() {
		if (transcriptNr== 0) {	// two transcripts involved: alternate donor/acceptor
			if (ss1.isDonor()&& ss2.isDonor())
				return TYPE_ALTERNATE_DONOR;
			if (ss1.isAcceptor()&& ss2.isAcceptor())
				return TYPE_ALTERNATE_ACCEPTOR;
			System.err.println("Unknown splice type ("+ss1+","+ss2
					+") in "+variation.getTranscript1()+","+variation.getTranscript2());
		} else {	// only one transcript involved: exon skipping/ intron retention
			if (ss1.isDonor()&& ss2.isAcceptor())
				return TYPE_INTRON_RETENTION;
			if (ss1.isAcceptor()&& ss2.isDonor())
				return TYPE_EXON_SKIPPING;
			System.err.println("Unknown splice type ("+ss1+","+ss2
					+") in "+((transcriptNr== 1)?variation.getTranscript1():variation.getTranscript2()));
		}
		return null;
	}

	public boolean isAlternateAcceptor() {
		return type.equals(TYPE_ALTERNATE_ACCEPTOR);
	}

	public boolean isAlternateDonor() {
		return type.equals(TYPE_ALTERNATE_DONOR);
	}

	public boolean isExonSkipping() {
		return type.equals(TYPE_EXON_SKIPPING);
	}

	public boolean isMutuallyExclusive() {
		return type.equals(TYPE_MUTUALLY_EXCLUSIVE);
	}

	public boolean isIntronRetention() {
		return type.equals(TYPE_INTRON_RETENTION);
	}
	public int getTranscriptNr() {
		return transcriptNr;
	}
	
	public void setType(String newType) {
		if (!newType.equals(TYPE_MUTUALLY_EXCLUSIVE))
			return;		// allow only mutually exclusive to be set
						// (uses 2 exon skipping splice events)
		this.type = newType;
	}
	
	public String toString() {
		return type;
	}
}
