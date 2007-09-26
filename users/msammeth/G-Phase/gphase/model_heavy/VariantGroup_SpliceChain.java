package gphase.model_heavy;

import gphase.tools.Arrays;

import java.util.Vector;

public class VariantGroup_SpliceChain {

	SpliceSite[] spliceChain= null; 
	Vector transcriptIDV= null;
	Vector asEventV= null;
	
	public VariantGroup_SpliceChain(SpliceSite[] newSpliceChain) {
		this.spliceChain= newSpliceChain;
	}
	
	public void addTranscriptID(String newTranscriptID) {
		if (transcriptIDV== null)
			transcriptIDV= new Vector();
		for (int i = 0; i < transcriptIDV.size(); i++) 
			if (transcriptIDV.elementAt(i).equals(newTranscriptID))
				return;
		transcriptIDV.add(newTranscriptID);
	}
	
	public void addASevent(ASVariation newASevent) {
		if (asEventV== null)
			asEventV= new Vector();
		for (int i = 0; i < asEventV.size(); i++) 
			if (asEventV.elementAt(i).equals(newASevent))
				return;		
		asEventV.add(newASevent);
	}
	
	public String[] getTranscriptIDs() {
		return (String[]) Arrays.toField(transcriptIDV);
	}
	
	public ASVariation[] getASevents() {
		return (ASVariation[]) Arrays.toField(asEventV);
	}
	

	
	
	
}
