package gphase.model_heavy;

import gphase.tools.Arrays;

import java.util.Comparator;
import java.util.Vector;

import org.freehep.graphicsio.swf.SWFAction.Equals2;

public class VariantGroup {

	SpliceSite[] spliceChain= null; 
	Vector transcriptIDV= null;
	ASVariation var= null;
	
	public VariantGroup(ASVariation event) {
		this.var= event;
	}
	
	public void addTranscriptID(String newTranscriptID) {
		if (transcriptIDV== null)
			transcriptIDV= new Vector();
		for (int i = 0; i < transcriptIDV.size(); i++) 
			if (transcriptIDV.elementAt(i).equals(newTranscriptID))
				return;
		transcriptIDV.add(newTranscriptID);
	}
	
	public String[] getTranscriptIDs() {
		return (String[]) Arrays.toField(transcriptIDV);
	}
	
	public ASVariation getASevent() {
		return var;
	}
	
	@Override
	public int hashCode() {
		return var.getSpliceUniverse()[0].hashCode();
	}
	
	@Override
	public boolean equals(Object obj) {
		Comparator compi= new ASVariation.IdentityComparator();
		return (compi.compare(this, obj)== 0);
	}
	
}
