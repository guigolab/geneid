package gphase.model.intransparent.init;

import java.util.Comparator;

public class Intron extends Exon {

	public Intron(Transcript newTranscript, int start,
			int end) {
		
		setStrand(newTranscript.getStrand());
		setStart(start);
		setEnd(end);
		addTranscript(newTranscript);
	}
	
	@Override
	public SpliceSite getDonor() {
		return getTranscripts()[0].getUpstreamExon(this).getDonor();
	}
	
	@Override
	public SpliceSite getAcceptor() {
		return getTranscripts()[0].getDownstreamExon(this).getAcceptor();
	}
	
	@Override
	public boolean addTranscript(Transcript newTranscript) {
		if(transcripts== null) {
			transcripts= new Transcript[] {newTranscript};
			return true;
		}
		
			// add transcript		
		Transcript[] nTranscripts= new Transcript[transcripts.length+ 1];
		for (int i= 0; i < transcripts.length; i++) {
			if (transcripts[i].getTranscriptID().equals(newTranscript.getTranscriptID()))
				return false;
			nTranscripts[i]= transcripts[i];
		}
		nTranscripts[nTranscripts.length- 1]= newTranscript;
		transcripts= nTranscripts;
		
		return true;
	}

	public boolean isInternal() {
		for (int i = 0; i < getTranscripts().length; i++) {
			if (getTranscripts()[i].getUpstreamExon(this).isInternal()&&
					getTranscripts()[i].getDownstreamExon(this).isInternal())
				return true;
		}
		return false;
	}
	
}
