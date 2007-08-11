package gphase.graph;

import com.sun.org.apache.xpath.internal.operations.Equals;

import gphase.model.SpliceSite;
import gphase.model.Transcript;
import gphase.tools.Arrays;

public class SpliceEdge {
	SpliceNode head= null;
	SpliceNode tail= null;
	Transcript[] transcripts= null;
	
	SpliceEdge(SpliceNode newTail, SpliceNode newHead) {
		this.head= newHead;
		this.tail= newTail;
		head.addInEdge(this);
		tail.addOutEdge(this);
	}

	public SpliceEdge(SpliceNode newTail, SpliceNode newHead, Transcript trans) {
		this(newTail, newHead);
		addTranscript(trans);
	}
	
	public SpliceEdge(SpliceNode newTail, SpliceNode newHead, Transcript[] trans) {
		this(newTail, newHead);
		addTranscripts(trans);
	}
	
	public String toString() {
		return "["+ getTail()+ " -> "+ getHead()+"]";
	}

	public int getLength() {
		return (head.getSite().getPos()- tail.getSite().getPos()+ 1);	// hols also for neg strand
	}
	
	public boolean isExonic() {
		if (getTail().getSite() instanceof SpliceSite) {
			SpliceSite ss= (SpliceSite) getTail().getSite();
			if (ss.isAcceptor())
				return true;
			else
				return false;
		}
		
		if (getHead().getSite() instanceof SpliceSite) {
			SpliceSite ss= (SpliceSite) getHead().getSite();
			if (ss.isDonor())
				return true;
			else
				return false;
		}
		
			// two abstract sites
		return true;	// worx for tss->tes  
	}
	
	public void addTranscripts(Transcript[] newTranscripts) {
		for (int i = 0; i < newTranscripts.length; i++) 
			addTranscript(newTranscripts[i]);		
	}
	
	public void addTranscript(Transcript newTranscript) {
		for (int i = 0; transcripts!= null&& i < transcripts.length; i++) 
			if (transcripts[i]== newTranscript)
				return;
		transcripts= (Transcript[]) Arrays.extendField(transcripts, newTranscript);
	}
	
	public void removeTranscript(Transcript remTranscript) {
		transcripts= (Transcript[]) Arrays.remove(transcripts, remTranscript);
	}
	
	public boolean equals(Object obj) {
		
		SpliceEdge anotherEdge= (SpliceEdge) obj;
		if (anotherEdge.getHead()== head&& anotherEdge.getTail()== tail)
			return true;
		return false;
	}

	public SpliceNode getHead() {
		return head;
	}

	public SpliceNode getTail() {
		return tail;
	}

	public Transcript[] getTranscripts() {
		return transcripts;
	}
}
