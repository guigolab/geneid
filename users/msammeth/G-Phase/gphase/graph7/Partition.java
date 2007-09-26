package gphase.graph7;

import java.util.HashMap;

import com.sun.org.apache.xalan.internal.xsltc.runtime.Hashtable;

public class Partition {

	long[] transcripts;
	HashMap<Path,Partition> nextHash;
	
	Partition(long[] t) {
		this.transcripts= t;
		nextHash= new HashMap<Path, Partition>(5,1f);
	}
	
}
