package gphase.graph9;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

public class PartitionSet {
	HashMap<Partition,Partition> partitions= null;
	
	public PartitionSet() {
		partitions= new HashMap<Partition,Partition>(4,1f);
	}

	public String toString() {
		StringBuffer b= new StringBuffer("[");
		Iterator<Partition> iter= partitions.keySet().iterator();
		while (iter.hasNext()) {
			b.append(iter.next());
			b.append(",");
		}
		b.deleteCharAt(b.length()- 1);
		b.append("] ; ");
		
		return b.toString();
	}
	
}
