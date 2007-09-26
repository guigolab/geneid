package gphase.graph9;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;


public class Partition {
	HashMap<PartitionSet, PartitionSet> parents= null;
	long[] transcripts= null;

	public Partition() {
		parents= new HashMap<PartitionSet,PartitionSet>(4,1f);
	}
	
	public void addParent(PartitionSet newParent) {
		newParent.partitions.put(this,this);
		parents.put(newParent,newParent);
	}
	
	@Override
	public String toString() {
		return transcripts.toString();
	}
	
	public Partition clonePartition() {
		
		Partition p= new Partition();
		p.parents= (HashMap<PartitionSet, PartitionSet>) parents.clone();
		Iterator<PartitionSet> iter= parents.keySet().iterator();
		while (iter.hasNext())
			iter.next().partitions.put(p, p);
		return p;
	}
}


