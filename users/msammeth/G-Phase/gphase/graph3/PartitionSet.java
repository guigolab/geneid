package gphase.graph3;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import sun.security.util.PathList;

public class PartitionSet {
	
	Vector<long[]> partitions= null;
	HashMap<long[],long[]>[] splits= null;
	
	public PartitionSet(Vector<long[]> partitions) {
		setPartitions(partitions);
	}

	public Vector<long[]> getPartitions() {
		return partitions;
	}

	public void setPartitions(Vector<long[]> partitions) {
		this.partitions = partitions;
		splits= new HashMap[partitions.size()];
	}
	
	public void mergePartitions(Vector<long[]> oldPartitions) {
		long[] newPartition= Graph.unite(oldPartitions);
		for (int i = 0; i < oldPartitions.size(); i++) {
			partitions.remove(oldPartitions.elementAt(i));
		}
		partitions.add(newPartition);
	}
	
	public void splitPathes(PartitionSet ps) {
		Vector<long[]> inPart= ps.getPartitions();
		for (int i = 0; i < partitions.size(); i++) {
			for (int j = 0; j < inPart.size(); j++) {
				long[] inter= Graph.intersect(partitions.elementAt(i), inPart.elementAt(j));
				if (Graph.isNull(inter)||
						Graph.equalSet(inter, Graph.intersect(inter, partitions.elementAt(i))))
					continue;	// inPart does not split this partition
				if (splits[i]== null) 
					splits[i]= new HashMap<long[],long[]>(5);
				splits[i].put(inter,inter);
			}
		}
	}
	
	public long[][][] getSplitPartitions() {
		long[][][] sp= new long[partitions.size()][][];
		for (int i = 0; i < sp.length; i++) {
			if (splits[i]== null) {
				sp[i]= new long[][] {partitions.elementAt(i)};
				continue;
			}
			
			long[] unity= null;
			Iterator<long[]> iter= splits[i].keySet().iterator();
			while (iter.hasNext())  {
				long[] nxt= iter.next();
				if (unity== null) 
					unity= nxt.clone();
				else
					unity= Graph.unite(unity, nxt);
			}
			long[] rest= Graph.without(partitions.elementAt(i),unity);
			int c= 0;
			if (!Graph.isNull(rest)) {
				sp[i]= new long[splits[i].size()+1][];
				sp[i][c++]= rest;
			} else
				sp[i]= new long[splits[i].size()][];

			iter= splits[i].keySet().iterator();
			while (iter.hasNext())  
				sp[i][c++]= iter.next();
		}
		
		return sp;
	}
}
