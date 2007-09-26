package gphase.graph5;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import sun.security.util.PathList;

public class PartitionSet {
	
	Vector<Vector<Path>> partitions;
	
	public PartitionSet() {
	}

	public HashMap<Path,Path> getPartition(Edge srcEdge) {
		return partitionMap.get(srcEdge);
	}

	public void setPartitions(Set<Edge> edges) {
		Iterator<Edge> edgeIter= edges.iterator();
		while (edgeIter.hasNext()) {
			Edge e= edgeIter.next();
			Vector<Path> pathes= new Vector<Path>(e.getHead().getPathesFrom(e.getTail()));
			assert(pathes.size()== 1);
			HashMap<Path, Path> pMap= new HashMap<Path, Path>(2,1f);
			pMap.put(pathes.elementAt(0), pathes.elementAt(0));
		}
	}
	
	public void mergePartitions(Vector<Vector<Path>> part2merge) {
		
		
		long[] newPartition= Graph.unite(oldPartitions);
		for (int i = 0; i < oldPartitions.size(); i++) {
			partitions.remove(oldPartitions.elementAt(i));
		}
		partitions.add(newPartition);
	}

	public Vector<Vector<Path>> getPartitions() {
		return partitions;
	}

	public void setPartitions(Vector<Vector<Path>> partitions) {
		this.partitions = partitions;
	}
}
