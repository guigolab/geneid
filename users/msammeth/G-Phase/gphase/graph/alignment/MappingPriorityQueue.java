package gphase.graph.alignment;

import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.SortedSet;

public class MappingPriorityQueue extends PriorityQueue {

	HashMap alignedPosMap= new HashMap();
	
	@Override
	public Object poll() {
		Mapping map= (Mapping) super.poll();
		String key= map.getMaxI()+","+map.getMaxJ();
		alignedPosMap.remove(key);
		return map;
	}
	
	@Override
	public boolean add(Object o) {
		
		Mapping map= (Mapping) o;
		if (map.maxRel== 3) {
			String key= map.getMaxI()+","+map.getMaxJ();
			Mapping altMap= (Mapping) alignedPosMap.get(key);
			if (altMap== null) {
				alignedPosMap.put(key, map);
				return super.add(o);	// always true, see Collection policy
			}
			if (altMap.getCost()> map.getCost()) {	// else
				remove(altMap);
				alignedPosMap.put(key, map);
				return super.add(o);
			}
			return false;	// else dont add, cost>= altCost
		} 
		
		return super.add(o);	// else
	}
	
	public MappingPriorityQueue() {
		// TODO Auto-generated constructor stub
	}

	public MappingPriorityQueue(int initialCapacity) {
		super(initialCapacity);
		// TODO Auto-generated constructor stub
	}

	public MappingPriorityQueue(Collection c) {
		super(c);
		// TODO Auto-generated constructor stub
	}

	public MappingPriorityQueue(PriorityQueue c) {
		super(c);
		// TODO Auto-generated constructor stub
	}

	public MappingPriorityQueue(SortedSet c) {
		super(c);
		// TODO Auto-generated constructor stub
	}

	public MappingPriorityQueue(int initialCapacity, Comparator comparator) {
		super(initialCapacity, comparator);
		// TODO Auto-generated constructor stub
	}

}
