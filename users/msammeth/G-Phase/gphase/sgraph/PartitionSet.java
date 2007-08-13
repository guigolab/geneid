package gphase.sgraph;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import sun.security.util.PathList;

public class PartitionSet {
	Vector<HashMap<Path, Path>> bucketList;	// non-redundant
	HashMap<Path, HashMap<Path,Path>> pathBucketMap;	// for extension
	
	public PartitionSet(int initialSize) {
		bucketList= new Vector<HashMap<Path, Path>>(initialSize);
		pathBucketMap= new HashMap<Path, HashMap<Path,Path>>(initialSize, 1f);		
	}
	
	public void addPath(Path p) {
		p.setPartition(this);
		int trptCnt= Graph.getTranscriptNb(p.getTranscripts());
		HashMap<Path, Path> bucket= new HashMap<Path, Path>(trptCnt,1f);
		bucket.put(p,p);
		pathBucketMap.put(p,bucket);
		bucketList.add(bucket);
	}
	
	public void splitPath(Path oldPath, Vector<Path> newPathes) {
		HashMap<Path, Path> bucket= pathBucketMap.remove(oldPath);
		bucket.remove(oldPath);
		for (int i = 0; i < newPathes.size(); i++) {
			bucket.put(newPathes.elementAt(i), newPathes.elementAt(i));
			pathBucketMap.put(newPathes.elementAt(i),bucket);
		}
		//nothing to do about non-/mappable, bucket has same OID
	}
	
	public Vector<Vector<Path>> splitBucket(Vector<Path> bucket) {
		HashMap<HashMap<Path, Path>, Vector<Path>> map= new HashMap<HashMap<Path, Path>, Vector<Path>>(bucket.size(),1f);
		for (int i = 0; i < bucket.size(); i++) {
			HashMap<Path, Path> part= pathBucketMap.get(bucket.elementAt(i));
			Vector<Path> v= map.get(part);
			if (v== null)
				map.put(part, new Vector<Path>(bucket.size()));
			map.get(part).add(bucket.elementAt(i));
		}
		
		return new Vector<Vector<Path>>(map.values());
	}
	
	public boolean mergeBuckets(Vector<Path> pathes) {
		
		HashMap<HashMap<Path, Path>,HashMap<Path, Path>> oldBucketH= new HashMap<HashMap<Path, Path>,HashMap<Path, Path>>(pathes.size());
		int len= 0;
		for (int i = 0; i < pathes.size(); i++) {
			HashMap<Path, Path> bucket= pathBucketMap.get(pathes.elementAt(i));
			len+= bucket.size();
			oldBucketH.put(bucket,bucket);
		}
		
		if (oldBucketH.size()<= 1)
			return false;	// nothing to do
		if (oldBucketH.size()== bucketList.size())
			return true;
		
		Vector<HashMap<Path, Path>> oldBuckets= new Vector(oldBucketH.values());
		HashMap<Path, Path> newBucket= new HashMap<Path, Path>(len+1,1f);
		for (int i = 0; i < oldBuckets.size(); i++) {
			bucketList.remove(oldBuckets.elementAt(i));	// iterates bucketList !! inefficient
			Iterator<Path> iter= oldBuckets.elementAt(i).keySet().iterator();
			while(iter.hasNext()) {
				Path p= iter.next();
				newBucket.put(p,p);
				pathBucketMap.remove(p);
				pathBucketMap.put(p,newBucket);
			}
		}
		bucketList.add(newBucket);
		return false;
	}
	
	
	
}
