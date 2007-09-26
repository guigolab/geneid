package gphase.graph8;

import java.util.Vector;


public class Partition {

	Vector<long[]> pathes= null;
	long[] unity= null;
	
	
	Partition(Vector<long[]> newPathes, long[] newUnity) {
		this.pathes= newPathes;
		this.unity= newUnity;
	}
	
	public String toString() {
		StringBuffer b= new StringBuffer("{");
		for (int i = 0; i < pathes.size(); i++) {
			b.append("[");
			for (int j = 0; j < pathes.elementAt(i).length; j++) {
				b.append(pathes.elementAt(i)[j]);
				b.append(",");
			}
			b.deleteCharAt(b.length()- 1);
			b.append("] ; ");
		}
		b.delete(b.length()- 3, b.length());
		b.append("}");
		
		return b.toString();
	}
}
