package gphase.graph;

import java.util.Arrays;
import java.util.Vector;

public class Tuple {
	
	public static Tuple[] intersect(Tuple[] a, Tuple[] b) {
		Vector v= new Vector();
		for (int i = 0; i < a.length; i++) {
			int j;
			for (j = 0; j < b.length; j++) 
				if (a[i].equals(b[j]))
					break;
			if (j< b.length)
				v.add(a[i]);
		}
		
		return (Tuple[]) gphase.tools.Arrays.toField(v);
	}
	
	int[] tupleInt= null;
	Object[] tupleObj= null;
	
	public Tuple(int x, int y) {
		setTuple(new int[] {x, y});
	}
	
	public Tuple(Object x, Object y) {
		setTuple(new Object[] {x, y});
	}
	
	public void setTuple(int[] newTuple) {
		Arrays.sort(newTuple);
		this.tupleInt= newTuple;
	}
	
	public void setTuple(Object[] newTuple) {
		Arrays.sort(newTuple);
		this.tupleObj= newTuple;
	}
	
	public boolean equals(Object obj) {
		Tuple anotherTuple= (Tuple) obj;
		if (tupleObj== null) {
			int[] otherT= anotherTuple.getTupleInt();
			if ((otherT== null|| tupleInt== null)||
					(otherT.length!= tupleInt.length))
				return false;
			for (int i = 0; i < otherT.length; i++) 
				if (otherT[i]!= tupleInt[i])
					return false;
		} else {
			Object[] otherT= anotherTuple.getTupleObj();
			if ((otherT== null|| tupleObj== null)||
					(otherT.length!= tupleObj.length))
				return false;
			for (int i = 0; i < otherT.length; i++) 
				if (!otherT[i].equals(tupleObj[i]))
					return false;
		}
		
		return true;
	}

	public int[] getTupleInt() {
		return tupleInt;
	}
	public Object[] getTupleObj() {
		return tupleObj;
	}
}
