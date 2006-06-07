/*
 * Created on Mar 3, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.tools;

import gphase.model.ASVariation;
import gphase.model.SuperExon;

import java.lang.reflect.Array;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Collection;
import java.util.Comparator;
import java.util.Vector;


/**
 * 
 * 
 * @author msammeth
 */
public class Arrays {

	public static class FieldSizeRevComparator extends FieldSizeComparator{
		public int compare(Object arg0, Object arg1) {
			return -super.compare(arg0, arg1);
		}
	}
	public static class FieldSizeComparator implements Comparator {

		public int compare(Object arg0, Object arg1) {
			
			int size1= -1, size2= -1;
			
			if (arg0 instanceof Vector) 
				size1= ((Vector) arg0).size();
			else if (arg0 instanceof Object[])
				size1= ((Object[]) arg0).length;

			if (arg1 instanceof Vector) 
				size2= ((Vector) arg1).size();
			else if (arg1 instanceof Object[])
				size2= ((Object[]) arg1).length;
			
			if (size1< size2)
				return (-1);
			if (size2< size1)
				return 1;			
			return 0;
		}
	}
	
	public static void main(String[] args) {
		Integer[] a= new Integer[] {new Integer(2),new Integer(3),new Integer(5)};
		a= (Integer[]) insert(a, new Integer(4), 2);
		for (int i = 0; i < a.length; i++) 
			System.out.println(a[i].intValue());
	}

	public static Object[] sort2DField(Object o) {
		
		Object[] oo;
		if (o instanceof Collection)
			oo= (Object[]) ((Collection) o).toArray();
		else
			oo= (Object[]) o;
		
		java.util.Arrays.sort(oo, new Arrays.FieldSizeComparator());
		return oo;
	}
	
	
	public static Object[] sort2DFieldRev(Object o) {
		
		if (o== null)
			return null;
		
		Object[] oo;
		if (o instanceof Collection)
			oo= (Object[]) ((Collection) o).toArray();
		else
			oo= (Object[]) o;
		
		java.util.Arrays.sort(oo, new Arrays.FieldSizeRevComparator());
		return oo;
	}	
	
	
	public static Object[] add(Object[] a, Object o) {
		return insert(a, o, a.length);
	}
	
	public static Object[] extendField(Object[] a, Object o) {
//		if (a== null) {
//			Object[] newA= (Object[]) Array.newInstance(o.getClass(), 1);
//			newA[0]= o;
//			return newA;
//		}
		
		return insert(a, o, a.length);	
	}
	
	/**
	 * Assuming that a and o share the same class and some other things, this
	 * inserts an object in an array at the specified position. Automatically 
	 * converts (negative) insertion points as proposed by for instance 
	 * <code>Arrays.binarySearch</code>.
	 * 
	 * @param a
	 * @param o
	 * @param p
	 * @return
	 */
	public static Object[] insert(Object[] a, Object o, int p) {
		
		if (p< 0)
			p= (p+1)* (-1);
		
		Object[] newA= (Object[]) Array.newInstance(o.getClass(), a.length+ 1);
		for (int i = 0; i < p; i++) 
			newA[i]= a[i];
		newA[p]= o;
		for (int i = p+1; i < newA.length; i++) 
			newA[i]= a[i-1];
		
		return newA;
	}
	
	public static Object[] remove(Object[] a, Object o) {
		
		Vector v= new Vector();
		for (int i = 0; i < a.length; i++) 
			if (a[i]!= o)
				v.add(a[i]);
		
		return (Object[]) toField(v);
	}

	/**
	 * Converts an eventually highdimensional Array or Vector to an Array.
	 * @param base
	 * @return
	 */
	public static Object toField(Object base) {

		if (!(base instanceof Collection)&& !(base instanceof Object[])) {
			return base;
		}
		
		Object[] o= null;
		if (base instanceof Collection) 
			o= ((Collection) base).toArray();
		else 						// (base instanceof Object[]) 
			o= ((Object[]) base);

		if (o.length< 1)	// empty array
			return null; //cannot guess the basic class for vectors anyway
		
		Object r= null;
		int x= 0;
		while (r== null&& x< o.length) {
			r= toField(o[x]);
		}
		Object[] result= null;
		if (r!= null) {
			result= (Object[]) Array.newInstance(r.getClass(), o.length);
			result[0]= r;
			for (int i = 1; i < o.length; i++) 
				result[i]= toField(o[i]);
		}
		
		return result;
	}
		
	public static Collection merge(Object o1, Object o2) {
		Object[] oa1= null;
		if (o1 instanceof Collection)
			oa1= ((Collection) o1).toArray();
		else
			oa1= (Object[]) o1;
		
		Object[] oa2= null;
		if (o2 instanceof Collection)
			oa2= ((Collection) o2).toArray();
		else
			oa2= (Object[]) o2;
		
		//Object r= toField(oa1[0]);
		Object[] result= (Object[]) Array.newInstance(oa1[0].getClass(), oa1.length+ oa2.length);
		Vector resVec= new Vector(result.length);
		for (int i = 0; i < oa1.length; i++) 
			resVec.add(oa1[i]);
		for (int i = 0; i < oa2.length; i++) 
			resVec.add(oa2[i]);
		
		return resVec;
	}
}
