package gphase.algo;

import gphase.tools.Arrays;

import java.lang.reflect.Method;
import java.util.Vector;

public class FilterFactory {
	
	static public final String AND= "AND";
	final static public String OR= "OR";
	public final static String NOT= "NOT"; 
	
	
	String[] mNames= null;
	Class c= null;
	int nBefore= 0, nAfter= 0;
	int[][] nMethod= null;
	
	
	public FilterFactory(String className, String methodsExec) {
		try {
			c= Class.forName(className);
		} catch (ClassNotFoundException e) {
			System.err.println("Class "+className+" not found!");;
		}
		addMethodString(methodsExec);
	}
	
	public void addMethodString(String mString){
		mNames= (String[]) Arrays.add(mNames, mString);
		nMethod= new int[mNames.length][];
//		for (int i = 0; i < nMethod.length; i++) 
//			nMethod[i]= 0;
	}
	
	/**
	 * recursive version
	 * @param o
	 * @return
	 */
	public Object filterRek(Object o) {
		if (mNames== null|| mNames.length== 0|| o== null)
			return o;
		
		if (o instanceof Object[]) {
			Object[] oo = (Object[]) o;
			Vector ooV= new Vector(oo.length);
			for (int i = 0; i < oo.length; i++) {
				Object filtOO= filterRek(oo[i]);
				if (filtOO!= null)
					ooV.add(filtOO);
			}
			if (ooV.size()> 0)
				return Arrays.toField(ooV.toArray());
		} else {
			++nBefore;
			++nAfter;
			boolean bool= true;
			for (int j = 0; j < mNames.length; j++) {				
				String[] tokens= mNames[j].split("_");
				if (nMethod[j]== null) {
					nMethod[j]= new int[tokens.length];
					for (int i = 0; i < nMethod[j].length; i++) 
						nMethod[j][i]= 0;
					
				}
				boolean andActive= false;
				boolean orActive= false;
				boolean notActive= false;
				for (int k = 0; k < tokens.length; k++) {
					if (tokens[k].equalsIgnoreCase(OR)) {
						nMethod[j][k]= -1;
						orActive= true;
						continue;
					}
					if (tokens[k].equalsIgnoreCase(AND)) {
						nMethod[j][k]= -1;
						andActive= true;
						continue;
					}
					if (tokens[k].equalsIgnoreCase(NOT)) {
						nMethod[j][k]= -1;
						notActive= true;
						continue;
					}
					try {
						Method m= c.getMethod(tokens[k], null);
						Boolean newBool= (Boolean) m.invoke(o, null);
						if (notActive)
							newBool= new Boolean(!newBool.booleanValue());
						boolean oldBool= bool;
						if (!orActive)	// AND or nothing (first)
							bool&= newBool.booleanValue();
						else if (orActive)
							bool|= newBool.booleanValue();
						if ((!bool)&& oldBool)
							++nMethod[j][k];
					} catch (Exception e) {
						e.printStackTrace();
					} finally {
						andActive= false;
						orActive= false;
						notActive= false;
					}
				}
				if (!bool) {
					//++nMethod[j];
					break;
				}
			}
			if (bool)
				return o;
			else
				--nAfter;
		}
		
		return null;
	}

	public Object[] filter(Object[] o) {
		if (mNames== null|| mNames.length== 0|| o== null|| o.length== 0)
			return o;
		
		nBefore+= o.length;
		Vector v= new Vector(o.length);
		for (int i = 0; i < o.length; i++) {
			boolean bool= true;
			for (int j = 0; j < mNames.length; j++) {				
				String[] tokens= mNames[j].split("_");
				if (nMethod[j]== null) {
					nMethod[j]= new int[tokens.length];
					for (int x = 0; x < nMethod[j].length; x++) 
						nMethod[j][x]= 0;
					
				}
				boolean andActive= false;
				boolean orActive= false;
				boolean notActive= false;
				for (int k = 0; k < tokens.length; k++) {
					if (tokens[k].equalsIgnoreCase(OR)) {
						nMethod[j][k]= -1;
						orActive= true;
						continue;
					}
					if (tokens[k].equalsIgnoreCase(AND)) {
						nMethod[j][k]= -1;
						andActive= true;
						continue;
					}
					if (tokens[k].equalsIgnoreCase(NOT)) {
						nMethod[j][k]= -1;
						notActive= true;
						continue;
					}
					try {
						Method m= c.getMethod(tokens[k], null);
						Boolean newBool= (Boolean) m.invoke(o[i], null);
						boolean oldBool= bool;
						if (notActive)
							newBool= new Boolean(!newBool.booleanValue());
						if (!orActive)	// AND or nothing (first)
							bool&= newBool.booleanValue();
						else if (orActive)
							bool|= newBool.booleanValue();
						if ((!bool)&& oldBool)
							++nMethod[j][k];
					} catch (Exception e) {
						System.err.println("Method "+tokens[k]+" not found!");
					} finally {
						andActive= false;
						orActive= false;
						notActive= false;
					}
				}
				if (!bool) {
					//++nMethod[j];
					break;
				}
			}
			if (bool)
				v.add(o[i]);
		}
		nAfter+= v.size();
		
		return (Object[]) Arrays.toField(v);
	}

	public int getNAfter() {
		return nAfter;
	}

	public int getNBefore() {
		return nBefore;
	}

	public int[][] getNMethod() {
		return nMethod;
	}
	
	public String toStringStats() {
		String result= getNBefore()+ " -> "+getNAfter()+ " (";
		for (int i = 0; i < getNMethod().length; i++) {
			for (int j = 0; j < getNMethod()[i].length; j++) 
				if (getNMethod()[i][j]!= -1)
					result+= getNMethod()[i][j]+ " ";
			result+= mNames[i];
			if (i< getNMethod().length- 1)
				result+= ",";
		}
		result+= ")";
		return result;
	}

	public String[] getMNames() {
		return mNames;
	}

	public Class getC() {
		return c;
	}
}
