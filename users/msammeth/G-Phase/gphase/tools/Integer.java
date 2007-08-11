package gphase.tools;

public class Integer implements Comparable {
	
	int val= 0;
	
	public Integer(int newVal) {
		val= newVal;
	}
	
	public int intValue() {
		return val;
	}
	
	public int compareTo(Object o) {
		int val2= ((Integer) o).intValue();
		
		if (val< val2)
			return -1;
		if (val> val2)
			return 1;
		
		return 0;
	}
}
