package gphase.tools;

public class Double implements Comparable {
	
	double val= 0;
	
	public Double(double newVal) {
		val= newVal;
	}
	
	public double doubleValue() {
		return val;
	}
	
	public int compareTo(Object o) {
		double val2= ((Double) o).doubleValue();
		
		if (val< val2)
			return -1;
		if (val> val2)
			return 1;
		
		return 0;
	}
}
