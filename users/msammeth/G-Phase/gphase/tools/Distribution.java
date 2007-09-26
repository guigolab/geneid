package gphase.tools;

import gphase.model.DirectedRegion;

import java.util.Arrays;
import java.util.HashMap;

public class Distribution {

	double[] arrayD= null;	// sorted
	int[] arrayI= null;	// sorted
	HashMap histogram= null;
	       
	public Distribution(double[] surface) {
		Arrays.sort(surface);
		this.arrayD= surface;
	}
	
	public Distribution(int[] surface) {
		Arrays.sort(surface);
		this.arrayI= surface;
	}
	
	public double getMin() {
		if (arrayD!= null) {
			double min= java.lang.Double.MAX_VALUE;
			for (int i = 0; i < arrayD.length; i++) 
				if (arrayD[i]< min)
					min= arrayD[i];
			return min;
		} else {
			int min= java.lang.Integer.MAX_VALUE;
			for (int i = 0; i < arrayI.length; i++) 
				if (arrayI[i]< min)
					min= arrayI[i];
			return min;
		}
	}
	
	public double getMax() {
		if (arrayD!= null) {
			double max= java.lang.Double.MIN_VALUE;
			for (int i = 0; i < arrayD.length; i++) 
				if (arrayD[i]> max)
					max= arrayD[i];
			return max;
		} else {
			int max= java.lang.Integer.MIN_VALUE;
			for (int i = 0; i < arrayI.length; i++) 
				if (arrayI[i]> max)
					max= arrayI[i];
			return max;
		}
	}
	
	public HashMap getHistogram() {
		if (histogram == null) {
			histogram= new HashMap();
			for (int i = 0; arrayD!= null&& i < arrayD.length; i++) {
				Double key= new Double(arrayD[i]);
				Integer val= (Integer) histogram.get(key);
				if (val== null) 
					val= new Integer(1);
				else
					val= new Integer(val.intValue()+ 1);
				histogram.put(key, val);
			}
			
			for (int i = 0; arrayI!= null&& i < arrayI.length; i++) {
				Integer key= new Integer(arrayI[i]);
				Integer val= (Integer) histogram.get(key);
				if (val== null) 
					val= new Integer(1);
				else
					val= new Integer(val.intValue()+ 1);
				histogram.put(key, val);
			}
		}

		return histogram;
	}
	
	public double getMean() {
		
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return (sum/ arrayD.length);
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			double sum= 0d;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return (sum/ arrayI.length);
		}
		
		return 0d;
	}
	
	public double getTotal() {
		if (arrayD!= null) {
			double sum= 0;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return sum;
		} else if (arrayI!= null) {
			int sum= 0;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return sum;
		}
		
		return 0d;
	}
	
	double getMedian(int i) {
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			int medPos= arrayD.length/ i;
			if (arrayD.length% i== 0&& (medPos+1)< arrayD.length)
				return ((arrayD[medPos]+ arrayD[medPos+ 1])/ i);
			else
				return arrayD[medPos];
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			int medPos= arrayI.length/ i;
			if (arrayI.length% i== 0&& (medPos+1)< arrayI.length)
				return ((arrayI[medPos]+ arrayI[medPos+ 1])/ i);
			else
				return arrayI[medPos];
		}
		return 0d;
	}
	
	public double getMedian() {
		return getMedian(2);
	}
	public double get1stQuart() {
		return getMedian(4);
	}

	public double getSum() {		
		if (arrayD!= null) {
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) 
				sum+= arrayD[i];
			return sum;
		} else if (arrayI!= null) {
			int sum= 0;
			for (int i = 0; i < arrayI.length; i++) 
				sum+= arrayI[i];
			return sum;
		}
		return 0d;
	}
	
	public double getStandardDeviation() {
		if (arrayD!= null) {
			if (arrayD.length== 0)
				return 0d;
			double med= getMedian();		
			double sum= 0d;
			for (int i = 0; i < arrayD.length; i++) {
				double val= arrayD[i]- med;
				val*= val;
				sum+= val;
			}
			
			sum/= (arrayD.length- 1);
			return java.lang.Math.sqrt(sum);
			
			
		} else if (arrayI!= null) {
			if (arrayI.length== 0)
				return 0d;
			double med= getMedian();		
			double sum= 0d;
			for (int i = 0; i < arrayI.length; i++) {
				double val= arrayI[i]- med;
				val*= val;
				sum+= val;
			}
			
			sum/= (arrayI.length- 1);
			return java.lang.Math.sqrt(sum);
		}
		
		return 0d;
	}
	
	/**
	 * 
	 * @return mean, median, stdDev
	 */
	public String[] toStatString() {
		String[] statStr= new String[3];
		statStr[0]= new Double(getMean()).toString();
		statStr[1]= new Double(getMedian()).toString();
		statStr[2]= new Double(getStandardDeviation()).toString();
		for (int i = 0; i < statStr.length; i++) 
			statStr[i]= statStr[i].substring(0, statStr[i].indexOf('.')+ 2);

		return statStr;
	}
}
