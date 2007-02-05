package gphase.tools;

import gphase.model.DirectedRegion;

import java.util.Arrays;
import java.util.HashMap;

public class Distribution {

	double[] surface= null;	// sorted
	double[][] histogram= null;
	       
	public Distribution(double[] surface) {
		this.surface= surface;
		setSurface(this.surface);
	}
	
	public Distribution(int[] surface) {
		this.surface= new double[surface.length];
		for (int i = 0; i < surface.length; i++) 
			this.surface[i]= (double) surface[i];
		setSurface(this.surface);
	}
	
	public void setSurface(double[] surface) {
		Arrays.sort(surface);
		this.surface= surface;
	}
	
	public double[][] getHistogram() {
		if (histogram == null) {
			HashMap distrTable= new HashMap();
			for (int i = 0; i < surface.length; i++) {
				Double key= new Double(surface[i]);
				Integer val= (Integer) distrTable.get(key);
				if (val== null) 
					val= new Integer(1);
				else
					val= new Integer(val.intValue()+ 1);
				distrTable.put(key, val);
			}
			
			Object[] keys= distrTable.keySet().toArray();
			Arrays.sort(keys);
			histogram= new double[keys.length][];
			for (int i = 0; i < keys.length; i++) {
				histogram[i]= new double[2];
				histogram[i][0]= ((Double) keys[i]).doubleValue();
				histogram[i][1]= (double) ((Integer) distrTable.get(keys[i])).intValue();
			}
		}

		return histogram;
	}
	
	public double getMean() {
		if (surface.length== 0)
			return 0d;
		double sum= 0d;
		for (int i = 0; i < surface.length; i++) 
			sum+= surface[i];
		return (sum/ surface.length);
	}
	
	public double getMedian() {
		if (surface.length== 0)
			return 0d;
		int medPos= surface.length/ 2;
		if (surface.length% 2== 0&& (medPos+1)< surface.length)
			return ((surface[medPos]+ surface[medPos+ 1])/ 2);
		else
			return surface[medPos];
	}
	
	public double getStandardDeviation() {
		if (surface.length== 0)
			return 0d;
		double med= getMedian();		
		double sum= 0d;
		for (int i = 0; i < surface.length; i++) {
			double val= surface[i]- med;
			val*= val;
			sum+= val;
		}
		
		sum/= (surface.length- 1);
		return Math.sqrt(sum);
	}
}
