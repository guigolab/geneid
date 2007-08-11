package gphase.tools;

import com.sun.mail.iap.Argument;

public class Math {	
	
	public static long pow(int a, int b) {
		if (b== 0)
			return 1;
		if (b== 1)
			return a;
		long c= a;
		for (int i = 0; i < b-1; i++) 
			c*= a;
		return c;
	}
	
	public static long powSafe(int a, int b) throws IllegalArgumentException {
		if (b== 0)
			return 1;
		if (b== 1)
			return a;
		long c= a;
		for (int i = 0; i < b-1; i++) {
			if (c*a> Long.MAX_VALUE)
				throw new IllegalArgumentException("Exceeding long range!");
			c*= a;
		}
		return c;
	}

}
