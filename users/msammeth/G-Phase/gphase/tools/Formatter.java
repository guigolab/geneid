package gphase.tools;

public class Formatter {
	public static String fprint(double fp, int dec) {
		String s= Double.toString(fp);
		int p= s.lastIndexOf(".");
		if (p< 0) {
			s+= ".";
			for (int i = 0; i < dec; i++) 
				s+= "0";
		} else {
			int end= p+ dec+ 1;
			if (end< s.length())
				s= s.substring(0, end);
			else
				for (int i = s.length(); i < end; i++) 
					s+= "0";
		}
		
		return s;
	}
}
