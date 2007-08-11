package gphase.tools;

public class Formatter {
	public static String fprint(double fp, int dec) {
		String s= java.lang.Double.toString(fp); 
		int p= s.lastIndexOf(".");
		if (p< 0) {
			s+= ".";
			for (int i = 0; i < dec; i++) 
				s+= "0";
		} else {
			int q= s.indexOf("E");
			String exp= "";
			if (q>= 0)
				exp= s.substring(q);
			int end= p+ dec+ 1;
			if (end< s.length())
				s= s.substring(0, end);
			else
				for (int i = s.length(); i < end; i++) 
					s+= "0";
			s+= exp;
		}
		
		
		return s;
	}
}
