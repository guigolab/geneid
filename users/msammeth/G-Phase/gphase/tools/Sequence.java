package gphase.tools;

public class Sequence {

	public static int countGC(String seq) {
		seq= seq.toUpperCase();
		int cnt= 0;
		for (int i = 0; i < seq.length(); i++) {
			char c= seq.charAt(i);
			if (c== 'G'|| c== 'C')
				++cnt;
		}
		return cnt;
	}
}
