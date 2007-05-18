package gphase.tools;

/*
 The class KMP manages a string matching for pattern in text with 
 an algorithm due to Knuth, Morris and Pratt.
 Its attributes are Strings for target and pattern.
 Its methods are a constructor, match(), findPrefix() and showMatch().
 */
public class KMP {
	private String target, pattern;

	public static void main(String[] args) {
		KMP myKMP= new KMP("abCDefAbab", "ab", false);
		System.out.println(myKMP.countMatches());
	}
	/*
	 * Pad index 0 of target and pattern with a blank to force initial character
	 * index of 1.
	 */
	public KMP(String newTarget, String newPattern, boolean caseSensitive) {
		if (!caseSensitive) {
			newTarget= newTarget.toUpperCase();
			newPattern= newPattern.toUpperCase();
		}
		target = " " + newTarget;
		pattern = " " + newPattern;
	}

	/*
	 * Implement Knuth Morris Pratt string matching algorithm.
	 */
	public void match() {
		int n = target.length(), m = pattern.length() - 1, q = 0;
		int[] pi = findPrefix();
	
		System.out.print("\tTarget  " + target.substring(1) + "\n\tPattern "
				+ pattern.substring(1) + "\n\t");
		for (int i = 1; i < n; i++) {
			while ((q > 0) && (pattern.charAt(q + 1) != target.charAt(i)))
				q = pi[q];
			if (pattern.charAt(q + 1) == target.charAt(i))
				q++;
			if (q == m) {
				showMatch(i - m);
				q = pi[q];
			}
		}
		System.out.println();
	}

	/*
	 * Implement Knuth Morris Pratt string matching algorithm.
	 */
	public int countMatches() {
		int n = target.length(), m = pattern.length() - 1, q = 0;
		int[] pi = findPrefix();
	
		int cnt= 0;
		for (int i = 1; i < n; i++) {
			while ((q > 0) && (pattern.charAt(q + 1) != target.charAt(i)))
				q = pi[q];
			if (pattern.charAt(q + 1) == target.charAt(i))
				q++;
			if (q == m) {
				++cnt;
				q = pi[q];
			}
		}
		
		//System.out.println(cnt+" "+ pattern+"in\n"+target);
		return cnt;
	}
	/*
	 * Compute prefix function for pattern.
	 */
	private int[] findPrefix() {
		int m = pattern.length(), k = 0;
		int[] pi = new int[m];

		pi[0] = pi[1] = 0;
		for (int q = 2; q < m; q++) {
			while ((k > 0) && (pattern.charAt(k + 1) != pattern.charAt(q)))
				k = pi[k];
			if (pattern.charAt(k + 1) == pattern.charAt(q))
				k++;
			pi[q] = k;
		}

		return pi;
	}

	/*
	 * Display an occurrence of pattern at index i in text.
	 */
	private void showMatch(int i) {
		System.out.println("Match at ndx " + i);
		for (int k = 0; k < i; k++)
			System.out.print(" ");
		System.out.println(pattern + "\n" + target);
	}
}