package gphase.graph7;

import java.util.Arrays;
import java.util.Comparator;

public class Combination {
	public static class CombiSorter implements Comparator<long[]> {
		public int compare(long[] o1, long[] o2) {
			for (int i = 0; i < o1.length; i++) {
				if (o1[i]< o2[i])
					return -1;
				else if (o2[i]< o1[i])
					return 1;
			}
			return 0;
		}
	}
	
	static CombiSorter defaultCombiSorter= new CombiSorter();
	
	long[][] combi;
	
	public Combination(long[][] newCombi) {
		this.combi= newCombi;
		Arrays.sort(this.combi, defaultCombiSorter);
	}
	
	
}
