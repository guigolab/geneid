import java.lang.reflect.Array;

public class Test {

	public static void main(String[] args) {
		Object[] oo= (Object[]) Array.newInstance(Test[].class, 3);
		for (int i = 0; i < oo.length; i++) 
			oo[i]= Array.newInstance(Test.class, 2); 
		
		Test[][] t= (Test[][]) oo;
		for (int i = 0; i < t.length; i++) 
			for (int j = 0; j < t[i].length; j++) 
				System.out.println(i+","+j);
	}
}
