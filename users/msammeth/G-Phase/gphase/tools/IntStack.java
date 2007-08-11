package gphase.tools;

public class IntStack {
	IntVector vector;
	
	public IntStack() {
		vector= new IntVector();
	}
	
	public boolean isEmpty() {
		return vector.size()== 0;
	}
	
	public int top() {
		return vector.get(vector.size()- 1);
	}
	
	public int pop() {
		return vector.remove(vector.size()- 1);
	}
	
	public void push(int value) {
		vector.add(value);
	}
}
