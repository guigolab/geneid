package gphase.tools;

public class IntVector {

	static final int DEFAULT_LOAD_FACTOR= 13;
	int[] vector;
	int size= -1;
	int incrementSize;
	
	public IntVector() {
		this(DEFAULT_LOAD_FACTOR);
	}
	
	public IntVector(int initialCapacity) {
		this(initialCapacity, DEFAULT_LOAD_FACTOR);
	}
	
	public IntVector(int initialCapacity, int loadFactor) {
		vector= new int[initialCapacity];
		incrementSize= loadFactor;
		size= 0;
	}
	
	public void add(int x) {
		if (size== vector.length)
			extendVector();
		vector[size++]= x;
	}
	
	void extendVector() {
		int[] newVector= new int[vector.length+ incrementSize];
		for (int i = 0; i < vector.length; i++) 
			newVector[i]= vector[i];
		vector= newVector;
	}
	
	public int size() {
		return size;
	}
	
	public int remove(int pos) {
		int result= vector[pos];
		--size;
		for (int i = pos; i < size; i++) 
			vector[pos]= vector[pos+1];
		return result;
	}
	
	public int get(int pos) {
		return vector[pos];
	}
	
	public int[] toIntArray() {
		int[] result= new int[size];
		for (int i = 0; i < size; i++) 
			result[i]= vector[i];
		return result;
	}
}
