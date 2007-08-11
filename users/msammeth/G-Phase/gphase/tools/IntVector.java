package gphase.tools;

public class IntVector {

	static final int DEFAULT_LOAD_FACTOR= 13;
	int[] vector;
	public int length= -1;
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
		length= 0;
	}
	
	public void add(int x) {
		if (length== vector.length)
			extendVector();
		vector[length++]= x;
	}
	
	@Override
	protected Object clone() throws CloneNotSupportedException {
		IntVector newV= new IntVector();
		newV.vector= this.vector.clone();
		newV.length= this.length;
		newV.incrementSize= this.incrementSize;
		return newV;
	}
	
	public IntVector cloneIntVector() {
		try {
			return (IntVector) clone();
		} catch (CloneNotSupportedException e) {
			e.printStackTrace();
		}
		return null;
	}
	
	void extendVector() {
		int[] newVector= new int[vector.length+ incrementSize];
		for (int i = 0; i < vector.length; i++) 
			newVector[i]= vector[i];
		vector= newVector;
	}
	
	public int size() {
		return length;
	}
	
	public int remove(int pos) {
		int result= vector[pos];
		--length;
		for (int i = pos; i < length; i++) 
			vector[pos]= vector[pos+1];
		return result;
	}
	
	public int get(int pos) {
		return vector[pos];
	}
	
	public void set(int pos, int val) {
		if (pos< vector.length)
			vector[pos]= val;
	}
	public int[] toIntArray() {
		int[] result= new int[length];
		for (int i = 0; i < length; i++) 
			result[i]= vector[i];
		return result;
	}

	public void insert(int val, int p) {
		
		if (p< 0)
			p= (p+1)* (-1);
		
		int[] newA= new int[vector.length+ 1];
		for (int i = 0; i < p; i++) 
			newA[i]= vector[i];
		newA[p]= val;
		for (int i = p+1; i < newA.length; i++) 
			newA[i]= vector[i-1];
		
		vector= newA;
		length= newA.length;
	}
	
	public int getValue(int pos) {
		if (pos< length)
			return vector[pos];
		return 0;
	}
	
	public void putValue(int pos, int val) {
		if (pos>= length) {
			int rounds= (pos+ 1)- size();
			for (int i = 0; i < rounds; i++) 
				add(0);
		}
		vector[pos]= val;
	}
}
