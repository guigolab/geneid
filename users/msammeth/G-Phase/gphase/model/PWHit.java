package gphase.model;

import java.io.Serializable;

/**
 * 
 * 
 * @author micha
 */
public class PWHit implements Serializable {
	
	public static final int TYPE_NOINIT= -1;
	public static final int TYPE_NOBRH= 0;
	public static final int TYPE_MBRH= 1;
	public static final int TYPE_UBRH= 2;
//	public static final int TYPE_SYNTHENY= 4;
	
	Object o1, o2;
	String[] alignment= null;
	int score= -1;
	double percentID= -1d;
	int type= TYPE_NOINIT;
	
	public PWHit(Object o1, Object o2) {
		setObjects(o1, o2);
	}
	public void setObjects(Object o1, Object o2) {
		setObject1(o1);
		setObject2(o2);
	}
	public String[] getAlignment() {
		return alignment;
	}
	public void setAlignment(String[] alignment) {
		this.alignment = alignment;
	}
	public int getScore() {
		return score;
	}
	public void setScore(int cost) {
		this.score = cost;
	}
	public Object getObject1() {
		return o1;
	}
	public void setObject1(Object o1) {
		this.o1 = o1;
	}
	public Object getObject2() {
		return o2;
	}
	public void setObject2(Object o2) {
		this.o2 = o2;
	}
	public Object getOtherObject(Object o) {
		if (o== o1)
			return o2;
		if (o== o2)
			return o1;
		return null;	// not found
	}
	public int getType() {
		return type;
	}
	public void setType(int type) {
		this.type = type;
	}
	public void addType(int type) {
		this.type|= type;
	}
	public boolean isType(int type) {
		return ((this.type& type)== type);
	}
	public double getPercentID() {
		return percentID;
	}
	public void setPercentID(double percentID) {
		this.percentID = percentID;
	}
}
