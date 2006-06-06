/*
 * Created on May 4, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.algorithm;

import java.util.Hashtable;
import java.util.Vector;

/**
 * 
 * 
 * @author micha
 */
public class ResultHash extends Hashtable {

	protected Vector result= null;
	protected float score= -1;
	protected Vector path= null;
	
	/**
	 * @return
	 */
	public Vector getResult() {
		return result;
	}

	/**
	 * @param vector
	 */
	public void setResult(Vector vector) {
		result= vector;
	}

	/**
	 * @return
	 */
	public Vector getPath() {
		return path;
	}

	/**
	 * @param vector
	 */
	public void setPath(Vector vector) {
		path= vector;
	}

	/**
	 * @return
	 */
	public float getScore() {
		return score;
	}

	/**
	 * @param i
	 */
	public void setScore(float i) {
		score= i;
	}

}
