/*
 * Created on Oct 20, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

/**
 * 
 * 
 * @author micha
 */
public class RepAlignComparator extends PSComparator {

	/**
	 * @param newFileBase
	 * @param newFileExt
	 */
	public RepAlignComparator(String newFileBase, String newFileExt) {
		super(newFileBase, newFileExt);
	}

	public static void main(String[] args) {
		start(DialignComparator.class, ".tfa.ms ");
	}
}
