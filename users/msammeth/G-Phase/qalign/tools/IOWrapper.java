package qalign.tools;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public interface IOWrapper {

	public SequenceWrapper[] getWrappedSequences();
	
	public void setWrappedSequences(SequenceWrapper[] x);
	public void read() throws Exception;
	public void write() throws Exception;	
}
