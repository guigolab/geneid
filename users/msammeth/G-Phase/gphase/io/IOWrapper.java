package gphase.io;


/**
 * Interface defining the default IO Wrapper.
 * 
 * @author micha, Thasso 
 */
public interface IOWrapper {

	public void read() throws Exception;
	public void write() throws Exception;
	public boolean isApplicable() throws Exception;
}
