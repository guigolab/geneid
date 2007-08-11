package gphase.ext;

import java.io.*;
/**  
 * Ein Thread, der eine InputStream ausliest und die erhaltenen Daten ins Nirvana schickt   
 * (bzw. nach /dev/null f??r den Linux-user ;-)). Eigentlich werden sie nur dem GarbageCollector  
 * ?berlassen.  
 * <br><br>  
 * Wird ein Name f??r den Stream ??bergeben (z.B. "standard_out"), so wird der Inhalt des Streams  
 * unter Angabe des Namens nach System.out geschrieben.<br>  
 * Der Thread beendet sich selbst, sobald der stream geschlossen wird oder ein .interrupt()   
 * aufgerufen wird.<br>  
 * Erstellungsdatum: (03.04.01 21:55:05)  
 * @author Tobias Vogele   
 */
public class DevPipeReaderThread extends Thread {

	// needed for run()
	// inherited from interface Constants
	protected boolean DEBUG = false;
	boolean finished= false;

	/**  
	 * The Inputstream, to be read.  
	 */
	protected InputStream in;
	/**  
	 * The secondary Inputstream, to be read.  
	 */
	private InputStream in2;
	/**  
	 * The Outputstream, to be written.  
	 */
	private OutputStream out;
	
	
	/**  
	 * DevNullThread - constructs a new DevPipeReader with the given streams.  
	 */
	public DevPipeReaderThread(InputStream in, OutputStream out) {
		this.in = in;
		this.out= out;
	}
	
	public DevPipeReaderThread() {
	}
	

	/**  
	 * DevNullThread - constructs a new DevPipeReader with the given streams:
	 * both input streams are written to the output stream.
	 */
	public DevPipeReaderThread(InputStream in, InputStream in2, OutputStream out) {
		this(in, out);
		this.in2= in2;
	}

	public InputStream getIn() {
		return in;
	}
	public void setIn(InputStream in) {
		this.in= in;
	}
	public InputStream getIn2() {
		return in2;
	}
	public OutputStream getOut() {
		return out;
	}

	/**  
	 * The main method.  
	 * Reads continously from the input stream(s) as long as there is something to read and
	 * interrupt() is not called
	 */
	public void run() {

		
		int w, x;
		String tst= "";
		finished= false;
		boolean firstFinished= false;
		
		do {
			
				// 
			if (isInterrupted()) {
				try {
					sleep(100);
					System.err.println("Not planned interrupt");
				} catch (InterruptedException e) {
					; //nothing
				}
				continue;
			}
				
				// read from first stream
			if (!firstFinished) {
				try {				
					w = in.read();
//					System.err.print(w+"\n");
					out.write(w);
					tst+= Character.toString(((char) w));	// strange aborting of process output
				} catch (IOException e) {
					w = -1;
				}
				if (w == -1)
					firstFinished= true;
			}

		} while (!firstFinished);

//		System.err.println("finished");		
		try {
			out.flush();
			finished= true;
//			out.close();	do not close System.out
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	public boolean hasFinished() {
		return finished;
	}
	/**
	 * Returns the finished.
	 * @return boolean
	 */
	public boolean isFinished() {
		return finished;
	}

	/**
	 * Sets the finished.
	 * @param finished The finished to set
	 */
	public void setFinished(boolean finished) {
		
		
		this.finished= finished;
	}

}
