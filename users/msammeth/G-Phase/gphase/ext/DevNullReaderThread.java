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
public class DevNullReaderThread extends Thread {  

// needed for run()
// inherited from interface Constants
protected boolean DEBUG= false;

/**  
 * Der Inputstream, der ausgelesen wird.  
 */  
private java.io.InputStream in;  
/**  
 * Falls ein Name f?r den Stream angegeben wird, wird er hier gespeichert und mit dem Inhalt  
 * des Streams ausgegeben.  
 */  
private java.lang.String streamName = null;  
/**  
 * DevNullThread - leerer Konstruktor.  
 */  
public DevNullReaderThread() {  
        super("DevNullReader");  
}  
/**  
 * DevNullThread - erstellt einen neunen DevNullReader mit dem ?bergebenen Stream ohne Name.  
 */  
public DevNullReaderThread(InputStream in) {  
        this(in, null);  
}  
/**  
 * DevNullThread - erstellt einen neuen DevNullReader mit den ?bergebenen Stream und Name.  
 */  
public DevNullReaderThread(InputStream in, String name) {  
        super("DevNullReader");  
        this.in = in;  
        streamName = name;  
}  
/**  
 * Gibt den Inputstream zur??ck, von dem gelesen werden soll.<br>  
 * Erstellungsdatum: (03.04.01 21:55:55)  
 * @return java.io.InputStream  
 */  
public java.io.InputStream getIn() {  
        return in;  
}  
/**  
 * Gibt den Name des Inputstreams zur??ck.<br>  
 * Erstellungsdatum: (03.04.01 21:57:24)  
 * @return java.lang.String  
 */  
public java.lang.String getStreamName() {  
        return streamName;  
}  
/**  
 * Die Hauptmethode.  
 * Hier wird immerfort der Stream ausgelesen, solange was da ist zum lesen und solange  
 * keine interrupt() aufgerufen wurde.<br>  
 * Erstellungsdatum: (03.04.01 21:56:13)  
 */  
public void run() {  
                if (DEBUG && streamName != null)
					System.out.println(streamName+"+ >");   

// evtl. std-out und std-err auslesen.  
                int w;  
                do {  
                        try {  
                                w = in.read();  
                        } catch (IOException e) {  
                                w = -1;  
                        }  
                        if (w != -1 && DEBUG && streamName != null)  
System.out.print((char)   
w);  
                } while (w != -1 && ! isInterrupted());  
                if (DEBUG && streamName != null) System.out.println("\n<"); 

  
}  
/**  
 * Setzt den Inputstream.<br>  
 * Erstellungsdatum: (03.04.01 21:55:55)  
 * @param newIn java.io.InputStream  
 */  
public void setIn(java.io.InputStream newIn) {  
        in = newIn;  
}  
/**  
 * Setzt den Name des Streams.<br>  
 * Erstellungsdatum: (03.04.01 21:57:24)  
 * @param newStreamName java.lang.String  
 */  
public void setStreamName(java.lang.String newStreamName) {  
        streamName = newStreamName;  
}  
}  

