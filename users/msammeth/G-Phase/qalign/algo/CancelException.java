package qalign.algo;

/**
 * Instances of this class are thrown if the user aborts the multiple alignment process.
 *
 * @author Michael Sammeth
 * @version 1.0 rc1
 * @since 1.0 rc1
 */

public class CancelException extends Exception {
  public CancelException() {
          super();
  }
  public CancelException(String message) {
          super(message);
  }
}
