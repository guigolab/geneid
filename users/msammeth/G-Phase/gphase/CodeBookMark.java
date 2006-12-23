package gphase;
	/**
	 * A class to guess "where am I" at runtime (why?)
	 * using StackTraceElement (JSDK1.4)
	 * 
	 * @author Giorgio Maone
	 * @version 1.0
	 */

public class CodeBookMark {
	 
	 StackTraceElement ste;
	  public CodeBookMark() {
	   try {
	    throw new Throwable();
	  } catch(Throwable t) {
	    ste=t.getStackTrace()[1];
	  }
	  }
	 public String getMethodName() {
	   return ste.getMethodName();
	 }
	 public String getClassName() {
	   return ste.getClassName();
	 }
	 public int getLineNumber() {
	   return ste.getLineNumber();
	 }
	 public String getFileName() {
	   return ste.getFileName();
	 }
	 
	  public static void main(String[] args) {
	    CodeBookMark bookMark = new CodeBookMark();
	  System.out.println(bookMark.getMethodName());
	  }
}
