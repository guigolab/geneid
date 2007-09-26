/*
 * Created on Nov 26, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import repeat.data.Repeat;

/**
 * 
 * 
 * @author micha
 */
public class BAliParserSchema {

	protected String fileName= null;
	protected Repeat firstRepeat= null; 
	protected Repeat[] repeats= null;
	
	public static void main(String[] args) {
		
		BAliParserSchema testParser= 
			new BAliParserSchema("D:\\pfeiffer\\BAliBASE2.01\\ref6\\test_4\\sh3_ref6_schema.html");
		Repeat result= testParser.getFirstRepeat();
		while(result!= null) {
			System.out.println(result);
			result= result.getNext();
		}
	}
	
	public BAliParserSchema(String newFileName) {
		
		this.fileName= newFileName;
	}
	
	public Repeat getFirstRepeat() {
		
		if (firstRepeat== null) {
			firstRepeat= decodeFile(readFile());
			
//			Repeat tmp= firstRepeat;
//			while (tmp!= null) {
//				System.out.println(tmp);
//				tmp= tmp.getNext();
//			}
		}
		
		return firstRepeat;
	}
	
	public Repeat[] getRepeats() {
		
		if (repeats == null) {
			if (getFirstRepeat()== null)
				return null;
				
			Repeat tmpRepeat= firstRepeat;			// count
			int counter= 0;
			while (tmpRepeat!= null) {
				++counter;
				tmpRepeat= tmpRepeat.getNext();
			}
			
			repeats= new Repeat[counter];			// copy
			tmpRepeat= firstRepeat;
			counter= 0;
			while (tmpRepeat!= null) {
				repeats[counter++]= tmpRepeat;
				tmpRepeat= tmpRepeat.getNext();
			}
		}

		return repeats;
	}
	
	protected String readFile() {
		
		if (fileName== null)
			return null;
			
		StringBuffer html= new StringBuffer(4000);
		
		try {		
			BufferedReader buffy= new BufferedReader(new FileReader(fileName));
			while (buffy.ready())
				html.append(buffy.readLine());
		} catch (FileNotFoundException e) {
			System.err.println(e); // :)
		} catch (IOException e) {
			System.err.println(e); // :)
		}
		
		return html.toString();
	}
	
	protected static Repeat decodeFile(String rawHtml) {
		
		return decodeTable(findTable(rawHtml));
	}
	
	protected static String findTable(String rawHtml) {
		
		int pos1= 0;
		int pos2= rawHtml.indexOf("<table");
		String subString= null;
		while(pos2!= (-1)) {
			 subString= rawHtml.substring(pos1, pos2);
			 if ((subString.indexOf("Sequence Name")!= (-1))&&
			 		(subString.indexOf("# Repeats")!= (-1)))
			 	break;
			 
			 pos1= pos2;	// find next
			 pos2= rawHtml.indexOf("<table", pos1+ 1);
			 subString= null;
		}
		
		return subString;
	}
	
	/* now with two auto-switches for markers, hope now its ok */
	public static Repeat decodeTable(String table) {
		
		int startCounter= -1, lengthCounter= -1, IDCounter= -1;
		boolean continu= true;
		
			// fill in start, length
		String posMarker= "<font size=2> <pre>";		// group 4c, or group test: "<font size=2> <pre>";
														//	group1: "<pre><font size=\"2\">" 		

		int pos= table.indexOf(posMarker);
		if (pos== -1) {
			posMarker= "<font size=2> <pre>";			// group 4c, or group test
			pos= table.indexOf(posMarker);
		}
		Repeat repeat= new Repeat();
		String seqID= null;
		while(pos!= (-1)) {
			repeat.setStart(Integer.parseInt(table.substring(pos+ posMarker.length(),	 
				table.indexOf("</font>", pos+ posMarker.length())).trim())- 1);					// start 0-based
				
			pos= table.indexOf(posMarker, pos+ 1);
			repeat.setLength(Integer.parseInt(table.substring(pos+ posMarker.length(), 
				table.indexOf("</font>", pos+ posMarker.length())).trim())- repeat.getStart());	// last pos is INCLUDED
				
				// prepare next
			repeat.setNext(new Repeat());
			repeat.getNext().setPrev(repeat);
			repeat= repeat.getNext();
			pos= table.indexOf(posMarker, pos+1);
		}
		repeat.getPrev().setNext(null);
		
			// fill in ID, seqID
		while (repeat.getPrev()!= null)
			repeat= repeat.getPrev();
		String nameMarker= "<td><pre>";

		posMarker= "colspan=2 align=center>";	// group4c: "colspan=\"2\" align=\"center\">"
												// for group1: "colspan=2 align=center>";
		if (table.indexOf(posMarker)== -1)
			posMarker= "colspan=\"2\" align=\"center\">";	// group4c
		pos= 0;
		for (int i= 0; i< 3; ++i)
			pos= table.indexOf(nameMarker, pos+ 1);
		String currentName= table.substring(pos+ nameMarker.length(), table.indexOf("/pre></td>", pos)).trim();
		currentName= currentName.substring(0, currentName.length()- 1);
		pos= table.indexOf(posMarker, pos);			
		while(pos!= (-1)) {
			String tmpType= table.substring(pos+ posMarker.length(), 
								table.indexOf("</td>", pos+ posMarker.length())).trim();
			while (tmpType.indexOf('<')>= 0) {		// purge
				if (tmpType.indexOf('<')== 0) 		// tag at start
					tmpType= tmpType.substring(
						tmpType.indexOf('>')+ 1,
						tmpType.length()
					).trim();
				else								// tag at end
					tmpType= tmpType.substring(
						0,
						tmpType.indexOf('<')
					).trim();
			}
			
			if (tmpType.trim().length()== 0) {
				tmpType= "NOID";		// test/lrr_ref6
				System.err.println("No ID for "+currentName+"! assinged NOID");	
			}
			repeat.setType(tmpType);
			repeat.setSeqName(currentName);
			
			int tmp= table.indexOf(nameMarker, pos+1);
			if ((tmp!= (-1))&& (tmp< table.indexOf(posMarker, pos+1))) {
				currentName= table.substring(tmp+ nameMarker.length(), table.indexOf("/pre></td>", tmp+1)).trim();
				currentName= currentName.substring(0, currentName.length()- 1);		// kill trailing '>'				
			}
										
			pos= table.indexOf(posMarker, pos+1);
			if (pos!= (-1))
				repeat= repeat.getNext();
		}
		
		while (repeat.getPrev()!= null)
			repeat= repeat.getPrev();
		return repeat;
	}
}
