/*
 * Created on Nov 27, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import repeat.data.Alignment;
import repeat.data.Block;

/**
 * 
 * 
 * @author micha
 */
public class BAliParserCentre {

	public static final String GAP_CHARS= "-.~";
	
	protected String fileName= null;
	protected Alignment alignment= null;

	public static void main(String[] args) {
		BAliParserCentre myParser= 
			new BAliParserCentre("D:\\pfeiffer\\BAliBASE2.01\\ref6\\test\\kringle_ref6_centre.html");
		Alignment ali= myParser.getAlignment();
		System.out.println(ali);
	}
	
	public BAliParserCentre(String newFileName) {
		
		this.fileName= newFileName;
	}
	
	protected String readFile() {
		
		if (fileName== null)
			return null;
			
		StringBuffer html= new StringBuffer(4000);
		
		try {		
			BufferedReader buffy= new BufferedReader(new FileReader(fileName));
			while (buffy.ready())
				html.append(buffy.readLine()+ "\n");
		} catch (FileNotFoundException e) {
			; // :)
		} catch (IOException e) {
			; // :)
		}
		
		return html.toString();
	}
	
	public Alignment getAlignment() {
		
		if (alignment== null) {
			alignment= decodeFile(readFile());
		}
		
		return alignment;
	}

	protected static Alignment decodeFile(String rawHtml) {
		
		return decodeAlignment(findAlignment(rawHtml));
	}
	
	protected static String findAlignment(String rawHtml) {
		
//		Pattern patti= Pattern.compile(".*<pr>(.*)<pr>.*", Pattern.DOTALL);
		try {
//			Matcher matti= patti.matcher(rawHtml);
//			String res= matti.group(1);
			int pos1= rawHtml.indexOf("<pr>");
			int pos2= rawHtml.indexOf("<pr>", (pos1+ 1));			
			return rawHtml.substring(pos1+ 4, pos2);
		} catch (Exception e) {
			return null;
		}
	}
	
	public static Alignment decodeAlignment(String htmlCode) {
		
		String repeatBlocks= null;
		StringTokenizer toki= new StringTokenizer(htmlCode, "\n");
//		Pattern pattern= Pattern.compile("^(.+)(\\d+)(.+)(\\d+)");
		Pattern pattern= Pattern.compile("(.*)\\s+(\\d+)\\s+(.*)\\s+(\\d+)\\D*");
		Matcher matcher;
		String line;
		int start, end;
		String name= null, seq= null;
		StringBuffer seqPurged= null;
		Integer repeatOffset= null;
		Vector nameVec= new Vector(), seqVec= new Vector(), startVec= new Vector(), endVec= new Vector();
		
		Vector blocksDetected= new Vector();
		Integer aHelix= null, bStrand= null, cBlock= null;
		while (toki.hasMoreTokens()) {
			
				// close eventually open blocks (which ran into border of upper row)
			if (aHelix!= null) {									// close all open Blocks
				Block tmpB= new Block(new String[] {name}, new int[] {aHelix.intValue()}, 
					seqPurged.length()- aHelix.intValue());
				tmpB.setID(Block.ID_AHELIX);
				blocksDetected.add(tmpB);
				aHelix= null;
			}
			if (bStrand!= null) {
				Block tmpB= new Block(new String[] {name}, new int[] {bStrand.intValue()}, 
					seqPurged.length()- bStrand.intValue());
				tmpB.setID(Block.ID_BSTRAND);
				blocksDetected.add(tmpB);
				bStrand= null;
			}
			if (cBlock!= null) {
				Block tmpB= new Block(new String[] {name}, new int[] {cBlock.intValue()}, 
					seqPurged.length()- cBlock.intValue());
				tmpB.setID(Block.ID_CBLOCK);
				blocksDetected.add(tmpB);
				cBlock= null;
			}


			line= new String(toki.nextToken()); 
			matcher= pattern.matcher(line);
			if (!matcher.matches()) {				// try to find repeat blocks
				if ((repeatBlocks== null)&& (line.indexOf("*")!= (-1)))
					repeatBlocks= line;					
				continue;
			}
			
				// match
			name= matcher.group(1).trim();
			nameVec.add(name);
			start= Integer.parseInt(matcher.group(2).trim())- 1;	// convert to 0-based
			startVec.add(new Integer(start));
			seq= matcher.group(3).trim(); 
			end= Integer.parseInt(matcher.group(4).trim())- 1;		// convert to 0-based
			endVec.add(new Integer(end));
			
				// gotta get the f* offset for the f* repeats
			if ((repeatBlocks!= null)&& (repeatOffset== null)) {
				Pattern patti= Pattern.compile("^.*\\s+\\d+(\\s+).*\\s+\\d+\\s+.*");
				Matcher matti= patti.matcher(line);
				if (matti.matches())
					repeatOffset= new Integer(matti.end(1));
			}
			
				// parse sequence
			seqPurged= new StringBuffer();
			int pointer= 0;
			for (;pointer< seq.length();) {

				if (seq.charAt(pointer)!= '<') {
					if (aHelix!= null) {									// close all open Blocks
						Block tmpB= new Block(new String[] {name}, new int[] {aHelix.intValue()}, 
							seqPurged.length()- aHelix.intValue());
						tmpB.setID(Block.ID_AHELIX);
						blocksDetected.add(tmpB);
						aHelix= null;
					}
					if (bStrand!= null) {
						Block tmpB= new Block(new String[] {name}, new int[] {bStrand.intValue()}, 
							seqPurged.length()- bStrand.intValue());
						tmpB.setID(Block.ID_BSTRAND);
						blocksDetected.add(tmpB);
						bStrand= null;
					}
					if (cBlock!= null) {
						Block tmpB= new Block(new String[] {name}, new int[] {cBlock.intValue()}, 
							seqPurged.length()- cBlock.intValue());
						tmpB.setID(Block.ID_CBLOCK);
						blocksDetected.add(tmpB);
						cBlock= null;
					}
					seqPurged.append(seq.charAt(pointer++));	// lazily add
					continue;
				}
				
					// else ...
				int addOffset= 0;					
				if (seq.charAt(pointer+1)== 'u') {	// underlined
					if (cBlock== null)				// open new Block
						cBlock= new Integer(seqPurged.length());
					String tmp= seq.substring(pointer+ 3, seq.indexOf("</u>", pointer));
					addOffset= 8;	// for the char, <u>, </u>
					
					if (tmp.indexOf("font")!= (-1)) {	// colored block spotted
						if (tmp.indexOf("RED")!= (-1)) 
							if (aHelix== null)
								aHelix= new Integer(seqPurged.length());	// open
								
						if (tmp.indexOf("#33cc00")!= (-1))
							if (bStrand== null)
								bStrand= new Integer(seqPurged.length());	// open
						
						int lenB4= tmp.length();
						tmp= tmp.substring(tmp.indexOf(">")+ 1, tmp.indexOf("</font>"));
						addOffset+= lenB4- tmp.length();
					
					} else {	// close open structure blocks
						if (aHelix!= null) {
							Block tmpB= new Block(new String[] {name}, new int[] {aHelix.intValue()}, 
								seqPurged.length()- aHelix.intValue());	// close
							tmpB.setID(Block.ID_AHELIX);
							blocksDetected.add(tmpB);
							aHelix= null;
						}
						if (bStrand!= null) {
							Block tmpB= new Block(new String[] {name}, new int[] {bStrand.intValue()}, 
								seqPurged.length()- bStrand.intValue());	// close
							tmpB.setID(Block.ID_BSTRAND);
							blocksDetected.add(tmpB);
							bStrand= null;
						}
					}
					
						// add character
					seqPurged.append(tmp.trim());
					pointer+= addOffset;
					continue;
					
				} else {						// not underlined
					if (cBlock!= null) {			// close cblock and add to result			
						Block tmp= new Block(new String[] {name}, new int[] {cBlock.intValue()}, 
							seqPurged.length()- cBlock.intValue());
						tmp.setID(Block.ID_CBLOCK);
						blocksDetected.add(tmp);
						cBlock= null;
					}
						
						// must contain <font> if reaching here...
					String tmp= seq.substring(pointer, seq.indexOf("</font>", pointer));
					if (tmp.indexOf("font")!= (-1)) {	// colored block spotted

						addOffset= 8;

						if (tmp.indexOf("RED")!= (-1))
							if (aHelix== null)
								aHelix= new Integer(seqPurged.length());	// open
								
						if (tmp.indexOf("#33cc00")!= (-1))
							if (bStrand== null)
								bStrand= new Integer(seqPurged.length());	// open
						
						int lenB4= tmp.length();
						tmp= tmp.substring(tmp.indexOf(">")+ 1, tmp.length());
						addOffset+= lenB4- tmp.length();
						
						seqPurged.append(tmp.trim());
						pointer+= addOffset;
					
					} else {	// close open structure blocks
						if (aHelix!= null) {
							Block tmpB= new Block(new String[] {name}, new int[] {aHelix.intValue()}, 
								seqPurged.length()- aHelix.intValue()- 1);	// close
							tmpB.setID(Block.ID_AHELIX);
							blocksDetected.add(tmpB);
							aHelix= null;
						}
						if (bStrand!= null) {
							Block tmpB= new Block(new String[] {name}, new int[] {bStrand.intValue()}, 
								seqPurged.length()- bStrand.intValue()- 1);	// close
							tmpB.setID(Block.ID_BSTRAND);
							blocksDetected.add(tmpB);
							bStrand= null;
						}
					}

				}
			} // end for seq parsing
			
			seqVec.add(seqPurged.toString());
		}

			// close all remaning blocks
		if (aHelix!= null) {
			Block tmpB= new Block(new String[] {name}, new int[] {aHelix.intValue()}, 
				seqPurged.length()- aHelix.intValue());	// close
			tmpB.setID(Block.ID_AHELIX);
			blocksDetected.add(tmpB);
			aHelix= null;
		}
		if (bStrand!= null) {
			Block tmpB= new Block(new String[] {name}, new int[] {bStrand.intValue()}, 
				seqPurged.length()- bStrand.intValue());	// close
			tmpB.setID(Block.ID_BSTRAND);
			blocksDetected.add(tmpB);
			bStrand= null;
		}
		if (cBlock!= null) {
			Block tmpB= new Block(new String[] {name}, new int[] {cBlock.intValue()}, 
				seqPurged.length()- cBlock.intValue());	// close
			tmpB.setID(Block.ID_CBLOCK);
			blocksDetected.add(tmpB);
			cBlock= null;
		}
		
			// create alignment
		String[] names= new String[nameVec.size()];
		String[] seqs= new String[seqVec.size()];
		int[] starts= new int[startVec.size()];
		int[] ends= new int[endVec.size()];
		for (int i= 0; i< nameVec.size(); ++i) {
			
			names[i]= (String) nameVec.elementAt(i);
			seqs[i]= (String) seqVec.elementAt(i);
			starts[i]= ((Integer) startVec.elementAt(i)).intValue();
			ends[i]= ((Integer) endVec.elementAt(i)).intValue();
		}
			
			// bug in Balibase: 
			// end position is 1 too high if last char is a gap
		for (int i= 0; i < ends.length; i++) {
			int pos= seqs[i].length()- 1;
			char c= seqs[i].charAt(pos);
			if (GAP_CHARS.indexOf(c)>= 0)
				ends[i]= ends[i]- 1;
		}
		
		Alignment ali= new Alignment();
		ali.setStartPos(starts);
		ali.setEndPos(ends);
		ali.setSequences(seqs);
		ali.setSeqIDs(names);
		
			// close repeat block
		int pointer= 0;
		while (pointer< repeatBlocks.length()) {
			
			while((pointer< repeatBlocks.length())&& (repeatBlocks.charAt(pointer)!= '*'))
					++pointer;	// find first star
			int startRep= pointer- repeatOffset.intValue();
			if (pointer>= seq.length())
				break;
			
			while((pointer< repeatBlocks.length())&& (repeatBlocks.charAt(pointer)== '*'))
				++pointer;	// find end
			
			int[] repStarts= new int[names.length];
			for (int i= 0; i< repStarts.length; ++i)
				repStarts[i]= startRep;
			Block tmpB= new Block(names, repStarts, 
				pointer- startRep- repeatOffset.intValue());	// close
			tmpB.setID(Block.ID_RAREA);
			blocksDetected.add(tmpB);
		}

		for (int i= 0; i < blocksDetected.size(); ++i) 
			ali.addBlock((Block) blocksDetected.elementAt(i));

//		for (int i= 0; i < blocksDetected.size(); i++) {
//			if (((Block) blocksDetected.elementAt(i)).isCBlock())
//				System.out.println(blocksDetected.elementAt(i)+"\n");
//		}
		return ali;
	}
	
}
