/*
 * Created on Oct 25, 2004
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Method;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import repeat.data.Repeat;

/**
 * 
 * 
 * @author micha
 */
public class RadekRepeatPoundWrapper extends RepeatPoundWrapper {
	
	protected String seqName= null;
	protected String type= null;
	
	public static void conquerBAliBase(String ext) {
		
		Class[] pars= {String.class, String.class};			// get redirection method
		Method m= null;
		Constructor c= null;
		try {
			c= RadekRepeatPoundWrapper.class.getConstructor(pars);
			m= RadekRepeatPoundWrapper.class.getMethod("convert",null);
		} catch (NoSuchMethodException e) {
			; // :)
		}
		
		BaliConqueror conqueror= new BaliConqueror(c,m,ext);	// iterate BaliBase
		conqueror.conquer();
	}
	
	public void convert() {
		
		setRepeats(getRepeats());
		fileBase= fileBase.substring(0, fileBase.lastIndexOf('.'))+ ".ttt";
		write();
	}
	
	
	protected Repeat[] read() {
		
			// create reader
		BufferedReader buffy= null;
		try {
			buffy= new BufferedReader(new FileReader(fileBase));
		} catch (FileNotFoundException e) {System.err.println(e);;}	// :)
		
			// read in
		Vector repVec= new Vector();
		try {
			while (buffy.ready()) {

					// find start of sequence or tor
				String line= buffy.readLine();
				while(buffy.ready() && 		// look for new protein or repeat type
						(line.indexOf(MARKER_TOR)< 0 
						&& line.indexOf(MARKER_PROTEIN_SEPARATOR)< 0)
						&& line.indexOf(">")< 0)
					line= buffy.readLine();
				if (!buffy.ready())
					break;
				if(line.indexOf(MARKER_PROTEIN_SEPARATOR)>= 0) {		// eop -> >newFasta
					while (buffy.ready()&& line.indexOf(">")< 0)
						line= buffy.readLine();				
				}
				if (!buffy.ready())
					break;
				if(line.indexOf(MARKER_TOR)< 0) {			// >newFasta -> TOR / EOP
					seqName= line.substring(1, line.length()).trim();
					while (buffy.ready()&& line.indexOf(MARKER_TOR)< 0
										&& line.indexOf(MARKER_EOP)< 0)
						line= buffy.readLine();
				}
				if (!buffy.ready())
					break;
				if (line.indexOf(MARKER_EOP)>= 0) {
					String confirmEnd= line.substring(line.indexOf(MARKER_EOP)+ MARKER_EOP.length(), 
												line.length()).trim();
					if (!confirmEnd.equals(seqName))
						System.err.println("Error in file "+fileBase+": partially processed @ "+seqName);
					continue;
				}
								
				type= buffy.readLine().trim();	// else: (line.indexOf(MARKER_TOR)>= 0)
		
				Vector subVec= readRepeatType(buffy);
				for (int i= 0; (subVec!= null)&& i < subVec.size(); i++) 
					repVec.add(subVec.elementAt(i));
			}
		} catch (IOException e) {;}	// :)
	
		
			// convert
		Repeat[] rep= new Repeat[repVec.size()];
		for (int i= 0; i < rep.length; i++) 
			rep[i]= (Repeat) repVec.elementAt(i);
			
		return rep;
	}

	public static final String MARKER_TOR= "type_of_the_repeat";
	public static final String MARKER_EOP= "end of the protein";
	public static final String MARKER_PROTEIN_SEPARATOR= "//";
	
	public static final String MARKER_LIST= "START LENGTH [PVALUE [SCORE]]";
	public static final String MARKER_ALI= "Profile pattern, \"X\": profile column, \"-\": a gap";

	
	/**
	 * @param newFileBase
	 */
	public RadekRepeatPoundWrapper(String newFileBase) {
		super(newFileBase);
		// TODO Auto-generated constructor stub
	}
	
		// compatibility with BaliConqueror
	public RadekRepeatPoundWrapper(String newFileBase, String newFileExt) {
		this(newFileBase+ newFileExt);
	}

	/**
	 * @param newFileBase
	 * @param newRepeats
	 */
	public RadekRepeatPoundWrapper(String newFileBase, Repeat[] newRepeats) {
		super(newFileBase, newRepeats);
		// TODO Auto-generated constructor stub
	}

	public static void main(String[] args) {
		
			// redirect to bali conqueror
		conquerBAliBase(".trs");		
		if (1== 1)
			System.exit(0);			// cutoff normal main


		RadekRepeatPoundWrapper rrepWrap= 
			new RadekRepeatPoundWrapper("D:\\workspace\\Trust\\test4_ank_ref6.tst");
//		Repeat[] repeats= rrepWrap.getRepeats();
//		for (int i= 0; i < repeats.length; i++) 
//			System.out.println(repeats[i]);
		RepeatPoundWrapper repWrap= new RepeatPoundWrapper("D:\\workspace\\Trust\\test4_ank_ref6.ttt");
		repWrap.setRepeats(rrepWrap.getRepeats());
		repWrap.write();
	}

	protected Vector readRepeatType(BufferedReader buffy) throws IOException {

		
			// list 
		String marker= MARKER_POUND+" "+MARKER_LIST;
		String line= buffy.readLine();
		while(line.indexOf(marker)< 0)
			line= buffy.readLine();
		line= buffy.readLine();
		Pattern patty= Pattern.compile("^(\\d+)\\s+(\\d+)\\s+#\\s+\\w+\\s+(\\d+).*");	
		Matcher matty= patty.matcher(line);
		Vector repVec= new Vector();
		while (matty.matches()) {
			
			Repeat tmpRep= new Repeat();
			tmpRep.setSeqName(seqName);
			tmpRep.setStart(Integer.parseInt(matty.group(1)));
			tmpRep.setLength(Integer.parseInt(matty.group(2)));
			tmpRep.setNb(Integer.parseInt(matty.group(3)));
			tmpRep.setType(type);
			repVec.add(tmpRep);

			line= buffy.readLine();			// next line
			matty= patty.matcher(line);
		}
		
			// check for alignment
		while (buffy.ready()) {
			if (line.indexOf(MARKER_ALI)> 0) 		// alignment
				break;
			line= buffy.readLine();
		}
				
		line= buffy.readLine();				// profile
		int counter= repVec.size();
		while(counter> 0) {
			line= buffy.readLine().trim();			// > Name
//			patty= Pattern.compile(".*(\\d+)$");	
//			matty= patty.matcher(line);
			line= line.substring(line.lastIndexOf(" ")+ 1, line.length()).trim();
			int repNb= Integer.parseInt(line);
			for (int i= 0; i < repVec.size(); i++) 
				if (((Repeat) repVec.elementAt(i)).getNb()== repNb) {
					((Repeat) repVec.elementAt(i)).setAlignedSeq(buffy.readLine());
					((Repeat) repVec.elementAt(i)).setAlignmentID(seqName+":"+type);		// redundant,ok?
				}
			--counter;
		}
		
		return repVec;
	}
}
