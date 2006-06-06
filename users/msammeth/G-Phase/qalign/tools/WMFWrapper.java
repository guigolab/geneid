package qalign.tools;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Vector;

import sun.misc.HexDumpEncoder;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class WMFWrapper {
	
	public final static int FW_DONTCARE= 0;
	public final static int FW_THIN= 100;
	public final static int FW_EXTRALIGHT= 200;
	public final static int FW_ULTRALIGHT= FW_EXTRALIGHT;
	public final static int FW_LIGHT= 300;
	public final static int FW_NORMAL= 400;
	public final static int FW_REGULAR= FW_NORMAL;
	public final static int FW_MEDIUM= 500;
	public final static int FW_SEMIBOLD= 600;
	public final static int FW_DEMIBOLD= FW_SEMIBOLD;
	public final static int FW_BOLD= 700;
	public final static int FW_EXTRABOLD= 800;
	public final static int FW_ULTRABOLD= FW_EXTRABOLD;
	public final static int FW_HEAVY= 900;
	public final static int FW_BLACK= FW_HEAVY;

	public final static int WMF_SHDR_KEY= 0x9AC6CDD7;
	public final static int META_SETBKCOLOR= 0x0201; //
	public final static int META_SETBKMODE= 0x0102; //
	public final static int META_SETMAPMODE= 0x0103; //
	public final static int META_SETROP2= 0x0104; //
	public final static int META_SETRELABS= 0x0105;
	public final static int META_SETPOLYFILLMODE= 0x0106; //
	public final static int META_SETSTRETCHBLTMODE= 0x0107;
	public final static int META_SETTEXTCHAREXTRA= 0x0108;
	public final static int META_SETTEXTCOLOR= 0x0209; //
	public final static int META_SETTEXTJUSTIFICATION= 0x020A; //
	public final static int META_SETWINDOWORG= 0x020B; //
	public final static int META_SETWINDOWEXT= 0x020C; //
	public final static int META_SETVIEWPORTORG= 0x020D;
	public final static int META_SETVIEWPORTEXT= 0x020E;
	public final static int META_OFFSETWINDOWORG= 0x020F;
	public final static int META_SCALEWINDOWEXT= 0x0410;
	public final static int META_OFFSETVIEWPORTORG= 0x0211;
	public final static int META_SCALEVIEWPORTEXT= 0x0412;
	public final static int META_LINETO= 0x0213; //
	public final static int META_MOVETO= 0x0214; //
	public final static int META_EXCLUDECLIPRECT= 0x0415;
	public final static int META_INTERSECTCLIPRECT= 0x0416;
	public final static int META_ARC= 0x0817; //
	public final static int META_ELLIPSE= 0x0418; //
	public final static int META_FLOODFILL= 0x0419;
	public final static int META_PIE= 0x081A; //
	public final static int META_RECTANGLE= 0x041B; //
	public final static int META_ROUNDRECT= 0x061C; //
	public final static int META_PATBLT= 0x061D; //
	public final static int META_SAVEDC= 0x001E;
	public final static int META_SETPIXEL= 0x041F; //
	public final static int META_OFFSETCLIPRGN= 0x0220;
	public final static int META_TEXTOUT= 0x0521; //
	public final static int META_BITBLT= 0x0922;
	public final static int META_STRETCHBLT= 0x0B23; //
	public final static int META_POLYGON= 0x0324; //
	public final static int META_POLYLINE= 0x0325; //
	public final static int META_ESCAPE= 0x0626; //
	public final static int META_RESTOREDC= 0x0127; //
	public final static int META_FILLREGION= 0x0228;
	public final static int META_FRAMEREGION= 0x0429;
	public final static int META_INVERTREGION= 0x012A;
	public final static int META_PAINTREGION= 0x012B;
	public final static int META_SELECTCLIPREGION= 0x012C; //
	public final static int META_SELECTOBJECT= 0x012D; //
	public final static int META_SETTEXTALIGN= 0x012E; //
	public final static int META_DRAWTEXT= 0x062F;
	public final static int META_CHORD= 0x0830; //
	public final static int META_SETMAPPERFLAGS= 0x0231;
	public final static int META_EXTTEXTOUT= 0x0A32; //
	public final static int META_SETDIBTODEV= 0x0D33;
	public final static int META_SELECTPALETTE= 0x0234; //
	public final static int META_REALIZEPALETTE= 0x0035; //
	public final static int META_ANIMATEPALETTE= 0x0436;
	public final static int META_SETPALENTRIES= 0x0037;
	public final static int META_POLYPOLYGON= 0x0538; //
	public final static int META_RESIZEPALETTE= 0x0139;
	public final static int META_DIBBITBLT= 0x0940; //
	public final static int META_DIBSTRETCHBLT= 0x0B41; //
	public final static int META_DIBCREATEPATTERNBRUSH= 0x0142; //
	public final static int META_STRETCHDIB= 0x0F43; //
	public final static int META_EXTFLOODFILL= 0x0548;
	public final static int META_RESETDC= 0x014C;
	public final static int META_STARTDOC= 0x014D;
	public final static int META_STARTPAGE= 0x004F;
	public final static int META_ENDPAGE= 0x0050;
	public final static int META_ABORTDOC= 0x0052;
	public final static int META_ENDDOC= 0x005E;
	public final static int META_DELETEOBJECT= 0x01F0; //
	public final static int META_CREATEPALETTE= 0x00F7; //
	public final static int META_CREATEBRUSH= 0x00F8;
	public final static int META_CREATEPATTERNBRUSH= 0x01F9;
	public final static int META_CREATEPENINDIRECT= 0x02FA; //
	public final static int META_CREATEFONTINDIRECT= 0x02FB; //
	public final static int META_CREATEBRUSHINDIRECT= 0x02FC; //
	public final static int META_CREATEBITMAPINDIRECT= 0x02FD;
	public final static int META_CREATEBITMAP= 0x06FE;
	public final static int META_CREATEREGION= 0x06FF; //

	public final static int MFCOMMENT= 15;
	public final static int SRCCOPY= 0xCC0020;
	public final static int PATCOPY= 0xF00021;
	public final static int PATINVERT= 0x5A0049;
	public final static int DSTINVERT= 0x550009;
	public final static int BLACKNESS= 0x000042;
	public final static int WHITENESS= 0xFF0062;
	public final static int BI_RLE8= 1;
	public final static int BI_RLE4= 2;

	public final static int TA_BASELINE= 24; // TextAlign options
	public final static int TA_BOTTOM= 8;
	public final static int TA_CENTER= 6;
	public final static int TA_UPDATECP= 1; // FIXME: update current postion
	public final static int TA_TOP= 0;
	public final static int OPAQUE= 2;
	public final static int TRANSPARENT= 1;
	public final static int ETO_GRAYED= 1;
	public final static int ETO_OPAQUE= 2;
	public final static int ETO_CLIPPED= 4;
	public final static int PS_SOLID= 0;
	public final static int PS_DASH= 1;
	public final static int PS_DOT= 2;
	public final static int PS_DASHDOT= 3;
	public final static int PS_DASHDOTDOT= 4;
	public final static int PS_NULL= 5;
	public final static int PS_INSIDEFRAME= 6;

	public static void main(String[] args) {

		WMFWrapper my= new WMFWrapper("test.wmf");
		my.setImgBounds(new Rectangle(0,0,500,500));
//		my.setColor(Color.black);
		my.drawLine(3,3,400,400);
		my.setFont(new Font("Arial", Font.PLAIN, 12), 30, -30, true, true);
//		my.selectObject();	// must be in there for not missing font
//		my.setBKMode(TRANSPARENT);
//		my.setTextColor(Color.black);
//		my.setTextAlign();
		my.drawString("Hallo", 20, 20);
		my.setFont(new Font("Script", Font.PLAIN, 24), 0, 0, false, false);
		my.drawString("Welt", 50, 50);
//		my.drawStringExt("Welt", 30, 30, 0, null);
		my.deleteObject(my.objCounter- 1);
		my.setTerminator();
		my.writeOut();
	}
	
	
	/**
	 * @author micha
	 */
	public class WmfRecord {
	
		protected int size= -1;
		protected int function= -1;
		protected Vector parameters= new Vector();
	
		public WmfRecord(int func) {
			this.function= func;
		}
	
		public int getSize() {
		
			return (3+ parameters.size());
		}
	
		public int[] getWords() {
			
			int[] par= new int[parameters.size()+ 3];
			int size= getSize();
			par[0]= size;
			par[1]= size >> 16;
			par[2]= function;
			for (int i= 0; i< parameters.size(); ++i)
				par[i+3]= ((Integer) parameters.elementAt(i)).intValue();
		
			return par;
		}
	
		public void addParameter(int par) {
		
			parameters.add(new Integer(par));
		}
	}

	public WMFWrapper(String fName) {
		
		this.fileName= fName;
	}

	protected String fileName= null;
		
	/**
	 * Coordinates of upper-left and lower-right corner in twips (metafile units).
	 */
	protected Rectangle imgRect= null;

	/**
	 * Number of twips (metafile units) per inch used to scale the image:<br>
	 * Normally, there are 1440 twips per inch (1:1 scale ratio), which may be changed
	 * to 720 (2:1 double size), 360 (4:1), ...
	 */
	protected int imgTPI= 1440;
	
	protected int scrRes= Toolkit.getDefaultToolkit().getScreenResolution();
	protected Vector wmfRecords= new Vector();
	protected int objCounter= 0;
	
	protected void write(OutputStream out, int data, int bytes) throws IOException {
		
		for (int i= 0; i< bytes; ++i) 
			out.write(data >> (i* 8));
	}
	
	protected void writeWmfSpecialHeader(OutputStream out) throws IOException {

		int chkSum= 0;
		write(out, WMF_SHDR_KEY, 4);	// key
		chkSum^= (WMF_SHDR_KEY& 0x0000FFFF);
		chkSum^= ((WMF_SHDR_KEY& 0xFFFF0000)>> 16);
		write(out, 0, 2);				// handle
		chkSum^= (0& 0x0000FFFF);
		
		write(out, imgRect.x, 2);			// left
		chkSum^= (imgRect.x& 0x0000FFFF);
		write(out, imgRect.y, 2);			// top
		chkSum^= (imgRect.y& 0x0000FFFF);
		int imgRight= imgRect.x+ imgRect.width;
		write(out, imgRight, 2);		// right
		chkSum^= (imgRight& 0x0000FFFF);
		int imgBottom= imgRect.y+ imgRect.height;
		write(out, imgBottom, 2);		// bottom
		chkSum^= (imgBottom& 0x0000FFFF);

		write(out, imgTPI, 2);			// twips per inch, scale
		chkSum^= (imgTPI& 0x0000FFFF);
		
		write(out, 0, 4);				// reserved
		chkSum^= (0& 0x0000FFFF);
		write(out, chkSum, 2);			// checksum
	}
	
	public void writeOut() {
		
		try {
			
			FileOutputStream f= new FileOutputStream(fileName);
			writeWmfHeader(f);
			
			writeRecords(f);
			
			f.flush();
			f.close();
		} catch (Exception e) {
			System.err.println(e);
			e.printStackTrace();
		}
	}
	
	public void writeRecords(OutputStream out) throws IOException {
		
		int[] tmpRec;
		for (int i= 0; i< wmfRecords.size(); ++i) {
			
			tmpRec= ((WmfRecord) wmfRecords.elementAt(i)).getWords();
			for (int j= 0; j< tmpRec.length; ++j)
				write(out, tmpRec[j], 2);
		}
	}
	
	protected void writeWmfHeader(OutputStream out) throws IOException {
		
			// the prepending header for portable wmf
		writeWmfSpecialHeader(out);
		
			// wmf-header
		write(out, 1, 2);				// file type (1=memory, 2=disk), but memory normal...
		write(out, 9, 2);				// header size (always 9 words) 
		write(out, 0x300, 2);			// version of MS windows used (Windows 3.0..)

		write(out, getSize(), 4);			// total file size (in words)
		write(out, (objCounter+ 1), 2);	// number of objects in file
		write(out, getMaxRecordSize(), 4);	// size of largest record (in words)
		write(out, 0, 2);					// no parameters (not used, always 0)
	}
	
	public int getSize() {
		
		int size= 
			11		// special header
			+ 9
			+ getRecordsSize();	// header
		return size;
	}
	
	public int getRecordsSize() {
		
		int sz= 0;
		for (int i= 0; i< wmfRecords.size(); ++i) 
			sz+= ((WmfRecord) wmfRecords.elementAt(i)).getSize();
		
		return sz;
	}
	
	public int getMaxRecordSize() {
		
		int sz= 0;
		int max= 0;
		for (int i= 0; i< wmfRecords.size(); ++i) {
			sz= ((WmfRecord) wmfRecords.elementAt(i)).getSize();
			if (sz> max)
				max= sz;
		}
		
		return max;
	}	
	
	public int getNumOfObjects() {
		
		return (objCounter+ 1);
	}
	

	
	public void setImgBounds(Rectangle rect) {
		
//		System.out.println("set img bounds");
		rect.x= correctValue(rect.x);
		rect.y= correctValue(rect.y);
		rect.width= correctValue(rect.width);
		rect.height= correctValue(rect.height);
		
		this.imgRect= rect;

		WmfRecord tmpRec;
		tmpRec= new WmfRecord(META_SETWINDOWORG);
		tmpRec.addParameter(rect.x);
		tmpRec.addParameter(rect.y);
		wmfRecords.add(tmpRec);

		tmpRec= new WmfRecord(META_SETWINDOWEXT);
		tmpRec.addParameter(rect.height);
		tmpRec.addParameter(rect.width);
		wmfRecords.add(tmpRec);
	}
	
	public void setImgTPI(int newTPI) {
		
		this.imgTPI= newTPI;
	}
	
	public void drawLine(int x1, int y1, int x2, int y2) {
		
//		System.out.println("drw line");
		x1= correctValue(x1);
		y1= correctValue(y1);
		x2= correctValue(x2);
		y2= correctValue(y2);

		WmfRecord rec;
		rec= new WmfRecord(META_MOVETO);
		rec.addParameter(y1);
		rec.addParameter(x1);
		wmfRecords.add(rec);

		rec= new WmfRecord(META_LINETO);
		rec.addParameter(y2);
		rec.addParameter(x2);
		wmfRecords.add(rec);
	}
	
	public void setColor(Color c) {
		
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_CREATEPENINDIRECT);
		tmpRecord.addParameter(PS_SOLID);
		tmpRecord.addParameter(0);
		tmpRecord.addParameter(0);
		tmpRecord.addParameter((c.getGreen() >> 8)+ c.getRed());
		tmpRecord.addParameter(c.getBlue());
		wmfRecords.add(tmpRecord);
		
/*		tmpRecord= new WmfRecord(META_SELECTOBJECT);
		tmpRecord.addParameter(0);
		wmfRecords.add(tmpRecord);
*/		
	}
	
	public void setTextColor(Color c) {
		
//		System.out.println("set tcol");
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_SETTEXTCOLOR);
		tmpRecord.addParameter(0);
		tmpRecord.addParameter(0);
		wmfRecords.add(tmpRecord);
	}
	
	/**
	 * see msdn-library logfont
	 */
	public void setFont(Font font, int esc, int ori, boolean uline, boolean sout) {

//		System.out.println("set font");
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_CREATEFONTINDIRECT);

			// 1. size:
			// Specifies the height, in logical units, of the font's character cell or character. 
			// The character height value (also known as the em height) is the character cell height 
			// value minus the internal-leading value.
			// 
			// The font mapper interprets the value specified here in the following manner.
			// Value 	Meaning
			//	> 0 	The font mapper transforms this value into device units and matches it against 
			//			the CELL HEIGHT of the available fonts.
			//	0 		The font mapper uses a default height value when it searches for a match.
			//	< 0 	The font mapper transforms this value into device units and matches its absolute value against 
			//			the CHARACTER HEIGHT of the available fonts.
			// For all height comparisons, the font mapper looks for the largest font that does not exceed 
			// the requested size. This mapping occurs when the font is used for the first time.
			//
			// Thus, correct the value again with screen resolution (96) against fix value for MM_TEXT mapping (72).
		tmpRecord.addParameter((-1)* ((correctValue(font.getSize())* scrRes)/ 72)); 
			
			// 2. width
			// Specifies the average width, in logical units, of characters in the font. 
			// If zero, the aspect ratio of the device is matched against the digitization aspect ratio of the 
			// available fonts to find the closest match, determined by the absolute value of the difference. 		tmpRecord.addParameter(0);
		tmpRecord.addParameter(0);

			// 3. escapement
			// Escapement: Specifies the angle, in tenths of degrees, between the escapement vector and the x-axis 
			// of the device. The escapement vector is parallel to the base line of a row of text.
		tmpRecord.addParameter(esc);

			// 4. orientation
			// Orientation: Specifies the angle, in tenths of degrees, between each character's base line and 
			// the x-axis of the device.
		tmpRecord.addParameter(ori);

			// 5. weight (bold/ plain)
			// Specifies the weight of the font in the range 0 through 1000. 
			// For example, 400 is normal and 700 is bold. If this value is zero, 
			// a default weight is used.		
		int weight= FW_DONTCARE;
		if (font.isBold())
			weight= FW_BOLD;
		if (font.isPlain())
			weight= FW_NORMAL;
		tmpRecord.addParameter(weight);
		
			// 6. italic + underline
		tmpRecord.addParameter(							
			((uline? 1: 0)<< 8)
			+ (font.isItalic()?1:0)
		);	

			// 7. strikeout + charset
		int charSet= 0;
		tmpRecord.addParameter(							
			((sout? 1: 0)<< 8)
			+ charSet
		);	
	
			// 8. outprecision + clipprecision
		tmpRecord.addParameter(0);	
		
			// 9. quality + pitch
		tmpRecord.addParameter(0);	
	
		
			// font name
		byte[] fName= font.getName().getBytes();
		for (int i= 0; (i+1)< fName.length; i+= 2)
			tmpRecord.addParameter((fName[i+1] << 8)+ fName[i]);
		if ((fName.length% 2)!= 0)
			tmpRecord.addParameter(fName[fName.length- 1]& 0x000000FF);
			
		wmfRecords.add(tmpRecord);
		
		selectObject(objCounter);
		if (objCounter> 0)
			deleteObject(objCounter- 1);
		setTextColor(Color.black);
		setBKMode(TRANSPARENT);
		++objCounter;
	}
	
	public void selectObject(int handle) {

//		System.out.println("sel obj "+handle);
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_SELECTOBJECT);
		tmpRecord.addParameter(handle);
		
		wmfRecords.add(tmpRecord);
	}
	
	public void deleteObject(int handle) {

//		System.out.println("del obj "+handle);
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_DELETEOBJECT);
		tmpRecord.addParameter(handle);
		
		wmfRecords.add(tmpRecord);
	}	
	
	public void setTerminator() {
		
//		System.out.println("set term");
		WmfRecord tmpRecord= new WmfRecord(0x0000);
		wmfRecords.add(tmpRecord);
	}
	

	/**
	 * The SetBkMode function sets the background mix mode of the specified device context. 
	 * The background mix mode is used with text, hatched brushes, and pen styles that are not solid lines. 
	 */
	public void setBKMode(int iBkMode) {
		
//		System.out.println("set bk");
		WmfRecord tmpRecord= new WmfRecord(META_SETBKMODE);
		tmpRecord.addParameter(iBkMode);
		tmpRecord.addParameter(iBkMode);
		
		wmfRecords.add(tmpRecord);
	}	
	
	public void setTextAlign() {
		
		WmfRecord tmpRecord= new WmfRecord(META_SETTEXTALIGN);
		tmpRecord.addParameter(0x0018);
		tmpRecord.addParameter(0);
		
		wmfRecords.add(tmpRecord);
	}		
	
	public void drawStringExt(String s, int x, int y, int option, Rectangle clipRect) {
		
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_EXTTEXTOUT);

		tmpRecord.addParameter(correctValue(y));
		tmpRecord.addParameter(correctValue(x));
		tmpRecord.addParameter(s.length());	// text length

		tmpRecord.addParameter(option);		// option: ETO_OPAQUE, ETO_GRAYED, ETO_CLIPPED

			// clipping
		if (((option& ETO_CLIPPED)!= 0)&& (clipRect!= null)) {
			tmpRecord.addParameter(correctValue(clipRect.y));
			tmpRecord.addParameter(correctValue(clipRect.x));
			tmpRecord.addParameter(correctValue(clipRect.height));
			tmpRecord.addParameter(correctValue(clipRect.width));
		}

			// write string
		byte[] sb= s.getBytes();
		for (int i= 0; (i+1)< sb.length; i+= 2)
			tmpRecord.addParameter((sb[i+1] << 8)+ sb[i]);
		if ((sb.length% 2)!= 0)
			tmpRecord.addParameter(sb[sb.length- 1]& 0x000000FF);

		wmfRecords.add(tmpRecord);
	}
	
	public void drawString(String s, int x, int y) {
		
//		System.out.println("drw str");
		WmfRecord tmpRecord;
		tmpRecord= new WmfRecord(META_TEXTOUT);

		tmpRecord.addParameter(s.length());	// text length
		byte[] sb= s.getBytes();
		for (int i= 0; (i+1)< sb.length; i+= 2)
			tmpRecord.addParameter((sb[i+1] << 8)+ sb[i]);
		if ((sb.length% 2)!= 0)
			tmpRecord.addParameter(sb[sb.length- 1]& 0x000000FF);

		tmpRecord.addParameter(correctValue(y));
		tmpRecord.addParameter(correctValue(x));

		wmfRecords.add(tmpRecord);
	}
	
	protected int correctValue(int orig) {
		
		return (int) (((double) orig / (double) scrRes) * (double) imgTPI);
	}
}
