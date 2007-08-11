/*
 * Created on Mar 3, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Vector;

import gphase.io.TabDelimitedFormatWrapper;
import gphase.model.ASVariation;
import gphase.model.ASVariationDomain;
import gphase.model.AbstractRegion;
import gphase.model.AbstractSite;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.TSSite;
import gphase.model.Transcript;
import gphase.model.ASVariation.SpliceChainComparator;
import gphase.tools.File;


import javax.swing.JFrame;
import javax.swing.JPanel;

import sun.java2d.loops.DrawLine;

/**
 * 
 * 
 * @author msammeth
 */
public class DomainOSigner extends JPanel {
	static final String DEFAULT_COLOR_MAP= "color.map";
	static HashMap<String,Color> colorMap= new HashMap<String,Color>();
	static {
		readInDomainColorMap(new File(DEFAULT_COLOR_MAP));
	}
	static int delme= 0;
	static final int BORDER= 5;
	static final int VBORDER= 2;
	static final int SPLICE_HEIGHT= 15;
	static final Dimension INTRON_DIMENSION= new Dimension(15, 3);
	static final Dimension EXON_DIMENSION= new Dimension(30, 10);
	static final Color EXON_COLOR= Color.green.darker().darker();
	static final Color[] DOMAIN_COLORS= new Color[] {
		Color.red, Color.cyan, Color.pink, Color.yellow, Color.orange, Color.blue, Color.magenta,
		Color.red.darker(), Color.cyan.darker(), Color.pink.darker(), Color.yellow.darker(), Color.orange.darker(), Color.blue.darker(), Color.magenta.darker(),
		Color.red.brighter(), Color.cyan.brighter(), Color.pink.brighter(), Color.yellow.brighter(), Color.orange.brighter(), Color.blue.brighter(), Color.magenta.brighter()
	};
	AbstractRegion[] highlightRegions= null;
	
	ASVariationDomain variation;

	public static HashMap<String,Color> readInDomainColorMap(File inFile) {
		if (!inFile.exists()) {
			System.out.println("No color configuration found for domains.");
			return colorMap;
		}
		try {
			TabDelimitedFormatWrapper reader= new TabDelimitedFormatWrapper(inFile.getAbsolutePath());
			reader.read();
			String[][] tab= reader.getTable();
			for (int i = 0; tab!= null&& i < tab.length; i++) {
				if (tab[i].length!= 4)
					continue;
				Color c= new Color(Integer.parseInt(tab[i][1]), 
						Integer.parseInt(tab[i][2]), Integer.parseInt(tab[i][3]));
				colorMap.put(tab[i][0], c);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return colorMap;
	}
	
	public static void writeOutDomainColorMap() {
		writeOutDomainColorMap(colorMap, new File(DEFAULT_COLOR_MAP));
	}
	public static void writeOutDomainColorMap(HashMap<String,Color> map, File outFile) {
		try {
			String[][] tab= new String[map.size()][];
			Object[] keys= map.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				Color c= map.get(keys[i]);
				tab[i]= new String[4];
				tab[i][0]= (String) keys[i];
				tab[i][1]= Integer.toString(c.getRed());
				tab[i][2]= Integer.toString(c.getGreen());
				tab[i][3]= Integer.toString(c.getBlue());
			}
			TabDelimitedFormatWrapper writer= new TabDelimitedFormatWrapper(outFile.getAbsolutePath());
			writer.setTable(tab);
			writer.write();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		
	}
	
	public static Color getNewColor(Color[] cols) {
		// http://www.compuphase.com/cmetric.htm
		//System.out.println("cols used: "+(delme+1));
		return DOMAIN_COLORS[delme++%DOMAIN_COLORS.length];
	}
	
	public DomainOSigner(ASVariationDomain var) {
		variation= var;
		setOpaque(false);
	}
	
	private Color getDomainColor(int counter) {
		return DOMAIN_COLORS[counter% DOMAIN_COLORS.length];
	}
	
	public static Color getDomainColor(String domID) {
		String id= ASVariationDomain.stripOrderSuffix(domID);
		Color col= colorMap.get(id);
		if (col== null) {
			col= getNewColor((Color[]) gphase.tools.Arrays.toField(colorMap.values()));
			colorMap.put(id, col);
		}
		return col;
	}
	
	public static Color[] getDomainColors(String[] domIDs) {
		Color[] cols= new Color[domIDs.length];
		for (int i = 0; i < domIDs.length; i++) 
			cols[i]= getDomainColor(domIDs[i]);
		return cols;
	}
	
	private int paintDomain(Graphics g, int xPos, int delta, int height, boolean up, int i, int domCtr1, int[] dom1Starts, int[] dom1Ends, DirectedRegion[] dom1) {
		
		for (int j = 0; j < dom1.length; j++) {
			for (int k = (j+1); k < dom1.length; k++) {
				if (dom1[j].overlaps(dom1[k]))
					System.currentTimeMillis();
			}
		}
		
			// countSegments
		int segCnt= 0, domCnt= 0;
		if (dom1Starts[domCtr1]== i)	// starts with spacer
			++segCnt;
		for (int j = domCtr1; j < dom1.length; j++) {
			if (dom1Starts[j]> i)
				break;
			if (j> 0&& dom1Ends[j-1]== i&& dom1Starts[j]== i) // adjacent in exon
				++segCnt;
			if (dom1Starts[j]== i&& dom1Ends[j]== i)	// starts and ends in exon
				segCnt+= 2;
			else //if (dom1Starts[j]== i^ dom1Ends[j]== i)	// starts or ends
				segCnt+= 1;
			// starts and ends outside of exon
			if (dom1Ends[j]== i)	// close domain
				++domCnt;
		}
		if (domCnt> 0&& (dom1Ends[domCtr1+ domCnt- 1]== i&& 
				(domCtr1+ domCnt== dom1.length|| dom1Starts[domCtr1+ domCnt]> i)))	// last domain needs spacer
			++segCnt;
		int xDom= xPos, deltaSeg= delta/ segCnt;
		boolean rest= (delta% segCnt)!= 0;
		if (dom1Starts[domCtr1]== i) {	// starts with spacer
			xDom+= deltaSeg;
			if (rest)
				++xDom;
		}
		
		int yPos= height/ 2;
		if (up)
			yPos-= EXON_DIMENSION.height/ 2;
		for (int j = domCtr1; j < dom1.length; j++) {
			if (dom1Starts[j]> i)
				break;
			int deltaDom= deltaSeg;
			if (dom1Starts[j]== i&& dom1Ends[j]== i)	// starts and ends in exon
				deltaDom+= deltaSeg;
			g.setColor(getDomainColor(dom1[j].getID()));
			g.fillRect(
					xDom, yPos,
					deltaDom, EXON_DIMENSION.height/ 2
			);
			if (dom1Ends[j]> i) {
				xDom+= deltaDom;
				break;
			}
			xDom+= deltaDom+ deltaSeg;
			if (rest)
				++xDom;
		}
		if (xDom< xPos+ delta&& domCtr1+ domCnt< dom1.length&& 
				dom1Starts[domCtr1+ domCnt]<= i&&
				dom1Ends[domCtr1+ domCnt]> i) {
			g.setColor(getDomainColor(dom1[domCtr1+ domCnt].getID()));
			g.fillRect(
					xDom, yPos,
					xPos+ delta- xDom, EXON_DIMENSION.height/ 2
			);
		}
			
		domCtr1+= domCnt;
		return domCtr1;
	}

	private int paintDomain_save(Graphics g, int xPos, int delta, int height, boolean up, int i, int domCtr1, int[] dom1Starts, int[] dom1Ends, DirectedRegion[] dom1) {
			// countSegments
		int segCnt= 0, domCnt= 0;
		if (dom1Starts[domCtr1]== i)	// starts with spacer
			++segCnt;
		for (int j = domCtr1; j < dom1.length; j++) {
			if (dom1Starts[j]> i)
				break;
			if (j> 0&& dom1Ends[j-1]== i&& dom1Starts[j]== i) // adjacent in exon
				++segCnt;
			if (dom1Starts[j]== i&& dom1Ends[j]== i)	// starts and ends in exon
				segCnt+= 2;
			else //if (dom1Starts[j]== i^ dom1Ends[j]== i)	// starts or ends
				segCnt+= 1;
			// starts and ends outside of exon
			if (dom1Ends[j]== i)	// close domain
				++domCnt;
		}
		if (domCnt> 0&& (dom1Ends[domCtr1+ domCnt- 1]== i&& 
				(domCtr1+ domCnt== dom1.length|| dom1Starts[domCtr1+ domCnt]> i)))	// last domain needs spacer
			++segCnt;
		int xDom= xPos, deltaSeg= delta/ segCnt;
		boolean rest= (delta% segCnt)!= 0;
		if (dom1Starts[domCtr1]== i) {	// starts with spacer
			xDom+= deltaSeg;
			if (rest)
				++xDom;
		}
		
		int yPos= height/ 2;
		if (up)
			yPos-= EXON_DIMENSION.height/ 2;
		for (int j = domCtr1; j < dom1.length; j++) {
			if (dom1Starts[j]> i)
				break;
			int deltaDom= deltaSeg;
			if (dom1Starts[j]== i&& dom1Ends[j]== i)	// starts and ends in exon
				deltaDom+= deltaSeg;
			g.setColor(getDomainColor(dom1[j].getID()));
			g.fillRect(
					xDom, yPos,
					deltaDom, EXON_DIMENSION.height/ 2
			);
			if (dom1Ends[j]> i) {
				xDom+= deltaDom;
				break;
			}
			xDom+= deltaDom+ deltaSeg;
			if (rest)
				++xDom;
		}
		if (xDom< xPos+ delta&& domCtr1+ domCnt< dom1.length&& 
				dom1Starts[domCtr1+ domCnt]<= i&&
				dom1Ends[domCtr1+ domCnt]> i) {
			g.setColor(getDomainColor(dom1[domCtr1+ domCnt].getID()));
			g.fillRect(
					xDom, yPos,
					xPos+ delta- xDom, EXON_DIMENSION.height/ 2
			);
		}
			
		domCtr1+= domCnt;
		return domCtr1;
	}
	
	protected void paintComponent(Graphics g) {
	
			setBackground(Color.white);
			super.paintComponent(g);
			
			Transcript trpt1= variation.getTranscript1();
//			if (variation.getSpliceChain2().length== 0|| 
//					(variation.getSpliceChain1().length> 0&& variation.getSpliceChain1()[0].getPos()< variation.getSpliceChain2()[0].getPos()))
//				trpt1= variation.getTranscript1();
//			else
//				trpt1= variation.getTranscript2();
			
			//super.paintComponent(g);
			SpliceSite[] su= variation.getSpliceUniverse();
			DirectedRegion[] dom1= variation.getDom1(), dom2= variation.getDom2();
			int[] dom1Starts= variation.getDom1Starts();
			int[] dom1Ends= variation.getDom1Ends();
			int[] dom2Starts= variation.getDom2Starts();
			int[] dom2Ends= variation.getDom2Ends();
			int domCtr1= 0, domCtr2= 0;
	
			
				// left border
			int height= 2* SPLICE_HEIGHT+ EXON_DIMENSION.height;
			g.setColor(EXON_COLOR);
			g.fillRect(
					0, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
			if (domCtr1< dom1Starts.length&& dom1Starts[domCtr1]== 0) {	// domain starts before variable area
				g.setColor(getDomainColor(dom1[domCtr1].getID()));
				g.fillRect(
						0, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}
			if (domCtr2< dom2Starts.length&& dom2Starts[domCtr2]== 0) {
				g.setColor(getDomainColor(dom2[domCtr2].getID()));
				g.fillRect(
						0, height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}
			
			g.setColor(Color.black);
	//		g.drawRect(
	//				0, height/ 2- EXON_DIMENSION.height/ 2,
	//				BORDER, EXON_DIMENSION.height
	//		);
			g.drawLine(0, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, height/ 2- EXON_DIMENSION.height/ 2);
			g.drawLine(BORDER, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, height/ 2+ EXON_DIMENSION.height/ 2);
			g.drawLine(0, height/ 2+ EXON_DIMENSION.height/ 2,
					BORDER, height/ 2+ EXON_DIMENSION.height/ 2);
			
	
				// paint all
			int lastX1= -1;
			int lastX2= -1;
			boolean exonic1, exonic2;
			if (su[0].isDonor()) {
				exonic1= true; exonic2= true;
			} else {
				exonic1= false; exonic2= false;
				lastX1= BORDER; lastX2= BORDER;
			}
			int xPos= BORDER;
			for (int i = 0; i < su.length; i++) {
				if (i== 0&& su[i].isDonor()) {
					if (trpt1.containsSS(su[i])) {
						exonic1= false;
						lastX1= xPos;
					} else {
						exonic2= false;
						lastX2= xPos;
					}
					continue;
				}
				int delta;
				if (exonic1|| exonic2) {	// exonic: i is the end of the end of the exonic range
					delta= EXON_DIMENSION.width;
					g.setColor(EXON_COLOR);
					g.fillRect(
							xPos, height/ 2- EXON_DIMENSION.height/ 2,
							delta, EXON_DIMENSION.height
					);
					int i1= variation.getSplicePos1(variation.getSpliceUniverse()[i]);
					if (exonic1&& domCtr1< dom1Starts.length&& dom1Starts[domCtr1]<= i) {	// domain starts before end of exon

						domCtr1= paintDomain(g, xPos, delta, height, true, i, domCtr1, dom1Starts, dom1Ends, dom1);
					}
					int i2= variation.getSplicePos2(variation.getSpliceUniverse()[i]);
					if (exonic2&& domCtr2< dom2Starts.length&& dom2Starts[domCtr2]<= i) {	// domain starts before end of exonic
						domCtr2= paintDomain(g, xPos, delta, height, false, i, domCtr2, dom2Starts, dom2Ends, dom2);
					}
					
					g.setColor(Color.black);
					g.drawRect(
							xPos, height/ 2- EXON_DIMENSION.height/ 2,
							delta, EXON_DIMENSION.height
					);
				} else {
					delta= INTRON_DIMENSION.width;
					g.setColor(Color.black);
					g.fillRect(
							xPos, height/ 2- INTRON_DIMENSION.height/ 2,
							delta, INTRON_DIMENSION.height
					);
				}
				
				xPos+= delta;
				if (trpt1.containsSS(su[i])) {
					if (su[i].isDonor()) {
						exonic1= false;
						lastX1= xPos;
					} else {
						exonic1= true;
						int diff= xPos- lastX1;
						g.setColor(Color.black);
						g.drawLine(lastX1, height/2- EXON_DIMENSION.height/ 2,
								lastX1+ diff/2, VBORDER);
						g.drawLine(lastX1+ diff/2, VBORDER, 
								xPos, height/2- EXON_DIMENSION.height/ 2);
					}
				} else {
					if (su[i].isDonor()) {
						exonic2= false;
						lastX2= xPos;
					} else {
						exonic2= true;
						int diff= xPos- lastX2;
						g.setColor(Color.black);
						g.drawLine(lastX2, height/2+ EXON_DIMENSION.height/ 2,
								lastX2+ diff/2, height- VBORDER);
						g.drawLine(lastX2+ diff/2, height- VBORDER, 
								xPos, height/2+ EXON_DIMENSION.height/ 2);
					}
				}
	
			}
					
				// paint last one
			if (su[su.length- 1].isAcceptor()) {	// following is a donor
				int delta= EXON_DIMENSION.width;
				g.setColor(EXON_COLOR);
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height
				);
				if (domCtr1< dom1Starts.length) {	// domain left
					g.setColor(getDomainColor(dom1[domCtr1].getID()));
					g.fillRect(
							xPos, height/ 2- EXON_DIMENSION.height/ 2,
							BORDER, EXON_DIMENSION.height/ 2
					);
				}			
				if (domCtr2< dom2Starts.length) {	// domain left
					g.setColor(getDomainColor(dom2[domCtr2].getID()));
					g.fillRect(
							xPos, height/ 2,
							BORDER, EXON_DIMENSION.height/ 2
					);
				}
				
				g.setColor(Color.black);
				g.drawRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height
				);
				
			} else {
				int delta= INTRON_DIMENSION.width;
				g.setColor(Color.black);
				g.fillRect(
						xPos, height/ 2- INTRON_DIMENSION.height/ 2,
						delta, INTRON_DIMENSION.height
				);
	
				xPos+= delta;
				if (!exonic1) {
					exonic1= true;
					int diff= xPos- lastX1;
					g.setColor(Color.black);
					g.drawLine(lastX1, height/2- EXON_DIMENSION.height/ 2,
							lastX1+ diff/2, VBORDER);
					g.drawLine(lastX1+ diff/2, VBORDER, 
							xPos, height/2- EXON_DIMENSION.height/ 2);
				} 
				if (!exonic2) {
					exonic2= true;
					int diff= xPos- lastX2;
					g.setColor(Color.black);
					g.drawLine(lastX2, height/2+ EXON_DIMENSION.height/ 2,
							lastX2+ diff/2, height- VBORDER);
					g.drawLine(lastX2+ diff/2, height- VBORDER, 
							xPos, height/2+ EXON_DIMENSION.height/ 2);
				}
	
				g.setColor(EXON_COLOR);
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height
				);
				if (domCtr1< dom1Starts.length) {	// domain left
					g.setColor(getDomainColor(dom1[domCtr1].getID()));
					g.fillRect(
							xPos, height/ 2- EXON_DIMENSION.height/ 2,
							BORDER, EXON_DIMENSION.height/ 2
					);
				}			
				if (domCtr2< dom2Starts.length) {	// domain left
					g.setColor(getDomainColor(dom2[domCtr2].getID()));
					g.fillRect(
							xPos, height/ 2,
							BORDER, EXON_DIMENSION.height/ 2
					);
				}			
	
				
				g.setColor(Color.black);
				g.drawRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height
				);
			}
		}

	protected void paintComponent_save(Graphics g) {

		setBackground(Color.white);
		super.paintComponent(g);
		
		Transcript trpt1= null;
		if (variation.getSpliceChain2().length== 0|| 
				(variation.getSpliceChain1().length> 0&& variation.getSpliceChain1()[0].getPos()< variation.getSpliceChain2()[0].getPos()))
			trpt1= variation.getTranscript1();
		else
			trpt1= variation.getTranscript2();
		
		//super.paintComponent(g);
		SpliceSite[] su= variation.getSpliceUniverse();
		DirectedRegion[] dom1= variation.getDom1(), dom2= variation.getDom2();
		int[] dom1Starts= variation.getDom1Starts();
		int[] dom1Ends= variation.getDom1Ends();
		int[] dom2Starts= variation.getDom2Starts();
		int[] dom2Ends= variation.getDom2Ends();
		int domCtr1= 0, domCtr2= 0, domColor1= -1, domColor2= -1;
		if (dom1Starts.length> 0)
			domColor1= Math.max(domColor1,domColor2)+ 1;
		if (dom2Starts.length> 0)
			domColor2= Math.max(domColor1,domColor2)+ 1;

		
			// left border
		int height= 2* SPLICE_HEIGHT+ EXON_DIMENSION.height;
		g.setColor(EXON_COLOR);
		g.fillRect(
				0, height/ 2- EXON_DIMENSION.height/ 2,
				BORDER, EXON_DIMENSION.height
		);
		if (domCtr1< dom1Starts.length&& dom1Starts[domCtr1]== 0) {	// domain starts before variable area
			g.setColor(getDomainColor(domColor1));
			g.fillRect(
					0, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height/ 2
			);
		}
		if (domCtr2< dom2Starts.length&& dom2Starts[domCtr2]== 0) {
			g.setColor(getDomainColor(domColor2));
			g.fillRect(
					0, height/ 2,
					BORDER, EXON_DIMENSION.height/ 2
			);
		}
		
		g.setColor(Color.black);
//		g.drawRect(
//				0, height/ 2- EXON_DIMENSION.height/ 2,
//				BORDER, EXON_DIMENSION.height
//		);
		g.drawLine(0, height/ 2- EXON_DIMENSION.height/ 2,
				BORDER, height/ 2- EXON_DIMENSION.height/ 2);
		g.drawLine(BORDER, height/ 2- EXON_DIMENSION.height/ 2,
				BORDER, height/ 2+ EXON_DIMENSION.height/ 2);
		g.drawLine(0, height/ 2+ EXON_DIMENSION.height/ 2,
				BORDER, height/ 2+ EXON_DIMENSION.height/ 2);
		

			// paint all
		int lastX1= -1;
		int lastX2= -1;
		boolean exonic1, exonic2;
		if (su[0].isDonor()) {
			exonic1= true; exonic2= true;
		} else {
			exonic1= false; exonic2= false;
			lastX1= BORDER; lastX2= BORDER;
		}
		int xPos= BORDER;
		for (int i = 0; i < su.length; i++) {
			if (i== 0&& su[i].isDonor()) {
				if (trpt1.containsSS(su[i])) {
					exonic1= false;
					lastX1= xPos;
				} else {
					exonic2= false;
					lastX2= xPos;
				}
				continue;
			}
			int delta;
			if (exonic1|| exonic2) {	// exonic: i is the end of the end of the exonic range
				delta= EXON_DIMENSION.width;
				g.setColor(EXON_COLOR);
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						delta, EXON_DIMENSION.height
				);
				if (exonic1&& domCtr1< dom1Starts.length&& dom1Starts[domCtr1]<= i) {	// domain starts before end of exonic
					g.setColor(getDomainColor(domColor1));
					int xDom= xPos, deltaDom= delta;
					if (dom1Starts[domCtr1]<= i-1) {	// dom starts before start of exonic region
						xDom= xPos;
						if (dom1Ends[domCtr1]<= i) { 	// dom ends before end of exon
							deltaDom= delta/2;
							++domCtr1;	// next domain
							domColor1= Math.max(domColor1,domColor2)+ 1;
						} else 	// dom starts before exon and ends after exon
							deltaDom= delta;
					} else {	// dom starts after start of exon
						deltaDom= delta/ 2;
						if (dom1Ends[domCtr1]<= i) { 	// dom ends before end of exon
							xDom= xPos+ (delta/4);
							++domCtr1;
							domColor1= Math.max(domColor1,domColor2)+ 1;
						} else 	// dom starts in exon and ends after exon
							xDom= xPos+ (delta/2);
					}
					g.fillRect(
							xDom, height/ 2- EXON_DIMENSION.height/ 2,
							deltaDom, EXON_DIMENSION.height/ 2
					);
				}
				if (exonic2&& domCtr2< dom2Starts.length&& dom2Starts[domCtr2]<= i) {	// domain starts before end of exonic
					g.setColor(getDomainColor(domColor2));
					int xDom= xPos, deltaDom= delta;
					if (dom2Starts[domCtr2]<= i-1) {	// dom starts before start of exonic region
						xDom= xPos;
						if (dom2Ends[domCtr2]<= i) { 	// dom ends before end of exon
							deltaDom= delta/2;
							++domCtr2;	// next domain
							domColor2= Math.max(domColor1,domColor2)+ 1;
						} else 	// dom starts before exon and ends after exon
							deltaDom= delta;
					} else {	// dom starts after start of exon
						deltaDom= delta/ 2;
						if (dom2Ends[domCtr2]<= i) { 	// dom ends before end of exon
							xDom= xPos+ (delta/4);
							++domCtr2;
							domColor2= Math.max(domColor1,domColor2)+ 1;
						} else 	// dom starts in exon and ends after exon
							xDom= xPos+ (delta/2);
					}
					g.fillRect(
							xDom, height/ 2,
							deltaDom, EXON_DIMENSION.height/ 2
					);
				}
				
				g.setColor(Color.black);
				g.drawRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						delta, EXON_DIMENSION.height
				);
			} else {
				delta= INTRON_DIMENSION.width;
				g.setColor(Color.black);
				g.fillRect(
						xPos, height/ 2- INTRON_DIMENSION.height/ 2,
						delta, INTRON_DIMENSION.height
				);
			}
			
			xPos+= delta;
			if (trpt1.containsSS(su[i])) {
				if (su[i].isDonor()) {
					exonic1= false;
					lastX1= xPos;
				} else {
					exonic1= true;
					int diff= xPos- lastX1;
					g.setColor(Color.black);
					g.drawLine(lastX1, height/2- EXON_DIMENSION.height/ 2,
							lastX1+ diff/2, VBORDER);
					g.drawLine(lastX1+ diff/2, VBORDER, 
							xPos, height/2- EXON_DIMENSION.height/ 2);
				}
			} else {
				if (su[i].isDonor()) {
					exonic2= false;
					lastX2= xPos;
				} else {
					exonic2= true;
					int diff= xPos- lastX2;
					g.setColor(Color.black);
					g.drawLine(lastX2, height/2+ EXON_DIMENSION.height/ 2,
							lastX2+ diff/2, height- VBORDER);
					g.drawLine(lastX2+ diff/2, height- VBORDER, 
							xPos, height/2+ EXON_DIMENSION.height/ 2);
				}
			}

		}
				
			// paint last one
		if (su[su.length- 1].isAcceptor()) {	// following is a donor
			int delta= EXON_DIMENSION.width;
			g.setColor(EXON_COLOR);
			g.fillRect(
					xPos, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
			if (domCtr1< dom1Starts.length) {	// domain left
				g.setColor(getDomainColor(domColor1));
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}			
			if (domCtr2< dom2Starts.length) {	// domain left
				g.setColor(getDomainColor(domColor2));
				g.fillRect(
						xPos, height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}
			
			g.setColor(Color.black);
			g.drawRect(
					xPos, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
			
		} else {
			int delta= INTRON_DIMENSION.width;
			g.setColor(Color.black);
			g.fillRect(
					xPos, height/ 2- INTRON_DIMENSION.height/ 2,
					delta, INTRON_DIMENSION.height
			);

			xPos+= delta;
			if (!exonic1) {
				exonic1= true;
				int diff= xPos- lastX1;
				g.setColor(Color.black);
				g.drawLine(lastX1, height/2- EXON_DIMENSION.height/ 2,
						lastX1+ diff/2, VBORDER);
				g.drawLine(lastX1+ diff/2, VBORDER, 
						xPos, height/2- EXON_DIMENSION.height/ 2);
			} 
			if (!exonic2) {
				exonic2= true;
				int diff= xPos- lastX2;
				g.setColor(Color.black);
				g.drawLine(lastX2, height/2+ EXON_DIMENSION.height/ 2,
						lastX2+ diff/2, height- VBORDER);
				g.drawLine(lastX2+ diff/2, height- VBORDER, 
						xPos, height/2+ EXON_DIMENSION.height/ 2);
			}

			g.setColor(EXON_COLOR);
			g.fillRect(
					xPos, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
			if (domCtr1< dom1Starts.length) {	// domain left
				g.setColor(getDomainColor(domColor1));
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}			
			if (domCtr2< dom2Starts.length) {	// domain left
				g.setColor(getDomainColor(domColor2));
				g.fillRect(
						xPos, height/ 2,
						BORDER, EXON_DIMENSION.height/ 2
				);
			}			

			
			g.setColor(Color.black);
			g.drawRect(
					xPos, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
		}
	}
	
	/* (non-Javadoc)
	 * @see javax.swing.JComponent#getPreferredSize()
	 */
	public Dimension getPreferredSize() {
		
		SpliceSite[] su= variation.getSpliceUniverse();
		int w= 0;
		boolean exonic1, exonic2;
		if (su[0].isDonor()) {
			exonic1= true; exonic2= true;
		} else {
			exonic1= false; exonic2= false;
		}
		for (int i = 0; i < su.length; i++) {
			if (i== 0&& su[i].isDonor()) 
				w+= 0;
			else {
				if (exonic1|| exonic2) 
					w+= EXON_DIMENSION.width;
				else
					w+= INTRON_DIMENSION.width;
			}
			
			if (variation.getTranscript1().containsSS(su[i])) {
				if (su[i].isDonor()) 
					exonic1= false;
				else 
					exonic1= true;
			} else {
				if (su[i].isDonor()) 
					exonic2= false;
				else
					exonic2= true;
			}

		}
				
			// last one
		if (su[su.length- 1].isAcceptor()) 	// following is a donor
			w+= 0;
		else
			w+= INTRON_DIMENSION.width;
		
		return new Dimension(
				2* BORDER+ w,
				2* SPLICE_HEIGHT+ EXON_DIMENSION.height
		);
	}

	public static HashMap<String, Color> getColorMap() {
		return colorMap;
	}
}
