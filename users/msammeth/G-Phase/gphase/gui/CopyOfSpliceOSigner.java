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
import java.util.Vector;

import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.AbstractSite;
import gphase.model.DefaultRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.TSSite;
import gphase.model.Transcript;
import gphase.model.ASVariation.SpliceChainComparator;

import javax.swing.JFrame;
import javax.swing.JPanel;

import sun.java2d.loops.DrawLine;

/**
 * 
 * 
 * @author msammeth
 */
public class CopyOfSpliceOSigner extends JPanel {
	
	static final int BORDER= 5;
	static final int VBORDER= 2;
	static final int SPLICE_HEIGHT= 15;
	static final Dimension INTRON_DIMENSION= new Dimension(15, 3);
	static final Dimension EXON_DIMENSION= new Dimension(30, 10);
	static final Color EXON_COLOR= Color.green.darker().darker();
	AbstractRegion[] highlightRegions= null;
	
	ASVariation variation;

	public static void main(String[] args) {
		
		CopyOfSpliceOSigner mySplicOSigner= new CopyOfSpliceOSigner(getAGene());
		AbstractRegion[] regions= new AbstractRegion[] {
			new DefaultRegion(170,180)	
		};
		mySplicOSigner.setHighlightRegions(regions);
		JFrame myFrame= new JFrame();
		myFrame.getContentPane().add(mySplicOSigner);
		myFrame.pack();
		myFrame.setVisible(true);
	}
	
	public CopyOfSpliceOSigner(ASVariation var) {
		variation= var;
		setOpaque(false);
	}
	
	protected void paintComponent(Graphics g) {

		Transcript trpt1= null;
		if (variation.getSpliceChain2().length== 0|| 
				(variation.getSpliceChain1().length> 0&& variation.getSpliceChain1()[0].getPos()< variation.getSpliceChain2()[0].getPos()))
			trpt1= variation.getTranscript1();
		else
			trpt1= variation.getTranscript2();
		
		//super.paintComponent(g);
		SpliceSite[] su= variation.getSpliceUniverse();
		
			// left border
		int height= 2* SPLICE_HEIGHT+ EXON_DIMENSION.height;
		if (su[0].isAcceptor()) {	// -> conserved site is donor and paint exon left of it
			g.setColor(EXON_COLOR);
			g.fillRect(
					0, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
			g.setColor(Color.black);
			g.drawRect(
					0, height/ 2- EXON_DIMENSION.height/ 2,
					BORDER, EXON_DIMENSION.height
			);
		} else {
			g.setColor(Color.black);
			g.fillRect(
					0, height/ 2- INTRON_DIMENSION.height/ 2,
					BORDER, INTRON_DIMENSION.height
			);
		}

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
			int delta;
			if (exonic1|| exonic2) {
				delta= EXON_DIMENSION.width;
				g.setColor(EXON_COLOR);
				g.fillRect(
						xPos, height/ 2- EXON_DIMENSION.height/ 2,
						delta, EXON_DIMENSION.height
				);
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
					delta, EXON_DIMENSION.height
			);
			g.setColor(Color.black);
			g.drawRect(
					xPos, height/ 2- EXON_DIMENSION.height/ 2,
					delta, EXON_DIMENSION.height
			);
			
			xPos+= delta;
			
			g.setColor(Color.black);
			g.fillRect(
					xPos, height/ 2- INTRON_DIMENSION.height/ 2,
					BORDER, INTRON_DIMENSION.height
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
			int delta;
			if (exonic1|| exonic2) 
				w+= EXON_DIMENSION.width;
			else
				w+= INTRON_DIMENSION.width;
			
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
			w+= EXON_DIMENSION.width;
		else
			w+= INTRON_DIMENSION.width;
		
		return new Dimension(
				2* BORDER+ w,
				2* SPLICE_HEIGHT+ EXON_DIMENSION.height
		);
	}
}
