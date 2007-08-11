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
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Vector;

import gphase.model.AbstractRegion;
import gphase.model.AbstractSite;
import gphase.model.DefaultRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Species;
import gphase.model.SpliceSite;
import gphase.model.TSSite;
import gphase.model.Transcript;

import javax.swing.JFrame;
import javax.swing.JPanel;

import sun.java2d.loops.DrawLine;

/**
 * 
 * 
 * @author msammeth
 */
public class SpliceOSigner_firstVersion extends JPanel {
	
	static final int SPLICE_ADVANCE= 20;	// x
	static final int TRANSCRIPT_ADVANCE= 50;	// y
	static final int LEFT_BORDER= 200;
	static final Color TSS_TSE_COLOR= Color.black;
	static final Color DONOR_COLOR= Color.blue;
	static final Color ACCEPTOR_COLOR= Color.red;
	static final Color EXON_COLOR= Color.green;
	AbstractRegion[] highlightRegions= null;
	
	Gene gene= null;

	static Gene getAGene() {
		Gene newGene= new Gene(new Species("homo_sapiens"), "ENSG00001000001");
		Transcript trans1= new Transcript(newGene, "ENST00000100001");
		newGene.addTranscript(trans1);
		Transcript trans2= new Transcript(newGene, "ENST00000100002");
		newGene.addTranscript(trans2);
		Exon e= new Exon(trans1, "ENSE00000100001", 50,100);
		trans1.addExon(e);
		trans2.addExon(e);
		e= new Exon(trans1, "ENSE00000100002", 150,200);
		trans1.addExon(e);
		e= new Exon(trans1, "ENSE00000100003", 250,350);
		trans1.addExon(e);
		e= new Exon(trans2, "ENSE00000100004", 300,350);
		trans2.addExon(e);
		
		return newGene;
	}
	
	
	public static void main(String[] args) {
		
		SpliceOSigner mySplicOSigner= new SpliceOSigner(getAGene());
		AbstractRegion[] regions= new AbstractRegion[] {
			new DefaultRegion(170,180)	
		};
		mySplicOSigner.setHighlightRegions(regions);
		JFrame myFrame= new JFrame();
		myFrame.getContentPane().add(mySplicOSigner);
		myFrame.pack();
		myFrame.setVisible(true);
	}
	
	public SpliceOSigner_firstVersion(Gene g) {
		gene= g;
	}
	
	protected void paintHighlights(Graphics g) {

		if (highlightRegions== null)
			return;
		
		AbstractSite[] sUniverse= gene.getSites();
		for (int i = 0; i < highlightRegions.length; i++) {
			
			int start= -1, end= -1;	// highlight region starts before the start-ths splice site 
									// (and ends before the end-ths splice site)
			
			for (int j = 0; j < sUniverse.length; j++) {
				if (start< 0&& highlightRegions[i].getStart()< sUniverse[j].getPos()) 
					start= j;
				if (highlightRegions[i].getEnd()< sUniverse[j].getPos()) { 
					end= j;
					break;
				}
			}
			
			if (start< 0)
				start= sUniverse.length;
			if (end< 0)
				end= sUniverse.length;
			
			g.setColor(new Color(255,255,255,200));
			Rectangle r= null;
			if (end== start)
				r= new Rectangle(
						LEFT_BORDER+ SPLICE_ADVANCE* start+ SPLICE_ADVANCE/ 4, 
						TRANSCRIPT_ADVANCE- TRANSCRIPT_ADVANCE/ 4, 
						SPLICE_ADVANCE/ 2, 
						TRANSCRIPT_ADVANCE* gene.getTranscripts().length+ TRANSCRIPT_ADVANCE/ 2
				);
			else
				r= new Rectangle(
						LEFT_BORDER+ SPLICE_ADVANCE* (start+1)- SPLICE_ADVANCE/ 2, 
						TRANSCRIPT_ADVANCE- TRANSCRIPT_ADVANCE/ 4, 
						SPLICE_ADVANCE* (end-start), 
						TRANSCRIPT_ADVANCE* gene.getTranscripts().length+ TRANSCRIPT_ADVANCE/ 2
				);
			
			g.fillRect(r.x, r.y, r.width, r.height);
		}
	}
	
	protected void paintComponent(Graphics g) {

		super.paintComponent(g);
		AbstractSite[] sUniverse= gene.getSites();
		
			// paint gene with splice sites
		g.setColor(Color.black);
		g.drawString(gene.getStableID(), 0, TRANSCRIPT_ADVANCE);
		g.drawLine(LEFT_BORDER, TRANSCRIPT_ADVANCE, LEFT_BORDER+ (sUniverse.length+ 1)* SPLICE_ADVANCE, TRANSCRIPT_ADVANCE);
		for (int i = 0; i < sUniverse.length; i++) {
			// dotted line
			g.setColor(Color.lightGray);
			g.drawLine(LEFT_BORDER+ SPLICE_ADVANCE* (i+1), TRANSCRIPT_ADVANCE, 
					LEFT_BORDER+ SPLICE_ADVANCE* (i+1), TRANSCRIPT_ADVANCE* (gene.getTranscripts().length+ 1));
			
			if (sUniverse[i] instanceof TSSite) {	// tss/tse
				g.setColor(TSS_TSE_COLOR);
				g.drawLine(LEFT_BORDER+ SPLICE_ADVANCE* (i+1), 
						TRANSCRIPT_ADVANCE- TRANSCRIPT_ADVANCE/ 10,
						LEFT_BORDER+ SPLICE_ADVANCE* (i+1),
						TRANSCRIPT_ADVANCE+ TRANSCRIPT_ADVANCE/ 10);
						
			} else if (sUniverse[i] instanceof SpliceSite&& ((SpliceSite) sUniverse[i]).isDonor()) {
				g.setColor(DONOR_COLOR);
				g.fillPolygon(
						new int[] {LEFT_BORDER+ SPLICE_ADVANCE* (i+1), 
								LEFT_BORDER+ SPLICE_ADVANCE* (i+1), 
								LEFT_BORDER+ SPLICE_ADVANCE* (i+1)+ SPLICE_ADVANCE/ 4},
						new int[] {TRANSCRIPT_ADVANCE- TRANSCRIPT_ADVANCE/ 10, 
								TRANSCRIPT_ADVANCE+ TRANSCRIPT_ADVANCE/ 10, 
								TRANSCRIPT_ADVANCE},
						3);

			} else  {		// acceptor
				g.setColor(ACCEPTOR_COLOR);
				g.fillPolygon(
						new int[] {LEFT_BORDER+ SPLICE_ADVANCE* (i+1), 
								LEFT_BORDER+ SPLICE_ADVANCE* (i+1), 
								LEFT_BORDER+ SPLICE_ADVANCE* (i+1)- SPLICE_ADVANCE/ 4},
						new int[] {TRANSCRIPT_ADVANCE- TRANSCRIPT_ADVANCE/ 10, 
								TRANSCRIPT_ADVANCE+ TRANSCRIPT_ADVANCE/ 10, 
								TRANSCRIPT_ADVANCE},
						3);
			}
		}
		
			// transcripts and exons
		for (int i = 0; i < gene.getTranscripts().length; i++) {
			g.setColor(Color.black);
			g.drawString(gene.getTranscripts()[i].getStableID(), 0, TRANSCRIPT_ADVANCE* (i+2));
			g.drawLine(LEFT_BORDER, TRANSCRIPT_ADVANCE* (i+2), LEFT_BORDER+ (sUniverse.length+ 1)* SPLICE_ADVANCE, TRANSCRIPT_ADVANCE* (i+2));
			Exon[] exons= gene.getTranscripts()[i].getExons();
			Comparator compi= new SpliceSite.PositionComparator();
			for (int j = 0; j < exons.length; j++) {
				int dPos= Arrays.binarySearch(sUniverse, exons[j].getStartSite(), compi);
				int aPos= Arrays.binarySearch(sUniverse, exons[j].getEndSite(), compi);
				g.setColor(EXON_COLOR);
				Rectangle r= new Rectangle(LEFT_BORDER+ (dPos+ 1)* SPLICE_ADVANCE,
						TRANSCRIPT_ADVANCE* (i+ 2)- TRANSCRIPT_ADVANCE/ 8,
						(aPos- dPos)* SPLICE_ADVANCE,
						TRANSCRIPT_ADVANCE/ 4);
				g.fillRect(r.x, r.y, r.width, r.height);
				g.setColor(Color.black);
				g.drawRect(r.x, r.y, r.width, r.height);
			}
		}
		paintHighlights(g);
		
			// draw footer
		g.setColor(Color.lightGray);
		String sign= "SpliceOSigner 1.0";
		int signWidth= g.getFontMetrics().stringWidth(sign);
		g.drawString(sign,
				getPreferredSize().width- signWidth, getPreferredSize().height- g.getFontMetrics().getHeight());
	}
	
	/* (non-Javadoc)
	 * @see javax.swing.JComponent#getPreferredSize()
	 */
	public Dimension getPreferredSize() {
		return new Dimension(LEFT_BORDER+ SPLICE_ADVANCE* (gene.getSites().length+ 2), TRANSCRIPT_ADVANCE* (gene.getTranscripts().length+ 2));
	}
	public AbstractRegion[] getHighlightRegions() {
		return highlightRegions;
	}
	public void setHighlightRegions(AbstractRegion[] highlightRegions) {
		this.highlightRegions = highlightRegions;
	}
}
