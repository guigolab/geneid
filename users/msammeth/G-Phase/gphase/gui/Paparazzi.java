/*
 * Created on Mar 8, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.gui;

import gphase.db.EnsemblDBAdaptor;
import gphase.model.AbstractRegion;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;

import java.awt.Component;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Properties;

import javax.swing.JLabel;
import javax.swing.JOptionPane;

import org.ensembl.datamodel.SnapShot;
import org.freehep.util.export.ExportFileType;


/**
 * 
 * 
 * @author msammeth
 */
public class Paparazzi {
	
		public static void main(String[] args) {
			Paparazzi papa= new Paparazzi();
			File f= new File("test.gif");
			
			EnsemblDBAdaptor adaptor= new EnsemblDBAdaptor();
			Graph g= adaptor.getGraphAllGenes(new Species(EnsemblDBAdaptor.SPECIES_SMALL[0]));
			Gene[] tst= new Gene[100];
			Gene[] ge= g.getGenes();
			for (int i = 0; i < tst.length; i++) 
				tst[i]= ge[i];
			
			burst(tst, new File("tyler"));
		}
		
	   public static void burst(Gene[] genes, File dir) {
	   
	   		for (int i = 0; i < genes.length; i++) {
				try {
			   		if (genes[i]!= null)
			   			writeGIF(snapShot(genes[i]), new File(dir.getPath()+ File.separator+ genes[i].getStableID()+".gif"));
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
	   }
	   
	   public static void burst(Gene[] genes, AbstractRegion[] highlightAreas, File dir) {
		   
		   		for (int i = 0; i < genes.length; i++) {
					try {
				   		if (genes[i]!= null)
				   			writeGIF(snapShot(genes[i], highlightAreas[i]), 
				   					new File(dir.getPath()+ File.separator+ genes[i].getStableID()+".gif"));
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
		   }
		   
	   public static Component snapShot(Gene gene) {
	   		SpliceOSigner pic= new SpliceOSigner(gene);
			pic.setSize(pic.getPreferredSize());
			return pic;
	   }

	   public static Component snapShot(Gene gene, AbstractRegion highlightArea) {
   		SpliceOSigner pic= new SpliceOSigner(gene);
   		pic.setHighlightRegions(new AbstractRegion[] {highlightArea});
		pic.setSize(pic.getPreferredSize());
		return pic;
	   }
	   
	   public static boolean writeGIF(Component component, File f) throws IOException
	   {
	      ExportFileType t= (ExportFileType) ExportFileType.getExportFileTypes("gif").get(0);
	      Properties props= null;
	      String creator= null;
	      t.exportToFile(f,component,component.getParent(),props,creator);
//	      props.put(SAVE_AS_FILE,file.getText());
//	      props.put(SAVE_AS_TYPE,currentType().getFileFilter().getDescription());
	      return true;
	   }
}
