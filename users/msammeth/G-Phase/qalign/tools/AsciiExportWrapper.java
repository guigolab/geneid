package qalign.tools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import javax.swing.table.TableModel;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class AsciiExportWrapper {

	String fileName= null;
	TableModel tableModel= null;

	/**
	 * Constructor for AsciiExportWrapper.
	 */
	public AsciiExportWrapper(String fName, TableModel model) {
		
		this.fileName= fName;
		this.tableModel= model;
	}

	public static void main(String[] args) {
	}
	
	public void writeOut() throws Exception {
		
		BufferedWriter buffy= new BufferedWriter(new FileWriter(fileName));
		
		for (int i= 0; i< tableModel.getRowCount(); ++i) {
			for (int j= 0; j< (tableModel.getColumnCount()-1); ++j)
				if (tableModel.getValueAt(i,j)!= null)
					buffy.write(tableModel.getValueAt(i,j)+ "\t");
				else 
					buffy.write("\t");
			buffy.write(tableModel.getValueAt(i,(tableModel.getColumnCount()- 1))+ "\n");
		}
		
		buffy.flush();
		buffy.close();
	}
}
