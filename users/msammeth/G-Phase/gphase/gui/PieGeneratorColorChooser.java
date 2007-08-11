package gphase.gui;

import gphase.tools.File;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.EventObject;
import java.util.HashMap;
import java.util.Vector;

import javax.annotation.processing.Filer;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JColorChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JWindow;
import javax.swing.ListSelectionModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

public class PieGeneratorColorChooser {
	
	public static class ColorComparator implements Comparator {
		public int compare(Object o1, Object o2) {
			int rgb1= ((Circle) ((Object[]) o1)[0]).getColor().getRGB();
			int rgb2= ((Circle) ((Object[]) o2)[0]).getColor().getRGB();
			
			if (rgb1< rgb2)
				return -1;
			if (rgb1> rgb2)
				return 1;
			return 0;
		}
	}
	
	static class CellComponentRenderer implements TableCellRenderer { 

		  public static final DefaultTableCellRenderer DEFAULT_RENDERER = new DefaultTableCellRenderer();

		  public Component getTableCellRendererComponent(JTable table, Object value,
		      boolean isSelected, boolean hasFocus, int row, int column) {
		    Component renderer = DEFAULT_RENDERER.getTableCellRendererComponent(
		        table, value, isSelected, hasFocus, row, column);
		    if (value instanceof Component)
		    	return ((Component) value);
		    return renderer;
		  }
		
	}
	class SelectionListener implements ListSelectionListener {
    
        // It is necessary to keep the table since it is not possible
        // to determine the table from the event's source
        SelectionListener() {
        }
        public void valueChanged(ListSelectionEvent e) {
            // If cell selection is enabled, both row and column change events are fired
            if (e.getSource() == PieGeneratorColorChooser.this.getTable().getSelectionModel()
                  && PieGeneratorColorChooser.this.getTable().getRowSelectionAllowed()) {
                // Column selection changed
                int first = e.getFirstIndex();
                int last = e.getLastIndex();	// works better :)
                PieGeneratorColorChooser.this.getColChooser().setColor(((Circle) PieGeneratorColorChooser.this.getTable().getValueAt(getTable().getSelectedRow(), 0)).getColor());
            } else if (e.getSource() == PieGeneratorColorChooser.this.getTable().getColumnModel().getSelectionModel()
                   && PieGeneratorColorChooser.this.getTable().getColumnSelectionAllowed() ){
                // Row selection changed
                int first = e.getFirstIndex();
                int last = e.getLastIndex();
            }
    
            if (e.getValueIsAdjusting()) {
            }
        }
        
        
    }
	
	public static void main(String[] args) {
		if (args== null|| args.length!= 1) {
			System.err.println("Hi! I need a input file. Bye.");
			System.exit(-1);
		}
			
		PieGeneratorColorChooser myCM= new PieGeneratorColorChooser(args[0]);
	}
	
	JFrame window= null;
	JColorChooser colChooser= null;
	String fName= null;
	JTable table= null;
	JScrollPane scroller= null;
	Object[][] rowData= null;
	
	public PieGeneratorColorChooser(String newFName) {
		this.fName= newFName;
		window= new JFrame();
		window.getContentPane().setLayout(new BorderLayout());
		window.getContentPane().add(getColChooser(), BorderLayout.NORTH);
		addScroller();
		window.getContentPane().add(getButtons(), BorderLayout.CENTER);
		
		window.pack();
		window.setVisible(true);
		window.addWindowListener(new WindowAdapter() {
			
			public void windowClosing(WindowEvent e) {
				super.windowClosing(e);
				
				writeTable(new File(PieGeneratorColorChooser.this.fName));
				System.out.println("wrote color map.");
				System.exit(0);
			}
		});
	}
	
	private void addScroller() {
		if (scroller!= null)
			window.getContentPane().remove(scroller);
		scroller= new JScrollPane(getTable());
		window.getContentPane().add(scroller, BorderLayout.SOUTH);		
	}

	private Component getButtons() {
		JPanel panel= new JPanel();
		panel.setLayout(new FlowLayout());
		JButton butAdd= new JButton("Add");
		panel.add(butAdd);
		butAdd.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				Object[][] newRows= new Object[rowData.length+ 1][];
				for (int i = 0; i < rowData.length; i++) 
					newRows[i]= rowData[i];
				
				Object[] o= new Object[2];
				o[0]= new Circle(Color.cyan, 10);
				o[1]= "1-2^ , 3-4^";
				newRows[newRows.length- 1]= o;
				
				rowData= newRows;
				table= null;
				addScroller();
				window.setVisible(true);
				
			}
		});
		
		JButton butDel= new JButton("Remove");
		panel.add(butDel);
		butDel.addActionListener(new ActionListener(){
			public void actionPerformed(ActionEvent e) {
				if (getTable().getSelectedRow()< 0|| getTable().getSelectedRow()>= rowData.length)
					return;
				Object[][] newRows= new Object[rowData.length- 1][];
				int ctr= 0;
				for (int i = 0; i < rowData.length; i++) {
					if (i!= getTable().getSelectedRow())
						newRows[ctr++]= rowData[i];
				}
				
				rowData= newRows;
				table= null;
				addScroller();
				window.setVisible(true);
			}
		});

		return panel;
	}

	JTable getTable() {
		if (table== null) {
			if (rowData== null)
				rowData= readTable(new gphase.tools.File(fName));
			Object[] colNames= new String[] {"color", "ID"};
			table= new JTable(rowData, colNames) {
				public boolean isCellEditable(int row, int column) {
					if (column!= 1)
						return false;
					return true;
				}
			};
			table.setDefaultRenderer(Object.class, new CellComponentRenderer());
			
				// listener
			table.getSelectionModel().addListSelectionListener(new SelectionListener());
			table.setSelectionMode(ListSelectionModel.SINGLE_INTERVAL_SELECTION);
		}
		return table;
	} 

	
	
	private Object[][] readTable(gphase.tools.File file) {
		
		if (file== null|| !file.exists())
			return new Object[0][];
		Vector v= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(file));
			while (buffy.ready()) {
				String line= buffy.readLine().trim();
				String[] tokens= line.split("\t");
				if (tokens.length!= 2)
					System.err.println("Invalid line "+line);		
				tokens[1]= tokens[1].substring(2, tokens[1].length()- 1);
				int r= Integer.valueOf(tokens[1].substring(0,2), 16);
				int g= Integer.valueOf(tokens[1].substring(2,4), 16);
				int b= Integer.valueOf(tokens[1].substring(4,6), 16);
				int a= Integer.valueOf(tokens[1].substring(6,8), 16);
				Color c= new Color(r,g,b,a);
				Object[] o= new Object[2];
				o[1]= tokens[0];
				o[0]= new Circle(c, 10);
				v.add(o);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		Object[][] o= new Object[v.size()][];
		for (int i = 0; i < o.length; i++) 
			o[i]= (Object[]) v.elementAt(i);
		
		
		return o;
	}

	private void writeTable(gphase.tools.File file) {
		try {
			BufferedWriter buffy= new BufferedWriter(new FileWriter(file));
			for (int i = 0; rowData!= null&& i < rowData.length; i++) {
				Color col= ((Circle) rowData[i][0]).getColor();
				String colStr= getHexString(col.getRed());
				colStr+= getHexString(col.getGreen());
				colStr+= getHexString(col.getBlue());
				colStr+= getHexString(col.getAlpha());
				String line= rowData[i][1]+"\t\"#"+ colStr+ "\"";
				buffy.write(line+"\n");
			}
			buffy.flush(); buffy.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	private String getHexString(int x) {
		String s= Integer.toHexString(x);
		if (s.length()== 1)
			s= "0"+s;
		return s;
	}

	public JColorChooser getColChooser() {
		if (colChooser == null) {
			colChooser= new JColorChooser();
			colChooser.getSelectionModel().addChangeListener(new ChangeListener() {
				public void stateChanged(ChangeEvent e) {
				    Color newColor = PieGeneratorColorChooser.this.getColChooser().getColor();
				    for (int i = 0; i < getTable().getSelectedRows().length; i++) 
					    getTable().setValueAt(new Circle(newColor, 10), 
					    		getTable().getSelectedRows()[i], 0);
				    getTable().repaint();
				}
			});
			
		}

		return colChooser;
	}
	
	
}
