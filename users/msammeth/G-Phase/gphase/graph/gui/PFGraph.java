package gphase.graph.gui;

import java.util.Arrays;
import java.util.Comparator;

import gphase.graph.SpliceEdge;
import gphase.graph.SpliceGraph;
import gphase.graph.SpliceNode;
import gphase.model.AbstractSite;
import gphase.model.SpliceSite;
import prefuse.data.Graph;
import prefuse.data.Table;
import prefuse.data.column.Column;
import prefuse.data.column.IntColumn;
import prefuse.data.column.ObjectColumn;

public class PFGraph extends Graph {
	public static Table getTable(SpliceNode[] nodes) {
		Table t= new Table(nodes.length, 1);
//		t.addColumn(null, String.class);
//		ObjectColumn c= (ObjectColumn) t.getColumn(0);
//		for (int i = 0; i < nodes.length; i++) 
//			c.set(nodes[i].getSite().toString(), i);	// nodes[i].getSite().getPos()

		return t;
	}
	
	public static Table getTable(SpliceNode[] nodes, SpliceEdge[] edges) {
		Table t= new Table(edges.length, 2);
		t.addColumn(Graph.DEFAULT_SOURCE_KEY, int.class);
		t.addColumn(Graph.DEFAULT_TARGET_KEY, int.class);
		IntColumn cs= (IntColumn) t.getColumn(0);
		IntColumn ct= (IntColumn) t.getColumn(1);
		Comparator compi= new SpliceGraph.NodeOrderComparator();
		for (int i = 0; i < edges.length; i++) {			
			cs.setInt(Arrays.binarySearch(nodes, edges[i].getTail(), compi), i);
			ct.setInt(Arrays.binarySearch(nodes, edges[i].getHead(), compi), i);
		}

		return t;
	}
	
	public PFGraph(SpliceGraph g) {
		super(getTable(g.getNodeList()), getTable(g.getNodeList(), g.getEdgeList()), true);
	}
}
