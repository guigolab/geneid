/*
 * Created on Nov 8, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.db;

import gphase.model.Graph;
import junit.framework.AssertionFailedError;
import junit.framework.TestCase;

/**
 * 
 * 
 * @author msammeth
 */
public class GraphTester extends TestCase {

	Graph graph;
	
	public GraphTester(Graph newGraph) {		
		super(newGraph.toString());	// no standard constructor
		this.graph= newGraph;
	}

	
	public void testGene(String stableID, String sequence) throws AssertionFailedError{
		if (!graph.readSequence(graph.getGene(stableID)).equals(sequence))
			throw new AssertionFailedError("Test gene sequences does not match: "+stableID);				
	}
}
