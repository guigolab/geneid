package gphase.sgraph;

import gphase.model.ASEvent;
import gphase.model.ASVariation;
import gphase.model.DirectedRegion;
import gphase.model.SpliceSite;
import gphase.model.Transcript;

import java.io.FileWriter;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import org.freehep.graphicsio.emf.IntersectClipRect;

public class Bubble {
	static int bubbleID= 0;
	
	Vector<Path> pathes;
	Vector<Path[]> tuples;
	Graph g;
	
	public Bubble(Graph newGraph, Vector<Path> newPathes, Vector<Path[]> newTuples) {
		pathes= newPathes;
		tuples= newTuples;
		g= newGraph;
		for (int i = 0; i < pathes.size(); i++) 
			pathes.elementAt(i).addBubbleID(bubbleID);
		++bubbleID;
		
	}
	
	public void extractEvents(HashMap<String, Integer> map) {
		for (int i = 0; i < tuples.size(); i++) {
			Transcript[][] t= new Transcript[tuples.elementAt(i).length][];
			Vector[] v= new Vector[t.length];
			
			for (int j = 0; j < v.length; j++) {
				v[j]= tuples.elementAt(i)[j].getAllNodes();
				v[j].remove(0);
				v[j].remove(v[j].size()-1);
			}

			SpliceSite[][] ss= new SpliceSite[v.length][];
			for (int j = 0; j < ss.length; j++) {
				ss[j]= new SpliceSite[v[j].size()];
				for (int m = 0; m < ss[j].length; m++) 
					ss[j][m]= ((Node) v[j].elementAt(m)).getSite();
			}
			
			for (int j = 0; j < t.length; j++) 
				t[j]= getGraph().decodeTset(tuples.elementAt(i)[j].getTranscripts());
			
			ASEvent ev= new ASEvent(t, ss);
			try {
				FileWriter f= new FileWriter("new.asta",true);
				f.write(ev.toStringASTA());
				f.write("\n");
				f.flush();
				f.close();
			} catch (Exception e) {
				e.printStackTrace();
			}
			//Vector<ASEvent> vec= map.get(ev.toString());
			Integer cnt= map.get(ev.toString());
			if (cnt== null) 
				cnt= new Integer(0);
			map.put(ev.toString(), new Integer(cnt.intValue()+1));
		}

	}
	
	Graph getGraph() {
		return g;
	}
	
}
