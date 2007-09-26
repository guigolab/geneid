package gphase.io;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Vector;

import gphase.tools.Arrays;
import gphase.tools.File;

public class ASTAcomparator {
	public static void main(String[] args) {
		ASTAcomparator compi= new ASTAcomparator();
		parseArguments(compi, args);
		if (compi.ready())
			compi.run(System.out);		
	}
	
	static void parseArguments(ASTAcomparator compi, String[] args) {
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("-s")|| args[i].equalsIgnoreCase("--struct")) {
				compi.setStruct(args[++i].split("&")); 
				continue;
			}
			if (args[i].equalsIgnoreCase("-n")|| args[i].equalsIgnoreCase("--notFound")) {
				compi.setOutputNotFound(true);
				continue;
			}

			if (args[i].equalsIgnoreCase("-in1")|| args[i].equalsIgnoreCase("--input1")) {
				File f= new File(args[++i]);
				if (!f.exists())
					System.err.println("File "+f+" not found!");
				else
					compi.setF1(f);
				continue;
			}
			if (args[i].equalsIgnoreCase("-in2")|| args[i].equalsIgnoreCase("--input2")) {
				File f= new File(args[++i]);
				if (!f.exists())
					System.err.println("File "+f+" not found!");
				else
					compi.setF2(f);
				continue;
			}
			if (args[i].equalsIgnoreCase("-st")|| args[i].equalsIgnoreCase("--strict")) {
				compi.strictChk= true;
				continue;
			}
			if (args[i].equalsIgnoreCase("-ro")|| args[i].equalsIgnoreCase("--refOnly")) {
				compi.diffInNew= false;
			}
			if (args[i].equalsIgnoreCase("-to")|| args[i].equalsIgnoreCase("--targetOnly")) {
				compi.diffInRef= false;
			}
		}
	}
	
	File f1= null, f2= null;
	String[] struct= null;
	boolean outputNotFound= false, strictChk= false, diffInRef= true, diffInNew= true;
	
	public ASTAcomparator() {
	}
	
	public boolean ready() {
		if (getF1()!= null|| getF2()!= null&& getStruct()!= null)
			return true;
		return false;
	}
	
	public void run(PrintStream p) {

		String[] structs= readStructures(getStruct(), getF1());
		HashMap foundMap= findStructures(structs, getF2());
		removeStructuresFromInput(foundMap, getF1());
		p.println("Investigated structures");
		for (int i = 0; i < getStruct().length; i++) 
			p.print(getStruct()[i]+"&");
		p.println("\n=======================\n\n");
		
		if (outputNotFound) {
			String[] nFound1= getNotFoundStructs2(foundMap, getF1(), true);
			String[] nFound2= getNotFoundStructs2(foundMap, getF2(), false);
			
			if (diffInRef) {
				p.println("NOT found from REFERENCE file "+getF1());
				p.println("===========================");
				for (int i = 0; nFound1!= null&& i < nFound1.length; i++) 
					p.println(nFound1[i]);
			}
			if (diffInNew)  {
				p.println();
				p.println("NOT found in COMPARISON file "+getF2());
				p.println("===========================");
				for (int i = 0; nFound2!= null&& i < nFound2.length; i++) 
					p.println(nFound2[i]);
			}			
		} else {
			Object[] keys= foundMap.keySet().toArray();
			if (diffInRef) {
				p.println("FOUND from REFERENCE file "+getF1());
				p.println("========================");
				for (int i = 0; i < keys.length; i++) 
					p.println(keys[i]);
			}
			if (diffInNew) {
				p.println();
				p.println("FOUND in COMPARISON file "+getF2());
				p.println("===========================");
				for (int i = 0; i < keys.length; i++) 
					p.println(foundMap.get(keys[i]));
			}
		}
		
	}
	
	private void removeStructuresFromInput(HashMap foundMap, File f) {
		
		HashMap coordMap= new HashMap();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(f));
			while (buffy.ready()) {
				String line= buffy.readLine();
				String[] tokens= line.split("\t");				
				int i;
				for (i = 0; i < getStruct().length; i++) {
					if (tokens[0].equals(getStruct()[i]))
						break;
				}
				if (i< getStruct().length)
					continue;
				String code= tokens[3]+"|";
				if (tokens.length> 5)
					code+= tokens[5];
				coordMap.put(code, code);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		
		Object[] keys= foundMap.keySet().toArray();
		for (int i = 0; i < keys.length; i++) {
			Vector v= (Vector) foundMap.get(keys[i]);
			for (int j = 0; j < v.size(); j++) {
				String[] tokens= ((String) v.elementAt(j)).split("\t");
				String code= tokens[3]+"|";
				if (tokens.length> 5)
					code+= tokens[5];
				if (coordMap.get(code)!= null)
					v.remove(j--);				
			}
			if (v.size()== 0)
				foundMap.remove(keys[i]);
		}
		
	}

	private String[] getNotFoundStructs(HashMap foundMap, File f, boolean searchKeys) {
			
			HashMap<String,Vector> structMap= new HashMap<String,Vector>();
			Object[] keys= foundMap.keySet().toArray();
			for (int i = 0; i < keys.length; i++) {
				String[] line= null;
				if (searchKeys) 
					line= new String[] {(String) keys[i]};
				else {
					Vector v= (Vector) foundMap.get(keys[i]);
					// doesnt help
	//				String baseStruct= ((String) keys[i]).split("\t")[0];
	//				for (int j = 0; j < v.size(); j++) {	// exact hits preferred
	//					String s= (String) v.elementAt(j);
	//					if (s.split("\t").equals(baseStruct)) {
	//						v= new Vector();
	//						v.add(s);
	//						break;
	//					}
	//				}
					line= (String[]) Arrays.toField(v);
				}
				for (int j = 0; j < line.length; j++) {
					String[] tokens= line[j].split("\t");
					Vector v= structMap.get(tokens[0]);
					if (v== null)
						v= new Vector();
					v.add(line[j]);
					structMap.put(tokens[0],v);
				}
				
			}
			
			Vector v= new Vector();
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				while (buffy.ready()) {
					String line= buffy.readLine();
					String[] tokens= line.split("\t");
					if (structMap.get(tokens[0])== null)
						continue;
					Vector foundV= (Vector) structMap.get(tokens[0]);
					int i;
					for (i = 0; i < foundV.size(); i++) {
						if (line.equalsIgnoreCase((String) foundV.elementAt(i)))
							break;
					}
					if (i== foundV.size())
						v.add(line);
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return (String[]) Arrays.toField(v);
		}

	private String[] getNotFoundStructs2(HashMap foundMap, File f, boolean searchKeys) {
			
			HashMap<String,Vector> structMap= new HashMap<String,Vector>();
			Object[] keys= foundMap.keySet().toArray();
			Vector v= new Vector();
			
			try {
				BufferedReader buffy= new BufferedReader(new FileReader(f));
				while (buffy.ready()) {
					String line= buffy.readLine();
					boolean found= false;
					for (int i = 0; i < keys.length; i++) {
						String[] token= null;
						if (searchKeys) 
							token= new String[] {(String) keys[i]};
						else {
							Vector vv= (Vector) foundMap.get(keys[i]);
							token= (String[]) Arrays.toField(vv);
						}
						for (int j = 0; j < token.length; j++) {
							if (token[j].equalsIgnoreCase(line)) {
								found= true;
								break;
							}
						}
						if (found)
							break;
						
					}
					if (!found) {
						if ((!strictChk)|| (strictChk&& line.startsWith(struct[0])))
							v.add(line);
					}
				}
			} catch (Exception e) {
				e.printStackTrace();
			}
			
			return (String[]) Arrays.toField(v);
		}

	private HashMap findStructures(String[] structs, File f) {
		
		Vector<String> coordV1= new Vector<String>(structs.length),
			coordV2= new Vector<String>(structs.length);
		for (int i = 0; i < structs.length; i++) {
			String[] tokens= structs[i].split("\t");
			coordV1.add(tokens[3]);
			if (tokens.length> 5)
				coordV2.add(tokens[5]);
			else
				coordV2.add(null);
		}
		
		HashMap foundMap= new HashMap();
		try {
			BufferedReader buffy = new BufferedReader(new FileReader(f));
			while (buffy.ready()) {
				String line= buffy.readLine();
				String[] tokens= line.split("\t");
				
				for (int i = 0; i < coordV1.size(); i++) {
					if (strictChk) {
						if (tokens[3].equals(coordV1.elementAt(i))) {
							if (((tokens.length< 6|| tokens[5].length()== 0)&& coordV2.elementAt(i)== null)
									||tokens[5].equals(coordV2.elementAt(i))) {
								Vector v= (Vector) foundMap.get(structs[i]);
								if (v== null)
									v= new Vector();
								v.add(line);
								foundMap.put(structs[i], v);
							}
						}
						continue;
					}
					
//					if (tokens[3].indexOf(coordV1.elementAt(i))>= 0|| (coordV2.elementAt(i)!= null&& tokens[3].indexOf(coordV2.elementAt(i))>= 0) ||
//							tokens.length> 5&& (tokens[5].indexOf(coordV1.elementAt(i))>= 0|| (coordV2.elementAt(i)!= null&& tokens[5].indexOf(coordV2.elementAt(i))>= 0))) {
					if (tokens[3].indexOf(coordV1.elementAt(i))>= 0|| tokens.length> 5&& tokens[5].indexOf(coordV1.elementAt(i))>= 0) {						
						if (coordV2.elementAt(i)!= null) {	// check second
							String otherToken= tokens[3];
							if (tokens[3].indexOf(coordV1.elementAt(i))>= 0) {
								if (tokens.length< 6)
									continue;
								otherToken= tokens[5];
							}
							if (otherToken.indexOf(coordV2.elementAt(i))< 0)
								continue;
						}
						Vector v= (Vector) foundMap.get(structs[i]);
						if (v== null)
							v= new Vector();
						v.add(line);
						foundMap.put(structs[i], v);
//						coordV1.removeElementAt(i);
//						coordV2.removeElementAt(i);
//						break;
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return foundMap;
	}

	private String[] readStructures(String[] struct, File f) {
		Vector v= new Vector();
		try {
			BufferedReader buffy= new BufferedReader(new FileReader(f));			
			while (buffy.ready()) {
				String line= buffy.readLine();
				for (int i = 0; i < struct.length; i++) {
					if (line.startsWith(struct[i]))
						v.add(line);
				}				
			}
		} catch (Exception e) {			
			e.printStackTrace();
		}
		
		return (String[]) Arrays.toField(v);
	}

	private String insertSpaces(String string) {
		String[] tokens= string.split(",");
		
		return tokens[0]+" , "+tokens[1];
	}

	public File getF1() {
		return f1;
	}

	public void setF1(File f1) {
		this.f1 = f1;
	}

	public File getF2() {
		return f2;
	}

	public void setF2(File f2) {
		this.f2 = f2;
	}

	public String[] getStruct() {
		return struct;
	}

	public void setStruct(String[] struct1) {
		this.struct = struct1;
		for (int i = 0; i < struct.length; i++) {
			struct[i]= struct[i];//insertSpaces(struct[i]);
		}
	}

	public boolean isOutputNotFound() {
		return outputNotFound;
	}

	public void setOutputNotFound(boolean outputNotFound) {
		this.outputNotFound = outputNotFound;
	}
}
