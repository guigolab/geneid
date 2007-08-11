package gphase.io.bed;

import gphase.tools.Arrays;

import java.awt.Color;
import java.util.Vector;

public class BEDobject {
	int start= -1, end= -1, strand= 0, thickStart= -1, thickEnd= -1, blockCount= 0;
	int score= -1;
	String name= null, chrom= null;
	Color col= null;
	int[] blockSizes= null, blockStarts= null;
	
	public BEDobject(String chromName, int newStrand) {
		setChrom(chromName);
		setStrand(newStrand);
	}
	public BEDobject() {
		
	}
	public BEDobject(String chromName, int newStrand, int newStart, int newEnd) {
		this(chromName, newStrand);
		setStart(newStart);
		setEnd(newEnd);
	}

	public int getBlockCount() {
		return blockCount;
	}

	public void setBlockCount(int blockCount) {
		this.blockCount = blockCount;
		if ((blockStarts!= null&& blockStarts.length!= blockCount)||
				(blockSizes!= null&& blockSizes.length!= blockCount))
			System.out.println("WARNING: block count doesnt match.");
	}

	public int[] getBlockSizes() {
		return blockSizes;
	}

	public void setBlockSizes(int[] blockSizes) {
		this.blockSizes = blockSizes;
		if ((blockStarts!= null&& blockStarts.length!= blockSizes.length)||
				(blockCount> 0&& blockSizes.length!= blockCount))
			System.out.println("WARNING: block sizes doesnt match.");
	}
	
	public void addBlockSize(int blockSize) {
		blockSizes= Arrays.add(blockSizes, blockSize);
	}
	
	public void addBlockStart(int blockStart) {
		blockStarts= Arrays.add(blockStarts, Math.abs(blockStart));
	}
	
	public void addBlockStart(int blockStart, int blockSize) {
		addBlockStart(blockStart);
		addBlockSize(blockSize);
		Vector v= new Vector();
		v.add(blockSizes);
		Arrays.synchroneousSort((Object) blockStarts, v);
	}
	

	public int[] getBlockStarts() {
		return blockStarts;
	}

	public void setBlockStarts(int[] blockStarts) {
		for (int i = 0; i < blockStarts.length; i++) {
			blockStarts[i]= Math.abs(blockStarts[i]);
		}
		this.blockStarts = blockStarts;
		if ((blockSizes!= null&& blockSizes.length!= blockStarts.length)||
				(blockCount> 0&& blockStarts.length!= blockCount))
			System.out.println("WARNING: block starts doesnt match.");
	}

	public String getChrom() {
		return chrom;
	}

	public void setChrom(String chrom) {
		this.chrom = chrom;
	}

	public Color getCol() {
		return col;
	}

	public void setCol(Color col) {
		this.col = col;
	}
	
	public void setCol(int red, int green, int blue) {
		setCol(new Color(red, green, blue));
	}
	
	public void setCol(String rgbVal) {
		int[] tokens= parseCommaSeparatedInts(rgbVal);
		if (tokens.length!= 3) {
			if (!rgbVal.equals("0"))
				System.out.println("WARNING: invalid color "+rgbVal);
			return;
		}
		setCol(tokens[0],tokens[1],tokens[2]);
	}
	
	int[] parseCommaSeparatedInts(String in) {
		String[] tokens= in.split(",");
		int[] out= new int[tokens.length];
		for (int i = 0; i < out.length; i++) 
			out[i]= Integer.parseInt(tokens[i]);
		return out;
	}
	
	public void setBlockStarts(String in) {
		setBlockStarts(parseCommaSeparatedInts(in));
	}
	
	public void setBlockSizes(String in) {
		setBlockSizes(parseCommaSeparatedInts(in));
	}
	

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = Math.abs(end);
		if (start>= 0&& start> this.end) {
			int h= start;
			start= this.end;
			this.end= h;
		}
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public int getScore() {
		return score;
	}

	public void setScore(int score) {
		this.score = score;
	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = Math.abs(start);
		if (end>= 0&& this.start> end) {
			int h= this.start;
			this.start= end;
			end= h;
		}
			
	}

	public int getStrand() {
		return strand;
	}

	public void setStrand(int strand) {
		if (strand!= 1&& strand!= -1)
			System.out.println("WARNING: No strand assignment for "+this);
		
		this.strand = strand;
		
		if (strand< 0) {	// check for swap
			start= -1* Math.abs(getStart());
			end= -1* Math.abs(getEnd());
			if (start> end) {
				int h= start;
				start= end;
				end= h;
			}
		}
	}
	
	public void setStrand(String str) {
		if (str.equals("+"))
			setStrand(1);
		else if (str.equals("-"))
			setStrand(-1);
		else
			System.out.println("WARNING: invalid strand tag "+str);
	}

	public int getThickEnd() {
		return thickEnd;
	}

	public void setThickEnd(int thickEnd) {
		thickEnd= Math.abs(thickEnd);
		this.thickEnd = thickEnd;
		if (thickStart>= 0&& thickStart> this.thickEnd) {
			int h= thickStart;
			thickStart= this.thickEnd;
			this.thickEnd= h;
		}
	}

	public int getThickStart() {
		return thickStart;
	}

	public void setThickStart(int thickStart) {
		thickStart= Math.abs(thickStart);
		this.thickStart = thickStart;
		if (thickEnd>= 0&& this.thickStart> thickEnd) {
			int h= this.thickStart;
			this.thickStart= thickEnd;
			thickEnd= h;
		}
	}
	
	@Override
	public String toString() {
		String s= getChrom()+"\t"+getStart()+"\t"+getEnd();
		String extra= "";
		if (getName()!= null)
			extra+= getName()+"\t";
		else
			extra+= "\t";
		if (getScore()>= 0)
			extra+= getScore()+"\t";
		else
			extra+= ".\t";
		if (getStrand()== 1)
			extra+="+\t";
		else if (getStrand()== -1)
			extra+= "-\t";
		else
			extra+= ".\t";
		if (thickStart>= 0)
			extra+= thickStart+ "\t";
		else
			extra+= "\t";
		if (thickEnd>= 0)
			extra+= thickEnd+ "\t";
		else
			extra+= "\t";
		if (getCol()!= null)
			extra+= getCol().getRed()+","+getCol().getGreen()+","+getCol().getBlue()+"\t";
		else
			extra+= "\t";
		extra+= blockCount+ "\t";
		for (int i = 0; blockSizes!= null&& i < blockSizes.length; i++) {
			extra+= blockSizes[i];
			if (i< blockSizes.length-1)
				extra+= ",";
		}
		extra+= "\t";
		for (int i = 0; blockStarts!= null&& i < blockStarts.length; i++) {
			extra+= blockStarts[i];
			if (i< blockStarts.length-1)
				extra+= ",";
		}
			
		return s+"\t"+extra;
	}
}
