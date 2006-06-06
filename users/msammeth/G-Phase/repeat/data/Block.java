/*
 * Created on Nov 27, 2003
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package repeat.data;

/**
 * 
 * 
 * @author micha
 */
public class Block {

	public static final byte ID_AHELIX= 1;
	public static final byte ID_BSTRAND= 2;
	public static final byte ID_CBLOCK= 3;
	public static final byte ID_RAREA= 4;
	
	protected byte ID= 0;
	protected String[] seqID= null;
	protected int[] startPos= {-1};
	protected int length= -1;

	public static void main(String[] args) {
		Block oneBlock= new Block(new String[] {"bla","blubb"}, new int[] {1,2},1);
		Block eqBlock= new Block(new String[] {"bla","blubb"}, new int[] {1,2},1);
		Block uneqBlock= new Block(new String[] {"bla"}, new int[] {1},1);
		System.out.println(oneBlock.equals(eqBlock));
		System.out.println(oneBlock.equals(uneqBlock));
	}

	public Block() {	
	}
	
	public Block(String[] newSeqID, int[] newStartPos, int newLength) {
		
		this();
		this.seqID= newSeqID;
		this.startPos= newStartPos;
		this.length= newLength;
	}
	
	public boolean isAHelix() {
		return isAHelix(this);
	}
	
	public static boolean isAHelix(Block aBlock) {
		return (aBlock.getID()== ID_AHELIX);
	}

	public boolean isBStrand() {
		return isBStrand(this);
	}
	
	public static boolean isBStrand(Block aBlock) {
		return (aBlock.getID()== ID_BSTRAND);
	}
	public boolean isCBlock() {
		return isCBlock(this);
	}
	
	public static boolean isCBlock(Block aBlock) {
		return (aBlock.getID()== ID_CBLOCK);
	}
	public boolean isRArea() {
		return isRArea(this);
	}
	
	public static boolean isRArea(Block aBlock) {
		return (aBlock.getID()== ID_RAREA);
	}
	
	
	/**
	 * @return
	 */
	public byte getID() {
		return ID;
	}

	/**
	 * @return
	 */
	public int getLength() {
		return length;
	}

	/**
	 * @return
	 */
	public String[] getSeqID() {
		return seqID;
	}

	/**
	 * @return
	 */
	public int[] getStartPos() {
		return startPos;
	}

	/**
	 * @param b
	 */
	public void setID(byte b) {
		ID= b;
	}

	/**
	 * @param i
	 */
	public void setLength(int i) {
		length= i;
	}

	/**
	 * @param string
	 */
	public void setSeqID(String[] string) {
		seqID= string;
	}

	/**
	 * @param i
	 */
	public void setStartPos(int[] i) {
		startPos= i;
	}
	
	public boolean equals(Object anObject) {
		
		Block anotherBlock= null;
		try {
			anotherBlock= (Block) anObject;
		} catch (ClassCastException e) {
			return false;
		}
		
		try {
			
				// compare IDs
			if (getID()!= anotherBlock.getID())
				return false;

				// compare seqID
			int[] mapping= null;
			if (seqID!= null) {
				mapping= new int[seqID.length];
				for (int i= 0; i< seqID.length; ++i)
					mapping[i]= i;
					
				for (int i= 0; i< seqID.length; ++i) {
					
					int j;
					for (j= 0; j< seqID.length; ++j) 
						if (seqID[i].equals(anotherBlock.getSeqID()[j]))
							break;
					if (j== seqID.length)
						return false;
					mapping[i]= j;
				}
			} else
				if (anotherBlock.getSeqID()!= null)
					return false;
			
				// compare Positions
			if (startPos!= null) 
				for (int i= 0; i< startPos.length; ++i)
					if (startPos[i]!= anotherBlock.getStartPos()[mapping[i]])
						return false;
			
		} catch (NullPointerException e) {
			return false;
		} catch (ArrayIndexOutOfBoundsException e) {
			return false;
		}
		
		return true;
	}
	
	public String toString() {
		
		String result= "Block: ";
		switch (ID) {
			case ID_AHELIX:
				result+= "aHelix";
				break;
			case(ID_BSTRAND):
				result+= "bStrand";
				break;
			case(ID_CBLOCK):
				result+= "cBlock";
				break;
			case(ID_RAREA):
				result+= "rArea";
				break;
		}
		result+= " {";
		for (int i= 0; i < seqID.length; i++) {
			result+= seqID[i]+ " ("+ startPos[i]+ ","+ (startPos[i]+ length)+ "), ";
		}
		result= result.substring(0, result.length()- 2);
		
		return result;
	}

}
