package qalign.tools;

import java.util.Comparator;
import java.util.Date;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class SequenceWrapper {

	protected boolean selected= false;
	protected String name= null;
	protected String comment= null;
	protected String fName= null;
	protected Date fDate= null;
	protected String sequence= null;
	protected IOWrapper wrapper= null;
	
	public static Comparator getNameComparator() {
		
		Comparator comp= new Comparator() {
			
			public boolean equals(Object o1, Object o2) {
				
				if (((SequenceWrapper) o1).getName().equals(((SequenceWrapper) o2).getName()))
					return true;
				return false;
			}
			
			public int compare(Object o1, Object o2) {
			
				return (((SequenceWrapper) o1).getName().compareTo(((SequenceWrapper) o2).getName()));
			}
		};
		
		return comp;
	}

	public static Comparator getFNameComparator() {
		
		Comparator comp= new Comparator() {
			
			public boolean equals(Object o1, Object o2) {
				
				if (((SequenceWrapper) o1).getFName().equals(((SequenceWrapper) o2).getFName()))
					return true;
				return false;
			}
			
			public int compare(Object o1, Object o2) {
			
				return (((SequenceWrapper) o1).getFName().compareTo(((SequenceWrapper) o2).getFName()));
			}
		};
		
		return comp;
	}

	public static Comparator getFDateComparator() {
		
		Comparator comp= new Comparator() {
			
			public boolean equals(Object o1, Object o2) {
				
				if (((SequenceWrapper) o1).getFDate().equals(((SequenceWrapper) o2).getFDate()))
					return true;
				return false;
			}
			
			public int compare(Object o1, Object o2) {
			
				return (((SequenceWrapper) o1).getFDate().compareTo(((SequenceWrapper) o2).getFDate()));
			}
		};
		
		return comp;
	}
	
	public SequenceWrapper(String newName, String newComment, String newFName, Date newFDate) {
		
		setName(newName);
		setComment(newComment);
		setFName(newFName);
		setFDate(newFDate);
	}
	
	public SequenceWrapper() {
	}
	
	public int compareTo(Object o2) throws ClassCastException {
		
		if (getSequence()== null)
			System.out.println(getName());
			
		if ((o2== null)|| (! (o2 instanceof SequenceWrapper))
			|| (((SequenceWrapper) o2).getSequence()== null)
			|| (getSequence()== null))
			throw new ClassCastException();
		
		return getSequence().compareTo(((SequenceWrapper) o2).getSequence());
	}
	
	public boolean equals(Object o2) {
		
		try {
			if (compareTo(o2)== 0)
				return true;
		} catch (ClassCastException e) {
			return false;
		}
		
		return false;
	}
	
	public static SequenceWrapper[] extend(SequenceWrapper[] sw1, SequenceWrapper[] sw2) {
		
			// nothing to do
		if (sw1== null)
			return sw2;
		if (sw2== null)
			return sw1;
		
			// merge
		Vector result= new Vector(sw1.length+ sw2.length);
		for (int i= 0; i< sw1.length; ++i)
			result.add(sw1[i]);
		for (int i= 0; i< sw2.length; ++i)
			result.add(sw2[i]);
			
			// kill doubles
		for (int i= 0; i< result.size(); ++i)
			for (int j= i+1; j< result.size(); ++j)
				if (result.elementAt(i).equals(result.elementAt(j)))
					result.remove(j--);
		
		SequenceWrapper[] res= new SequenceWrapper[result.size()];
		for (int i= 0; i< res.length; ++i)
			res[i]= (SequenceWrapper) result.elementAt(i);
			
		return res;
	}
	
	public String[] getAll() {
		
		String[] result= new String[4];
		result[0]= getName();
		result[1]= getComment();
		result[2]= getFName();

		String bDate= getFDate().toString();
		StringTokenizer st= new StringTokenizer(bDate, " ");
		String date= "";
		String time= null;
		st.nextToken();	// dow
		date+= st.nextToken()+"/";
		date+= st.nextToken()+"/";
		time= st.nextToken();
		time= time.substring(0, time.lastIndexOf(':'));
		st.nextToken();	// time zone
		date+= st.nextToken();
		
		result[3]= date+ " "+ time;
		
		return result;
	}
	
	/**
	 * Returns the comment.
	 * @return String
	 */
	public String getComment() {
		return comment;
	}

	/**
	 * Returns the name.
	 * @return String
	 */
	public String getName() {
		return name;
	}

	/**
	 * Returns the wrapper.
	 * @return IOWrapper
	 */
	public IOWrapper getWrapper() {
		return wrapper;
	}

	/**
	 * Sets the comment.
	 * @param comment The comment to set
	 */
	public void setComment(String comment) {
		this.comment = comment;
	}

	/**
	 * Sets the name.
	 * @param name The name to set
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Sets the wrapper.
	 * @param wrapper The wrapper to set
	 */
	public void setWrapper(IOWrapper wrapper) {
		this.wrapper = wrapper;
	}

	/**
	 * Returns the fDate.
	 * @return String
	 */
	public Date getFDate() {
		return fDate;
	}

	/**
	 * Returns the fName.
	 * @return String
	 */
	public String getFName() {
		return fName;
	}

	/**
	 * Sets the fDate.
	 * @param fDate The fDate to set
	 */
	public void setFDate(Date fDate) {
		this.fDate = fDate;
	}

	/**
	 * Sets the fName.
	 * @param fName The fName to set
	 */
	public void setFName(String fName) {
		this.fName = fName;
	}

	/**
	 * Returns the sequence.
	 * @return String
	 */
	public String getSequence() {
		return sequence;
	}

	/**
	 * Sets the sequence.
	 * @param sequence The sequence to set
	 */
	public void setSequence(String sequence) {
		this.sequence = sequence.toUpperCase();
	}

	/**
	 * Returns the selected.
	 * @return boolean
	 */
	public boolean isSelected() {
		return selected;
	}

	/**
	 * Sets the selected.
	 * @param selected The selected to set
	 */
	public void setSelected(boolean selected) {
		this.selected = selected;
	}

}
