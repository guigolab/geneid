/*
 * Created on Nov 26, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.model.copy;

import java.util.Date;
import java.util.Vector;

/**
+--------------------------+---------------+------------------+----------------+--------------+--------------------+------------------------+
| glook_encode_region_name | filt_chr_name | filt_chrom_start | filt_chrom_end | silent_type  | encode_description | glook_encode_region    |
+--------------------------+---------------+------------------+----------------+--------------+--------------------+------------------------+
| ENm001                   | 7             |        115404471 |      117281897 | manual_picks | CFTR               | 7:115404471:117281897  |
| ENm002                   | 5             |        131284313 |      132284313 | manual_picks | Interleukin        | 5:131284313:132284313  |
| ENm003                   | 11            |        115962315 |      116462315 | manual_picks | Apo                | 11:115962315:116462315 |
| ENm004                   | 22            |         30128507 |       31828507 | manual_picks | Chr22              | 22:30128507:31828507   |
| ENm005                   | 21            |         32668236 |       34364221 | manual_picks | Chr21              | 21:32668236:34364221   |
| ENm006                   | X             |        152635144 |      153973591 | manual_picks | ChrX               | X:152635144:153973591  |
| ENm007                   | 19            |         59023584 |       60024460 | manual_picks | Chr19              | 19:59023584:60024460   |
| ENm008                   | 16            |                1 |         500000 | manual_picks | Alpha              | 16:1:500000            |
| ENm009                   | 11            |          4730995 |        5732587 | manual_picks | Beta               | 11:4730995:5732587     |
| ENm010                   | 7             |         26730760 |       27230760 | manual_picks | HOXA               | 7:26730760:27230760    |
| ENm011                   | 11            |          1699991 |        2306039 | manual_picks | 1GF2/H19           | 11:1699991:2306039     |
| ENm012                   | 7             |        113527083 |      114527083 | manual_picks | FOXP2              | 7:113527083:114527083  |
| ENm013                   | 7             |         89428339 |       90542763 | manual_picks | NULL               | 7:89428339:90542763    |
| ENm014                   | 7             |        125672606 |      126835803 | manual_picks | NULL               | 7:125672606:126835803  |
| ENr111                   | 13            |         29418015 |       29918015 | random_picks | NULL               | 13:29418015:29918015   |
| ENr112                   | 2             |         51570355 |       52070355 | random_picks | NULL               | 2:51570355:52070355    |
| ENr113                   | 4             |        118604258 |      119104258 | random_picks | NULL               | 4:118604258:119104258  |
| ENr114                   | 10            |         55153818 |       55653818 | random_picks | NULL               | 10:55153818:55653818   |
| ENr121                   | 2             |        118010803 |      118510803 | random_picks | NULL               | 2:118010803:118510803  |
| ENr122                   | 18            |         59412300 |       59912300 | random_picks | NULL               | 18:59412300:59912300   |
| ENr123                   | 12            |         38626476 |       39126476 | random_picks | NULL               | 12:38626476:39126476   |
| ENr131                   | 2             |        234273824 |      234773888 | random_picks | NULL               | 2:234273824:234773888  |
| ENr132                   | 13            |        112338064 |      112838064 | random_picks | NULL               | 13:112338064:112838064 |
| ENr133                   | 21            |         39244466 |       39744466 | random_picks | NULL               | 21:39244466:39744466   |
| ENr211                   | 16            |         25780427 |       26280428 | random_picks | NULL               | 16:25780427:26280428   |
| ENr212                   | 5             |        141880150 |      142380150 | random_picks | NULL               | 5:141880150:142380150  |
| ENr213                   | 18            |         23719231 |       24219231 | random_picks | NULL               | 18:23719231:24219231   |
| ENr221                   | 5             |         55871006 |       56371006 | random_picks | NULL               | 5:55871006:56371006    |
| ENr222                   | 6             |        132218539 |      132718539 | random_picks | NULL               | 6:132218539:132718539  |
| ENr223                   | 6             |         73789952 |       74289952 | random_picks | NULL               | 6:73789952:74289952    |
| ENr231                   | 1             |        147971133 |      148471133 | random_picks | NULL               | 1:147971133:148471133  |
| ENr232                   | 9             |        128764855 |      129264855 | random_picks | NULL               | 9:128764855:129264855  |
| ENr233                   | 15            |         41520088 |       42020088 | random_picks | NULL               | 15:41520088:42020088   |
| ENr311                   | 14            |         52947075 |       53447075 | random_picks | NULL               | 14:52947075:53447075   |
| ENr312                   | 11            |        130604797 |      131104797 | random_picks | NULL               | 11:130604797:131104797 |
| ENr313                   | 16            |         60833949 |       61333949 | random_picks | NULL               | 16:60833949:61333949   |
| ENr321                   | 8             |        118882220 |      119382220 | random_picks | NULL               | 8:118882220:119382220  |
| ENr322                   | 14            |         98458223 |       98958223 | random_picks | NULL               | 14:98458223:98958223   |
| ENr323                   | 6             |        108371396 |      108871396 | random_picks | NULL               | 6:108371396:108871396  |
| ENr324                   | X             |        122507849 |      123007849 | random_picks | NULL               | X:122507849:123007849  |
| ENr331                   | 2             |        220102850 |      220602850 | random_picks | NULL               | 2:220102850:220602850  |
| ENr332                   | 11            |         63940888 |       64440888 | random_picks | NULL               | 11:63940888:64440888   |
| ENr333                   | 20            |         33304928 |       33804928 | random_picks | NULL               | 20:33304928:33804928   |
| ENr334                   | 6             |         41405894 |       41905894 | random_picks | NULL               | 6:41405894:41905894    |
+--------------------------+---------------+------------------+----------------+--------------+--------------------+------------------------+
 * 
 * @author msammeth
 */
public class EncodeRegion extends DefaultRegion {
	
	public static final String ENC_REG_HUMAN= "regions_coords.txt";
	
	int number= -1;
	String description= null;
	boolean silentTypeManual= false;	
	
	EncodeFragment[] fragments= null;
	
	public static final EncodeRegion[] toArray(Vector v) {
		
		EncodeRegion[] regions= new EncodeRegion[v.size()];
		for (int i = 0; i < v.size(); i++) 
			regions[i]= (EncodeRegion) v.elementAt(i);
		
		return regions;
	}
	
	public EncodeRegion(Species spe) {
		this.species= spe;
	}
	
	public EncodeRegion(Species spe, int newNr) {
		this(spe);
		setNumber(newNr);
	}	

	/**
	 * checks if a new fragment is already in the <code>framgents</code> or not and 
	 * in case, adds it.
	 * 
	 * @param newFragment
	 * @return <code>true</code> if the new fragment was successfully added
	 */
	public boolean addFragment(EncodeFragment newFragment) {
		
		if(newFragment== null)
			return false;
		if (fragments== null) {
			fragments= new EncodeFragment[] {newFragment};
		}
		
		for (int i = 0; i < fragments.length; i++) 		// check if already in fragment[]
			if (fragments[i].equals(newFragment))
				return false;
		
		
		EncodeFragment[] newFragments= new EncodeFragment[fragments.length+ 1];
		for (int i = 0; i < fragments.length; i++) 
			newFragments[i]= fragments[i];
		newFragments[fragments.length]= newFragment;
		
		fragments= newFragments;
		return true;
	}
	public String getName() {
		
		String digits= Integer.toString(number);
		for (; digits.length() < 3; ) 
			digits= "0"+ digits;
		return "EN"+ (silentTypeManual?"m":"r")+digits;
	}
		
	public String toString() {
		return getName()+ ":"+ super.toString();
	}
	/**
	 * @return Returns the description.
	 */
	public String getDescription() {
		return description;
	}
	/**
	 * @param description The description to set.
	 */
	public void setDescription(String description) {
		this.description = description;
	}
	/**
	 * @return Returns the number.
	 */
	public int getNumber() {
		return number;
	}
	/**
	 * @param number The number to set.
	 */
	public void setNumber(int number) {
		this.number = number;
	}
	/**
	 * @return Returns the silent_type_manual.
	 */
	public boolean isSilentTypeManual() {
		return silentTypeManual;
	}
	/**
	 * @param silent_type_manual The silent_type_manual to set.
	 */
	public void setSilentTypeManual(boolean silent_type_manual) {
		this.silentTypeManual = silent_type_manual;
	}
	/**
	 * takes input of the form "ENm014", "ENr111", ...
	 * @param name
	 */
	public void setSilentAndNumber(String name) {
		
		name.trim();
		if (name.startsWith("EN")|| name.startsWith("en"))
			name= name.substring(2, name.length());
		else
			System.err.println("Incomplete/invalid name "+ name);
		
		if (name.charAt(0)== 'm'|| name.charAt(0)== 'M')
			silentTypeManual= true;
		else if (name.charAt(0)== 'r'|| name.charAt(0)== 'R')
			silentTypeManual= false;
		else
			System.err.println("Illegal silent type "+ name.charAt(0));
		
		number= Integer.parseInt(name.substring(1));
	}
	public EncodeFragment[] getFragments() {
		return fragments;
	}
	Date assemblyDate = null;

	String assemblyId = null;

	String assemblyProvider;

	Date freezeDate = null;

	int ncbiTaxonId = -1;
	public Date getAssemblyDate() {
		return assemblyDate;
	}
	public void setAssemblyDate(Date assemblyDate) {
		this.assemblyDate = assemblyDate;
	}
	public String getAssemblyId() {
		return assemblyId;
	}
	public void setAssemblyId(String assemblyId) {
		this.assemblyId = assemblyId;
	}
	public String getAssemblyProvider() {
		return assemblyProvider;
	}
	public void setAssemblyProvider(String assemblyProvider) {
		this.assemblyProvider = assemblyProvider;
	}
	public Date getFreezeDate() {
		return freezeDate;
	}
	public void setFreezeDate(Date freezeDate) {
		this.freezeDate = freezeDate;
	}
	public int getNcbiTaxonId() {
		return ncbiTaxonId;
	}
	public void setNcbiTaxonId(int ncbiTaxonId) {
		this.ncbiTaxonId = ncbiTaxonId;
	}
}
