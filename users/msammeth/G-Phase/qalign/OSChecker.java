package qalign;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class OSChecker {

	private static boolean winXP= false;
	private static boolean win2K= false;
	private static boolean winME= false;
	private static boolean winNT= false;
	private static boolean win9x= false;
	private static boolean linux= false;
	private static boolean sunOS= false;
	private static boolean hpOS= false;
	private static boolean macOSX= false;
	
	static {
		init();
	}
	

	private static void init() {
		
		String osName= System.getProperty("os.name").toLowerCase();
		
		if (osName.indexOf("windows xp") > -1) 
			winXP= true;
		if (osName.indexOf("windows 2000") > -1) 
			win2K= true;
		if (osName.indexOf("me") > -1) 
			winME= true;
		if (osName.indexOf("nt") > -1) 
			winNT= true;
		if (osName.indexOf("windows 9") > -1) 
			win9x= true;
		
		if (osName.indexOf("linux") > -1) 		// startsWith()
			linux= true;
		if ((osName.indexOf("sunos") > -1)||
			(osName.indexOf("solaris") > -1)) 	// startsWith()
			sunOS= true;
		if (osName.indexOf("hp") > -1) 		// startsWith()
			hpOS= true;
			
		if (osName.indexOf("mac os x") > -1) 		// startsWith()
			macOSX= true;
	}
	
	
	/**
	 * Returns the linux.
	 * @return boolean
	 */
	public static boolean isLinux() {
		return linux;
	}

	/**
	 * Returns the sunOS.
	 * @return boolean
	 */
	public static boolean isSunOS() {
		return sunOS;
	}

	/**
	 * Returns the win2K.
	 * @return boolean
	 */
	public static boolean isWin2K() {
		return win2K;
	}

	/**
	 * Returns the win9x.
	 * @return boolean
	 */
	public static boolean isWin9x() {
		return win9x;
	}

	/**
	 * Returns the winME.
	 * @return boolean
	 */
	public static boolean isWinME() {
		return winME;
	}

	/**
	 * Returns the winNT.
	 * @return boolean
	 */
	public static boolean isWinNT() {
		return winNT;
	}

	/**
	 * Returns the winXP.
	 * @return boolean
	 */
	public static boolean isWinXP() {
		return winXP;
	}
	
	/**
	 * Returns the winXP.
	 * @return boolean
	 */
	public static boolean isWindows() {
		return (isWinXP()|| isWin2K()|| isWin9x()|| isWinME()|| isWinNT());
	}
	
	/**
	 * Returns the macOSX.
	 * @return boolean
	 */
	public static boolean isMacOSX() {
		return macOSX;
	}	

	/**
	 * Returns the hpOS.
	 * @return boolean
	 */
	public static boolean isHpOS() {
		return hpOS;
	}

}
