package gphase.tools;

import java.util.Date;

public class Time {

	long time= -1l;

	public Time(long newTime) {
		this.time= newTime;
	}
	
	public int getHours() {
		return (int) (time/ 3600000);
	}
	
	private long getRestHours() {
		return time% 3600000;
	}
	
	public int getMinutes() {
		return (int) (getRestHours()/ 60000);
	}

	private long getRestMinutes() {
		return getRestHours()% 60000;
	}
	
	public int getSeconds() {
		return (int) (getRestMinutes()/ 1000);
	}

	public long getMillis() {
		return getRestMinutes()% 1000;
	}
	
	public String toString() {
		return getHours()+ ":"+ getMinutes()+":"+ getSeconds()+"."+getMillis();
	}
	
	public static String toTime(long temps) {
		Time t= new Time(temps);
		return t.toString();
	}
}
