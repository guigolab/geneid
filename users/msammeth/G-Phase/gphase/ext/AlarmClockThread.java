/*
 * Created on Sep 5, 2005
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */
package gphase.ext;

import java.util.Date;


/**
 * 
 * 
 * @author msammeth
 */
public class AlarmClockThread extends Thread {

	Thread anotherThread;
	Date alarm;
	long check= 60000l;
	
	public void run() {
		Date now= new Date(System.currentTimeMillis());
		while (true) {
			if (!anotherThread.isAlive())
				break;
			if (now.compareTo(alarm)< 0) {
				anotherThread.interrupt();
				break;
			}
			try {Thread.sleep(check);} 
			catch (InterruptedException e) {;}
		}
	}
	
	public AlarmClockThread(Thread aThread, Date newAlarm) {
		setAnotherThread(aThread);
		setAlarm(newAlarm);
	}
	/**
	 * @return Returns the alarm.
	 */
	public Date getAlarm() {
		return alarm;
	}
	/**
	 * @param alarm The alarm to set.
	 */
	public void setAlarm(Date alarm) {
		this.alarm = alarm;
	}
	/**
	 * @return Returns the anotherThread.
	 */
	public Thread getAnotherThread() {
		return anotherThread;
	}
	/**
	 * @param anotherThread The anotherThread to set.
	 */
	public void setAnotherThread(Thread anotherThread) {
		this.anotherThread = anotherThread;
	}
	/**
	 * @return Returns the check.
	 */
	public long getCheck() {
		return check;
	}
	/**
	 * @param check The check to set.
	 */
	public void setCheck(long check) {
		this.check = check;
	}
}
