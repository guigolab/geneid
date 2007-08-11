package gphase.model;

public class DirectedSite extends AbstractSite {

	public static class PositionComparator extends AbstractSite.PositionComparator {
	
		public int compare(Object arg0, Object arg1) {
			
			DirectedSite s1= null, s2= null;
			try {
				s1= (DirectedSite) arg0;
				s2= (DirectedSite) arg1;
			} catch (ClassCastException e) {
				return super.compare(arg0, arg1);
			}
			
			
			int pos1= s1.isForward()?s1.getPos():-s1.getPos();
			int pos2= s2.isForward()?s2.getPos():-s2.getPos();
			
			if (pos1< pos2)
				return -1;
			if (pos1> pos2)
				return 1;
			return 0;
		}
		
	}

	public boolean isForward() {
		if (pos< 0)
			return false;
		return true;
	}
}
