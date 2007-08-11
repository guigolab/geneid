package gphase.tools;

public class FilenameChecker {
	public static void main(String[] args) {
		checkFileNameLength();
	}
	
	public static void checkFileNameLength() {
		int x= 0;
		File dir= new File("delme");
		dir.mkdir();
		try {
			for (x = 1;; x++) {
				StringBuffer sb= new StringBuffer();
				for (int i = 0; i < x; i++) 
					sb.append('a');
				File f= new File(dir.getAbsolutePath()+File.separator+sb.toString());
				f.createNewFile();
				f.delete();
			}
		} catch (Exception e) {
			System.out.println(e.getMessage());
			System.out.println("max filename "+x);
		}
		dir.delete();
	}
}
