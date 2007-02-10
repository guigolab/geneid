import java.io.File;

/*
 * Created on Jul 24, 2006
 *
 * To change the template for this generated file go to
 * Window&gt;Preferences&gt;Java&gt;Code Generation&gt;Code and Comments
 */

public class Renamer {
	public static void main(String[] args) {
		
		if (args.length< 2) {
			System.out.println("Usage: Renamer [path] [outPath] pattern(-a for append to all files) substitute\ne.g. Renamer /home/ug/root .tmp .del");
			System.exit(0);
		}
		
		int cnt= 0;
		int contain= 0;
		try {
			String fPath= "./";
			String pattern= args[0];
			String subst= args[1];
			if (args.length== 3) {
				fPath= args[0];
				pattern= args[1];
				subst= args[2];
			}
			String outPath= fPath;
			if (args.length== 4) {
				fPath= args[0];
				outPath= args[1];
				pattern= args[2];
				subst= args[3];
			}

			String[] dir= new File(fPath).list();
			contain= dir.length;
			for (int i = 0; dir!= null&& i < dir.length; i++) {
				StringBuffer sb= new StringBuffer(dir[i]);
				int x= -1;
				//while (i= sb.indexOf(args[1], i+1)>= 0)	// can have recursive impliction if (pattern contains substitute)
				if (pattern.equals("-a")) {
					sb.append(subst);
				} else {
					x= sb.indexOf(pattern, x+1);
					if (x>= 0)
						sb.replace(x, x+ pattern.length(), subst);
				}
				
				File f= new File(fPath+ File.separator+ dir[i]);
				String outName= outPath+ File.separator+ sb.toString();
				if (!f.getAbsolutePath().equals(outName)) {
					f.renameTo(new File(outName));
					//System.out.println(f.getAbsolutePath()+ " -> "+ outName);
					++cnt;
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.out.println(cnt+ " files renamed, "+ (contain- cnt)+ " files unchanged");
	}
}
