package qalign.gui;

import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

//import org.eclipse.swt.widgets.Display;
//import org.eclipse.swt.widgets.ProgressBar;

/**
 * @author micha
 *
 * To change this generated comment edit the template variable "typecomment":
 * Window>Preferences>Java>Templates.
 * To enable and disable the creation of type comments go to
 * Window>Preferences>Java>Code Generation.
 */
public class GUIProxy {

        public Object display;

        /**
         * The GUI layer imposed:
         * 0 (not inited)
         * 1 for SWT
         * 2 for Swing
         */
        protected byte guiModeProgressBar= 0;
        /**
         * The GUI layer imposed:
         * 0 (not inited)
         * 1 for SWT
         * 2 for Swing
         */
        protected byte guiModeMessageLabel= 0;
        public static final byte GUI_NONE= 0;
        public static final byte GUI_SWT= 1;
        public static final byte GUI_SWING= 2;
        public static final String[] GUIS_SUPPORTED= {"<NONE>", "SWT", "Swing"};
        protected Object progressBar= null;
        protected Object messageLabel= null;

class ExecRunnable implements Runnable {

        Method method= null;
        Object target= null;
        Object[] args= null;
        Object result= null;

        public ExecRunnable(Method newMethod, Object newTarget, Object[] newArgs) {

                this.method= newMethod;
                this.target= newTarget;
                this.args= newArgs;
        }

        public void run() {

                try {
                        result= method.invoke(target, args);
                } catch (InvocationTargetException e) {
                        e.printStackTrace(); 	// nothing
                } catch (IllegalAccessException e) {
                        e.printStackTrace(); 	// nothing
                }
        }

        public Object getResult() {

                return result;
        }

}

        public void setProgressBar(Object newProgressBar) throws ClassCastException {

                progressBar= newProgressBar;

                try {
                        if (Class.forName("org.eclipse.swt.widgets.ProgressBar").isInstance(newProgressBar))
                                guiModeProgressBar= GUI_SWT;
                        else						// assuming...
                                guiModeProgressBar= GUI_SWING;
                } catch (ClassNotFoundException e) {
                                                                                // assuming...
                        guiModeProgressBar= GUI_SWING;
                }
        }

        public void setMessageLabel(Object newMessageLabel) {

                messageLabel= newMessageLabel;

                guiModeMessageLabel= GUI_NONE;
                try {
                        if (Class.forName("org.eclipse.swt.widgets.Label").isInstance(newMessageLabel))
                                guiModeMessageLabel= GUI_SWT;
                        else						// assuming...
                                guiModeMessageLabel= GUI_SWING;
                } catch (ClassNotFoundException e) {
                                                                                // assuming...
                        guiModeMessageLabel= GUI_SWING;
                }
        }

        public void setMessage(String newMessage) {

                        // error
                if ((messageLabel== null)|| (newMessage== null))
                        return;

                Class messageLabelClass= null;
                Class[] parameters= new Class[1];
                parameters[0]= String.class;
                try {
                        if (guiModeMessageLabel== GUI_SWT)
                                messageLabelClass= Class.forName("org.eclipse.swt.widgets.Label");
                        else
                                messageLabelClass= Class.forName("javax.swing.JLabel");
                } catch (ClassNotFoundException e) {
                        ; // cannot happen, checked in the set methods
                }

                Method method= null;
                try {
                        method= messageLabelClass.getMethod("setText", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }
                Object[] args= new Object[1];
                args[0]= newMessage;
                asyncGUIExec(new ExecRunnable(method, messageLabel, args), guiModeMessageLabel);
        }

        public void setMaximum(int newMaximum) {


                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return;

                Method method= null;
                Class[] parameters= new Class[1];
                parameters[0]= int.class;
                try {
                        method= tclass.getMethod("setMaximum", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }

                Object[] args= new Object[1];
                args[0]= new Integer(newMaximum);
                asyncGUIExec(new ExecRunnable(method, progressBar, args), guiModeProgressBar);
        }

        public void setMinimum(int newMinimum) {

                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return;

                Method method= null;
                Class[] parameters= new Class[1];
                parameters[0]= int.class;
                try {
                        method= tclass.getMethod("setMinimum", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }

                Object[] args= new Object[1];
                args[0]= new Integer(newMinimum);
                asyncGUIExec(new ExecRunnable(method, progressBar, args), guiModeProgressBar);
        }

        public void setValue(int newValue) {

                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return;

                Method method= null;
                Class[] parameters= new Class[1];
                parameters[0]= int.class;
                try {
                        if (guiModeProgressBar== GUI_SWT)
                                method= tclass.getMethod("setSelection", parameters);
                        else					// assuming...
                                method= tclass.getMethod("setValue", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }

                Object[] args= new Object[1];
                args[0]= new Integer(newValue);
                asyncGUIExec(new ExecRunnable(method, progressBar, args), guiModeProgressBar);
        }

        public Object asyncGUIExec(ExecRunnable run, byte modeUI) {

                if ((modeUI== GUI_SWT)&& (display!= null)) {

                        Class displayClass= null;
                        try {
                                displayClass= Class.forName("org.eclipse.swt.widgets.Display");
                                if (!displayClass.isInstance(display))
                                        return null;
                        } catch (ClassNotFoundException e) {
                                return null; // nothing
                        }

                        Class[] args= new Class[1];
                        args[0]= Runnable.class;		// why does 'run.getClass();' not work ?
                        Object[] parms= new Object[1];
                        parms[0]= run;
                        try {
                                Method asyncExec= displayClass.getMethod("asyncExec", args);
                                asyncExec.invoke(display, parms);
                                return run.getResult();
                        } catch (NoSuchMethodException e) {
                                ; // nothing
                        } catch (InvocationTargetException e) {
                                ; // nothing
                        } catch (IllegalAccessException e) {
                                ; // nothing
                        }
                }

                return null;
        }

        public Object syncGUIExec(ExecRunnable run, byte modeUI) {

                if ((modeUI== GUI_SWT)&& (display!= null)) {

                        Class displayClass= null;
                        try {
                                displayClass= Class.forName("org.eclipse.swt.widgets.Display");
                                if (!displayClass.isInstance(display))
                                        return null;
                        } catch (ClassNotFoundException e) {
                                return null; // nothing
                        }

                        Class[] args= new Class[1];
                        args[0]= Runnable.class;		// why does 'run.getClass();' not work ?
                        Object[] parms= new Object[1];
                        parms[0]= run;
                        try {
                                Method syncExec= displayClass.getMethod("syncExec", args);
                                syncExec.invoke(display, parms);
                                return run.getResult();
                        } catch (NoSuchMethodException e) {
                                ; // nothing
                        } catch (InvocationTargetException e) {
                                ; // nothing
                        } catch (IllegalAccessException e) {
                                ; // nothing
                        }

                }

                return null;

        }

        public void setDisplay(Object newDisplay) {

                display= newDisplay;
        }

        public int getValue() {

                int result= (-1);
                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return result;

                Method method= null;
                Class[] parameters= null;
                try {
                        if (guiModeProgressBar== GUI_SWT)
                                method= tclass.getMethod("getSelection", parameters);
                        else					// assuming...
                                method= tclass.getMethod("getValue", parameters);
                } catch (NoSuchMethodException e) {
                        e.printStackTrace(); // does not happen
                }


                Object res= syncGUIExec(new ExecRunnable(method, progressBar, null), guiModeProgressBar);
                if (res!= null)
                        result= ((Integer) res).intValue();

                return result;
        }


        public int getMaximum() {

                int result= (-1);
                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return result;

                Method method= null;
                Class[] parameters= null;
                try {
                        method= tclass.getMethod("getMaximum", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }

                Object res= syncGUIExec(new ExecRunnable(method, progressBar, null), guiModeProgressBar);
                if (res!= null)
                        result= ((Integer) res).intValue();

                return result;
        }


        public int getMinimum() {

                int result= (-1);
                Class tclass= getClassProgressBar();
                        // error
                if (tclass== null)
                        return result;

                Method method= null;
                Class[] parameters= null;
                try {
                        method= tclass.getMethod("getMinimum", parameters);
                } catch (NoSuchMethodException e) {
                        ; // does not happen
                }

                Object res= syncGUIExec(new ExecRunnable(method, progressBar, null), guiModeProgressBar);
                if (res!= null)
                        result= ((Integer) res).intValue();

                return result;
        }


        public Class getClassProgressBar() {
                                // error
                if (progressBar== null)
                        return null;

                Class tclass= null;
                try {
                        if (guiModeMessageLabel== GUI_SWT)
                                tclass= Class.forName("org.eclipse.swt.widgets.ProgressBar");
                        else
                                tclass= Class.forName("javax.swing.JProgressBar");
                } catch (ClassNotFoundException e) {
                        ; // cannot happen, checked in the set methods
                }

                return tclass;
        }

        public void setOutput(String output) {
                ;
        }

        public void increase() {

                int val= getValue();
                int max= getMaximum();

                if ((val+1)>= max)
                        setValue(max);
                else
                        setValue(val+1);

        }

}
