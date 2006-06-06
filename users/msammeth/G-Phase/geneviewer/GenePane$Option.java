// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.*;

class GenePane$Option extends GenePane$Action {

    JFrame optionPane;

    public void actionPerformed(ActionEvent actionevent) {
        showOptionPane();
    }

    void init() {
        JPanel jpanel = new JPanel();
        jpanel.setLayout(new GridLayout(0, 2));
        jpanel.add(new JLabel("ORF Color"));
        JButton jbutton = new JButton("change");
        jbutton.setBackground(GenePane.ORFCOLOR1);
        jbutton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                Color color = JColorChooser.showDialog(null, "Choose ORF Color", GenePane.ORFCOLOR1);
                if(color == null) {
                    return;
                } else {
                    GenePane.ORFCOLOR1 = color;
                    optionPane.setVisible(false);
                    repaint();
                    return;
                }
            }

             {
                super();
            }
        });
        jpanel.add(jbutton);
        jpanel.add(new JLabel("UTR Color"));
        jbutton = new JButton("change");
        jbutton.setBackground(GenePane.UTRCOLOR1);
        jbutton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                Color color = JColorChooser.showDialog(null, "Choose UTR Color", GenePane.UTRCOLOR1);
                if(color == null) {
                    return;
                } else {
                    GenePane.UTRCOLOR1 = color;
                    optionPane.setVisible(false);
                    repaint();
                    return;
                }
            }

             {
                super();
            }
        });
        jpanel.add(jbutton);
        jpanel.add(new JLabel("Background Color"));
        jbutton = new JButton("change");
        jbutton.setBackground(GenePane.BGCOLOR);
        jbutton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                Color color = JColorChooser.showDialog(null, "Choose Background Color", GenePane.BGCOLOR);
                if(color == null) {
                    return;
                } else {
                    GenePane.BGCOLOR = color;
                    optionPane.setVisible(false);
                    GenePane.access$700(this$0).setBackground(GenePane.BGCOLOR);
                    GenePane.access$800(this$0).setBackground(GenePane.BGCOLOR);
                    GenePane.access$800(this$0).resetAllAlternativeContainersColor();
                    GenePane.access$900(this$0).setBackground(GenePane.BGCOLOR.darker());
                    repaint();
                    return;
                }
            }

             {
                super();
            }
        });
        jpanel.add(jbutton);
        jpanel.add(new JLabel("Range Color"));
        jbutton = new JButton("change");
        jbutton.setBackground(GenePane.ALTERNATIVERANGE_COLOR);
        jbutton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent actionevent) {
                Color color = JColorChooser.showDialog(null, "Choose Range Color", GenePane.ALTERNATIVERANGE_COLOR);
                if(color == null) {
                    return;
                } else {
                    GenePane.ALTERNATIVERANGE_COLOR = color;
                    optionPane.setVisible(false);
                    repaint();
                    return;
                }
            }

             {
                super();
            }
        });
        jpanel.add(jbutton);
        optionPane = new JFrame();
        optionPane.getContentPane().add(jpanel);
    }

    void showOptionPane() {
        optionPane.pack();
        optionPane.setVisible(true);
    }


    public GenePane$Option() {
        super(GenePane.this, "Option", GenePane.access$500());
        optionPane = null;
        mnemonic = 'o';
        toolTipText = "Show Options";
        init();
    }
}