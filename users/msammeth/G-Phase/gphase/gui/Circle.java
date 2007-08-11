package gphase.gui;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import javax.swing.JPanel;

import javax.swing.JPanel;;

public class Circle extends JPanel {
	int diameter= 0;
	Color color= null;
	public Circle(Color c, int dia) {
		super();
		this.diameter= dia;
		this.color= c;
		setOpaque(false);
	}

	@Override
	public Dimension getPreferredSize() {
		return new Dimension(diameter, diameter);
	}
	
	@Override
	protected void paintComponent(Graphics g) {
		super.paintComponent(g);
		g.setColor(color);
		g.fillOval(0, 0, diameter, diameter);
		g.setColor(Color.black);
		g.drawOval(0, 0, diameter, diameter);
	}

	public Color getColor() {
		return color;
	}
}
