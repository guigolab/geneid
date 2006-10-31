package gphase.gui.pie;

// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 11/07/2006 08:44:21 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) 
// Source File Name:   Pie2D.java

import java.applet.Applet;
import java.awt.*;
import java.awt.geom.Arc2D;
import java.awt.geom.Point2D;
import java.text.DecimalFormat;
import java.util.StringTokenizer;

public class Pie2D extends Applet
{

    public Pie2D()
    {
    }

    public void init()
    {
        demo = 0;
        legend_border_off = 0;
        show_values_on_slices = 0;
        max = 0;
        distance = 30;
        max_description = 0;
        legend_rows = 3;
        show_legend = 1;
        init = 0;
        indicator_lines = 0;
        show_title = 1;
        show_border = 1;
        show_legend_on_right = 0;
        title_pos_y = 30;
        title_pos_x = -1;
        pie_x = -1D;
        pie_y = -1D;
        radius = -1D;
        legend_x = -1D;
        legend_y = -1D;
        legend_w = -1D;
        legend_h = -1D;
        setLayout(null);
        screen_w = size().width;
        screen_h = size().height;
        resize((int)screen_w, (int)screen_h);
        max_num_slices = 40;
        show_percents_on_slices = 0;
        show_percents_beside_slices = 0;
        show_percents_on_legend = 0;
        show_percents = 0;
        font1 = new Font("SansSerif", 0, 12);
        font2 = new Font("SansSerif", 1, 14);
        title = "Add here your title";
        numSlices = 0;
        total = 0.0D;
        value = new double[max_num_slices];
        color = new Color[max_num_slices];
        desc = new String[max_num_slices];
        startAngle = -45;
        pie_x = -1D;
        pie_y = -1D;
        radius = -1D;
        InitColors();
        getParam();
        if(demo == 1)
            InitDemoParams();
    }

    public void InitDemoParams()
    {
        demo_font = new Font("SansSerif", 1, 12);
        int d_red = backg_color.getRed();
        if(d_red + 30 < 255)
            d_red += 30;
        else
            d_red -= 30;
        int d_green = backg_color.getGreen();
        if(d_green + 30 < 255)
            d_green += 30;
        else
            d_green -= 30;
        int d_blue = backg_color.getBlue();
        if(d_blue + 30 < 255)
            d_blue += 30;
        else
            d_blue -= 30;
        demo_color = new Color(d_red, d_green, d_blue);
    }

    public void DrawDemoText()
    {
        g2b.setPaint(demo_color);
        g2b.setFont(demo_font);
        for(int y = 20; (double)y < screen_h; y += 60)
        {
            for(int x = 10; (double)x < screen_w; x += 40)
                g2b.drawString("DEMO", x, y);

        }

        for(int y = 50; (double)y < screen_h; y += 60)
        {
            for(int x = 30; (double)x < screen_w; x += 40)
                g2b.drawString("DEMO", x, y);

        }

    }

    public void InitColors()
    {
        color[0] = new Color(103, 105, 168);
        color[1] = new Color(99, 156, 255);
        color[2] = new Color(129, 126, 157);
        color[3] = new Color(140, 140, 140);
        color[4] = new Color(129, 135, 181);
        color[5] = new Color(52, 184, 222);
        color[6] = new Color(79, 237, 224);
        color[7] = new Color(202, 134, 177);
        color[8] = new Color(198, 99, 165);
        color[9] = new Color(99, 90, 255);
        color[10] = new Color(199, 174, 145);
        color[11] = new Color(219, 208, 165);
        color[12] = new Color(173, 198, 148);
        color[13] = new Color(137, 238, 151);
        color[14] = new Color(247, 189, 132);
        color[15] = new Color(206, 173, 156);
        color[16] = new Color(219, 208, 165);
        color[17] = new Color(99, 165, 156);
        color[18] = new Color(99, 204, 213);
        color[19] = new Color(244, 129, 114);
        color[20] = new Color(103, 105, 168);
        color[21] = new Color(99, 156, 255);
        color[22] = new Color(129, 126, 157);
        color[23] = new Color(140, 140, 140);
        color[24] = new Color(129, 135, 181);
        color[25] = new Color(52, 184, 222);
        color[26] = new Color(79, 237, 224);
        color[27] = new Color(202, 134, 177);
        color[28] = new Color(198, 99, 165);
        color[29] = new Color(99, 90, 255);
        color[30] = new Color(199, 174, 145);
        color[31] = new Color(219, 208, 165);
        color[32] = new Color(173, 198, 148);
        color[33] = new Color(137, 238, 151);
        color[34] = new Color(247, 189, 132);
        color[35] = new Color(206, 173, 156);
        color[36] = new Color(219, 208, 165);
        color[37] = new Color(99, 165, 156);
        color[38] = new Color(99, 204, 213);
        color[39] = new Color(244, 129, 114);
        backg_color = Color.white;
        title_color = Color.black;
        legend_color = Color.gray;
        border_color = Color.gray;
    }

    public void getParam()
    {
        String param = "";
        param = getParameter("legend_border_off");
        if(param != null)
            SetLegendBorderOff();
        param = getParameter("show_values_on_slices");
        if(param != null)
            SetValuesOnSlices();
        param = getParameter("set_legend_off");
        if(param != null)
            SetLegendOff();
        param = getParameter("set_title_off");
        if(param != null)
            ShowTitleOff();
        param = getParameter("title");
        if(param != null)
            SetTitle(param);
        param = getParameter("show_legend_on_right");
        if(param != null)
            ShowLegendOnRight();
        if(show_legend == 1 && show_legend_on_right == 0)
        {
            param = getParameter("legend_rows");
            if(param != null)
                SetLegendRows(Integer.parseInt(param));
            param = getParameter("legend_distance");
            if(param != null)
                SetDistance(Integer.parseInt(param));
        }
        param = getParameter("title_position_y");
        if(param != null)
        {
            param = param.trim();
            SetTitleYPos(Integer.parseInt(param));
        }
        param = getParameter("title_position_x");
        if(param != null)
        {
            param = param.trim();
            SetTitleXPos(Integer.parseInt(param));
        }
        param = getParameter("pie_x");
        if(param != null)
        {
            param = param.trim();
            SetPieX(Double.parseDouble(param));
        }
        param = getParameter("pie_y");
        if(param != null)
        {
            param = param.trim();
            SetPieY(Double.parseDouble(param));
        }
        param = getParameter("legend_x");
        if(param != null)
        {
            param = param.trim();
            SetLegendX(Double.parseDouble(param));
        }
        param = getParameter("legend_y");
        if(param != null)
        {
            param = param.trim();
            SetLegendY(Double.parseDouble(param));
        }
        param = getParameter("legend_w");
        if(param != null)
        {
            param = param.trim();
            SetLegendW(Double.parseDouble(param));
        }
        param = getParameter("legend_h");
        if(param != null)
        {
            param = param.trim();
            SetLegendH(Double.parseDouble(param));
        }
        param = getParameter("show_border_off");
        if(param != null)
            ShowBorderOff();
        param = getParameter("radius");
        if(param != null)
        {
            param = param.trim();
            SetRadius(Double.parseDouble(param));
        }
        param = getParameter("show_percents_on_slices");
        if(param != null)
            SetPercentsOnSlices();
        param = getParameter("show_percents_beside_slices");
        if(param != null)
            SetPercentsBesideSlices();
        if(show_percents_beside_slices == 1)
        {
            param = getParameter("show_indicator_lines");
            if(param != null)
                ShowIndicatorLines();
        }
        param = getParameter("show_percents_on_legend");
        if(param != null)
            SetPercentsOnLegend();
        param = getParameter("backg_color");
        Color c = null;
        if(param != null)
        {
            StringTokenizer t = new StringTokenizer(param, ",");
            int green;
            int blue;
            int red = green = blue = -1;
            String s = t.nextToken();
            s = s.trim();
            if(s != null)
                red = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                green = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                blue = Integer.parseInt(s);
            if(red != -1 && green != -1 && blue != -1)
            {
                c = new Color(red, green, blue);
                SetBackgColor(c);
            }
        }
        param = getParameter("title_color");
        c = null;
        if(param != null)
        {
            StringTokenizer t = new StringTokenizer(param, ",");
            int green;
            int blue;
            int red = green = blue = -1;
            String s = t.nextToken();
            s = s.trim();
            if(s != null)
                red = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                green = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                blue = Integer.parseInt(s);
            if(red != -1 && green != -1 && blue != -1)
            {
                c = new Color(red, green, blue);
                SetTitleColor(c);
            }
        }
        param = getParameter("legend_color");
        c = null;
        if(param != null)
        {
            StringTokenizer t = new StringTokenizer(param, ",");
            int green;
            int blue;
            int red = green = blue = -1;
            String s = t.nextToken();
            s = s.trim();
            if(s != null)
                red = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                green = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                blue = Integer.parseInt(s);
            if(red != -1 && green != -1 && blue != -1)
            {
                c = new Color(red, green, blue);
                SetLegendColor(c);
            }
        }
        param = getParameter("border_color");
        c = null;
        if(param != null)
        {
            StringTokenizer t = new StringTokenizer(param, ",");
            int green;
            int blue;
            int red = green = blue = -1;
            String s = t.nextToken();
            s = s.trim();
            if(s != null)
                red = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                green = Integer.parseInt(s);
            s = t.nextToken();
            s = s.trim();
            if(s != null)
                blue = Integer.parseInt(s);
            if(red != -1 && green != -1 && blue != -1)
            {
                c = new Color(red, green, blue);
                SetBorderColor(c);
            }
        }
        double val = 0.0D;
        String descr = "";
        for(int keep_reading = 1; keep_reading == 1;)
        {
            val = 0.0D;
            String nume_param = "val_" + (numSlices + 1);
            param = getParameter(nume_param);
            int has_val = 0;
            if(param != null)
            {
                param = param.trim();
                val = Double.parseDouble(param);
                has_val = 1;
            }
            nume_param = "description_" + (numSlices + 1);
            param = getParameter(nume_param);
            if(param != null)
            {
                param = param.trim();
                descr = param;
            }
            nume_param = "color_" + (numSlices + 1);
            param = getParameter(nume_param);
            c = null;
            if(param != null)
            {
                StringTokenizer t = new StringTokenizer(param, ",");
                int green;
                int blue;
                int red = green = blue = -1;
                String s = t.nextToken();
                s = s.trim();
                if(s != null)
                    red = Integer.parseInt(s);
                s = t.nextToken();
                s = s.trim();
                if(s != null)
                    green = Integer.parseInt(s);
                s = t.nextToken();
                s = s.trim();
                if(s != null)
                    blue = Integer.parseInt(s);
                if(red != -1 && green != -1 && blue != -1)
                    c = new Color(red, green, blue);
            }
            if(has_val == 1 && descr != null)
            {
                if(c == null)
                    AddPieSlice(val, descr);
                else
                    AddPieSlice(val, descr, c);
            } else
            {
                keep_reading = 0;
            }
        }

    }

    public void SetLegendBorderOff()
    {
        legend_border_off = 1;
    }

    public void SetValuesOnSlices()
    {
        show_values_on_slices = 1;
        show_percents_on_slices = 0;
    }

    public void SetDistance(int d)
    {
        distance = d;
    }

    public void SetLegendRows(int nr)
    {
        legend_rows = nr;
    }

    public void SetLegendOff()
    {
        show_legend = 0;
    }

    public void SetRadius(double r)
    {
        radius = r;
    }

    public void SetPieX(double x)
    {
        pie_x = x;
    }

    public void SetPieY(double y)
    {
        pie_y = y;
    }

    public void SetTitle(String t)
    {
        title = t;
    }

    public void SetTitleColor(Color c)
    {
        title_color = c;
    }

    public void ShowLegendOnRight()
    {
        show_legend_on_right = 1;
    }

    public void SetTitleYPos(int pos)
    {
        title_pos_y = pos;
    }

    public void SetTitleXPos(int pos)
    {
        title_pos_x = pos;
    }

    public void ShowTitleOff()
    {
        show_title = 0;
    }

    public void ShowBorderOff()
    {
        show_border = 0;
    }

    public void SetLegendX(double x)
    {
        legend_x = x;
    }

    public void SetLegendY(double y)
    {
        legend_y = y;
    }

    public void SetLegendW(double w)
    {
        legend_w = w;
    }

    public void SetLegendH(double h)
    {
        legend_h = h;
    }

    public void SetBackgColor(Color c)
    {
        backg_color = c;
    }

    public void SetLegendColor(Color c)
    {
        legend_color = c;
    }

    public void SetBorderColor(Color c)
    {
        border_color = c;
    }

    public void SetPercentsOnSlices()
    {
        show_percents_on_slices = 1;
        show_values_on_slices = 0;
    }

    public void SetPercentsBesideSlices()
    {
        show_percents_beside_slices = 1;
    }

    public void ShowIndicatorLines()
    {
        indicator_lines = 1;
    }

    public void SetPercentsOnLegend()
    {
        show_percents_on_legend = 1;
    }

    public void paint(Graphics g)
    {
        g2 = (Graphics2D)g;
        if(pie_image == null)
            pie_image = createImage((int)screen_w, (int)screen_h);
        backg = (Graphics2D)pie_image.getGraphics();
        g2b = (Graphics2D)backg;
        g2b.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2b.setPaint(backg_color);
        g2b.fillRect(0, 0, (int)screen_w, (int)screen_h);
        g2b.setColor(border_color);
        if(show_border == 1)
            g2b.drawRect(0, 0, (int)screen_w - 1, (int)screen_h - 1);
        if(init == 0)
        {
            CalcDimensions();
            init = 1;
        }
        if(demo == 1)
            DrawDemoText();
        paint_pie(backg);
        if(show_percents_on_slices == 1)
            DrawPercentsOnSlices(backg);
        if(show_values_on_slices == 1)
            DrawValuesOnSlices(backg);
        if(show_percents_beside_slices == 1)
            DrawPercentsBesideSlices(backg);
        if(show_legend == 1)
            paintLegend(backg);
        g2.drawImage(pie_image, 0, 0, null);
    }

    public void PaintTitle(Graphics backg)
    {
        g2b.setFont(font2);
        g2b.setColor(title_color);
        double start_at;
        if(title_pos_x == -1)
        {
            double l = g2b.getFontMetrics().stringWidth(title);
            start_at = screen_w / 2D - l / 2D;
        } else
        {
            start_at = title_pos_x;
        }
        g2b.drawString(title, (int)start_at, title_pos_y);
    }

    public void CalcDimensions()
    {
        g2b.setFont(font1);
        int max_len = 0;
        for(int i = 0; i < numSlices; i++)
            if(desc[i].length() > max_len)
            {
                max_len = desc[i].length();
                if(show_percents_on_legend == 1)
                    max_description = g2b.getFontMetrics().stringWidth(desc[i] + "(100.00%)");
                else
                    max_description = g2b.getFontMetrics().stringWidth(desc[i]);
            }

        if(legend_rows > numSlices)
            legend_rows = numSlices;
        g2b.setFont(font1);
        max = 0;
        for(int i = 0; i < numSlices; i++)
            if((int)value[i] > max)
                max = (int)value[i];

        if(show_legend == 0)
        {
            if(radius == -1D)
                radius = (int)((4D * screen_w) / 6D);
            if(radius + pie_y > screen_h - 20D)
                radius = (6D * (screen_h - pie_y)) / 8D;
            if(pie_x == -1D)
                pie_x = (screen_w - radius) / 2D;
            if(pie_y == -1D)
                pie_y = (double)title_pos_y + (screen_h - radius - (double)title_pos_y) / 2D;
        } else
        if(show_legend_on_right == 0)
        {
            if(legend_x == -1D)
                legend_x = 30D;
            if(legend_h == -1D)
                legend_h = legend_rows * 20;
            if(legend_w == -1D)
                legend_w = screen_w - 2D * legend_x;
            if(legend_y == -1D)
                legend_y = screen_h - legend_h - 10D;
            if(radius == -1D)
                radius = (int)((4D * screen_w) / 6D);
            if(radius + 30D + legend_h > screen_h - 20D)
                radius = (6D * (screen_h - (double)title_pos_y - legend_h)) / 8D;
            if(pie_x == -1D)
                pie_x = (screen_w - radius) / 2D;
            if(pie_y == -1D)
                pie_y = (double)title_pos_y + (screen_h - legend_h - radius - (double)title_pos_y) / 2D;
        } else
        {
            if(legend_h == -1D)
                legend_h = numSlices * 20;
            if(legend_w == -1D)
                legend_w = max_description + latura + 10;
            if(legend_y == -1D)
                legend_y = (screen_h - legend_h) / 2D;
            if(legend_x == -1D)
                legend_x = screen_w - legend_w - 20D;
            if(radius == -1D)
                radius = (int)((4D * (screen_w - legend_w - (double)title_pos_y)) / 6D);
            if(radius + (double)title_pos_y > screen_h - 20D)
                radius = (6D * (screen_h - (double)title_pos_y)) / 8D;
            if(pie_x == -1D)
                pie_x = (screen_w - legend_w - 20D - radius) / 2D;
            if(pie_y == -1D)
                pie_y = (double)title_pos_y + (screen_h - radius - (double)title_pos_y) / 2D;
        }
    }

    public void paintLegend(Graphics backg)
    {
        int vert_dist = 20;
        int a = 0;
        int b = 0;
        g2b.setFont(font1);
        if(show_legend_on_right == 0)
        {
            double dist = max_description + latura + distance + vert_dist;
            for(int i = 0; i < numSlices; i++)
                if(desc[i] != null)
                {
                    g2b.setColor(color[i]);
                    g2b.fillRect((int)(legend_x + (double)b * dist), (int)(legend_y + (double)(a * vert_dist)), latura, latura);
                    g2b.setColor(Color.gray);
                    g2b.draw3DRect((int)(legend_x + (double)b * dist), (int)(legend_y + (double)(a * vert_dist)), latura, latura, true);
                    g2b.setColor(legend_color);
                    String str = desc[i];
                    if(show_percents_on_legend == 1)
                    {
                        double percent = (value[a] / total) * 100D;
                        DecimalFormat fmt = new DecimalFormat("0.##");
                        String s = " " + fmt.format(percent) + "%";
                        str = str + '(' + s + ')';
                    }
                    g2b.drawString(str, (int)(legend_x + (double)b * dist + 20D), (int)(legend_y + (double)(a * vert_dist) + 10D));
                    if(i % legend_rows == legend_rows - 1)
                    {
                        b++;
                        a = 0;
                    } else
                    {
                        a++;
                    }
                }

        } else
        {
            for(int i = 0; i < numSlices; i++)
                if(desc[i] != null)
                {
                    g2b.setColor(color[i]);
                    g2b.fillRect((int)legend_x, (int)(legend_y + (double)(a * vert_dist)), latura, latura);
                    g2b.setColor(Color.gray);
                    g2b.draw3DRect((int)legend_x, (int)(legend_y + (double)(a * vert_dist)), latura, latura, true);
                    String str = desc[i];
                    if(show_percents_on_legend == 1)
                    {
                        double percent = (value[a] / total) * 100D;
                        DecimalFormat fmt = new DecimalFormat("0.##");
                        String s = " " + fmt.format(percent) + "%";
                        str = str + '(' + s + ')';
                    }
                    g2b.setColor(legend_color);
                    g2b.drawString(str, (int)(legend_x + 20D), (int)(legend_y + (double)(a * vert_dist) + 10D));
                    a++;
                }

            if(legend_border_off == 0)
                g2b.drawRect((int)legend_x - 5, (int)legend_y - 5, (int)legend_w + 10, (int)legend_h + 5);
        }
    }

    public void DrawPercentsOnSlices(Graphics backg)
    {
        g2b.setFont(font1);
        for(int i = 0; i < numSlices; i++)
        {
            g2b.setPaint(Color.black);
            angle = (int)Math.round(360D * (value[i] / total));
            Arc2D arc = new java.awt.geom.Arc2D.Double(pie_x, pie_y, radius, radius, startAngle, angle / 2, 2);
            Point2D p = arc.getEndPoint();
            Point2D c = new java.awt.geom.Point2D.Double();
            c.setLocation(arc.getCenterX(), arc.getCenterY());
            int x = (int)((p.getX() + (c.getX() - p.getX()) / 3D) - 10D);
            int y = (int)(p.getY() + (c.getY() - p.getY()) / 3D) + 7;
            double percent = (value[i] / total) * 100D;
            DecimalFormat fmt = new DecimalFormat("0.##");
            String str = fmt.format(percent) + "%";
            g2b.drawString(str, x, y);
            startAngle += angle;
        }

    }

    public void DrawValuesOnSlices(Graphics backg)
    {
        g2b.setFont(font1);
        for(int i = 0; i < numSlices; i++)
        {
            g2b.setPaint(Color.black);
            angle = (int)Math.round(360D * (value[i] / total));
            Arc2D arc = new java.awt.geom.Arc2D.Double(pie_x, pie_y, radius, radius, startAngle, angle / 2, 2);
            Point2D p = arc.getEndPoint();
            Point2D c = new java.awt.geom.Point2D.Double();
            c.setLocation(arc.getCenterX(), arc.getCenterY());
            int x = (int)((p.getX() + (c.getX() - p.getX()) / 3D) - 10D);
            int y = (int)(p.getY() + (c.getY() - p.getY()) / 3D) + 7;
            String str = Integer.toString((int) value[i]);
            g2b.drawString(str, x, y);
            startAngle += angle;
        }

    }

    public void DrawPercentsBesideSlices(Graphics backg)
    {
        g2b.setFont(font1);
        for(int i = 0; i < numSlices; i++)
        {
            g2b.setPaint(Color.black);
            angle = (int)Math.round(360D * (value[i] / total));
            Arc2D arc = new java.awt.geom.Arc2D.Double(pie_x, pie_y, radius, radius, startAngle, angle / 2, 2);
            Point2D p = arc.getEndPoint();
            Point2D c = new java.awt.geom.Point2D.Double();
            c.setLocation(arc.getCenterX(), arc.getCenterY());
            double percent = (value[i] / total) * 100D;
            DecimalFormat fmt = new DecimalFormat("0.##");
            String str = fmt.format(percent) + "%";
            double dist1 = 20D;
            double dist2 = 5D;
            double px = p.getX();
            double py = p.getY();
            double cx = c.getX();
            double cy = c.getY();
            double x1;
            double x_str;
            if(px > cx)
            {
                x1 = px + dist1;
                x_str = x1;
            } else
            {
                x1 = px - dist1;
                x_str = x1 - 40D;
            }
            double y1;
            double y_str;
            if(py > cy)
            {
                y1 = py + dist1;
                y_str = y1 + 10D;
            } else
            {
                y1 = py - dist2;
                y_str = y1;
            }
            if(indicator_lines == 1)
                if(py > cy)
                    g2b.drawLine((int)x1, (int)y1, (int)px, (int)py);
                else
                    g2b.drawLine((int)x1, (int)y1, (int)px, (int)py);
            g2b.drawString(str, (int)x_str, (int)y_str);
            startAngle += angle;
        }

    }

    public void paint_pie(Graphics backg)
    {
        g2b = (Graphics2D)backg;
        startAngle = -45;
        if(show_title == 1)
            PaintTitle(backg);
        for(int i = 0; i < numSlices; i++)
        {
            g2b.setPaint(color[i]);
            angle = (int)Math.round(360D * (value[i] / total));
            g2b.fill(new java.awt.geom.Arc2D.Double(pie_x, pie_y, radius, radius, startAngle, angle, 2));
            startAngle += angle;
        }

        g2b.fill(new java.awt.geom.Arc2D.Double(pie_x, pie_y, radius, radius, startAngle, 315 - startAngle, 2));
    }

    public void AddPieSlice(double value, String description)
    {
        this.value[numSlices] = value;
        desc[numSlices] = description;
        numSlices++;
        total += value;
    }

    public void AddPieSlice(double value, String description, Color c)
    {
        this.value[numSlices] = value;
        desc[numSlices] = description;
        color[numSlices] = c;
        numSlices++;
        total += value;
    }

    public void Destroy()
    {
        backg.dispose();
        g2.dispose();
        g2b.dispose();
        pie_image.flush();
        System.gc();
    }

    int demo;
    Font demo_font;
    Color demo_color;
    int legend_border_off;
    int show_values_on_slices;
    int max;
    int distance;
    int max_description;
    int legend_rows;
    int show_legend;
    int init;
    static int latura = 12;
    String title;
    int title_pos_y;
    int title_pos_x;
    int show_title;
    int show_border;
    int show_legend_on_right;
    int indicator_lines;
    double legend_x;
    double legend_y;
    double legend_w;
    double legend_h;
    int numSlices;
    double total;
    double value[];
    double screen_h;
    double screen_w;
    int startAngle;
    int angle;
    double pie_x;
    double pie_y;
    double radius;
    static int max_num_slices;
    int show_percents_on_slices;
    int show_percents_beside_slices;
    int show_percents_on_legend;
    int show_percents;
    Color color[];
    Color backg_color;
    Color title_color;
    Color legend_color;
    Color border_color;
    String desc[];
    Image pie_image;
    Graphics backg;
    Graphics2D g2;
    Graphics2D g2b;
    Font font1;
    Font font2;

}