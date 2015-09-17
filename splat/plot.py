# Splat - splat/interpol.py
#
# Copyright (C) 2014 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from tkinter import *

def spline(spline, n=100, w=400, h=400):
    if (n < 2) or (w < 2) or (h < 2):
        raise ValueError("Invalid view parameters")

    step = (spline.end - spline.start) / float(n - 1)
    pts = list()
    ymin = ymax = None
    for i in range(n):
        x = spline.start + (i * step)
        y = spline.value(x)
        if ymin is None or y < ymin:
            ymin = y
        elif ymax is None or y > ymax:
            ymax = y
        pts.append(spline.value(x))

    master = Tk()
    c = Canvas(master, width=w, height=h)
    c.pack()
    c.create_rectangle(0, 0, w, h, fill="white")
    xratio = float(w) / float(n - 1)
    yratio = float(h) / (ymax - ymin)
    y0 = (pts[0] - ymin) * yratio
    for x, y in enumerate(pts[1:]):
        y = (y - ymin) * yratio
        x1 = (x + 1) * xratio
        x = x * xratio
        c.create_line(x, (h - y0), x1, (h - y))
        y0 = y
    dim = w / 100
    xratio = w / (spline.end - spline.start)
    for x0, y0 in spline.points():
        x = (x0 - spline.start) * xratio
        y = h - ((y0 - ymin) * yratio)
        c.create_oval((x - dim), (y - dim), (x + dim), (y + dim))
    mainloop()
