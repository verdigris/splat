from Tkinter import *

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
    mainloop()
