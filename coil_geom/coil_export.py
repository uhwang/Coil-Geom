'''
    04/08/2026  Export coil data to pptx
    04/09/2026  Export coil data to SVG
    04/11/2026  Debug Plot
'''
from abc import ABC, abstractmethod
    
import collections 
import collections.abc
import pptx
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
import numpy as np

if __name__ == "__main__":
    import coil_color
    import coil_util as cu
else:
    from . import coil_color
    from . import coil_util as cu

def get_color(c):
    cc = coil_color.color_collection
    
    if c in cc:
        c_ = cc[c]
    else:
        c_ = cc['k']
        
    return (c_[0], c_[1], c_[2])
    
def Polyline(slide, x, y, lcol, lthk, fcol, closed):
    lcol = get_color(lcol)
    fcol = get_color(fcol)

    free_form = slide.shapes.build_freeform(Inches(x[0]), Inches(y[0]))
    free_form.add_line_segments(
        [(Inches(x1), Inches(y1)) for x1, y1 in zip(x[1:],y[1:])],
        close=closed
    )
    
    line_shape = free_form.convert_to_shape()
    
    if closed == False: 
        line_shape.fill.background()
    else:
        if isinstance(fcol, tuple):
            line_shape.fill.solid()
            line_shape.fill.fore_color.rgb = RGBColor(fcol[0], fcol[1], fcol[2])
            # lcol is None, set lcol as fcol if not the default colot works
            if not isinstance(lcol, tuple):
                line_shape.line.color.rgb = RGBColor(fcol[0], fcol[1], fcol[2])
                line_shape.line.width = Inches(0.001) # dummy thinkness
        else:
            line_shape.fill.background()    
    
    if lcol:
        if fcol and lcol == fcol:
            line_shape.line.color.rgb = RGBColor(fcol[0], fcol[1], fcol[2])
            line_shape.line.width = Inches(lthk)
        else:
            line_shape.line.color.rgb = RGBColor(lcol[0], lcol[1], lcol[2])
            line_shape.line.width = Inches(lthk)
        
    line_shape.shadow.inherit     = False
    line_shape.shadow.blur_radius = 0
    line_shape.shadow.distance    = 0
 
    
def draw_transition_ellipse(dev, co):
    ex, ey, ts = co.c1x, co.c1y, co.t_star
    aa, bb = co.axlen, co.bxlen
    xt = aa*np.cos(ts)
    yt = bb*np.sin(ts)
    dd = np.linspace(0, 2*np.pi, 100)
    
    px, py = co.P[2], co.P[3]
    xx1 = co.c1x+aa*np.cos(dd)
    yy1 = co.c1y+bb*np.sin(dd)
    
    xx2 = xx1+co.r_dist
    yy2 = yy1
    
    xs_c = px + co.r_fillet*np.cos(dd)
    ys_c = py + co.r_fillet*np.sin(dd)
    dev.polyline(xs_c, ys_c, lcol="gold", lthk=0.01)
    
    if co.coil_type == cu._coil_type_ellipse_curvature or \
       co.coil_type == cu._coil_type_ellipse_shape:
        xx3 = co.xc_s+co.a_s*np.cos(dd)
        yy3 = co.yc_s+co.b_s*np.sin(dd)
        dev.polyline(xx3, yy3, lcol="blueviolet", lthk=0.01)

    p1x, p1y = ex+xt, ey+yt
    p2x, p2y = px-(ex+xt-px), ey+yt
    dev.polyline([px, p1x], [py, p1y], lcol="cadetblue", lthk=0.01)
    dev.polyline([px, p2x], [py, p2y], lcol="cadetblue", lthk=0.01)
        
class Device():
    def __init__(self, xx, yy):
        self.xmin, self.xmax = min(xx), max(xx)
        self.ymin, self.ymax = min(yy), max(yy)
        self.max_range = max(self.xmax-self.xmin, self.ymax-self.ymin)
        
    @abstractmethod
    def xs_(self, xs): pass
    @abstractmethod
    def x_(self, x): pass
    
    @abstractmethod
    def ys_(self, ys): pass
    @abstractmethod
    def y_(self, y): pass

    @abstractmethod
    def close(self): pass
    
    @abstractmethod
    def polyline(self, xs, ys, lcol, lthk): pass
    
class DevicePPT(Device):
    def __init__(self, fname, xx, yy):
        super().__init__(xx,yy)
        self.fname = fname
        self.ppt = pptx.Presentation()
        self.blank_slide_layout = self.ppt.slide_layouts[6]
        self.slide = self.ppt.slides.add_slide(self.blank_slide_layout)
        
        # Get width and height in EMUs
        width_emu = self.ppt.slide_width
        height_emu = self.ppt.slide_height
        
        # Convert to inches for readability
        width_inches = width_emu / 914400
        height_inches = height_emu / 914400
        min_edge = min(width_inches, height_inches)
        self.scale = min_edge/self.max_range
        
    def xs_(self, xs): return (xs-self.xmin)*self.scale
    def ys_(self, ys): return (ys-self.ymin)*self.scale
    
    def x_(self, x): return (x-self.xmin)*self.scale
    def y_(self, y): return (y-self.ymin)*self.scale
    
    def polyline(self, xs, ys, lcol, lthk):
        Polyline(self.slide, self.xs_(xs), self.ys_(ys), lcol, lthk, 'w', False)
    
    def close(self):
        self.ppt.save(self.fname)
     

def save_ppt(coil, fname, trace=False, lcol='b', lthk=0.005, debug=False):
    fc = 'w'
    if trace:
        c1, c2, c3 = 'r', 'g', 'b'
        x, y, seg = coil.create_geom(trace)
        dev = DevicePPT(fname, x, y)
        
        bit=True
        for i,s in enumerate(seg):
            ii=i+1
            if ii%2 == 0:
                Polyline(dev.slide, dev.xs_(s[0]), dev.ys_(s[1]), c2, lthk, fc, False)
            else:
                Polyline(dev.slide, dev.xs_(s[0]), dev.ys_(s[1]), c1 if bit else c3, lthk, fc, False)
                bit = not bit
    else:
        x, y = coil.create_geom(trace)
        dev = DevicePPT(fname, x, y)
        Polyline(dev.slide, dev.xs_(x), dev.ys_(y), lcol, lthk, fc, False)
        
    if debug:
        draw_transition_ellipse(dev, coil)
        
    dev.close()
    
'''
    Device SVG
'''
    
_line_format_begin = "<line x1=\"%3.3f\" y1=\"%3.3f\" x2=\"%3.3f\" y2=\"%3.3f\" "
_line_format_end = " style=\"fill:none;stroke:rgb(%d,%d,%d);stroke-width:%3.3f\" />\n"
_polygon_format_end = "style=\"stroke:rgb(%d,%d,%d);stroke-width:%d;fill:rgb(%d,%d,%d);\"/>\n"
_polygon_format_end_nostroke = "style=\"stroke:none;fill:rgb(%d,%d,%d);\"/>\n"
_polygon_format_end_nofill = "style=\"stroke:rgb(%d,%d,%d);stroke-width:%d;fill:none;\"/>\n"
_next_line = "\n"              
_points_per_line = 5
_move_to_cmd = 'M'
_line_to_cmd = "L"
               
class DeviceSVG(Device):
    def __init__(self, fname, xx, yy, hgt=400, ar=1.33, xgap=3, ygap=3):
        super().__init__(xx, yy)
        wid = int(hgt*ar)
        self.fp = open(fname, "wt")
        self.fp.write("<svg version=\"1.1\"\n"\
                      "width=\"%d\" height=\"%d\"\n"
                      "xmlns=\"http://www.w3.org/2000/svg\">\n"\
                      %(int(wid),int(hgt)))
        self.scale = min(hgt,wid)/self.max_range
        self.xx = xx
        self.yy = yy
        self.xgap = xgap
        self.ygap = ygap
        
    def polyline(self, xx, yy, lcol, lthk):
        lc = get_color(lcol)
        lt = self.scale*lthk
        self.fp.write("<polyline points=\"\n")
        self.create_pnt_list(xx, yy)
        self.fp.write(_line_format_end%(\
                      lc[0], 
                      lc[1], 
                      lc[2],
                      1 if lt < 0.001 else lt))
                      
    def create_pnt_list(self, xx, yy):
        xs = self.xs_(xx)
        ys = self.ys_(yy)
        
        for i, (x1, y1) in enumerate(zip(xs, ys)):
            self.fp.write("%3.3f %3.3f, "%(x1, y1))
            if (i+1)%_points_per_line == 0:
                self.fp.write(_next_line)
        self.fp.write("\"\n")      
        
    def xs_(self, xs): return self.xgap+(xs-self.xmin)*self.scale
    def ys_(self, ys): return self.ygap+(ys-self.ymin)*self.scale

    def x_(self, x): return self.xgap+(x-self.xmin)*self.scale
    def y_(self, y): return self.ygap+(y-self.ymin)*self.scale
    
    def close(self):
        self.fp.write("</svg>")
        self.fp.close()     

def save_svg(coil, fname, trace=False, lcol='b', lthk=0.005, debug=False):
    
    if trace:
        c1, c2, c3 = 'r', 'g', 'b'
        x, y, seg = coil.create_geom(trace)
        dev = DeviceSVG(fname, x, y)
        
        bit=True
        for i,s in enumerate(seg):
            ii=i+1
            if ii%2 == 0:
                dev.polyline(s[0], s[1], c2, lthk)
            else:
                dev.polyline(s[0], s[1], c1 if bit else c3, lthk)
                bit = not bit
    else:
        x, y = coil.create_geom(trace)
        dev = DeviceSVG(fname, x, y)
        dev.polyline(x, y, lcol, lthk)
        
    if debug:
        draw_transition_ellipse(dev, coil)
        
    dev.close()        