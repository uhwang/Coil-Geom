'''
    04/08/2026  Export coil data to pptx
'''
if __name__ == "__main__":
    import coil_color
else:
    from . import coil_color

import collections 
import collections.abc
import pptx
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
import numpy as np

def Polyline(slide, x, y, lcol, lthk, fcol, closed):

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
    
class DevicePPT():
    def __init__(self, fname, xx, yy):
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
        
        xmin, xmax = min(xx), max(xx)
        ymin, ymax = min(yy), max(yy)
        
        max_range = max(xmax-xmin, ymax-ymin)
        min_edge = min(width_inches, height_inches)
        
        self.xmin = xmin
        self.ymin = ymin
        self.scale = min_edge/max_range
        
    def x_(self, xs): return (xs-self.xmin)*self.scale
    def y_(self, ys): return (ys-self.ymin)*self.scale
    
    def close(self):
        self.ppt.save(self.fname)
     
def get_color(c):
    cc = coil_color.color_collection
    
    if c in cc:
        c_ = cc[c]
    else:
        c_ = cc['k']
        
    return (c_[0], c_[1], c_[2])
    
def save_ppt(coil, fname, trace=False, lcol='b', lthk=0.005, fcol='w', closed=False):
    fc = get_color(fcol)
    
    if trace:
        c1 = get_color('r')
        c2 = get_color('g')
        c3 = get_color('b')
        x, y, seg = coil.create_geom(trace)
        dev = DevicePPT(fname, x, y)
        
        bit=True
        for i,s in enumerate(seg):
            ii=i+1
            if ii%2 == 0:
                Polyline(dev.slide, dev.x_(s[0]), dev.y_(s[1]), c2, lthk, fc, closed)
            else:
                Polyline(dev.slide, dev.x_(s[0]), dev.y_(s[1]), c1 if bit else c3, lthk, fc, closed)
                bit = not bit
    else:
        lc = get_color(lcol)
        x, y = coil.create_geom(trace)
        
        xmin, xmax = min(x), max(x)
        ymin, ymax = min(y), max(y)
        dev = DevicePPT(fname, x, y)
        
        Polyline(dev.slide, dev.x_(x), dev.y_(y), lc, lthk, fc, closed)
    dev.close()