'''
    coil_3phase.py
    
    generate 3 phase delta & Y configuration

    04/25/2026
    
'''
import numpy as np

if __name__ == "__main__":
    import coil_geom as cg
    import coil_export
else:
    from . import coil_geom as cg
    from . import coil_export as cxp
    
# coil_t: 
#   cc : Circle Coil Circle
#   ec : Ellipse Coil Circle
#   ecc: Ellipse Coil Curvature
#   ecs: Ellipse Coil Shape

coil_t_list = ["cc", "ec", "ecs", "ecc"]
    
def three_delta(coil, lead_l, lead_r):

    xs, ys = [], []
    c = coil.create_geom(False, lead_l, lead_r)
    c_= c.copy()
    x1, y1 = c.copy_data()
    xs.append(x1)
    ys.append(y1)
    
    c.flipud().rotate(-60, axis=0)
    x2, y2 = c.copy_data()
    xs.append(x2)
    ys.append(y2)
    
    c_.flipud().rotate(60, axis=2)
    x3, y3 = c_.copy_data()
    xs.append(x3)
    ys.append(y3)
        
    return xs, ys
        
def save_3phase_delta(fname, p_dist=-0.4, ncoil=5, lead_l=4, lead_r=4, lcol='b', lthk=0.01, coil_t="ecc"):

    c_t = coil_t.lower()
    if c_t in coil_t_list:
        if c_t == coil_t_list[0]:
            coil = cg.CircleCoil(p_dist=p_dist, ncoil=ncoil)
        elif c_t == coil_t_list[1]:
            coil = cg.EllipseCoil(p_dist=p_dist, ncoil=ncoil)
        elif c_t == coil_t_list[2]:
            coil = cg.EllipseCoilShape(p_dist=p_dist, ncoil=ncoil)
        elif c_t == coil_t_list[3]:
            coil = cg.EllipseCoilCurvature(p_dist=p_dist, ncoil=ncoil)
   
        xs, ys = three_delta(coil, lead_l, lead_r)
        dev = cxp.DevicePPT(fname, xs, ys)
        for x_, y_ in zip(xs, ys):
            dev.polyline(dev.xs_(x_), dev.ys_(y_), lcol, lthk)

        dev.close()        
        
def three_y(cg):
    pass