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
    c1 = coil.create_geom(False, lead_l, lead_r)
    x1, y1 = c1.get()
    xs.append(x1)
    ys.append(y1)
    
    c2 = c1.flipud().rotate(60, axis=0)
    xs.append(c2.x)
    ys.append(c2.y)
    
    c3 = c1.flipud().rotate(-60, axis=2)
    xs.append(c3.x)
    ys.append(c3.y)
        
    return xs, ys

def three_y(coil, lead_l, lead_r):

    xs, ys = [], []
    c = coil.create_geom(False, lead_l, lead_r)
    c1 = c.rotate(90, axis=0)
    xs.append(c1.x)
    ys.append(c1.y)
    
    c2 = c1.rotate(120, axis=0)
    xs.append(c2.x)
    ys.append(c2.y)
    
    c3 = c2.rotate(120, axis=0)
    xs.append(c3.x)
    ys.append(c3.y)        
    return xs, ys
    
def create_coil(c_t, p_dist, ncoil):
        
    if c_t == coil_t_list[0]:
        coil = cg.CircleCoil(p_dist=p_dist, ncoil=ncoil)
    elif c_t == coil_t_list[1]:
        coil = cg.EllipseCoil(p_dist=p_dist, ncoil=ncoil)
    elif c_t == coil_t_list[2]:
        coil = cg.EllipseCoilShape(p_dist=p_dist, ncoil=ncoil)
    elif c_t == coil_t_list[3]:
        coil = cg.EllipseCoilCurvature(p_dist=p_dist, ncoil=ncoil)
 
    return coil
    
def save_3phase_delta(fname, p_dist=0.4, ncoil=5, lead_l=4, lead_r=4, lcol='b', lthk=0.01, coil_t="ecc"):

    c_t = coil_t.lower()
    if c_t not in coil_t_list:
        return 
    
    coil = create_coil(c_t, p_dist, ncoil)
    xs, ys = three_delta(coil, lead_l, lead_r)
    
    dev = cxp.DevicePPT(fname, xs, ys)
    for x_, y_ in zip(xs, ys):
        dev.polyline(dev.xs_(x_), dev.ys_(y_), lcol, lthk)

    dev.close()        

def plot_3phase_delta(p_dist=0.4, ncoil=5, lead_l=5, lead_r=5, lcol='b', lthk=0.01, coil_t="ecc"):
    import matplotlib.pyplot as plt
    
    c_t = coil_t.lower()
    if c_t not in coil_t_list:
        return 
    
    coil = create_coil(c_t, p_dist, ncoil)
    xs, ys = three_delta(coil, lead_l, lead_r)
        
    plt.plot(xs[0], ys[0])
    plt.plot(xs[1], ys[1])
    plt.plot(xs[2], ys[2])
    plt.axis('equal')
    plt.show()
        
def save_3phase_y(fname, p_dist=0.4, ncoil=5, lead_l=6, lead_r=6, lcol='b', lthk=0.01, coil_t="ecc"):

    c_t = coil_t.lower()
    if c_t not in coil_t_list:
        return 
    
    coil = create_coil(c_t, p_dist, ncoil)
    xs, ys = three_y(coil, lead_l, lead_r)
    
    dev = cxp.DevicePPT(fname, xs, ys)
    for x_, y_ in zip(xs, ys):
        dev.polyline(dev.xs_(x_), dev.ys_(y_), lcol, lthk)

    dev.close()        

def plot_3phase_y(p_dist=0.4, ncoil=5, lead_l=5, lead_r=5, lcol='b', lthk=0.01, coil_t="ecc"):
    import matplotlib.pyplot as plt
    
    c_t = coil_t.lower()
    if c_t not in coil_t_list:
        return 
    
    coil = create_coil(c_t, p_dist, ncoil)
    xs, ys = three_y(coil, lead_l, lead_r)
        
    plt.plot(xs[0], ys[0])
    plt.plot(xs[1], ys[1])
    plt.plot(xs[2], ys[2])
    plt.axis('equal')
    plt.show()