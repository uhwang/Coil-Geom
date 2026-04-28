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
phase_t_list = ["delta", "y"]
phase_t_delta = phase_t_list[0]
phase_t_y = phase_t_list[1]
    
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

def create_3phase(coil_t, p_dist, ncoil, lead_l, lead_r, phase_t):
    c_t = coil_t.lower()
    if c_t not in coil_t_list:
        raise ValueError(f"Error(create_3phase): invalid coil type->{c_t}")
    
    coil = create_coil(c_t, p_dist, ncoil)
    if phase_t == phase_t_delta: 
        xs, ys = three_delta(coil, lead_l, lead_r)
    elif phase_t == phase_t_y: 
        xs, ys = three_y(coil, lead_l, lead_r)
    else:
        raise ValueError(f"Error(create_3phase): invalid phase type->{phase_t}")
        
    return xs, ys
    
def save_3phase_delta(fname, p_dist=0.4, ncoil=5, lead_l=4, lead_r=4, lcol='b', lthk=0.01, coil_t="ecc"):
    try:
        xs, ys = create_3phase(coil_t, p_dist, ncoil, lead_l, lead_r, phase_t_delta)
    except Exception as e:
        print(f"Error(save_3phase_delta->create_3phase): {e}")
        return
        
    dev = cxp.DevicePPT(fname, xs, ys)
    for x_, y_ in zip(xs, ys):
        dev.polyline(dev.xs_(x_), dev.ys_(y_), lcol, lthk)

    dev.close()  

    return xs, ys    

def plot_3phase_delta(p_dist=0.4, ncoil=5, lead_l=5, lead_r=5, lcol='b', lthk=0.01, coil_t="ecc"):
    import matplotlib.pyplot as plt
    
    try:
        xs, ys = create_3phase(coil_t, p_dist, ncoil, lead_l, lead_r, phase_t_y)
    except Exception as e:
        print(f"Error(plot_3phase_delta->create_3phase): {e}")
        return
        
    plt.plot(xs[0], ys[0], color=lcol)
    plt.plot(xs[1], ys[1], color=lcol)
    plt.plot(xs[2], ys[2], color=lcol)
    plt.axis('equal')
    plt.show()
        
def save_3phase_y(fname, p_dist=0.4, ncoil=5, lead_l=6, lead_r=6, lcol='b', lthk=0.01, coil_t="ecc"):

    try:
        xs, ys = create_3phase(coil_t, p_dist, ncoil, lead_l, lead_r, phase_t_y)
    except Exception as e:
        print(f"Error(save_3phase_y->create_3phase): {e}")
        return
    
    dev = cxp.DevicePPT(fname, xs, ys)
    for x_, y_ in zip(xs, ys):
        dev.polyline(dev.xs_(x_), dev.ys_(y_), lcol, lthk)

    dev.close()

    return xs, ys

def plot_3phase_y(p_dist=0.4, ncoil=5, lead_l=5, lead_r=5, lcol='b', lthk=0.01, coil_t="ecc"):
    import matplotlib.pyplot as plt
    
    try:
        xs, ys = create_3phase(coil_t, p_dist, ncoil, lead_l, lead_r, phase_t_y)
    except Exception as e:
        print(f"Error(plot_3phase_y->create_3phase): {e}")
        return
        
    plt.plot(xs[0], ys[0], color=lcol)
    plt.plot(xs[1], ys[1], color=lcol)
    plt.plot(xs[2], ys[2], color=lcol)
    plt.axis('equal')
    plt.show()