'''
    coil_print.py
    
    03/06/2026
'''
if __name__ == "__main__":
    import coil_util as cu
else:
    from . import coil_util as cu

_deg = 180/3.14159265358979

def print_coil(c):
    direction = "DN" if c.p_dist >= 0 else "UP"
    
    info_1 = f"====================================\n" \
             f"Type : {c.coil_type}\n" \
             f"Dir  : {direction}\n" \
             f"====================================\n" \
             f"c1x  : {c.c1x}\n" \
             f"c1y  : {c.c1y}\n" \
             f"c2x  : {c.c2x}\n" \
             f"c2y  : {c.c2y}\n" \
             f"_v1  : {', '.join(f'{x:.3f}' for x in c._v1)}\n" \
             f"_v2  : {', '.join(f'{x:.3f}' for x in c._v2)}\n" \
             f"_v3  : {', '.join(f'{x:.3f}' for x in c._v3)}\n" \
             f"_v4  : {', '.join(f'{x:.3f}' for x in c._v4)}\n" \
             f"_v5  : {', '.join(f'{x:.3f}' for x in c._v5)}\n" \
             f" v1  : {', '.join(f'{x:.3f}' for x in c.v1)}\n" \
             f" v2  : {', '.join(f'{x:.3f}' for x in c._v2)}\n" \
             f" v3  : {', '.join(f'{x:.3f}' for x in c._v3)}\n" \
             f" v4  : {', '.join(f'{x:.3f}' for x in c._v4)}\n" \
             f" v5  : {', '.join(f'{x:.3f}' for x in c._v5)}\n" \
             f"_P   : {', '.join(f'{x:.3f}' for x in c._P)}\n" \
             f" P   : {', '.join(f'{x:.3f}' for x in c.P)}\n" \
             f"P1   : {', '.join(f'{x:.3f}' for x in c.P1)}\n" \
             f"P2   : {', '.join(f'{x:.3f}' for x in c.P2)}\n" \
             f"p_d  : {c.p_dist:2.3f}\n" \
             f"r_f  : {c.r_fillet:2.3f}\n" \
             f"t_*  : {c.t_star:2.3f}({c.t_star*_deg:2.3f})\n" \
             f"Npnt : {c.npnt}\n" \
             f"Ax   : {c.axlen}\n" \
             f"Bx   : {c.bxlen}\n" \
             f"====================================\n"
             
    if hasattr(c, "a_s"):
        info_2 = f"Ys Optimize\n" \
                 f"====================================\n" \
                 f"a_s  : {c.a_s:2.3f}\n" \
                 f"b_s  : {c.b_s:2.3f}\n" \
                 f"y_s  : {c.y_s:2.3f}\n" \
                 f"xc_s : {c.xc_s:2.3f}\n" \
                 f"yc_s : {c.yc_s:2.3f}\n" \
                 f"phi_s: {c.phi_s:2.3f}({c.phi_s*_deg:2.3f})\n" \
                 f"phi_e: {c.phi_e:2.3f}({c.phi_e*_deg:2.3f})\n" \
                 f"Delta: {c.delta:2.3f}\n" 
                 
        if c.coil_type == cu._coil_type_ellipse_shape:
            info_3 = f"b/a  : {c.bxlen/c.axlen:2.3f}\n" \
                     f"bs/as: {c.b_s/c.a_s:2.3f}\n" \
                     "====================================\n" 
        else:
            info_3 = "====================================\n" 
        return info_1+info_2+info_3
    return info_1    
    