'''
    gen_coil.py

'''
import numpy as np

if __name__ == "__main__":
    import coil_util
else:
    from . import coil_util

def create_coil_geom(coil, trace, lead_l, lead_r):
    _pi2 = np.pi*2
    _pi  = np.pi
    
    if trace:
        coil_segment = list()

    p_dist = coil.p_dist
    ncoil = coil.ncoil
    npnt = coil.npnt
    npnt_sub = coil.npnt_sub
    sx, sy = coil.v1[2], coil.v1[3]
    sx_sub, sy_sub = coil.P[2], coil.P[3]
    cx_sub = sx_sub
    r_fillet = coil.r_fillet
    
    p1 = coil.P1[2:]-coil.P1[0:2]
    p2 = coil.P2[2:]-coil.P2[0:2]
    angle_fillet1 = np.atan2(p1[1],p1[0])
    angle_fillet2 = np.atan2(p2[1],p2[0])
    
    if coil.coil_type == coil_util._coil_type_circle:
        axlen = coil.r
        bxlen = coil.r
        t_star1 = angle_fillet1
        t_star2 = angle_fillet2
        t_star3 = angle_fillet1
        t_star4 = angle_fillet2
    else:
        axlen = coil.axlen
        bxlen = coil.bxlen
        t_star1 = coil.t_star
        #if p_dist >= 0:
        if p_dist < 0:
            t_star2 = _pi+(_pi2-t_star1)
        else:
            t_star2 = _pi-t_star1
        t_star3 = angle_fillet1
        t_star4 = angle_fillet2
    
    dist_coil = coil_util._norm_vec(coil._v3)
    dist_sub_coil = coil.r_dist
    
    ncoil_sub = ncoil-1
    total_pnt = npnt*ncoil+npnt_sub*ncoil_sub
    
    if lead_l > 0: total_pnt += 1
    if lead_r > 0: total_pnt += 1
    
    xx = np.zeros(total_pnt) 
    yy = np.zeros(total_pnt)
    
    def create_segment(cx,cy,axlen,bxlen,deg1, deg2, npnt, index):
        ddeg = np.linspace(deg1, deg2, npnt, endpoint=True)
        xx[index:index+npnt] = cx+axlen*np.cos(ddeg)
        yy[index:index+npnt] = cy+bxlen*np.sin(ddeg)
            
        return index+npnt
    
    # First coil
    if coil.coil_type == coil_util._coil_type_circle:
        #start_angle = _pi if p_dist >= 0 else -_pi
        start_angle = _pi if p_dist < 0 else -_pi
        end_angle = t_star1
    else:
        #start_angle = _pi if p_dist >= 0 else -_pi
        start_angle = _pi if p_dist < 0 else -_pi
        #end_angle = t_star1-_pi2 if p_dist >= 0 else t_star1
        end_angle = t_star1-_pi2 if p_dist < 0 else t_star1
        
    if lead_l > 0:
        xx[0] = sx+axlen*np.cos(_pi)-lead_l
        yy[0] = sy
        if trace:
            coil_segment.append((xx[0], yy[0]))
        prv_index = 1
    else:
        prv_index = 0
        
    cur_index = create_segment(sx, 
                               sy, 
                               axlen, 
                               bxlen, 
                               start_angle, 
                               end_angle, 
                               npnt, 
                               prv_index)
    if trace:
        coil_segment.append((xx[prv_index:cur_index], 
                             yy[prv_index:cur_index]))
    prv_index = cur_index
    
    for i in range(1, ncoil-1):
        # Next sub coil
        cx_sub = sx_sub+dist_sub_coil*(i-1)
        if coil.coil_type == coil_util._coil_type_circle or \
           coil.coil_type == coil_util._coil_type_ellipse:
            axlen = r_fillet
            bxlen = r_fillet
            sy_sub = coil.P[3]
            start_angle = t_star3
            end_angle = t_star4
        elif hasattr(coil, "a_s"):
            axlen = coil.a_s
            bxlen = coil.b_s
            sy_sub= coil.yc_s
            start_angle = coil.phi_s
            end_angle = coil.phi_e
            
        cur_index = create_segment(cx_sub, 
                                   sy_sub, 
                                   axlen, 
                                   bxlen, 
                                   start_angle,
                                   end_angle,
                                   npnt_sub, 
                                   prv_index)
        if trace:
            coil_segment.append((xx[prv_index:cur_index], 
                                 yy[prv_index:cur_index]))
        prv_index = cur_index
        
        # Next Coil
        sx += dist_coil
        if coil.coil_type == coil_util._coil_type_circle:
            start_angle = _pi2+angle_fillet2 if p_dist < 0 \
                         else angle_fillet2
            end_angle = angle_fillet1 if p_dist < 0 else _pi2+angle_fillet1

            axlen = coil.r
            bxlen = coil.r
        else:
            start_angle = t_star2 #if p_dist >= 0 else angle_fillet2
            #end_angle = coil.t_star-_pi2 if p_dist >= 0 \
            end_angle = coil.t_star-_pi2 if p_dist < 0 \
                        else _pi2+coil.t_star
            axlen = coil.axlen
            bxlen = coil.bxlen
            
        cur_index = create_segment(sx, 
                                   sy, 
                                   axlen, 
                                   bxlen, 
                                   start_angle, 
                                   end_angle, 
                                   npnt,
                                   prv_index)
        
        if trace:
            coil_segment.append((xx[prv_index:cur_index], 
                                 yy[prv_index:cur_index]))
        prv_index = cur_index
    
    if ncoil_sub > 1:
        # Next sub coil
        cx_sub += dist_sub_coil
    
    if ncoil > 1:
        if coil.coil_type == coil_util._coil_type_circle or \
           coil.coil_type == coil_util._coil_type_ellipse:
            axlen = r_fillet
            bxlen = r_fillet
            sy_sub = coil.P[3]
            start_angle = t_star3
            end_angle = t_star4
        elif hasattr(coil, "a_s"):
            axlen = coil.a_s
            bxlen = coil.b_s
            sy_sub= coil.yc_s
            start_angle = coil.phi_s
            end_angle = coil.phi_e
   
        cur_index = create_segment(cx_sub, 
                                   sy_sub, 
                                   axlen, 
                                   bxlen, 
                                   start_angle,
                                   end_angle,
                                   npnt_sub, 
                                   prv_index)
        
        if trace:
            coil_segment.append((xx[prv_index:cur_index], 
                                 yy[prv_index:cur_index]))
        prv_index = cur_index
                
        sx += dist_coil
        if coil.coil_type == coil_util._coil_type_circle:
            #start_angle = _pi2+angle_fillet2 if p_dist >= 0 \
            start_angle = _pi2+angle_fillet2 if p_dist < 0 \
                         else angle_fillet2
            #end_angle = 0 if p_dist >= 0 else _pi2
            end_angle = 0 if p_dist < 0 else _pi2
            axlen = coil.r
            bxlen = coil.r
        else:
            start_angle = t_star2 
            #end_angle = 0 if p_dist >= 0 else _pi2
            end_angle = 0 if p_dist < 0 else _pi2
            axlen = coil.axlen
            bxlen = coil.bxlen

        cur_index = create_segment( sx, 
                                    sy, 
                                    axlen, 
                                    bxlen, 
                                    start_angle, 
                                    end_angle, 
                                    npnt, 
                                    prv_index)
        if trace:
            coil_segment.append((xx[prv_index:cur_index], 
                                 yy[prv_index:cur_index]))

    if lead_r > 0:
        xx[-1] = xx[-2]+lead_r
        yy[-1] = yy[-2]
        if trace:
            coil_segment.append((xx[-1],yy[-1]))
                                 
    return (xx, yy, coil_segment) if trace else (xx, yy)
                 
                 