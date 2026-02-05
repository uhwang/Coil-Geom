'''
    Drawing Coil 
    02/04/2026
    Uisang hwang
'''

import numpy as np

_norm_vec = lambda v: np.sqrt(v[0]**2+v[1]**2)
_norm_pnt = lambda v1,v2: np.sqrt(v1**2+v2**2)

class CoilGeom():
    def __init__(self, cx=1, cy=1, r=1, 
                       r_dist = 1.3,
                       p_dist = 0.7,
                       ncoil = 5,
                       ):
        
        self.r = r
        self.cx = cx
        self.cy = cy
        self.c1x = cx
        self.c1y = cy
        self.c2x = r+r*r_dist
        self.c2y = cy
        self.r_dist = r_dist
        self.p_dist = p_dist
        self.ncoil = ncoil
        
        self._create()
        
    def _create(self):
        self._v1 = np.array([self.c1x, self.c1y], dtype='float32')
        self._v2 = np.array([self.c2x, self.c2y], dtype='float32')
        self._v3 = self._v2-self._v1
        
        self._v4 = self._v3*0.5
        h = np.sqrt(self.r**2-self._v4[0]**2)
        
        _c = np.cos(np.deg2rad(-90))
        _s = np.sin(np.deg2rad(-90))
        _r = np.array([[_c, -_s], [_s, _c]])
        
        self._v5 = np.dot(_r, self._v4)
        self._v5 = self._v5/np.linalg.norm(self._v5)*h
        
        self.v1 = self._v1
        self.v2 = self._v2
        self.v3 = np.array([self.v1[0],
                            self.v1[1], 
                            self.v1[0]+self._v3[0],
                            self.v1[1]+self._v3[1]])
        self.v4 = np.array([self.v1[0],
                            self.v1[1],
                            self.v1[0]+self._v4[0],
                            self.v1[1]+self._v4[1]])
        self.v5 = np.array([self.v1[0]+self._v4[0],
                            self.v1[1]+self._v4[1],
                            self.v1[0]+self._v4[0]+self._v5[0],
                            self.v1[1]+self._v4[1]+self._v5[1]])
        self.vt = self.v4+self.v5
        self.theta = np.atan2(self.vt[3]-self.vt[1], self.vt[2]-self.vt[0])
                            
        # Find fillet circle
        self._P = self._v5*self.p_dist
        L = _norm_pnt(_norm_vec(self._P),_norm_vec(self._v4))
        R = self.r
        self.r_fillet = np.abs(R-L)
        self.P = np.array([self.v4[2],
                           self.v4[3],
                           self.v4[2]+self._P[0], 
                           self.v4[3]+self._P[1]])
        V = self.P[2:]-self._v1
        uv=V/_norm_vec(V)
        f1 = self.r_fillet*uv
        
        V = self.P[2:]-self._v2
        uv=V/_norm_vec(V)
        f2 = self.r_fillet*uv
        
        self.P1 = np.array([self.P[2], 
                            self.P[3], 
                            self.P[2]+f1[0], 
                            self.P[3]+f1[1]])
        self.P2 = np.array([self.P[2], 
                            self.P[3], 
                            self.P[2]+f2[0], 
                            self.P[3]+f2[1]])

    def create_pnt(self, ncoil = 5, npnt=50, cx=1, cy=1, r=1, r_dist=1.3, p_dist=0.7, trace=False):
    
        check_ncoil = [ncoil < 2, self.ncoil < 2]
        if any(check_ncoil):
            print(f"=> Invalid number of coil: ncoil{ncoil} or self.ncoil{self.ncoil}")
            return
            
        ncoil = ncoil if ncoil != self.ncoil else self.ncoil

        updates = [
            ("cx", cx, cx != self.cx),
            ("cy", cy, cy != self.cy),
            ("r", r, r != self.r),
            ("r_dist", r_dist, r_dist != self.r_dist),
            ("p_dist", p_dist, p_dist != self.p_dist)
        ]
        
        to_update = [item for item in updates if item[2]]

        if to_update:
            for name, value, _ in to_update:
                setattr(self, name, value)
            self._create()
        
        if trace:
            coil_segment = list()
    
        sub_npnt = int(npnt/2)
        
        p1 = self.P1[2:]-self.P1[0:2]
        p2 = self.P2[2:]-self.P2[0:2]
        deg_fillet1 = np.atan2(p1[1],p1[0])
        deg_fillet2 = np.atan2(p2[1],p2[0])
        
        deg_fillet  = deg_fillet2 - deg_fillet1
        dist_coil = _norm_vec(self._v3)
        dist_sub_coil = _norm_vec(self._v4)*2
        
        ncoil_sub = ncoil-1
        total_pnt = npnt*ncoil+sub_npnt*ncoil_sub
        xx = np.zeros(total_pnt) 
        yy = np.zeros(total_pnt)
        
        def circle_pnt(cx,cy,r,deg1, deg2, npnt, index):
            ddeg = (deg2-deg1)/(npnt-1)
            for i in range(npnt):
                deg = deg1+ddeg*i
                x = cx+r*np.cos(deg)
                y = cy+r*np.sin(deg)
                xx[i+index], yy[i+index] = x, y
                
            return i+index+1

        # First coil
        sx, sy = self.v1[0], self.v1[1]
        sx_sub, sy_sub = self.P[2], self.P[3]
        r, r_sub = self.r, self.r_fillet
        cx_sub = sx_sub
        start_angle = np.pi if self.p_dist > 0 else -np.pi
        end_angle = deg_fillet1
        prv_index = 0
        cur_index = circle_pnt(sx, sy, r, start_angle, end_angle, npnt, prv_index)
        
        if trace:
            coil_segment.append((xx[prv_index:cur_index], 
                              yy[prv_index:cur_index]))
        prv_index = cur_index
        
        for i in range(1, ncoil-1):
            # Next sub coil
            cx_sub = sx_sub+dist_sub_coil*(i-1)
            cur_index = circle_pnt(cx_sub, sy_sub, r_sub, 
                                   deg_fillet1, deg_fillet2, 
                                   sub_npnt, prv_index)
            if trace:
                coil_segment.append((xx[prv_index:cur_index], 
                                  yy[prv_index:cur_index]))
            prv_index = cur_index
            
            # Next Coil
            sx += dist_coil
            start_angle = np.pi*2+deg_fillet2 if self.p_dist > 0 else deg_fillet2
            end_angle = deg_fillet1 if self.p_dist > 0 else np.pi*2+deg_fillet1
            cur_index = circle_pnt(sx, sy, r, 
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
            cur_index = circle_pnt(cx_sub, sy_sub, r_sub, 
                                   deg_fillet1, 
                                   deg_fillet2, 
                                   sub_npnt, 
                                   prv_index)
            if trace:
                coil_segment.append((xx[prv_index:cur_index], 
                                  yy[prv_index:cur_index]))
            prv_index = cur_index
                    
            sx += dist_coil
            start_angle = np.pi*2+deg_fillet2 if self.p_dist > 0 else \
                         deg_fillet2
            end_angle = 0 if self.p_dist > 0 else np.pi*2
            cur_index = circle_pnt(sx, sy, r, 
                               start_angle, 
                               end_angle, 
                               npnt, 
                               prv_index)
            if trace:
                coil_segment.append((xx[prv_index:cur_index], 
                                  yy[prv_index:cur_index]))
                                  
        return (xx, yy, coil_segment) if trace else (xx, yy)
        
    def circle_coil(self, ncoil = 5, npnt=50, cx=1, cy=1, r=1, 
                    r_dist = 1.3, p_dist = 0.7, trace = False):    
        return self.create_pnt(ncoil, npnt, cx, cy, r, r_dist, p_dist, trace)
        
    def ellipse_coil(self, ncoil = 5, npnt=50, cx=1, cy=1, r=1, 
                     r_dist = 1.3, p_dist = 0.7, compress=0.5, trace = False):
        if trace:
            xx, yy, seg = self.create_pnt(ncoil, npnt, cx, cy, r, r_dist, p_dist, trace)
            for s in seg:
                s[0] *= compress
        else:
            xx, yy = self.create_pnt(ncoil, npnt, cx, cy, r, r_dist, p_dist, trace)
            
        return (xx*compress, yy, seg) if trace else (xx*compress, yy)
        
    def __str__(self):
        return f"c1x: {self.c1x}\n" \
               f"c1y: {self.c1y}\n" \
               f"c2x: {self.c2x}\n" \
               f"c2y: {self.c2y}\n" \
               f"v1 : {self._v1[0]:2.3f}, {self._v1[1]:2.3f}\n" \
               f"v2 : {self._v2[0]:2.3f}, {self._v2[1]:2.3f}\n" \
               f"v3 : {self._v3[0]:2.3f}, {self._v3[1]:2.3f}\n" \
               f"v4 : {self._v4[0]:2.3f}, {self._v4[1]:2.3f}\n" \
               f"v5 : {self._v5[0]:2.3f}, {self._v5[1]:2.3f}\n" \
               f"V5 : {self.v5[0]:2.3f}, {self.v5[1]:2.3f}\n" \
               f"Th : {self.theta/np.pi*180:2.3f}\n" \
               f"_P : {self._P[0]:2.3f}, {self._P[1]:2.3f}\n" \
               f" P : {self.P[0]:2.3f}, {self.P[1]:2.3f}, {self.P[2]:2.3f}, {self.P[3]:2.3f}\n" \
               f"P1 : {self.P1[0]:2.3f}, {self.P1[1]:2.3f}, {self.P1[2]:2.3f}, {self.P1[3]:2.3f}\n" \
               f"P2 : {self.P2[0]:2.3f}, {self.P2[1]:2.3f}, {self.P2[2]:2.3f}, {self.P2[3]:2.3f}\n" \
               f"R  : {self.r_fillet:2.3f}\n"
               
