'''
    Drawing Coil 
    
    02/04/2026 Initial Version
    02/18/2026 Refactor CoilGeom class
    04/25/2026 Added geometric transformation methods 
               (rotate, flipud, fliplr) for enhanced 
               coil design flexibility.
               
               Implemented Deferred Transformation in CoilGeom 
               to support .rotate() and .flip() 
               before .create_geom() is called.
    Uisang hwang
'''

import numpy as np
import copy
from scipy.optimize import minimize_scalar, fsolve

_pi = np.pi
_pi2= np.pi*2
_d2r= _pi/180

if __name__ == "__main__":
    import coil_gen
    import coil_util
    import coil_print
else:
    from . import coil_gen
    from . import coil_util
    from . import coil_print

def get_coil_type(coil):
    if isinstance(coil, CircleCoil): 
        return coil_util._coil_type_circle
    elif isinstance(coil, EllipseCoilCurvature): 
        return coil_util._coil_type_ellipse_curvature
    elif isinstance(coil, EllipseCoilShape):
        return coil_util._coil_type_ellipse_shape
    elif isinstance(coil, EllipseCoil): 
        return coil_util._coil_type_ellipse
    else: 
        return "Invalid Coil Type"

class CoilData:
    def __init__(self, *args, tasks=None):
        self.x_ = args[0]
        self.y_ = args[1]
        """
        args[0]: x 좌표 데이터
        args[1]: y 좌표 데이터
        args[2]: seg (trace 데이터, 선택사항)
        tasks: 지연 계산을 위한 명령 리스트
        """
        if tasks is None:
            # 최초 생성 시에만 데이터를 복사하여 원본을 안전하게 보관 (Deep Copy)
            self._raw_x = np.array(args[0], copy=True)
            self._raw_y = np.array(args[1], copy=True)
            if len(args) > 2 and args[2] is not None:
                self._raw_seg = copy.deepcopy(args[2])
            else:
                self._raw_seg = None
            self.tasks = []
        else:
            # 변환 체이닝 중에는 복사 없이 원본 데이터의 참조만 전달 (메모리 효율 극대화)
            self._raw_x = args[0]
            self._raw_y = args[1]
            self._raw_seg = args[2] if len(args) > 2 else None
            self.tasks = tasks

        # 계산 결과 캐싱 (필요할 때 한 번만 계산)
        self._cached_data = None

    def get(self):
        return self.x_, self.y_
        
    def rotate(self, deg, axis=1):
        new_tasks = self.tasks + [('rotate', deg, axis)]
        return CoilData(self._raw_x, self._raw_y, self._raw_seg, tasks=new_tasks)

    def flipud(self):
        new_tasks = self.tasks + [('flipud',)]
        return CoilData(self._raw_x, self._raw_y, self._raw_seg, tasks=new_tasks)

    def fliplr(self, axis=1):
        new_tasks = self.tasks + [('fliplr', axis)]
        return CoilData(self._raw_x, self._raw_y, self._raw_seg, tasks=new_tasks)

    @property
    def data(self):
        if self._cached_data is None:
            self._cached_data = self._apply_all_tasks()
        return self._cached_data

    @property
    def x(self): return self.data[0]

    @property
    def y(self): return self.data[1]

    @property
    def seg(self): return self.data[2]

    @property
    def x_c(self):
        return (np.max(self.x) + np.min(self.x)) * 0.5

    @property
    def y_c(self):
        return (np.max(self.y) + np.min(self.y)) * 0.5

    def __iter__(self):
        return zip(self.x, self.y)

    def __len__(self):
        return len(self.x)

    def _update_center(self):
        pass

    def _apply_all_tasks(self):
        curr_x = np.copy(self._raw_x)
        curr_y = np.copy(self._raw_y)
        curr_seg = copy.deepcopy(self._raw_seg)

        for task in self.tasks:
            op = task[0]
            if op == 'rotate':
                curr_x, curr_y, curr_seg = self._exec_rotate(curr_x, curr_y, curr_seg, *task[1:])
            elif op == 'flipud':
                curr_x, curr_y, curr_seg = self._exec_flipud(curr_x, curr_y, curr_seg)
            elif op == 'fliplr':
                curr_x, curr_y, curr_seg = self._exec_fliplr(curr_x, curr_y, curr_seg, *task[1:])
        
        return curr_x, curr_y, curr_seg

    def _exec_rotate(self, x, y, seg, deg, axis):
        theta = deg * _d2r
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        
        xc = (np.max(x) + np.min(x)) * 0.5
        yc = (np.max(y) + np.min(y)) * 0.5

        if axis == 0:
            ref_x, ref_y = x[0], y[0]
        elif axis == 2:
            ref_x, ref_y = x[-1], y[-1]
        else:
            ref_x, ref_y = xc, yc

        new_x = (x - ref_x) * cos_t - (y - ref_y) * sin_t + ref_x
        new_y = (x - ref_x) * sin_t + (y - ref_y) * cos_t + ref_y
        
        new_seg = None
        if seg:
            new_seg = [
                ((sx - ref_x) * cos_t - (sy - ref_y) * sin_t + ref_x,
                 (sx - ref_x) * sin_t + (sy - ref_y) * cos_t + ref_y) 
                for sx, sy in seg
            ]
        return new_x, new_y, new_seg

    def _exec_flipud(self, x, y, seg):
        new_y = 2 * y[0] - y
        new_seg = [(sx, 2 * yc - sy) for sx, sy in seg] if seg else None
        return x, new_y, new_seg

    def _exec_fliplr(self, x, y, seg, axis):
        xc = (np.max(x) + np.min(x)) * 0.5
        
        if axis == 0:
            ref_x = x[0]
        elif axis == 2:
            ref_x = x[-1]
        else: # axis == 1
            ref_x = xc

        new_x = 2 * ref_x - x
        
        new_seg = None
        if seg:
            new_seg = [(2 * ref_x - sx, sy) for sx, sy in seg]
            
        return new_x, y, new_seg
            
class VectorDiagram():
    def __init__(self):
        self._v1 = np.zeros(2)
        self._v2 = np.zeros(2)
        self._v3 = np.zeros(2)
        self._v4 = np.zeros(2)
        self._v5 = np.zeros(2)
        self._P  = np.zeros(2)
        self.v1  = np.zeros(4)
        self.v2  = np.zeros(4)
        self.v3  = np.zeros(4)
        self.v4  = np.zeros(4)
        self.v5  = np.zeros(4)
        self.P   = np.zeros(4)
        self.P1  = np.zeros(4)
        self.P2  = np.zeros(4)
        
class CoilGeom():
    def __init__( self, 
                  cx, 
                  cy,
                  r_dist,
                  p_dist,
                  ncoil,
                  npnt,
                  npnt_sub
                ):
                
        self.cx = cx
        self.cy = cy
        self.c1x = cx
        self.c1y = cy
        self.c2x = cx+r_dist
        self.c2y = cy
        self.r_dist = r_dist
        self.p_dist = p_dist
        self.ncoil = ncoil
        self.npnt = npnt
        self.npnt_sub = npnt_sub
        self._pending_transform = {
            'rot': 0, 
            'axis': 1, 
            'flip': None,
            'flip_axis': 1  # fliplr을 위한 축 설정 추가
        }

    def rotate(self, deg, axis=1):
        self._pending_transform['rot'] = deg
        self._pending_transform['axis'] = axis
        return self
    
    def flipud(self):
        self._pending_transform['flip'] = 'ud'
        return self
    
    def fliplr(self, axis=1):
        self._pending_transform['flip'] = 'lr'
        self._pending_transform['flip_axis'] = axis
        return self
    
    def apply_pending(self, coil_data):
        """CoilData 객체에 예약된 작업을 MoviePy 스타일로 이관"""
        out = coil_data
        
        # 회전 예약 적용 (실제 계산은 지연됨)
        if self._pending_transform['rot'] != 0:
            out = out.rotate(
                self._pending_transform['rot'], 
                axis=self._pending_transform['axis']
            )
        
        # 대칭 이동 예약 적용 (실제 계산은 지연됨)
        if self._pending_transform['flip'] == 'ud':
            out = out.flipud()
        elif self._pending_transform['flip'] == 'lr':
            out = out.fliplr(axis=self._pending_transform['flip_axis'])
            
        return out
        
    def get_geom(self, trace=False, lead_l=0, lead_r=0):
        return coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
        
    def create_geom(self, trace=False, lead_l=0, lead_r=0):
        result = coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
        
        data = CoilData(*result)

        return self.apply_pending(data)  
        
class EllipseCoilOptimize():
    def __init__(self):
        self.a_s = 0.0
        self.b_s = 0.0
        self.y_s = 0.0
        self.xc_s = 0.0
        self.yc_s = 0.0
        self.phi_s = 0.0
        self.phi_e = 0.0
        self.delta = 0.0
        
class EllipseCoil(CoilGeom, VectorDiagram):
    def __init__(self, cx = 0, 
                       cy = 0, 
                       axlen = 1, 
                       bxlen = 2, 
                       r_dist = 1.3, # distance from cx
                       p_dist = -0.4,
                       ncoil = 2,
                       npnt = 50,
                       npnt_sub = 25
                       ):
        VectorDiagram.__init__(self)
        CoilGeom.__init__( self, 
                           cx, 
                           cy, 
                           r_dist, 
                           p_dist, 
                           ncoil, 
                           npnt, 
                           npnt_sub )
                           
        self.coil_type = get_coil_type(self)
        self.axlen = axlen
        self.bxlen = bxlen
        self.create_vector()
        
    def create_vector(self):

        xmid = (self.c2x-self.c1x)*0.5
        cx, cy = self.c1x, self.c1y
        a,b = self.axlen, self.bxlen
        hh = b*np.sqrt(1-xmid**2/a**2)

        self._v1[:] = (self.c1x, self.c1y)
        self._v2[:] = (self.c2x, self.c2y)
        self._v3[:] = self._v2-self._v1
        self._v4[:] = self._v3*0.5
        if self.p_dist < 0:
            self._v5[:] = (self._v4[1], -self._v4[0])
        else:
            self._v5[:] = (-self._v4[1], self._v4[0])
        self._v5[:] = self._v5/coil_util._norm_vec(self._v5)*hh
        self._P [:] = self._v5*abs(self.p_dist)

        self.v1[:] = ( 0,0,self._v1[0], self._v1[1] ) 
        self.v2[:] = ( 0,0,self._v2[0], self._v2[1] ) 
        self.v3[:] = ( self.v1[2],
                       self.v1[3], 
                       self.v1[2]+self._v3[0],
                       self.v1[3]+self._v3[1] )
        self.v4[:] = ( self.v1[2], 
                       self.v1[3],
                       self.v1[2]+self._v4[0],
                       self.v1[3]+self._v4[1])
        self.v5[:] = ( self.v1[2]+self._v4[0],
                       self.v1[3]+self._v4[1],
                       self.v1[2]+self._v4[0]+self._v5[0],
                       self.v1[3]+self._v4[1]+self._v5[1] )
        self.P[:]  = ( self.v4[2],
                       self.v4[3],
                       self.v4[2]+self._P[0], 
                       self.v4[3]+self._P[1] )
        
        # Objective: Minimize squared distance from P to Ellipse
        def dist_sq(t):
            tx = self.c1x+self.axlen * np.cos(t)
            ty = self.c1y+self.bxlen * np.sin(t)
            return (tx - self.P[2])**2 + (ty - self.P[3])**2
        
        # Search bottom arc (pi to 2*pi) for the fillet connection
        # Search upper arc (0 to pi) for the fillet connection
        search_range=(_pi, _pi2) if self.p_dist < 0 else (0, _pi) 
        res = minimize_scalar(dist_sq, bounds=search_range, method='bounded')

        self.t_star = res.x
        self.r_fillet = np.sqrt(res.fun)
        tx = self.c1x+self.axlen*np.cos(self.t_star)-self.P[2]
        ty = self.c1x+self.bxlen*np.sin(self.t_star)-self.P[3]
    
        self.P1[:] = ( self.P[2], 
                       self.P[3], 
                       self.P[2]+tx, 
                       self.P[3]+ty )
        self.P2[:] = ( self.P[2], 
                       self.P[3], 
                       self.P[2]-tx, 
                       self.P[3]+ty )
    
    def set_param(self, cx = 0,
                        cy = 0,
                        axlen=1, 
                        bxlen=2, 
                        r_dist = 1.3, 
                        p_dist = -0.4, 
                        ncoil = 5, 
                        npnt = 50, 
                        npnt_sub = 25):
                        
        check_ncoil = [ncoil < 2, self.ncoil < 2]
        if any(check_ncoil):
            print(f"=> EllipseCoil: Invalid number of coil (ncoil:{ncoil} or self.ncoil:{self.ncoil})")
            return
            
        updates = [
            ("cx", cx, cx != self.cx),
            ("cy", cy, cy != self.cy),
            ("axlen", axlen, axlen != self.axlen),
            ("bxlen", bxlen, bxlen != self.bxlen),
            ("r_dist", r_dist, r_dist != self.r_dist),
            ("p_dist", p_dist, p_dist != self.p_dist)
        ]
        
        to_update = [item for item in updates if item[2]]

        if to_update:
            for name, value, _ in to_update:
                setattr(self, name, value)
            self.create_vector()
            if isinstance(self, EllipseCoilOptimize):
                optimize_vertex(self)
        
        self.ncoil = ncoil if ncoil != self.ncoil else self.ncoil
        self.npnt = npnt if npnt != self.npnt else self.npnt
        self.npnt_sub = npnt_sub if npnt_sub != self.npnt_sub else self.npnt_sub    
       
    '''
    def create_geom(self, trace = False, lead_l=0, lead_r=0):
        #return coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
        result = coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
        data = CoilData(*result)
        return self.apply_pending(data)
    '''
    
    def __str__(self):
        return coil_print.print_coil(self)

# find Optimum Y(ys) coord of transition ellipse 
def optimize_vertex(self):
    """
    target: 0 - 1.0
    mode: curvature or shape
    """
    # 1. Basic Geometric Parameters
    a, b, t_s = self.axlen, self.bxlen, self.t_star
    xt, yt = a*np.cos(t_s), b*np.sin(t_s)
    m = -(b**2*xt)/(a**2*yt)
    xs, yc_c = self.P[2], self.P[3]
    r_f = self.r_fillet 
    is_up = -1 if self.p_dist >= 0 else 1
    dx = xt-xs
    rho_main = ((a*np.sin(t_s))**2+(b*np.cos(t_s))**2)**1.5/(a*b)
    
    def get_asbs_c(ys):
        dy = yt - ys
        b_sq = dy**2-m*dy*dx
        a_sq = -(b_sq*dx)/(m*dy) if abs(m*dy) > 1e-6 else -1
        return np.sqrt(a_sq), np.sqrt(max(0,b_sq))

    def get_asbs_s(ys):
        dy = yt - ys
        b_s = (m*dx*dy-dy**2)/(2*dy-m*dx)
        a_s = np.sqrt(max(0,-(dx**2)/((dy**2/b_s**2)+(2*dy/b_s))))
        return a_s, b_s
        
    asbs_func = get_asbs_c if self.target_mode == coil_util._optimze_curvature \
               else get_asbs_s
    
    def obj_curv(delta):
        #ys_tmp = yc_c + r_f + delta
        ys_tmp = yc_c + delta
        a_s, b_s = asbs_func(ys_tmp)
        dy = yt-ys_tmp
        phi_t = np.arctan2(dy/b_s, dx/a_s)
        rho_b = ((a_s*np.sin(phi_t))**2 + (b_s*np.cos(phi_t))**2)**1.5 / (a_s*b_s)
        return (min(rho_main, rho_b) / max(rho_main, rho_b)) - self.target
        
    def obj_shape(delta):
        ys_tmp = yc_c + r_f + delta
        a_s, b_s = asbs_func(ys_tmp)
        cur_shape = (min(b/a, b_s/a_s) / max(b/a, b_s/a_s))
        return cur_shape-self.target
        
    # 2. solve 
    obj_func = obj_curv if self.target_mode == coil_util._optimze_curvature \
               else obj_shape
               
    delta_found = False
    
    if self.target_mode == coil_util._optimze_curvature:
        if self.p_dist >= 0: 
            bounds=(-self.right_bound, -self.left_bound)
        else:
            bounds=(self.left_bound, self.right_bound)
        res = minimize_scalar(lambda d: obj_func(d)**2, bounds=bounds, method='bounded')
        if res.success:
            self.delta = res.x    
            delta_found = True
    else:
        initial_guess = self.initial_guess
        sol, info, ier, msg = fsolve(obj_func, initial_guess, full_output=True)
        if ier == 1: # converge OK
            self.delta = sol[0]   
            delta_found = True
            
    if delta_found:
        ys = (yc_c + r_f) + self.delta
        a_s, b_s = get_asbs_s(ys)
        #a_s, b_s = asbs_func(ys)
        self.a_s = a_s
        self.b_s = b_s
        self.y_s = ys
        self.xc_s = (self.c2x-self.c1x)*0.5
        self.yc_s = ys-self.b_s
        # properties of transition ellipse
        dy = yt - ys + self.b_s
        self.phi_s = np.arctan2(dy/self.b_s, dx/self.a_s)
        self.phi_e = np.arctan2(dy/self.b_s, (self.r_dist-xt-xs)/self.a_s)
    else:
        # Fail
        print(f"[Warning] fsolve failed to converge: {msg}")
        print("Fallback to default circular fillet (offset=0).")
        self.applied_offset = 0.0 
        self.a_s = self.r_fillet
        self.b_s = self.r_fillet
        self.xc_s = (self.c2x-self.c1x)*0.5
        self.yc_s = self.P[3]
        self.y_s = yc_c + self.r_fillet

class EllipseCoilCurvature(EllipseCoil, EllipseCoilOptimize):
    def __init__(self,cx = 0, 
                      cy = 0, 
                      axlen = 1, 
                      bxlen = 2, 
                      r_dist = 1.3,
                      p_dist = -0.4,
                      ncoil = 2,
                      npnt = 50,
                      npnt_sub = 25,
                      target = 0.5,
                      left_bound=0.01,
                      right_bound=3):
        
        EllipseCoilOptimize.__init__(self)
        EllipseCoil.__init__(self, cx,cy,axlen,bxlen,r_dist,p_dist,ncoil,npnt,npnt_sub)
        
        self.coil_type = get_coil_type(self)
        self.target_mode = coil_util._optimze_curvature
        self.target = target
        self.left_bound = left_bound
        self.right_bound = right_bound
        self.create_vector()
        optimize_vertex(self)

class EllipseCoilShape(EllipseCoil, EllipseCoilOptimize):
    def __init__(self,cx = 0, 
                      cy = 0, 
                      axlen = 1, 
                      bxlen = 2, 
                      r_dist = 1.3,
                      p_dist = -0.4,
                      ncoil = 2,
                      npnt = 50,
                      npnt_sub = 25,
                      target = 0.9,
                      initial_guess=0.01):
        
        EllipseCoilOptimize.__init__(self)              
        EllipseCoil.__init__(self, cx,cy,axlen,bxlen,r_dist,p_dist,ncoil,npnt,npnt_sub)
        
        self.coil_type = get_coil_type(self)
        self.target_mode = coil_util._optimze_shape
        self.target = target
        self.initial_guess = initial_guess
        self.create_vector()
        optimize_vertex(self)

class CircleCoil(CoilGeom, VectorDiagram):
    def __init__(self, cx=1, 
                       cy=1, 
                       r=1, 
                       r_dist = 1.3,
                       p_dist = -0.7,
                       ncoil = 2,
                       npnt = 50,
                       npnt_sub = 25
                       ):
     
        VectorDiagram.__init__(self)
        CoilGeom.__init__( self, 
                           cx, cy, 
                           r_dist, 
                           p_dist, 
                           ncoil, 
                           npnt, 
                           npnt_sub )
                                
        self.coil_type = get_coil_type(self)
        self.r = r
        self.axlen = r
        self.bxlen = r
        self.create_vector()
        
    def create_vector(self):
        self._v1[:] = (self.c1x, self.c1y)
        self._v2[:] = (self.c2x, self.c2y)
        self._v3[:] = self._v2-self._v1
        self._v4[:] = self._v3*0.5
        
        hh = np.sqrt(self.r**2-coil_util._norm_vec(self._v4)**2)
        
        if self.p_dist < 0:
            self._v5[:] = (self._v4[1], -self._v4[0])
        else:
            self._v5[:] = (-self._v4[1], self._v4[0])        
        #self._v5[:] = (self._v4[1], -self._v4[0])
        self._v5[:] = self._v5/coil_util._norm_vec(self._v5)*hh
        self._P [:] = self._v5*abs(self.p_dist)
        
        self.v1[:] = ( 0,0,self._v1[0], self._v1[1] ) 
        self.v2[:] = ( 0,0,self._v2[0], self._v2[1] ) 
        self.v3[:] = ( self.v1[2],
                       self.v1[3], 
                       self.v1[2]+self._v3[0],
                       self.v1[3]+self._v3[1] )
        self.v4[:] = ( self.v1[2], 
                       self.v1[3],
                       self.v1[2]+self._v4[0],
                       self.v1[3]+self._v4[1])
        self.v5[:] = ( self.v1[2]+self._v4[0],
                       self.v1[3]+self._v4[1],
                       self.v1[2]+self._v4[0]+self._v5[0],
                       self.v1[3]+self._v4[1]+self._v5[1] )
        self.P[:]  = ( self.v4[2],
                       self.v4[3],
                       self.v4[2]+self._P[0], 
                       self.v4[3]+self._P[1] )
        
        P, V1, V2, R = self.P[2:], self.v1[2:], self.v2[2:], self.r
        
        L = P - V1
        L_len = coil_util._norm_vec(L)
        self.r_fillet = np.abs(R-L_len)
        uv=L/L_len
        f1 = self.r_fillet*uv
        
        L = P - V2
        L_len = coil_util._norm_vec(L)
        uv=L/L_len
        f2 = self.r_fillet*uv
    
        self.P1[:] = ( self.P[2], 
                       self.P[3], 
                       self.P[2]+f1[0], 
                       self.P[3]+f1[1] )
        self.P2[:] = ( self.P[2], 
                       self.P[3], 
                       self.P[2]+f2[0], 
                       self.P[3]+f2[1])
        p1 = self.P1[2:]-self.P1[0:2]
        self.t_star = np.atan2(p1[1],p1[0])
    
    def set_param(self, cx = 1,
                        cy = 1,
                        r = 1,
                        r_dist = 1.3, 
                        p_dist = -0.7, 
                        ncoil = 2, 
                        npnt = 50, 
                        npnt_sub = 25,
                        trace = False):
    
            check_ncoil = [ncoil < 2, self.ncoil < 2]
            if any(check_ncoil):
                print(f"=> Invalid number of coil: ncoil{ncoil} or self.ncoil{self.ncoil}")
                return
                
            updates = [
                ("cx", cx, cx != self.cx),
                ("cy", cy, cy != self.cy),
                ("r", r, r != self.r),
                #("ncoil", ncoil, ncoil != self.ncoil),
                ("r_dist", r_dist, r_dist != self.r_dist),
                ("p_dist", p_dist, p_dist != self.p_dist)
            ]
            
            to_update = [item for item in updates if item[2]]
    
            if to_update:
                for name, value, _ in to_update:
                    setattr(self, name, value)
                self.create_vector()
                
            self.ncoil = ncoil if ncoil != self.ncoil else self.ncoil
            self.npnt = npnt if npnt != self.npnt else self.npnt
            self.npnt_sub = npnt_sub if npnt_sub != self.npnt_sub else self.npnt_sub        
            #return coil_gen.create_coil_geom(self, trace)

    '''
    def create_geom(self, trace = False, lead_l=0, lead_r=0):
        result = coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
        data = CoilData(*result)
        return self.apply_pending(data)
        #return coil_gen.create_coil_geom(self, trace, lead_l, lead_r)
    '''
    
    def __str__(self):
        return coil_print.print_coil(self)

if __name__ == "__main__":
    cc = EllipseCoilCurvature(axlen=2.5,bxlen=6.0,r_dist=3.527)
    #cc.create_geom(axlen=2.5,bxlen=6.0,r_dist=3.527)
    #cc = EllipseCoilShape()