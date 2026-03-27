'''

    coil_const.py

'''
import numpy as np
    
_norm_vec = lambda v: np.sqrt(v[0]**2+v[1]**2)
_norm_pnt = lambda v1,v2: np.sqrt(v1**2+v2**2)
_dist_pnt = lambda x1, y1, x2, y2: np.sqrt((x2-x1)**2+(y2-y1)**2)

_circle_coil_name = "Circular Coil"
_ellipse_coil_name = "Ellipic Coil"
_default_npnt = 50
_pi2 = np.pi*2

_optimze_curvature = "curvature"
_optimze_shape = "shape"

_coil_type_name = [
    "Circle Coil",
    "Ellipse Coil(circle)",
    "Ellipse Coil(curvature optimize)",
    "Ellipse Coil(shape optimize)"
]

_circle_coil = 0
_ellipse_coil_circle = 1
_ellipse_coil_optimze_curvature = 2
_ellipse_coil_optimze_shape = 3

_coil_type_circle = _coil_type_name[_circle_coil]
_coil_type_ellipse = _coil_type_name[_ellipse_coil_circle]
_coil_type_ellipse_curvature = _coil_type_name[_ellipse_coil_optimze_curvature]
_coil_type_ellipse_shape = _coil_type_name[_ellipse_coil_optimze_shape]

