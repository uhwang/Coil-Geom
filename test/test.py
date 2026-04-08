import matplotlib.pyplot as plt
import coil_geom as cg

c_c_up = cg.CircleCoil(p_dist=0.7)
c_c_dn = cg.CircleCoil()

c_e_up = cg.EllipseCoil(p_dist=0.4)
c_e_dn = cg.EllipseCoil()

c_es_up = cg.EllipseCoilShape(p_dist=0.4)
c_es_dn = cg.EllipseCoilShape()

c_ec_up = cg.EllipseCoilCurvature(p_dist=0.4)
c_ec_dn = cg.EllipseCoilCurvature()

cg.save_ppt(c_c_up , "c_c_up.pptx" )
cg.save_ppt(c_c_dn , "c_c_dn.pptx" )

cg.save_ppt(c_e_up , "c_e_up.pptx" )
cg.save_ppt(c_e_dn , "c_e_dn.pptx" )

cg.save_ppt(c_es_up, "c_es_up.pptx")
cg.save_ppt(c_es_dn, "c_es_dn.pptx")

cg.save_ppt(c_ec_up, "c_ec_up.pptx")
cg.save_ppt(c_ec_dn, "c_ec_dn.pptx")