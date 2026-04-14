import coil_geom as cg

c_c_up = cg.CircleCoil(p_dist=0.7, ncoil=5)
c_c_dn = cg.CircleCoil(ncoil=5)

c_e_up = cg.EllipseCoil(p_dist=0.4, ncoil=5)
c_e_dn = cg.EllipseCoil(ncoil=5)

c_es_up = cg.EllipseCoilShape(p_dist=0.4, target=0.8, ncoil=5)
c_es_dn = cg.EllipseCoilShape(ncoil=5)

c_ec_up = cg.EllipseCoilCurvature(p_dist=0.4, ncoil=5)
c_ec_dn = cg.EllipseCoilCurvature(ncoil=5)

cg.save_ppt(c_c_up ,  "c_c_up.pptx" )
cg.save_ppt(c_c_dn ,  "c_c_dn.pptx" )
cg.save_ppt(c_e_up ,  "c_e_up.pptx" )
cg.save_ppt(c_e_dn ,  "c_e_dn.pptx" )
cg.save_ppt(c_es_up, "c_es_up.pptx", debug=True)
cg.save_ppt(c_es_dn, "c_es_dn.pptx", debug=True)
cg.save_ppt(c_ec_up, "c_ec_up.pptx", lead_l=2, lead_r=2)
cg.save_ppt(c_ec_dn, "c_ec_dn.pptx")

cg.save_svg(c_c_up ,  "c_c_up.svg" )
cg.save_svg(c_c_dn ,  "c_c_dn.svg" )
cg.save_svg(c_e_up ,  "c_e_up.svg" )
cg.save_svg(c_e_dn ,  "c_e_dn.svg" )
cg.save_svg(c_es_up, "c_es_up.svg", debug=True)
cg.save_svg(c_es_dn, "c_es_dn.svg", debug=True)
cg.save_svg(c_ec_up, "c_ec_up.svg", lead_l=2, lead_r=2)
cg.save_svg(c_ec_dn, "c_ec_dn.svg")

cg.save_pdf(c_c_up ,  "c_c_up.pdf" )
cg.save_pdf(c_c_dn ,  "c_c_dn.pdf" )
cg.save_pdf(c_e_up ,  "c_e_up.pdf" )
cg.save_pdf(c_e_dn ,  "c_e_dn.pdf" )
cg.save_pdf(c_es_up, "c_es_up.pdf", debug=True)
cg.save_pdf(c_es_dn, "c_es_dn.pdf", debug=True)
cg.save_pdf(c_ec_up, "c_ec_up.pdf")
cg.save_pdf(c_ec_dn, "c_ec_dn.pdf")
