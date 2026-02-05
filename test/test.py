import matplotlib.pyplot as plt
import coil_geom as cg

coil = cg.CoilGeom()

xc, yc = coil.circle_coil()
xe, ye = coil.ellipse_coil()

fig, axs = plt.subplots(2)
fig.suptitle('Circle & Ellipse Coil Geometry')
axs[0].plot(xc, yc)
axs[1].plot(xe, ye)

xmax = max(xc)
axs[0].set_xlim(0,xmax)
axs[0].set_aspect('equal')
axs[1].set_xlim(0,xmax)
plt.show()

