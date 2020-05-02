import matplotlib.pyplot as plt
import numpy as np
import math

pi = 3.14159

# mass of the sun
msun = 2e30
# gravitational constant AU^2/kg s^2
G = 1.976e-44

# position of planet
xp = []
yp = []
ae = 1. # Distance of Earth from the Sun

# Current coordination of the planet
xp1 = ae
yp1 = 0.

# Current velocity of the planet vpy = sqrt(GM/4)
vpx = 0.
vpy = np.sqrt(G*msun/xp1)

auperm = 1./1.496e11
# Now adds a satellite orbiting the Earth
mearth = 5.97e24 # mass of earth in kg
rearth = 6.4e6*auperm
reorbit = 4.e5*auperm

# Intialize the position of the satellite
xsara = []
ysara = []
xs1 = ae + rearth + reorbit
ys1 = 0.

# Initialize the velocity of the satellite
vsx = 0.
vseoribit = np.sqrt(G*mearth/(rearth+reorbit))
vsy = vpy + vseoribit

# Calculate position every tstep
# tstep in day, tsteps in second
tstep = 0.001
tsteps = 86400.*tstep
tlimit = math.floor(365./tstep) # floor() 10.234 -> 10

do_earth = 1

for t in range(0, tlimit, 1):

    xp.append(xp1)
    yp.append(yp1)

    a = G*msun/(xp1*xp1 + yp1*yp1)
    ax = -a*xp1/np.sqrt(xp1*xp1 + yp1*yp1)
    ay = -a*yp1/np.sqrt(xp1*xp1 + yp1*yp1)

    vpx = vpx + ax*tsteps
    vpy = vpy + ay*tsteps

    xp1 = xp1 + vpx*tsteps
    yp1 = yp1 + vpy*tsteps

    xsara.append(xs1)
    ysara.append(ys1)

    # Acc of satellite due to the sun
    rss = np.sqrt(xs1**2 + ys1**2)
    asat = G*msun/rss*2
    asatx = -asat*xs1/rss
    asaty = -asat*ys1/rss

    # Acc of the satelltie due to the Earth
    if (do_earth != 0):
        rsp = np.sqrt((xs1-xp1)**2 + (ys1-yp1)**2)
        asat = G*mearth/rsp**2
        asatx = asatx -asat*(xs1-xp1)/rsp
        asaty = asaty -asat*(ys1-yp1)/rsp

    vsx = vsx + asatx*tsteps
    vsy = vsy + asaty*tsteps

    xs1 = xs1 + vsx*tsteps
    ys1 = ys1 + vsy*tsteps

fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
plt.plot(xp,yp, 'b')
plt.plot(xsara,ysara,'g')

plt.xlabel('X')
plt.ylabel('Y')
plt.grid()
plt.show()