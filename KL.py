import matplotlib.pyplot as plt
import numpy as np
import math

'''Constants'''
# pi
pi = 3.14159
# gravitational constant AU^2/kg s^2
G = 1.976e-44
# Astronomical Unit per meter
auperm = 1. / 1.496e11
# mass of common stars in kg
mSun = 1.98e30
# mass of common planets in kg
mJupiter = 1.90e27
mEarth = 5.97e24
# orbital distance of common planets in AU
rJupiter = 778330000000 * auperm  # 5.202 AU
# rJupiter = 1.1
rEarth = 1.

'''Masses'''
# mass of the star
ms = mSun
# mass of planet p
mp = 0.
# mass of planet q
mq = mJupiter

'''Orbit radius'''
# radius of planet p
rp = rEarth
# radius of planet q
rq = rJupiter
# Period of planet P
Tp = 1.
# Period of planet Q
Tq = 1.

# Current coordination of the planet p
xp1 = rp * 0.99
yp1 = rp * 0.01
# xp1 = rEarth
# yp1 = 0.
zp1 = np.sqrt(rp ** 2 - xp1 ** 2 - yp1 ** 2)

# Current coordination of the planet q
xq1 = rq * 0.98
yq1 = rq * 0.01
# xq1 = rJupiter
# yq1 = 0.
zq1 = np.sqrt(rq ** 2 - xq1 ** 2 - yq1 ** 2)

# Current velocity of the planet p         vpy = sqrt(GM/4)
vpx = 0.
vpy = np.sqrt((G * ms) / rp) * 0.99
vpz = np.sqrt((G * ms) / rp - vpy ** 2)
# vpy = np.sqrt((G * ms) / rp)
# vpz = 0.
vp = np.sqrt(vpx ** 2 + vpy ** 2 + vpz ** 2)

# Current velocity of the planet q         vpy = sqrt(GM/4)
vqx = 0.
vqy = np.sqrt((G * ms) / rq) * 0.98
vqz = np.sqrt((G * ms) / rq - vqy ** 2)
# vqy = np.sqrt((G * ms) / rq)
# vqz = 0.
vq = np.sqrt(vqx ** 2 + vqy ** 2 + vqz ** 2)

# Specific Angular Momentum of planet P and Q
hp = np.sqrt((yp1 * vpz - vpy * zp1) ** 2 + (vpx * zp1 - xp1 * vpz) ** 2 + (xp1 * vpy - vpx * yp1) ** 2)
hq = np.sqrt((yq1 * vqz - vqy * zq1) ** 2 + (vqx * zq1 - xq1 * vqz) ** 2 + (xq1 * vqy - vqx * yq1) ** 2)

# Total energy of planet P and Q
Ep = 0.5 * mp * vp ** 2 - G * ms * mp / rp
Eq = 0.5 * mq * vq ** 2 - G * ms * mq / rq

# Time list
time = []
t = 0.

# position of planet p list
xp = []
yp = []
zp = []

# position of planet q list
xq = []
yq = []
zq = []

# inclination of planet p list
ip = []

# inclination of planet q list
iq = []

# eccentricity of planet P list
ep = []

# eccentricity of planet Q list
eq = []

# /sqrt (1-e^2)
edp = []

# 1/cosi
idelta = []

# Lz
Lz = []
SumLz = 0.

# Calculate position every tstep
# tstep in day, tsteps in second
tstep = 1.
tsteps = 86400. * tstep
tlimit = math.floor(25*365. / tstep)  # floor() 10.234 -> 10

for t in range(0, tlimit, 1):
    # Add updated time to list
    time.append(t)

    # Add updated position of planet P to list
    xp.append(xp1)
    yp.append(yp1)
    zp.append(zp1)

    # Add updated position of planet Q to list
    xq.append(xq1)
    yq.append(yq1)
    zq.append(zq1)

    # Add updated inclination to list
    ip1 = np.arctan(zp1 / np.sqrt(xp1 * xp1 + yp1 * yp1))
    iq1 = np.arctan(zq1 / np.sqrt(xq1 * xq1 + yq1 * yq1))

    ip.append(ip1 * 180 / pi)
    iq.append(iq1 * 180 / pi)

    # Add updated eccentricity to list
    ep1sqr = 1 + (2 * (0.5 * vp ** 2 - G * ms / rp)) * (hp / (G * ms)) ** 2
    eq1sqr = 1 + (2 * (0.5 * vq ** 2 - G * ms / rq)) * (hq / (G * ms)) ** 2

    ep.append(np.sqrt(ep1sqr))
    eq.append(np.sqrt(eq1sqr))

    # Calculate and add updated Lz to list
    Lz1 = np.sqrt(1 - ep1sqr) * (np.cos(np.absolute(ip1)+np.absolute(iq1)))
    Lz.append(Lz1)
    SumLz += Lz1

    # Distance from P to S
    Dps = np.sqrt(xp1 ** 2 + yp1 ** 2 + zp1 ** 2)
    # Distance from P to Q
    Dpq = np.sqrt((xp1 - xq1) ** 2 + (yp1 - yq1) ** 2 + (zp1 - zq1) ** 2)
    # Distance from Q to S
    Dqs = np.sqrt(xq1 ** 2 + yq1 ** 2 + zq1 ** 2)

    # Update acc of planet P due to S and Q
    aPs = G * ms / Dps ** 2
    aPq = G * mq / Dpq ** 2

    # Update acc components of planet P
    aPx = -aPs * xp1 / Dps - aPq * (xp1 - xq1) / Dpq
    aPy = -aPs * yp1 / Dps - aPq * (yp1 - yq1) / Dpq
    aPz = -aPs * zp1 / Dps - aPq * (zp1 - zq1) / Dpq

    # Update acc of planet Q due to S and P
    aQs = G * ms / Dqs ** 2
    aQp = G * mp / Dpq ** 2

    # Update acc components of planet Q
    aQx = -aQs * xq1 / Dqs - aQp * (xq1 - xp1) / Dpq
    aQy = -aQs * yq1 / Dqs - aQp * (yq1 - yp1) / Dpq
    aQz = -aQs * zq1 / Dqs - aQp * (zq1 - zp1) / Dpq

    # Update velocity components of planet P
    vpx = vpx + aPx * tsteps
    vpy = vpy + aPy * tsteps
    vpz = vpz + aPz * tsteps

    # Update velocity components of planet Q
    vqx = vqx + aQx * tsteps
    vqy = vqy + aQy * tsteps
    vqz = vqz + aQz * tsteps

    # Update position components of planet P
    xp1 = xp1 + vpx * tsteps
    yp1 = yp1 + vpy * tsteps
    zp1 = zp1 + vpz * tsteps

    # Update position components of planet Q
    xq1 = xq1 + vqx * tsteps
    yq1 = yq1 + vqy * tsteps
    zq1 = zq1 + vqz * tsteps

    # Update combine velocity of planet P and Q
    vp = np.sqrt(vpx ** 2 + vpy ** 2 + vpz ** 2)
    vq = np.sqrt(vqx ** 2 + vqy ** 2 + vqz ** 2)

    # Update radius of orbit of planet P and Q
    rp = np.sqrt(xp1 ** 2 + yp1 ** 2 + zp1 ** 2)
    rq = np.sqrt(xq1 ** 2 + yq1 ** 2 + zq1 ** 2)

    # Update total energy of planet P and Q
    Ep = 0.5 * mp * vp ** 2 - G * ms * mp / rp
    Eq = 0.5 * mq * vq ** 2 - G * ms * mq / rq

    # Update specific angular momentum of planet P and Q , h = v x r
    hp = np.sqrt((yp1 * vpz - vpy * zp1) ** 2 + (vpx * zp1 - xp1 * vpz) ** 2 + (xp1 * vpy - vpx * yp1) ** 2)
    hq = np.sqrt((yq1 * vqz - vqy * zq1) ** 2 + (vqx * zq1 - xq1 * vqz) ** 2 + (xq1 * vqy - vqx * yq1) ** 2)

    # Update Time
    t += tsteps
    print(t)

LzAverage = SumLz / tlimit
print(LzAverage)

plt.figure()

ax1 = plt.subplot2grid((2, 2), (0, 0), colspan = 1, rowspan = 1)
plt.plot(time, ip, 'b')
plt.plot(time, iq, 'r')
ax1.set_title('Inclination versus Time')
plt.xlabel('Time')
plt.ylabel('Inclination')
plt.grid()

ax2 = plt.subplot2grid((2, 2), (0, 1), colspan = 1, rowspan = 1)
plt.plot(time, ep, 'b')
plt.plot(time, eq, 'r')
ax2.set_title('Eccentricity versus Time')
plt.xlabel('Time')
plt.ylabel('Eccentricity')
plt.grid()

'''
ax3 = plt.subplot2grid((3, 2), (1, 0), colspan = 1, rowspan = 1)
plt.plot(xp, yp, 'b')
plt.plot(xq, yq, 'r')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid()

ax4 = plt.subplot2grid((3, 2), (1, 1), colspan = 1, rowspan = 1)
plt.plot(xp, zp, 'b')
plt.plot(xq, zq, 'r')
plt.xlabel('X')
plt.ylabel('Z')
plt.grid()

ax5 = plt.subplot2grid((3, 2), (2, 0), colspan = 1, rowspan = 1)
plt.plot(yp, zp, 'b')
plt.plot(yq, zq, 'r')
plt.xlabel('Y')
plt.ylabel('Z')
plt.grid()
'''

ax6 = plt.subplot2grid((2, 2), (1, 0), colspan = 2, rowspan = 1)
plt.plot(time, Lz, 'g')
ax6.axhline(LzAverage, ls = '--')
ax6.set_title("Lz versus Time")
ax6.annotate('Average Lz: ' + str(round(LzAverage,3)), xy=(9300,0.96))
plt.xlabel('time')
plt.ylabel('Lz')
plt.grid()



plt.show()

figx = plt.figure()
axx = figx.add_subplot()
plt.plot(time, Lz, 'g')
axx.axhline(LzAverage, ls = '--')
axx.set_title("Lz versus Time")
axx.annotate('Average Lz: ' + str(round(LzAverage,3)), xy=(9300,0.96))
plt.xlabel('time')
plt.ylabel('Lz')
plt.grid()

plt.show()

'''
fig = plt.figure()

# Inclination versus Time
plt.plot(time, ip, 'b')
plt.plot(time, iq, 'r')

plt.xlabel('Time')
plt.ylabel('Inclination')
plt.grid()
plt.show()

fig2 = plt.figure()

# Eccentricity versus Time
plt.plot(time, ep, 'b')
plt.plot(time, eq, 'r')

plt.xlabel('Time')
plt.ylabel('Eccentricity')
plt.grid()
plt.show()

fig3 = plt.figure()
plt.plot(xp, yp, 'b')
plt.plot(xq, yq, 'r')

plt.xlabel('X')
plt.ylabel('Y')
plt.grid()
plt.show()

fig4 = plt.figure()
plt.plot(xp, zp, 'b')
plt.plot(xq, zq, 'r')

plt.xlabel('X')
plt.ylabel('Z')
plt.grid()
plt.show()

fig5 = plt.figure()
plt.plot(yp, zp, 'b')
plt.plot(yq, zq, 'r')

plt.xlabel('Y')
plt.ylabel('Z')
plt.grid()
plt.show()

figx = plt.figure()
plt.plot(idp, edp, 'b')
plt.plot(idq, edq, 'r')

plt.xlabel('1/cosi')
plt.ylabel('/sqrt (1-e^2)')
plt.grid()
plt.show()

'''
