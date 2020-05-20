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
mp = mEarth
# mass of planet q
mq = mJupiter

# Period of planet P
Tp = 1.

# Period of planet Q
Tq = 1.

# Current coordination of the planet p
xp1 = rEarth * 0.8
yp1 = rEarth * 0.1
# xp1 = rEarth
# yp1 = 0.
zp1 = np.sqrt(rEarth ** 2 - xp1 ** 2 - yp1 ** 2)
rp = rEarth

# Current coordination of the planet q
xq1 = rJupiter * 0.7
yq1 = rJupiter * 0.1
# xq1 = rJupiter
# yq1 = 0.
zq1 = np.sqrt(rJupiter ** 2 - xq1 ** 2 - yq1 ** 2)
rq = rJupiter

# Current velocity of the planet p         vpy = sqrt(GM/4)
vpx = 0.
vpy = np.sqrt((G * ms) / rEarth) * 0.95
vpz = np.sqrt((G * ms) / rEarth - vpy ** 2)
# vpy = np.sqrt((G * ms) / rEarth)
# vpz = 0.
vp = np.sqrt(vpx ** 2 + vpy ** 2 + vpz ** 2)

# Current velocity of the planet q         vpy = sqrt(GM/4)
vqx = 0.
vqy = np.sqrt((G * ms) / rJupiter) * 0.92
vqz = np.sqrt((G * ms) / rJupiter - vqy ** 2)
# vqy = np.sqrt((G * ms) / rJupiter)
# vqz = 0.
vq = np.sqrt(vqx ** 2 + vqy ** 2 + vqz ** 2)

# Specific Angular Momentum of planet P and Q
hp = np.sqrt((vpy * zp1 - vpz * yp1) ** 2 + (vpx * zp1 - vpz * xp1) ** 2 + (vpx * yp1 - vpy * zp1) ** 2)
hq = np.sqrt((vqy * zq1 - vqz * yq1) ** 2 + (vqx * zq1 - vqz * xq1) ** 2 + (vqx * yq1 - vqy * zq1) ** 2)

# hp = np.sqrt((yp1 * vpz - vpy * zp1) ** 2 + (vpx * zp1 - xp1 * vpz) ** 2 + (xp1 * vpy - vpx * yp1) ** 2)
# hq = np.sqrt((yq1 * vqz - vqy * zq1) ** 2 + (vqx * zq1 - xq1 * vqz) ** 2 + (xq1 * vqy - vqx * yq1) ** 2)

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

# Calculate position every tstep
# tstep in day, tsteps in second
tstep = 1.
tsteps = 86400. * tstep
tlimit = math.floor(36500. / tstep)  # floor() 10.234 -> 10

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
    ip.append(np.arctan(zp1 / np.sqrt(xp1 * xp1 + yp1 * yp1)) * 180 / pi)
    iq.append(np.arctan(zq1 / np.sqrt(xq1 * xq1 + yq1 * yq1)) * 180 / pi)

    # Add updated eccentricity to list \sqrt{1+E(2h^2/mGM^2))

    # print(np.sqrt(1 + 2*Eq * (hq ** 2 / (mq * (G * ms) ** 2))))
    # print(np.sqrt(1 + 2*Ep * (hp ** 2 / (mp * (G * ms) ** 2))))

    '''
    print("vp:"+str(vp)+", rp:"+str(rp))
    print("vq:"+str(rp)+", rq:"+str(rq))
    print("hp:"+str(hp))
    print("hq:"+str(hq))
    '''

    ep.append(np.sqrt(1 + 2 * Ep * (hp ** 2 / ((G * ms) ** 2))))
    eq.append(np.sqrt(1 + 2 * Eq * (hq ** 2 / ((G * ms) ** 2))))

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
    aPx = -aPs * xp1 / Dps - aPq * xp1 / Dpq
    aPy = -aPs * yp1 / Dps - aPq * yp1 / Dpq
    aPz = -aPs * zp1 / Dps - aPq * zp1 / Dpq

    # Update acc of planet Q due to S and P
    aQs = G * ms / Dqs * Dqs
    aQp = G * mp / Dpq * Dpq

    # Update acc components of planet Q
    aQx = -aQs * xq1 / Dqs - aQp * xq1 / Dpq
    aQy = -aQs * yq1 / Dqs - aQp * yq1 / Dpq
    aQz = -aQs * zq1 / Dqs - aQp * zq1 / Dpq

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
    hp = np.sqrt((vpy * zp1 - vpz * yp1) ** 2 + (vpx * zp1 - vpz * xp1) ** 2 + (vpx * yp1 - vpy * zp1) ** 2)
    hq = np.sqrt((vqy * zq1 - vqz * yq1) ** 2 + (vqx * zq1 - vqz * xq1) ** 2 + (vqx * yq1 - vqy * zq1) ** 2)

    # Update Time
    t += tsteps

fig = plt.figure()
# ax = fig.add_subplot(111, aspect='equal')

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

# KE + PE = Energy of the orbit --> instantaneous eccentricity
