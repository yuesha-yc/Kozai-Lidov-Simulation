import matplotlib.pyplot as plt
import numpy as np
import math

SimulationIntervalInYears = int(input())
SimulationPlanetQ = input()

def getPlanetQMass(stringinput):
    if stringinput == "jupiter":
        return mJupiter
    else:
        return mSaturn

def getPlanetQRadius(stringinput):
    if stringinput == "jupiter":
        return rJupiter
    else:
        return rSaturn

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
mNeptune = 1.0247e26
mSaturn = 5.6846e26
mMercury = 3.3022e23
# orbital distance of common planets in AU
rJupiter = 5.203
rNeptune = 30.10
rSaturn = 9.582
rEarth = 1.000
rMercury = 0.460

'''Intialization of Object Masses'''
# mass of the star
ms = mSun
# mass of planet p
mp = mEarth
# mass of planet q
mq = getPlanetQMass(SimulationPlanetQ)

'''Intialization of Object Orbit radius'''
# radius of planet p
rp = rEarth
# radius of planet q
rq = getPlanetQRadius(SimulationPlanetQ)
# Period of planet P
Tp = 1.
# Period of planet Q
Tq = 1.

# Initial Inclination to X-Y plane of planet P orbit
eccP = 10 / 180 * pi
# Initial Inclination to X-Y plane planet Q orbit
eccQ = -10 / 180 * pi

# Current coordination of the planet p [10 degree to X-Y plane initially]
xp1 = rp * 0.98475
yp1 = rp * 0.01000
zp1 = rp * 0.17368

# Current coordination of the planet q [10 degree to X-Y plane initially]
xq1 = rq * 0.98475
yq1 = rq * 0.01000
zq1 = rq * -0.17368


# Current velocity of the planet p         vpy = sqrt(GM/4)
vpx = 0.
vpy = np.sqrt((G * ms) / rp) * 0.99
vpz = np.sqrt((G * ms) / rp - vpy ** 2)
vp = np.sqrt(vpx ** 2 + vpy ** 2 + vpz ** 2)

# Current velocity of the planet q         vpy = sqrt(GM/4)
vqx = 0.
vqy = np.sqrt((G * ms) / rq) * 0.98
vqz = np.sqrt((G * ms) / rq - vqy ** 2)
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

# mutal inclination between p and q
ipq = []

# eccentricity of planet P list
ep = []

# eccentricity of planet Q list
eq = []

# Lz
Lz = []
SumLz = 0.

# Calculate position every tstep
# tstep in day, tsteps in second
tstep = 1.
tsteps = 86400. * tstep
tlimit = math.floor(SimulationIntervalInYears*365. / tstep)  # floor() 10.234 -> 10

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
    ipq.append(np.absolute(ip1-iq1) * 180 / pi)

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

LzAverage = SumLz / tlimit


def plotTripile():
    plt.figure()

    ax1 = plt.subplot2grid((2, 2), (0, 0), colspan = 1, rowspan = 1)
    plt.plot(time, ipq, 'b', linewidth=0.5)
    # plt.plot(time, ip, 'r')
    # plt.plot(time, iq, 'r')
    ax1.set_title('Inclination versus Time'+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    plt.xlabel('Time / Day')
    plt.ylabel('Inclination (i-total) / Degree')
    plt.grid()

    ax2 = plt.subplot2grid((2, 2), (0, 1), colspan = 1, rowspan = 1)
    plt.plot(time, ep, 'b', linewidth=0.5)
    plt.plot(time, eq, 'r', linewidth=0.5)
    ax2.set_title('Eccentricities versus Time'+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    plt.xlabel('Time / Day')
    plt.ylabel('Eccentricity')
    plt.grid()

    ax6 = plt.subplot2grid((2, 2), (1, 0), colspan = 2, rowspan = 1)
    plt.plot(time, Lz, 'g', linewidth=0.5)
    ax6.axhline(LzAverage, ls = '--')
    ax6.set_title("Lz versus Time"+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    ax6.annotate('Average Lz: ' + str(round(LzAverage,3)), xy=(SimulationIntervalInYears*365,0.96))
    plt.xlabel('Time / Day')
    plt.ylabel('Lz')
    plt.grid()

    plt.show()


def plotEccentricity():
    plt.figure()
    ax2 = plt.subplot2grid((1, 1), (0, 0), colspan = 1, rowspan = 1)
    plt.plot(time, ep, 'b', linewidth=0.5)
    plt.plot(time, eq, 'r', linewidth=0.5)
    ax2.set_title('Eccentricities of planetary orbits versus Time'+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    plt.xlabel('Time / Day')
    plt.ylabel('Eccentricity')
    plt.grid()
    plt.show()

def plotMutalInc():
    plt.figure()
    ax1 = plt.subplot2grid((1, 1), (0, 0), colspan = 1, rowspan = 1)
    plt.plot(time, ipq, 'g', linewidth=0.5)
    ax1.set_title('Mutual Inclination between planetary orbits versus Time'+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    plt.xlabel('Time / Day')
    plt.ylabel('Inclination (i-total) / Degree')
    plt.grid()
    plt.show()

def plotLz():
    plt.figure()
    ax6 = plt.subplot2grid((1, 1), (0, 0), colspan = 1, rowspan = 1)
    plt.plot(time, Lz, 'g', linewidth=0.5)
    ax6.axhline(LzAverage, ls = '--')
    ax6.set_title("Kozai-Lidov Constant (Lz) versus Time"+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    ax6.annotate('Average Lz: ' + str(round(LzAverage,3)), xy=(SimulationIntervalInYears*365,0.96))
    plt.xlabel('Time / Day')
    plt.ylabel('Lz')
    plt.grid()
    plt.show()

def plotSeparateInclinations():
    plt.figure()
    ax1 = plt.subplot2grid((1, 1), (0, 0), colspan = 1, rowspan = 1)
    plt.plot(time, ip, 'b', linewidth=0.5)
    plt.plot(time, iq, 'r', linewidth=0.5)
    ax1.set_title('Inclinations between planetary orbits and X-Y Plane versus Time'+' - Simulation of '+str(SimulationIntervalInYears)+ ' Years')
    plt.xlabel('Time / Day')
    plt.ylabel('Inclination (planetP=blue, planetQ=red) / Degree')
    plt.grid()
    plt.show()



# plotTripile()
plotEccentricity()
plotMutalInc()
plotLz()
plotSeparateInclinations()