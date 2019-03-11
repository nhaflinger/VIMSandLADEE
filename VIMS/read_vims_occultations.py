#
# read_vims_occultations.py
#
# author: Doug Creel
#         University of Idaho, Physics Dept.
# created: October 2017
# description: process Cassini VIMS ring occultation data
#

import os
import argparse
import sys
from numpy import sqrt, pi, exp, linspace, random
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import Circle, PathPatch, Wedge, FancyArrowPatch
import mpl_toolkits.mplot3d.art3d as art3d
import spiceypy as spice
from matplotlib.backends.backend_pdf import PdfPages
from astropy import units as u
import astropy.coordinates as coord

sys.path.extend(['/Users/nhaflinger/research/libs/pds-tools-master'])
import pdsparser as pds


class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


# orthographic projection
def orthographic_proj(zfront, zback):
    a = (zfront+zback)/(zfront-zback)
    b = -2*(zfront*zback)/(zfront-zback)

    return np.array([[1,0,0,0],
                        [0,1,0,0],
                        [0,0,a,b],
                        [0,0,-0.0001,zback]])


# convert Julian date to ephemeris time
def cspice_utc2et(kernelfile, ctime):
    spice_ver = spice.tkvrsn('TOOLKIT')
    spice.furnsh(kernelfile)
    kernels_loaded = spice.ktotal("ALL")
    #print(kernels_loaded)
    ntime = spice.utc2et(ctime.value)
    #print(spice_ver)
    return ntime


# read PDS label file
def read_label(labelfile):
    # general data
    global stime, ttime, rstime, rttime, starname, occultdir, profiledir, wavemax, wavemin, ringradiusmax, ringradiusmin, ringlongmax, ringlongmin, ringelev, ligthinc, ringazmax, ringazmin, opaclowest, opachighest, radres

    result = pds.PdsLabel.from_file(labelfile)
    #print(result["START_TIME"])
    stime = result["START_TIME"]
    ttime = result["STOP_TIME"]
    rstime = result["RING_EVENT_START_TIME"]
    rttime = result["RING_EVENT_STOP_TIME"]
    starname = result["STAR_NAME"]
    occultdir = result["RING_OCCULTATION_DIRECTION"]
    #profiledir = result["RING_PROFILE_DIRECTION"]
    wavemax = result["MAXIMUM_WAVELENGTH"]
    wavemin = result["MINIMUM_WAVELENGTH"]
    radres = result["RADIAL_RESOLUTION"]
    ringradiusmax = result["MAXIMUM_RING_RADIUS"]
    ringradiusmin = result["MINIMUM_RING_RADIUS"]
    ringlongmax = result["MAXIMUM_RING_LONGITUDE"]
    ringlongmin = result["MINIMUM_RING_LONGITUDE"]
    ringelev = result["OBSERVED_RING_ELEVATION"]
    lightinc = result["LIGHT_SOURCE_INCIDENCE_ANGLE"]
    ringazmax = result["MAXIMUM_OBSERVED_RING_AZIMUTH"]
    ringazmin = result["MINIMUM_OBSERVED_RING_AZIMUTH"]
    opachighest = result["HIGHEST_DETECTABLE_OPACITY"]
    opaclowest = result["LOWEST_DETECTABLE_OPACITY"]

    # table data (column names)
    global numcol, numrow, minsamp, maxsamp

    subnode = result["SERIES"]
    numcol = subnode["COLUMNS"]
    numrow = subnode["ROWS"]
    minsamp = subnode["MINIMUM_SAMPLING_PARAMETER"]
    maxsamp = subnode["MAXIMUM_SAMPLING_PARAMETER"]


# read PDS table file
def read_table(tablefile):
    # table data (column names)
    global ringradius, ringlong, ringaz, normoptdepth, maxoptdepth, minoptdepth, meansignal, meansigunc, obsevent, ringevent, backgroundmodel, unoccultstar, numsamp, noteflag

    ringradius = np.zeros(numrow.value)
    ringlong = np.zeros(numrow.value)
    ringaz = np.zeros(numrow.value)
    normoptdepth = np.zeros(numrow.value)
    maxoptdepth = np.zeros(numrow.value)
    minoptdepth = np.zeros(numrow.value)
    meansignal = np.zeros(numrow.value)
    meansigunc = np.zeros(numrow.value)
    obsevent = np.zeros(numrow.value)
    ringevent = np.zeros(numrow.value)
    backgroundmodel = np.zeros(numrow.value)
    unoccultstar = np.zeros(numrow.value)
    numsamp = np.zeros(numrow.value)
    noteflag = np.zeros(numrow.value)

    filein = open(tablefile, "r")

    nr = 0
    for line in filein:
        newline = line.strip()
        newline = newline.replace(" ", "")
        rowdata = newline.split(",")

        if (len(rowdata[0]) == 0):
            break

        ringradius[nr] = float(rowdata[0])
        ringlong[nr] = float(rowdata[1])
        ringaz[nr] = float(rowdata[2])
        normoptdepth[nr] = float(rowdata[3])
        maxoptdepth[nr] = float(rowdata[4])
        minoptdepth[nr] = float(rowdata[5])
        if (rowdata[6] != "*******"):
            meansignal[nr] = float(rowdata[6])
        meansigunc[nr] = float(rowdata[7])
        obsevent[nr] = float(rowdata[8])
        ringevent[nr] = float(rowdata[9])
        backgroundmodel[nr] = float(rowdata[10])
        unoccultstar[nr] = float(rowdata[11])
        numsamp[nr] = float(rowdata[12])
        noteflag[nr] = float(rowdata[13])

        nr += 1

    filein.close()


# read star database file
def read_stardb(filename, star):
    global star_ra, star_dec

    filein = open(filename, "r")

    for line in filein:
        newline = line.strip()
        rowdata = newline.split(" ")
        rowdata[0] = rowdata[0].upper()
        if (rowdata[0] == str(star)):
            ra_str = rowdata[1] + ":" + rowdata[2] + ":" + rowdata[3]
            dec_str = rowdata[4] + ":" + rowdata[5] + ":" + rowdata[6]
            ra = coord.Angle(ra_str, unit=u.hour)
            dec = coord.Angle(dec_str, unit=u.degree)
            print(rowdata[0] + " ra and dec: ", ra, dec)
            star_ra, star_dec = ra.degree, dec.degree
            break

    filein.close()


# setup saturn plot
def setup_saturn_plot(ax3, elev, azim, drawz):
    ax3.set_aspect('equal','box')
    ax3.view_init(elev=elev, azim=azim)
    proj3d.persp_transformation = orthographic_proj

    # hide grid and background
    ax3.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax3.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax3.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax3.grid(False)

    # hide z axis in orthographic top view
    if (drawz == False):
        ax3.w_zaxis.line.set_lw(0.)
        ax3.set_zticks([])
        ax3.set_xlabel('X (1000 km)')
        ax3.set_ylabel('Y (1000 km)')

    if (drawz == True):
        ax3.set_zlabel('Z (1000 km)')
        ax3.set_xlabel('X (1000 km)')
        ax3.set_ylabel('Y (1000 km)')
        #ax3.w_zaxis.line.set_lw(0.)
        #ax3.set_zticks([])
        #ax3.w_xaxis.line.set_lw(0.)
        #ax3.set_xticks([])
        #ax3.w_yaxis.line.set_lw(0.)
        #ax3.set_yticks([])



# draw Saturn rings
def draw_rings(ax3, elev, azim, draw_mode):
    # Saturn dimensions
    radius = 60268. / 1000.

    # Saturn rings
    dringmin = 1.110 * radius 
    dringmax = 1.236 * radius 
    cringmin = 1.239 * radius 
    titanringlet = 1.292 * radius 
    maxwellgap = 1.452 * radius 
    cringmax = 1.526 * radius 
    bringmin = 1.526 * radius 
    bringmax = 1.950 * radius 
    aringmin = 2.030 * radius 
    enckegap = 2.214 * radius 
    keelergap = 2.265 * radius 
    aringmax = 2.270 * radius 
    fringmin = 2.320 * radius 
    gringmin = 2.754 * radius 
    gringmax = 2.874 * radius 
    eringmin = 2.987 * radius 
    eringmax = 7.964 * radius 

    if (draw_mode == 'back'):
        offset = -azim*np.pi/180. - 0.5*np.pi
    if (draw_mode == 'front'):
        offset = -azim*np.pi/180. + 0.5*np.pi

    rad, theta = np.mgrid[dringmin:dringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    z = 0. * rad
    line1 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    rad, theta = np.mgrid[cringmin:cringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    z = 0. * rad
    line2 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    rad, theta = np.mgrid[bringmin:bringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    z = 0. * rad
    line3 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    rad, theta = np.mgrid[aringmin:aringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    z = 0. * rad
    line4 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    #rad, theta = np.mgrid[gringmin:gringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    #x = rad * np.cos(theta)
    #y = rad * np.sin(theta)
    #z = 0. * rad
    #line5 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    #rad, theta = np.mgrid[eringmin:eringmax:4j, 0.0-offset:1.0*np.pi-offset:100j]
    #x = rad * np.cos(theta)
    #y = rad * np.sin(theta)
    #z = 0. * rad
    #line6 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.25)

    rad, theta = np.mgrid[fringmin:1.005*fringmin:2j, 0.0-offset:1.0*np.pi-offset:100j]
    x = rad * np.cos(theta)
    y = rad * np.sin(theta)
    z = 0. * rad
    line7 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=25, shade=False, lw=0.1)


# draw Saturn (no rings)
def draw_saturn(ax3, elev, azim):
    # Saturn dimensions
    radius = 60268. / 1000.
    radius_pole = 54364. / 1000.

    # draw Saturn
    phi, theta = np.mgrid[0.0:np.pi:100j, 0.0:2.0*np.pi:100j]
    x = radius*np.sin(phi)*np.cos(theta)
    y = radius*np.sin(phi)*np.sin(theta)
    z = radius_pole*np.cos(phi)

    line3 = ax3.plot_surface(x, y, z, color="w", edgecolor='b', rstride = 8, cstride=5, shade=False, lw=0.25)
    #line3 = ax3.plot_wireframe(x, y, z, color="w", edgecolor='b', rstride = 5, cstride=5, lw=0.25)
    
    ax3.tick_params(labelsize=8)


# draw vector for sun direction from Saturn barycenter 
def draw_sundir(ax4, sundir, start, length, updir):
    a = Arrow3D([start*sundir[0],(start+length)*sundir[0]], [start*sundir[1],(start+length)*sundir[1]],[start*sundir[2],(start+length)*sundir[2]], mutation_scale=10, lw=1, arrowstyle="-|>", color="k", zorder=101)
    ax4.add_artist(a)
    ax4.text((1.00*start+length)*sundir[0], (1.00*start+length)*sundir[1], (1.00*start+length)*sundir[2], "Sun Direction", size=8, zdir=updir, zorder=102)


# draw terminator for shadowed side of planet 
def draw_terminator(ax4, radius, radius_pole, sundir, updir):
    ydir = np.cross(sundir, updir)
    zdir = np.cross(ydir, sundir)

    dotp = np.dot(sundir, updir)
    ang = np.arccos(dotp)

    theta = np.linspace(0., 2.*np.pi, 100)
    phi = np.linspace(0., np.pi, 100)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = 0. * np.sin(phi)

    ax4.plot(x, y, z, color='r', linewidth=1.0, zorder=102)


# draw path of occultation
def draw_occult_track(ax3):
    otrackx = np.zeros(len(ringradius))
    otracky = np.zeros(len(ringradius))
    otrackz = np.zeros(len(ringradius))
    for i in range(len(ringradius)):
        # when angle "theta" is zero this gives same result as "ringlong"
        #otrackx[i] = ringradius[i] * np.cos((ringaz[i])*np.pi/180.) / 1000.;
        #otracky[i] = ringradius[i] * np.sin((ringaz[i])*np.pi/180.) / 1000.;
        otrackx[i] = ringradius[i] * np.cos((ringlong[i])*np.pi/180.) / 1000.;
        otracky[i] = ringradius[i] * np.sin((ringlong[i])*np.pi/180.) / 1000.;
        otrackz[i] = 0.;

    line5 = ax3.plot(otrackx, otracky, otrackz, zorder=100, linewidth=1.0, color='r', linestyle='--')


# finish plot setings
def wrapup_saturn_plot(ax3, datafile):
    if (float(ringelev) < 0.):
        yloc = 0.1
    else:
        yloc = 0.9

    ax3.set_title('Occultation Track for ' + datafile, fontsize=8, y=yloc)

    ax3.set_xlim([-200, 200])
    ax3.set_ylim([-200, 200])
    ax3.set_zlim([-200, 200])

    scaling = np.array([getattr(ax3, 'get_{}lim'.format(dim))() for dim in 'xyz']); ax3.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3) 


# plot occultation data
def plot_results(datafile, plot_names):
    colors = ['grey','magenta','violet','blue','cyan','green','orange','red','maroon','brown']

    fig1, (ax1, ax2) = plt.subplots(2, figsize=(9,9), sharex=True, sharey=False)

    fig1.suptitle(datafile, fontsize=12)

    plot_names.append("norm_optical_depth")
    #xr = [105000, 145500]
    line1 = ax1.plot(ringradius/1000., meansignal, lw=0.5)
    #ax1.fill_between(ringradius, meansignal-meansigunc, meansignal+meansigunc, alpha=0.5)
    line2 = ax2.plot(ringradius/1000., normoptdepth, lw=0.5)
    #ax2.invert_yaxis()
    ax1.grid(linestyle='--', linewidth=1)
    ax2.grid(linestyle='--', linewidth=1)
    #ax2.fill_between(ringradius, minoptdepth, maxoptdepth, alpha=0.5)
    #ax1.set_xlim(xr)
    #ax1.set_xlabel('Ring Radius (km)', fontsize=8)
    ax1.set_ylabel('Mean Signal', fontsize=8)
    ax2.set_xlabel(r'Ring Radius, $\rho$ (1000 km)', fontsize=8)
    ax2.set_ylabel(r'Normalized Optical Depth, $\tau$', fontsize=8)

    #ax1.yaxis.set_major_locator(ticker.MultipleLocator(200))
    #ax1.yaxis.set_minor_locator(ticker.MultipleLocator(40))
    #ax2.yaxis.set_major_locator(ticker.MultipleLocator(0.1))
    #ax2.yaxis.set_minor_locator(ticker.MultipleLocator(0.02))
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(1))

    # get ephemeris times during occultation
    step = int(numrow)
    stimes = [x*(ett-ets)/step + ets for x in range(step)]
    rtimes = [x*(rtt-rts)/step + rts for x in range(step)]
    rstimes = []
    for i in range(len(stimes)):
        rstimes.append(stimes[i] - rtimes[i])

    #Run spkpos as a vectorized function
    #scstate, lightTimes = spice.spkpos('Cassini', stimes, 'J2000', 'NONE', 'SATURN BARYCENTER')

    # determine direction to sun from Saturn barycenter
    sunstate, sLightTimes = spice.spkpos('SUN', stimes, 'J2000', 'NONE', 'SATURN BARYCENTER')
    sidx = 0
    sundirmag = sqrt(sunstate[sidx][0]*sunstate[sidx][0] + sunstate[sidx][1]*sunstate[sidx][1] + sunstate[sidx][2]*sunstate[sidx][2])
    sundir = np.zeros(3)
    sundir[0] = sunstate[sidx][0] / sundirmag
    sundir[1] = sunstate[sidx][1] / sundirmag
    sundir[2] = sunstate[sidx][2] / sundirmag

    scstate = np.zeros(6*len(stimes)).reshape(len(stimes), 6)
    lightTimes = np.zeros(6*len(stimes)).reshape(len(stimes), 6)

    # set spice error handling
    spice.erract("SET", 100, action='report')

    for i in range(len(stimes)):
        scstate[i], lightTimes[i] = spice.spkezr('Cassini', stimes[i], 'J2000', 'NONE', 'SATURN BARYCENTER')
    
    # Clean up the kernels
    spice.kclear()

    # determine if ingress or egress pass
    file_parts = datafile.split("_")
    if (file_parts[4] == "I"):
        viewpos = scstate[step-1]
    elif (file_parts[4] == "E"):
        viewpos = scstate[0]

    #RA and dec of Saturn N pole from Jacobson 2008 SAT 287
    # compute ra and dec of Saturn
    et_center = 0.5*(ets + ett)
    global polera, poledec
    if (file_parts[4] == "I"):
        et_center = ett
    elif (file_parts[4] == "E"):
        et_center = ets

    polera = 40.58238-.0275*et_center/24./3600./365.25/100.
    poledec = 83.53762-.0048*et_center/24./3600./365.25/100.

    # determine view orientation from Cassini
    #viewpos = spice.rotvec(viewpos[0:3], (90. - (poledec+star_dec)), 1)
    #viewpos = spice.rotvec(viewpos[0:3], (180. - (polera+star_ra)), 3)

    # determine view orientation from star
    rad = 1.0*np.cos(star_dec * np.pi/180.)
    rho = np.zeros(3)
    rho[0] = rad*np.cos(star_ra * np.pi/180.)
    rho[1] = rad*np.sin(star_ra * np.pi/180.)
    rho[2] = 1.0*np.sin(star_dec * np.pi/180.)
    capN = (polera + 90.) * np.pi/180.
    capJ = (90. - poledec) * np.pi/180.

    rot = np.zeros(3*3).reshape(3,3)
    rot[0,0] = np.cos(capN)
    rot[1,0] = -np.sin(capN)*np.cos(capJ)
    rot[2,0] = np.sin(capN)*np.sin(capJ)
    rot[0,1] = np.sin(capN)
    rot[1,1] = np.cos(capN)*np.cos(capJ)
    rot[2,1] = -np.cos(capN)*np.sin(capJ)
    rot[0,2] = 0.0
    rot[1,2] = np.sin(capJ)
    rot[2,2] = np.cos(capJ)

    #rotT = rot.transpose()
    rho = rot.dot(rho)

    sundir = rot.dot(sundir)

    updir = (0,0,1)
    updir = rot.dot(updir)

    theta = np.arctan2(rho[1], rho[0])*180./np.pi + 00.
    phi = np.arctan2(rho[2], sqrt(rho[0]*rho[0] + rho[1]*rho[1]))*180./np.pi
    #theta = np.arctan2(viewpos[1], viewpos[0])*180./np.pi
    #phi = np.arctan2(viewpos[2], sqrt(viewpos[0]*viewpos[0] + viewpos[1]*viewpos[1]))*180./np.pi
    print(theta, phi)

    radius = 60268. / 1000.
    radius_pole = 54364. / 1000.
    fringmin = 2.320 * radius 

    # view from star
    plot_names.append("occultation_track_" + starname)
    fig2 = plt.figure(figsize=(9,9))
    ax3 = fig2.add_subplot(111, projection='3d')

    setup_saturn_plot(ax3, phi, theta, True)
    draw_rings(ax3, phi, theta, 'back')
    draw_saturn(ax3, phi, theta)
    draw_rings(ax3, phi, theta, 'front')
    draw_occult_track(ax3)
    zdir = (sundir[0], sundir[1], sundir[2])
    draw_sundir(ax3, sundir, fringmin+20, 100, zdir)
    draw_terminator(ax3, radius, radius_pole, sundir, updir)
    wrapup_saturn_plot(ax3, datafile)

    # view from above
    plot_names.append("occultation_track_top")
    fig3 = plt.figure(figsize=(9,9))
    ax4 = fig3.add_subplot(111, projection='3d')

    setup_saturn_plot(ax4, 90, 90, False)
    draw_rings(ax4, 90, 90, 'back')
    draw_saturn(ax4, 90, 90)
    draw_rings(ax4, 90, 90, 'front')
    draw_occult_track(ax4)
    zdir = (0,0,1)
    draw_sundir(ax4, sundir, fringmin+20, 100, zdir)
    draw_terminator(ax4, radius, radius_pole, sundir, updir)
    wrapup_saturn_plot(ax4, datafile)

    # additional plot data
    plot_names.append("occultation_data")
    tdiff = obsevent - ringevent
    fig4, (ax5, ax6, ax7, ax8) = plt.subplots(4, figsize=(9,9), sharex=True, sharey=False)
    line5 = ax5.plot(obsevent, tdiff, lw=1.0)
    line6 = ax6.plot(obsevent, ringradius/1000., lw=1.0)
    line7 = ax7.plot(obsevent, ringlong, lw=1.0)
    line8 = ax8.plot(obsevent, ringaz, lw=1.0)

    ax8.set_xlabel(r'Observed Event Time, $t_{OET}$ (sec)', fontsize=10)
    ax5.set_ylabel(r'$t_{OET}$ - $t_{RET}$ (sec)', fontsize=10)
    ax6.set_ylabel(r'$\rho$ (1000 km)', fontsize=10)
    ax7.set_ylabel(r'$\phi_{RL}$ (deg)', fontsize=10)
    ax8.set_ylabel(r'$\phi_{ORA}$ (deg)', fontsize=10)

    ax8.xaxis.set_major_locator(ticker.MultipleLocator(1000))
    ax8.xaxis.set_minor_locator(ticker.MultipleLocator(200))

    ax5.xaxis.set_tick_params(which='both', direction='in')
    ax6.xaxis.set_tick_params(which='both', direction='in')
    ax7.xaxis.set_tick_params(which='both', direction='in')
    ax8.xaxis.set_tick_params(which='both', direction='in')

    ax5.yaxis.set_label_coords(-0.075, 0.5)
    ax6.yaxis.set_label_coords(-0.075, 0.5)
    ax7.yaxis.set_label_coords(-0.075, 0.5)
    ax8.yaxis.set_label_coords(-0.075, 0.5)

    fig4.subplots_adjust(hspace=0.00)

    plot_names.append("spacecraft_statevec")
    fig5, (ax10, ax11, ax12, ax13, ax14, ax15) = plt.subplots(6, figsize=(9,9), sharex=True, sharey=False)
    line10 = ax10.plot(obsevent, scstate.T[0]/1000., lw=1.0)
    line11 = ax11.plot(obsevent, scstate.T[1]/1000., lw=1.0)
    line12 = ax12.plot(obsevent, scstate.T[2]/1000., lw=1.0)
    line13 = ax13.plot(obsevent, scstate.T[3], lw=1.0)
    line14 = ax14.plot(obsevent, scstate.T[4], lw=1.0)
    line15 = ax15.plot(obsevent, scstate.T[5], lw=1.0)

    ax15.set_xlabel(r'Observed Event Time, $t_{OET}$ (sec)', fontsize=10)
    ax10.set_ylabel(r'$r_x$ (1000 km)', fontsize=10)
    ax11.set_ylabel(r'$r_y$ (1000 km)', fontsize=10)
    ax12.set_ylabel(r'$r_z$ (1000 km)', fontsize=10)
    ax13.set_ylabel(r'$v_x$ (km/s)', fontsize=10)
    ax14.set_ylabel(r'$v_y$ (km/s)', fontsize=10)
    ax15.set_ylabel(r'$v_z$ (km/s)', fontsize=10)

    ax10.yaxis.set_label_coords(-0.085, 0.5)
    ax11.yaxis.set_label_coords(-0.085, 0.5)
    ax12.yaxis.set_label_coords(-0.085, 0.5)
    ax13.yaxis.set_label_coords(-0.085, 0.5)
    ax14.yaxis.set_label_coords(-0.085, 0.5)
    ax15.yaxis.set_label_coords(-0.085, 0.5)

    ax12.xaxis.set_major_locator(ticker.MultipleLocator(1000))
    ax12.xaxis.set_minor_locator(ticker.MultipleLocator(200))

    ax10.xaxis.set_tick_params(which='both', direction='in')
    ax11.xaxis.set_tick_params(which='both', direction='in')
    ax12.xaxis.set_tick_params(which='both', direction='in')
    ax13.xaxis.set_tick_params(which='both', direction='in')
    ax14.xaxis.set_tick_params(which='both', direction='in')
    ax15.xaxis.set_tick_params(which='both', direction='in')

    fig5.subplots_adjust(hspace=0.00)


def main():
    global ets, ett, rts, rtt

    # define command line arguments
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("-f", "--filelist", help="text file with list of files")
    parser.add_argument("-s", "--stardb", help="star database")
    parser.add_argument("-d", "--datadir", help="data directory)")
    parser.add_argument("-r", "--results", help="destination directory")
    parser.add_argument("-k", "--kernels", help="kernel meta file")
    parser.add_argument("-p", "--saveplot", help="save plot file", action="store_true")
    args = parser.parse_args()

    # get command line arguments
    filename = args.filelist
    datadir = args.datadir
    stardb = args.stardb
    results = args.results
    kernelfile = args.kernels
    saveplot = args.saveplot

    # open file with list of label files
    filein = open(filename, "r")

    occ = []
    for line in filein:
        occ.append(line.strip())

    nocc = len(occ)
    for i in range(nocc):
        labelfile = datadir + "/" + occ[i] + ".LBL"
        tablefile = datadir + "/" + occ[i] + ".TAB"
        if (os.path.exists(labelfile)):
            print("Processing: " + labelfile)
            read_label(labelfile)
            
            # get times for occultation
            ets = cspice_utc2et(kernelfile, stime)
            ett = cspice_utc2et(kernelfile, ttime)
            rts = cspice_utc2et(kernelfile, rstime)
            rtt = cspice_utc2et(kernelfile, rttime)

            if (os.path.exists(tablefile)):
                read_table(tablefile)
                read_stardb(stardb, starname)
                plot_names = []
                plot_results(occ[i], plot_names)
    
                if (saveplot):
                    print("Plotting results.")
                    newdir = results + "/" + occ[i] 
                    if (os.path.exists(newdir) == False):
                        os.mkdir(newdir)
                
                    for i in plt.get_fignums():
                        plotfile = newdir + "/" + plot_names[i-1] + ".pdf"
                        plt.figure(i)
                        plt.savefig(plotfile)
                        plt.close()

                else:
                    plt.show()

                    #with PdfPages(plotfile) as pp:
                        #for i in plt.get_fignums():
                            #plt.figure(i)
                            #pp.savefig()
                            #plt.close()
        else :
            print(labelfile + " does not exist!")

    print("All tasks completed!")
    filein.close()



if __name__ == "__main__":
    main()


