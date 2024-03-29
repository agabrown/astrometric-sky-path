"""
This module provides functions for calculating the astrometric paths of stars on the sky as seen by an
observer in orbit around the solar system barycentre. That is, the topocentric coordinate directions as a
function of time are calculated.

Anthony Brown Nov 2018 - Nov 2018
"""

import numpy as np

from astropy import constants
from astropy import units
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
from erfa import pmpx, c2s, epj2jd, pmsafe, s2c

_radtomas = (180*3600*1000)/np.pi
_mastorad = np.pi/(180*3600*1000)
_kmps_to_aupyr = (units.year.to(units.s)*units.km.to(units.m))/constants.au.value


def ephemeris_earth_astropy(t):
    """
    Calculate the ephemeris for "earth" in the BCRS using astropy tools.    
    
    NOTE: There are several versions of the solar system ephemeris available in astropy and others can be
    provided through a URL (see the documentation of astropy.coordinates.solar_system_ephemeris).
    Depending on the accuracy needed, employing  an ephemeris different from the default may be better.

    Parameters
    ----------
    
    t : float array
        Times at which to calculate the ephemeris in Julian years TCB.
        
    Returns
    -------
    
    Array of shape (3,t.size) representing the xyz components of the ephemeris at times t.
    """
    times = Time(t, format='jyear', scale='tcb')
    earthEph = get_body_barycentric('earth', times)
    return np.vstack((earthEph.x.value, earthEph.y.value, earthEph.z.value))


def epoch_topocentric_coordinates(alpha, delta, parallax, mura, mudec, vrad, t, refepoch, ephem):
    """
    For each observation epoch calculate the topocentric coordinate directions (alpha(t), delta(t)) given
    the astrometric parameters of a source, the observation times, and the ephemeris (in the BCRS) for
    the observer. Also calculate the local plane coordinates xi(t) and eta(t).

    The code uses the Pyhton ERFA routines (https://pyerfa.readthedocs.io/)

    Parameters
    ----------
    
    alpha : float
        Right ascension at reference epoch (radians)
    delta : float
        Declination at reference epoch (radians)
    parallax : float
        Parallax (mas), negative values allowed
    mura : float
        Proper motion in right ascension, including cos(delta) factor (mas/yr)
    mudec : float
        Proper motion in declination (mas/yr)
    vrad : float
        Radial velocity (km/s)
    t : float array
        Observation times (Julian year TCB)
    refepoch : float
        Reference epoch (Julian year TCB)
    ephem : function
        Funtion providing the observer's ephemeris in BCRS at times t (units of AU)
                
    Returns
    -------
    
    alpha, delta, xi, eta : float arrays
        The coordinates directions to the sources for each time t and the corresponsing local plane coordinates (xi,
        eta), referred to the source direction at the reference epoch.
        Units are radians for (alpha, delta) and mas for (xi, eta).
    """
    # Normal triad, defined at the reference epoch.
    p = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
    q = np.array([-np.sin(delta)*np.cos(alpha), -np.sin(delta)*np.sin(alpha), np.cos(delta)])
    r = np.array([np.cos(delta)*np.cos(alpha), np.cos(delta)*np.sin(alpha), np.sin(delta)])

    # Calculate observer's ephemeris.
    bO_bcrs = ephem(t)

    # Unit conversions
    plx = parallax/1000.0
    pmra = mura*_mastorad/np.cos(delta)
    pmdec = mudec*_mastorad

    uO = pmpx(alpha, delta, pmra, pmdec, plx, vrad, t-refepoch, bO_bcrs.T)

    # Local plane coordinates which approximately equal delta_alpha*cos(delta) and delta_delta
    xi = np.dot(p,uO.T)/np.dot(r,uO.T)*_radtomas
    eta = np.dot(q,uO.T)/np.dot(r,uO.T)*_radtomas

    alpha_obs, delta_obs = c2s(uO)

    return alpha_obs, delta_obs, xi, eta


def epoch_barycentric_coordinates(alpha, delta, parallax, mura, mudec, vrad, t, refepoch):
    """
    For each observation epoch calculate the barycentric coordinate directions (alpha(t), delta(t)) given
    the astrometric parameters of a source, the observation times. Also calculate the local plane
    coordinates xi(t) and eta(t). 
    
    HERE THE PARALLACTIC MOTION IS THUS NOT INCLUDED and one obtains a sky path as seen by an observer at
    the solar system barycentre.

    The code uses the Pyhton ERFA routines (https://pyerfa.readthedocs.io/)

    Parameters
    ----------
    
    alpha : float
        Right ascension at reference epoch (radians)
    delta : float
        Declination at reference epoch (radians)
    parallax : float
        Parallax (mas), negative values allowed
    mura : float
        Proper motion in right ascension, including cos(delta) factor (mas/yr)
    mudec : float
        Proper motion in declination (mas/yr)
    vrad : float
        Radial velocity (km/s)
    t : float array
        Observation times (Julian year TCB)
    refepoch : float
        Reference epoch (Julian year TCB)
                
    Returns
    -------
    
    alpha, delta, xi, eta : float arrays
        The coordinates directions to the sources for each time t and the corresponsing local plane coordinates (xi,
        eta), referred to the source direction at the reference epoch.
        Units are radians for (alpha, delta) and mas for (xi, eta).
    """
    # Normal triad, defined at the reference epoch.
    p = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
    q = np.array([-np.sin(delta)*np.cos(alpha), -np.sin(delta)*np.sin(alpha), np.cos(delta)])
    r = np.array([np.cos(delta)*np.cos(alpha), np.cos(delta)*np.sin(alpha), np.sin(delta)])

    # Unit conversions
    plx = parallax/1000.0
    pmra = mura*_mastorad/np.cos(delta)
    pmdec = mudec*_mastorad
    zpjdref, mjdred = epj2jd(refepoch)
    zpjd, mjd = epj2jd(t)

    alpha_obs, delta_obs, _, _, _, _ = pmsafe(alpha, delta, pmra, pmdec, plx, vrad, zpjdref, mjdred, zpjd, mjd)
    uO = s2c(alpha_obs, delta_obs)

    # Local plane coordinates which approximately equal delta_alpha*cos(delta) and delta_delta
    xi = np.dot(p,uO.T)/np.dot(r,uO.T)*_radtomas
    eta = np.dot(q,uO.T)/np.dot(r,uO.T)*_radtomas

    return alpha_obs, delta_obs, xi, eta

def epoch_topocentric_coordinates_noerfa(alpha, delta, parallax, mura, mudec, vrad, t, refepoch, ephem):
    """
    For each observation epoch calculate the topocentric coordinate directions (alpha(t), delta(t)) given
    the astrometric parameters of a source, the observation times, and the ephemeris (in the BCRS) for
    the observer. Also calculate the local plane coordinates xi(t) and eta(t).

    The code is partly based on the SOFA library (http://www.iausofa.org/) pmpx.c code.

    Parameters
    ----------
    
    alpha : float
        Right ascension at reference epoch (radians)
    delta : float
        Declination at reference epoch (radians)
    parallax : float
        Parallax (mas), negative values allowed
    mura : float
        Proper motion in right ascension, including cos(delta) factor (mas/yr)
    mudec : float
        Proper motion in declination (mas/yr)
    vrad : float
        Radial velocity (km/s)
    t : float array
        Observation times (Julian year TCB)
    refepoch : float
        Reference epoch (Julian year TCB)
    ephem : function
        Funtion providing the observer's ephemeris in BCRS at times t (units of AU)
                
    Returns
    -------
    
    Arrays alpha, delta, xi, eta. Units are radians for (alpha, delta) and mas for (xi, eta).
    """
    
    # Normal triad, defined at the reference epoch.
    p = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
    q = np.array([-np.sin(delta)*np.cos(alpha), -np.sin(delta)*np.sin(alpha), np.cos(delta)])
    r = np.array([np.cos(delta)*np.cos(alpha), np.cos(delta)*np.sin(alpha), np.sin(delta)])
   
    # Calculate observer's ephemeris.
    bO_bcrs = ephem(t)

    # Calculate the Roemer delay, take units into account.
    tB = t + np.dot(r, bO_bcrs)*constants.au.value/constants.c.value/units.year.to(units.s)
    
    plxrad = parallax*_mastorad
    murarad = mura*_mastorad
    mudecrad = mudec*_mastorad
    mur = vrad*_kmps_to_aupyr*np.abs(plxrad)

    uO = np.repeat(r, t.size).reshape((r.size, t.size))
    uO = uO + np.tensordot((p*murarad + q*mudecrad + r*mur),(tB-refepoch),axes=0) - plxrad*bO_bcrs
    
    # Local plane coordinates which approximately equal delta_alpha*cos(delta) and delta_delta
    xi = np.dot(p,uO)/np.dot(r,uO)*_radtomas
    eta = np.dot(q,uO)/np.dot(r,uO)*_radtomas

    alpha_obs = np.arctan2(uO[1,:], uO[0,:])
    delta_obs = np.arctan2(uO[2,:],np.sqrt(uO[0,:]**2+uO[1,:]**2))
                 
    return alpha_obs, delta_obs, xi, eta


def epoch_barycentric_coordinates_noerfa(alpha, delta, parallax, mura, mudec, vrad, t, refepoch):
    """
    For each observation epoch calculate the barycentric coordinate directions (alpha(t), delta(t)) given
    the astrometric parameters of a source, the observation times. Also calculate the local plane
    coordinates xi(t) and eta(t). 
    
    HERE THE PARALLACTIC MOTION IS THUS NOT INCLUDED and one obtains a sky path as seen by an observer at
    the solar system barycentre.

    The code is partly based on the SOFA library (http://www.iausofa.org/) pmpx.c code.

    Parameters
    ----------
    
    alpha : float
        Right ascension at reference epoch (radians)
    delta : float
        Declination at reference epoch (radians)
    parallax : float
        Parallax (mas), negative values allowed
    mura : float
        Proper motion in right ascension, including cos(delta) factor (mas/yr)
    mudec : float
        Proper motion in declination (mas/yr)
    vrad : float
        Radial velocity (km/s)
    t : float array
        Observation times (Julian year TCB)
    refepoch : float
        Reference epoch (Julian year TCB)
                
    Returns
    -------
    
    Arrays alpha, delta, xi, eta. Units are radians for (alpha, delta) and mas for (xi, eta).
    """
    
    # Normal triad, defined at the reference epoch.
    p = np.array([-np.sin(alpha), np.cos(alpha), 0.0])
    q = np.array([-np.sin(delta)*np.cos(alpha), -np.sin(delta)*np.sin(alpha), np.cos(delta)])
    r = np.array([np.cos(delta)*np.cos(alpha), np.cos(delta)*np.sin(alpha), np.sin(delta)])
   
    plxrad = parallax*_mastorad
    murarad = mura*_mastorad
    mudecrad = mudec*_mastorad
    mur = vrad*_kmps_to_aupyr*np.abs(plxrad)

    uO = np.repeat(r, t.size).reshape((r.size, t.size))
    uO = uO + np.tensordot((p*murarad + q*mudecrad + r*mur),(t-refepoch),axes=0)
    
    # Local plane coordinates which approximately equal delta_alpha*cos(delta) and delta_delta
    xi = np.dot(p,uO)/np.dot(r,uO)*_radtomas
    eta = np.dot(q,uO)/np.dot(r,uO)*_radtomas

    alpha_obs = np.arctan2(uO[1,:], uO[0,:])
    delta_obs = np.arctan2(uO[2,:],np.sqrt(uO[0,:]**2+uO[1,:]**2))
                 
    return alpha_obs, delta_obs, xi, eta
