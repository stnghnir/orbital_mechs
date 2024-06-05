""" Analytic formulas of conic sections
"""
from math import radians, degrees, sin, cos, acos, sqrt, pi

def orbital_height(a, e, f):
    """ Orbital height for orbital position f
    Input:
        - a [km]
        - e [-]
        - f [DEG]
    """
    return ( a*(1-e**2) )/( 1+e*cos(radians(f)) )

def orbital_position(a, e, r):
    """ Orbital position for orbital height r
    """
    return degrees(acos( (a*(1-e**2))/(e*r) - e**-1 ))

def orbital_velocity(mu, a, r):
    """ Vis Viva Equation
    """
    return sqrt( mu*(2/r-1/a) )

def circular_velocity(mu, r):
    """ Circular Velocity
    """
    return sqrt(mu/r)

def orbital_period(mu, a):
    """ Orbital Period in seconds
    """
    return 2*pi*sqrt(a**3/mu)


def plane_change(vi, alpha):
    """ Plane Change maneuver
    Input: 
        - Initial velocity: Vi [km/s]
        - Inclination change: alpha [deg]
        """
    return 2*vi*sin(radians(alpha)/2)

def general_coplanar(vi, vf, alpha):
    """ General Coplanar maneuver
    Input: 
        - Initial velocity: Vi [km/s]
        - Final velocity: Vf [km/s]
        - Inclination change: alpha [deg]
        """
    dv2 = vi**2 + vf**2 - 2*vi*vf*cos(radians(alpha))
    return sqrt(dv2)