""" Analytic formulas
"""
# Conic sections
from math import radians, degrees, sin, cos, tan, asin, acos, atan, sqrt, pi

def orbital_height(a, e, f):
    """ Orbital height for orbital position f
    Input:
    :a [km]
    :e [-]
    :f [DEG]
    """
    f = radians(f)
    return ( a*(1-e**2) )/( 1+e*cos(f) )

def orbital_position(a, e, r):
    """ Orbital position for orbital height r
    """
    return degrees( acos( (a*(1-e**2))/(e*r) - e**-1 ) )

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

def getTimeAtf(mu, a, e, f):
    """ Application of the kepler equation
    """
    f = radians(f)
    tanE2 = sqrt((1-e)/(1+e)) * tan(f/2)
    E = 2*atan(tanE2)
    if E<0:E+=2*pi
    M = E - e*sin(E)
    return M*sqrt(a**3/mu)

# Maneuvers
def plane_change(vi, alpha):
    """ Plane Change maneuver
    Input: 
        - Initial velocity: Vi [km/s]
        - Inclination change: alpha [deg]
        """
    alpha = radians(alpha)
    return 2*vi*sin(alpha/2)

def general_coplanar(vi, vf, alpha):
    """ General Coplanar maneuver
    Input: 
        - Initial velocity: Vi [km/s]
        - Final velocity: Vf [km/s]
        - Inclination change: alpha [deg]
        """
    alpha = radians(alpha)
    dv2 = vi**2 + vf**2 - 2*vi*vf*cos(alpha)
    return sqrt(dv2)

def flight_path_angle(e, f):
    """ Flight path angle
    Inputs and outputs in Degree
    """
    f = radians(f)
    tany = (e*sin(f))/(1+e*cos(f))
    return degrees(atan(tany))

# Hyperbolic
def hyperbolic_excess_velocity(mu, r, C3):
    return sqrt(2*mu/r + C3)