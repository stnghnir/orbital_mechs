""" Analytic formulas
"""
# Conic sections
from math import radians, degrees, sin, cos, tan, asin, acos, atan, sqrt, pi
import numpy as np
norm = np.linalg.norm
head = lambda vec: vec/norm(vec)
ary = np.array

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

def orbital_angular_momentum(mu, a, e):
    """ Orbital Angular Momentum
    """
    return sqrt(mu*a*(1-e**2))

def orbital_period(mu, a):
    """ Orbital Period in seconds
    """
    return 2*pi*sqrt(a**3/mu)

def circular_velocity(mu, r):
    """ Circular Velocity
    """
    return sqrt(mu/r)

def getTimeAtf(mu, a, e, f):
    """ Application of the kepler equation
    """
    f = radians(f)
    tanE2 = sqrt((1-e)/(1+e)) * tan(f/2)
    if E:=(2*atan(tanE2))<0:E+=2*pi
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

# Conversion of orbit representations
def state2kep(mu, vr, vv, angles="RAD"):
    """ Converts state vector to kepler elements
    a: Semi-Major Axis [km]
    e: Eccentricity [-]
    i: Inclination [rad]
    RAAN: Longitude / Right Ascension of the Ascending Node [rad]
    w: Argument of Periapsis [rad]
    f: True Anomaly [rad]

    With kwarg angles="DEG", angles are converted to degree
    """
    # Ensures np.array
    vr = ary(vr)
    vv = ary(vv)
    r = norm(vr)
    v = norm(vv)
    # Semi-Major-Axis
    a = r / (2-r*v**2/mu)
    # Eccentricity (and its vector )
    ve = (v**2/mu - 1/r) * vr - 1/mu * (vr @ vv) * vv
    e = norm(ve)
    # Specific Orbital Angular Momentum
    vh = np.cross(vr, vv)
    # Inclination
    i = np.arccos(vh[2] / norm(vh))
    # Nodal Vector
    n = np.cross(ary([0,0,1]), head(vh))
    # Right Ascension of the Ascending Node
    RAAN = np.arccos(n[0] / norm(n))
    if n[1] < 0: RAAN = 2*np.pi - RAAN
    # Argument of Periapsis
    w = np.arccos( n@ve / norm(n) / e)
    if ve[2] < 0: w = 2*np.pi - w
    # True Anomaly
    f = np.arccos( ve@vr / e / r)
    if np.dot(ve, vr) < 0: f = 2*np.pi - f
    # RAD 2 DEG
    if angles=="DEG":
        i = degrees(i)
        RAAN = degrees(RAAN)
        w = degrees(w)
        f = degrees(f)

    return a, e, i, RAAN, w, f

def kep2state(mu, a, e, i, RAAN, w, f, angles="RAD"):
    """ Converts kepler elements to state vector
    vr: Orbital Position as vector
    vv: Orbital Velocity as vector

    With kwarg angles="DEG", angles are treated as input in degree
    """
    if angles=="DEG":
        i = radians(i)
        RAAN = radians(RAAN)
        w = radians(w)
        f = radians(f)

    r = a*(1-e**2)/(1+e*cos(f))
    v = orbital_velocity(mu, a, r)
    h = orbital_angular_momentum(mu, a, e)

    theta = w+f
    vr = r * ary([
        cos(RAAN)*cos(theta) - sin(RAAN)*sin(theta)*cos(i),
        sin(RAAN)*cos(theta) + cos(RAAN)*sin(theta)*cos(i),
        sin(theta)*sin(i)
    ])
    vv = -mu/h * ary([
        cos(RAAN)*( sin(theta) + e*sin(w) ) + sin(RAAN)*( cos(theta) + e*cos(w) )*cos(i),
        sin(RAAN)*( sin(theta) + e*sin(w) ) - cos(RAAN)*( cos(theta) + e*cos(w) )*cos(i),
        -( cos(theta) + e*cos(w) )*sin(i)
    ])

    return vr, vv