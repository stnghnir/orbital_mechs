import os
from math import sqrt
from functools import wraps

import numpy as np
ary = np.array
norm = np.linalg.norm

import spiceypy as spice
from lamberthub import izzo2015 as izzo

from .constants import AU
from .analytic import (
    circular_velocity as circ_v, 
    orbital_period, 
    orbital_velocity as visviva,
    hyperbolic_excess_velocity as hypr
    )
from .bodies import Celestial

BASE_MISSION = {"UTC", "Planet", "Orbit"}

KERNELS = [os.path.dirname(__file__)+"/kernels/de405.bsp",
           os.path.dirname(__file__)+"/kernels/leapseconds.tls"]
# SPICE Wrapper
def SPICE(fun):
    @wraps(fun)
    def wrapper(*args, **kwargs):
        spice.furnsh(KERNELS)
        try:
            result = fun(*args, **kwargs)
        finally:
            spice.kclear()
        return result
    return wrapper

@SPICE
def _getPlanetaryState(target, observer, et, frame="J2000") -> tuple:
    state, lt = spice.spkezr(target, et, frame, "None", observer)
    r = state[0:3]
    v = state[3:6]
    return r, v

@SPICE
def _UTC2ET(UTC):
    return spice.utc2et(UTC)

class Interplanetary2D:
    def __init__(self, central_body: Celestial, 
                 departure_body: Celestial, 
                 arrival_body: Celestial) -> None:
        
        self.mu: float = central_body.mu
        self.departure: Celestial = departure_body
        self.arrival: Celestial = arrival_body

        self.a: float
        self.P: float
        self.departure_velocity: float
        self.arrival_velocity: float
        self.departure_c3: float
        self.arrival_c3: float
        self.departure_deltaV: float
        self.arrival_deltaV: float
    
    def __str__(self) -> str:
        return (f"Departure ΔV: {self.departure_deltaV:.4f} km/s\n"
            f"Arrival ΔV: {self.arrival_deltaV:.4f} km/s\n"
            f"Combined ΔV: {self.departure_deltaV + self.arrival_deltaV:.4f} km/s\n"
            f"Travel time is {self.P/365.25/86400:.4f} years")

    def calc_maneuver(self, departure_height, arrival_height) -> float:
        """ Calculate interplanetary maneuver
        Departure orbit
        Arrival orbit
        Hohmann and C3s
        Delta V
        """
        self.departure.r = self.departure.orbit.a*AU
        self.arrival.r = self.arrival.orbit.a*AU

        self.departure.v = circ_v(self.mu, self.departure.r)
        self.arrival.v = circ_v(self.mu, self.arrival.r)

        self.a = (self.departure.r + self.arrival.r) / 2
        self.P = orbital_period(self.mu, self.a) / 2

        self.departure_velocity = visviva(self.mu, self.a, self.departure.r)
        self.arrival_velocity = visviva(self.mu, self.a, self.arrival.r)

        self.departure_c3 = (self.departure_velocity - self.departure.v)**2
        self.arrival_c3 = (self.arrival_velocity - self.arrival.v)**2

        self.departure_deltaV = sqrt(2*self.departure.mu/departure_height + self.departure_c3) - circ_v(self.departure.mu, departure_height)
        self.arrival_deltaV = sqrt(2*self.arrival.mu/arrival_height + self.arrival_c3) - circ_v(self.arrival.mu, arrival_height)

        return self.departure_deltaV + self.arrival_deltaV

class MissionState:
    def __init__(self) -> None:
        pass

class Transfer:
    def __init__(self, center, departure, arrival) -> None:
        self.center: Celestial = center
        self.departure: MissionState = departure
        self.arrival: MissionState = arrival

        self.C3 = None
        self.deltaV = None
        self.V_SC = None
        self.TOF = None

        self._calculate()

    def _calculate(self):
        departure_et = _UTC2ET(self.departure["UTC"])
        arrival_et = _UTC2ET(self.arrival["UTC"])

        R_Planet_Departure, V_Planet_Departure = _getPlanetaryState(self.departure["Planet"].name, self.center.name, departure_et)
        R_Planet_Arrival, V_Planet_Arrival = _getPlanetaryState(self.arrival["Planet"].name, self.center.name, arrival_et)

        self.TOF = arrival_et - departure_et

        V_SC_Departure, V_SC_Arrival = izzo(self.center.mu, R_Planet_Departure, R_Planet_Arrival, self.TOF)

        self.V_SC = {"Departure": V_SC_Departure,
                     "Arrival": V_SC_Arrival}

        self.C3 = {"Departure": norm(V_Planet_Departure - V_SC_Departure)**2,
              "Arrival": norm(V_Planet_Arrival - V_SC_Arrival)**2}
        
        self.deltaV = {"Departure": self._getInsertionDeltaV(self.departure, self.C3["Departure"]),
                "Arrival": self._getInsertionDeltaV(self.arrival, self.C3["Arrival"])}
    
    def _getInsertionDeltaV(self, mission, C3):
        mu = mission["Planet"].mu
        r = mission["Orbit"]
        if r:
            return hypr(mu, r, C3) - circ_v(mu, r)
        else:
            return None
    
    def __str__(self) -> str:
        return (
            "Transfer details:\n"
            f"C3 Departure: {self.C3["Departure"]} km2/s2\n" +
            f"C3 Arrival: {self.C3["Arrival"]} km2/s2\n" + 
            f"Departure delta V: {self.deltaV["Departure"]} km/s\n" + 
            f"Arrival delta V: {self.deltaV["Arrival"]} km/s"
            )