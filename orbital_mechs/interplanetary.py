from math import sqrt
from .constants import AU
from .analytic import (
    circular_velocity as circ_v, 
    orbital_period, 
    orbital_velocity as visviva
    )
from .bodies import Celestial

class Interplanetary:
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