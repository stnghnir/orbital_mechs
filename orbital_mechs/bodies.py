class Celestial:
    def __init__(self, mu, R) -> None:
        self.mu = mu # km3/s2: gravitational parameter
        self.R = R # km: radius
        self.orbit: Kepler = None

class Kepler:
    def __init__(self, a, e, i, L, cw, W) -> None:
        self.a = a # au
        self.e = e # -
        self.i = i # deg
        self.L = L # deg
        self.cw = cw # deg
        self.W = W # deg


sun = Celestial( 132_712_439_935.5, 696_000.0  )

mercury = Celestial( 22_032.1, 2439.7 )

venus = Celestial( 324_858.8, 6051.8 )
venus.orbit = Kepler( 0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255 )

earth = Celestial( 398_600.4, 6378.14 )
earth.orbit = Kepler(1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.)
moon =  Celestial(    4902.8, 1737.4  )

mars = Celestial( 42_828.3, 3397. )
mars.orbit = Kepler( 1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891 )

jupiter = Celestial( 126_711_995.4, 71_492. )

saturn = Celestial( 37_939_519.7, 60_268. )
saturn.orbit = Kepler(9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448)

uranus = Celestial( 5_780_158.5, 25_559. )

neptune = Celestial( 6_871_307.8, 24_764. )