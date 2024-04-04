# -*- coding: utf-8 -*-
# pylint: disable=missing-module-docstring, missing-class-docstring, missing-function-docstring
# pylint: disable=non-ascii-name

# sunpos - Compute the position of the sun relative to
# Earth's orbit, according to the current date and time and
# the observer's position on Earth.
#
# Aaron Mansheim
# created: 2011-06-22
#
# Based closely on: http://stjarnhimlen.se/comp/tutorial.html 2009-07-02
# The similar code sunriset.c by the same author is in the public domain.
#
# Ultimately my purpose is to make this as clear as possible to
# anyone who understands a) basics of the Python programming language
# and b) how the sine and cosine functions relate to a point on
# the circumference of the unit circle.

import dataclasses
from math import sqrt, hypot, degrees, fsum, fmod
from degrees import Degrees, sin, cos

#   Regions of this file:
# 1              coordinate_systems
# 2              time_and_angle
# 3 astronomical_time_and_angle (functions)
# 4 astronomical_coordinate_systems
# 5 astronomical_coordinate_transformations (functions)
# 6 astronomical_predictions (functions)
# 7              orbits
# 8              the_Moon (functions)

# 1 -------------------- coordinate_systems

@dataclasses.dataclass
class Cartesian2d:
    x: float
    y: float

    # Is not under test
    def tuple(self) -> tuple[float, float]:
        return self.x, self.y

    def polar(self) -> 'Polar':
        return Polar.from_cartesian2d(self)

    def rotate(self, θ: Degrees) -> 'Cartesian2d':
        return self.polar().rotate(θ).cartesian2d()

    @classmethod
    def from_polar(cls, p: 'Polar') -> 'Cartesian2d':
        x = p.r * p.θ.cos()
        y = p.r * p.θ.sin()

        self = cls.__new__(cls)
        self.__class__.__init__(self, x, y)
        return self

    @classmethod
    def from_eccentric(cls, E: Degrees, e) -> 'Cartesian2d':
        # Would this be easy to understand from a diagram?
        x = E.cos() - e
            # == radius * cos(θ)
        y = E.sin() * sqrt(1 - e ** 2)
            # == radius * sin(θ)

        self = cls.__new__(cls)
        self.__class__.__init__(self, x, y)
        return self

@dataclasses.dataclass
class Polar:
    r: float
    θ: Degrees

    def tuple(self) -> tuple[float, Degrees]:
        return self.r, self.θ

    def cartesian2d(self) -> Cartesian2d:
        return Cartesian2d.from_polar(self)

    def rotate(self, Δθ: Degrees) -> 'Polar':
        return Polar(self.r, Degrees(self.θ + Δθ).rev())

    @classmethod
    def from_cartesian2d(cls, p: Cartesian2d) -> 'Polar':
        r = hypot(p.x, p.y)
        θ = Degrees.atan2(p.y, p.x).rev()

        self = cls.__new__(cls)
        self.__class__.__init__(self, r, θ)
        return self

@dataclasses.dataclass
class Cartesian3d:
    x: float
    y: float
    z: float

    # Is not under test
    def tuple(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z

    def spherical(self) -> 'Spherical':
        return Spherical.from_cartesian3d(self)

    def cylindrical(self) -> 'Cylindrical':
        return Cylindrical.from_cartesian3d(self)

    def onto_xy_plane(self) -> Cartesian2d:
        return Cartesian2d(self.x, self.y)

    def rotate_about_x(self, θ: Degrees) -> 'Cartesian3d':
        p = Cartesian2d(self.y, self.z).rotate(θ)

        return self.__class__(self.x, p.x, p.y)

    def decline_about_y(self, θ: Degrees) -> 'Cartesian3d':
        # Equivalent to:
        # (lambda neg_y, x, z: (x, -neg_y, z)) (
        #   *Cartesian3d(-y, x, z).rotate_about_x(90 - θ).tuple()
        # )
        # This (-y, x, z) is just (x, y, z) rotated around z.

        point = Cartesian2d(self.x, self.z).rotate(Degrees(90 - θ))

        return self.__class__(point.x, self.y, point.y)

    @classmethod
    def from_cylindrical(cls, p) -> 'Cartesian3d':
        q = Polar(p.cylindrical_radius, p.longitude).cartesian2d()
        z = p.h

        self = cls.__new__(cls)
        self.__class__.__init__(self, q.x, q.y, z)
        return self

    @classmethod
    def from_spherical(cls, p) -> 'Cartesian3d':
        q = p.cylindrical().cartesian3d()

        self = cls.__new__(cls)
        self.__class__.__init__(self, q.x, q.y, q.z)
        return self

@dataclasses.dataclass
class Cylindrical:
    cylindrical_radius: float
    h: float
    longitude: Degrees

    # Is not under test
    def tuple(self) -> tuple[float, float, Degrees]:
        return self.cylindrical_radius, self.h, self.longitude

    def spherical(self) -> 'Spherical':
        return Spherical.from_cylindrical(self)

    def cartesian3d(self) -> Cartesian3d:
        return Cartesian3d.from_cylindrical(self)

    @classmethod
    def from_spherical(cls, p) -> 'Cylindrical':
        q = p.onto_meridian_plane().cartesian2d()

        self = cls.__new__(cls)
        self.__class__.__init__(self, q.x, q.y, p.longitude)
        return self

    @classmethod
    def from_cartesian3d(cls, p) -> 'Cylindrical':
        q = p.onto_xy_plane().polar()
        h = p.z
        return Cylindrical(q.r, h, q.θ)

@dataclasses.dataclass
class Spherical:
    radius: float
    latitude: Degrees
    longitude: Degrees

    def tuple(self) -> tuple[float, Degrees, Degrees]:
        return self.radius, self.latitude, self.longitude

    def cylindrical(self) -> Cylindrical:
        return Cylindrical.from_spherical(self)

    def cartesian3d(self) -> Cartesian3d:
        return Cartesian3d.from_spherical(self)

    # Is not under test
    def onto_equator_plane(self) -> Polar:
        return Polar(self.radius, self.longitude)

    def onto_meridian_plane(self) -> Polar:
        return Polar(self.radius, self.latitude)

    @classmethod
    def from_cylindrical(cls, p: Cylindrical) -> 'Spherical':
        q = Cartesian2d(p.cylindrical_radius, p.h).polar()

        self = cls.__new__(cls)
        self.__class__.__init__(self, q.r, q.θ, p.longitude)
        return self

    @classmethod
    def from_cartesian3d(cls, p: Cartesian3d) -> 'Spherical':
        return p.cylindrical().spherical()

# 2 -------------------- time_and_angle

@dataclasses.dataclass
class Date:
    year: int
    month: int
    day: int
    def day_number(self) -> float:
        """
        This formula models the common astronomical calendar,
        which coincides with the Gregorian calendar for years 1901-2099.
        On the other hand it is actually a Julian calendar in that it has
        no exceptions to the rule that every fourth year is a leap year.
        """
        return (
            367 * self.year
            - (7 * (self.year + (self.month + 9) // 12)) // 4
            + (275 * self.month) // 9
            + self.day
            - 730530)

class Seconds(float):
    def minutes(self) -> float:
        return self / 60.0

    @classmethod
    def fromhours(cls, hours: float) -> 'Seconds':
        return cls(3600 * hours)

    @classmethod
    def fromminutes(cls, minutes: float) -> 'Seconds':
        return cls(60 * minutes)

class Arcdegrees(Degrees):
    def hours(self) -> float:
        return self * 24 / 360.0

    def arcminutes(self) -> float:
        return 60 * self

    def arcseconds(self) -> Seconds:
        return Seconds(3600 * self)

    @classmethod
    def fromhours(cls, hours: float) -> 'Arcdegrees':
        return cls(hours * 360 / 24.0).rev()

# 3 -------------------- astronomical_time_and_angle

def toGMST0(Ls: Degrees) -> Arcdegrees:
    """
    Sidereal time (in other words, right ascension in hours)
    at the 00:00 meridian at Greenwich right now.
    """
    return Arcdegrees(Degrees(Ls + 180).rev()).hours()

def to_sidereal_time(GMST0, UT, terrestrial_longitude: Degrees) -> float:
    """
    Locally adjusted sidereal time: more info at function toGMST0.
    """
    return fmod(GMST0 + UT + Arcdegrees(terrestrial_longitude).hours(), 24.0)

def hour_angle(sidereal_time, RA) -> float:
    return fmod(sidereal_time - RA, 24.0)

def obliquity_of_the_ecliptic(day_number) -> float:
    """
    Inclination of Earth's axis of rotation
    to its plane of revolution.
    """
    return 23.4393 - 3.563E-7 * day_number

def ecliptic_anomaly_to_longitude(
        longitude_of_periapsis: Degrees,
        angle_from_periapsis: Degrees,
) -> Degrees:
    """
    Ecliptic longitude, sometimes written simply as
    "longitude", has its zero at the ascending node.
    For the Sun-Earth system, the ascending node is
    the vernal point (if I'm not mistaken!?).

    In contrast to longitude, the mean anomaly and true
    anomaly have their zero at periapsis. For the
    Sun-Earth system, the periapsis is perihelion.

    When we say "mean longitude" of the Sun, we mean
    "mean ecliptic longitude" with zero longitude at the
    vernal point. We are not referring to the terrestrial
    longitude or to the celestial longitude (also known
    as right ascension) that are zero at a meridian.

    We should also explain "mean longitude" vs. "longitude".
    The mean longitude progresses quite cyclically,
    while the (true) longitude goes elliptically.
    Both complete one revolution in the same time.
    True longitude measures the actual accelerating
    and decelerating position of a body. Mean longitude
    averages out the position over the whole revolution
    to measure out a fictitious steady motion.
    """
    return Degrees(longitude_of_periapsis + angle_from_periapsis).rev()

def sun_earth_elements(day_number) -> tuple[Degrees, float, float, Degrees]:
    """
    Why this function is named "sun_earth_elements":
    We're primarily concerned with the Earth-based observer,
    so we'll use these orbital elements to determine how the Sun
    moves through Earth's sky. Of course, outside of that frame
    it is more sensible to interpret these elements vice versa,
    with Earth moving around the Sun.
    """

    ω = Degrees(282.9404 + 4.70935E-5 * day_number).rev()   # longitude of perihelion
    a = 1.000000                                   # mean_distance, in a.u.
    e = 0.016709 - 1.151E-9 * day_number           # eccentricity
    Ms = Degrees(356.0470 + 0.9856002585 * day_number).rev() # mean anomaly
    return (ω, a, e, Ms)

def eccentric_anomaly_first_approximation(M: Degrees, e) -> Degrees:
    """
    This truncated Taylor series is claimed accurate enough for a small
    eccentricity such as that of the Sun-Earth orbit (0.017).
    Properly the eccentric anomaly is the solution E
    of Kepler's equation in mean anomaly M and eccentricity e:
        M = E - e sin(E).
    Thanks to Paul Schlyter himself for clarifying this in a
    private email which, frankly, I have not yet fully analyzed.
    """
    return Degrees(
        M
        + degrees(e * M.sin())
        * (1 + e * M.cos())).rev()

    # Note if we write E(n+1) = M + e sin (E(n)) and iterate,
    # we get this, which is supposedly equivalent!?
    # E(1) = M + e sin (E(0))
    # E(2) = M + e sin ( M + e sin (E(0)) )
    #      = M + e sin (M) cos (e sin (E(0))) + e cos (M) sin (e sin (E(0)) )
    # E(3) = M + e sin ( M + e sin ( M + e sin ( E(0) ) ) )

def iterate_for_eccentric_anomaly(M, e, E0: Degrees) -> Degrees:
    # Approximately invert Kepler's equation M = E - e sin E
    # by Newton's method, E1 = E0 - f(E0) / f'(E0),
    # where f(E) = E - e sin(E) - M,
    # so that Newton's method approaches E such that 0 = E - e sin(E) - M.

    def iterate_once_for_eccentric_anomaly(M, e, E0: Degrees) -> Degrees:
        # Take one step of Newton's method, E1 = E0 - f(E0) / f'(E0).
        # Here f(E) = E - e sin(E) - M, so that f'(E) = 1 - e cos(E). That is,
        # E1 = E0 - (E - e sin(E) - M) / (1 - e cos(E)).
        radians = E0 - (E0 - degrees(e * E0.sin()) - M) / (1 - e * E0.cos())
        return Degrees(radians)

    E1: Degrees = E0
    E0 = Degrees(0)
    k = 0
    while abs(E1 - E0) >= 0.005 and k < 30:
        k += 1
        E0 = E1
        E1 = iterate_once_for_eccentric_anomaly(M, e, E0)
    return E1

# 4 -------------------- astronomical_coordinate_systems

@dataclasses.dataclass
class Celestial:
    RA: float
    Decl: Degrees

@dataclasses.dataclass
class Geographic:
    latitude: Degrees
    longitude: Degrees

@dataclasses.dataclass
class Horizontal:
    altitude: float
    azimuth: float

    def tuple(self) -> tuple[float, float]:
        return (self.altitude, self.azimuth)

    @classmethod
    def from_celestial(cls, mlon: Degrees,
                       celestial: Celestial,
                       hours_UT,
                       geographic: Geographic) -> 'Horizontal':
        _distance, altitude, azimuth = Spherical(
                1.0,
                celestial.Decl,
                Arcdegrees.fromhours(hour_angle(
                    to_sidereal_time(toGMST0(mlon), hours_UT, geographic.longitude),
                    Arcdegrees(celestial.RA).hours(),
                )),
            ).cartesian3d().decline_about_y(geographic.latitude).spherical().tuple()

        azimuth = Degrees(azimuth + 180).rev()
        return cls(altitude, azimuth)

# 5 -------------------- astronomical_coordinate_transformations

def ecliptic_to_celestial(distance, latitude, longitude: Degrees, obliquity):
    (distance1, Decl, RA) = Spherical(
            distance,
            latitude,
            longitude
        ).cartesian3d().rotate_about_x(obliquity).spherical().tuple()
    return (distance1, Decl, RA)

def position_from_plane_of_orbit_to_ecliptic(polar: Polar, Ω: Degrees, i: Degrees, ω):
    # TODO: break this down to easily understood operations
    #   using angle sum formulas etc.
    #                           _ _     _ _     _
    # cos(-x) = cos(x)        1 ┆╳ ╲   ╱┆╳ ╲   ╱┆ cos
    # sin(-x) = -sin(x)       0 ╱┄╲┄╲┄╱┄╱┄╲┄╲┄╱┄╱ sin
    # sin(x) = cos(x - 90)   -1 ┆  ╲_╳_╱┆  ╲_╳_╱┆
    # sin(x - 90) = -cos(x)    -2pi     0       2pi
    #
    # a - b = a + (-b)
    # cos(a + b) = cos(a)cos(b) - sin(a)sin(b)
    # cos(a - b) = cos(a)cos(b) + sin(a)sin(b)
    # sin(a + b) = sin(a)cos(b) + cos(a)sin(b)
    # sin(a - b) = sin(a)cos(b) - cos(a)sin(b)

    return Cartesian3d(
        polar.r * (
            Ω.cos() * cos(polar.θ + ω)
            - Ω.sin() * sin(polar.θ + ω) * i.cos()
        ),
        polar.r * (
            Ω.sin() * cos(polar.θ + ω)
            + Ω.cos() * sin(polar.θ + ω) * i.cos()
        ),
        polar.r * sin(polar.θ + ω) * i.sin(),
    )

def sun_earth_ecliptic_to_celestial(distance, true_lon: Degrees, obliquity: Degrees):
    # equivalent to:
    # ecliptic_to_celestial(distance, 0.0, true_lon, obliquity)

    p = Polar(distance, true_lon).cartesian2d()
    q = Cartesian3d(p.x, p.y, 0.0).rotate_about_x(obliquity).spherical()
    return Celestial(RA=q.longitude, Decl=q.latitude)

# 6 -------------------- astronomical_predictions

# return: (distance, true_longitude, Ls, obliquity, )
def date_to_sun_earth_ecliptic(date: Date):

    # return: (ω, a, e, M, obliquity, )
    def sun_earth_elements_and_obliquity(
            date: Date) -> tuple[Degrees, float, float, Degrees, float]:
        d = date.day_number()
        return sun_earth_elements(d) + (obliquity_of_the_ecliptic(d), )

    # return: Polar(distance, ecliptic_angle)
    def ecliptic_polar(Ms: Degrees, e) -> Polar:
        return Cartesian2d.from_eccentric(
            eccentric_anomaly_first_approximation(Ms, e),
            e,
        ).polar()

    # return: (distance, true_longitude, Ls, obliquity, )
    def longitudes(Ms: Degrees, e, ω: Degrees, obliquity: Degrees):
        p = ecliptic_polar(Ms, e)
        (distance, ν) = p.r, p.θ
        return (
            distance,
            ecliptic_anomaly_to_longitude(ω, ν),
            ecliptic_anomaly_to_longitude(ω, Ms),
            obliquity,
        )

    # return: (distance, true_longitude, Ls, obliquity, )
    def date_to_sun_earth_ecliptic_1(date: Date):
        (ω, _a, e, Ms, obliquity, ) = sun_earth_elements_and_obliquity(
            date
        )
        return longitudes(Ms, e, ω, obliquity)

    return date_to_sun_earth_ecliptic_1(date)

def time_and_location_to_sun_horizontal(date: Date,
                                         hours_UT, geographic: Geographic):

    def sun_earth_ecliptic_to_celestial_for_horizontal(
        distance0,
        true_longitude,
        Ls,
        obliquity: Degrees,
    ):
        return (Ls, sun_earth_ecliptic_to_celestial(
            distance0,
            true_longitude,
            obliquity,
        ), hours_UT, geographic)

    (altitude, azimuth) = Horizontal.from_celestial(
        *sun_earth_ecliptic_to_celestial_for_horizontal(
            *date_to_sun_earth_ecliptic(date)
        )
    ).tuple()

    return (altitude, azimuth)

# 7 -------------------- orbits

@dataclasses.dataclass
class OrbitalElements:
    '''
    (☉,♈︎) Let ☉∈R³. Without loss of generality, let ☉ be the origin.
    Let ♈︎∈R³ and ☉≠♈︎. Call ♈︎ the vernal point.

    (☊) Let N∈R³ and N⟂♈︎. Then there exists a unique plane P such that ☉∈P and P⟂N.
    We will use P as a reference plane. Let ☊∈P and ☊≠☉.
    Call ☊ the ascending node.

    (m∠aob) Recall that for all a,b,o∈R³ such that a≠o and b≠o,
    m∠aob=arccos(â₀⋅b̂₀), where a₀=a-o, b₀=b-o,
    â₀=a₀/‖a₀‖, and b̂₀=b₀/‖b₀‖. Recall that m∠aob∈[0,π].

    (Ω,i) Define the inclination i as follows.
    Let Ω=m∠♈︎☉☊ . Call Ω the longitude of the ascending node.
    Let I∈R³, such that I⟂☊ and I⋅N>0.
    Let i=m∠(N×☊)☉I so that i∈[0,π).

    Define the orbital plane as the plane containing points ☉, ☊, and I.

    (e,a,ω) Define an orbit as the ellipse in the orbital plane having
    a focus at ☉, eccentricity e, measure a of semimajor axis,
    and angle ω from ☊ to periapsis.

    Define the eccentricity of the orbit as e = sqrt(1 - (b/a)²),
    where b is measure of semiminor axis.

    Define the periapsis as the vector from ☉ to the nearest point on the orbit.

    (ν) Specify the location of the body in the orbit by the true anomaly ν.
    The true anomaly is the angle at ☉ between periapsis and the body.

    (E) Alternatively, specify the location of the body in the orbit by the
    eccentric anomaly E, which is related to ν as follows:

    tan E = sqrt(1-e²) sin ν / (e + cos ν).

    (M) Alternatively, specify the location of the body in the orbit by the
    mean anomaly M, which is related to E as follows:

    M = E - e sin E.
    '''
    Ω: Degrees
    i: Degrees
    ω: Degrees
    a: float
    e: float
    M: Degrees

# 8 -------------------- the_Moon

def moon_elements(d) -> OrbitalElements:
    return OrbitalElements(
        Ω = Degrees(125.1228 - 0.0529538083 * d).rev(),  # Long asc. node
        i = Degrees(5.1454),                             # Inclination
        ω = Degrees(318.0634 + 0.1643573223 * d).rev(),  # Arg. of perigee
        a = 60.2666,                           # Mean distance,
                                               #   in Earth equatorial radii
        e = 0.054900,                          # Eccentricity
        M = Degrees(115.3654 + 13.0649929509 * d).rev(), # Mean anomaly
    )

def moon_perturbation_arguments(Ls, moon_N, moon_ω, Mm, Ms) -> tuple[
                Degrees, Degrees, Degrees, Degrees, Degrees, Degrees]:
    Lm = moon_N + moon_ω + Mm # moon mean longitude
    (
        D,  # moon mean elongation
        F,  # moon argument of latitude
    ) = (
        Degrees(Lm - Ls).rev(),
        Degrees(Lm - moon_N).rev(),
    )
    return (
            Ls,
            Degrees(Lm).rev(),
            Ms, # sun mean anomaly
            Mm, # moon mean anomaly
            D,
            F)

def moon_2_distance_perturbation_terms( _Ls, _Lm, _Ms, Mm, D, _F) -> tuple[float, float]:
    # all & only the terms larger than 0.1 Earth radii
    return (
        -0.58 * cos(Mm - 2 * D),
        -0.46 * cos(2 * D),
    )

def moon_5_latitude_perturbation_terms(
        _Ls, _Lm, _Ms, Mm, D, F) -> tuple[ float, float, float, float, float]:
    # All & only the terms larger than 0.01 degrees
    # aiming for 1-2 arcmin accuracy.
    # Even with the first term the error is usually under 0.15 degree.
    return (
        -0.173 * sin(F - 2 * D),
        -0.055 * sin(Mm - F - 2 * D),
        -0.046 * sin(Mm + F - 2 * D),
        +0.033 * sin(F + 2 * D),
        +0.017 * sin(2 * Mm + F),
    )

def moon_12_longitude_perturbation_terms( _Ls, _Lm, Ms, Mm, D, F) -> tuple[
        float, float, float, float, float, float,
        float, float, float, float, float, float]:
    # All & only the terms larger than 0.01 degrees,
    # aiming for 1-2 arcmin accuracy.
    # Even with the first two terms the error is usually under 0.25 degree.
    # These might come from EPL.
    return (
        -1.274 * sin(Mm - 2 * D),      # Evection (Ptolemy)
        +0.658 * sin(2 * D),           # Variation (Tycho Brahe)
        -0.186 * sin(Ms),              # Yearly equation (Tycho Brahe)
        -0.059 * sin(2 * Mm - 2 * D),
        -0.057 * sin(Mm - 2 * D + Ms),
        +0.053 * sin(Mm + 2 * D),
        +0.046 * sin(2 * D - Ms),
        +0.041 * sin(Mm - Ms),
        -0.035 * sin(D),               # Parallactic equation
        -0.031 * sin(Mm + Ms),
        -0.015 * sin(2 * F - 2 * D),
        +0.011 * sin(Mm - 4 * D),
    )

def cartesian2d_in_plane_of_orbit(E: Degrees, e, distance):
    return Cartesian2d(
        distance * (E.cos() - e),
        distance * sqrt(1 - e ** 2) * E.sin(),
    )

def date_and_sun_mean_to_moon_ecliptic(date: Date,
                                       Ls: Degrees, M: Degrees):
    # TODO: Generalize. Much of this is not specific to the moon.

    def moon_elliptic_to_polar(a, e, M):
        return cartesian2d_in_plane_of_orbit(
                    iterate_for_eccentric_anomaly(
                        M,
                        e,
                        eccentric_anomaly_first_approximation(M, e)
                    ),
                    e,
                    a
                ).polar()

    def moon_elements_to_spherical(elems: OrbitalElements) -> Spherical:
        return Spherical.from_cartesian3d(position_from_plane_of_orbit_to_ecliptic(
                moon_elliptic_to_polar(elems.a, elems.e, elems.M),
                elems.Ω,
                elems.i,
                elems.ω
                ))

    def moon_elements_to_spherical_perturbation(elems: OrbitalElements) -> map:
        args = moon_perturbation_arguments(Ls, elems.Ω, elems.ω, elems.M, M)
        return map (
            fsum,
            (
                moon_2_distance_perturbation_terms(*args),
                moon_5_latitude_perturbation_terms(*args),
                moon_12_longitude_perturbation_terms(*args),
            )
        )

    def date_and_sun_mean_to_moon_ecliptic_1(date: Date):
        elems = moon_elements(date.day_number())
        return map(sum, zip(
            *(
                moon_elements_to_spherical(elems).tuple(),
                moon_elements_to_spherical_perturbation(elems),
            )
        ))

    (pdist, plat, plon) = date_and_sun_mean_to_moon_ecliptic_1(date)
    return (pdist, plat, plon)
