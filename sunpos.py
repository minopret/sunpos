# -*- coding: utf-8 -*-
# pylint: disable=invalid-name, non-ascii-name

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

from math import sqrt, hypot, pi, fsum
from degrees import sin, cos, atan2, rev

class Cartesian2d:
    def __init__(self, x, y):
        self.x, self.y = x, y

    # Is not under test
    def as_tuple(self):
        return self.x, self.y

    def to_polar(self):
        return Polar.from_cartesian2d(self)

    @classmethod
    def from_polar(cls, p):
        x = p.r * cos(p.θ)
        y = p.r * sin(p.θ)
        self = cls.__new__(cls)
        self.__init__(x, y)
        return self

    @classmethod
    def from_eccentric(cls, eccentric_anomaly, eccentricity):
        # This would be easy to understand from a simple diagram.
        x = cos(eccentric_anomaly) - eccentricity
            # == radius * cos(θ)
        y = sin(eccentric_anomaly) * sqrt(1 - eccentricity ** 2)
            # == radius * sin(θ)
        self = cls.__new__(cls)
        self.__init__(x, y)
        return self

class Cartesian3d:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z

    # Is not under test
    def as_tuple(self):
        return self.x, self.y, self.z

    def to_spherical(self):
        return Spherical.from_cartesian3d(self)

    def to_cylindrical(self):
        return Cylindrical.from_cartesian3d(self)

    def onto_xy_plane(self):
        return Cartesian2d(self.x, self.y)

    def rotate_about_x(self, θ):
        p = Cartesian2d(self.y, self.z).to_polar().rotate(θ).to_cartesian2d()
        return Cartesian3d(self.x, p.x, p.y)

    def decline_about_y(self, θ):
        # Equivalent to:
        # (lambda neg_y, x, z: (x, -neg_y, z)) (
        #   *Cartesian3d(-y, x, z).rotate_about_x(90 - θ).as_tuple()
        # )
        # We could have written (y, x, z), but (-y, x, z) has the
        # nice property that it is just (x, y, z) rotated around z.

        point = Cartesian2d(self.x, self.z).to_polar().rotate(90 - θ).to_cartesian2d()
        return Cartesian3d(point.x, self.y, point.y)

    @classmethod
    def from_cylindrical(cls, p):
        q = Polar(p.cylindrical_radius, p.longitude).to_cartesian2d()
        z = p.height
        self = cls.__new__(cls)
        self.__init__(q.x, q.y, z)
        return self

    @classmethod
    def from_spherical(cls, p):
        q = p.to_cylindrical().to_cartesian3d()
        self = cls.__new__(cls)
        self.__init__(q.x, q.y, q.z)
        return self

class Cylindrical:
    def __init__(self, cylindrical_radius, height, longitude):
        self.cylindrical_radius, self.height, self.longitude = cylindrical_radius, height, longitude

    # Is not under test
    def as_tuple(self):
        return self.cylindrical_radius, self.height, self.longitude

    def to_spherical(self):
        return Spherical.from_cylindrical(self)

    def to_cartesian3d(self):
        return Cartesian3d.from_cylindrical(self)

    @classmethod
    def from_spherical(cls, p):
        q = p.onto_meridian_plane().to_cartesian2d()
        self = cls.__new__(cls)
        self.__init__(q.x, q.y, p.longitude)
        return self

    @classmethod
    def from_cartesian3d(cls, p):
        q = p.onto_xy_plane().to_polar()
        height = p.z
        return Cylindrical(q.r, height, q.θ)

class OrbitalElements:
    '''
    Assume a point O in R³. Assume a nonzero vector ♈︎in R³.
    Let N be a vector perpendicular to ♈︎.
    Define P as the plane through O and normal to N.

    Let ☊ be a nonzero vector in P.
    Let Ω ∈ [0, 2π) be the angle from ♈︎to ☊.

    Let I be a vector perpendicular to ☊, with non-negative projection on N.
    Let i ∈ [0, π) be the angle from N × ☊ to I.

    Define the orbital plane as the plane containing points O, O + ☊, and O + I.

    Define the orbit as the ellipse in the orbital plane having
    a focus at O, eccentricity e, measure a of semimajor axis,
    and angle ω from ☊ to periapsis.

    Periapsis is the vector from O to the nearest point on the ellipse.

    e = sqrt(1 - (b/a)²), where b is measure of semiminor axis.

    Specify the location of the body in the orbit by the true anomaly ν.
    The true anomaly is the angle at O between periapsis and the body.

    Alternatively, specify the location of the body in the orbit by the
    eccentric anomaly E, or by the mean anomaly M, which quantities are
    related to ν as follows:

    tan E = sqrt(1-e²) sin ν / (e + cos ν).

    M = E - e sin E.
    '''
    def __init__(self, Ω, i, ω, a, e, M):
        self.Ω, self.i, self.ω, self.a, self.e, self.M = Ω, i, ω, a, e, M

class Polar:
    def __init__(self, r, θ):
        self.r, self.θ = r, θ

    def as_tuple(self):
        return self.r, self.θ

    def to_cartesian2d(self):
        return Cartesian2d.from_polar(self)

    def rotate(self, Δθ):
        return Polar(self.r, rev(self.θ + Δθ))

    @classmethod
    def from_cartesian2d(cls, p):
        r = hypot(p.x, p.y)
        θ = rev(atan2(p.y, p.x))
        self = cls.__new__(cls)
        self.__init__(r, θ)
        return self

class Spherical:
    def __init__(self, radius, latitude, longitude):
        self.radius, self.latitude, self.longitude = radius, latitude, longitude

    def as_tuple(self):
        return self.radius, self.latitude, self.longitude

    def to_cylindrical(self):
        return Cylindrical.from_spherical(self)

    def to_cartesian3d(self):
        return Cartesian3d.from_spherical(self)

    # Is not under test
    def onto_equator_plane(self):
        return Polar(self.radius, self.longitude)

    def onto_meridian_plane(self):
        return Polar(self.radius, self.latitude)

    @classmethod
    def from_cylindrical(cls, p):
        q = Cartesian2d(p.cylindrical_radius, p.height).to_polar()
        self = cls.__new__(cls)
        self.__init__(q.r, q.θ, p.longitude)
        return self

    @classmethod
    def from_cartesian3d(cls, p):
        return p.to_cylindrical().to_spherical()



def arcdegrees_to_hours(arcdegrees):
    return arcdegrees * 24 / 360.0

def hours_to_arcdegrees(hours):
    return rev(hours * 360 / 24.0)

def arcdegrees_to_arcminutes(arcdegrees):
    return 60 * arcdegrees

# Is not under test
def arcdegrees_to_arcseconds(arcdegrees):
    return 3600 * arcdegrees

def hours_to_seconds(hours):
    return 3600 * hours

def minutes_to_seconds(minutes):
    return 60 * minutes

def seconds_to_minutes(seconds):
    return seconds / 60.0

def day_number(year, month, day):
    """
    This formula models the common astronomical calendar,
    which coincides with the Gregorian calendar for years 1901-2099.
    On the other hand it is actually a Julian calendar in that it has
    no exceptions to the rule that every fourth year is a leap year.
    """
    return (
        367 * year
        - (7 * (year + (month + 9) // 12)) // 4
        + (275 * month) // 9
        + day
        - 730530)



def GMST0(mean_longitude):
    """
    Sidereal time (in other words, right ascension in hours)
    at the 00:00 meridian at Greenwich right now.
    """
    return arcdegrees_to_hours(rev(mean_longitude + 180))

def sidereal_time(GMST0, UT, terrestrial_longitude):
    """
    Locally adjusted sidereal time: more info at function GMST0.
    """
    # TODO: normalize to 24 hours
    return GMST0 + UT + arcdegrees_to_hours(terrestrial_longitude)

def hour_angle(sidereal_time, RA):
    # TODO: normalize to 24 hours
    return sidereal_time - RA



def cartesian2d_in_plane_of_orbit(eccentric_anomaly, eccentricity, distance):
    (x, y) = (
        distance * (cos(eccentric_anomaly) - eccentricity),
        distance * sqrt(1 - eccentricity ** 2) * sin(eccentric_anomaly),
    )
    return (x, y)

def eccentric_anomaly_first_approximation(mean_anomaly, eccentricity):
    """
    This truncated Taylor series is claimed accurate enough for a small
    eccentricity such as that of the Sun-Earth orbit (0.017).
    Properly the eccentric anomaly is the solution E
    of Kepler's equation in mean anomaly M and eccentricity e:
        M = E - e sin(E).
    Thanks to Paul Schlyter himself for clarifying this in a
    private email which, frankly, I have not yet fully analyzed.
    """
    return rev(
        mean_anomaly
        + (180 / pi) * eccentricity * sin(mean_anomaly)
        * (1 + eccentricity * cos(mean_anomaly)))

    # Note if we write E(n+1) = M + e sin (E(n)) and iterate,
    # we get this, which is supposedly equivalent!?
    # E(1) = M + e sin (E(0))
    # E(2) = M + e sin ( M + e sin (E(0)) )
    #      = M + e sin (M) cos (e sin (E(0))) + e cos (M) sin (e sin (E(0)) )
    # E(3) = M + e sin ( M + e sin ( M + e sin ( E(0) ) ) )

def ecliptic_anomaly_to_longitude(
        longitude_of_periapsis,
        angle_from_periapsis,
):
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
    return rev(longitude_of_periapsis + angle_from_periapsis)

def ecliptic_to_celestial(distance, latitude, longitude, oblecl):
    (distance1, Decl, RA) = Spherical(
            distance,
            latitude,
            longitude
        ).to_cartesian3d().rotate_about_x(oblecl).to_spherical().as_tuple()
    return (distance1, Decl, RA)

def iterate_for_eccentric_anomaly(M, e, E0):

    def iterate_once_for_eccentric_anomaly(M, e, E0):
        return E0 - (E0 - (180 / pi) * e * sin(E0) - M) / (1 - e * cos(E0))

    """
    Approximately invert Kepler's equation M = E - e sin E
    by Newton's method.
    """
    E1 = E0
    E0 = 0
    k = 0
    while abs(E1 - E0) >= 0.005 and k < 30:
        k = k + 1
        E0 = E1
        E1 = iterate_once_for_eccentric_anomaly(M, e, E0)
    return E1

def position_from_plane_of_orbit_to_ecliptic(r, θ, Ω, i, ω):
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

    (xeclip, yeclip, zeclip) = (
        r * (
            cos(Ω) * cos(θ + ω)
            - sin(Ω) * sin(θ + ω) * cos(i)
        ),
        r * (
            sin(Ω) * cos(θ + ω)
            + cos(Ω) * sin(θ + ω) * cos(i)
        ),
        r * sin(θ + ω) * sin(i),
    )
    return (xeclip, yeclip, zeclip)

def sun_earth_celestial_to_alt_azimuth(mlon, Decl, RA, hours_UT, lat, lon):
    distance, altitude, azimuth = Spherical(
            1.0,
            Decl,
            hours_to_arcdegrees(hour_angle(
                sidereal_time(GMST0(mlon), hours_UT, lon),
                arcdegrees_to_hours(RA),
            )),
        ).to_cartesian3d().decline_about_y(lat).to_spherical().as_tuple()

    azimuth = rev(azimuth + 180)
    return (altitude, azimuth)

def sun_earth_ecliptic_to_celestial(distance, true_lon, oblecl):
    # equivalent to:
    # ecliptic_to_celestial(distance, 0.0, true_lon, oblecl)

    p = Polar(distance, true_lon).to_cartesian2d()
    q = Cartesian3d(p.x, p.y, 0.0).rotate_about_x(oblecl).to_spherical()
    distance1, Decl, RA = q.radius, q.latitude, q.longitude
    return (distance1, Decl, RA)



# return: (distance, true_longitude, mean_longitude, oblecl, )
def date_to_sun_earth_ecliptic(year, month, day):

    # return: (ω, a, e, M, oblecl, )
    def sun_earth_elements_and_oblecl(year, month, day):
        d = day_number(year, month, day)
        return sun_earth_elements(d) + (obliquity_of_the_ecliptic(d), )

    # return: Polar(distance, ecliptic_angle)
    def ecliptic_polar(M, e):
        return Cartesian2d.from_eccentric(
            eccentric_anomaly_first_approximation(M, e),
            e,
        ).to_polar()

    # return: (distance, true_longitude, mean_longitude, oblecl, )
    def longitudes(M, e, ω, oblecl):
        p = ecliptic_polar(M, e)
        (distance, true_anomaly) = p.r, p.θ
        return (
            distance,
            ecliptic_anomaly_to_longitude(ω, true_anomaly),
            ecliptic_anomaly_to_longitude(ω, M),
            oblecl,
        )

    # return: (distance, true_longitude, mean_longitude, oblecl, )
    def date_to_sun_earth_ecliptic_1(year, month, day):
        (ω, a_ignored, e, M, oblecl, ) = sun_earth_elements_and_oblecl(
            year,
            month,
            day,
        )
        return longitudes(M, e, ω, oblecl)

    return date_to_sun_earth_ecliptic_1(year, month, day)

def obliquity_of_the_ecliptic(day_number):
    """
    Inclination of Earth's axis of rotation
    to its plane of revolution.
    """
    return 23.4393 - 3.563E-7 * day_number

def sun_earth_elements(day_number):
    """
    Why this function is named "sun_earth_elements":
    We're primarily concerned with the Earth-based observer,
    so we'll use these orbital elements to determine how the Sun
    moves through Earth's sky. Of course, outside of that frame
    it is more sensible to interpret these elements vice versa,
    with Earth moving around the Sun.
    """

    ω = rev(282.9404 + 4.70935E-5 * day_number)    # longitude of perihelion
    a = 1.000000                                   # mean_distance, in a.u.
    e = 0.016709 - 1.151E-9 * day_number           # eccentricity
    M = rev(356.0470 + 0.9856002585 * day_number)  # mean anomaly
    return (ω, a, e, M)

def time_and_location_to_sun_alt_azimuth(year, month, day, hours_UT, lat, lon):

    def sun_earth_ecliptic_to_celestial_for_alt_azimuth(
        distance0,
        true_longitude,
        mean_longitude,
        oblecl,
    ):
        return (mean_longitude, ) + sun_earth_ecliptic_to_celestial(
            distance0,
            true_longitude,
            oblecl,
        )[1:3] + (hours_UT, lat, lon, )

    (altitude, azimuth) = sun_earth_celestial_to_alt_azimuth(
        *sun_earth_ecliptic_to_celestial_for_alt_azimuth(
            *date_to_sun_earth_ecliptic(
                year,
                month,
                day,
            )
        )
    )

    return (altitude, azimuth)



def moon_elements(d):
    Ω = rev(125.1228 - 0.0529538083 * d)   # Long asc. node
    i = 5.1454                             # Inclination
    ω = rev(318.0634 + 0.1643573223 * d)   # Arg. of perigee
    a = 60.2666                            # Mean distance,
                                           #   in Earth equatorial radii
    e = 0.054900                           # Eccentricity
    M = rev(115.3654 + 13.0649929509 * d)  # Mean anomaly
    return OrbitalElements(Ω, i, ω, a, e, M)

def moon_perturbation_arguments(mlon, moon_N, moon_ω, moon_M, M):
    Lm = moon_N + moon_ω + moon_M
    (
        Ls,  # Sun's mean longitude
        Ms,  # Sun's mean anomaly
        Mm,  # Moon's mean anomaly
        D,   # Moon's mean elongation
        F,   # Moon's argument of latitude
    ) = (
        mlon,
        M,
        moon_M,
        rev(Lm - mlon),
        rev(Lm - moon_N),
    )
    return (Ls, rev(Lm), Ms, Mm, D, F)

def moon_2_distance_perturbation_terms(Ls, Lm, Ms, Mm, D, F):
    # all & only the terms larger than 0.1 Earth radii
    return (
        -0.58 * cos(Mm - 2 * D),
        -0.46 * cos(2 * D),
    )

def moon_5_latitude_perturbation_terms(Ls, Lm, Ms, Mm, D, F):
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

def moon_12_longitude_perturbation_terms(Ls, Lm, Ms, Mm, D, F):
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

def date_and_sun_mean_to_moon_ecliptic(year, month, day, mean_longitude, mean_anomaly):
    # TODO: Generalize. Much of this is not specific to the moon.

    def moon_elliptic_to_polar(a, e, M):
        return Cartesian2d(*cartesian2d_in_plane_of_orbit(
                    iterate_for_eccentric_anomaly(
                        M,
                        e,
                        eccentric_anomaly_first_approximation(M, e)
                    ),
                    e,
                    a
                )
            ).to_polar().as_tuple()

    def moon_elements_to_spherical(elems):
        return Cartesian3d(*position_from_plane_of_orbit_to_ecliptic(
                *(moon_elliptic_to_polar(elems.a, elems.e, elems.M) + (elems.Ω, elems.i, elems.ω, ))
            )).to_cylindrical().to_spherical()

    def moon_elements_to_spherical_perturbation(elems):

        return map (
            fsum,
            (lambda args : (
                moon_2_distance_perturbation_terms(*args),
                moon_5_latitude_perturbation_terms(*args),
                moon_12_longitude_perturbation_terms(*args),
            )) (moon_perturbation_arguments(
                mean_longitude,
                elems.Ω,
                elems.ω,
                elems.M,
                mean_anomaly,
            )),
        )

    def date_and_sun_mean_to_moon_ecliptic_1(year, month, day):
        return map(sum, zip(
            *(lambda elems: (
                moon_elements_to_spherical(elems).as_tuple(),
                moon_elements_to_spherical_perturbation(elems),
            )) (moon_elements(day_number(year, month, day)))
        ))

    (pdist, plat, plon) = date_and_sun_mean_to_moon_ecliptic_1(year, month, day)
    return (pdist, plat, plon)
