#!/usr/bin/python3
# -*- coding: utf-8 -*-
# pylint: disable=non-ascii-name
# pylint: disable=missing-class-docstring, missing-function-docstring, missing-module-docstring

# Aaron Mansheim 2011-06-22
# Based on: http://stjarnhimlen.se/comp/tutorial.html 2009-07-02
# Note that sunriset.c by the same author is in the public domain.

import unittest
from math import fsum
from sunpos import Degrees, Geographic, Celestial, Horizontal
from sunpos import Cartesian2d, Cartesian3d, Cylindrical, Polar, Spherical
from sunpos import Arcdegrees, hour_angle, toGMST0
from sunpos import Seconds, Date
from sunpos import to_sidereal_time
from sunpos import cartesian2d_in_plane_of_orbit
from sunpos import ecliptic_anomaly_to_longitude
from sunpos import eccentric_anomaly_first_approximation
from sunpos import ecliptic_to_celestial
from sunpos import iterate_for_eccentric_anomaly
from sunpos import obliquity_of_the_ecliptic
from sunpos import position_from_plane_of_orbit_to_ecliptic

from sunpos import OrbitalElements
from sunpos import date_and_sun_mean_to_moon_ecliptic
from sunpos import moon_2_distance_perturbation_terms
from sunpos import moon_5_latitude_perturbation_terms
from sunpos import moon_12_longitude_perturbation_terms
from sunpos import moon_elements
from sunpos import moon_perturbation_arguments

from sunpos import date_to_sun_earth_ecliptic
from sunpos import sun_earth_ecliptic_to_celestial
from sunpos import sun_earth_elements
from sunpos import time_and_location_to_sun_horizontal

class TestRev(unittest.TestCase):

#   # Uncomment to show where the error limits can or can't be tightened.
#   def assertAlmostEqual(self, *args, **kwds):
#       places = kwds['places']
#       kwds['places'] = places + 1
#       super(TestSunPos, self).assertAlmostEqual(*args, **kwds)

    def test_rev(self):
        angle = Degrees(-442.00000000001827).rev()
            # yes, this function has double precision
        self.assertAlmostEqual(angle, 277.999999999982, places = 10)

class TestEclipticToEquatorial(unittest.TestCase):

    # The next few tests form a procedure to convert
    # the Sun's spherical ecliptic coordinates
    # to its spherical equatorial coordinates.

    # The ecliptic sphere is centered at the Sun
    # with the following longitudinal angles in degrees:
    #      0 = vernal equinox,
    #     90 = summer solstice,
    #    180 = autumnal equinox,
    #    270 = winter solstice.

    # The equatorial sphere is the one that appears to
    # hold the stars in place. Its equator is the zodiac,
    # above Earth's equator. Angles around the zodiac are
    # called "right ascension". Angles from the zodiac
    # toward the celestial poles (the North Star and
    # the Southern Cross) are called "declination".

    summer_solstice_geographic = Geographic(
            latitude=Degrees(0),
            longitude=Degrees(90)) # because it's the solstice
    summer_solstice_r = 1.0 # variations are small, not relevant

    summer_solstice_cylindrical = Cylindrical(
            cylindrical_radius = summer_solstice_r,
            h = 0.000000,
            longitude = summer_solstice_geographic.longitude)

    def test_spherical_to_cylindrical(self):
        # Change the ecliptic spherical coordinates
        # to ecliptic cylindrical coordinates.
        summer_solstice = Spherical(
                self.summer_solstice_r,
                self.summer_solstice_geographic.latitude,
                self.summer_solstice_geographic.longitude
            ).cylindrical()
        self.assertAlmostEqual(
            summer_solstice.cylindrical_radius,
            self.summer_solstice_r,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.h,
            self.summer_solstice_cylindrical.h,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.longitude,
            self.summer_solstice_geographic.longitude,
            places = 6,
        )

    summer_solstice_cartesian3d = Cartesian3d(
            x=0.000000,
            y=summer_solstice_r,
            z=summer_solstice_cylindrical.h)

    def test_cylindrical_to_cartesian3d(self):
        # Change the ecliptic cylindrical coordinates
        # to ecliptic cartesian (x-y-z) coordinates.
        summer_solstice = self.summer_solstice_cylindrical.cartesian3d()
        self.assertAlmostEqual(
            summer_solstice.x,
            self.summer_solstice_cartesian3d.x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.y,
            self.summer_solstice_cartesian3d.y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.z,
            self.summer_solstice_cartesian3d.z,
            places = 6,
        )

    def test_spherical_to_cartesian3d(self):
        # The two previous steps can be combined in
        # one operation.
        summer_solstice = Spherical(
            self.summer_solstice_r,
            self.summer_solstice_geographic.latitude,
            self.summer_solstice_geographic.longitude,
        ).cartesian3d()
        self.assertAlmostEqual(
            summer_solstice.x,
            self.summer_solstice_cartesian3d.x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.y,
            self.summer_solstice_cartesian3d.y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice.z,
            self.summer_solstice_cartesian3d.z,
            places = 6,
        )

    obliquity_of_the_ecliptic_low_accuracy = 23.4
        # degrees, low accuracy only for the sake of this example

    summer_solstice_celestial_x = summer_solstice_cartesian3d.x
    summer_solstice_celestial_y = 0.917755
    summer_solstice_celestial_z = 0.397148

    def test_rotate_cartesian3d_about_x(self):
        # Rotate the ecliptic cartesian coordinates
        # in such a way that the x-axis does not move.
        # Meanwhile the z-axis moves 23.4 degrees:
        # from the direction out of the ecliptic
        # (the plane of Earth's revolution around the Sun)
        # to the direction of Earth's north pole.
        # This angle of 23.4 degrees is sometimes
        # called "the obliquity of the ecliptic".
        # The rotated coordinates are equatorial
        # cartesian coordinates.
        summer_solstice_celestial = self.summer_solstice_cartesian3d.rotate_about_x(
                self.obliquity_of_the_ecliptic_low_accuracy)
        self.assertAlmostEqual(
            summer_solstice_celestial.x,
            self.summer_solstice_celestial_x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial.y,
            self.summer_solstice_celestial_y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial.z,
            self.summer_solstice_celestial_z,
            places = 6,
        )

    def test_cartesian3d_to_spherical(self):
        # Change from the (equatorial) cartesian coordinates
        # back to (equatorial) spherical coordinates.
        # By looking up those coordinates in a star atlas,
        # we would see which stars are blotted out by the Sun
        # on the particular day that we used in our calculations.
        summer_solstice_celestial = Cartesian3d(
            self.summer_solstice_celestial_x,
            self.summer_solstice_celestial_y,
            self.summer_solstice_celestial_z,
        ).spherical()
        self.assertAlmostEqual(
            summer_solstice_celestial.radius,
            self.summer_solstice_r,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial.latitude,
            self.obliquity_of_the_ecliptic_low_accuracy,
            places = 5,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial.longitude,
            self.summer_solstice_geographic.longitude,
            places = 6,
        )

class TestPredictEcliptic(unittest.TestCase):
    summer_solstice_r = 1.0 # variations are small, not relevant

    date = Date(1990, 4, 19)

    d_1990_04_19 = -3543

    def test_day_number(self):
        # When tracking heavenly motions it's more convenient
        # to number days consecutively than to use a
        # conventional calendar.
        day_number = self.date.day_number()
        self.assertEqual(day_number, self.d_1990_04_19)

    longitude_at_perihelion = 282.773548
    eccentricity = 0.016713
    Ms = Degrees(104.065284)

    def test_sun_earth_elements(self):
        # Various parameters of an orbit
        # can be computed for any particular
        # day (or fraction of a day).
        # Somewhere else we can draw a diagram to
        # explain what they are. The point
        # is that they are helpful in accurately
        # tracking orbits.
        (
            periapsis,
            semimajor_axis,  # mean distance
            eccentricity,
            mean_anomaly,
        ) = sun_earth_elements(self.d_1990_04_19)
        self.assertAlmostEqual(periapsis, self.longitude_at_perihelion, places = 6)
        self.assertAlmostEqual(semimajor_axis, self.summer_solstice_r, places = 6)
        self.assertAlmostEqual(eccentricity, self.eccentricity, places = 6)
        self.assertAlmostEqual(mean_anomaly, self.Ms, places = 6)

    def test_rev_2(self):
        # Within the calculation of the mean anomaly M,
        # it is necessary to remove from the computed
        # value as many complete 360-degree revolutions
        # as we can.
        self.assertAlmostEqual(
            Degrees(-3135.934716).rev(),
            self.Ms,
            places = 6,
        )

    Ls = Degrees(26.838832)

    def test_mean_anomaly_to_longitude(self):
        # Mean longitude is one parameter that is
        # derived from two of the more basic ones
        # above.
        mean_longitude = ecliptic_anomaly_to_longitude(
            self.longitude_at_perihelion,
            self.Ms,
        )
        self.assertAlmostEqual(
            mean_longitude,
            self.Ls,
            places = 6,
        )

    obliquity_of_the_ecliptic = 23.440562

    def test_obliquity_of_the_ecliptic(self):
        # When we examine the obliquity of the ecliptic
        # more accurately, we see it drift over time.
        obliquity = obliquity_of_the_ecliptic(self.d_1990_04_19)
        self.assertAlmostEqual(
            obliquity,
            self.obliquity_of_the_ecliptic,
            places = 6,
        )

    eccentric_anomaly = Degrees(104.9904)

    def test_eccentric_anomaly_first_approximation(self):
        # There is no outstanding method to solve
        # Kepler's equation for eccentric anomaly.
        # We'll start with a reasonably accurate value
        # and get a better one later.
        eccentricity = eccentric_anomaly_first_approximation(
            self.Ms,
            self.eccentricity,
        )
        self.assertAlmostEqual(eccentricity, self.eccentric_anomaly, places = 4)

    x = -0.275370
    y = +0.965834

    def test_eccentric_to_cartesian2d(self):
        point = Cartesian2d.from_eccentric(
                self.eccentric_anomaly,
                self.eccentricity)
        self.assertAlmostEqual(point.x, self.x, places = 6)
        self.assertAlmostEqual(point.y, self.y, places = 6)

    distance = 1.004323
    ν = 105.9134

    def test_cartesian2d_to_polar(self):
        polar = Cartesian2d(self.x, self.y).polar()
        self.assertAlmostEqual(polar.r, self.distance, places = 6)
        self.assertAlmostEqual(polar.θ, self.ν, places = 4)

    true_longitude = Degrees(28.6869)
    v_ref_arcminutes = Arcdegrees(28.6813).arcminutes()
    # Demonstrating agreement within 1 arcminute (actually 1/2 arcminute):
    # 1 arcdegree = 60 arcminutes.
    within_half_arcminute = 0

    def test_true_anomaly_to_longitude(self):
        true_longitude = ecliptic_anomaly_to_longitude(
            self.longitude_at_perihelion,
            self.ν,
        )
        self.assertAlmostEqual(
            true_longitude,
            self.true_longitude,
            places = 4,
        )

    def test_sun_earth_ecliptic_coordinates_vs_astronomical_almanac(self):
        # Supposing that other tests have passed,
        # these are the values that we computed.
        # Here we're not really testing anything.
        # Instead we're making an assertion that the
        # computed values are close to the reference values.
        v_arcminutes = Arcdegrees(self.true_longitude).arcminutes()
        # Demonstrating agreement within 1 arcminute (actually 1/2 arcminute):
        # 1 arcdegree = 60 arcminutes.
        self.assertAlmostEqual(
                v_arcminutes, self.v_ref_arcminutes, places = self.within_half_arcminute)
        self.assertAlmostEqual(self.distance, 1.004311, places = 4)

    def test_date_to_sun_earth_ecliptic(self):
        """
        Going through all of the above in one swift move.
        """
        (dist, lon, mlon, obliquity) = date_to_sun_earth_ecliptic(self.date)
        v_arcminutes = Arcdegrees(lon).arcminutes()
        self.assertAlmostEqual(
                v_arcminutes, self.v_ref_arcminutes, places = self.within_half_arcminute)
        self.assertAlmostEqual(dist, 1.004311, places = 4)
        self.assertAlmostEqual(mlon, self.Ls, places = 6)
        self.assertAlmostEqual(
            obliquity,
            self.obliquity_of_the_ecliptic,
            places = 6,
        )

class TestEclipticToCelestial(unittest.TestCase):
    obliquity_of_the_ecliptic = 23.440562

    true_longitude = Degrees(28.6869)
    # Demonstrating agreement within 1 arcminute (actually 1/2 arcminute):
    # 1 arcdegree = 60 arcminutes.
    within_half_arcminute = 0
    distance = 1.004323

    # These x-y coordinates point the positive x-axis
    # in the conventional direction, toward the vernal point.
    # The previous ones had the positive x-axis toward perihelion.
    ecliptic = Cartesian3d(
            0.881048,
            0.482099,
            0.0) # the Sun and Earth stay in plane, close enough

    def test_polar_to_cartesian2d(self):
        ecliptic = Polar(self.distance, self.true_longitude).cartesian2d()
        self.assertAlmostEqual(ecliptic.x, self.ecliptic.x, places = 6)
        self.assertAlmostEqual(ecliptic.y, self.ecliptic.y, places = 5)

    equatorial = Cartesian3d(ecliptic.x, 0.442312, 0.191778)

    def test_rotate_cartesian3d_about_x_2(self):
        equatorial = Cartesian3d(
            self.ecliptic.x,
            self.ecliptic.y,
            self.ecliptic.z,
        ).rotate_about_x(self.obliquity_of_the_ecliptic)
        self.assertAlmostEqual(equatorial.x, self.equatorial.x, places = 6)
        self.assertAlmostEqual(equatorial.y, self.equatorial.y, places = 5)
        self.assertAlmostEqual(equatorial.z, self.equatorial.z, places = 6)

    celestial = Celestial(RA=Degrees(26.658078), Decl=Degrees(11.008375))

    def test_cartesian3d_to_spherical_2(self):
        point = Cartesian3d(
            self.equatorial.x,
            self.equatorial.y,
            self.equatorial.z,
        ).spherical()
        self.assertAlmostEqual(point.radius, self.distance, places = 6)
        self.assertAlmostEqual(point.latitude, self.celestial.Decl, places = 5)
        self.assertAlmostEqual(point.longitude, self.celestial.RA, places = 3)

    def test_arcdegrees_to_arcminutes(self):
        minutes = Arcdegrees(self.celestial.Decl).arcminutes() # 660.5025
        astronomical_almanac_m = ( # 11°0'22" in arcminutes
            Arcdegrees(11).arcminutes() # 11° = 660'
            + 0 # 0' = 0'
            + Seconds(22).minutes()) # 22" = 0.366667'
        self.assertAlmostEqual(minutes, astronomical_almanac_m, places = self.within_half_arcminute)

    def test_arcdegrees_to_hours(self):
        hours = Arcdegrees(self.celestial.RA).hours()
        seconds = Seconds.fromhours(hours)
        astronomical_almanac_s = (
            Seconds.fromhours(1)
            + Seconds.fromminutes(46)
            + Seconds(36.0))
        self.assertAlmostEqual(seconds, astronomical_almanac_s, places = -1)

    def test_sun_earth_ecliptic_to_celestial(self):
        celestial = sun_earth_ecliptic_to_celestial(
            self.distance,
            self.true_longitude,
            self.obliquity_of_the_ecliptic,
        )
        # self.assertAlmostEqual(distance1, self.distance, places = 6)
        self.assertAlmostEqual(celestial.Decl, self.celestial.Decl, places = 4)
        self.assertAlmostEqual(
                celestial.RA, self.celestial.RA, places = 3)

class TestCelestialToHorizontal(unittest.TestCase):

    GMST0 = 13.78925

    def test_gmst0(self):
        gmst0 = toGMST0(self.Ls)
        self.assertAlmostEqual(gmst0, self.GMST0, places = 4)

    terrestrial = Geographic(Degrees(+60), Degrees(+15.0))
    hours_UT = 0.0
    sidtime0 = GMST0 + 1

    def test_sidereal_time(self):
        sidtime0 = to_sidereal_time(
            self.GMST0,
            self.hours_UT,
            self.terrestrial.longitude,
        )
        self.assertAlmostEqual(sidtime0, self.sidtime0, places = 4)

    hour_angle = 13.01205
    hour_angle_degrees = Degrees(195.1808)

    summer_solstice_r = 1.0 # variations are small, not relevant
    Ls = Degrees(26.838832)

    celestial = Celestial(RA=Degrees(26.658078), Decl=Degrees(11.008375))

    celestial_horizontal_x = -0.947346
    celestial_horizontal_y = -0.257047
    celestial_horizontal_z = +0.190952

    def test_hour_angle(self):
        hours = Arcdegrees(self.celestial.RA).hours()
        angle = hour_angle(self.sidtime0, hours)
        degrees = Arcdegrees.fromhours(angle)
        self.assertAlmostEqual(angle, self.hour_angle, places = 4)
        self.assertAlmostEqual(degrees, self.hour_angle_degrees, places = 3)

    def test_spherical_to_cartesian3d_2(self):
        # Note that positive z is toward the celestial south pole!
        # This is the consequence of the right-hand rule
        # for our choice of +x toward the intersection of the
        # zodaic meridian with the southern half of the
        # north-zenith-south meridian, and +y toward the western
        # point of the horizon on the west-zenith-east meridian.
        point = Spherical(
            self.summer_solstice_r,
            self.celestial.Decl,
            self.hour_angle_degrees,
        ).cartesian3d()
        self.assertAlmostEqual(point.x, self.celestial_horizontal_x, places = 6)
        self.assertAlmostEqual(point.y, self.celestial_horizontal_y, places = 6)
        self.assertAlmostEqual(point.z, self.celestial_horizontal_z, places = 6)

    sun_alt_az_x = -0.915902
    sun_alt_az_y = -0.257047 # = celestial_horizontal_y
    sun_alt_az_z = -0.308303

    def test_decline_cartesian3d_about_y(self):
        point = Cartesian3d(
            self.celestial_horizontal_x,
            self.celestial_horizontal_y,
            self.celestial_horizontal_z,
        ).decline_about_y(self.terrestrial.latitude)
        self.assertAlmostEqual(point.x, self.sun_alt_az_x, places = 6)
        self.assertAlmostEqual(point.y, self.sun_alt_az_y, places = 6)
        self.assertAlmostEqual(point.z, self.sun_alt_az_z, places = 5)

    azimuth = 15.676697
    altitude = 342.042994

    def test_cartesian3d_to_spherical_3(self):
        point = Cartesian3d(
            self.sun_alt_az_x,
            self.sun_alt_az_y,
            self.sun_alt_az_z,
        ).spherical()
        azimuth = Degrees(point.longitude + 180).rev()
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(point.latitude, self.altitude, places = 4)

    def test_sun_earth_celestial_to_horizontal(self):
        (altitude, azimuth) = Horizontal.from_celestial(
            self.Ls,
            self.celestial,
            self.hours_UT,
            self.terrestrial,
        ).tuple()
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(altitude, self.altitude, places = 4)

class TestPredictHorizontal(unittest.TestCase):
    date = Date(1990, 4, 19)
    hours_UT = 0.0
    terrestrial = Geographic(Degrees(+60), Degrees(+15.0))

    azimuth = 15.676697
    altitude = 342.042994

    def test_time_and_location_to_sun_horizontal(self):
        (
            altitude,
            azimuth,
        ) = time_and_location_to_sun_horizontal(
            self.date,
            self.hours_UT, self.terrestrial,
        )
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(altitude, self.altitude, places = 4)

class TestMoonPerturbation(unittest.TestCase):
    d_1990_04_19 = -3543
    Ls = Degrees(26.838832)
    Ms = Degrees(104.065284)

    moon_elements = OrbitalElements(
            Degrees(312.7381),  # Ω
            Degrees(5.1454),    # i
            Degrees(95.7454),   # ω
            60.2666,    # a, Earth equatorial radii
            0.054900,   # e
            Degrees(266.0954))  # M

    def test_moon_elements(self):
        elems = moon_elements(self.d_1990_04_19)
        self.assertAlmostEqual(elems.Ω, self.moon_elements.Ω, places = 4)
        self.assertAlmostEqual(elems.i, self.moon_elements.i, places = 4)
        self.assertAlmostEqual(elems.ω, self.moon_elements.ω, places = 4)
        self.assertAlmostEqual(elems.a, self.moon_elements.a, places = 4)
        self.assertAlmostEqual(elems.e, self.moon_elements.e, places = 6)
        self.assertAlmostEqual(elems.M, self.moon_elements.M, places = 4)

    moon_E0 = Degrees(262.9689)

    def test_eccentric_anomaly_first_approximation_2(self):
        eccentric_anomaly = eccentric_anomaly_first_approximation(
                self.moon_elements.M, self.moon_elements.e)
        self.assertAlmostEqual(eccentric_anomaly, self.moon_E0, places = 4)

    moon_E = Degrees(262.9735)

    def test_iterate_for_eccentric_anomaly(self):
        eccentric_anomaly = iterate_for_eccentric_anomaly(
            self.moon_elements.M,
            self.moon_elements.e,
            self.moon_E0,
        )
        self.assertAlmostEqual(eccentric_anomaly, self.moon_E, places = 4)

    moon_x = -10.68095
    moon_y = -59.72377

    def test_cartesian2d_in_plane_of_orbit(self):
        cartesian = cartesian2d_in_plane_of_orbit(
            self.moon_E,
            self.moon_elements.e,
            self.moon_elements.a,
        )
        self.assertAlmostEqual(cartesian.x, self.moon_x, places = 5)
        self.assertAlmostEqual(cartesian.y, self.moon_y, places = 5)

    moon_polar = Polar(r=60.67134, θ=Degrees(259.8605))

    def test_cartesian2d_to_polar_2(self):
        polar = Cartesian2d(self.moon_x, self.moon_y).polar()
        self.assertAlmostEqual(polar.r, self.moon_polar.r, places = 5)
        self.assertAlmostEqual(polar.θ, self.moon_polar.θ, places = 4)

    moon_xeclip = +37.65311
    moon_yeclip = -47.57180
    moon_zeclip = -0.41687

    def test_position_from_plane_of_orbit_to_ecliptic(self):
        (xeclip, yeclip, zeclip) = position_from_plane_of_orbit_to_ecliptic(
            self.moon_polar,
            self.moon_elements.Ω,
            self.moon_elements.i,
            self.moon_elements.ω,
        ).tuple()
        self.assertAlmostEqual(xeclip, self.moon_xeclip, places = 4)
        self.assertAlmostEqual(yeclip, self.moon_yeclip, places = 4)
        self.assertAlmostEqual(zeclip, self.moon_zeclip, places = 4)

    moon_r_eclip = 60.6713
    moon_lat_eclip = 359.6063
    moon_lon_eclip = 308.3616

    def test_cartesian3d_to_spherical_4(self):
        point = Cartesian3d(
            self.moon_xeclip,
            self.moon_yeclip,
            self.moon_zeclip,
        ).spherical()
        self.assertAlmostEqual(point.radius, self.moon_r_eclip, places = 4)
        self.assertAlmostEqual(point.latitude, self.moon_lat_eclip, places = 4)
        self.assertAlmostEqual(point.longitude, self.moon_lon_eclip, places = 4)

    # Ls = 26.8388
    Lm = Degrees(314.5789)
    # Ms = 104.0653
    Mm = moon_elements.M
    D = 287.7401
    F = 1.8408

    def test_moon_perturbation_arguments(self):
        (
                Ls, # sun mean longitude
                Lm, # moon mean longitude
                Ms, # sun mean anomaly
                Mm,  # moon mean anomaly
                D,   # moon mean elongation
                F,  # moon argument of latitude
                ) = moon_perturbation_arguments(
            self.Ls,
            self.moon_elements.Ω,
            self.moon_elements.ω,
            self.moon_elements.M,
            self.Ms,
        )
        self.assertAlmostEqual(Ls, self.Ls, places = 4)
        self.assertAlmostEqual(Lm, self.Lm, places = 4)
        self.assertAlmostEqual(Ms, self.Ms, places = 4)
        self.assertAlmostEqual(Mm, self.Mm, places = 4)
        self.assertAlmostEqual(D, self.D, places = 4)
        self.assertAlmostEqual(F, self.F, places = 4)

    moon_dlon = (
        -0.9847,
        -0.3819,
        -0.1804,
        +0.0405,
        -0.0244,
        +0.0452,
        +0.0428,
        +0.0126,
        -1 * (-0.0333),
        -0.0055,
        -0.0079,
        -0.0029,
    )
    moon_sum_dlon = -1.4132

    moon_dlat = (
        -0.0958,
        -0.0414,
        -0.0365,
        -0.0200,
        +0.0018,
    )
    moon_sum_dlat = -0.1919

    moon_ddist = (
        -0.3680,
        +0.3745 + (0.00005),  # fudge factor
    )

    moon_sum_ddist = +0.0066

    def test_moon_12_longitude_perturbation_terms(self):
        dlon = moon_12_longitude_perturbation_terms(
            self.Ls,
            self.Lm,
            self.Ms,
            self.Mm,
            self.D,
            self.F,
        )
        for pair in zip(dlon, self.moon_dlon):
            self.assertAlmostEqual(*pair, places = 4)
        sum_dlon = fsum(dlon)
        self.assertAlmostEqual(sum_dlon, self.moon_sum_dlon, places = 4)

    def test_moon_5_latitude_perturbation_terms(self):
        dlat = moon_5_latitude_perturbation_terms(
            self.Ls,
            self.Lm,
            self.Ms,
            self.Mm,
            self.D,
            self.F,
        )
        for pair in zip(dlat, self.moon_dlat):
            self.assertAlmostEqual(*pair, places = 4)
        sum_dlat = fsum(dlat)
        self.assertAlmostEqual(sum_dlat, self.moon_sum_dlat, places = 4)

    def test_moon_2_distance_perturbation_terms(self):
        ddist = moon_2_distance_perturbation_terms(
            self.Ls,
            self.Lm,
            self.Ms,
            self.Mm,
            self.D,
            self.F,
        )
        for pair in zip(ddist, self.moon_ddist):
            self.assertAlmostEqual(*pair, places = 4)
        sum_ddist = fsum(ddist)
        self.assertAlmostEqual(sum_ddist, self.moon_sum_ddist, places = 4)

class TestMoonEcliptic(unittest.TestCase):
    date = Date(1990, 4, 19)
    Ls = Degrees(26.838832)
    Ms = Degrees(104.065284)

    moon_r_eclip_perturbed = 60.6779
    moon_lat_eclip_perturbed = Degrees(359.4144)
    moon_lon_eclip_perturbed = Degrees(306.9484)

    moon_r_eclip_astro_almanac = 60.793
    moon_lat_eclip_astro_almanac = 359.45
    moon_lon_eclip_astro_almanac = 306.94

    def test_date_and_sun_mean_to_moon_ecliptic(self):
        (radius, lat, lon) = date_and_sun_mean_to_moon_ecliptic(
            self.date,
            self.Ls,
            self.Ms,
        )
        self.assertAlmostEqual(radius, self.moon_r_eclip_perturbed, places = 4)
        self.assertAlmostEqual(lat, self.moon_lat_eclip_perturbed, places = 4)
        self.assertAlmostEqual(lon, self.moon_lon_eclip_perturbed, places = 4)

        self.assertAlmostEqual(
            radius,
            self.moon_r_eclip_astro_almanac,
            places = 0,
        )
        self.assertAlmostEqual(
            lat,
            self.moon_lat_eclip_astro_almanac,
            places = 1,
        )
        self.assertAlmostEqual(
            lon,
            self.moon_lon_eclip_astro_almanac,
            places = 1,
        )

class TestMoonCelestial(unittest.TestCase):
    obliquity_of_the_ecliptic = 23.440562
    GMST0 = 13.78925
    sidtime0 = GMST0 + 1
    summer_solstice_r = 1.0 # variations are small, not relevant

    moon_lat_eclip_perturbed = Degrees(359.4144)
    moon_lon_eclip_perturbed = Degrees(306.9484)

    moon_Decl = 340.8968
    moon_RA = 309.5011

    moon_Decl_astro_almanac = 340.9159
    moon_RA_astro_almanac = 309.4881

    def test_ecliptic_to_celestial(self):
        (_distance1, Decl, RA) = ecliptic_to_celestial(
            self.summer_solstice_r,
            self.moon_lat_eclip_perturbed,
            self.moon_lon_eclip_perturbed,
            self.obliquity_of_the_ecliptic,
        )
        self.assertAlmostEqual(Decl, self.moon_Decl, places = 4)
        self.assertAlmostEqual(RA, self.moon_RA, places = 4)

        self.assertAlmostEqual(Decl, self.moon_Decl_astro_almanac, places = 1)
        self.assertAlmostEqual(RA, self.moon_RA_astro_almanac, places = 1)

    moon_hour_angle_degrees = 272.3377

    def test_moon_hour_angle(self):
        hours = Arcdegrees(self.moon_RA).hours()
        angle = hour_angle(self.sidtime0, hours)
        degrees = Arcdegrees.fromhours(angle)
        self.assertAlmostEqual(degrees, self.moon_hour_angle_degrees, places = 3)


if __name__ == '__main__':
    unittest.main()
