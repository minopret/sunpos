#!/usr/bin/python

# Aaron Mansheim 2011-06-22
# Based on: http://stjarnhimlen.se/comp/tutorial.html 2009-07-02
# Note that sunriset.c by the same author is in the public domain.

from sunpos import *
import unittest
from math import fabs, floor, fsum


class TestSunPos(unittest.TestCase):

#   # Uncomment to show where the error limits can or can't be tightened.
#   def assertAlmostEqual(self, *args, **kwds):
#       places = kwds['places']
#       kwds['places'] = places + 1
#       super(TestSunPos, self).assertAlmostEqual(*args, **kwds)

    def test_rev(self):
        r = rev(-442.00000000001827)
            # yes, this function has double precision
        self.assertAlmostEqual(r, 277.999999999982, places = 10)

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

    summer_solstice_long = 90
    summer_solstice_lat = 0  # because it's the solstice
    summer_solstice_r = 1.0  # variations are small, not relevant

    summer_solstice_cylinder_radius = 1.000000
    summer_solstice_cylinder_height = 0.000000
    summer_solstice_cylinder_longitude = 90.000000

    def test_spherical_to_cylindrical(self):
        # Change the ecliptic spherical coordinates
        # to ecliptic cylindrical coordinates.
        (
            summer_solstice_cylinder_radius,
            summer_solstice_cylinder_height,
            summer_solstice_cylinder_longitude,
        ) = spherical_to_cylindrical(
            self.summer_solstice_r,
            self.summer_solstice_lat,
            self.summer_solstice_long,
        )
        self.assertAlmostEqual(
            summer_solstice_cylinder_radius,
            1.000000,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_cylinder_height,
            0.000000,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_cylinder_longitude,
            90.000000,
            places = 6,
        )

    summer_solstice_x = 0.000000
    summer_solstice_y = 1.000000
    summer_solstice_z = 0.000000

    def test_cylindrical_to_cartesian3d(self):
        # Change the ecliptic cylindrical coordinates
        # to ecliptic cartesian (x-y-z) coordinates.
        (
            summer_solstice_x,
            summer_solstice_y,
            summer_solstice_z,
        ) = cylindrical_to_cartesian3d(
            self.summer_solstice_cylinder_radius,
            self.summer_solstice_cylinder_height,
            self.summer_solstice_cylinder_longitude,
        )
        self.assertAlmostEqual(
            summer_solstice_x,
            self.summer_solstice_x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_y,
            self.summer_solstice_y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_z,
            self.summer_solstice_z,
            places = 6,
        )

    def test_spherical_to_cartesian3d(self):
        # The two previous steps can be combined in
        # one operation.
        (
            summer_solstice_x,
            summer_solstice_y,
            summer_solstice_z,
        ) = spherical_to_cartesian3d(
            self.summer_solstice_r,
            self.summer_solstice_lat,
            self.summer_solstice_long,
        )
        self.assertAlmostEqual(
            summer_solstice_x,
            self.summer_solstice_x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_y,
            self.summer_solstice_y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_z,
            self.summer_solstice_z,
            places = 6,
        )

    oblecl = 23.4
        # degrees, low accuracy only for the sake of this example

    summer_solstice_celestial_x = 0.000000
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
        (
            summer_solstice_celestial_x,
            summer_solstice_celestial_y,
            summer_solstice_celestial_z,
        ) = rotate_cartesian3d_about_x(
            self.summer_solstice_x,
            self.summer_solstice_y,
            self.summer_solstice_z,
            self.oblecl,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial_x,
            self.summer_solstice_celestial_x,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial_y,
            self.summer_solstice_celestial_y,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial_z,
            self.summer_solstice_celestial_z,
            places = 6,
        )

    def test_cartesian3d_to_spherical(self):
        # Change from the (equatorial) cartesian coordinates
        # back to (equatorial) spherical coordinates.
        # By looking up those coordinates in a star atlas,
        # we would see which stars are blotted out by the Sun
        # on the particular day that we used in our calculations.
        (
            summer_solstice_celestial_r,
            summer_solstice_Decl,
            summer_solstice_RA,
        ) = cartesian3d_to_spherical(
            self.summer_solstice_celestial_x,
            self.summer_solstice_celestial_y,
            self.summer_solstice_celestial_z,
        )
        self.assertAlmostEqual(
            summer_solstice_celestial_r,
            1.000000,
            places = 6,
        )
        self.assertAlmostEqual(
            summer_solstice_Decl,
            23.400000,
            places = 5,
        )
        self.assertAlmostEqual(
            summer_solstice_RA,
            90.000000,
            places = 6,
        )

    Y = 1990
    M = 4
    D = 19

    d_1990_04_19 = -3543

    def test_day_number(self):
        # When tracking heavenly motions it's more convenient
        # to number days consecutively than to use a
        # conventional calendar.
        d = day_number(self.Y, self.M, self.D)
        self.assertEqual(d, self.d_1990_04_19)

    longitude_at_perihelion = 282.773548
    eccentricity = 0.016713
    mean_anomaly = 104.065284

    def test_sun_earth_elements(self):
        # Various parameters of an orbit
        # can be computed for any particular
        # day (or fraction of a day).
        # Somewhere else we can draw a diagram to
        # explain what they are. The point
        # is that they are helpful in accurately
        # tracking orbits.
        (
            w,
            a,  # mean distance
            e,
            M,
        ) = sun_earth_elements(self.d_1990_04_19)
        self.assertAlmostEqual(w, self.longitude_at_perihelion, places = 6)
        self.assertAlmostEqual(a, 1.000000, places = 6)
        self.assertAlmostEqual(e, self.eccentricity, places = 6)
        self.assertAlmostEqual(M, self.mean_anomaly, places = 6)

    def test_rev_2(self):
        # Within the calculation of the mean anomaly M,
        # it is necessary to remove from the computed
        # value as many complete 360-degree revolutions
        # as we can.
        self.assertAlmostEqual(
            rev(-3135.934716),
            self.mean_anomaly,
            places = 6,
        )

    mean_longitude = 26.838832

    def test_mean_anomaly_to_longitude(self):
        # Mean longitude is one parameter that is
        # derived from two of the more basic ones
        # above.
        mean_longitude = ecliptic_anomaly_to_longitude(
            self.longitude_at_perihelion,
            self.mean_anomaly,
        )
        self.assertAlmostEqual(
            mean_longitude,
            self.mean_longitude,
            places = 6,
        )

    obliquity_of_the_ecliptic = 23.440562

    def test_obliquity_of_the_ecliptic(self):
        # When we examine the obliquity of the ecliptic
        # more accurately, we see it drift over time.
        oblecl = obliquity_of_the_ecliptic(self.d_1990_04_19)
        self.assertAlmostEqual(
            oblecl,
            self.obliquity_of_the_ecliptic,
            places = 6,
        )

    eccentric_anomaly = 104.9904

    def test_eccentric_anomaly_first_approximation(self):
        # There is no outstanding method to solve
        # Kepler's equation for eccentric anomaly.
        # We'll start with a reasonably accurate value
        # and get a better one later.
        e = eccentric_anomaly_first_approximation(
            self.mean_anomaly,
            self.eccentricity,
        )
        self.assertAlmostEqual(e, self.eccentric_anomaly, places = 4)

    x = -0.275370
    y = +0.965834

    def test_eccentric_to_cartesian2d(self):
        (x, y) = eccentric_to_cartesian2d(
            self.eccentric_anomaly,
            self.eccentricity,
        )
        self.assertAlmostEqual(x, self.x, places = 6)
        self.assertAlmostEqual(y, self.y, places = 6)

    distance = 1.004323
    true_anomaly = 105.9134

    def test_cartesian2d_to_polar(self):
        (distance, true_anomaly) = cartesian2d_to_polar(self.x, self.y)
        self.assertAlmostEqual(distance, self.distance, places = 6)
        self.assertAlmostEqual(true_anomaly, self.true_anomaly, places = 4)

    true_longitude = 28.6869

    def test_true_anomaly_to_longitude(self):
        true_longitude = ecliptic_anomaly_to_longitude(
            self.longitude_at_perihelion,
            self.true_anomaly,
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
        v_arcminutes = arcdegrees_to_arcminutes(self.true_longitude)
        v_ref_arcminutes = arcdegrees_to_arcminutes(28.6813)
        # Demonstrating agreement within 1 arcminute (actually 1/2 arcminute):
        # 1 arcdegree = 60 arcminutes.
        self.assertAlmostEqual(v_arcminutes, v_ref_arcminutes, places = 0)
        self.assertAlmostEqual(self.distance, 1.004311, places = 4)

    def test_date_to_sun_earth_ecliptic(self):
        """
        Going through all of the above in one swift move.
        """
        (dist, lon, mlon, oblecl) = date_to_sun_earth_ecliptic(
            self.Y,
            self.M,
            self.D,
        )
        v_arcminutes = arcdegrees_to_arcminutes(lon)
        v_ref_arcminutes = arcdegrees_to_arcminutes(28.6813)
        # Demonstrating agreement within 1 arcminute (actually 1/2 arcminute):
        # 1 arcdegree = 60 arcminutes.
        self.assertAlmostEqual(v_arcminutes, v_ref_arcminutes, places = 0)
        self.assertAlmostEqual(dist, 1.004311, places = 4)
        self.assertAlmostEqual(mlon, self.mean_longitude, places = 6)
        self.assertAlmostEqual(
            oblecl,
            self.obliquity_of_the_ecliptic,
            places = 6,
        )

    # These x-y coordinates point the positive x-axis
    # in the conventional direction, toward the vernal point.
    # The previous ones had the positive x-axis toward perihelion.
    ecliptic_x = 0.881048
    ecliptic_y = 0.482099

    def test_polar_to_cartesian2d(self):
        (
            ecliptic_x,
            ecliptic_y,
        ) = polar_to_cartesian2d(self.distance, self.true_longitude)
        self.assertAlmostEqual(ecliptic_x, self.ecliptic_x, places = 6)
        self.assertAlmostEqual(ecliptic_y, self.ecliptic_y, places = 5)

    ecliptic_z = 0.0  # the Sun and Earth stay in plane, close enough

    equatorial_x = 0.881048
    equatorial_y = 0.442312
    equatorial_z = 0.191778

    def test_rotate_cartesian3d_about_x_2(self):
        (
            equatorial_x,
            equatorial_y,
            equatorial_z,
        ) = rotate_cartesian3d_about_x(
            self.ecliptic_x,
            self.ecliptic_y,
            self.ecliptic_z,
            self.obliquity_of_the_ecliptic,
        )
        self.assertAlmostEqual(equatorial_x, self.equatorial_x, places = 6)
        self.assertAlmostEqual(equatorial_y, self.equatorial_y, places = 5)
        self.assertAlmostEqual(equatorial_z, self.equatorial_z, places = 6)

    Decl = 11.008375
    RA = 26.658078

    def test_cartesian3d_to_spherical_2(self):
        (
            distance,
            Decl,
            RA,
        ) = cartesian3d_to_spherical(
            self.equatorial_x,
            self.equatorial_y,
            self.equatorial_z,
        )
        self.assertAlmostEqual(distance, self.distance, places = 6)
        self.assertAlmostEqual(Decl, self.Decl, places = 5)
        self.assertAlmostEqual(RA, self.RA, places = 3)

    def test_arcdegrees_to_arcminutes(self):
        m = arcdegrees_to_arcminutes(self.Decl)
        astronomical_almanac_m = (
            arcdegrees_to_arcminutes(11)
            + 0
            + seconds_to_minutes(22))
        self.assertAlmostEqual(m, astronomical_almanac_m, places = 0)

    def test_arcdegrees_to_hours(self):
        h = arcdegrees_to_hours(self.RA)
        s = hours_to_seconds(h)
        astronomical_almanac_s = (
            hours_to_seconds(1)
            + minutes_to_seconds(46)
            + 36.0)
        self.assertAlmostEqual(s, astronomical_almanac_s, places = -1)

    def test_sun_earth_ecliptic_to_celestial(self):
        (distance1, Decl, RA) = sun_earth_ecliptic_to_celestial(
            self.distance,
            self.true_longitude,
            self.obliquity_of_the_ecliptic,
        )
        self.assertAlmostEqual(distance1, self.distance, places = 6)
        self.assertAlmostEqual(Decl, self.Decl, places = 4)
        self.assertAlmostEqual(RA, self.RA, places = 3)

    GMST0 = 13.78925

    def test_GMST0(self):
        gmst0 = GMST0(self.mean_longitude)
        self.assertAlmostEqual(gmst0, self.GMST0, places = 4)

    terrestrial_longitude = +15.0
    hours_UT = 0.0
    sidtime0 = 14.78925

    def test_sidereal_time(self):
        sidtime0 = sidereal_time(
            self.GMST0,
            self.hours_UT,
            self.terrestrial_longitude,
        )
        self.assertAlmostEqual(sidtime0, self.sidtime0, places = 4)

    hour_angle = 13.01205
    hour_angle_degrees = 195.1808

    def test_hour_angle(self):
        h = arcdegrees_to_hours(self.RA)
        ha = hour_angle(self.sidtime0, h)
        d = hours_to_arcdegrees(ha)
        self.assertAlmostEqual(ha, self.hour_angle, places = 4)
        self.assertAlmostEqual(d, self.hour_angle_degrees, places = 3)

    celestial_horizontal_x = -0.947346
    celestial_horizontal_y = -0.257047
    celestial_horizontal_z = +0.190952

    def test_spherical_to_cartesian3d_2(self):
        # Note that positive z is toward the celestial south pole!
        # This is the consequence of the right-hand rule
        # for our choice of +x toward the intersection of the
        # zodaic meridian with the southern half of the
        # north-zenith-south meridian, and +y toward the western
        # point of the horizon on the west-zenith-east meridian.
        (
            x,
            y,
            z,
        ) = spherical_to_cartesian3d(
            1.0,
            self.Decl,
            self.hour_angle_degrees,
        )
        self.assertAlmostEqual(x, self.celestial_horizontal_x, places = 6)
        self.assertAlmostEqual(y, self.celestial_horizontal_y, places = 6)
        self.assertAlmostEqual(z, self.celestial_horizontal_z, places = 6)

    terrestrial_latitude = +60
    sun_alt_az_x = -0.915902
    sun_alt_az_y = -0.257047
    sun_alt_az_z = -0.308303

    def test_decline_cartesian3d_about_y(self):
        (
            x,
            y,
            z,
        ) = decline_cartesian3d_about_y(
            self.celestial_horizontal_x,
            self.celestial_horizontal_y,
            self.celestial_horizontal_z,
            self.terrestrial_latitude,
        )
        self.assertAlmostEqual(x, self.sun_alt_az_x, places = 6)
        self.assertAlmostEqual(y, self.sun_alt_az_y, places = 6)
        self.assertAlmostEqual(z, self.sun_alt_az_z, places = 5)

    azimuth = 15.676697
    altitude = 342.042994

    def test_cartesian3d_to_spherical_3(self):
        (
            distance,
            altitude,
            azimuth,
        ) = cartesian3d_to_spherical(
            self.sun_alt_az_x,
            self.sun_alt_az_y,
            self.sun_alt_az_z,
        )
        azimuth = rev(azimuth + 180)
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(altitude, self.altitude, places = 4)

    def test_sun_earth_celestial_to_alt_azimuth(self):
        (altitude, azimuth) = sun_earth_celestial_to_alt_azimuth(
            self.mean_longitude,
            self.Decl,
            self.RA,
            self.hours_UT,
            self.terrestrial_latitude,
            self.terrestrial_longitude,
        )
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(altitude, self.altitude, places = 4)

    def test_time_and_location_to_sun_alt_azimuth(self):
        (
            altitude,
            azimuth,
        ) = time_and_location_to_sun_alt_azimuth(
            self.Y, self.M, self.D, self.hours_UT,
            self.terrestrial_latitude, self.terrestrial_longitude,
        )
        self.assertAlmostEqual(azimuth, self.azimuth, places = 4)
        self.assertAlmostEqual(altitude, self.altitude, places = 4)

    moon_N = 312.7381
    moon_i = 5.1454
    moon_w = 95.7454
    moon_a = 60.2666  # Earth equatorial radii
    moon_e = 0.054900
    moon_M = 266.0954

    def test_moon_elements(self):
        (N, i, w, a, e, M) = moon_elements(self.d_1990_04_19)
        self.assertAlmostEqual(N, self.moon_N, places = 4)
        self.assertAlmostEqual(i, self.moon_i, places = 4)
        self.assertAlmostEqual(w, self.moon_w, places = 4)
        self.assertAlmostEqual(a, self.moon_a, places = 4)
        self.assertAlmostEqual(e, self.moon_e, places = 6)
        self.assertAlmostEqual(M, self.moon_M, places = 4)

    moon_E0 = 262.9689

    def test_eccentric_anomaly_first_approximation_2(self):
        E0 = eccentric_anomaly_first_approximation(self.moon_M, self.moon_e)
        self.assertAlmostEqual(E0, self.moon_E0, places = 4)

    moon_E = 262.9735

    def test_iterate_for_eccentric_anomaly(self):
        E = iterate_for_eccentric_anomaly(
            self.moon_M,
            self.moon_e,
            self.moon_E0,
        )
        self.assertAlmostEqual(E, self.moon_E, places = 4)

    moon_x = -10.68095
    moon_y = -59.72377

    def test_cartesian2d_in_plane_of_orbit(self):
        (x, y) = cartesian2d_in_plane_of_orbit(
            self.moon_E,
            self.moon_e,
            self.moon_a,
        )
        self.assertAlmostEqual(x, self.moon_x, places = 5)
        self.assertAlmostEqual(y, self.moon_y, places = 5)

    moon_r = 60.67134
    moon_theta = 259.8605

    def test_cartesian2d_to_polar_2(self):
        (r, theta) = cartesian2d_to_polar(self.moon_x, self.moon_y)
        self.assertAlmostEqual(r, self.moon_r, places = 5)
        self.assertAlmostEqual(theta, self.moon_theta, places = 4)

    moon_xeclip = +37.65311
    moon_yeclip = -47.57180
    moon_zeclip = -0.41687

    def test_position_from_plane_of_orbit_to_ecliptic(self):
        (xeclip, yeclip, zeclip) = position_from_plane_of_orbit_to_ecliptic(
            self.moon_r,
            self.moon_theta,
            self.moon_N,
            self.moon_i,
            self.moon_w,
        )
        self.assertAlmostEqual(xeclip, self.moon_xeclip, places = 4)
        self.assertAlmostEqual(yeclip, self.moon_yeclip, places = 4)
        self.assertAlmostEqual(zeclip, self.moon_zeclip, places = 4)

    moon_r_eclip = 60.6713
    moon_lat_eclip = 359.6063
    moon_lon_eclip = 308.3616

    def test_cartesian3d_to_spherical_4(self):
        (
            r,
            latitude,
            longitude,
        ) = cartesian3d_to_spherical(
            self.moon_xeclip,
            self.moon_yeclip,
            self.moon_zeclip,
        )
        self.assertAlmostEqual(r, self.moon_r_eclip, places = 4)
        self.assertAlmostEqual(latitude, self.moon_lat_eclip, places = 4)
        self.assertAlmostEqual(longitude, self.moon_lon_eclip, places = 4)

    # mean_longitude = 26.8388
    moon_mean_longitude = 314.5789
    # mean_anomaly = 104.0653
    moon_mean_anomaly = 266.0954
    moon_mean_elongation = 287.7401
    moon_argument_of_latitude = 1.8408

    def test_moon_perturbation_arguments(self):
        (Ls, Lm, Ms, Mm, D, F) = moon_perturbation_arguments(
            self.mean_longitude,
            self.moon_N,
            self.moon_w,
            self.moon_M,
            self.mean_anomaly,
        )
        self.assertAlmostEqual(Ls, self.mean_longitude, places = 4)
        self.assertAlmostEqual(Lm, self.moon_mean_longitude, places = 4)
        self.assertAlmostEqual(Ms, self.mean_anomaly, places = 4)
        self.assertAlmostEqual(Mm, self.moon_mean_anomaly, places = 4)
        self.assertAlmostEqual(D, self.moon_mean_elongation, places = 4)
        self.assertAlmostEqual(F, self.moon_argument_of_latitude, places = 4)

    moon_dlon = (
        -0.9847,
        -0.3819,
        -0.1804,
        +0.0405,
        -0.0244,
        +0.0452,
        +0.0428,
        +0.0126,
        -(-0.0333),
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
            self.mean_longitude,
            self.moon_mean_longitude,
            self.mean_anomaly,
            self.moon_mean_anomaly,
            self.moon_mean_elongation,
            self.moon_argument_of_latitude,
        )
        for p in zip(dlon, self.moon_dlon):
            self.assertAlmostEqual(p[0], p[1], places = 4)
        sum_dlon = fsum(dlon)
        self.assertAlmostEqual(sum_dlon, self.moon_sum_dlon, places = 4)

    def test_moon_5_latitude_perturbation_terms(self):
        dlat = moon_5_latitude_perturbation_terms(
            self.mean_longitude,
            self.moon_mean_longitude,
            self.mean_anomaly,
            self.moon_mean_anomaly,
            self.moon_mean_elongation,
            self.moon_argument_of_latitude,
        )
        for p in zip(dlat, self.moon_dlat):
            self.assertAlmostEqual(p[0], p[1], places = 4)
        sum_dlat = fsum(dlat)
        self.assertAlmostEqual(sum_dlat, self.moon_sum_dlat, places = 4)

    def test_moon_2_distance_perturbation_terms(self):
        ddist = moon_2_distance_perturbation_terms(
            self.mean_longitude,
            self.moon_mean_longitude,
            self.mean_anomaly,
            self.moon_mean_anomaly,
            self.moon_mean_elongation,
            self.moon_argument_of_latitude,
        )
        for p in zip(ddist, self.moon_ddist):
            self.assertAlmostEqual(p[0], p[1], places = 4)
        sum_ddist = fsum(ddist)
        self.assertAlmostEqual(sum_ddist, self.moon_sum_ddist, places = 4)

    moon_r_eclip_perturbed = 60.6779
    moon_lat_eclip_perturbed = 359.4144
    moon_lon_eclip_perturbed = 306.9484

    moon_r_eclip_astro_almanac = 60.793
    moon_lat_eclip_astro_almanac = 359.45
    moon_lon_eclip_astro_almanac = 306.94

    def test_date_and_sun_mean_to_moon_ecliptic(self):
        (r, lat, lon) = date_and_sun_mean_to_moon_ecliptic(
            self.Y,
            self.M,
            self.D,
            self.mean_longitude,
            self.mean_anomaly,
        )
        self.assertAlmostEqual(r, self.moon_r_eclip_perturbed, places = 4)
        self.assertAlmostEqual(lat, self.moon_lat_eclip_perturbed, places = 4)
        self.assertAlmostEqual(lon, self.moon_lon_eclip_perturbed, places = 4)

        self.assertAlmostEqual(
            r,
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

    moon_Decl = 340.8968
    moon_RA = 309.5011
    
    moon_Decl_astro_almanac = 340.9159
    moon_RA_astro_almanac = 309.4881

    def test_ecliptic_to_celestial(self):
        (distance1, Decl, RA) = ecliptic_to_celestial(
            1.0,
            self.moon_lat_eclip_perturbed,
            self.moon_lon_eclip_perturbed,
            self.obliquity_of_the_ecliptic,
        )
        self.assertAlmostEqual(Decl, self.moon_Decl, places = 4)
        self.assertAlmostEqual(RA, self.moon_RA, places = 4)
        
        self.assertAlmostEqual(Decl, self.moon_Decl_astro_almanac, places = 1)
        self.assertAlmostEqual(RA, self.moon_RA_astro_almanac, places = 1)

    moon_hour_angle_degrees = 272.3377

    def test_hour_angle(self):
        h = arcdegrees_to_hours(self.moon_RA)
        ha = hour_angle(self.sidtime0, h)
        d = hours_to_arcdegrees(ha)
        self.assertAlmostEqual(d, self.moon_hour_angle_degrees, places = 3)


if __name__ == '__main__':
    unittest.main()
