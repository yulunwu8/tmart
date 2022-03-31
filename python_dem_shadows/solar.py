import numpy as np
from numpy import sin, cos, tan, arcsin


def to_juliandate(d):
    """Convert a datetime object to a julian date.

    A Julian date is the decimal number of days since January 1, 4713 BCE."""
    seconds_per_day = 86400
    return d.timestamp() / seconds_per_day + 2440587.5


# TODO TdR 07/10/16: test
def sun_vector(julian_date, latitude, longitude, timezone):
    """Calculate a unit vector in the direction of the sun.

    :param date_time: utc timestamp, time of observation
    :param latitude: latitude of observation
    :param longitude: longitude of observation
    :param timezone: timezone hour offset relative to UTC
    :return: 3-dimensional unit vector.
    """
    # TODO TdR 07/10/16: verify computations
    omega_r = _hour_angle(julian_date, longitude, timezone)
    delta_r = np.deg2rad(sun_declination(julian_date))
    lambda_r = np.deg2rad(latitude)

    # TODO TdR 27/09/16: factor out common computation
    svx = - sin(omega_r) * cos(delta_r)
    svy = sin(lambda_r) * cos(omega_r) * cos(delta_r) \
        - cos(lambda_r) * sin(delta_r)
    svz = cos(lambda_r) * cos(omega_r) * cos(delta_r) \
        + sin(lambda_r) * sin(delta_r)
    return np.array([svx, svy, svz])


def sun_declination(julian_date):
    """Compute the declination of the sun on a given day."""
    # TODO TdR 27/09/16: verify calculations.
    jdc = (julian_date - 2451545.0) / 36525.0
    sec = 21.448 - jdc * (46.8150 + jdc * (0.00059 - jdc * .001813))
    e0 = 23.0 + (26.0 + (sec / 60.0)) / 60.0
    oblcorr = e0 + 0.00256 * cos(np.deg2rad(125.04 - 1934.136 * jdc))
    l0 = 280.46646 + jdc * (36000.76983 + jdc * 0.0003032)
    l0 = (l0 - 360 * (l0 // 360)) % 360
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmas = np.deg2rad(gmas)
    seqcent = sin(gmas) * (1.914602 - jdc * (0.004817 + 0.000014 * jdc)) + \
        sin(2 * gmas) * (0.019993 - 0.000101 * jdc) + sin(3 * gmas) * 0.000289

    suntl = l0 + seqcent
    sal = suntl - 0.00569 - 0.00478 * sin(np.deg2rad(125.04 - 1934.136 * jdc))
    delta = arcsin(sin(np.deg2rad(oblcorr)) * sin(np.deg2rad(sal)))
    return np.rad2deg(delta)


def _equation_of_time(julian_date):
    """Calculate the equation of time.

    See https://en.wikipedia.org/wiki/Equation_of_time.
    """
    # TODO TdR 07/10/16: verify computations.
    jdc = (julian_date - 2451545.0) / 36525.0
    sec = 21.448 - jdc * (46.8150 + jdc * (0.00059 - jdc * 0.001813))
    e0 = 23.0 + (26.0 + (sec / 60.0)) / 60.0
    oblcorr = e0 + 0.00256 * cos(np.deg2rad(125.04 - 1934.136 * jdc))
    l0 = 280.46646 + jdc * (36000.76983 + jdc * 0.0003032)
    l0 = (l0 - 360 * (l0 // 360)) % 360
    gmas = 357.52911 + jdc * (35999.05029 - 0.0001537 * jdc)
    gmas = np.deg2rad(gmas)

    ecc = 0.016708634 - jdc * (0.000042037 + 0.0000001267 * jdc)
    y = (tan(np.deg2rad(oblcorr) / 2)) ** 2
    rl0 = np.deg2rad(l0)
    EqTime = y * sin(2 * rl0) \
        - 2.0 * ecc * sin(gmas) \
        + 4.0 * ecc * y * sin(gmas) * cos(2 * rl0)\
        - 0.5 * y * y * sin(4 * rl0) \
        - 1.25 * ecc * ecc * sin(2 * gmas)
    return np.rad2deg(EqTime) * 4


# TODO TdR 07/10/16: test
def _hour_angle(julian_date, longitude, timezone):
    """Internal function for solar position calculation."""
    # TODO TdR 07/10/16: verify computations
    hour = ((julian_date - np.floor(julian_date)) * 24 + 12) % 24
    time_offset = _equation_of_time(julian_date)
    standard_meridian = timezone * 15
    delta_longitude_time = (longitude - standard_meridian) * 24.0 / 360.0
    omega_r = np.pi * (
        ((hour + delta_longitude_time + time_offset / 60) / 12.0) - 1.0)
    return omega_r



