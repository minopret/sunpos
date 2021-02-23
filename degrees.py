import math


def sin(x):
    return math.sin(x * math.pi / 180)


def cos(x):
    return math.cos(x * math.pi / 180)


def atan2(y, x):
    return ((180 / math.pi) * math.atan2(y, x))


def asin(x):
    return ((180 / math.pi) * math.asin(x))


def rev(x):
    return (x - math.floor(x / 360.0) * 360.0)
