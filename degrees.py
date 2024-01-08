'''
Provide versions of math functions for angles in degrees.
'''
# pylint: disable=invalid-name
import math
import numbers

class Degrees(float): # numbers.Real):
    '''
    An angle in degrees.
    '''
    def sin(self):
        'Return the sine of ``self`` (measured in degrees).'
        return math.sin(math.radians(self))

    def cos(self):
        'Return the cosine of ``self`` (measured in degrees).'
        return math.cos(math.radians(self))

    def rev(self):
        '''
        The remainder of ``self/360``.
        Note that the Python expression ``self % 360`` may not return the same
        result.  The intent is that ``x.rev()`` be exactly (mathematically; to infinite precision)
        equal to ``x - n*360.0`` for some integer ``n`` such that the result is non-negative and
        less than 360.
        '''
        return self.__class__(math.fmod(self, 360.0) + (self < 0) * 360.0)

    @classmethod
    def fromradians(cls, radians):
        'Get an angle in degrees equivalent to the given radians.'
        return cls(math.degrees(radians))

    @classmethod
    def atan2(cls, y, x):
        '''
        Return the arc tangent (measured in degrees) of ``y/x``.
        
        Unlike ``atan(y/x)``, the signs of both ``x`` and ``y`` are considered.
        '''
        return cls(math.degrees(math.atan2(y, x)))

    @classmethod
    def asin(cls, x):
        '''
        Return the arc sine (measured in degrees) of ``x``.
        
        The result is between 0 and 360.
        '''
        return cls(math.degrees(math.asin(x)))

def sin(x):
    'Return the sine of ``x`` (measured in degrees).'
    return math.sin(math.radians(x))

def cos(x):
    'Return the cosine of ``x`` (measured in degrees).'
    return math.cos(math.radians(x))

def atan2(y, x):
    '''
    Return the arc tangent (measured in degrees) of ``y/x``.
    
    Unlike ``atan(y/x)``, the signs of both ``x`` and ``y`` are considered.
    '''
    return math.degrees(math.atan2(y, x))

def asin(x):
    '''
    Return the arc sine (measured in degrees) of x.
    
    The result is between 0 and 360.
    '''
    return math.degrees(math.asin(x))

def rev(x):
    '''
    The remainder of ``x/360``. Note that the Python expression ``x % 360`` may not return the same
    result.  The intent is that ``rev(x)`` be exactly (mathematically; to infinite precision) equal
    to ``x - n*360.0`` for some integer ``n`` such that the result is non-negative and less than
    360.
    '''
    return math.fmod(x, 360.0) + (x < 0) * 360.0
