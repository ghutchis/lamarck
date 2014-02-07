#!/usr/bin/env python

import math

def gaussianFraction (x, center, width):
    # FWHM = 2(sqrt(2 ln 2)) * c
    c = width / (2.0* sqrt(2.0 * math.log(2.0)))
    return math.exp(-1.0*(x - center)^2 / (2.0 * c^2))

