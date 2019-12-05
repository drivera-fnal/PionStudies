from ROOT import *

def is_cosmic( nUpper, nLower, upperLimit = 5, lowerLimit = 5 ):
  if nUpper >= upperLimit and nLower > lowerLimit:
    return True
  else: return False
  
