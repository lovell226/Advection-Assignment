import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell, squareWave
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2
from advection import compareSchemeSnapshots

#want to test courrant numbers 0.2, 0.4, 0.8, 1.0, 1.5, 2.0, 5.0
#if we keep u fixed as 1? and change nt instead?



schemes = [FTBS, CTCS, SemiLagrangian,SemiLagrangian2]

    
compareSchemeSnapshots(cosBell, schemes, "basicSchemesCosMass.pdf",
                    nx = 40, nt = 1200, u = 1, endTime = 12.0, xLim = [0,1], a=0, b=0.5)

compareSchemeSnapshots(cosBell, schemes, "basicSchemesCosMass2.pdf",
                    nx = 40, nt = 12000, u = 1, endTime = 100.0, xLim = [0,1], a=0, b=0.5)

compareSchemeSnapshots(squareWave, schemes, "basicSchemesSquareMass.pdf",
                    nx = 40, nt = 1200, u = 1, endTime = 12.0, xLim = [0,1], a=0, b=0.5)
