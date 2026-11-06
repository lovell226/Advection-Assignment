import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import squareWave
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2
from advection import compareSchemeSnapshots
from growingErrors import growingErrors

#want to test for different sized waves - how do I do this while keeping a single wave?
#currently have a funciton that calls just x but this would need aditional params unless make 
#multiple funcitons

schemes = [FTBS, CTCS, SemiLagrangian,SemiLagrangian2]

    
compareSchemeSnapshots(squareWave, schemes, "SquareWave.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.3)


growingErrors(squareWave, schemes, "SchemesCosErrors_Square.pdf",
              nx = 40, nt = 40, u = 1, endTime = 0.4,xLim = [0,1], a = 0, b = 0.3)


compareSchemeSnapshots(squareWave, schemes, "SquareWave_nt=8.pdf",
                    nx = 40, nt = 12, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.3)