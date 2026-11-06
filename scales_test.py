import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2
from advection import compareSchemeSnapshots
from growingErrors import growingErrors

#want to test for different sized waves - how do I do this while keeping a single wave?
#currently have a funciton that calls just x but this would need aditional params unless make 
#multiple funcitons


schemes = [FTBS, CTCS, SemiLagrangian,SemiLagrangian2]

    
compareSchemeSnapshots(cosBell, schemes, "CosScale_b=0.5.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)

compareSchemeSnapshots(cosBell, schemes, "CosScale2_b=0.3.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.3)

compareSchemeSnapshots(cosBell, schemes, "CosScale_b=1.0.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 1.0)

compareSchemeSnapshots(cosBell, schemes, "CosScale4_b=0.15.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.15)

growingErrors(cosBell, schemes, "SchemesCosErrors_scales_b=0.3.pdf",
              nx = 40, nt = 40, u = 1, endTime = 0.4,xLim = [0,1], a = 0, b = 0.3)

compareSchemeSnapshots(cosBell, schemes, "CosScale4_b=0.05.pdf",
                    nx = 40, nt = 40, u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.05)