import matplotlib.pyplot as plt
import numpy as np
import os
from advectionHelpers import cosBell
from schemesBasic import FTBS, CTCS, FTCS, SemiLagrangian, SemiLagrangian2
from advection import compareSchemeSnapshots

#want to test courrant numbers 0.2, 0.4, 0.8, 1.0, 1.5, 2.0, 5.0
#if we keep u fixed as 1? and change nt instead?



n=[80, 40, 20, 16, 10, 8, 2]
schemes = [FTBS, CTCS, FTCS, SemiLagrangian,SemiLagrangian2]

#Is there a way to fix the scales for this? so more obvious compared to analytical sol
    
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=80.pdf",
                    nx = 40, nt = n[0], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)

compareSchemeSnapshots(cosBell, schemes, "Cos_nt=40.pdf",
                    nx = 40, nt = n[1], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=20.pdf",
                    nx = 40, nt = n[2], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=16.pdf",
                    nx = 40, nt = n[3], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=10.pdf",
                    nx = 40, nt = n[4], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=8.pdf",
                    nx = 40, nt = n[5], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)
compareSchemeSnapshots(cosBell, schemes, "Cos_nt=2.pdf",
                    nx = 40, nt = n[6], u = 1, endTime = 0.4, xLim = [0,1], a = 0, b = 0.5)