Advection assignment. Hilary Weller <h.weller@reading.ac.uk> & Becks Lovell <r.lovell@pgr.reading.ac.uk>
Code to solve the advection equation.
Report_Run.py         : Run this to save plots of snapshots of the results
                        for schemes and for all metric plots.
function_setup.py     : Functions for evolving the schemes, initialising experiment
                        parameters and plotting function.                        
advection.py          : Function to save plots of snapshots of the results
                        for schemes.
advectionHelpers.py   : Functions for initialising and calculating error metrics
schemesBasic.py       : FTCS, CTCS, FTBS, Semi-Lagrangian (with linear and cubic interpolation)
                        advection schemes.
