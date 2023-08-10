
import ROOT
import ostap.fixes.fixes
from ostap.core.core import cpp, Ostap
from ostap.core.core import pwd, cwd, ROOTCWD
from ostap.core.core import rootID, funcID, funID, fID, histoID, hID, dsID
from ostap.core.core import VE
from ostap.histos.histos import h1_axis, h2_axes, h3_axes
from ostap.histos.graphs import makeGraph, hToGraph, hToGraph2, hToGraph3, lw_graph
import ostap.trees.trees
import ostap.trees.cuts
import ostap.histos.param
import ostap.histos.compare
import ostap.io.root_file
import ostap.math.models
import ostap.fitting.roofit
import ostap.fitting.models as Models


# ===============================================================================
# | 2023-08-08, CERN, Aleksei DZIUBA (PNPI), 20 C                               |
# ===============================================================================
# | Amplf. | Ch |    Integral      |        Input, mV       |   Amplitude, V    |
# -------------------------------------------------------------------------------
# | R23-P2 | 01 | 34.779 +/- 0.045 | 0.0049657 +/- 0.000002 | 0.0966 +/- 0.0032 |
# | R23-P2 | 01 | 34.875 +/- 0.024 | 0.0099487 +/- 0.000004 | 0.1925 +/- 0.0031 |
# | R23-P2 | 01 | 34.964 +/- 0.013 | 0.0199447 +/- 0.000004 | 0.3850 +/- 0.0036 |
# | R23-P2 | 01 | 35.025 +/- 0.006 | 0.0399431 +/- 0.000007 | 0.7710 +/- 0.0046 |
# | R23-P2 | 01 | 35.073 +/- 0.004 | 0.0797188 +/- 0.000014 | 1.5408 +/- 0.0069 |
# ===============================================================================


x = [ VE(4.500, 0.002**2), VE(9.949,0.004**2), VE(19.945,0.004**2), VE(39.943,0.007**2), VE(79.719,0.014**2) ]
y = [ VE(34.779, 0.045**2), VE(34.875, 0.024**2), VE(34.964, 0.013**2), VE(35.025, 0.006**2), VE(35.073, 0.004**2) ]

a = [ VE(0.096626, 0.003186**2), VE(0.192465, 0.003176**2), VE(0.385012, 0.003612**2), VE(0.770968, 0.004575**2), VE(1.540838, 0.006938**2)]

#gr = makeGraph(x,y)
#gr.SetMarkerStyle(24)
#gr.GetYaxis().SetRangeUser(33.9,36.1)
#gr.Draw("AP")
#tf = ROOT.TF1("tf","[0]+[1]*x/([2]+x)",0,100)
#r = gr.Fit( tf ,"S")
#print(r)

gr = makeGraph(x,a)
gr.Draw("AP")
r = gr.Fit("pol1","S")
print(r)
