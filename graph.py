
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


# ================================================================================================
# | R23-P2 | 01 | 41.477 +/- 0.010 | 41.465 | =    | \33 | 41.4769 +/- 0.0053 | 67.409 +/- 0.378 |
# | R23-P2 | 01 | 41.451 +/- 0.011 | 41.460 | *    | \32 | 41.4507 +/- 0.0074 | 33.751 +/- 0.006 |
# | R23-P2 | 01 | 41.369 +/- 0.016 | 41.374 | **   | \31 | 41.3688 +/- 0.0141 | 16.857 +/- 0.003 |
# | R23-P2 | 01 | 41.255 +/- 0.028 | 41.241 | ***  | \29 | 41.2554 +/- 0.0272 |  8.410 +/- 0.003 |
# | R23-P2 | 01 | 41.107 +/- 0.053 | 41.089 | **** | \30 | 41.1070 +/- 0.0530 |  4.200 +/- 0.001 |
# ================================================================================================



x = [ VE(4.200, 0.001**2), VE(8.410,0.003**2), VE(16.857,0.003**2), VE(33.751,0.006**2), VE(67.409,0.378**2) ]
y = [ VE(41.1070,0.0530**2), VE(41.2554, 0.0272**2), VE(41.3688, 0.0141**2), VE(41.4507,0.0074**2), VE( 41.4769,0.0053**2) ]

a = [ VE(0.096626, 0.003186**2), VE(0.192465, 0.003176**2), VE(0.385012, 0.003612**2), VE(0.770968, 0.004575**2), VE(1.540838, 0.006938**2)]

gr = makeGraph(x,y)
gr.SetMarkerStyle(24)
gr.GetYaxis().SetRangeUser(39.9,42.1)
gr.Draw("AP")
#tf = ROOT.TF1("tf","([0]+[1]*x)/([2]+x)",0,100)
tf = ROOT.TF1("tf","[0]+[1]*x/([2]+x)",0,100)
r = gr.Fit( tf ,"S")
print(r)

#gr = makeGraph(x,a)
#r = gr.Fit("pol1","S")
#print(r)
