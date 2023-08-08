#===============================================================================
path_to_folder = "../../data//20230808-0033"
#DIR_PATH     = "../20230605-0005/"
LEFT_BORDER  = 850
RIGHT_BORDER = 1300
Y_LOW        = 0.15
Y_HIGH       = 0.85
SHOW_SPEC    = 1
#===============================================================================
last_res = None # global var with recent results
#===============================================================================

import os
import numpy as np
#import matplotlib.pyplot as plt
import statistics as stat
from statistics import mean, median, stdev
from scipy.stats import pearsonr as rho

import tkinter as tk
from tkinter.ttk import Combobox
from tkinter.filedialog import askopenfilename, asksaveasfilename


from datetime import datetime





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




def get_files(dir_path,postfix=".txt"):
    res = []
    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path,path)):
            if path.endswith(postfix):
                res.append(os.path.join(dir_path,path))
    return res

def process_file(file_name):
    x=[]; y_in=[]; y_out=[]
    with open(file_name, encoding="utf-8") as fl:
        cnt=0
        for line in fl:
            if cnt>2 and len(line)>2:
                w=line[:-1].split("\t")
                x    .append( float( w[0] ) )
                y_in .append( float( w[1] ) )
                y_out.append( float( w[2] ) )
            cnt+=1
    return ( file_name, np.array(x), ( np.array(y_in), np.array(y_out) ) )


def get_data(dir_path,report=False):
    r = get_files(dir_path)
    super_r = []
    for f in r:
        if report:
            print(f)
        super_r.append( process_file(f) )
    if report:
        print( str(len(super_r)) + " files proceed")
    return super_r

def calc_energy_baseline(spec,chan,max,left_base,right_base):
    energy = 0.
    i = chan
    while spec[2][1][i-1]>left_base:
        energy += ( (spec[1][i]-spec[1][i-1])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        i-=1
    i_min = i
    i = chan
    while spec[2][1][i-1]>right_base:
        energy += ( (spec[1][i+1]-spec[1][i])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        i+=1
    i_max = i
    return (energy, i_min, i_max, spec[1][i_min], spec[1][i_max])

def find_x(y,x1,y1,x2,y2):
    if x1==x2:
        print("ERROR: find_x() : x1==x2 ==> same x")
        return 0.
    if y1==y2:
        print("ERROR: find_x() : y1==y2 ==> any x possible")
        return 0.
    return x1 + (y-y1)*(x2-x1)/(y2-y1)

def calc_energy_slopes(spec,chan,max,left_base,right_base,y_low=Y_LOW,y_high=Y_HIGH):
    y_low = spec[2][1][chan]*y_low
    y_high= spec[2][1][chan]*y_high
    energy = 0.

    # go left
    i = chan
    i_low = chan
    i_high = chan
    while spec[2][1][i]>y_high:
        energy += ( (spec[1][i]-spec[1][i-1])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        if spec[2][1][i-1]<y_high:
            i_high = i-1
        i-=1
    while spec[2][1][i]>y_low:
        energy += ( (spec[1][i]-spec[1][i-1])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        if spec[2][1][i-1]<y_low:
            i_low = i-1
        i-=1
    i_min = int( find_x(left_base,i_low,spec[2][1][i_low],i_high,spec[2][1][i_high]) )
    x_min = find_x(left_base,i_low,spec[2][1][i_low],i_high,spec[2][1][i_high])
    while i>i_min:
        energy += ( (spec[1][i]-spec[1][i-1])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        i-=1

    # go right
    i = chan
    while spec[2][1][i]>y_high:
        energy += ( (spec[1][i+1]-spec[1][i])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        if spec[2][1][i+1]<y_high:
            i_high = i+1
        i+=1
    while spec[2][1][i]>y_low:
        energy += ( (spec[1][i+1]-spec[1][i])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        if spec[2][1][i+1]<y_low:
            i_low = i+1
        i+=1
    i_max = int( find_x(left_base,i_low,spec[2][1][i_low],i_high,spec[2][1][i_high]) )
    x_max = find_x(left_base,i_low,spec[2][1][i_low],i_high,spec[2][1][i_high])
    while i<i_max:
        energy += ( (spec[1][i+1]-spec[1][i])*0.5*( spec[2][1][i]+spec[2][1][i+1] ) )
        i+=1

    return (energy, i_min, i_max, x_min, x_max )

def calc_energy(spec,chan,max,left_base,right_base):
    method = "slopes"
    if method=="slopes":
        return calc_energy_slopes(spec,chan,max,left_base,right_base)
    if method=="baseline":
        return calc_energy_baseline(spec,chan,max,left_base,right_base)
    return calc_energy_slopes(spec,chan,max,left_base,right_base)

def analyze_spectrum(spec,left_border=LEFT_BORDER, right_border=RIGHT_BORDER):
    chan       = np.argmax ( spec[2][1]        )
    val        = np.max    ( spec[2][1]        )
    left_in    = np.average( spec[2][0][0:left_border] )
    right_in   = np.average( spec[2][0][right_border:] )
    left_base  = np.average( spec[2][1][0:left_border] )
    right_base = np.average( spec[2][1][right_border:] )
#     energy     = calc_energy(spec,chan,max,left_base,right_base,method="baseline")
    energy     = calc_energy(spec,chan,max,left_base,right_base)
#     energy     = calc_energy_baseline(spec,chan,max,left_base,right_base)
#     return ( np.argmax(spec[2][1]) , np.max(spec[2][1]), np.average( spec[2][1][0:850] ) , np.average( spec[2][1][1300:] ) )
    return ( chan , val, left_base, right_base , energy, left_in, right_in, left_in - right_in )

def ratio_unc(A,sA,B,sB,rho):
    if B!=0 and A!=0:
        return A/B, abs(A/B)*np.sqrt(pow(sA/A,2)+pow(sB/B,2))
    return 0.

def process_dir( path_to_dir ):
    r = get_data(path_to_dir)

    energies  = []
    starts    = []
    stops     = []
    durations = []
    gen       = []
    amp       = []

    for s in r:
        m = analyze_spectrum(s)
        energies.append(m[4][0])
        starts.append(m[4][3])
        stops.append(m[4][4])
        durations.append( m[4][4]-m[4][3])
        gen.append( m[7]*0.001)
        amp.append( m[1])

    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")

    ss  = "================================================================================="
    ss += "\n  Date and time : " + dt_string
    ss += "\n------------------------------ CONDITIONS ---------------------------------------"
    ss += "\n  Folder  : " + path_to_dir
    ss += "\n------------------------------ RESULTS ------------------------------------------"
    ss += "\n  Mean  Input             : " + str(stat.mean(gen))
    ss += "\n  StDev of Mean Input     : " + str(stat.stdev(gen)/np.sqrt(float(len(gen))))
    ss += "\n  Median Input            : " + str(stat.median(gen))
    ss += "\n  StDev of Input          : " + str(stat.stdev(gen))
    ss += "\n  Mean  Energy            : " + str(stat.mean(energies))
    ss += "\n  StDev of Mean Energy    : " + str(stat.stdev(energies)/np.sqrt(float(len(energies))))
    ss += "\n  Median Energy           : " + str(stat.median(energies))
    ss += "\n  StDev of Energy         : " + str(stat.stdev(energies))
    ss += "\n  Mean of Amplitude       : " + str(stat.mean(amp))
    ss += "\n  StDev of Mean Amplitude : " + str(stat.stdev(amp)/np.sqrt(float(len(amp))))
    ss += "\n  Median Amplitude        : " + str(stat.median(amp))
    ss += "\n  StDev of Amplitude      : " + str(stat.stdev(amp))
    ss += "\n  Mean of Start Time      : " + str(stat.mean(starts))
    ss += "\n  StDev of Start Time     : " + str(stat.stdev(starts))
    ss += "\n  Mean of Stop Time       : " + str(stat.mean(stops))
    ss += "\n  StDev of Stop Time      : " + str(stat.stdev(stops))
    ss += "\n  Mean of Duration        : " + str(stat.mean(durations))
    ss += "\n  StDev of Duration       : " + str(stat.stdev(durations))
    ss += "\n  Correlation             : " + str(rho(energies,gen)[0])
    ss += "\n=================================================================================\n\n"


    rat = ratio_unc( stat.mean(energies),
                     stat.stdev(energies)/np.sqrt(float(len(energies))),
                     stat.mean(gen),
                     stat.stdev(gen)/np.sqrt(float(len(gen))),
                     rho(energies,gen)[0])

    log_ss  = " | "   + "{:.3f}".format( rat[0] )
    log_ss += " +/- " + "{:.3f}".format( rat[1] )
    log_ss += " | "   + "{:.3f}".format( stat.median(energies)/stat.median(gen) )
    log_ss += " | "
    print(log_ss)

    for idx in range(len(r)):
        if idx==0:
            av  = r[idx][2][1]*float(1/len(r))
        else:
            av += r[idx][2][1]*float(1/len(r))

    rr = {"energies":energies,"starts":starts,"stops":stops,"durations":durations,"data":r,"av":av,"gen":gen,"amp":amp}
    return (ss, rr)

#===============================================================================
# GUI
#===============================================================================

sss, rrr = process_dir( path_to_folder )

x = rrr["gen"]
y = rrr["energies"]



xx = ROOT.RooRealVar("xx","xx",0.0,0.1)
yy = ROOT.RooRealVar("yy","yy",0.0,4.0)



varset = ROOT.RooArgSet(xx,yy)
ds = ROOT.RooDataSet("ds","ds",varset)



for i in range(len(x)):
    xx.setVal(x[i])
    yy.setVal(y[i])
    ds.add(varset)



model = Models.Gauss2D_pdf("model",xvar=xx,yvar=yy)
ro, wo = model.fitTo( ds, draw= True,silent=True)
ro, wo = model.fitTo( ds, draw= True,silent=True)
ro, wo = model.fitTo( ds, draw= True,silent=True)
ro, wo = model.fitTo( ds, draw= True,silent=True)
ro, wo = model.fitTo( ds, draw= True,silent=True)
ro, wo = model.fitTo( ds, draw= True,silent=True)

print(ro)
print( ro.divide( 'mu_y_model' , 'mu_x_model' ) )
