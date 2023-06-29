#===============================================================================
DIR_PATH = "20230530-0001/"
LEFT_BORDER  = 850
RIGHT_BORDER = 1300
#===============================================================================
import os
import numpy as np
import matplotlib.pyplot as plt
import statistics as stat

def get_files(dir_path,postfix=".txt"):
    res = []
    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(dir_path,path)):
            if path.endswith(postfix):
                res.append(os.path.join(dir_path,path))
    return res

def process_file(file_name):
    x=[]; y_in=[]; y_out=[]
    with open(file_name) as fl:
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
    return (energy, i_min, i_max )

def analyze_spectrum(spec,left_border=LEFT_BORDER, right_border=RIGHT_BORDER):
     chan       = np.argmax ( spec[2][1]        )
     val        = np.max    ( spec[2][1]        )
     left_base  = np.average( spec[2][1][0:left_border] )
     right_base = np.average( spec[2][1][right_border:] )
     energy     = calc_energy_baseline(spec,chan,max,left_base,right_base)
#     return ( np.argmax(spec[2][1]) , np.max(spec[2][1]), np.average( spec[2][1][0:850] ) , np.average( spec[2][1][1300:] ) )
     return ( chan , val, left_base, right_base , energy )

r = get_data(DIR_PATH)

energies = []

for s in r:
    m = analyze_spectrum(s)
    ss  = s[0]
    ss += "\t" + str( m[0]          )
    ss += "\t" + str( s[1][m[0]]    )
    ss += "\t" + str( s[2][1][m[0]] )
#    ss += "\t" + str( m[1] )
    ss += "\t" + str( m[2] )
    ss += "\t" + str( m[3] )
    ss += "\t" + str( m[4][0] )
    ss += "\t" + str( m[4][1] )
    ss += "\t" + str( m[4][2] )
    energies.append(m[4][0])
    print(ss)

print("Mean  Energy = " +str(stat.mean(energies))+" +/- "+str(stat.stdev(energies)/np.sqrt(1000.)))
print("StDev Energy = " +str(stat.stdev(energies)))
plt.hist( energies )
plt.show()
