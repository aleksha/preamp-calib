USE_TIME = False

import numpy as np
import matplotlib.pyplot as plt

x=[]; y_in=[]; y_out=[]

with open("test_file.txt") as fl:
    cnt=0
    for line in fl:
        if cnt>2 and len(line)>2:
            w=line[:-1].split("\t")
#            print( w )
            x    .append( float( w[0] ) )
            y_in .append( float( w[1] ) )
            y_out.append( float( w[2] ) )
        cnt+=1

c = []
for i in range(len(x)):
    c.append(i)

if USE_TIME:
    plt.plot( np.array(x), np.array(y_in),"r-")
    plt.plot( np.array(x), np.array(y_out)*100.,"b-")
else:
    plt.plot( np.array(c), np.array(y_in),"r-")
    plt.plot( np.array(c), np.array(y_out)*100.,"b-")
plt.show()
