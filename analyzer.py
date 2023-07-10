#===============================================================================
#DIR_PATH     = "../20230530-0001/"
LEFT_BORDER  = 850
RIGHT_BORDER = 1300
Y_LOW        = 0.15
Y_HIGH       = 0.85
SHOW_SPEC    = 1
#===============================================================================
import os
import numpy as np
#import matplotlib.pyplot as plt
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

def calc_energy(spec,chan,max,left_base,right_base,method="slopes"):
    if method=="slopes":
        return calc_energy_slopes(spec,chan,max,left_base,right_base)
    if method=="baseline":
        return calc_energy_baseline(spec,chan,max,left_base,right_base)
    return calc_energy_slopes(spec,chan,max,left_base,right_base)

def analyze_spectrum(spec,left_border=LEFT_BORDER, right_border=RIGHT_BORDER):
     chan       = np.argmax ( spec[2][1]        )
     val        = np.max    ( spec[2][1]        )
     left_base  = np.average( spec[2][1][0:left_border] )
     right_base = np.average( spec[2][1][right_border:] )
#     energy     = calc_energy(spec,chan,max,left_base,right_base,method="baseline")
     energy     = calc_energy(spec,chan,max,left_base,right_base,method="slopes")
#     energy     = calc_energy_baseline(spec,chan,max,left_base,right_base)
#     return ( np.argmax(spec[2][1]) , np.max(spec[2][1]), np.average( spec[2][1][0:850] ) , np.average( spec[2][1][1300:] ) )
     return ( chan , val, left_base, right_base , energy )


def process_dir( path_to_dir , fig):
    r = get_data(path_to_dir)

    energies = []
    starts   = []
    stops    = []

    for s in r:
        m = analyze_spectrum(s)
        energies.append(m[4][0])
        starts.append(m[4][3])
        stops.append(m[4][4])

    ss  = path_to_dir
    ss += "\nMean  Energy = " +str(stat.mean(energies))+" +/- "+str(stat.stdev(energies)/np.sqrt(1000.))
    ss += "\nStDev Energy = " +str(stat.stdev(energies)) + "\n\n"

    for idx in range(len(r)):
        if idx==0:
            av  = r[idx][2][1]*float(1/len(r))
        else:
            av += r[idx][2][1]*float(1/len(r))

    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()

    plt.hist(energies)

#    plt.hist(starts)

#    plt = fig.add_subplot(111)
#    plt.clear()
#    plt.plot(r[0][1],av,"b-")
#    plt.plot(r[SHOW_SPEC][1],r[SHOW_SPEC][2][1],"r-")

    return ss

#===============================================================================
# GUI
#===============================================================================


import tkinter as tk
from tkinter.filedialog import askopenfilename, asksaveasfilename

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure



def open_file():
    """Open a file for editing."""
    filepath = askopenfilename(
        filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")]
    )
    if not filepath:
        return
    txt_edit.delete("1.0", tk.END)
    with open(filepath, mode="r", encoding="utf-8") as input_file:
        text = input_file.read()
        txt_edit.insert(tk.END, text)
    window.title(f"Simple Text Editor - {filepath}")

def save_file():
    """Save the current file as a new file."""
    filepath = asksaveasfilename(
        defaultextension=".txt",
        filetypes=[("Text Files", "*.txt"), ("All Files", "*.*")],
    )
    if not filepath:
        return
    with open(filepath, mode="w", encoding="utf-8") as output_file:
        text = txt_edit.get("1.0", tk.END)
        output_file.write(text)
    window.title(f"Simple Text Editor - {filepath}")


window = tk.Tk()
window.title("PicoScope's Spectra Analysis")

#window.rowconfigure(0, minsize=50, weight=1)
#window.columnconfigure(1, minsize=100, weight=1)

start_path_to_folder = "/home/adzyuba/dev/AMBER/20230530-0001/"
path_to_folder = tk.StringVar(value=start_path_to_folder)
txt_edit = tk.Text(window, width=100)

#scroll_bar = tk.Scrollbar(window)

fig_draw = Figure(figsize=(7,5), dpi=100)

canvas = FigureCanvasTkAgg(fig_draw, window)
canvas.draw()


def upd_result(txt=txt_edit,fig=fig_draw):
    folder=path_to_folder.get()
    print(folder)
    res = process_dir(folder,fig)
    canvas.draw()
    txt.insert(tk.END, res)


frame_input = tk.Frame(window)
#frame_input = tk.Frame(window, relief=tk.RAISED, bd=2)

#lbl_open = tk.Label(window, text="Open", command=open_file)
lbl_open = tk.Label(frame_input, text="Folder with data: ")
lbl_open.grid(row=0,column=0)
txt_input = tk.Entry(frame_input, textvariable=path_to_folder, width=80)
txt_input.grid(row=0,column=1)
btn_open = tk.Button(frame_input, text="Analyse", command=upd_result )
#btn_save = tk.Button(frame_input, text="Analyse", command=upd_result( txt_edit, path_to_folder.get() ) )
#btn_save = tk.Button(frame_input, text="Analyse", command=dummy )
btn_open.grid(row=0,column=2)

#lbl_open.pack()
#txt_input.pack()
#btn_save.pack()

btn_save = tk.Button(window, text="Save log as...", command=save_file )

frame_input.grid(row=0,column=0)

#canvas.get_tk_widget().grid(row=2,column=0,side=tk.BOTTOM, fill=tk.BOTH, expand=True)
canvas.get_tk_widget().grid(row=1,column=0)


txt_edit.grid(row=2,column=0)
btn_save.grid(row=3,column=0)


window.mainloop()
