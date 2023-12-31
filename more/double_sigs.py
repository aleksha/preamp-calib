#===============================================================================
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

import tkinter as tk
from tkinter.ttk import Combobox
from tkinter.filedialog import askopenfilename, asksaveasfilename

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

from datetime import datetime

window = tk.Tk()

start_path_to_folder = "./20230530-0001/"
path_to_folder = tk.StringVar(value=start_path_to_folder)

methods = ["slopes", "baseline"]
method_var = tk.StringVar(value=methods[0]) 
pdf_var = tk.StringVar(value="output") 

left_border_var  = tk.IntVar(value=LEFT_BORDER ) 
right_border_var = tk.IntVar(value=RIGHT_BORDER) 
low_threshold_var  = tk.IntVar(value=int(Y_LOW*100.) ) 
high_threshold_var = tk.IntVar(value=int(Y_HIGH*100.)) 

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
    y_low = low_threshold_var.get()/100.
    y_high= high_threshold_var.get()/100.
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
    method = method_var.get()
    if method=="slopes":
        return calc_energy_slopes(spec,chan,max,left_base,right_base)
    if method=="baseline":
        return calc_energy_baseline(spec,chan,max,left_base,right_base)
    return calc_energy_slopes(spec,chan,max,left_base,right_base)

def analyze_spectrum(spec,left_border=LEFT_BORDER, right_border=RIGHT_BORDER):
    left_border = left_border_var.get()
    right_border = right_border_var.get()
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


def process_dir( path_to_dir , fig):
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
        gen.append( m[7])
        amp.append( m[1])

    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d %H:%M:%S")

    ss  = "================================================================================="
    ss += "\n  Date and time : " + dt_string
    ss += "\n------------------------------ CONDITIONS ---------------------------------------"
    ss += "\n  Folder  : " + path_to_dir
    ss += "\n  Method  : " + str(method_var.get())
    if str(method_var.get())=="slopes":
        ss += "    (linear part in [" + str(low_threshold_var.get()) + ","
        ss += str(high_threshold_var.get())+ "]% of spectrum height)"
    ss += "\n  Borders : [" + str(left_border_var.get()) + "," + str(right_border_var.get())+ "]"
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
    ss += "\n=================================================================================\n\n"

    for idx in range(len(r)):
        if idx==0:
            av  = r[idx][2][1]*float(1/len(r))
        else:
            av += r[idx][2][1]*float(1/len(r))

    cntL = 0.; cntH = 0.;
    for idx in range(len(r)):
        if durations[idx]<513:
            cntL+=1.
        else:
            cntH+=1.
    fst = True
    for idx in range(len(r)):
        if durations[idx]<513:
            if fst:
                avL = r[idx][2][1] / cntL
                fst = False
            else:
                avL += r[idx][2][1] / cntL

    fst = True
    for idx in range(len(r)):
        if durations[idx]>513:
            if fst:
                avH = r[idx][2][1] / cntH
                fst = False
            else:
                avH += r[idx][2][1] / cntH

    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()

#    plt.hist(energies)

#    plt.hist(starts)

#    plt = fig.add_subplot(111)
#    plt.clear()
    plt.plot(r[0][1],av,"b-")
    plt.plot(r[0][1],avL,"r-")
    plt.plot(r[0][1],avH,"g-")
#    plt.plot(r[SHOW_SPEC][1],r[SHOW_SPEC][2][1],"r-")
    plt.set_xlabel("Time [us]")
    plt.set_ylabel("Signal [a.u.]")
    rr = {"energies":energies,"starts":starts,"stops":stops,"durations":durations,"data":r,"av":av,"gen":gen,"amp":amp}
    return (ss, rr)

#===============================================================================
# GUI
#===============================================================================



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


window.title("PicoScope's Spectra Analysis")


txt_edit = tk.Text(window, width=100)


def browse_button():
    # Allow user to select a directory and store it in global var
    # called folder_path
    filename = tk.filedialog.askdirectory()
    path_to_folder.set(filename)

#scroll_bar = tk.Scrollbar(window)

fig_draw = Figure(figsize=(7,5), dpi=100)

canvas = FigureCanvasTkAgg(fig_draw, window)
canvas.draw()


def upd_result(txt=txt_edit,fig=fig_draw):
    folder=path_to_folder.get()
    print(folder)
    res = process_dir(folder,fig)
    last_ras = res
    canvas.draw()
    txt.insert(tk.END, res[0])
    global last_res
    last_res = res
    return res

xlab = {  "energies":"Energy [a.u.]",
          "starts":"Signal start [ch.]","stops":"Signal stop [ch.]",
          "durations":"Signal duration [ch.]","gen":"Input signal [a.u.]",
          "amp":"Amplitude [a.u.]"}

def draw_canv(what="energies",txt=txt_edit,fig=fig_draw):
    global last_res
    if last_res==None:
        last_res = upd_result(txt,fig)
    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
    plt.hist( last_res[1][what] )
    plt.set_xlabel(xlab[what])
    plt.set_ylabel("Entries")
    canvas.draw()

def draw_energies(txt=txt_edit,fig=fig_draw):
    draw_canv("energies",txt,fig)

def draw_starts(txt=txt_edit,fig=fig_draw):
    draw_canv("starts",txt,fig)

def draw_stops(txt=txt_edit,fig=fig_draw):
    draw_canv("stops",txt,fig)

def draw_durations(txt=txt_edit,fig=fig_draw):
    draw_canv("durations",txt,fig)

def draw_input(txt=txt_edit,fig=fig_draw):
    draw_canv("gen",txt,fig)

def draw_amp(txt=txt_edit,fig=fig_draw):
    draw_canv("amp",txt,fig)

def draw_spectra(txt=txt_edit,fig=fig_draw):
    global last_res
    if last_res==None:
        last_res = upd_result(txt,fig)
    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
    plt.plot(last_res[1]["data"][0][1],last_res[1]["av"],"b-")
    plt.plot(last_res[1]["data"][SHOW_SPEC][1],last_res[1]["data"][SHOW_SPEC][2][1],"r-")
    plt.set_xlabel("Time [us]")
    plt.set_ylabel("Signal [a.u.]")
    canvas.draw()

def draw_scatter(txt=txt_edit,fig=fig_draw):
    global last_res
    if last_res==None:
        last_res = upd_result(txt,fig)
    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
    plt.scatter(last_res[1]["gen"],last_res[1]["energies"])
    plt.set_xlabel("Input signal [a.u]")
    plt.set_ylabel('Energy [a.u.]')
    canvas.draw()


def create_pdf(pdf_name="output.pdf",txt=txt_edit,fig=fig_draw):
    global last_res
    if last_res==None:
        last_res = upd_result(txt,fig)

    pdf = PdfPages( str(pdf_var.get()) +".pdf")

    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
#    plt.plot([0,1],[0,1],"b")
    plt.text(0.0, 0.0, last_res[0], fontsize = 6)
    plt.axis('off')
    pdf.savefig( fig )

    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
    plt.scatter(last_res[1]["gen"],last_res[1]["energies"])
    plt.set_xlabel("Input signal [a.u]")
    plt.set_ylabel('Energy [a.u.]')
    pdf.savefig( fig )

    fig.clear()
    plt = fig.add_subplot(111)
    plt.clear()
    plt.plot(last_res[1]["data"][0][1],last_res[1]["av"],"b-")
    plt.plot(last_res[1]["data"][SHOW_SPEC][1],last_res[1]["data"][SHOW_SPEC][2][1],"r-")
    plt.set_xlabel("Time [us]")
    pdf.savefig( fig )


    lst_pdf = ["energies", "starts", "stops", "durations", "gen","amp"]
    for what in lst_pdf:
        fig.clear()
        plt = fig.add_subplot(111)
        plt.clear()
        plt.hist( last_res[1][what] )
        plt.set_xlabel(xlab[what])
        plt.set_ylabel("Entries")
        pdf.savefig( fig )

    pdf.close()


    
def dummy():
    print("Dummy function has been called")
    
frame_input = tk.Frame(window)
#frame_input = tk.Frame(window, relief=tk.RAISED, bd=2)

#lbl_open = tk.Label(window, text="Open", command=open_file)
lbl_open = tk.Label(frame_input, text="Folder with data: ")
lbl_open.grid(row=0,column=0)
txt_input = tk.Entry(frame_input, textvariable=path_to_folder, width=50)
txt_input.grid(row=0,column=1)
btn_browse = tk.Button(frame_input, text="Browse", command=browse_button)
btn_open = tk.Button(frame_input, text="Analyse", command=upd_result )
#btn_save = tk.Button(frame_input, text="Analyse", command=upd_result( txt_edit, path_to_folder.get() ) )
#btn_save = tk.Button(frame_input, text="Analyse", command=dummy )
btn_save = tk.Button(frame_input, text="Save log as...", command=save_file )
lbl_pdf = tk.Label(frame_input, text="PDF file: ")
ent_pdf = tk.Entry(frame_input, textvariable=pdf_var, width=15)
btn_pdf  = tk.Button(frame_input, text="PDF", command=create_pdf )
btn_browse.grid(row=0, column=2)
btn_open.grid(row=0,column=3)
btn_save.grid(row=0,column=4)
lbl_pdf.grid(row=0,column=5)
ent_pdf.grid(row=0,column=6)
btn_pdf.grid(row=0,column=7)

frame_input.grid(row=0,column=0)



frame_control = tk.Frame(window)
lbl_method = tk.Label(frame_control, text="Method: ")
cbx_method = Combobox(frame_control, textvariable=method_var, values=methods)

lbl_left = tk.Label(frame_control, text="Left border [ch]: ")
ent_left = tk.Entry(frame_control, textvariable=left_border_var, width=5)
lbl_right= tk.Label(frame_control, text="Right border [ch]: ")
ent_right= tk.Entry(frame_control, textvariable=right_border_var, width=5)

lbl_low  = tk.Label(frame_control, text="Low threshold [%]: ")
ent_low  = tk.Entry(frame_control, textvariable=low_threshold_var, width=4)
lbl_high = tk.Label(frame_control, text="High threshold [%]: ")
ent_high = tk.Entry(frame_control, textvariable=high_threshold_var, width=4)

lbl_method.grid(row=0,column=0)
cbx_method.grid(row=0,column=1)
lbl_left.grid(row=0,column=2)
ent_left.grid(row=0,column=3)
lbl_right.grid(row=0,column=4)
ent_right.grid(row=0,column=5)
lbl_low.grid(row=0,column=6)
ent_low.grid(row=0,column=7)
lbl_high.grid(row=0,column=8)
ent_high.grid(row=0,column=9)

frame_control.grid(row=1,column=0)

frame_draw = tk.Frame(window)
lbl_draw = tk.Label(frame_draw, text="Draw buttons: ")
btn_energies  = tk.Button(frame_draw, text="Energies", command=draw_energies )
btn_starts    = tk.Button(frame_draw, text="Starts", command=draw_starts )
btn_stops     = tk.Button(frame_draw, text="Stops", command=draw_stops  )
btn_durations = tk.Button(frame_draw, text="Durations", command=draw_durations  )
btn_spectra   = tk.Button(frame_draw, text="Spectra", command=draw_spectra  )
btn_input     = tk.Button(frame_draw, text="Input", command=draw_input )
btn_scatter   = tk.Button(frame_draw, text="Scatter", command=draw_scatter )
btn_amp       = tk.Button(frame_draw, text="Amplitude", command=draw_amp )
lbl_draw.grid(row=0,column=0)
btn_energies.grid(row=0,column=1)
btn_starts.grid(row=0,column=2)
btn_stops.grid(row=0,column=3)
btn_durations.grid(row=0,column=4)
btn_spectra.grid(row=0,column=5)
btn_input.grid(row=0,column=6)
btn_scatter.grid(row=0,column=7)
btn_amp.grid(row=0,column=8)
frame_draw.grid(row=2,column=0)


#canvas.get_tk_widget().grid(row=2,column=0,side=tk.BOTTOM, fill=tk.BOTH, expand=True)
canvas.get_tk_widget().grid(row=3,column=0)


txt_edit.grid(row=4,column=0)
v=tk.Scrollbar(window, orient='vertical', command=txt_edit.yview)
v.grid(row=4,column=1,sticky='nsew')

txt_edit['yscrollcommand'] = v.set
#v.config(command=txt_edit.yview)

window.mainloop()
