# preamp-calib
A software to perform calibration of pre-amplifiers for the time-projection chamber (TPC) readout.

## Application

`analyzer.py` - application with GUI (Tkinter, Matplotlib, NumPy) 

To test unzip `20230605-0005.zip`

## Interactive scripts

`read_file.py` - draw pair PicoScope's spectra (channels A ans B) from one measurement

`process_dir.py` - analyse sample of measurements

## Energy calculation

Two mwthods are used. In both simple energy summation is used, i.e.
histogram approximation: `energy += (x[i]+x[i+1])*(0.5*(y[i]+y[i+1]))`

### Baseline based method

Go left (right) from the peak position and sum energies 
untill it is higher than left (right) baseline.

### Slopes based method

Find left (right) point as an intersept of left (right) slopes with left (right) baseline.
Slopes treated as lenear functions between 15% and 85% of maximum.

### Comparison of method

Slopes based mathod is found to be more effective extimate as it results with
energy distribution with smaller standard deviation wrt. baseline based method.

 * StDev Energy (slopes    based) = 0.0187
 * StDev Energy (baselines based) = 0.0275

