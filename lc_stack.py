#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-01-17
Version: 1.0
Description: Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022.

Contact: nathedd@unc.edu
"""

import math
import matplotlib.pyplot as plt
import sys

# storing variable data as dictionaries: key = index, value = column of file_name
jd = {}  # julian day; not used in these methods but may be passed to other files
forcediffimflux = {}
forcediffimfluxunc = {}
forceddiffimchisq = {}
zpdiff = {}  # zeropoint
filter = {}

forcediffimflux_new = {}  # not used in these methods but may be passed to other files
flux_by_filter = {}  # {filter, list of flux}
forcediffimfluxunc_new = {}  # not used in these methods but may be passed to other files
unc_by_filter = {}  # {filter, list of unc}

# used for plotting
combined_flux = {}
combined_unc = {}
combined_start = {}
combined_end = {}

windows = {}  # for helper function get_indices


def fill_vars(data):
    "Takes a string of ascii data and converts it to dictionary variables. Note: If file used is of a different format than ZTF file, column indices will need to be changed to match your file."
    for line in data:
        line = line.strip()
        columns = line.split()
        if (columns[24] != 'null' and columns[25] != 'null'): 
            forcediffimflux[int(columns[0])] = (float(columns[24]))  # creates dictionary of {index, flux}
            forcediffimfluxunc[int(columns[0])] = (float(columns[25]))  # creates dictionary of {index, flux uncertainty}
            zpdiff[int(columns[0])] = (float(columns[20]))  # creates dictionary of {index, zero points}
            filter[int(columns[0])] = (columns[4])  # creates dictionary of {index, filter}
            jd[int(columns[0])] = (float(columns[22]))  # creates dictionary of {index, julian day}
        if (columns[27] != 'null'):
            forceddiffimchisq[int(columns[0])] = (float(columns[27]))  # creates dictionary of {index, chi^2}


def correct_baseline():
    """Following an estimate of the baseline level, subtract estimate from differential flux measurements."""
    for index in forcediffimflux:
        forcediffimflux[index] = forcediffimflux[index] - baseline


def validate_uncertainties():
    """Check the distribution of PSF-fit reduced chi squared values for a given filter."""
    g_list = []
    r_list = []
    i_list = []
    # placeholders
    g_avg = 1
    r_avg = 1
    i_avg = 1
    for index in forceddiffimchisq: # sort chi squared by filter
        new = forceddiffimchisq[index]
        if filter[index] == 'ZTF_g':
            g_list.append(new)
        elif filter[index] == 'ZTF_r':
            r_list.append(new)
        elif filter[index] == 'ZTF_i':
            i_list.append(new)
    if len(g_list) != 0:
        g_avg = sum(g_list) / len(g_list)
    if len(r_list) != 0:
        r_avg = sum(r_list) / len(r_list)
    if len(i_list) != 0:
        i_avg = sum(i_list) / len(i_list)
# ensure pixel uncertainties are consistent
    if round(g_avg) != 1:
        for index in forcediffimfluxunc:
            new = forcediffimfluxunc[index] * math.sqrt(forceddiffimchisq[index])
            if filter[index] == 'ZTF_g':
                forcediffimfluxunc[index] = new
    if round(r_avg) != 1:
        if filter[index] == 'ZTF_r':
                forcediffimfluxunc[index] = new
    if round(i_avg) != 1:
        if filter[index] == 'ZTF_i':
                forcediffimfluxunc[index] = new


def rescale(start, end):  # start, end inclusive
    """Make new columns for forcediffimfluxi and forcediffimfluxunci with rescaled input fluxes and uncertainties and sort new fluxes and uncertainties by filter."""
    zpavg =  min(zpdiff.values())  # flux / uncertainty ratio is consistent for any value; picking min for simplicity
    g_list = []  # list of rescaled fluxes in the g band
    r_list = []  # list of rescaled fluxes in the r band
    i_list = []  # list of rescaled fluxes in the i band
    g_unc_list = []  # list of rescaled uncertainties in the g band
    r_unc_list = []  # list of rescaled uncertainties in the r band
    i_unc_list = []  # list of rescaled uncertainties in the i band
    for index in range(start, end+1): 
        if index in forcediffimflux:
            new = forcediffimflux[index]*10**(0.4*(zpavg-zpdiff[index]))  # place fluxes on the same photometric zeropoint
            forcediffimflux_new[index] = new
            # sorts new flux by filter and matches dictionary of jd to flux by filter by index
            if filter[index] == 'ZTF_g':
                g_list.append(new)
            elif filter[index] == 'ZTF_r':
                r_list.append(new)
            elif filter[index] == 'ZTF_i':
                i_list.append(new)
            new = forcediffimfluxunc[index]*10**(0.4*(zpavg-zpdiff[index]))  # place uncertainties on the same photometric zeropoint
            forcediffimfluxunc_new[index] = new
            # sorts new uncertainty by filter
            if filter[index] == 'ZTF_g':
                g_unc_list.append(new)
            elif filter[index] == 'ZTF_r':
                r_unc_list.append(new)
            elif filter[index] == 'ZTF_i':
                i_unc_list.append(new)
    flux_by_filter['ZTF_g'] = g_list
    flux_by_filter['ZTF_r'] = r_list
    flux_by_filter['ZTF_i'] = i_list
    unc_by_filter['ZTF_g'] = g_unc_list
    unc_by_filter['ZTF_r'] = r_unc_list
    unc_by_filter['ZTF_i'] = i_unc_list


def collapse_flux_by_filter(start, end):
    """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average. Bin Separately by Filter"""
    for filter in flux_by_filter:
        flux = 0
        flux_unc = 0
        w_tot = 0
        idx = 0
        for point in flux_by_filter[filter]:
            w = 1/((unc_by_filter[filter][idx])**2)
            flux += w*point
            w_tot += w
            idx += 1
        if w_tot != 0:
            flux = flux / w_tot
            flux_unc = w_tot**(-1/2) 
        if len(flux_by_filter[filter]) != 0: 
            if filter in combined_flux:
                combined_flux[filter].append(flux)
                combined_unc[filter].append(flux_unc)
                combined_start[filter].append(jd[start])
                combined_end[filter].append(jd[end])
            else:
                combined_flux[filter] = [flux]
                combined_unc[filter] = [flux_unc]
                combined_start[filter] = [jd[start]]
                combined_end[filter] = [jd[end]]


def cal_mag(flux, flux_unc, jd_start, jd_end, ra, dec, num_days):
    """Obtaining calibrated magnitudes (for transients)."""
    zpavg = min(zpdiff.values())
    mag = 0
    sigma = 0
    for filter in flux:
        i = 0
        while (i < len(flux[filter])):
            if ((flux[filter][i] / flux_unc[filter][i]) > 5):  # 5 is the signal to noise threshold for declaring a measurement a "non-detection", so that it can be assigned an upper-limit (see Masci et. al)
                # confident detection, plot magnitude with error bars
                if flux[filter][i] < 0:
                    mag = -(zpavg-2.5*math.log10(-flux[filter][i]))  # negative flux cannot be plotted using log10
                else:
                    mag = (zpavg - 2.5*math.log10(flux[filter][i]))  # plotted as points
                    sigma = 1.0857 * flux_unc[filter][i] / flux[filter][i]
                if filter == 'ZTF_g':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='ZTF_g', c='blue')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='blue')
                elif filter == 'ZTF_r':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='ZTF_r', c='red')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='red')
                else:  # ZTF_i
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='ZTF_i', c='green')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='green')
            else:
                # compute upper flux limits and plot as arrow
                mag = (zpavg - 2.5*math.log10(3*flux_unc[filter][i]))  # 3 is the actual signal to noise ratio to use when computing SNU-sigma upper-limit
                if filter == 'ZTF_g':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='blue', label='Single upper-epoch limits')  # plot as arrow
                elif filter == "ZTF_r":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='red', label='Single upper-epoch limits')
                else:  # ZTF_i
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='green', label='Single upper-epoch limits')
            i += 1
    plt.xlabel('jd')
    plt.ylabel('magnitude')
    plt.title("RA: " + ra + "\nDEC: " + dec + "\nDays Binned: " + str(num_days))  # will add title 
    plt.legend(labels=['ZTF_g', 'ZTF_r', 'ZTF_i'])
    leg = plt.gca().get_legend()
    leg.legend_handles
    leg.legend_handles[0].set_color('blue')
    leg.legend_handles[1].set_color('red')
    leg.legend_handles[2].set_color('green')
    plt.gca().invert_yaxis()
    plt.show()


def get_indices(start, i, num_days):
    """Creates a dictionary of indices {start, end} by day."""
    end = start
    if i >= len(jd):
        return None
    while ((start in jd and i+1 in jd) and (round(jd[start] + num_days, 7) >= round(jd[i+1], 7))):
        end = i+1
        i += 1
        if i == len(jd) - 1:
            break
    windows[start] = end
    get_indices(end + 1, i+1, num_days)
    

def main():
    # use of user inputs to select desired binning window
    file_name = str(input("Enter filename: "))
    start = 57  # starting line 57 appears to be consistent across ZTF files, but may need to be changed in future updates; note: if using a file that is not of ZTF format, indices on lines 20 and 21 of this code will need to be edited to match your file
    end = 58 + int(input("Enter ending index: "))  # readlines method appears to be end exclusive, hence 57 + 1; 
    baseline = float(input("If there is any residual baseline (nonzero), input it here. Else, input 0: "))  # you will need to have examined a plot of forcediffimflux to jd to determine if there is any residual offset in the baseline (Masci et. al section 10)
    num_days = float(input("Enter the number of days to be binned at a time: "))
    with open(file_name) as f: # opens .txt file
        data = f.readlines()[start: end]
    sys.setrecursionlimit(len(data))
    fill_vars(data)
    f.close()
    with open(file_name) as f:
        data = f.readlines()[3:4]
        for line in data:
            ra = str(line.strip())[25:]
    f.close()
    with open(file_name) as f:
        data = f.readlines()[4:5]
        for line in data:
            dec = str(line.strip())[25:]
    if baseline != 0:
        correct_baseline()
    validate_uncertainties()  # if file used is already uncertainty validated, you may remove this call
    get_indices(0, 0, num_days)
    for start in windows:
        rescale(start, windows[start])
        collapse_flux_by_filter(start, windows[start])
    cal_mag(combined_flux, combined_unc, combined_start, combined_end, ra, dec, num_days)
    
if __name__ == '__main__':  # invoke python main.py to run
    main()   