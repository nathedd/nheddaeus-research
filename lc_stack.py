#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-01-05
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

combined_meas = {}  # {filter, [flux, unc, jd start, jd end]}

windows = {}  # for helper function get_indices


def fill_vars(data):
    "Takes a string of ascii data and converts it to dictionary variables. Note: If file used is of a different format than ZTF file, column indices will need to be changed to match your file."
    for line in data:
        line = line.strip()
        columns = line.split()
        jd[int(columns[0])] = (float(columns[22]))  # creates dictionary of {index, julian day}
        if (columns[24] != 'null'): 
            forcediffimflux[int(columns[0])] = (float(columns[24]))  # creates dictionary of {index, flux}
        else:
            forcediffimflux[int(columns[0])] = (None)  # catches unusable data points; fills with None for consistency between variables 
        if (columns[25] != 'null'): 
            forcediffimfluxunc[int(columns[0])] = (float(columns[25]))  # creates dictionary of {index, flux uncertainty}
        else:
            forcediffimflux[int(columns[0])] = (None)  # catches unusable data points; fills with None for consistency between variables
        if (columns[27] != 'null'):
            forceddiffimchisq[int(columns[0])] = (float(columns[27]))  # creates dictionary of {index, chi^2}
        else:
            forceddiffimchisq[int(columns[0])] = (None)  # catches unusable data points; fills with None for consistency between variables
        zpdiff[int(columns[0])] = (float(columns[20]))  # creates dictionary of {index, zero points}
        filter[int(columns[0])] = (columns[4])  # creates dictionary of {index, filter}


def correct_baseline():
    """Following an estimate of the baseline level, subtract estimate from differential flux measurements."""
    for index in forcediffimflux:
        if forcediffimflux[index] != None:
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
        if forceddiffimchisq[index] != None:
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
            if filter[index] == 'ZTF_g' and forcediffimfluxunc[index] != None:
                forcediffimfluxunc[index] = new
    if round(r_avg) != 1:
        if filter[index] == 'ZTF_r' and forcediffimfluxunc[index] != None:
                forcediffimfluxunc[index] = new
    if round(i_avg) != 1:
        if filter[index] == 'ZTF_i' and forcediffimfluxunc[index] != None:
                forcediffimfluxunc[index] = new


def rescale():  # start, end inclusive
    """Make new columns for forcediffimfluxi and forcediffimfluxunci with rescaled input fluxes and uncertainties and sort new fluxes and uncertainties by filter."""
    zpavg =  min(zpdiff.values())  # flux / uncertainty ratio is consistent for any value; picking min for simplicity
    g_list = []  # list of rescaled fluxes in the g band
    r_list = []  # list of rescaled fluxes in the r band
    i_list = []  # list of rescaled fluxes in the i band
    g_unc_list = []  # list of rescaled uncertainties in the g band
    r_unc_list = []  # list of rescaled uncertainties in the r band
    i_unc_list = []  # list of rescaled uncertainties in the i band
    for index in forcediffimflux:  
        if forcediffimflux[index] != None and forcediffimfluxunc[index] != None:  # skips unusable data points
            new = forcediffimflux[index]*10**(0.4*(zpavg-zpdiff[index]))  # place fluxes on the same photometric zeropoint
            forcediffimflux_new[index] = new
            # sorts new flux by filter and matches dictionary of jd to flux by filter by index
            if filter[index] == 'ZTF_g':
                g_list.append(new)
            elif filter[index] == 'ZTF_r':
                r_list.append(new)
            elif filter[index] == 'ZTF_i':
                i_list.append(new)
    flux_by_filter['ZTF_g'] = g_list
    flux_by_filter['ZTF_r'] = r_list
    flux_by_filter['ZTF_i'] = i_list
    for index in forcediffimfluxunc:
        if forcediffimfluxunc[index] != None and forcediffimflux[index] != None:  # skips unusable data points
            new = forcediffimfluxunc[index]*10**(0.4*(zpavg-zpdiff[index]))  # place uncertainties on the same photometric zeropoint
            forcediffimfluxunc_new[index] = new
            # sorts new uncertainty by filter
            if filter[index] == 'ZTF_g':
                g_unc_list.append(new)
            elif filter[index] == 'ZTF_r':
                r_unc_list.append(new)
            elif filter[index] == 'ZTF_i':
                i_unc_list.append(new)
    unc_by_filter['ZTF_g'] = g_unc_list
    unc_by_filter['ZTF_r'] = r_unc_list
    unc_by_filter['ZTF_i'] = i_unc_list


def collapse_flux_by_filter(start, end):
    """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average. Bin Separately by Filter"""
    for filter in flux_by_filter:
        flux = 0
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
            print(filter + ' flux: ' + str(flux) + ' flux_unc: ' + str(flux_unc)) 
        combined_meas[filter] = [flux, flux_unc, jd[start - 57], jd[end - 58]]  # if not using a ZTF file, have these match the indices used in jd dictionary


def cal_mag(flux, flux_unc, filter, jd_start, jd_end):
    """Obtaining calibrated magnitudes (for transients). Run in main.py file."""
    zpavg = min(zpdiff.values())
    mag = 0
    sigma = 0
    i = 0
    while (i < len(flux)):
        if ((flux[i] / flux_unc[i]) > 5):  # 5 is the signal to noise threshold for declaring a measurement a "non-detection", so that it can be assigned an upper-limit (see Masci et. al)
            # confident detection, plot magnitude with error bars
            if flux[i] < 0:
                mag = (-(zpavg-2.5*math.log10(-flux[i])))  # negative flux cannot be plotted using log10
            else:
                mag = (zpavg - 2.5*math.log10(flux[i]))  # plotted as points
                sigma = 1.0857 * flux_unc[i] / flux[i]
            plt.scatter((jd_end[i]+jd_start[i])/2, mag, label='Magnitude', c='blue')
            plt.errorbar((jd_end[i]+jd_start[i])/2, mag, yerr=sigma, ls='none', c='blue')
        else:
            # compute upper flux limits and plot as arrow
            mag = (zpavg - 2.5*math.log10(3*flux_unc[i]))  # 3 is the actual signal to noise ratio to use when computing SNU-sigma upper-limit
            plt.scatter((jd_end[i]+jd_start[i])/2, mag, marker='v', c='red', label='Single upper-epoch limits')  # plot as arrow
        i += 1
    plt.xlabel('jd')
    plt.ylabel('magnitude in '+ str(filter)[0:len(str(filter))])
    # plt.legend(loc = 'upper left')
    plt.show()


def get_indices(start, i):
    """Creates a dictionary of indices {start, end} by day."""
    end = start
    if i >= len(jd):
        print(windows)
        return None
    while int(jd[i] % 10) == int(jd[i+1] % 10):
        end = i+1
        i += 1
        if i == len(jd) - 1:
            break
    windows[start] = end
    get_indices(end + 1, i+1)
    

def main():
    # use of user inputs to select desired binning window
    file_name = str(input("Enter filename: "))
    start = 57 + int(input("Enter starting index (Index will be found in column one of your ascii file, in the rows that do not start with '#'): "))  # starting line 57 appears to be consistent across ZTF files, but may need to be changed in future updates; note: if using a file that is not of ZTF format, indices on lines 20 and 21 of this code will need to be edited to match your file
    end = 58 + int(input("Enter ending index: "))  # readlines method appears to be end exclusive, hence 57 + 1; 
    baseline = float(input("If there is any residual baseline (nonzero), input it here. Else, input 0: "))  # you will need to have examined a plot of forcediffimflux to jd to determine if there is any residual offset in the baseline (Masci et. al section 10)
    with open(file_name) as f: # opens .txt file
        data = f.readlines()[start: end]
    sys.setrecursionlimit(len(data))
    fill_vars(data)
    if baseline != 0:
        correct_baseline()
    validate_uncertainties()  # if file used is already uncertainty validated, you may remove this call
    rescale()
    collapse_flux_by_filter(start, end)


if __name__ == '__main__':  # invoke python main.py to run
    main()   