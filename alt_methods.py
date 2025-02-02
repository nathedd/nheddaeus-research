#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: alt_methods.py
Author: Nat Heddaeus
Date: 2024-01-21
Version: 1.0
Description: Alternate methods to the ones designed for ZTF-formatted files in lc_stack.py. Meant for use with data used in "The Final Season Reimagined: 30 Tidal Disruption Events from the ZTF-I Survey" (Hammerstein et. al)

Contact: nathedd@unc.edu
"""

import sys
import math
import matplotlib.pyplot as plt

file_name = str(input("Enter filename: "))
num_days = int(input("Enter number of days to bin: "))
with open(file_name) as f:
    data = f.readlines()[1:]
sys.setrecursionlimit(len(data))

jd = []
filter = []
forcediffimflux = []
forcediffimfluxunc = []
windows = {}
flux_by_filter = {}
unc_by_filter = {}

# used for plotting
combined_flux = {}
combined_unc = {}
combined_start = {}
combined_end = {}


def hammerstein_vars():
    """Alternate method to fill_vars() in lc_stack.py."""
    for line in data:
        line = line.strip()
        columns = line.split()
        jd.append(float(columns[0]))
        filter.append(str(columns[1]))
        forcediffimflux.append(float(columns[2]))
        forcediffimfluxunc.append(float(columns[3]))


def hammerstein_windows(start, i, num_days):
    """Alternate method to get_indices() in lc_stack.py"""
    end = start
    if i >= len(jd):
        return None
    while ((i + 1 <= len(jd)) and (round(jd[start] + num_days, 12) >= round(jd[i+1], 7))) and (jd[i-1] - jd[i] < 5):
        end = i+1
        i += 1
        if i == len(jd) - 1:
            break
    windows[start] = end
    hammerstein_windows(end + 1, i+1, num_days)


def hammerstein_by_filter(start, end):
    """Sort flux by filters."""
    r_list = []
    g_list = []
    UVW2_list = []  # UVW2.uvot
    UVW1_list = []  # UVW1.uvot
    U_list = []  # U.uvot
    B_list = []  # B.uvot
    V_list = []  # V.uvot
    o_list = []  # o.atlas
    c_list = []  # c.atlas

    r_unc = []
    g_unc = []
    UVW2_unc = []  # UVW2.uvot
    UVW1_unc = []  # UVW1.uvot
    U_unc = []  # U.uvot
    B_unc = []  # B.uvot
    V_unc = []  # V.uvot
    o_unc = []  # o.atlas
    c_unc = []  # c.atlas

    for index in range(start, end+1):
        if filter[index] == "r.ztf":
            r_list.append(forcediffimflux[index])
            r_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "g.ztf":
            g_list.append(forcediffimflux[index])
            g_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "UVW2.uvot":
            UVW2_list.append(forcediffimflux[index])
            UVW2_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "UVW1.uvot":
            UVW1_list.append(forcediffimflux[index])
            UVW1_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "U.uvot":
            U_list.append(forcediffimflux[index])
            U_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "B.uvot":
            B_list.append(forcediffimflux[index])
            B_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "V.uvot":
            V_list.append(forcediffimflux[index])
            V_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "o.atlas":
            o_list.append(forcediffimflux[index])
            o_unc.append(forcediffimfluxunc[index])
        elif filter[index] == "c.atlas":
            c_list.append(forcediffimflux[index])
            c_unc.append(forcediffimfluxunc[index])

    flux_by_filter['g.ztf'] = g_list
    flux_by_filter['r.ztf'] = r_list
    flux_by_filter['UVW2.uvot'] = UVW2_list
    flux_by_filter['UVW1.uvot'] = UVW1_list
    flux_by_filter['U.uvot'] = U_list
    flux_by_filter['B.uvot'] = B_list
    flux_by_filter['V.uvot'] = V_list
    flux_by_filter['o.atlas'] = o_list
    flux_by_filter['c.atlas'] = c_list

    unc_by_filter['g.ztf'] = g_unc
    unc_by_filter['r.ztf'] = r_unc
    unc_by_filter['UVW2.uvot'] = UVW2_unc
    unc_by_filter['UVW1.uvot'] = UVW1_unc
    unc_by_filter['U.uvot'] = U_unc
    unc_by_filter['B.uvot'] = B_unc
    unc_by_filter['V.uvot'] = V_unc
    unc_by_filter['o.atlas'] = o_unc
    unc_by_filter['c.atlas'] = c_unc



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


def hammerstein_cal_mag(flux, flux_unc, jd_start, jd_end, num_days):
    """Alternate method to cal_mag() in lc_stack.py."""
    zpavg = 25
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
                if filter == 'g.ztf':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='ZTF_g', c='blue')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='blue')
                elif filter == 'r.ztf':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='ZTF_r', c='red')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='red')
                elif filter == "UVW2.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='UVW2.uvot', c='green')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='green')
                elif filter == "UVW1.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='UVW1.uvot', c='gold')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='gold')
                elif filter == "U.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='U.uvot', c='skyblue')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='skyblue')
                elif filter == "B.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='B.uvot', c='purple')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='purple')
                elif filter == "V.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, label='V.uvot', c='darkolivegreen')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='darkolivegreen')
                elif filter =="o.atlas":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='chocolate', label='Single upper-epoch limits')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='chocolate')
                elif filter =="c.atlas":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='tan', label='Single upper-epoch limits')
                    plt.errorbar((jd_end[filter][i]+jd_start[filter][i])/2, mag, yerr=sigma, ls='none', c='tan')
            else:
                # compute upper flux limits and plot as arrow
                mag = (zpavg - 2.5*math.log10(3*flux_unc[filter][i]))  # 3 is the actual signal to noise ratio to use when computing SNU-sigma upper-limit
                if filter == 'g.ztf':
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='blue', label='Single upper-epoch limits')  # plot as arrow
                elif filter == "r.ztf":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='red', label='Single upper-epoch limits')
                elif filter =="UVW2.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='green', label='Single upper-epoch limits')
                elif filter =="UVW1.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='gold', label='Single upper-epoch limits')
                elif filter =="U.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='skyblue', label='Single upper-epoch limits')
                elif filter =="B.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='purple', label='Single upper-epoch limits')
                elif filter =="V.uvot":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='darkolivegreen', label='Single upper-epoch limits')
                elif filter =="o.atlas":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='chocolate', label='Single upper-epoch limits')
                elif filter =="c.atlas":
                    plt.scatter((jd_end[filter][i]+jd_start[filter][i])/2, mag, marker='v', c='tan', label='Single upper-epoch limits')
            i += 1
    plt.xlabel('jd')
    plt.ylabel('magnitude')
    plt.title("Days Binned: " + str(num_days))  # will add title 
    plt.legend(labels=['ZTF_g', 'ZTF_r', 'ZTF_i', 'UVW2.uvot', 'UVW1.uvot', 'U.uvot', 'B.uvot', 'V.uvot', 'o.atlas', 'c.atlas'])
    leg = plt.gca().get_legend()
    leg.legend_handles[1].set_color('red')
    leg.legend_handles[2].set_color('green')
    leg.legend_handles[3].set_color('gold')
    leg.legend_handles[4].set_color('skyblue')
    leg.legend_handles[5].set_color('purple')
    leg.legend_handles[6].set_color('darkolivegreen')
    leg.legend_handles[7].set_color('chocolate')
    leg.legend_handles[8].set_color('tan')
    plt.gca().invert_yaxis()
    plt.show()



hammerstein_vars()
hammerstein_windows(0, 0, num_days)
for start in windows:
    hammerstein_by_filter(start, windows[start])
    collapse_flux_by_filter(start, windows[start])
hammerstein_cal_mag(combined_flux, combined_unc, combined_start, combined_end, num_days)
