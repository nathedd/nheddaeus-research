#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack_int.py
Author: Nat Heddaeus
Date: 2024-05-11
Version: 1.0
Description: An optimized version of the code from lc_stack.py made for integration into the ztfrest pipeline. Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022. 

Contact: nathedd@unc.edu
"""

from astropy.io import ascii  # for reading files
from astropy.table import Table  # for outputting results
import numpy as np
import matplotlib.pyplot as plt


def stack_lc(tbl, days_stack): 
    """Given a dataframe with a maxlike light curve, stack the flux."""
    snt_det=3  # signal to noise threshold for declaring a measurement a "non-detection"
    snt_ul=5  # actual signal to noise ratio for computing a sigma upper limit

    zpavg = np.mean(tbl['zpdiff,'])  # fiducial photometric zero point for rescaling fluxes

    length = len(tbl) - 1 # number of rows of data (zero indexed)

    t_out = Table([[],[],[],[],[],[],[],[]],
                  names=('jd', 'flux', 'flux_unc', 'zp',
                         'mag', 'mag_unc', 'flux_ul', 'filter'),
                  dtype=('double', 'f', 'f', 'f', 'f', 'f', 'f', 'S'))  # used to output the stacked flux under new windows
    
    # make bins for stacking within inputted time windows
    jd = np.array(tbl['jd,'])
    fil = tbl['filter,']  # filter by index
    bins = np.arange(jd[0] + days_stack, jd[length] + days_stack, days_stack)  # creates a bin for every day between the start and end date with mesh size of days_stack
    bin_len = len(bins)


    # place the fluxes on the same photometric zeropoint
    rs_flux = np.zeros(length+1)  # variable for rescaled fluxes
    rs_unc = np.zeros(length+1)  # variable for rescaled flux uncertainties

    for index in range (0, length + 1):
        if (tbl['forcediffimflux,'][index]) != 'null':
            flux_i = np.copy(tbl['forcediffimflux,'][index]).astype(np.float64)  # old forcediffimflux value
            unc_i = np.copy(tbl['forcediffimfluxunc,'][index]).astype(np.float64)  # old forcediffimfluxunc value
            zpdiff_i = np.copy(tbl['zpdiff,'][index]).astype(np.float64)
            rs_flux[index] = flux_i*10**(0.4 * (zpavg-zpdiff_i))  # rescale the flux
            rs_unc[index] = unc_i*10**(0.4 * (zpavg-zpdiff_i))  # rescale the corresponding uncertainty
        else:
            rs_flux[index] = np.nan
            rs_unc[index] = np.nan

    # combine flux measurements by filter
    filters = list(set(fil))  # fetches each unique filter
    bin_flux = np.zeros((bin_len, 2))
    bin_unc = np.zeros((bin_len, 2))
    for i in range(0, len(filters)):
        f_idx = np.where(fil==filters[i])[0]  # get all indices for a particular filter
        for j in range(0, bin_len):
            new_idx = np.copy(f_idx)
            bin_idx = np.where(jd[new_idx] <= bins[j])[0]  # get valid lower limit indices in jd_bin
            if bin_idx.size != 0:
                new_idx = new_idx[bin_idx]
            if j != 0 and new_idx.size != 0:
                bin_idx = np.where(jd[new_idx] > bins[j-1])[0]  # get valid upper limit indices in jd_bin
                new_idx = new_idx[bin_idx]
            w = 1/(rs_unc[new_idx])**2
            w_sum = np.nansum(w)  # prevents RuntimeWarning with scalar divide/scalar power
            if w.size > 0 and w_sum != 0:
                bin_flux[j, 0] = np.nansum(w*rs_flux[new_idx])/w_sum
                bin_unc[j, 0] = np.nansum(w)**(-1/2)  # combined unc
                bin_flux[j, 1] = i  # filter
                bin_unc[j, 1] = i  # filter
    
    # calculate calibrated magnitudes
    mag = np.zeros(bin_len)
    sigma = np.zeros(bin_len)
    flux_ul = np.zeros(bin_len)
    for i in range (0, bin_len):
        flux_i = bin_flux[i, 0]
        unc_i = bin_unc[i, 0]
        filt_i = int(bin_flux[i, 1])  # index in filters
        if flux_i != 0:
            if (flux_i/unc_i) > snt_det:  # confident detection
                mag[i] = zpavg - 2.5*np.log10(flux_i)  # compute AB magnitude
                sigma[i] = 1.0857*unc_i/flux_i  # compute uncertainty in magnitude
                flux_ul[i] = np.nan
            else:
                flux_ul[i] = zpavg - 2.5*np.log10(snt_ul*unc_i)  # compute flux upper limit
                mag[i] = np.nan
                sigma[i] = np.nan
        else:
            flux_ul[i] = np.nan
            mag[i] = np.nan
            sigma[i] = np.nan
        
        # fill in output table
        t_out.add_row([bins[i], flux_i, unc_i, zpavg, mag[i], sigma[i], flux_ul[i], filters[filt_i]])
    return t_out


def plot_lc(t_out):
    """Make a light curve of the binned fluxes."""
    jd = np.array(t_out['jd'])
    mag = np.array(t_out['mag'])
    sigma = np.array(t_out['mag_unc'])
    fil = np.array(t_out['filter'])  # filters start with b for some reason
    flux_ul = np.array(t_out['flux_ul'])

    filters = list(set(fil))
    colors = ['red', 'blue', 'green', 'purple', 'gray', 'olive', 'cyan', 'pink', 'brown', 'orange']
    for i in range(0, len(filters)):
        idx = np.where(fil == filters[i])[0]
        mags = mag[idx]
        jds = jd[idx]
        sigmas = sigma[idx]
        uls = flux_ul[idx]

        plt.scatter(jds, mags, c=colors[i])
        plt.errorbar(jds, mags, yerr=sigmas, ls='none', c=colors[i])
        plt.scatter(jds, uls, marker='v', color=colors[i])

    plt.gca().invert_yaxis()
    plt.xlabel('jd')
    plt.ylabel('AB magnitude')
    plt.show()
                

if __name__ == "__main__":
    filename1 = str(input("Input filename: "))  # this line can be replaced with a hardcoded file
    tbl = ascii.read(filename1, delimiter=' ', header_start=0)  # the script is hardcoded to work with ZTF files, references to header names may need to be changed for outside files
    t_out = stack_lc(tbl, 1.5)
    plot_lc(t_out)