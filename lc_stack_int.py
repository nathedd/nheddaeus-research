#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack_int.py
Author: Nat Heddaeus
Date: 2024-04-23
Version: 1.0
Description: An optimized version of the code from lc_stack.py made for integration into the ztfrest pipeline. Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022. 

Contact: nathedd@unc.edu
"""

from astropy.io import ascii  # for reading files
from astropy.table import Table  # for outputting results
import numpy as np


def stack_lc(tbl, days_stack): 
    """Given a dataframe with a maxlike light curve, stack the flux."""
    snt_det=3  # signal to noise threshold for declaring a measurement a "non-detection"
    snt_ul=5  # actual signal to noise ratio for computing a sigma upper limit

    zpavg = np.mean(tbl['zpdiff,'])  # fiducial photometric zero point for rescaling fluxes

    length = len(tbl) - 1 # number of rows of data (zero indexed)

    t_out = Table([[],[],[],[],[],[],[],[],[],[]],
                  names=('jd', 'flux', 'flux_unc', 'zp', 'ezp',
                         'mag', 'mag_unc', 'limmag', 'filter', 'programid'),
                  dtype=('double', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'S', 'int'))  # used to output the stacked flux under new windows
    
    # make bins for stacking within inputted time windows
    bins = np.arange(tbl['jd,'][0], tbl['jd,'][length] + days_stack, days_stack)  # creates a bin for every day between the start and end date with mesh size of days_stack
    bin_len = len(bins)

    # correct baseline, validate uncertainties, rescale, collapse flux by filter, cal mag

    # place the fluxes on the same photometric zeropoint
    rs_flux = np.zeros(length+1)  # variable for rescaled fluxes
    rs_unc = np.zeros(length+1)  # variable for rescaled flux uncertainties

    for index in range (0, length):
        if (tbl['forcediffimflux,'][index]) != 'null':
            flux_i = np.copy(tbl['forcediffimflux,'][index]).astype(np.float64)  # old forcediffimflux value
            unc_i = np.copy(tbl['forcediffimfluxunc,'][index]).astype(np.float64)  # old forcediffimfluxunc value
            zpdiff_i = np.copy(tbl['zpdiff,'][index]).astype(np.float64)
            rs_flux[index] = flux_i*10**(0.4 * (zpavg-zpdiff_i))  # rescale the flux
            rs_unc[index] = unc_i*10**(0.4 * (zpavg-zpdiff_i))  # rescale the corresponding uncertainty
        else:
            continue

    # combine flux measurements by filter
    bin_flux = np.zeros(bin_len)
    filters = list(set(tbl['filter,']))  # fetches each unique filter





if __name__ == "__main__":
    filename1 = str(input("Input filename: "))  # this line can be replaced with a hardcoded file
    tbl = ascii.read(filename1, delimiter=' ', header_start=0)  # the script is hardcoded to work with ZTF files, references to header names may need to be changed for outside files
    stack_lc(tbl, 1.5)