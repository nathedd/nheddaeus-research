#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: main.py
Author: Nat Heddaeus
Date: 2024-12-30
Version: 1.0
Description: Script for running lc_stack for multiple windows and plotting them.
Contact: nathedd@unc.edu
"""


from lc_stack import main, cal_mag, zpdiff, combined_meas


g_flux = []
r_flux = []
i_flux = []

g_unc = []
r_unc = []
i_unc = []

g_jd_start = []
g_jd_end = []
r_jd_start = []
r_jd_end = []
i_jd_start = []
i_jd_end = []


runtime = int(input("Input number of windows to bin: "))
idx = 0
while (idx < runtime):
    idx += 1
    main()
    g_flux.append(combined_meas['ZTF_g'][0])
    g_unc.append(combined_meas['ZTF_g'][1])
    g_jd_start.append(combined_meas['ZTF_g'][2])
    g_jd_end.append(combined_meas['ZTF_g'][3])

    r_flux.append(combined_meas['ZTF_r'][0])
    r_unc.append(combined_meas['ZTF_r'][1])
    r_jd_start.append(combined_meas['ZTF_r'][2])
    r_jd_end.append(combined_meas['ZTF_r'][3])

    i_flux.append(combined_meas['ZTF_i'][0])
    i_unc.append(combined_meas['ZTF_i'][1])
    i_jd_start.append(combined_meas['ZTF_i'][2])
    i_jd_end.append(combined_meas['ZTF_i'][3])


cal_mag(g_flux, g_unc, 'ZTF_g', g_jd_start, g_jd_end)
cal_mag(r_flux, r_unc, 'ZTF_r', r_jd_start, r_jd_end)
cal_mag(i_flux, i_unc, 'ZTF_i', i_jd_start, i_jd_end)