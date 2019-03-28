from __future__ import division
from sys import path
path.append('modules/')
from _curses import raw
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker
import matplotlib.pyplot as plt
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
import scivis.units as ut # for tmerg
import statsmodels.formula.api as smf
import scipy.optimize as opt
from math import pi, sqrt
import matplotlib as mpl
from glob import glob
import pandas as pd
import numpy as np
import itertools
import os.path
import cPickle
import time
import copy
import h5py
import csv
import os
import gc

from scidata.utils import locate
import scidata.carpet.hdf5 as h5
import scidata.xgraph as xg
from scidata.carpet.interp import Interpolator

from scipy import interpolate
cmap = plt.get_cmap("viridis")
# from sklearn.linear_model import LinearRegression-
from scipy.optimize import fmin
from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
    MultipleLocator
from matplotlib.colors import LogNorm, Normalize

from general import *
from lists import *
from filework import *
from units import time_constant, volume_constant, energy_constant

from math import pi, log10
import time



# class READ_CYLH5:
#
#     def __init__(self, fname):
#
#         self.fname = fname
#         self.file = h5py.File(fname, "r")
#
#         self.iteration = self.fname.split('/')[-1].split('.')[0]
#
#         file = h5py.File(fname, 'r')
#         print("Data File: {} [".format(self.iteration)),
#         for i in file:
#             print(i),
#         print(']'),
#
#         self.phi_3d     = np.array(self.file["phi_cyl_3d"])
#         self.r_3d       = np.array(self.file["r_cyl_3d"])
#         self.z_3d       = np.array(self.file["z_cyl_3d"])
#         self.dphi_3d    = np.array(self.file["dphi_cyl_3d"])
#         self.dr_3d      = np.array(self.file["dr_cyl_3d"])
#         self.dz_3d      = np.array(self.file["dz_cyl_3d"])
#         self.J          = np.array(self.file["J"])
#         self.Jflux      = np.array(self.file["Jflux"])
#         self.D          = np.array(self.file["D"])
#         self.vphi       = np.array(self.file["vphi"])
#         self.vr         = np.array(self.file["vr"])
#         self.Dgeo       = np.array(self.file["Dgeo"])
#         self.Dbern      = np.array(self.file["Dbern"])
#         self.rho        = np.array(self.file["rho"])
#         self.temp       = np.array(self.file["temp"])
#         self.ye         = np.array(self.file["Ye"])
#
#         self.time = self.file.attrs["time"] * 0.004925794970773136 * 1e-3
#         # units.conv_time(units.cactus, units.cgs, self.file.attrs["time"])
#         print(" | loaded it:{} time: {} ms".format(int(self.iteration), "%.1f" % (self.time * 1e3)))
#
#     def __delete__(self, instance):
#         del instance.phi_3d
#         del instance.r_3d
#         del instance.z_3d
#         del instance.dphi_3d
#         del instance.dr_3d
#         del instance.dz_3d
#         del instance.J
#         del instance.Jflux
#         del instance.D
#         del instance.vphi
#         del instance.vr
#         del instance.Dgeo
#         del instance.Dbern
#         del instance.rho
#         del instance.temp
#
#         instance.file.close()
#         gc.collect()
#
# class ANALYZE_CYLINDRICAL(READ_CYLH5):
#
#     def __init__(self, fname):
#         '''
#         LOADS AND ANALYZES THE .cylh5 file,
#         that contatins the CYLINDRICAL grid and interpolated variables.
#         :param fname:
#         '''
#
#         READ_CYLH5.__init__(self, fname)
#
#         # self.fname = fname
#         # self.file = h5py.File(fname, "r")
#         #
#         # self.iteration = self.fname.split('.')[0]
#         #
#         # file = h5py.File(fname, 'r')
#         # print("Data File: {} [".format(self.iteration)),
#         # for i in file:
#         #     print(i),
#         # print(']'),
#         #
#         # self.phi_3d  = np.array(self.file["phi_cyl_3d"]) # 3D
#         # self.r_3d    = np.array(self.file["r_cyl_3d"])   # 3D
#         # self.z_3d    = np.array(self.file["z_cyl_3d"])   # 3D
#         # self.dphi_3d = np.array(self.file["dphi_cyl_3d"])# 3D
#         # self.dr_3d   = np.array(self.file["dr_cyl_3d"])  # 3D
#         # self.dz_3d   = np.array(self.file["dz_cyl_3d"])  # 3D
#         # self.J       = np.array(self.file["J"])       # 3D
#         # self.Jflux   = np.array(self.file["Jflux"])   # 3D
#         # self.D       = np.array(self.file["D"])       # 3D
#         # self.vphi    = np.array(self.file["vphi"])    # 3D
#         # self.vr      = np.array(self.file["vr"])      # 3D
#         # self.Dgeo    = np.array(self.file["Dgeo"])    # 3D
#         # self.Dbern   = np.array(self.file["Dbern"])   # 3D
#         # self.rho     = np.array(self.file["rho"])     # 3D
#         # self.temp    = np.array(self.file["temp"])
#         #
#         # self.time = self.file.attrs["time"] * 0.004925794970773136 * 1e-3# units.conv_time(units.cactus, units.cgs, self.file.attrs["time"])
#         # print(" | time: {} ms".format("%.1f"%(self.time * 1e2)))
#
#     def get_j_z_integ_ave(self):
#
#         # J: integrate over 'z'; average over 'phi' for every 'r'
#         J_z_summed = 2 * np.sum(self.J * self.dz_3d, axis=(2))   # for every r, phi, summ over z
#         # J_z_phi_summed = np.sum(J_z_summed, axis=(1))   # for every r, summ over all phi
#
#         J_z_summed_ave = np.zeros(len(J_z_summed[0, :]))
#         for i in range (len(J_z_summed[:, 0])):
#             J_min_J_ave = J_z_summed[i, :] - (np.sum(J_z_summed[i, :]) / len(J_z_summed[i, :]))
#             J_min_J_ave *= self.r_3d[i, 0, 0] ** (.5)
#             J_z_summed_ave = np.vstack((J_z_summed_ave, J_min_J_ave))
#         J_z_summed_ave = np.delete(J_z_summed_ave, 0, 0)
#
#         return J_z_summed_ave
#
#     def get_jf_z_integ(self):
#
#         # J: integrate over 'z'; average over 'phi' for every 'r'
#         Jf_z_summed = 2 * np.sum(self.Jflux * self.dz_3d, axis=(2))  # for every r, phi, summ over z
#         # J_z_phi_summed = np.sum(J_z_summed, axis=(1))   # for every r, summ over all phi
#
#         return Jf_z_summed
#
#     def get_d_z_integ(self):
#         return 2 * np.sum(self.D * self.dz_3d, axis=(2))
#
#     def get_dbern_z_integ(self):
#         return 2 * np.sum(self.Dbern * self.dz_3d, axis=(2))
#
#     def get_vphi_z0_ave(self):
#
#         vphi_z_ave = self.vphi[:, :, 0]  # np.sum(vphi_cyl, axis=2) / len(vphi_cyl[0, 0, :])
#         vphi_minus_phi_ave = np.zeros(len(self.vphi[0, :, 0]))
#
#         for i in range(len(self.vphi[:, 0, 0])):
#             vphi_minus_phi_ave_i = vphi_z_ave[i, :] - (np.sum(vphi_z_ave[i, :]) / len(vphi_z_ave[i, :]))
#             vphi_minus_phi_ave = np.vstack((vphi_minus_phi_ave, vphi_minus_phi_ave_i))
#
#         vphi_minus_phi_ave = np.delete(vphi_minus_phi_ave, 0, 0)
#         return vphi_minus_phi_ave
#
#     def get_vr_z0_ave(self):
#
#         # velocities: Get z=0 slice ; v - sum_phi(v)/len(v) for every 'r'
#         vr_z_ave = self.vr[:, :, 0]  # np.sum(vr_cyl, axis=2) / len(vr_cyl[0, 0, :])
#         vr_minus_phi_ave = np.zeros(len(self.vr[0, :, 0]))
#
#         for i in range(len(self.vphi[:, 0, 0])):
#             vr_minus_phi_ave_i = vr_z_ave[i, :] - (np.sum(vr_z_ave[i, :]) / len(vr_z_ave[i, :]))
#             vr_minus_phi_ave = np.vstack((vr_minus_phi_ave, vr_minus_phi_ave_i))
#
#         vr_minus_phi_ave = np.delete(vr_minus_phi_ave, 0, 0)
#         return vr_minus_phi_ave
#
#     def get_temperature(self):
#
#         return self.temp[:, :, 0]
#
#     def get_ye(self):
#
#         return self.ye[:, :, 0]
#
#     def get_dens_decomp(self, dens_2d, phi_2d, dphi_2d, dr_2d, m=1):
#         '''
#         Returns complex arrays [\int(d\phi)] and [\int(dr d\phi)]
#         '''
#
#         integ_over_phi = np.sum(dens_2d * np.exp(1j * m * phi_2d) * dphi_2d, axis=1)
#
#         integ_over_phi_r = np.sum(integ_over_phi * dr_2d[:, 0])
#
#         return integ_over_phi, integ_over_phi_r
#
# class PLOT_CYLINDRICAL:
#
#     def __init__(self, fnames, polar=True):
#
#         # settings
#         self.rmin = 0 # min R
#         self.rmax = 50 # max R
#         self.jmin = -5e-4 # min J (plot)
#         self.jmax = 5e-4 # max J (plot)
#         self.jfmin = -5e-6
#         self.jfmax = 5e-6
#         self.dmin = 1e-7 # min Dunb (plot)
#         self.dmax = 6e-7 # max Dunb (plot)
#         self.dmask_below = None#2e-7 # mask Dunb below this value
#         self.dmask_above = None # mask Dunb above this value
#         self.dscale = 'log'
#         self.tmin = 5.
#         self.tmax = 10.
#         self.tcmap = 'Oranges'
#         self.tmask_above = None
#         self.tmask_below = None
#         self.yemin = 0.1
#         self.yemax = 0.3
#         self.yemask_below = None
#         self.yemask_above = None
#         self.tscale = 'log'
#         self.jcmap = 'RdBu_r' # cmap fpr J plot
#         self.dcmap = 'Greens' # cmap for Dunb plot
#         self.yecmap = 'Oranges'
#         self.v_i = 6 # plot ever v_i velocity (arrow plot)
#         self.vrmin = 15
#         self.vphi_eq_vr = True # fix vr = vphi (no angular component plotted)
#         self.vwidth = 0.008
#
#         self.dens_levels = [1e-5, 1e-4]
#
#
#         self.set_save_format = '.pdf'
#
#         if len(fnames) < 1: raise IOError("No files selected for plotting")
#
#         self.polar=polar
#
#         self.cl_cyl = []
#         for fname in fnames:
#             self.cl_cyl.append(ANALYZE_CYLINDRICAL(fname))
#
#         if len(self.cl_cyl) < 1: raise IOError("No files loaded, no classes appended")
#
#     def plot_dens_decomp_phase(self, ax, cl_cyl, norm_to_m0=True, d_or_rho='d'):
#
#         # choose what to integrate
#         if d_or_rho == 'd':
#             y = cl_cyl.D[:, :, 0]
#         else:
#             y = cl_cyl.rho[:, :, 0]
#
#         # integrate and normalize
#         if norm_to_m0:
#             m1_integ_over_phi, m1_integ_over_phi_r \
#                 = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=1)
#             m0_integ_over_phi, m0_integ_over_phi_r \
#                 = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=0)
#
#             m1_integ_over_phi_r = m1_integ_over_phi_r / m0_integ_over_phi_r
#             m1_integ_over_phi = m1_integ_over_phi / m0_integ_over_phi
#         else:
#             m1_integ_over_phi, m1_integ_over_phi_r \
#                 = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=1)
#
#         m1_integ_over_phi_abs = np.absolute(m1_integ_over_phi)
#         m1_integ_over_phi_phase = np.angle(m1_integ_over_phi, deg=True)
#
#         m1_integ_over_phi_r_abs = np.absolute(m1_integ_over_phi_r)
#         m1_integ_over_phi_r_phase = np.angle(m1_integ_over_phi_r, deg=True)
#
#
#         # plotting
#         r = np.array(cl_cyl.r_3d[:, 0, 0])
#         phi = np.zeros(r.shape)
#         phi.fill(m1_integ_over_phi_r_phase)
#         # print(phi)
#         ax.plot(np.deg2rad(phi), r, '-', color='black') # plot integrated
#         # ax.plot(np.deg2rad(m1_integ_over_phi_phase), r, '-', color='black') # plot for every 'r'
#
#
#         # fig, ax = plt.subplots()
#         # # ax.scatter(arr.real, arr.imag)
#         # ax.set_xlabel('Re')
#         # ax.set_ylabel(r'Absolute ($C_{m}=\int_0 ^{2\pi} \rho W \sqrt{-g}\cdot \exp(i m \phi) d\phi$ )')
#         # ax.set_xlim(0, 50)
#         # # ax.set_ylim(-0.001, 0.001)
#         #
#         # plt.plot(cl_cyl.r_3d[:, 0, 0], np.absolute(arr), '.', color='black')
#         # plt.plot(cl_cyl.r_3d[:, 0, 0], np.absolute(arr), '-', color='black')
#         # plt.savefig('{}.png'.format('tst_reim'), bbox_inches='tight', dpi=128)
#         # plt.close()
#
#     def plot_dens_decomp_phase2(self, ax, cyl_cl, m_list=[], norm_to_m0=True, int_over_z=True, d_or_rho='d'):
#
#         for i_m, m in enumerate(m_list):
#             if d_or_rho == 'd':
#                 using = cyl_cl.D
#             else:
#                 using = cyl_cl.rho
#             if int_over_z:
#                 m_int_phi, m_int_phi_r = \
#                     PHYSICS.get_dens_decomp_3d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, cyl_cl.dz_3d, m)
#                 if norm_to_m0:
#                     _int_phi, _int_phi_r = \
#                         PHYSICS.get_dens_decomp_3d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, cyl_cl.dz_3d, 0)
#                     m_int_phi /= _int_phi_r
#                     m_int_phi_r /= _int_phi_r
#             else:
#                 m_int_phi, m_int_phi_r = \
#                     PHYSICS.get_dens_decomp_2d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, m)
#                 if norm_to_m0:
#                     _int_phi, _int_phi_r = \
#                         PHYSICS.get_dens_decomp_2d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, 0)
#                     m_int_phi /= _int_phi_r
#                     m_int_phi_r /= _int_phi_r
#
#             if m == 1:
#                 r = np.array(cyl_cl.r_3d[:, 0, 0])
#                 r_min, r_max = ax.get_ylim()
#                 # for every 'r'
#                 ax.plot(np.angle(m_int_phi)[r<r_max], r[r<r_max], '-.', color='black')
#                 # integrated over 'r'
#                 phi = np.zeros(r.shape)
#                 phi.fill(np.angle(m_int_phi_r))
#                 ax.plot(phi[r<r_max], r[r<r_max], '-', color='black')  # plot integrated
#
#         #
#         #
#         # for m in m_list:
#         #
#         #     group = dfile["rl=%d" % m]
#         #     int_phi2d = np.array(group["int_phi"])
#         #     int_phi_r1d = np.array(group["int_phi_r"])  # | time   <-> r
#         #     times = np.array(group["time"])
#         #     r = np.array(group["r"])
#         #
#         #
#         #
#         #     if m == 1: color = 'black'
#         #     if m == 2: color = 'gray'
#         #     ax.plot(np.angle(int_phi2d[0, :]), r, '-', color=color)
#         #
#         #
#         #
#         #
#         #
#         #
#         #
#         #
#         #
#         # # choose what to integrate
#         # if d_or_rho == 'd':
#         #     y = cl_cyl.D[:, :, 0]
#         # else:
#         #     y = cl_cyl.rho[:, :, 0]
#         #
#         # # integrate and normalize
#         # if norm_to_m0:
#         #     m1_integ_over_phi, m1_integ_over_phi_r \
#         #         = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=1)
#         #     m0_integ_over_phi, m0_integ_over_phi_r \
#         #         = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=0)
#         #
#         #     m1_integ_over_phi_r = m1_integ_over_phi_r / m0_integ_over_phi_r
#         #     m1_integ_over_phi = m1_integ_over_phi / m0_integ_over_phi
#         # else:
#         #     m1_integ_over_phi, m1_integ_over_phi_r \
#         #         = cl_cyl.get_dens_decomp(y, cl_cyl.phi_3d[:, :, 0], cl_cyl.dphi_3d[:, :, 0], cl_cyl.dr_3d[:, :, 0], m=1)
#         #
#         # m1_integ_over_phi_abs = np.absolute(m1_integ_over_phi)
#         # m1_integ_over_phi_phase = np.angle(m1_integ_over_phi, deg=True)
#         #
#         # m1_integ_over_phi_r_abs = np.absolute(m1_integ_over_phi_r)
#         # m1_integ_over_phi_r_phase = np.angle(m1_integ_over_phi_r, deg=True)
#         #
#         #
#         # # plotting
#         # r = np.array(cl_cyl.r_3d[:, 0, 0])
#         # phi = np.zeros(r.shape)
#         # phi.fill(m1_integ_over_phi_r_phase)
#         # # print(phi)
#         # ax.plot(np.deg2rad(phi), r, '-', color='black') # plot integrated
#         # ax.plot(np.deg2rad(m1_integ_over_phi_phase), r, '-', color='black') # plot for every 'r'
#
#
#         # fig, ax = plt.subplots()
#         # # ax.scatter(arr.real, arr.imag)
#         # ax.set_xlabel('Re')
#         # ax.set_ylabel(r'Absolute ($C_{m}=\int_0 ^{2\pi} \rho W \sqrt{-g}\cdot \exp(i m \phi) d\phi$ )')
#         # ax.set_xlim(0, 50)
#         # # ax.set_ylim(-0.001, 0.001)
#         #
#         # plt.plot(cl_cyl.r_3d[:, 0, 0], np.absolute(arr), '.', color='black')
#         # plt.plot(cl_cyl.r_3d[:, 0, 0], np.absolute(arr), '-', color='black')
#         # plt.savefig('{}.png'.format('tst_reim'), bbox_inches='tight', dpi=128)
#         # plt.close()
#
#     def plot_J(self, ax, cl_cyl):
#
#         J = cl_cyl.get_j_z_integ_ave()
#         im = ax.pcolormesh(cl_cyl.phi_3d[:, :, 0], cl_cyl.r_3d[:, :, 0], J, cmap=self.jcmap, vmin=self.jmin, vmax=self.jmax)
#         im.set_rasterized(True)
#         return im
#
#     def plot_Jflux(self, ax, cl_cyl):
#
#         Jflux = cl_cyl.get_jf_z_integ()
#         im = ax.pcolormesh(cl_cyl.phi_3d[:, :, 0], cl_cyl.r_3d[:, :, 0], Jflux, cmap=self.jcmap, vmin=self.jfmin, vmax=self.jfmax)
#         im.set_rasterized(True)
#         return im
#
#     def plot_Dunb(self, ax, cl_cyl):
#
#         Dunb = cl_cyl.get_dbern_z_integ()
#         phi = cl_cyl.phi_3d[:, :, 0]
#         r = cl_cyl.r_3d[:, :, 0]
#
#         if self.dmin == None:
#             self.dmin = Dunb.min()
#         if self.dmax == None:
#             self.dmax = Dunb.max()
#
#         norm = LogNorm(vmin=self.dmin, vmax=self.dmax)
#
#         if self.dscale == 'norm':
#             if self.dmask_above == None and self.dmask_below == None:
#                 im2 = ax.pcolormesh(phi, r, Dunb, vmin=self.dmin, vmax=self.dmax, cmap=self.dcmap)
#             elif self.dmask_below != None and self.dmask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(Dunb, Dunb < self.dmask_below),
#                                     vmin=self.dmin, vmax=self.dmax, cmap=self.dcmap)
#             elif self.dmask_below == None and self.dmask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(Dunb, Dunb > self.dmask_above),
#                                     vmin=self.dmin, vmax=self.dmax, cmap=self.dcmap)
#             else:
#                 im2 = ax.pcolormesh(cl_cyl.phi, cl_cyl.r,
#                                     np.ma.masked_array(Dunb, (Dunb > self.dmask_above) & (Dunb < self.dmask_below)),
#                                     vmin=self.dmax, vmax=self.dmax, cmap=self.dcmap)
#         else:
#             if self.dmask_above == None and self.dmask_below == None:
#                 im2 = ax.pcolormesh(phi, r, Dunb, norm=norm, cmap=self.dcmap)
#             elif self.dmask_below != None and self.dmask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(Dunb, Dunb < self.dmask_below),
#                                     norm=norm, cmap=self.dcmap)
#             elif self.dmask_below == None and self.dmask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(Dunb, Dunb > self.dmask_above),
#                                     norm=norm, cmap=self.dcmap)
#             else:
#                 im2 = ax.pcolormesh(phi, r,
#                                     np.ma.masked_array(Dunb, (Dunb > self.dmask_above) & (Dunb < self.dmask_below)),
#                                     norm=norm, cmap=self.dcmap)
#
#         im2.set_rasterized(True)  # dbern_int<1e-6
#         return im2
#
#     def plot_vel_polar(self, ax, cl_cyl):
#
#         # Phi, R = np.meshgrid(cl_cyl.phi_3d[:, :, 0], cl_cyl.r_3d[0])
#         Phi = cl_cyl.phi_3d[:, :, 0]
#         R = cl_cyl.r_3d[:, :, 0]
#
#         vphi = cl_cyl.get_vphi_z0_ave()
#         vr = cl_cyl.get_vr_z0_ave()
#
#         if self.vphi_eq_vr: vphi = vr
#
#         # for phi_, vphi_ in zip(Phi[:, 0], vphi[:, 0]):
#         #     for r_, vr_ in zip(R[0, :], vr[0, :]):
#         #         Q = ax.quiver(phi_,r_, vphi * np.cos(phi_),vr_ * np.sin(phi_),
#         #                       units='width', pivot='mid', scale=1.5,
#         #                       width=np.array(self.vwidth * (r_ / self.rmax)))
#         Phi = Phi[::self.v_i, ::self.v_i]
#         R = R[::self.v_i, ::self.v_i]
#         vphi = vphi[::self.v_i, ::self.v_i]
#         vr = vr[::self.v_i, ::self.v_i]
#         mask = R > self.vrmin
#
#         Q = ax.quiver(Phi[mask],
#                       R[mask],
#                       vphi[mask] * np.cos(Phi[mask]),
#                       vr[mask] * np.sin(Phi[mask]),
#                       units='width', pivot='mid', scale=1.5,
#                       width= self.vwidth)
#
#     def plot_temp(self, ax, cl_cyl):
#
#         temp = cl_cyl.get_temperature()
#         phi = cl_cyl.phi_3d[:, :, 0]
#         r = cl_cyl.r_3d[:, :, 0]
#
#         if self.tmin == None:
#             self.tmin = temp.min()
#         if self.tmax == None:
#             self.tmax = temp.max()
#         norm = LogNorm(vmin=self.tmin, vmax=self.tmax)
#
#         if self.tscale == 'norm':
#             if self.tmask_above == None and self.dmask_below == None:
#                 im2 = ax.pcolormesh(phi, r, temp, vmin=self.tmin, vmax=self.tmax, cmap=self.tcmap)
#             elif self.tmask_below != None and self.tmask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(temp, temp < self.tmask_below),
#                                     vmin=self.tmin, vmax=self.tmax, cmap=self.tcmap)
#             elif self.tmask_below == None and self.tmask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(temp, temp > self.tmask_above),
#                                     vmin=self.tmin, vmax=self.tmax, cmap=self.tcmap)
#             else:
#                 im2 = ax.pcolormesh(cl_cyl.phi, cl_cyl.r,
#                                     np.ma.masked_array(temp, (temp > self.tmask_above) & (temp < self.tmask_below)),
#                                     vmin=self.tmax, vmax=self.tmax, cmap=self.tcmap)
#         else:
#             if self.tmask_above == None and self.tmask_below == None:
#                 im2 = ax.pcolormesh(phi, r, temp, norm=norm, cmap=self.tcmap)
#             elif self.tmask_below != None and self.tmask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(temp, temp < self.tmask_below),
#                                     norm=norm, cmap=self.tcmap)
#             elif self.tmask_below == None and self.tmask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(temp, temp > self.tmask_above),
#                                     norm=norm, cmap=self.tcmap)
#             else:
#                 im2 = ax.pcolormesh(phi, r,
#                                     np.ma.masked_array(temp, (temp > self.tmask_above) & (temp < self.tmask_below)),
#                                     norm=norm, cmap=self.tcmap)
#
#         im2.set_rasterized(True)  # dbern_int<1e-6
#         return im2
#
#     def plot_ye(self, ax, cl_cyl):
#
#         ye = cl_cyl.get_ye()
#         phi = cl_cyl.phi_3d[:, :, 0]
#         r = cl_cyl.r_3d[:, :, 0]
#
#         if self.yemin == None:
#             self.yemin = ye.min()
#         if self.yemax == None:
#             self.yemax = ye.max()
#         norm = LogNorm(vmin=self.yemin, vmax=self.yemax)
#
#         if self.tscale == 'norm':
#             if self.yemask_above == None and self.dmask_below == None:
#                 im2 = ax.pcolormesh(phi, r, ye, vmin=self.yemin, vmax=self.yemax, cmap=self.yecmap)
#             elif self.yemask_below != None and self.yemask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(ye, ye < self.yemask_below),
#                                     vmin=self.yemin, vmax=self.yemax, cmap=self.yecmap)
#             elif self.yemask_below == None and self.yemask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(ye, ye > self.yemask_above),
#                                     vmin=self.yemin, vmax=self.yemax, cmap=self.yecmap)
#             else:
#                 im2 = ax.pcolormesh(cl_cyl.phi, cl_cyl.r,
#                                     np.ma.masked_array(ye, (ye > self.yemask_above) & (ye < self.yemask_below)),
#                                     vmin=self.yemax, vmax=self.yemax, cmap=self.yecmap)
#         else:
#             if self.yemask_above == None and self.yemask_below == None:
#                 im2 = ax.pcolormesh(phi, r, ye, norm=norm, cmap=self.yecmap)
#             elif self.yemask_below != None and self.yemask_above == None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(ye, ye < self.yemask_below),
#                                     norm=norm, cmap=self.yecmap)
#             elif self.yemask_below == None and self.yemask_above != None:
#                 im2 = ax.pcolormesh(phi, r, np.ma.masked_array(ye, ye > self.yemask_above),
#                                     norm=norm, cmap=self.yecmap)
#             else:
#                 im2 = ax.pcolormesh(phi, r,
#                                     np.ma.masked_array(ye, (ye > self.yemask_above) & (ye < self.yemask_below)),
#                                     norm=norm, cmap=self.yecmap)
#
#         im2.set_rasterized(True)  # dbern_int<1e-6
#         return im2
#
#     def plot_2_files(self, out_name='tst_js'):
#
#         from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
#             MultipleLocator
#         from mpl_toolkits.axes_grid1 import make_axes_locatable
#
#         rows = 2
#         cols = 2
#
#         fig, (ax_list) = plt.subplots(rows, cols, subplot_kw=dict(projection='polar'), figsize=(cols * 4.5, rows * 3.6))
#
#         if rows == 1 and cols == 1:
#             ax = ax_list
#
#             im1 = self.plot_J(ax, self.cl_cyl[0])
#             im2 = self.plot_Dunb(ax, self.cl_cyl[0])
#
#             self.plot_vel_polar(ax, self.cl_cyl[0])
#
#             ax.contour(self.cl_cyl[0].phi_3d[:,0,:], self.cl_cyl[0].r_3d[:,0,0], self.cl_cyl[0].get_d_z_integ(),
#                        [1e-5, 1e-4], colors="gray")
#
#             cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#             cbar = plt.colorbar(im1, cax=cax1)
#             cbar.ax.set_title(r"$J$")
#
#             cax = fig.add_axes([1.1, 0.1, 0.03, 0.8])
#             cbar2 = plt.colorbar(im2, cax=cax)
#             cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#
#             ax.set_rlim(self.rmin, self.rmax)
#
#         if rows == 1 and cols > 1:
#
#             im1_list = []
#             im2_list = []
#             for cl, ax in zip(self.cl_cyl, ax_list):
#                 im1 = self.plot_J(ax, cl)
#                 im2 = self.plot_Dunb(ax, cl)
#
#                 self.plot_vel_polar(ax, cl)
#
#                 ax.contour(cl.phi, cl.r, cl.get_d_z_integ(),
#                            [1e-5, 1e-4], colors="gray")
#                 ax.set_rlim(self.rmin, self.rmax)
#                 im1_list.append(im1)
#                 im2_list.append(im2)
#
#             cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#             cbar = plt.colorbar(im1_list[-1], cax=cax1)
#             cbar.ax.set_title(r"$J$")
#
#             cax = fig.add_axes([1.1, 0.1, 0.03, 0.8])
#             cbar2 = plt.colorbar(im2_list[-1], cax=cax)
#             cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#
#         if rows == 2 and cols == 2:
#
#             im1_list = []
#             im2_list = []
#             i = 0
#             for i1 in range(len(ax_list[:][0])):
#                 for i2 in range(len(ax_list[0][:])):
#                     ax = ax_list[i1][i2]
#                     cl = self.cl_cyl[i]
#
#                     im1 = self.plot_J(ax, cl)
#                     im2 = self.plot_Dunb(ax, cl)
#
#                     self.plot_vel_polar(ax, cl)
#
#                     ax.contour(cl.phi_3d[0,:,0], cl.r_3d[:,0,0], cl.get_d_z_integ(),
#                                [1e-5, 1e-4], colors="gray")
#                     ax.set_rlim(self.rmin, self.rmax)
#                     im1_list.append(im1)
#                     im2_list.append(im2)
#                     i+=1
#
#             cax1 = fig.add_axes([0.9, 0.5, 0.03, 0.35]) # [*left*, *bottom*, *width*, *height*]
#             cbar = plt.colorbar(im1_list[-1], cax=cax1, extend='both')
#             cbar.ax.set_title(r"$J$")
#
#             cax2 = fig.add_axes([0.9, 0.1, 0.03, 0.35]) # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im2_list[-1], cax=cax2, extend='both')
#             cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#
#         plt.subplots_adjust(hspace=0.2)
#         plt.subplots_adjust(wspace=-0.1)
#         # plt.tight_layout()
#         plt.savefig('{}.png'.format(out_name), bbox_inches='tight', dpi=128)
#         plt.close()
#
#     def test_plot_ang_mom_flux(self, out_name):
#
#         rows = 1
#         cols = 1
#
#         fig, (ax_list) = plt.subplots(rows, cols, subplot_kw=dict(projection='polar'), figsize=(cols * 4.5, rows * 3.6))
#
#         if rows == 1 and cols == 1:
#             ax = ax_list
#
#             im1 = self.plot_Jflux(ax, self.cl_cyl[0])
#             cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#             cbar = plt.colorbar(im1, cax=cax1, format='%.1e')
#             cbar.ax.set_title(r"$J$ flux")
#
#             im2 = self.plot_Dunb(ax, self.cl_cyl[0])
#             cax = fig.add_axes([1.1, 0.1, 0.03, 0.8])
#             cbar2 = plt.colorbar(im2, cax=cax)
#             cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#
#             # self.plot_vel_polar(ax, self.cl_cyl[0])
#
#             ax.contour(self.cl_cyl[0].phi_3d[:, :, 0], self.cl_cyl[0].r_3d[:, :, 0], self.cl_cyl[0].get_d_z_integ(),
#                        [1e-5, 1e-4], colors="gray")
#
#
#             self.plot_dens_decomp_phase(ax, self.cl_cyl[0], norm_to_m0=True, d_or_rho='d')
#
#             ax.set_rlim(self.rmin, self.rmax)
#
#
#         # plt.subplots_adjust(hspace=0.2)
#         # plt.subplots_adjust(wspace=-0.1)
#         # plt.tight_layout()
#         plt.savefig('{}.png'.format(out_name), bbox_inches='tight', dpi=128)
#         plt.close()
#
#     def plot_j_and_jflux_for_3_profiles(self, out_name):
#
#         rows = 4
#         cols = 4
#
#         fig, (ax_list) = plt.subplots(rows, cols, subplot_kw=dict(projection='polar'), figsize=(cols * 4.0, rows * 3.0))
#
#         # print(np.array(ax_list).shape)
#
#         # if rows == 1 and cols == 1:
#         #     ax = ax_list
#         #
#         #     im1 = self.plot_J(ax, self.cl_cyl[0])
#         #     im2 = self.plot_Dunb(ax, self.cl_cyl[0])
#         #
#         #     self.plot_vel_polar(ax, self.cl_cyl[0])
#         #
#         #     ax.contour(self.cl_cyl[0].phi_3d[:, 0, :], self.cl_cyl[0].r_3d[:, 0, 0], self.cl_cyl[0].get_d_z_integ(),
#         #                [1e-5, 1e-4], colors="gray")
#         #
#         #     cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#         #     cbar = plt.colorbar(im1, cax=cax1)
#         #     cbar.ax.set_title(r"$J$")
#         #
#         #     cax = fig.add_axes([1.1, 0.1, 0.03, 0.8])
#         #     cbar2 = plt.colorbar(im2, cax=cax)
#         #     cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#         #
#         #     ax.set_rlim(self.rmin, self.rmax)
#         #
#         # if rows == 1 and cols > 1:
#         #
#         #     im_j_list = []
#         #     im_dunb_list = []
#         #     for cl, ax in zip(self.cl_cyl, ax_list):
#         #         im1 = self.plot_J(ax, cl)
#         #         im2 = self.plot_Dunb(ax, cl)
#         #
#         #         self.plot_vel_polar(ax, cl)
#         #
#         #         ax.contour(cl.phi, cl.r, cl.get_d_z_integ(),
#         #                    [1e-5, 1e-4], colors="gray")
#         #         ax.set_rlim(self.rmin, self.rmax)
#         #         im_j_list.append(im1)
#         #         im_dunb_list.append(im2)
#         #
#         #     cax1 = fig.add_axes([0.9, 0.1, 0.03, 0.8])
#         #     cbar = plt.colorbar(im_j_list[-1], cax=cax1)
#         #     cbar.ax.set_title(r"$J$")
#         #
#         #     cax = fig.add_axes([1.1, 0.1, 0.03, 0.8])
#         #     cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax)
#         #     cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#         #
#         # if rows == 2 and cols == 2:
#         #
#         #     im_j_list = []
#         #     im_dunb_list = []
#         #     i_cl = 0
#         #     for i1 in range(len(ax_list[:][0])):
#         #         for i2 in range(len(ax_list[0][:])):
#         #             ax = ax_list[i1][i2]
#         #             cl = self.cl_cyl[i_cl]
#         #
#         #             im1 = self.plot_J(ax, cl)
#         #             im2 = self.plot_Dunb(ax, cl)
#         #
#         #             self.plot_vel_polar(ax, cl)
#         #
#         #             ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#         #                        [1e-5, 1e-4], colors="gray")
#         #             ax.set_rlim(self.rmin, self.rmax)
#         #             im_j_list.append(im1)
#         #             im_dunb_list.append(im2)
#         #             i_cl += 1
#         #
#         #     cax1 = fig.add_axes([0.9, 0.5, 0.03, 0.35])  # [*left*, *bottom*, *width*, *height*]
#         #     cbar = plt.colorbar(im_j_list[-1], cax=cax1, extend='both')
#         #     cbar.ax.set_title(r"$J$")
#         #
#         #     cax2 = fig.add_axes([0.9, 0.1, 0.03, 0.35])  # [*left*, *bottom*, *width*, *height*]
#         #     cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax2, extend='both')
#         #     cbar2.ax.set_title(r"$\rho_{unb;Bern}$")
#
#         if rows == 3 and cols == 2:
#
#             im_j_list = []
#             im_dunb_list = []
#             im_jf_list = []
#             i_cl = 0
#
#             for i1 in range(len(ax_list[:])):
#                 for i2 in range(len(ax_list[0][:])):
#                     ax = ax_list[i1][i2]
#                     cl = self.cl_cyl[i1] # here i2 loops over all profiles only.
#
#                     if i1 == 0 and i2 == 0: ax.title.set_text(r'Angular Momentum $\cdot r^{1/2}$')
#                     if i1 == 0 and i2 == 1: ax.title.set_text(r'Angular Momentum Flux')
#                     ax.annotate('t:{} ms'.format("%.1f"%(cl.time * 1e3)),
#                                 xy=(-150 / 180 * np.pi, 50),  # theta, radius
#                                 xytext=(0.03, 0.03),  # fraction, fraction
#                                 xycoords='axes fraction',
#                                 horizontalalignment='center',
#                                 verticalalignment='center'
#                                 )
#
#                     # ---
#                     if i2 == 0:
#                         im_j_list.append(self.plot_J(ax, cl))
#                     else:
#                         im_jf_list.append(self.plot_Jflux(ax, cl))
#
#                     # unbound density
#                     im_dunb_list.append(self.plot_Dunb(ax, cl))
#
#                     # self.plot_vel_polar(ax, cl)
#
#                     # density countours
#                     ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#                                [1e-5, 1e-4], colors="gray")
#
#                     self.plot_dens_decomp_phase(ax, cl, norm_to_m0=True, d_or_rho='d')
#
#                     ax.set_rlim(self.rmin, self.rmax)
#
#                     # ---
#                     i_cl += 1
#
#             cax1 = fig.add_axes([0.9, 0.63, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar = plt.colorbar(im_j_list[-1], cax=cax1, extend='both', format='%.1e')
#             cbar.ax.set_title(r"$Ang.Mom$ $\cdot r^{1/2}$")
#
#             cax2 = fig.add_axes([0.9, 0.37, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_jf_list[-1], cax=cax2, extend='both', format='%.1e')
#             cbar2.ax.set_title(r"$Ang.Mom.Flux$")
#
#             cax3 = fig.add_axes([0.9, 0.1, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax3, extend='both')
#             cbar2.ax.set_title(r"$D_{unb;Bern}$")
#         if rows == 3 and cols == 3:
#
#             im_j_list = []
#             im_dunb_list = []
#             im_jf_list = []
#             i_cl = 0
#
#             for i1 in range(len(ax_list[:])):
#                 for i2 in range(len(ax_list[0][:])):
#                     ax = ax_list[i1][i2]
#                     cl = self.cl_cyl[i1] # here i2 loops over all profiles only.
#
#                     if i1 == 0 and i2 == 0: ax.title.set_text(r'Angular Momentum $\cdot r^{1/2}$')
#                     if i1 == 0 and i2 == 1: ax.title.set_text(r'Angular Momentum Flux')
#                     if i1 == 0 and i2 == 2: ax.title.set_text(r'Unbound Matter Density')
#
#                     ax.annotate('t:{} ms'.format("%.1f"%(cl.time * 1e3)),
#                                 xy=(-150 / 180 * np.pi, 50),  # theta, radius
#                                 xytext=(0.03, 0.03),  # fraction, fraction
#                                 xycoords='axes fraction',
#                                 horizontalalignment='center',
#                                 verticalalignment='center'
#                                 )
#
#                     # ---
#                     if i2 == 0:
#                         im_j_list.append(self.plot_J(ax, cl))
#                     elif i2 == 1:
#                         im_jf_list.append(self.plot_Jflux(ax, cl))
#                     elif i2 == 2:
#                         im_dunb_list.append(self.plot_Dunb(ax, cl))
#                     else:
#                         raise ValueError("More cols of plots than tasks")
#
#                     # self.plot_vel_polar(ax, cl)
#
#                     # density countours
#                     ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#                                [1e-5, 1e-4], colors="gray")
#
#                     self.plot_dens_decomp_phase(ax, cl, norm_to_m0=True, d_or_rho='d')
#
#                     ax.set_rlim(self.rmin, self.rmax)
#
#                     # ---
#                     i_cl += 1
#
#             cax1 = fig.add_axes([0.9, 0.63, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar = plt.colorbar(im_j_list[-1], cax=cax1, extend='both', format='%.1e')
#             cbar.ax.set_title(r"$Ang.Mom$ $\cdot r^{1/2}$")
#
#             cax2 = fig.add_axes([0.9, 0.37, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_jf_list[-1], cax=cax2, extend='both', format='%.1e')
#             cbar2.ax.set_title(r"$Ang.Mom.Flux$")
#
#             cax3 = fig.add_axes([0.9, 0.1, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax3, extend='both')
#             cbar2.ax.set_title(r"$D_{unb;Bern}$")
#         if rows == 3 and cols == 4:
#
#             im_j_list = []
#             im_dunb_list = []
#             im_jf_list = []
#             i_cl = 0
#
#             for i1 in range(len(ax_list[:])):
#                 for i2 in range(len(ax_list[0][:])):
#                     ax = ax_list[i1][i2]
#                     cl = self.cl_cyl[i1] # here i2 loops over all profiles only.
#
#                     if i1 == 0 and i2 == 0: ax.title.set_text(r'Angular Momentum $\cdot r^{1/2}$')
#                     if i1 == 0 and i2 == 1: ax.title.set_text(r'Angular Momentum Flux')
#                     if i1 == 0 and i2 == 2: ax.title.set_text(r'Unbound Matter Density')
#                     if i1 == 0 and i2 == 3: ax.title.set_text(r'Temperature')
#
#                     ax.annotate('t:{} ms'.format("%.1f"%(cl.time * 1e2)),
#                                 xy=(-150 / 180 * np.pi, 50),  # theta, radius
#                                 xytext=(0.03, 0.03),  # fraction, fraction
#                                 xycoords='axes fraction',
#                                 horizontalalignment='center',
#                                 verticalalignment='center'
#                                 )
#
#                     # ---
#                     if i2 == 0:
#                         im_j_list.append(self.plot_J(ax, cl))
#                     elif i2 == 1:
#                         im_jf_list.append(self.plot_Jflux(ax, cl))
#                     elif i2 == 2:
#                         im_dunb_list.append(self.plot_Dunb(ax, cl))
#                     elif i2 == 3:
#                         self.plot_temp(ax, cl)
#                     else:
#                         raise ValueError("More cols of plots than tasks")
#
#                     # self.plot_vel_polar(ax, cl)
#
#                     # density countours
#                     ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#                                [1e-5, 1e-4], colors="gray")
#
#                     self.plot_dens_decomp_phase(ax, cl, norm_to_m0=True, d_or_rho='d')
#
#                     ax.set_rlim(self.rmin, self.rmax)
#
#                     # ---
#                     i_cl += 1
#
#             cax1 = fig.add_axes([0.9, 0.63, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar = plt.colorbar(im_j_list[-1], cax=cax1, extend='both', format='%.1e')
#             cbar.ax.set_title(r"$Ang.Mom$ $\cdot r^{1/2}$")
#
#             cax2 = fig.add_axes([0.9, 0.37, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_jf_list[-1], cax=cax2, extend='both', format='%.1e')
#             cbar2.ax.set_title(r"$Ang.Mom.Flux$")
#
#             cax3 = fig.add_axes([0.9, 0.1, 0.03, 0.23])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax3, extend='both')
#             cbar2.ax.set_title(r"$D_{unb;Bern}$")
#         if rows == 4 and cols == 4:
#
#             im_j_list = []
#             im_dunb_list = []
#             im_jf_list = []
#             im_temp_list = []
#             i_cl = 0
#
#             for i1 in range(len(ax_list[:])):
#                 for i2 in range(len(ax_list[0][:])):
#                     ax = ax_list[i1][i2]
#                     cl = self.cl_cyl[i1] # here i2 loops over all profiles only.
#
#                     if i1 == 0 and i2 == 0: ax.title.set_text(r'Angular Momentum $\cdot r^{1/2}$')
#                     if i1 == 0 and i2 == 1: ax.title.set_text(r'Angular Momentum Flux')
#                     if i1 == 0 and i2 == 2: ax.title.set_text(r'Unbound Matter Density')
#                     if i1 == 0 and i2 == 3: ax.title.set_text(r'Ye')
#
#                     ax.set_rlim(self.rmin, self.rmax)
#                     ax.annotate('t:{} ms'.format("%.1f"%(cl.time * 1e3)),
#                                 xy=(-150 / 180 * np.pi, 50),  # theta, radius
#                                 xytext=(0.03, 0.03),  # fraction, fraction
#                                 xycoords='axes fraction',
#                                 horizontalalignment='center',
#                                 verticalalignment='center'
#                                 )
#                     # ---
#                     if i2 == 0:
#                         im_j_list.append(self.plot_J(ax, cl))
#                         self.plot_dens_decomp_phase2(ax, cl, m_list=[1], norm_to_m0=True, int_over_z=True, d_or_rho='d')
#                     elif i2 == 1:
#                         im_jf_list.append(self.plot_Jflux(ax, cl))
#                         self.plot_dens_decomp_phase2(ax, cl, m_list=[1], norm_to_m0=True, int_over_z=True, d_or_rho='d')
#                     elif i2 == 2:
#                         im_dunb_list.append(self.plot_Dunb(ax, cl))
#                         ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#                                    self.dens_levels, colors="gray")
#                     elif i2 == 3:
#                         im_temp_list.append(self.plot_ye(ax, cl))
#                         ax.contour(cl.phi_3d[0, :, 0], cl.r_3d[:, 0, 0], cl.get_d_z_integ(),
#                                    self.dens_levels, colors="gray")
#                     else:
#                         raise ValueError("More cols of plots than tasks")
#
#                     # self.plot_vel_polar(ax, cl)
#
#                     # density countours
#
#
#                     # ---
#                     i_cl += 1
#
#             cax1 = fig.add_axes([0.9, 0.70, 0.02, 0.17])  # [*left*, *bottom*, *width*, *height*]
#             cbar = plt.colorbar(im_j_list[-1], cax=cax1, extend='both', format='%.1e')
#             cbar.ax.set_title(r"$Ang.Mom$ $\cdot r^{1/2}$")
#
#             cax2 = fig.add_axes([0.9, 0.50, 0.02, 0.17])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_jf_list[-1], cax=cax2, extend='both', format='%.1e')
#             cbar2.ax.set_title(r"$Ang.Mom.Flux$")
#
#             cax3 = fig.add_axes([0.9, 0.30, 0.02, 0.17])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_dunb_list[-1], cax=cax3, extend='both')
#             cbar2.ax.set_title(r"$D_{unb;Bern}$")
#
#             cax4 = fig.add_axes([0.9, 0.1, 0.02, 0.17])  # [*left*, *bottom*, *width*, *height*]
#             cbar2 = plt.colorbar(im_temp_list[-1], cax=cax4, extend='both')
#             cbar2.ax.set_title(r"Ye")
#
#
#         plt.subplots_adjust(hspace=0.2)
#         plt.subplots_adjust(wspace=-0.3)
#         # plt.tight_layout()
#         plt.savefig('{}{}'.format(out_name, self.set_save_format), bbox_inches='tight', dpi=128)
#         plt.close()




class FULL_DATA_ANALYSIS:

    def __init__(self, sim, files):

        self.files = files
        self.sim = sim
        # self.cyl_cls = []
        # for file in files:
        #     self.cyl_cls.append(READ_CYLH5(file))

    def save_density_modes(self, m_list=[], norm_to_m0=True, d_or_rho='d', int_over_z=True):
        '''
        Uses 2D slize
        :param m_list:
        :param norm_to_m0:
        :param d_or_rho:
        :return:
        '''

        Matrix1 = [[0 for x in range(len(self.files))] for y in range(len(m_list))]
        Matrix2 = [[0 for x in range(len(self.files))] for y in range(len(m_list))]
        Matrix3 = [[0 for x in range(len(self.files))] for y in range(len(m_list))]
        Matrix4 = [[0 for x in range(len(self.files))] for y in range(len(m_list))]

        for i_file, file in enumerate(self.files):
            cyl_cl = READ_CYLH5(file)

            # tst:

            # rho_edges = 10.0 ** np.linspace(4.0, 16.0, 120)
            # temp_edges = 10.0 ** np.linspace(-2, 2, 50)
            # Ye_edges = np.linspace(0, 0.5, 50)
            # dens = cyl_cl.D
            # delta = cyl_cl.dz_3d * cyl_cl.dr_3d * cyl_cl.dphi_3d

            # rho = cyl_cl.rho
            # temp = cyl_cl.temp
            # Ye = cyl_cl.ye
            # data = np.column_stack((rho.flatten(), temp.flatten(), Ye.flatten()))
            # H, _ = np.histogramdd(data,
            #                       bins=(rho_edges, temp_edges, Ye_edges),
            #                       weights=dens.flatten() * np.prod(delta) * 2)
            #
            # plt.plot(rho_edges[:-1], np.sum(H, axis=(1,2)), '.', color="black")
            # plt.savefig('{}{}.png'.format(LISTS.fig_dir, 'test_my_hist'), bbox_inches='tight', dpi=128)
            # plt.close()

            # print(H.shape)
            # exit(1)

            for i_m, m in enumerate(m_list):
                if d_or_rho == 'd': using = cyl_cl.D
                else: using = cyl_cl.rho
                if int_over_z:
                    m_int_phi, m_int_phi_r = \
                        PHYSICS.get_dens_decomp_3d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, cyl_cl.dz_3d, m)
                    if norm_to_m0:
                        _int_phi, _int_phi_r = \
                            PHYSICS.get_dens_decomp_3d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, cyl_cl.dz_3d, 0)
                        m_int_phi /= _int_phi_r
                        m_int_phi_r /= _int_phi_r
                else:
                    m_int_phi, m_int_phi_r = \
                        PHYSICS.get_dens_decomp_2d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, m)
                    if norm_to_m0:
                        _int_phi, _int_phi_r = \
                            PHYSICS.get_dens_decomp_2d(using, cyl_cl.phi_3d, cyl_cl.dphi_3d, cyl_cl.dr_3d, 0)
                        m_int_phi /= _int_phi_r
                        m_int_phi_r /= _int_phi_r

                Matrix1[i_m][i_file] = m_int_phi
                Matrix2[i_m][i_file] = m_int_phi_r
                Matrix3[i_m][i_file] = cyl_cl.time
                Matrix4[i_m][i_file] = cyl_cl.r_3d[:, 0, 0]

            cyl_cl.__delete__(cyl_cl) # clear the arrays and instances

        dfile = h5py.File(Paths.ppr_sims + self.sim + '/densitymodes_from_3d.h5', "w")

        start_t = time.time()
        print("\t Saving '{}' into sim. dir...".format('densitymodes_from_3d.h5')),
        for i_m, m in enumerate(m_list):

            r =  np.array(Matrix4[i_m][0]) # assuming r is the same

            int_phi_reshaped = np.reshape(np.array(Matrix1[i_m]), (len(self.files), len(r)))

            group = dfile.create_group("rl=%d" % m)
            group["int_phi"] = int_phi_reshaped
            group["int_phi_r"] = Matrix2[i_m] # | time   <-> r
            group["time"] = Matrix3[i_m]
            group["r"] = r

        dfile.close()
        print(" done! (%.2f sec)" % (time.time() - start_t))

    def plot_density_modes(self, sim, m_list=[], load_file='densitymodes_from_3d.h5'):

        dfile = h5py.File(LISTS.loc_of_sims + sim + '/'+ load_file, "r")

        rows = 2
        cols = 1

        fig = plt.figure(figsize=(6.5, 1.5 * 3.6))  # figsize=(4.5, 2.5 * 3.6)  # (<->; v)
        ax_list = []
        for n in range(1, rows + 1):
            if n == 1:
                ax_list.append(fig.add_subplot(rows, cols, n))
            else:
                ax_list.append(fig.add_subplot(rows, cols, n, sharex=ax_list[n - 2]))  # sharex=axs[n - 2]))

        for m in m_list:

            if m == 1: color = 'black'; ls = '-'; lw=1.
            elif m == 2: color = 'gray';  ls = '-.'; lw=1.
            elif m == 3: color = 'blue';  ls = '-.'; lw=0.4
            elif m == 4: color = 'orange';  ls = '-.'; lw=0.4
            elif m == 5: color = 'green';  ls = '-.'; lw=0.4
            elif m == 6: color = 'pink';  ls = '-.'; lw=0.4
            elif m == 7: color = 'purple';  ls = '-.'; lw=0.4
            else: raise ValueError('m is not in color/ls list')

            group = dfile["rl=%d" % m]
            int_phi2d = np.array(group["int_phi"])
            int_phi_r1d = np.array(group["int_phi_r"])  # | time   <-> r
            times = np.array(group["time"])
            r = np.array(group["r"])

            # phase plot
            ax_list[0].plot(times * 1e3, np.unwrap(np.angle(int_phi_r1d)), ls, lw=lw, color=color, label='m:{}'.format(m))

            ax_list[0].set_ylabel(r'$C_m/C_0$ Phase')
            # ax_list[0].annotate(r'$C_m = \int{\rho W \sqrt{-g}\cdot\exp(i m \phi) dz dr d\phi}$',
            #             xy=(-150 / 180 * np.pi, 50),  # theta, radius
            #             xytext=(0.65, 0.90),  # fraction, fraction
            #             xycoords='axes fraction',
            #             horizontalalignment='center',
            #             verticalalignment='center'
            #             )
            ax_list[0].legend()

            # magnitude plot
            ax_list[1].plot(times * 1e3, np.abs(int_phi_r1d), ls, lw=lw,  color=color, label='m:{}'.format(m))

            ax_list[1].set_yscale('log')
            ax_list[1].set_ylabel(r'$C_m/C_0$ Magnitude')
            ax_list[1].legend()

        plt.savefig('{}{}.png'.format(LISTS.fig_dir, 'desnity_modes'), bbox_inches='tight', dpi=128)
        plt.close()


# --- --- ---

class READ_CYLH5:

    def __init__(self, fname):

        self.fname = fname
        self.file = h5py.File(fname, "r")

        self.iteration = self.fname.split('/')[-1].split('.')[0]

        file = h5py.File(fname, 'r')
        print("Data File: {} [".format(self.iteration)),
        for i in file:
            print(i),
        print(']'),

        self.phi_3d     = np.array(self.file["phi_cyl_3d"])
        self.r_3d       = np.array(self.file["r_cyl_3d"])
        self.z_3d       = np.array(self.file["z_cyl_3d"])
        self.dphi_3d    = np.array(self.file["dphi_cyl_3d"])
        self.dr_3d      = np.array(self.file["dr_cyl_3d"])
        self.dz_3d      = np.array(self.file["dz_cyl_3d"])
        self.J          = np.array(self.file["J"])
        self.Jflux      = np.array(self.file["Jflux"])
        self.D          = np.array(self.file["D"])
        self.vphi       = np.array(self.file["vphi"])
        self.vr         = np.array(self.file["vr"])
        self.Dgeo       = np.array(self.file["Dgeo"])
        self.Dbern      = np.array(self.file["Dbern"])
        self.rho        = np.array(self.file["rho"])
        self.temp       = np.array(self.file["temp"])
        self.ye         = np.array(self.file["Ye"])

        self.time = self.file.attrs["time"] * 0.004925794970773136 * 1e-3
        # units.conv_time(units.cactus, units.cgs, self.file.attrs["time"])
        print(" | loaded it:{} time: {} ms".format(int(self.iteration), "%.1f" % (self.time * 1e3)))

    def __delete__(self, instance):
        del instance.phi_3d
        del instance.r_3d
        del instance.z_3d
        del instance.dphi_3d
        del instance.dr_3d
        del instance.dz_3d
        del instance.J
        del instance.Jflux
        del instance.D
        del instance.vphi
        del instance.vr
        del instance.Dgeo
        del instance.Dbern
        del instance.rho
        del instance.temp

        instance.file.close()
        gc.collect()

class ANALYZE_CYLINDRICAL:

    def __init__(self, fname):
        '''
        LOADS AND ANALYZES THE .cylh5 file,
        that contatins the CYLINDRICAL grid and interpolated variables.
        :param fname:
        '''

        self.fname = fname
        self.file = h5py.File(fname, "r")

        self.iteration = self.fname.split('/')[-1].split('.')[0]

        self.file = h5py.File(fname, 'r')
        print("Data File: {} [".format(self.iteration)),
        for i in self.file:
            print(i),
        print(']'),

        # self.time = self.file.attrs["time"] * 0.004925794970773136 * 1e-3
        # units.conv_time(units.cactus, units.cgs, self.file.attrs["time"])
        # print(" | loaded it:{} time: {} ms".format(int(self.iteration), "%.1f" % (self.time * 1e3)))

    def __delete__(self, instance):

        instance.file.close()
        gc.collect()


    def get_arr(self, v_n):

        if v_n not in self.file.keys():
            raise NameError("v_n: {} not in the available v_ns in the file")

        return np.array(self.file[v_n])

    def get_z_integ(self, v_n, multiplier=2.):

        return multiplier * np.sum(self.get_arr(v_n) * self.get_arr("dz_cyl"), axis=(2))

    def get_z_integ__phi_ave(self, v_n):
        # J integrate over 'z'; average over 'phi' for every 'r'

        # if v_n not in self.file.keys():
        #     raise NameError("v_n: {} not in the available v_ns in the file")
        # J_z_summed = 2 * np.sum(self.get_arr(v_n) * self.get_arr("dz_cyl"), axis=(2))   # for every r, phi, summ over z
        # J_z_phi_summed = np.sum(J_z_summed, axis=(1))   # for every r, summ over all phi

        J_z_summed = self.get_z_integ(v_n)

        J_z_summed_ave = np.zeros(len(J_z_summed[0, :]))
        for i in range (len(J_z_summed[:, 0])):
            J_min_J_ave = J_z_summed[i, :] - (np.sum(J_z_summed[i, :]) / len(J_z_summed[i, :]))
            # J_min_J_ave *= self.get_arr("r_cyl")[i, 0, 0] ** (.5)
            J_z_summed_ave = np.vstack((J_z_summed_ave, J_min_J_ave))
        J_z_summed_ave = np.delete(J_z_summed_ave, 0, 0)

        return J_z_summed_ave

    # def get_jf_z_integ(self):
    #
    #     # J: integrate over 'z'; average over 'phi' for every 'r'
    #     Jf_z_summed = 2 * np.sum(self.Jflux * self.dz_3d, axis=(2))  # for every r, phi, summ over z
    #     # J_z_phi_summed = np.sum(J_z_summed, axis=(1))   # for every r, summ over all phi
    #
    #     return Jf_z_summed

    # def get_d_z_integ(self):
    #     return 2 * np.sum(self.D * self.dz_3d, axis=(2))

    # def get_dbern_z_integ(self):
    #     return 2 * np.sum(self.Dbern * self.dz_3d, axis=(2))

    def get_slice_z0(self, v_n):

        return self.get_arr(v_n)[:, :, 0]

    def get_slice_z0__phi_ave(self, v_n):

        slice = self.get_slice_z0(v_n)

        vphi_minus_phi_ave = np.zeros(len(self.vphi[0, :, 0]))
        for i in range(len(slice[:, 0, 0])):
            vphi_minus_phi_ave_i = slice[i, :] - (np.sum(slice[i, :]) / len(slice[i, :]))
            vphi_minus_phi_ave = np.vstack((vphi_minus_phi_ave, vphi_minus_phi_ave_i))

        vphi_minus_phi_ave = np.delete(vphi_minus_phi_ave, 0, 0)
        return vphi_minus_phi_ave

    def get_arr_mod(self, v_n, mod='slice'):

        # available modes: z_int_phi_ave, slice,

        if mod == 'z_int_phi_ave':
            return self.get_z_integ__phi_ave(v_n)

        elif mod == 'slice':
            return self.get_slice_z0(v_n)

        elif mod == 'slice_phi_ave':
            return self.get_slice_z0__phi_ave(v_n)

        elif mod == 'z_int':
            return self.get_z_integ(v_n)

        elif mod == '-':
            return self.get_arr(v_n)

        else:
            raise NameError('mod: {} for v_n: {} not recognized'
                            .format(mod, v_n))

    # def get_vphi_z0_ave(self):
    #
    #     vphi_z_ave = self.vphi[:, :, 0]  # np.sum(vphi_cyl, axis=2) / len(vphi_cyl[0, 0, :])
    #     vphi_minus_phi_ave = np.zeros(len(self.vphi[0, :, 0]))
    #
    #     for i in range(len(self.vphi[:, 0, 0])):
    #         vphi_minus_phi_ave_i = vphi_z_ave[i, :] - (np.sum(vphi_z_ave[i, :]) / len(vphi_z_ave[i, :]))
    #         vphi_minus_phi_ave = np.vstack((vphi_minus_phi_ave, vphi_minus_phi_ave_i))
    #
    #     vphi_minus_phi_ave = np.delete(vphi_minus_phi_ave, 0, 0)
    #     return vphi_minus_phi_ave
    #
    # def get_vr_z0_ave(self):
    #
    #     # velocities: Get z=0 slice ; v - sum_phi(v)/len(v) for every 'r'
    #     vr_z_ave = self.vr[:, :, 0]  # np.sum(vr_cyl, axis=2) / len(vr_cyl[0, 0, :])
    #     vr_minus_phi_ave = np.zeros(len(self.vr[0, :, 0]))
    #
    #     for i in range(len(self.vphi[:, 0, 0])):
    #         vr_minus_phi_ave_i = vr_z_ave[i, :] - (np.sum(vr_z_ave[i, :]) / len(vr_z_ave[i, :]))
    #         vr_minus_phi_ave = np.vstack((vr_minus_phi_ave, vr_minus_phi_ave_i))
    #
    #     vr_minus_phi_ave = np.delete(vr_minus_phi_ave, 0, 0)
    #     return vr_minus_phi_ave

    # def get_temperature(self):
    #
    #     return self.temp[:, :, 0]
    #
    # def get_ye(self):
    #
    #     return self.ye[:, :, 0]

    def get_dens_decomp(self, dens_2d, phi_2d, dphi_2d, dr_2d, m=1):
        '''
        Returns complex arrays [\int(d\phi)] and [\int(dr d\phi)]
        '''

        integ_over_phi = np.sum(dens_2d * np.exp(1j * m * phi_2d) * dphi_2d, axis=1)

        integ_over_phi_r = np.sum(integ_over_phi * dr_2d[:, 0])

        return integ_over_phi, integ_over_phi_r

class PLOT_CYLINDRICAL:

    def __init__(self):
        self.gen_set_dic = {
            "style": "polar",
            "figir": "./",
            "figformat": 'tst.png',
            "figsize": (8.5, 6.5),
            "file_extension": ".cylh5", # file type is "12345.cylh5"
            "path": Paths.ppr_sims + "DD2_M13641364_M0_LK_LR_R04/3d/"
        }

        self.plot_dics = []

        ang_mom = {
            'it':        [254606],
            'v_n':       ['ang_mom'],
            'vmin':      [-5e-4],  # for plotting
            'vmax':      [5e-4],
            'rmin':      [0.],
            'rmax':      [50.],
            'mask_below':[None],
            'mask_above':[None],
            'cmap':      ['RdBu_r'],
            'norm':      ['log'],
            'todo':      ['z_int_phi_ave']
        }
        self.plot_dics.append(ang_mom)

        ang_mom_flux = {
            'it':        [254606],
            'v_n':       ['ang_mom_flux'],
            'vmin':      [-5e-6],  # for plotting
            'vmax':      [5e-6],
            'rmin':      [0.],
            'rmax':      [50.],
            'mask_below':[None],
            'mask_above':[None],
            'cmap':      ['RdBu_r'],
            'norm':      ['log'],
            'todo':      ['z_int_phi_ave']
        }
        self.plot_dics.append(ang_mom_flux)

        dens_unb_bern = {
            'it':        [254606],
            'v_n':       ['dens_unb_bern'],
            'vmin':      [1e-7],  # for plotting
            'vmax':      [6e-7],
            'rmin':      [0.],
            'rmax':      [50.],
            'mask_below':[None],
            'mask_above':[None],
            'cmap':      ['Greens'],
            'norm':      ['log'],
            'todo':      ['z_int_phi_ave']
        }
        self.plot_dics.append(dens_unb_bern)

        temp = {
            'it':        [254606],
            'v_n':       ['temp'],
            'vmin':      [5.],  # for plotting
            'vmax':      [10.],
            'rmin':      [0.],
            'rmax':      [50.],
            'mask_below':[None],
            'mask_above':[None],
            'cmap':      ['Oranges'],
            'norm':      ['log'],
            'todo':      ['slice']
        }
        self.plot_dics.append(temp)

        ye = {
            'it':        [254606],
            'v_n':       ['Ye'],
            'vmin':      [0.1],  # for plotting
            'vmax':      [0.4],
            'rmin':      [0.],
            'rmax':      [50.],
            'mask_below':[None],
            'mask_above':[None],
            'cmap':      ['Oranges'],
            'norm':      ['norm'],
            'todo':      ['slice']
        }
        self.plot_dics.append(ye)

    def plot_one_task(self, ax, data_cl, task_dic, position=0.):

        ax.set_rlim(task_dic["rmin"], task_dic["rmax"])

        data = data_cl.get_arr_mod(task_dic["v_n"], task_dic["todo"])
        phi = data_cl.get_slice_z0("phi_cyl")
        r = data_cl.get_slice_z0("r_cyl")

        # set maximum/minimum
        if task_dic["vmin"] == None: task_dic["vmin"] = data.min()
        if task_dic["vmax"] == None: task_dic["vmax"] = data.max()

        # set norm
        if task_dic["norm"] == "norm":
            norm = LogNorm(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        elif task_dic["norm"] == "log":
            norm = Normalize(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        else:
            raise NameError("unrecognized norm: {} in task {}".format(task_dic["norm"],
                                                                      task_dic["v_n"]))

        # set masks
        if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
            im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
            im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
                                norm=norm, cmap=task_dic["cmap"])
        elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
            im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
                                norm=norm, cmap=task_dic["cmap"])
        else:
            im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
                                np.ma.masked_array(data,
                                                   (data > task_dic["mask_above"])
                                                   & (data < task_dic["mask_below"])),
                                norm=norm, cmap=task_dic["cmap"])


        # # if the norm
        # if task_dic["norm"] == 'norm':
        #     if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
        #         im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        #     elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
        #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
        #                             norm=norm, cmap=task_dic["cmap"])
        #     elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
        #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
        #                             norm=norm, cmap=task_dic["cmap"])
        #     else:
        #         im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
        #                             np.ma.masked_array(data,
        #                                                (data > task_dic["mask_above"])
        #                                                & (data < task_dic["mask_below"])),
        #                             norm=norm, cmap=task_dic["cmap"])
        #
        # elif task_dic["norm"] == 'log':
        #     if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
        #         im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        #     elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
        #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
        #                             norm=norm, cmap=task_dic["cmap"])
        #     elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
        #         im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
        #                             norm=norm, cmap=task_dic["cmap"])
        #     else:
        #         im2 = ax.pcolormesh(phi, r,
        #                             np.ma.masked_array(data,
        #                                                (data > task_dic["mask_above"])
        #                                                & (data < task_dic["mask_below"])),
        #                             norm=norm, cmap=task_dic["cmap"])
        # else:
        #     raise NameError("norm:{} is not recognized (task:{})".format(task_dic["norm"],
        #                                                                  task_dic["v_n"]))

        im2.set_rasterized(True)
        return im2


        # if task_dic["v_n"] == "ang_mom":
        #     im = ax.pcolormesh(data_cl.get_slice_z0("phi_cyl"),
        #                        data_cl.get_slice_z0("r_cyl"),
        #                        data_cl.get_z_integ__phi_ave(task_dic["v_n"]), # int and ave
        #                        cmap=task_dic["cmap"],
        #                        vmin=task_dic["vmin"],
        #                        vmax=task_dic["vmax"])
        #     im.set_rasterized(True)
        #     return im
        #
        # elif task_dic["v_n"] == "ang_mom_flux":
        #     im = ax.pcolormesh(data_cl.get_slice_z0("phi_cyl"),
        #                        data_cl.get_slice_z0("r_cyl"),
        #                        data_cl.get_z_integ(task_dic["v_n"]), # int only
        #                        cmap=task_dic["cmap"],
        #                        vmin=task_dic["vmin"],
        #                        vmax=task_dic["vmax"])
        #     im.set_rasterized(True)
        #     return im
        #
        # elif task_dic["v_n"] in ["dens_unb_bern", "temp", "ye"]:
        #     data = data_cl.get_z_integ(task_dic["v_n"])
        #     phi = data_cl.get_slice_z0("phi_cyl")
        #     r = data_cl.get_slice_z0("r_cyl")
        #
        #     if task_dic["vmin"] == None: task_dic["vmin"] = data.min()
        #     if task_dic["vmax"] == None: task_dic["vmax"] = data.max()
        #
        #     if task_dic["norm"] == "norm":
        #         norm = LogNorm(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        #     elif task_dic["norm"] == "log":
        #         norm = Normalize(vmin=task_dic["vmin"], vmax=task_dic["vmax"])
        #     else:
        #         raise NameError("unrecognized norm: {} in task {}".format(task_dic["norm"],
        #                                                                   task_dic["v_n"]))
        #     # if the norm
        #     if task_dic["norm"] == 'norm':
        #         if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
        #             im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        #         elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
        #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #         elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
        #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #         else:
        #             im2 = ax.pcolormesh(data_cl.get_arr("phi_cyl"), data_cl.get_arr("r_cyl"),
        #                                 np.ma.masked_array(data,
        #                                                    (data > task_dic["mask_above"])
        #                                                    & (data < task_dic["mask_below"])),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #     elif task_dic["norm"] == 'log':
        #         if task_dic["mask_above"] == None and task_dic["mask_below"] == None:
        #             im2 = ax.pcolormesh(phi, r, data, norm=norm, cmap=task_dic["cmap"])
        #         elif task_dic["mask_below"] != None and task_dic["mask_above"] == None:
        #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data < task_dic["mask_below"]),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #         elif task_dic["mask_below"] == None and task_dic["mask_above"] != None:
        #             im2 = ax.pcolormesh(phi, r, np.ma.masked_array(data, data > task_dic["mask_above"]),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #         else:
        #             im2 = ax.pcolormesh(phi, r,
        #                                 np.ma.masked_array(data, (data > task_dic["mask_above"])
        #                                                    & (data < task_dic["mask_below"])),
        #                                 norm=norm, cmap=task_dic["cmap"])
        #     else:
        #         raise NameError("norm:{} is not recognized (task:{})".format(task_dic["norm"],
        #                                                                      task_dic["v_n"]))
        #
        #     im2.set_rasterized(True)
        #     return im2
        #
        # else:
        #     raise NameError("task: {} not recognized".format(task_dic["v_n"]))

    def plot_from_dics(self):
        """
        1 row per plot
        :return:
        """

        # check if all entries in the dictionary have the same length
        dum_n = len(self.plot_dics[0]['it'])
        for n_task in range(len(self.plot_dics)):
            if dum_n != len(self.plot_dics[n_task]['it']):
                raise NameError("plot_dic[{}] contains wrong n_of_elements "
                                "(should be the same for all)".format(self.plot_dics[n_task]['it']))

        n_rows = len(self.plot_dics) #  |
        n_cols = len(self.plot_dics[0]['it'])  # <->

        print("Plotting: {} rows {} columns".format(n_rows, n_cols))

        fig = plt.figure(figsize=self.gen_set_dic['figsize'])  # (<->; v)
        # initializing the matrix with dummy axis objects
        sbplot_matrix = [[fig.add_subplot(n_rows, n_cols, 1, projection='polar')
                          for x in range(n_rows)] for y in range(n_cols)]

        # loading data processing classes
        data_cl = []
        for it in self.plot_dics[0]['it']:
            # file type is "12345.cylh5"
            fname = self.gen_set_dic["path"] + str(int(it))+self.gen_set_dic["file_extension"]
            data_cl.append(ANALYZE_CYLINDRICAL(fname))

        i = 1
        for n_row in range(n_rows):
            for n_col in range(n_cols):

                if n_col == 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')

                elif n_col == 0 and n_row > 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')#,
                                                                  # sharex=sbplot_matrix[n_col][0])

                elif n_cols > 0 and n_row == 0:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')#,
                                                                  # sharey=sbplot_matrix[0][n_row])

                else:
                    sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i, projection='polar')#,
                                                                  # sharex=sbplot_matrix[n_col][0],
                                                                  # sharey=sbplot_matrix[0][n_row])

                    # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                i += 1

        # plotting for every point in the grid
        for n_row in range(n_rows):
            im = None
            for n_col in range(n_cols):
                ax = sbplot_matrix[n_col][n_row]

                task_dic = {}
                for v_n in self.plot_dics[n_row].keys():
                    task_dic[v_n] = self.plot_dics[n_row][v_n][n_col]
                im = self.plot_one_task(ax, data_cl[n_col], task_dic, n_cols-(n_col+1))


            # adding axis for colobar and plotting colobars for every row
            ax_for_row = sbplot_matrix[-1][n_row]
            pos1 = ax_for_row.get_position()  # get the original position
            pos2 = [pos1.x0 + pos1.width,
                    pos1.y0, #- 0.5 * pos1.height,
                    0.02, pos1.height]
            cax1 = fig.add_axes(pos2)
            cbar = plt.colorbar(im, cax=cax1, extend='both', format='%.1e')
            cbar.ax.set_title(r"{}".format(str(self.plot_dics[n_row]["v_n"][-1]).replace('_', '\_')))


        plt.subplots_adjust(hspace=0.2)
        plt.subplots_adjust(wspace=-0.3)
        # plt.tight_layout()
        plt.savefig('{}{}'.format(self.gen_set_dic["figir"], self.gen_set_dic["figformat"]),
                    bbox_inches='tight', dpi=128)
        plt.close()

# --- --- ---

# def get_profiles(sim, time_list = [], it_list = [], n_more=0, ftype = '.cylh5', time_units='s', add_path='profiles/'):
#     '''
#     Looks at all the available profiles in the sim.dir + '/profiles/'.
#     Using the .it_time_file, interpolates [times] for all [iterations]
#
#     '''
#
#     if len(time_list) > 0 and len(it_list) > 0: raise NameError("Only time_list or it_list can be specified")
#
#     files = sorted(glob(Paths.ppr_sims + sim + "/" + add_path + "*{}".format(ftype)),
#                    key=lambda x: int(x.split('/')[-1].replace("{}".format(ftype), "")))
#     iterations = [int(it.split('/')[-1].replace("{}".format(ftype), "")) for it in files]
#     get_times = interpoate_time_form_it(iterations, Paths.ppr_sims + sim + '/' + Files.it_time, time_units)
#
#     print('|------------------------------------------------------------------------|')
#     print("ALL TIMES: \n{}".format(get_times))
#     print("ALL ITERATIONS: \n{}".format(np.array(iterations, dtype=int)))
#     print('|------------------------------------------------------------------------|')
#
#     # if nothing specified, return all the files found and all the timestpes
#     if not any(time_list) and not any(it_list):
#         return get_times, files
#
#     # if times to find are given, but not iterations, - compute iterations for those times and find files
#     elif any(time_list) and not any(it_list):
#         get_iterations = interpoate_it_form_time(time_list, Paths.ppr_sims + sim + '/' + Files.it_time, time_units)
#
#     # if times are not given, but iterations to use are, - compute times for those iterations
#     elif not any(time_list) and any(it_list):
#         if np.array(it_list, dtype=int).min() < np.array(iterations, dtype=int).min():
#             raise ValueError("Given it_list it:{} < iterations.min()"
#                              .format(np.array(it_list, dtype=int).min(), np.array(iterations, dtype=int).min()))
#         if np.array(it_list, dtype=int).max() > np.array(iterations, dtype=int).max():
#             raise ValueError("Given it_list it:{} < iterations.min()"
#                              .format(np.array(it_list, dtype=int).max(), np.array(iterations, dtype=int).max()))
#         get_iterations = it_list
#     else:
#         raise IOError("Input is not recongized.")
#
#     # select those iterations that are present in the 'iterations' list of files.
#     available_iterations = []
#     files_of_available_iterations = []
#     if len(get_iterations) < 1:
#         raise ValueError("No iterations for times:{} were found".format(time_list))
#
#     for it in get_iterations:
#         available_iterations.append(iterations[find_nearest_index(np.array(iterations, dtype=int), it)])
#         fname = Paths.ppr_sims + sim + "/" + add_path + "{}{}".format(available_iterations[-1], ftype)
#         files_of_available_iterations.append(fname)
#
#     if len(available_iterations) < 1:
#         raise ValueError("No available iterations seleted from required list: {}".format(get_iterations))
#     if len(available_iterations) < len(get_iterations):
#         raise ValueError("N of available it:{} != N of required:{}".format(len(available_iterations), len(get_iterations)))
#
#     # if just one time/iteration is given and n_more > 0, get MORE iterations (consecutive)
#     if n_more > 0:
#         if len(time_list) == 1 or len(it_list) == 1:
#             idx = find_nearest_index(np.array(iterations, dtype=int), int(available_iterations[0]))
#             for i in range(1, n_more):
#                 fname = Paths.ppr_sims + sim + "/" + add_path + "{}{}".format(iterations[int(idx + i)], ftype)
#                 available_iterations.append(iterations[idx + i])
#                 files_of_available_iterations.append(fname)
#
#     available_times = interpoate_time_form_it(available_iterations,
#                                               Paths.ppr_sims + sim + '/' + Files.it_time, time_units)
#
#     return available_times, files_of_available_iterations

if __name__ == '__main__':

    ''' ANG. MOMENTUM '''
    times, files = get_profiles("DD2_M13641364_M0_LK_LR_R04",
                                time_list=[], it_list=[], n_more=0, ftype='.cylh5', time_units='s', add_path='3d/')
    # 0.050, 0.060, 0.070, 0.080
    # pl = PLOT_CYLINDRICAL(files)
    # pl.plot_j_and_jflux_for_3_profiles(LISTS.fig_dir + 'tst_j_jf')

    pl = PLOT_CYLINDRICAL()
    pl.plot_from_dics()

    exit(0)
    ''' FULL DATA ANALYSIS '''
    # times, files = get_profiles("DD2_M13641364_M0_SR",
    #                             time_list=[0.050], it_list=[], n_more=0, ftype='.cylh5', time_units='s')
    # full = FULL_DATA_ANALYSIS("DD2_M13641364_M0_SR", files)
    # full.save_density_modes(m_list=[1, 2, 3, 4, 5, 6, 7])
    # full.plot_density_modes("DD2_M13641364_M0_SR", m_list=[1, 2, 3, 4, 5, 6, 7], load_file="densitymodes_from_3d.h5")

    # tst_plot_dens_modes("DD2_M13641364_M0_SR", [1, 2], "")

    # LOAD_CYLINDRICAL('389962_cyl.h5');
    # exit(1)


    # am = ANALYZE_CYLINDRICAL('389962.cylh5')
    # am.test_plot_angular_momentum(0, 50)
    # am.plot_dens_decomp()
    #
    # pl = PLOT_CYLINDRICAL(['860718.cylh5', '880652.cylh5', '897036.cylh5'])#, '411290.cylh5', '431090.cylh5', '452608.cylh5'])
    # pl.plot_2_files()
    # pl.test_plot_ang_mom_flux('tst_ang_mom_flux')
    # pl.plot_j_and_jflux_for_3_profiles('tst_j_jf')
    # exit(0)






    # exit(0)
    #
    # files = sorted(glob("*.h5"), key=lambda x: int(x.replace(".h5", "")))
    # print("Total files loaded: {}".format(len(files)))
    #
    # it_dmass = []
    # for file in files:
    #
    #     # disk mass estimation
    #     it = np.int(file.split('.h5')[0])
    #     start_t = time.time()
    #     # print("Computing disk mass for it:{} ... ".format(it)),
    #
    #     save = SAVE_CYLINDRICAL(file, ['J', 'D', 'vphi', 'vr', 'Dgeo', 'Dbern', 'Dgerch'])
    #     save.main()
    #     exit(1)


        # disk = COMPUTE_DISK(file)
        # disk.set_print_add_info = False

        # disk_mass = disk.get_disk_mass(True) # Save plots
        # it_dmass = np.append(it_dmass, [it, disk_mass])
        # print("done! ({:.2f} sec)".format(time.time() - start_t))

        # disk angular momentum
        # ang = DISK_ANGULAR_MOMENTUM(file)
        # ang.compute_for_all_reflevels(False) # remove ghost points - False. (True does not work)
        # saves into MJ_encl_ITERATION.txt


    # it_dmass_reshaped = np.reshape(it_dmass, (len(files), 2))
    # np.savetxt('./disk_mass.asc', it_dmass_reshaped, '%.6f', '  ', '\n',
    #            '\n rho from {} to {} , reflevel:{}  '
    #            .format(disk.rho_mask_bottom, disk.rho_mask_up, disk.NLEVELS), '', '# 1:iteration 3:data')





    # prof = READ_PROFILE(files[-1])
    # grid_ = prof.grid
    # print(prof.iteration)
    # print(grid_[0].delta)
    # # print(len(grid_.mesh()))
    # x, y, z = prof.get_x_y_z(0)
    # print(prof.get_arr('dens',0,True))