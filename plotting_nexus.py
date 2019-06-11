from __future__ import division
from sys import path
import matplotlib
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
from matplotlib.colors import LogNorm, Normalize
from matplotlib.ticker import AutoMinorLocator, FixedLocator, NullFormatter, \
    MultipleLocator
import numpy as np


class BASIC_PARTS():

    def __init__(self):
        pass

    def set_plot_title(self, ax, plot_dic):
        if "title" in plot_dic.keys():
            if plot_dic["title"] != '' and plot_dic["title"] != None:

                title = plot_dic["title"]

                # data = plot_dic['data']
                #
                # if plot_dic["title"] == 'it':
                #     title = plot_dic["it"]
                # elif plot_dic["title"] == 'time [s]' or \
                #     plot_dic["title"] == 'time':
                #     title = "%.3f" % data.get_time(plot_dic["it"]) + " [s]"
                # elif plot_dic["title"] == 'time [ms]':
                #     title = "%.1f" % (data.get_time(plot_dic["it"]) * 1000) + " [ms]"
                # else:
                #     title = plot_dic["title"]
                ax.title.set_text(r'{}'.format(title))

    def set_min_max_scale(self, ax, dic):
        if dic["ptype"] == "cartesian":
            if "xmin" in dic.keys() and "xmax" in dic.keys():
                if dic["xmin"] != None and dic["xmax"] != None:
                    ax.set_xlim(dic["xmin"], dic["xmax"])
            if "ymin" in dic.keys() and "ymax" in dic.keys():
                if dic["ymin"] != None and dic["ymax"] != None:
                    ax.set_ylim(dic["ymin"], dic["ymax"])
            if "xscale" in dic.keys():
                if dic["xscale"] == 'log':
                    ax.set_xscale("log")
            if "yscale" in dic.keys():
                if dic["yscale"] == 'log':
                    ax.set_yscale("log")
        elif dic["ptype"] == "polar":

            if "phimin" in dic.keys() and "phimax" in dic.keys():
                if dic["phimin"] != None and dic["phimax"] != None:
                    ax.set_philim(dic["phimin"], dic["phimax"])

            if "rmin" in dic.keys() and "rmax" in dic.keys():
                if dic["rmin"] != None and dic["rmax"] != None:
                    ax.set_rlim(dic["rmin"], dic["rmax"])

            if "xscale" in dic.keys():
                if dic["xscale"] == 'log':
                    raise NameError("log scale is not available for x in polar")

            if "yscale" in dic.keys():
                if dic["yscale"] == 'log':
                    raise NameError("log scale is not available for y in polar")
        else:
            raise NameError("Unknown 'ptype' of the plot: {}".format(dic["ptype"]))

    def set_xy_labels(self, ax, dic):

        if dic["ptype"] == "cartesian":
            if "v_n_x" in dic.keys():
                if "xlabel" in dic.keys():
                    ax.set_xlabel(dic["xlabel"], fontsize=12)
                else:
                    ax.set_xlabel(dic["v_n_x"].replace('_', '\_'), fontsize=12)

            if "v_n_y" in dic.keys():
                if "ylabel" in dic.keys():
                    ax.set_ylabel(dic["ylabel"], fontsize=12)
                else:
                    ax.set_ylabel(dic["v_n_y"].replace('_', '\_'), fontsize=12)


        elif dic["ptype"] == "polar":
            pass
        else:
            raise NameError("Unknown 'ptype' of the plot: {}".format(dic["ptype"]))

    def set_legend(self, ax, dic):
        if "legend" in dic.keys():
            if dic["legend"]:
                print("legend")
                ax.legend(fancybox=False, loc='best', shadow=False, fontsize=8)

    def remover_some_ticks(self, ax, dic):

        if "rmxlbls" in dic.keys():
            if dic["rmxlbls"]:
                ax.set_xticklabels([])
                ax.axes.xaxis.set_ticklabels([])

        if "rmylbls" in dic.keys():
            if dic["rmylbls"]:
                ax.set_yticklabels([])
                ax.axes.yaxis.set_ticklabels([])


    @staticmethod
    def plot_colormesh(ax, dic, x_arr, y_arr, z_arr):

        """

        int_ang_mom_flux_dic = {
            'ptype': 'polar',
            'v_n_x': 'phi_cyl', 'v_n_y': 'r_cyl', 'v_n': 'ang_mom_flux',
            'phimin': None, 'phimax': None, 'rmin': 0, 'rmax': 50, 'vmin': 1e-8, 'vmax': 1e-5,
            'mask_below': None, 'mask_above': None, 'cmap': 'RdBu_r', 'norm': "log"
        }

        corr_dic_ang_mom_flux_dens_unb_bern = {
            'ptype': 'cartesian',
            'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
            'xmin': 1e-11, 'xmax': 1e-7, 'ymin': 1e-11, 'ymax': 1e-7, 'vmin': 1e-7, 'vmax': None,
            'xscale': 'log', 'yscale': 'log',
            'mask_below': None, 'mask_above': None, 'cmap': 'inferno_r', 'norm': 'log'
        }

        :param ax:
        :param dic:
        :param x_arr:
        :param y_arr:
        :param z_arr:
        :return:
        """


        if "mask" in dic.keys():
            if dic["mask"] == "negative":
                z_arr = np.ma.masked_array(z_arr, z_arr < 0) # [z_arr < 0] = np.nan#
                print(z_arr)
            elif dic["mask"] == "positive":
                z_arr = -1 * np.ma.masked_array(z_arr, z_arr > 0)

        if "mask_below" in dic.keys():
            z_arr = np.ma.masked_array(z_arr, z_arr < dic["mask_below"])
        elif "mask_above" in dic.keys():
            z_arr = np.ma.masked_array(z_arr, z_arr > dic["mask_below"])


        if dic["vmin"] == None: dic["vmin"] = z_arr.min()
        if dic["vmax"] == None: dic["vmax"] = z_arr.max()

        if dic["norm"] == "norm" or dic["norm"] == "linear" or dic["norm"] == None:
            norm = Normalize(vmin=dic["vmin"], vmax=dic["vmax"])
        elif dic["norm"] == "log":
            norm = LogNorm(vmin=dic["vmin"], vmax=dic["vmax"])
        else:
            raise NameError("unrecognized norm: {} in task {}"
                            .format(dic["norm"], dic["v_n"]))



        im = ax.pcolormesh(x_arr, y_arr, z_arr, norm=norm, cmap=dic["cmap"])#, vmin=dic["vmin"], vmax=dic["vmax"])
        im.set_rasterized(True)

        return im

    def fill_arr_with_vmin(self, arr, dic):

        if 'fill_vmin' in dic.keys():
            if dic['fill_vmin']:
                if 'vmin' in dic.keys():
                    if dic['vmin'] != None:
                        arr = np.maximum(arr, dic['vmin'])
        return arr

    def add_fancy_to_ax(self, ax, dic):

        if 'fancyticks' in dic.keys():
            if dic['fancyticks']:
                ax.tick_params(axis='both', which='both', labelleft=True,
                               labelright=False, tick1On=True, tick2On=True,
                               labelsize=12, direction='in')

        if 'minorticks' in dic.keys():
            if dic["minorticks"]:
                ax.minorticks_on()

        if 'v_n_x' in dic.keys() and 'v_n_y' in dic.keys():
            if dic['v_n_x'] == 'hist_theta' and dic['v_n_y'] == 'hist_theta_m':
                xmajorticks = np.arange(5) * 90. / 4.
                xminorticks = np.arange(17) * 90. / 16
                xmajorlabels = [r"$0^\circ$", r"$22.5^\circ$", r"$45^\circ$",
                                r"$67.5^\circ$", r"$90^\circ$"]
                ax.xaxis.set_major_locator(FixedLocator(xmajorticks))
                ax.xaxis.set_minor_locator(FixedLocator(xminorticks))
                ax.set_xticklabels(xmajorlabels)

            if dic['v_n_x'] == 'hist_ye' and dic['v_n_y'] == 'hist_ye_m':
                if 'sharey' in dic.keys() and dic['sharey']:
                    ax.set_xticks(np.arange(0.1, 0.6, .1))
                else:
                    ax.xaxis.set_major_locator(MultipleLocator(0.1))
                    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

            if dic['v_n_x'] == 'hist_vel_inf' and dic['v_n_y'] == 'hist_vel_inf_m':
                if 'sharey' in dic.keys() and dic['sharey']:
                    ax.set_xticks(np.arange(0.2, 1.2, .2))
                else:
                    ax.set_xticks(np.arange(0, 1.2, .2))

            if dic['v_n_x'] == 'A' and dic['v_n_y'] == 'Y_final':
                ax.set_xticks(np.arange(50, 200, 50))

        if 'sharey' in dic.keys():
            if dic['sharey']:
                ax.set_yticklabels([])
                ax.axes.yaxis.set_ticklabels([])
                ax.set_ylabel('')
                ax.tick_params(labelleft=False)

        if 'sharex' in dic.keys():
            if dic['sharex']:
                ax.set_xticklabels([])
                ax.axes.xaxis.set_ticklabels([])
                ax.set_xlabel('')
                ax.tick_params(labelbottom=False)

        # if 'invert_x' in dic.keys():
        #     if dic['invert_x']:
        #         ax.axes.invert_xaxis()
        #
        # if 'invert_y' in dic.keys():
        #     if dic['invert_y']:
        #         print("revertin!")
        #         ax = plt.gca()
        #         ax.set_ylim(ax.get_ylim()[::-1])
        #         ax.invert_yaxis()
        #         ax.axes.invert_yaxis()
        #         plt.gca().invert_yaxis()

    def plot_generic_vertical_line(self, ax, dic):

        value = dic['value']

        if 'label' in dic.keys():
            if dic['label'] == None:
                ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'])
            else:
                # exit(1)
                ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'], label=dic['label'])
        else:
            ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'])

        return 0

    def plot_generic_line(self, ax, dic, x_arr, y_arr):


        if 'alpha' in dic.keys():
            if 'marker' in dic.keys():
                if 'label' in dic.keys():
                    print("alpha marker label")
                    ax.plot(x_arr, y_arr, dic['marker'],  color=dic['color'], markersize=dic['ms'], alpha=dic['alpha'], label=dic['label'])
                else:
                    print("alpha marker")
                    ax.plot(x_arr, y_arr, dic['marker'], color=dic['color'], markersize=dic['ms'], alpha=dic['alpha'])

            elif 'ls' in dic.keys():
                if 'ds' in dic.keys():
                    if 'label' in dic.keys():
                        print("alpha ls ds label")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'],  drawstyle=dic['ds'], alpha=dic['alpha'], label=dic['label'])
                    else:
                        print("alpha ls ds")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], drawstyle=dic['ds'], alpha=dic['alpha'])
                else:
                    if 'label' in dic.keys():
                        print("alpha ls label")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], alpha=dic['alpha'], label=dic['label'])
                        print("ls")
                    else:
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], alpha=dic['alpha'])
            else:
                raise NameError("what am I: marker or a line?")
        else:
            if 'marker' in dic.keys():
                if 'label' in dic.keys():
                    print("marker label")
                    ax.plot(x_arr, y_arr, dic['marker'], color=dic['color'], markersize=dic['ms'], label=dic['label'])
                    print("marker")
                else:
                    ax.plot(x_arr, y_arr, dic['marker'], color=dic['color'], markersize=dic['ms'])

            elif 'ls' in dic.keys():
                if 'ds' in dic.keys():
                    if 'label' in dic.keys():
                        print("ls ds label")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], drawstyle=dic['ds'], label=dic['label'])
                    else:
                        print("ls ds")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], drawstyle=dic['ds'])
                else:
                    if 'label' in dic.keys():
                        print("ls label")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'], label=dic['label'])
                    else:
                        print("ls")
                        ax.plot(x_arr, y_arr, ls=dic['ls'], color=dic['color'])
            else:
                raise NameError("what am I: marker or a line?")

    def plot_generic_errorbar(self, ax, dic, x_arr, y_arr, yerr):


        if 'label' in dic.keys() and dic['label'] != None:
            ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'],
                        alpha=dic['alpha'], label=dic['label'])
        else:
            ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'],
                        alpha=dic['alpha'])

        #
        #
        # if 'alpha' in dic.keys():
        #     if 'marker' in dic.keys():
        #         if 'label' in dic.keys():
        #             print("alpha marker label")
        #             ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'], alpha=dic['alpha'], label=dic['label'])
        #         else:
        #             print("alpha marker")
        #             ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'], alpha=dic['alpha'])
        #
        # else:
        #     if 'marker' in dic.keys():
        #         if 'label' in dic.keys():
        #             print("marker label")
        #             ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'], label=dic['label'])
        #             print("marker")
        #         else:
        #             ax.errorbar(x_arr, y_arr, yerr=yerr, fmt=dic['marker'], color=dic['color'], markersize=dic['ms'])
        #     else:
        #         raise NameError("what am I: marker or a line?")

    def plot_generic_band(self, ax, dic, x_arr, y1_arr, y2_arr):

        assert len(y1_arr) == len(x_arr)
        assert len(y2_arr) == len(x_arr)

        if dic["label"] == None:
            ax.fill_between(x_arr, y1_arr, y2_arr, alpha=dic['alpha'], color=dic['color'])
        else:
            ax.fill_between(x_arr, y1_arr, y2_arr, alpha=dic['alpha'], color=dic['color'], label=dic['label'])

    def treat_time_acis(self, x_arr, dic):

        if '-t' in dic.keys():
            # print("tmerger: {}".format(dic['-t']))
            x_arr = x_arr - float(dic['-t'])

        # print("x[0]-tmrg:{}".format(x_arr[0]))
        if 'xunits' in dic.keys():
            if dic['xunits'] == 's':
                pass
            elif dic['xunits'] == 'ms':
                x_arr = x_arr * 1e3
            else:
                raise NameError("x_units {} not recognized".format(dic['xunits']))

        return x_arr

    def treat_mass_acis(self, y_arr, dic):

        if 'yunits' in dic.keys():
            if dic['yunits'] == '1e-2Msun':
                return y_arr * 1e2

class PLOT_TASK(BASIC_PARTS):

    def __init__(self):
        BASIC_PARTS.__init__(self)

    def plot_2d_generic_colormesh(self, ax, dic):

        if dic['dtype'] == "corr":

            o_data = dic["data"]
            table = o_data.get_res_corr(dic["it"], dic["v_n_x"], dic["v_n_y"])
            table = np.array(table)
            x_arr = np.array(table[0, 1:])  # * 6.176269145886162e+17
            y_arr = np.array(table[1:, 0])
            z_arr = np.array(table[1:, 1:])

            z_arr = z_arr / np.sum(z_arr)
            z_arr = np.maximum(z_arr, z_arr.min())

            im = self.plot_colormesh(ax, dic, x_arr, y_arr, z_arr)

        else:
            raise NameError("plot type dic['dtype'] is not recognised (given: {})".format(dic["dtype"]))

        return im

    def plot_2d_projection(self, ax, dic):

        """

        dic = {
            'dtype': 'corr',
            'it': it, 'v_n_x': 'dens_unb_bern', 'v_n_y': 'ang_mom_flux', 'v_n': 'mass',
        }

        :param ax:
        :param dic:
        :return:
        """


        if dic['dtype'] == 'int':

            o_data = dic["data"]

            x_arr, y_arr, z_arr = o_data.get_modified_2d_data(
                dic["it"], dic["v_n_x"], dic["v_n_y"], dic["v_n"], dic["mod"]
            )
            #
            #
            # y_arr = o_data.get_grid_data(dic["it"], dic["v_n_x"])  # phi
            # x_arr = o_data.get_grid_data(dic["it"], dic["v_n_y"])  # r
            # z_arr =
            #
            #
            # if 'mod' in dic.keys() and dic['mod'] != None:
            #     if dic['mod'] == 'integ_over_z':
            #         x_arr = np.array(x_arr[:, 0, 0])
            #         y_arr = np.array(y_arr[0, :, 0])
            #         z_arr = o_data.get_integ_over_z(dic['it'], dic['v_n'])
            #     elif dic['mod'] == 'integ_over_z int':
            #         x_arr = np.array(x_arr[:, 0, 0]) # r
            #         y_arr = np.array(y_arr[0, :, 0]) # phi
            #         z_arr = o_data.get_integ_over_z(dic['it'], dic['v_n'])
            #         # print(np.rad2deg(y_arr)); exit(1)
            #         print(x_arr[-1], y_arr[-1])
            #         print(x_arr.shape, y_arr.shape, z_arr.shape)
            #         y_arr = np.append(y_arr, 2*np.pi)
            #         z_arr = np.vstack((z_arr.T, z_arr[:, -1])).T
            #         y_arr = np.ins(y_arr, 2 * np.pi)
            #         z_arr = np.vstack((z_arr[:, 0], z_arr.T)).T
            #         print(x_arr.shape, y_arr.shape, z_arr.shape)
            #         # from scipy import interpolate
            #         # grid_x = np.linspace(0.0, 2.0 * np.pi, 360)
            #         # grid_y = np.linspace(dic['rmin'], dic['rmax'], 200)
            #         # print(x_arr.shape, y_arr.shape, z_arr.shape)
            #         # X, Y = np.meshgrid(x_arr, y_arr)
            #         # f_ = interpolate.interp2d(X, Y, z_arr)
            #         # z_arr = f_(grid_x, grid_y)
            #     else:
            #         raise NameError("Unknown 'mod' parameter:{} ".format(dic['mod']))
            # else:
            #     x_arr = np.array(x_arr[:, 0, 0])
            #     y_arr = np.array(y_arr[0, :, 0])
            #     z_arr = o_data.get_int_data(dic["it"], dic["v_n"])
            #     z_arr = np.array(z_arr[:, :, 0])  # take a slice
            #
            #
            #
            #
            #
            #
            # # z_arr = self.fill_arr_with_vmin(z_arr, dic)
            #
            # # print(x_arr.shape, y_arr.shape, z_arr.shape)
            # # print(y_arr)

            im = self.plot_colormesh(ax, dic, y_arr, x_arr, z_arr)  # phi, r, data

        elif dic['dtype'] == 'dm':

            o_data = dic['data']
            im = 0
            iterations = o_data.get_grid('iterations')
            r_cyl = o_data.get_grid(dic['v_n_y'])

            if dic['v_n'] == 'int_phi':
                int_phi2d = o_data.get_data(dic['mode'], dic['v_n'])
                int_phi2d_for_it = int_phi2d[int(np.where(iterations == dic['it'])[0]), :]
                x_arr = np.angle(int_phi2d_for_it)[r_cyl < dic['rmax']]
                y_arr = r_cyl[r_cyl < dic['rmax']]

                if 'int' in dic.keys():
                    if dic['int'] == 'spline':
                        from scipy import interpolate
                        y_grid = np.mgrid[y_arr[0]:y_arr[-1]:1000j]
                        x_grid = interpolate.interp1d(y_arr, x_arr, kind='linear')(y_grid)
                        ax.plot(x_grid, y_grid, dic['ls'], color=dic['color'])
                else:
                    ax.plot(x_arr, y_arr, dic['ls'], color=dic['color'])

            elif dic['v_n'] == 'int_phi_r':
                int_phi_r1d = o_data.get_data(dic['mode'], dic['v_n'])
                int_phi_r1d_for_it = int_phi_r1d[iterations == dic['it']]
                phi = np.zeros(r_cyl.shape)
                phi.fill(float(np.angle(int_phi_r1d_for_it)))
                ax.plot(phi[r_cyl < dic['rmax']], r_cyl[r_cyl < dic['rmax']], dic['ls'], color=dic['color'])

            else:
                raise NameError("dic['v_n'] is not recognized. Use 'int_phi' or 'int_phi_r' ")

        else:
            raise NameError("plot type dic['dtype'] is not recognised (given: {})".format(dic["dtype"]))

        return im

    def plot_density_mode_line(self, ax, dic):


        if dic['dtype'] == 'dm':
            o_data = dic['data']
            im = 0
            x_arr = o_data.get_grid(dic['v_n_x'])
            # print("x[0]:{}".format(x_arr[0]))

            x_arr = self.treat_time_acis(x_arr, dic)

            # print("x[0]-tmrg * 1e3:{}".format(x_arr[0]))

            if dic['v_n_y'] == 'int_phi_r abs':
                int_phi_r1d = o_data.get_data(dic['mode'], 'int_phi_r')
                y_arr = np.abs(int_phi_r1d)
                if 'norm_to_m' in dic.keys():
                    # print('Normalizing')
                    norm_int_phi_r1d = o_data.get_data(dic['norm_to_m'], 'int_phi_r')
                    # print(norm_int_phi_r1d); exit(1)
                    y_arr = y_arr / abs(norm_int_phi_r1d)[0]

            elif dic['v_n_y'] == 'int_phi_r phase':
                int_phi_r1d = o_data.get_data(dic['mode'], 'int_phi_r')
                y_arr = np.unwrap(np.angle(int_phi_r1d))
                if 'norm_to_m' in dic.keys() and dic['norm_to_m'] != None:
                    raise NameError("cannot normalize phase of the mode")
            else:
                raise NameError("v_n_y {} is not recognized".format(dic['v_n_y']))



            self.plot_generic_line(ax, dic, x_arr, y_arr)

            # if dic['label'] == 'mode':
            #     self.plot_generic_line(ax, dic, x_arr, y_arr)

                # ax.plot(x_arr, y_arr, ls=dic['ls'], lw=dic['lw'], color=dic['color'],
                #         label='m:{}'.format(dic['mode']))
            # elif dic['label'] == 'sim':
                # sim = o_data.sim
                # ax.plot(x_arr, y_arr, ls=dic['ls'], lw=dic['lw'], color=dic['color'],
                #                 label='{}'.format(sim.replace('_', '\_')))
            # else:
            #     ax.plot(x_arr, y_arr, ls=dic['ls'], lw=dic['lw'], color=dic['color'])
            del(x_arr)
        else:
            raise NameError("plot type dic['dtype'] is not recognised (given: {})".format(dic["dtype"]))

    def plot_tcoll_vert_line(self, ax, dic):

        o_data = dic['data']
        try:
            value = o_data.get_par("tcoll_gw")
            value = self.treat_time_acis(value, dic)

            # print(value); exit(1)
            if 'label' in dic.keys():
                if dic['label'] != None:
                    ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'])
                else:
                    ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'], label=dic['label'])
            else:
                ax.axvline(x=value, linestyle=dic['ls'], color=dic['color'], linewidth=dic['lw'])
        except:
            print("Warning! tcoll failed to be plotted")

    def plot_histogram_1d(self, ax, dic):

        o_data = dic['data']
        if dic['v_n_x'] == 'hist_theta' and dic['v_n_y'] == 'hist_theta_m':
            tht = o_data.get_arr('hist_theta', dic['criterion'])
            M = o_data.get_arr('hist_theta_m', dic['criterion'])

            if dic['norm']: M /= np.sum(M)

            # ax.step(90. - (tht / np.pi * 180.), M, color=dic['color'], where='mid', label=dic['label'])
            # ax.plot(90. - (tht / np.pi * 180.), M, color=dic['color'], ls=dic['ls'], drawstyle=dic['ds'], label=dic['label'])
            self.plot_generic_line(ax, dic, 90. - (tht / np.pi * 180.), M)
            dtht = tht[1] - tht[0]


            ax.set_xlim(xmin=0 - dtht / np.pi * 180, xmax=90.)
            # xmajorticks = np.arange(5) * 90. / 4.
            # xminorticks = np.arange(17) * 90. / 16
            # xmajorlabels = [r"$0^\circ$", r"$22.5^\circ$", r"$45^\circ$",
            #                 r"$67.5^\circ$", r"$90^\circ$"]
            # ax.xaxis.set_major_locator(FixedLocator(xmajorticks))
            # ax.xaxis.set_minor_locator(FixedLocator(xminorticks))
            # ax.set_xticklabels(xmajorlabels)
            #
            # ax.set_xlabel(r"Angle from orbital plane")

        elif dic['v_n_x'] == 'hist_ye' and dic['v_n_y'] == 'hist_ye_m':
            o_data = dic['data']
            ye = o_data.get_arr('hist_ye', dic['criterion'])
            M = o_data.get_arr('hist_ye_m', dic['criterion'])

            if dic['norm']: M /= np.sum(M)

            # ax.step(ye, M, color=dic['color'], where='mid', label=dic['label'])

            # ax.plot(ye, M, color=dic['color'], ls=dic['ls'], drawstyle=dic['ds'], label=dic['label'])
            self.plot_generic_line(ax, dic, ye, M)

        elif dic['v_n_x'] == 'hist_vel_inf' and dic['v_n_y'] == 'hist_vel_inf_m':

            o_data = dic['data']
            vel_inf = o_data.get_arr('hist_vel_inf', dic['criterion'])
            M = o_data.get_arr('hist_vel_inf_m', dic['criterion'])
            # print(M)
            if dic['norm']: M /= np.sum(M)

            # print(M)
            # ax.step(ye, M, color=dic['color'], where='mid', label=dic['label'])
            # ax.plot(vel_inf, M, color=dic['color'], ls=dic['ls'], drawstyle=dic['ds'], label=dic['label'])
            self.plot_generic_line(ax, dic, vel_inf, M)

    def plot_ejecta_profile(self, ax, dic):

        o_data = dic['data']

        x_arr = o_data.get_arr(dic['v_n_x'], dic['criterion'])
        y_arr = o_data.get_arr(dic['v_n_y'], dic['criterion'])

        x_arr = self.treat_time_acis(x_arr, dic)
        y_arr = self.treat_mass_acis(y_arr, dic)

        # if 'xunits' in dic.keys():
        #     if dic['xunits'] == 'ms':
        #         x_arr *= 1e3

        self.plot_generic_line(ax, dic, x_arr, y_arr)

        # ax.plot(x_arr, y_arr, ls = dic['ls'], color=dic['color'], lw = dic['lw'])

    def plot_nucleo_yeilds_line(self, ax, dic):

        o_data = dic['data']
        if dic['v_n_x'] in o_data.list_sims_v_ns:
            x_arr = o_data.get_normalized_sim_data(dic['v_n_x'], criterion=dic['criterion'], method=dic['method'])
        elif dic['v_n_x'] in o_data.list_sol_v_ns:
            x_arr = o_data.get_normalized_sol_data(dic['v_n_x'], method=dic['method'])
        else:
            raise NameError("v_n_x:{} is not in available v_ns lists"
                            .format(dic["v_n_x"]))

        o_data = dic['data']
        if dic['v_n_y'] in o_data.list_sims_v_ns:
            y_arr = o_data.get_normalized_sim_data(dic['v_n_y'], criterion=dic['criterion'], method=dic['method'])
        elif dic['v_n_y'] in o_data.list_sol_v_ns:
            y_arr = o_data.get_normalized_sol_data(dic['v_n_y'], method=dic['method'])
        else:
            raise NameError("v_n_y:{} is not in available v_ns lists"
                            .format(dic["v_n_y"]))

        self.plot_generic_line(ax, dic, x_arr, y_arr)

    def plot_mkn_lightcurve(self, ax, dic):

        data = dic['data']
        m_time, m_min, m_max = data.get_model_min_max(dic['band'], fname=dic['fname'])
        if dic["label"] == None:
            ax.fill_between(m_time, m_min, m_max, alpha=dic['alpha'], color=dic['color'])
        else:
            ax.fill_between(m_time, m_min, m_max, alpha=dic['alpha'], color=dic['color'], label=dic['label'])

    def plot_mkn_obs_data(self, ax, dic):

        data = dic['data']
        data_list = data.get_obs_data(dic["band"], fname="AT2017gfo.h5")

        for i_, arr in enumerate(data_list):
            if dic["label"] == None:
                self.plot_generic_errorbar(ax, dic, arr[:, 0], arr[:, 1], yerr=arr[:, 2])
            else:
                if i_ == 0:
                    self.plot_generic_errorbar(ax, dic, arr[:, 0], arr[:, 1], yerr=arr[:, 2])
                else:
                    dic['label'] = None
                    self.plot_generic_errorbar(ax, dic, arr[:, 0], arr[:, 1], yerr=arr[:, 2])

    def plot_ejecta_band_2_objects(self, ax, dic):

        o_data1 = dic["data1"]
        o_data2 = dic["data2"]

        tmerg1 = o_data1.get_par("tmerger_gw")
        time_arr1 = o_data1.get_arr(dic["v_n_x"], dic["criterion1"])
        time_arr1 = time_arr1 - tmerg1
        mass_flux_arr1 = o_data1.get_arr(dic["v_n_y"], dic["criterion1"])
        mass_flux_arr1 = self.treat_mass_acis(mass_flux_arr1, dic)

        tmerg2 = o_data2.get_par("tmerger_gw")
        time_arr2 = o_data2.get_arr(dic["v_n_x"], dic["criterion2"])
        time_arr2 = time_arr2 - tmerg2
        mass_flux_arr2 = o_data2.get_arr(dic["v_n_y"], dic["criterion2"])
        mass_flux_arr2 = self.treat_mass_acis(mass_flux_arr2, dic)

        from scipy import interpolate

        print(time_arr1[-1], time_arr2[-1])

        if time_arr2[-1] > time_arr1[-1]:

            mass_flux_arr2 = interpolate.interp1d(time_arr2, mass_flux_arr2, kind='linear', bounds_error=False)(time_arr1)
            # print(len(time_arr1))
            time_arr1 = self.treat_time_acis(time_arr1, dic)
            self.plot_generic_band(ax, dic, time_arr1, mass_flux_arr1, mass_flux_arr2)
        else:
            #  time_arr2[-1] < time_arr1[-1]
            mass_flux_arr1 = interpolate.interp1d(time_arr1, mass_flux_arr1, kind='linear', bounds_error=False)(time_arr2)
            time_arr2 = self.treat_time_acis(time_arr2, dic)
            self.plot_generic_band(ax, dic, time_arr2, mass_flux_arr1, mass_flux_arr2)


        return 0


    def plot_task(self, ax, dic):

        if dic["task"] == '2d projection':
            return self.plot_2d_projection(ax, dic)
        elif dic["task"] == '2d colormesh':
            return self.plot_2d_generic_colormesh(ax, dic)
        elif dic["task"] == 'line':
            return self.plot_density_mode_line(ax, dic)
        elif dic['task'] == 'vertline':
            return self.plot_generic_vertical_line(ax, dic)
        elif dic['task'] == 'hist1d':
            return self.plot_histogram_1d(ax, dic)
        elif dic['task'] == 'ejprof':
            return self.plot_ejecta_profile(ax, dic)
        elif dic['task'] == 'ejband':
            return self.plot_ejecta_band_2_objects(ax, dic)
        elif dic['task'] == 'nucleo':
            return self.plot_nucleo_yeilds_line(ax, dic)
        elif dic['task'] == 'mkn model':
            return self.plot_mkn_lightcurve(ax, dic)
        elif dic['task'] == 'mkn obs':
            return self.plot_mkn_obs_data(ax, dic)
        else:
            raise NameError("dic['task'] is not recognized ({})".format(dic["task"]))


class PLOT_MANY_TASKS(PLOT_TASK):

    def __init__(self):

        PLOT_TASK.__init__(self)

        self.gen_set = {
            "figdir": './',
            "figname": "rename_me.png",
            # "figsize": (13.5, 3.5), # <->, |
            "figsize": (3.8, 3.5),  # <->, |
            "type": "cartesian",
            "subplots_adjust_h": 0.2,
            "subplots_adjust_w": 0.3,
            "fancy_ticks": True,
            "minorticks_on": True,
            "invert_y": False,
            "invert_x": False
        }

        self.set_plot_dics = []

    def set_ncols_nrows(self):

        tmp_rows = []
        tmp_cols = []

        for dic in self.set_plot_dics:
            tmp_cols.append(dic['position'][1])
            tmp_rows.append(dic['position'][0])

        max_row = max(tmp_rows)
        max_col = max(tmp_cols)

        for row in range(1, max_row):
            if not row in tmp_rows:
                raise NameError("Please set vertical plot position in a subsequent order: 1,2,3... not 1,3...")

        for col in range(1, max_col):
            if not col in tmp_cols:
                raise NameError("Please set horizontal plot position in a subsequent order: 1,2,3... not 1,3...")

        print("\tSet {} rows {} columns (total {}) of plots".format(max_row, max_col, len(self.set_plot_dics)))

        return int(max_row), int(max_col)

    def set_plot_dics_matrix(self):

        plot_dic_matrix = [[0
                             for x in range(self.n_rows)]
                             for y in range(self.n_cols)]

        # get a matrix of dictionaries describing plots (for ease of representation)
        for dic in self.set_plot_dics:
            col, row = int(dic['position'][1]-1), int(dic['position'][0]-1) # -1 as position starts with 1
            # print(col, row)
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):
                    if int(col) == int(n_col) and int(row) == int(n_row):
                        plot_dic_matrix[n_col][n_row] = dic
                        # print('adding {} {}'.format(col, row))

            if isinstance(plot_dic_matrix[col][row], int):
                raise ValueError("Dictionary to found for n_row {} n_col {} in "
                                 "creating matrix of dictionaries".format(col, row))

        return plot_dic_matrix

    def set_plot_matrix(self):

        fig = plt.figure(figsize=self.gen_set['figsize'])  # (<->; v)

        if self.gen_set['type'] == 'cartesian':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1)
                              for x in range(self.n_rows)]
                             for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    elif n_col == 0 and n_row > 0:
                        if self.gen_set['sharex']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharex=sbplot_matrix[n_col][0])
                        else:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    elif n_col > 0 and n_row == 0:
                        if self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharey=sbplot_matrix[0][n_row])
                        else:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i)
                    else:
                        if self.gen_set['sharex'] and not self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharex=sbplot_matrix[n_col][0])
                        elif not self.gen_set['sharex'] and self.gen_set['sharey']:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharey=sbplot_matrix[0][n_row])
                        else:
                            sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i,
                                                                          sharex=sbplot_matrix[n_col][0],
                                                                          sharey=sbplot_matrix[0][n_row])

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
                    i += 1

        elif self.gen_set['type'] == 'polar':
            # initializing the matrix with dummy axis objects
            sbplot_matrix = [[fig.add_subplot(self.n_rows, self.n_cols, 1, projection='polar')
                                  for x in range(self.n_rows)]
                                  for y in range(self.n_cols)]

            i = 1
            for n_row in range(self.n_rows):
                for n_col in range(self.n_cols):

                    if n_col == 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                    elif n_col == 0 and n_row > 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0])
                    elif n_col > 0 and n_row == 0:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharey=self.sbplot_matrix[0][n_row])
                    else:
                        sbplot_matrix[n_col][n_row] = fig.add_subplot(self.n_rows, self.n_cols, i, projection='polar')
                                                                      # sharex=self.sbplot_matrix[n_col][0],
                                                                      # sharey=self.sbplot_matrix[0][n_row])
                    i += 1

                        # sbplot_matrix[n_col][n_row].axes.get_yaxis().set_visible(False)
                    # sbplot_matrix[n_col][n_row] = fig.add_subplot(n_rows, n_cols, i)
        else:
            raise NameError("type of the plot is not recognized. Use 'polar' or 'cartesian' ")

        return fig, sbplot_matrix

    def plot_images(self):

        # initializing the matrix of images for colorbars (if needed)
        image_matrix = [[0
                        for x in range(self.n_rows)]
                        for y in range(self.n_cols)]


        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                for dic in self.set_plot_dics:
                    if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
                        print("\tPlotting n_row:{} n_col:{}".format(n_row, n_col))
                        ax = self.sbplot_matrix[n_col][n_row]
                        # dic = self.plot_dic_matrix[n_col][n_row]
                        if isinstance(dic, int):
                            print("Warning: Dictionary for row:{} col:{} not set".format(n_row, n_col))
                            self.fig.delaxes(ax)  # delets the axis for empty plot
                        else:
                            dic = dict(dic)
                            im = self.plot_task(ax, dic)

                            if not isinstance(im, int):
                                image_matrix[n_col][n_row] = im

                            self.set_plot_title(ax, dic)
                            # self.account_for_shared(ax, n_row, n_col)
                            self.set_min_max_scale(ax, dic)
                            self.set_xy_labels(ax, dic)
                            self.set_legend(ax, dic)
                            self.remover_some_ticks(ax, dic)
                            self.add_fancy_to_ax(ax, dic)


        return image_matrix

    def account_for_shared(self, ax, n_row, n_col):

        ax.axes.xaxis.set_ticklabels([])
        ax.axes.yaxis.set_ticklabels([])

        if n_col > 0 and n_row < self.n_rows:
            # ax.tick_params(labelbottom=False)
            # ax.tick_params(labelleft=False)
            #
            # ax.set_yticklabels([])
            # ax.set_xticklabels([])
            #
            # ax.get_yaxis().set_ticks([])
            #
            # ax.set_yticks([])
            # ax.set_yticklabels(labels=[])
            # # ax.set_yticklabels([]).remove()

            # ax.axes.get_xaxis().set_visible(False)
            # ax.axes.get_yaxis().set_visible(False)

            ax.axes.xaxis.set_ticklabels([])
            ax.axes.yaxis.set_ticklabels([])

            # ax.tick_params(labelbottom=False)
            # ax.tick_params(labelleft=False)
            # ax.tick_params(labelright=False)

            # ax.get_yaxis().set_visible(False)

        if n_col > 0:
            ax.set_ylabel('')

        if n_row != self.n_rows-1:
            ax.set_xlabel('')

    def plot_one_cbar(self, im, dic, n_row, n_col):

        if "cbar" in dic.keys():
            if dic["cbar"] != None and dic["cbar"] != '':

                location = dic["cbar"].split(' ')[0]
                shift_h = float(dic["cbar"].split(' ')[1])
                shift_w = float(dic["cbar"].split(' ')[2])
                cbar_width = 0.02


                if location == 'right':
                    ax_to_use = self.sbplot_matrix[-1][n_row]
                    pos1 = ax_to_use.get_position()
                    pos2 = [pos1.x0 + pos1.width + shift_h,
                            pos1.y0,
                            cbar_width,
                            pos1.height]
                elif location == 'left':
                    ax_to_use = self.sbplot_matrix[-1][n_row]
                    pos1 = ax_to_use.get_position()
                    pos2 = [pos1.x0 - pos1.width - shift_h,
                            pos1.y0,
                            cbar_width,
                            pos1.height]
                elif location == 'bottom':
                    ax_to_use = self.sbplot_matrix[n_col][-1]
                    pos1 = ax_to_use.get_position()
                    pos2 = [pos1.x0,
                            pos1.y0 - pos1.height + shift_w,
                            cbar_width,
                            pos1.height]
                else:
                    raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                    .format(location))

                cax1 = self.fig.add_axes(pos2)
                if location == 'right':
                    if 'cbar fmt' in dic.keys() and dic['cbar fmt'] != None:
                        cbar = plt.colorbar(im, cax=cax1, extend='both', format=dic['cbar fmt'])
                    else:
                        cbar = plt.colorbar(im, cax=cax1, extend='both')  # , format='%.1e')
                elif location == 'left':
                    if 'cbar fmt' in dic.keys() and dic['cbar fmt'] != None:
                        cbar = plt.colorbar(im, cax=cax1, extend='both', format=dic['cbar fmt'])#, format='%.1e')
                    else:
                        cbar = plt.colorbar(im, cax=cax1, extend='both')
                    cax1.yaxis.set_ticks_position('left')
                    cax1.yaxis.set_label_position('left')
                else:
                    raise NameError("cbar location {} not recognized. Use 'right' or 'bottom' "
                                    .format(location))
                if 'cbar label' in dic.keys() and dic['cbar label'] != None:
                    cbar.ax.set_title(dic['cbar label'])

                else:
                    cbar.ax.set_title(r"{}".format(str(dic["v_n"]).replace('_', '\_')))

    def plot_colobars(self):

        for n_row in range(self.n_rows):
            for n_col in range(self.n_cols):
                for dic in self.set_plot_dics:
                    if n_col + 1 == int(dic['position'][1]) and n_row + 1 == int(dic['position'][0]):
                        print("\tColobar for n_row:{} n_col:{}".format(n_row, n_col))
                        # ax  = self.sbplot_matrix[n_col][n_row]
                        # dic = self.plot_dic_matrix[n_col][n_row]
                        im  = self.image_matrix[n_col][n_row]
                        if isinstance(dic, int):
                            print("Warning: Dictionary for row:{} col:{} not set".format(n_row, n_col))
                        else:
                            self.plot_one_cbar(im, dic, n_row, n_col)


        # for n_row in range(self.n_rows):
        #     for n_col in range(self.n_cols):
        #         print("Colobar for n_row:{} n_col:{}".format(n_row, n_col))
        #         # ax  = self.sbplot_matrix[n_col][n_row]
        #         dic = self.plot_dic_matrix[n_col][n_row]
        #         im  = self.image_matrix[n_col][n_row]
        #         if isinstance(dic, int):
        #             Printcolor.yellow("Dictionary for row:{} col:{} not set".format(n_row, n_col))
        #         else:
        #             self.plot_one_cbar(im, dic, n_row, n_col)

    def save_plot(self):

        if self.gen_set["invert_y"]:
            plt.gca().invert_yaxis()

        if self.gen_set["invert_x"]:
            plt.gca().invert_xaxis()

        plt.subplots_adjust(hspace=self.gen_set["subplots_adjust_h"])
        plt.subplots_adjust(wspace=self.gen_set["subplots_adjust_w"])
        # plt.tight_layout()
        plt.savefig('{}{}'.format(self.gen_set["figdir"], self.gen_set["figname"]),
                    bbox_inches='tight', dpi=128)
        plt.close()

    def main(self):

        if len(self.set_plot_dics) == 0:
            raise ValueError("No plot dics have been passed. Exiting")

        # initializing the n_cols, n_rows
        self.n_rows, self.n_cols = self.set_ncols_nrows()
        # initializing the matrix of dictionaries of the
        self.plot_dic_matrix = self.set_plot_dics_matrix()
        # initializing the axis matrix (for all subplots) and image matrix fo colorbars
        self.fig, self.sbplot_matrix = self.set_plot_matrix()
        # plotting
        self.image_matrix = self.plot_images()
        # adding colobars
        self.plot_colobars()

        # saving the result
        self.save_plot()