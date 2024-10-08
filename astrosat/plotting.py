    
import numpy
import matplotlib.pylab as plt
import os
import astropy.coordinates
import astropy.units as u
import scipy.interpolate
import os


class Plot:
    def __init__(self, AS, stars):
        self.AS = AS
        self.stars = stars

        # setup plots
        try:
            os.mkdir(AS.parameters.fdir)
        except FileExistsError:
            pass
        try:
            os.mkdir(AS.parameters.TLEdir)
        except FileExistsError:
            pass
        # plt.style.use('dark_background')
        if AS.parameters.radius == 180:
            ## full sky
            self.markersize = 1.
        else:
            self.markersize = 5

        self.fig = plt.figure('Sky field', figsize=(7, 6))
        self.axSky = self.fig.add_subplot(111)  # , projection='polar')

    def plot_stars(self):
        # add stars to plot
        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        for i_d in range(len(self.stars.bright_stars)):
            RAstar = self.stars.bright_stars[i_d][1]
            DECstar = self.stars.bright_stars[i_d][2]
            RAanglediff = (self.AS.parameters.RA - RAstar + 180 + 360) % 360 - 180
            DECanglediff = (self.AS.parameters.DEC - DECstar + 180 + 360) % 360 - 180

            # over sample as field is square
            if abs(RAanglediff) < self.AS.parameters.radius*1.5:
                if abs(DECanglediff) < self.AS.parameters.radius*1.5:
                    mag_star = self.stars.bright_stars[i_d][3]
                    # r = 1+1*(Mmin-mag_star)/float(Mmax-Mmin)
                    r = -1*(Mmin-mag_star)/float(Mmax-Mmin)
                    r = numpy.clip(r, 0, 1)
                    plt.plot(RAstar*12./180., DECstar, color='k', marker='o', markersize=(r*-1+1)*self.markersize)


    def plot_satellites(self, sat_dict):
        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        # Add satellites to plot
        for i_sat in range(len(sat_dict)):
            mag_sat = sat_dict[i_sat][3]
            if mag_sat is not None:
                r = abs(-1*(Mmin-mag_sat)/float(Mmax-Mmin))
                r = numpy.clip(r, 1e-3, 0.9)
                if len(sat_dict[i_sat][4])==1:
                    plt.plot(sat_dict[i_sat][4], sat_dict[i_sat][5], 'x', color='red',markersize=(r*-1+1) * self.markersize)
                else:
                    # find transit time through FoV
                    f = scipy.interpolate.interp1d(numpy.array(sat_dict[i_sat][4]), sat_dict[i_sat][5], fill_value='extrapolate')
                    RAtemp = numpy.arange(self.AS.parameters.RA-self.AS.parameters.radius,
                                            self.AS.parameters.RA+self.AS.parameters.radius, 1./3600.)
                    DECextrap = f(RAtemp)
                    plt.plot(RAtemp*12./180., DECextrap, linestyle='--', color='red', linewidth=(r*-1+1)* self.markersize)

    def plot_legend(self, labels):
        Mmin = self.AS.parameters.Mmin
        Mmax = self.AS.parameters.Mmax
        # add legend to plot
        label_list = []
        p1, = plt.plot([],[], 'o', color=str(0), markersize=self.markersize*(0), label='stars')
        label_list.append(p1)
        for label in labels:
            # r = 1+1*(Mmin-label)/float(Mmax-Mmin)
            r = abs(-1*(Mmin-label)/float(Mmax-Mmin))
            r = numpy.clip(r, 0, 1.)
            p1, = plt.plot([],[], 'o', color=str(r), markersize=self.markersize*(r*-1+1), label='mag %i' % label)

            label_list.append(p1)
        p1, = plt.plot([],[], 'o', color=str(0), markersize=self.markersize*(0), label='satellites')
        label_list.append(p1)
        for label in labels:
            # r = 1+1*(Mmin-label)/float(Mmax-Mmin)
            r = abs(-1*(Mmin-label)/float(Mmax-Mmin))
            r = numpy.clip(r, 0, 1.)
            p2, = plt.plot([],[], color='r', linestyle='--', linewidth=(r*-1+1)*self.markersize, label='mag %i' % label)
            label_list.append(p2)
        for label in labels:
            # r = 1+1*(Mmin-label)/float(Mmax-Mmin)
            r = abs(-1*(Mmin-label)/float(Mmax-Mmin))
            r = numpy.clip(r, 0, 1.)
            p2, = plt.plot([],[], 'x',color='r', markersize=(r*-1+1)*self.markersize, label='mag %i' % label)
            label_list.append(p2)

        plt.legend(handles=label_list, numpoints=1, bbox_to_anchor=(1.05, 1), loc='upper left')

    def make_plot(self):
        # make axis scale and labels
        if self.AS.parameters.radius == 180:
            # ful sky
            self.axSky.axis([0,24,-90,90])
        else:
            self.axSky.axis([self.AS.parameters.RA*12./180. - self.AS.parameters.radius*24./360.,
                            self.AS.parameters.RA*12./180. + self.AS.parameters.radius*24./360.,
                            self.AS.parameters.DEC-self.AS.parameters.radius, self.AS.parameters.DEC+self.AS.parameters.radius])
        labels = self.axSky.xaxis.get_majorticklocs()
        xlabel_list = []
        for label in labels:
            angle = astropy.coordinates.Angle(float(label), u.hour)
            xlabel_list.append( angle.to_string(unit=u.hour, precision=0))
        labels = self.axSky.yaxis.get_majorticklocs()
        ylabel_list = []
        for label in labels:
            angle = astropy.coordinates.Angle(float(label),u.degree)
            ylabel_list.append(angle.to_string(unit=u.degree, sep=':', precision=0))
        self.axSky.set_xticklabels(xlabel_list, rotation=30)
        self.axSky.set_yticklabels(ylabel_list)
        plt.subplots_adjust(left=0.15, right=0.75, top=0.9, bottom=0.2)
        self.axSky.set_xlabel('RA')
        self.axSky.set_ylabel('DEC')

    def save_plot(self, file_name):
        if not os.path.isdir(self.AS.parameters.plotDir):
            os.mkdir(self.AS.parameters.plotDir)

        self.fig.savefig(self.AS.parameters.plotDir + file_name, bbox_inches='tight')
