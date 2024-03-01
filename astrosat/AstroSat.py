# authors: James Osborn, james.osborn@durham.ac.uk
#          Matthew Townson, matthew.townson@northumbria.ac.uk

import numpy
import ephem
import datetime
import math
import scipy.interpolate
import os
import warnings
from rich.progress import track
from urllib.request import urlopen
import logging
warnings.filterwarnings("ignore")

from .parameters import Parameters


class AstroSat:
    def __init__(self, parameters):
        if type(parameters) is str:
            self.parameters = Parameters(parameters)
        else:
            self.parameters = parameters
        self.search_radius = numpy.max([15.,self.parameters.radius])
        self.sun = ephem.Sun()

        logging.basicConfig(filename=self.parameters.logfile, encoding='utf-8', level=self.parameters.log_level)

    def get_TLEs(self, satellite_type='ACTIVE', forceNew=0):
        '''
        Download TLEs
        check all times take closest one, if time now is closer to obs date download new

        '''
        timeStamp = self.parameters.date.timestamp()
        timeNow = datetime.datetime.now().timestamp()
        timeTempList=[]

        # Create list of archived TLE files
        fnlist = os.listdir(self.parameters.TLEdir)
        fnlist.sort()
        for fn in fnlist:
            if '.dat' in fn:
                timeTemp = int(fn.split('_')[0])
                typeTemp = fn.split('_')[1][:-4]
                if typeTemp == satellite_type:
                    timeTempList.append(timeTemp)

        # Apply cleaning rules to decide if new TLE file is needed
        if len(timeTempList) == 0:
            forceNew = 1
        elif len(timeTempList) == 1:
            timeTemp = timeTempList[0]
            if abs(timeStamp - timeTemp) > abs(timeStamp - timeNow):
                # if archived file is further in time to now download new TLE
                forceNew = 1
        else:
            # find closest match list in time
            timeTemp = timeTempList[numpy.argmin(abs(timeStamp-numpy.array(timeTempList)))]
            if abs(timeStamp - timeTemp) > abs(timeStamp - timeNow):
                # if archived file is further in time to now download new TLE
                forceNew = 1

        # Load TLEs from file if possible
        if forceNew == 0:
            fn = '%i_%s.dat' % (timeTemp,satellite_type)
            if self.parameters.verbose:
                print(f'Loading TLE:{satellite_type} from file: {fn}')
            satTLEs = []
            with open(f'{self.parameters.TLEdir}/{fn}', 'r') as f:
                satTLEs.extend(line.split(',') for line in f)
        else:
            satTLEs = self.download_TLEs(satellite_type)
        return satTLEs

    def download_TLEs(self, satellite_type, format='TLE'):
        if self.parameters.verbose:
            print(f'Downloading TLE:{satellite_type}')
        TLE_URL = f'https://celestrak.org/NORAD/elements/gp.php?GROUP={satellite_type}&FORMAT={format}'
        celestrak = urlopen(TLE_URL)
        TLEs = [item.strip() for item in celestrak]
        celestrak.close()

        # Clean raw TLE data
        result = [(TLEs[i].decode('utf-8'), TLEs[i+1].decode('utf-8'), TLEs[i+2].decode('utf-8')) for i in numpy.arange(0, len(TLEs)-2, 3)]

        # write TLE data to file
        fn = '%i_%s.dat' %(datetime.datetime.now().timestamp(), satellite_type)
        with open(f'{self.parameters.TLEdir}/{fn}', 'w') as f:
            for TLE in result:
                f.write('%s,%s,%s\n' %(TLE[0],TLE[1],TLE[2]))

        return result

    def get_satellites(self,satTLEs):
        self.sats={}
        for tle in satTLEs:
            self.sats[tle[0]] = ephem.readtle(tle[0], tle[1], tle[2])
        return self.sats

    def process_satellite(self, sat, date, Fmodel=None, satDict=None):
        '''
        Parse sat into satellite orbits and project onto observers RA/DEC
        '''

        if satDict is None:
            satDict = {}

        self.parameters.obs.date = date

        sat.compute(self.parameters.obs)

        RA_angle_diff = (self.parameters.RA - sat.ra*180/numpy.pi + 180 + 360) % 360 - 180
        DEC_angle_diff = (self.parameters.DEC - sat.dec*180/numpy.pi + 180 + 360) % 360 - 180

        if abs(RA_angle_diff) < self.search_radius and abs(DEC_angle_diff) < self.search_radius:
            self.sun.compute(self.parameters.obs)
            
            if Fmodel == 'diffuseSpherical':
                # solar phase angle
                a = self.sun.earth_distance * 1.496e+11  # distance sun from observer (Km)
                b = sat.range*1. 
                if self.parameters.obs.elevation < -1000.:
                    b -= ephem.earth_radius
                angle_c = ephem.separation((sat.az, sat.alt), (self.sun.az, self.sun.alt))
                c = math.sqrt(math.pow(a, 2) + math.pow(b, 2) - 2*a*b*math.cos(angle_c))
                angle_a = math.acos((math.pow(b, 2) + math.pow(c, 2) - math.pow(a, 2)) / (2 * b * c))
                phase_angle = angle_a 
        
                # "Optical Tracking and Spectral Characterization of Cubesats for Operational Missions"
                # Gasdia, Forrest, (2016). PhD Dissertations and Master's Theses. 212.
                # https://commons.erau.edu/edt/212
                F = (2/(3*numpy.pi**2))*(numpy.sin(phase_angle) + (numpy.pi-phase_angle)*numpy.cos(phase_angle))
        
            elif Fmodel is None:
                # no solar phase angle dependence
                F = 1
        
            # extinction due to airmass
            if self.parameters.obs.elevation < -1000:
                gamma = 0.
            else:
                # https://www.aanda.org/articles/aa/full_html/2020/04/aa37501-20/aa37501-20.html
                # include curvature of Earth in airmass (for low elevation angles)
                # gamma = 0.12*1./(numpy.sin(sat.alt)+0.15*(sat.alt*180./numpy.pi + 3.885)**(-1.253))
        
                # calculate directly from TLE
                gamma = 0.12 * sat.range/sat.elevation
                
            # "Optical Tracking and Spectral Characterization of Cubesats for Operational Missions"
            # Gasdia, Forrest, (2016). PhD Dissertations and Master's Theses. 212.
            # https://commons.erau.edu/edt/212
            mag_sat = -26.74 - 2.5*numpy.log10((self.parameters.cross_section*F)/sat.range**2) + gamma
        
            if sat.eclipsed:
                # in shadow
                mag_sat = None
            elif mag_sat<0.:
                # Too bright - geometry gone wrong (small angles)
                mag_sat = None
            elif sat.alt<0.:
                # below horizon
                if self.parameters.obs.elevation > -1000:
                    # not if centre of Earth observer
                    mag_sat = None
        
            if self.parameters.duration > 1:

                if sat.name not in satDict.keys():
                    satDict[sat.name] = {}
                    satDict[sat.name]['RA'] = []
                    satDict[sat.name]['DEC'] = []
                    satDict[sat.name]['MAG'] = []
                    satDict[sat.name]['Time'] = []
                    satDict[sat.name]['sunElev'] = []
                    satDict[sat.name]['ALT'] = []
                    satDict[sat.name]['AZ'] = []
                    satDict[sat.name]['RANGE'] = []
                    satDict[sat.name]['ELEV'] = []

                satDict[sat.name]['RA'].append(sat.ra*12./numpy.pi)
                satDict[sat.name]['DEC'].append(sat.dec*180/numpy.pi)
                satDict[sat.name]['MAG'].append(mag_sat)
                satDict[sat.name]['Time'].append(date)
                satDict[sat.name]['ALT'].append(sat.alt*180./numpy.pi)
                satDict[sat.name]['AZ'].append(sat.az*180./numpy.pi)
                satDict[sat.name]['sunElev'].append(self.sun.alt*180./numpy.pi)
                satDict[sat.name]['RANGE'].append(sat.range)
                satDict[sat.name]['ELEV'].append(sat.elevation)
            else:
                if sat.name not in satDict.keys():
                    satDict[sat.name] = {}
                    satDict[sat.name]['RA'] = []
                    satDict[sat.name]['DEC'] = []
                    satDict[sat.name]['MAG'] = []
                    satDict[sat.name]['Time'] = []
                    satDict[sat.name]['sunElev'] = []
                    satDict[sat.name]['ALT'] = []
                    satDict[sat.name]['AZ'] = []
                    satDict[sat.name]['RANGE'] = []
                    satDict[sat.name]['ELEV'] = []

                satDict[sat.name]['RA'] = (sat.ra*12./numpy.pi)
                satDict[sat.name]['DEC'] = (sat.dec*180/numpy.pi)
                satDict[sat.name]['MAG'] = (mag_sat)
                satDict[sat.name]['Time'] = (date)
                satDict[sat.name]['ALT'] = (sat.alt*180./numpy.pi)
                satDict[sat.name]['AZ'] = (sat.az*180./numpy.pi)
                satDict[sat.name]['sunElev'] = (self.sun.alt*180./numpy.pi)
                satDict[sat.name]['RANGE'] = (sat.range)
                satDict[sat.name]['ELEV'] = (sat.elevation)

        return satDict

    def print_satellite_dictionary(self, satellite_dictionary):
        '''
        Calculate and print satellites that intersect with field of view
        '''
        sat_table = []
        RAtemp = numpy.arange(self.parameters.RA-self.parameters.radius, self.parameters.RA+self.parameters.radius, 1./3600.)
        if self.parameters.verbose == 1:
            for i_sat in track(satellite_dictionary.keys(), description='printing satellites...'):
                sat_table = self.print_satellite_dictionary_loop(i_sat, satellite_dictionary, RAtemp, sat_table)

        else:
            for i_sat in satellite_dictionary.keys():
                sat_table = self.print_satellite_dictionary_loop(i_sat, satellite_dictionary, RAtemp, sat_table)

        if len(sat_table)>0:
            print("{: <30} {: <15} {: <15}{: <10}".format(*['Name', 'Time (UTC)', 'Duration (s)', 'Mag (V)']))
            for row in sat_table:
                print("{: <30} {: <15} {: <15.1f}{: <10.2f}".format(*row[:-2]))
        else:
            print("No satellite intercept predicted")
        
        return sat_table

    def print_satellite_dictionary_loop(self, i_sat, satellite_dictionary, RAtemp, sat_table):

        if type(satellite_dictionary[i_sat]['MAG']) == list:
            mag_sat = satellite_dictionary[i_sat]['MAG'][0]
        else:
            mag_sat = satellite_dictionary[i_sat]['MAG']
        if mag_sat is not None:
            RA_sat = satellite_dictionary[i_sat]['RA']
            DEC_sat = satellite_dictionary[i_sat]['DEC']
            if type(satellite_dictionary[i_sat]['Time']) == list:
                time_sat = []
                for timeTemp in satellite_dictionary[i_sat]['Time']:
                    time_sat.append(timeTemp.timestamp())

                # find transit time through FoV
                f = scipy.interpolate.interp1d(numpy.array(RA_sat)*180/12., DEC_sat, fill_value='extrapolate',kind='slinear')
                DECextrap = f(RAtemp)

                RA1 = RAtemp[numpy.argmin(abs(DECextrap-(self.parameters.DEC-self.parameters.radius)))]
                RA2 = RAtemp[numpy.argmin(abs(DECextrap-(self.parameters.DEC+self.parameters.radius)))]
                
                if RA1 < self.parameters.RA-self.parameters.radius:
                    RAmin = self.parameters.RA-self.parameters.radius
                elif RA1 > self.parameters.RA+self.parameters.radius:
                    RAmin = self.parameters.RA+self.parameters.radius
                else:
                    RAmin = RA1
                if RA2 < self.parameters.RA-self.parameters.radius:
                    RAmax = self.parameters.RA-self.parameters.radius
                elif RA2 > self.parameters.RA+self.parameters.radius:
                    RAmax = self.parameters.RA+self.parameters.radius
                else:
                    RAmax = RA2

                DEC1 = f(RAmin)
                DEC2 = f(RAmax)

                DECmin = numpy.max([DEC1,self.parameters.DEC-self.parameters.radius])
                DECmax = numpy.min([DEC2,self.parameters.DEC+self.parameters.radius])

                if DECmin != DECmax and RAmin != RAmax:
                    f2 = scipy.interpolate.interp1d(numpy.array(RA_sat)*180/12.,time_sat,fill_value='extrapolate')
                    timeRA = f2([RAmin,RAmax])
                    f3 = scipy.interpolate.interp1d(DEC_sat,time_sat,fill_value='extrapolate')
                    timeDEC = f3([DECmin,DECmax])

                    if numpy.array_equal(timeRA.astype(int), timeDEC.astype(int)):
                        # make sure the satellite is where we think it is!
                        try:
                            satDict = self.process_satellite(self.sats[i_sat],datetime.datetime.fromtimestamp(numpy.mean(timeRA)))
                            mag_sat=satDict[i_sat]['MAG'][0]

                            time0 = numpy.min(timeRA)
                            if self.parameters.date.timestamp()+self.parameters.duration > time0 > self.parameters.date.timestamp():
                                # intercept FOV
                                if satDict[i_sat]['RA'][0]*180/12.< RAmax and satDict[i_sat]['RA'][0]*180./12. > RAmin:
                                    if satDict[i_sat]['DEC'][0] < DECmax and satDict[i_sat]['DEC'][0] > DECmin:
                                        if mag_sat != None:
                                            sat_table.append([i_sat, datetime.datetime.strftime(datetime.datetime.fromtimestamp(time0),'%H:%M:%S'), abs(timeRA[1]-timeRA[0]), mag_sat, [RAmin,RAmax],[DECmin,DECmax]])
                        except Exception:
                            # satellite outside range
                            pass
        return sat_table

    def find_intercept_sats(self, Fmodel, timeStep=10):
        '''
        do recursive search -  default 10 sec steps to find small number potential sats and time window then 1 sec
        '''
        sats = self.sats
        cadence = 'low resolution'

        for timeStep in [timeStep, 1]:
            satDict = {}
            if self.parameters.verbose == 1:
                for istep in track(range(0, int(self.parameters.duration), timeStep), description=f'propagating satellites ({cadence})...'):
                    dateTemp = self.parameters.date+datetime.timedelta(seconds=istep)
                    for isat in sats.keys():
                        satDict = self.process_satellite(sats[isat], dateTemp, Fmodel, satDict)
            else:
                for istep in range(0, int(self.parameters.duration), timeStep):
                    dateTemp = self.parameters.date+datetime.timedelta(seconds=istep)
                    for isat in sats.keys():
                        satDict = self.process_satellite(sats[isat], dateTemp, Fmodel, satDict)
            sats = {sat_key:sats[sat_key] for sat_key in satDict.keys()}
            cadence = 'high resolution'

        return satDict