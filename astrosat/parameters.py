import numpy
import ephem
from astropy.coordinates import AltAz, SkyCoord
from astropy import units as u
from astropy.coordinates import EarthLocation
import datetime
import yaml
import logging


class Parameters:
    def __init__(self, parameter_file):

        self.parameter_file = parameter_file
        self.parameters_dict = self.load_parameters()
        logging.basicConfig(filename=self.parameters_dict['logfile'], encoding='utf-8', level=self.parameters_dict['log_level'])

        self.cross_section = float(self.parameters_dict['cross_section'])
        self.radius = float(self.parameters_dict['radius'])
        self.duration = float(self.parameters_dict['duration'])
        self.Mmin = float(self.parameters_dict['Mmin'])
        self.Mmax = float(self.parameters_dict['Mmax'])
        self.plotDir = self.parameters_dict['plotDir']
        self.fdir = self.parameters_dict['fdir']
        self.TLEdir = self.parameters_dict['TLEdir']

        self.obs = ephem.Observer()
        self.obs.lon = str(self.parameters_dict['lon'])
        self.obs.lat = str(self.parameters_dict['lat'])
        self.obs.elevation = self.parameters_dict['alt']

        self.date = datetime.datetime(year=self.parameters_dict['year'], month=self.parameters_dict['month'],
                                      day=self.parameters_dict['day'], hour=self.parameters_dict['hour'],
                                      minute=self.parameters_dict['minute'], second=self.parameters_dict['seconds'])
        
        self.location = EarthLocation.from_geodetic(self.obs.lon, self.obs.lat, height=self.obs.elevation)

        if 'RA' and 'DEC' in self.parameters_dict:
            logging.debug('Ra/Dec provided')
            ra = self.parameters_dict['RA'].split(':')
            dec = self.parameters_dict['DEC'].split(':')
            self.RA = (numpy.array(ra).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()*180./12.
            self.DEC = (numpy.array(dec).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
        elif 'ALT' and 'AZ' in self.parameters_dict:
            logging.debug('Alt/Az provided, converting to Ra/Dec')
            alt = self.parameters_dict['ALT'].split(':')
            az = self.parameters_dict['AZ'].split(':')
            alt = (numpy.array(alt).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
            az = (numpy.array(az).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()

            # Use astropy to convert co-ordinates
            c = SkyCoord(AltAz(alt=alt*u.degree, az=az*u.degree, location=self.location, obstime=self.date))
            pointing = c.transform_to('icrs')
            self.RA = pointing.ra.value
            self.DEC = pointing.dec.value

        else:
            raise TypeError('Observing pointing must be provided in either RA/DEC, or ALT/AZ format')

        self.verbose = self.parameters_dict['verbose']
        self.logfile = self.parameters_dict['logfile']
        self.log_level = self.parameters_dict['log_level']

    def load_parameters(self):
        with open(self.parameter_file, 'r') as stream:
            try:
                parameters_dict = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        return parameters_dict
    
    def update_time(self):
        if 'RA' and 'DEC' in self.parameters_dict:
            logging.debug('Ra/Dec provided')
            ra = self.parameters_dict['RA'].split(':')
            dec = self.parameters_dict['DEC'].split(':')
            self.RA = (numpy.array(ra).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()*180./12.
            self.DEC = (numpy.array(dec).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
        elif 'ALT' and 'AZ' in self.parameters_dict:
            logging.debug('Alt/Az provided, converting to Ra/Dec')
            alt = self.parameters_dict['ALT'].split(':')
            az = self.parameters_dict['AZ'].split(':')
            alt = (numpy.array(alt).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
            az = (numpy.array(az).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()

            # Use astropy to convert co-ordinates
            c = SkyCoord(AltAz(alt=alt*u.degree, az=az*u.degree, location=self.location, obstime=self.date))
            pointing = c.transform_to('icrs')
            self.RA = pointing.ra.value
            self.DEC = pointing.dec.value

        else:
            raise TypeError('Observing pointing must be provided in either RA/DEC, or ALT/AZ format')
