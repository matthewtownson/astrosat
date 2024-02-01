import ephem
import numpy
import astropy
import os
from astroquery.simbad import Simbad
import astropy.units as u

class Stars(object):
    def __init__(self, parameters) -> None:
        self.parameters = parameters
        self.search_radius = numpy.max([15.,self.parameters.radius])
        

        self.sun = ephem.Sun()

        self.bright_stars = self.get_bright_stars()
        

    def get_bright_stars(self):
            '''
            Extract bright star catalogue to python list
            '''
            stars=[]
            # find the stars
            if self.parameters.radius < 2.:
                # small FoV - use simbad
                if self.parameters.verbose:
                    print('Connecting to Simbad')
                useBSC=0
                try:
                    self.parameters.obs.date = self.parameters.date
                    star = ephem.FixedBody()
                    star._ra = self.parameters.RA*numpy.pi/180.
                    star._dec = self.parameters.DEC*numpy.pi/180
                    star.compute(self.parameters.obs)
                    RA_simbad, DEC_simbad = self.parameters.obs.radec_of(star.az, star.alt)
                    custom_simbad = Simbad()
                    custom_simbad.add_votable_fields('flux(V)')
                    pointing = astropy.coordinates.SkyCoord(ra=RA_simbad*u.rad, dec=DEC_simbad*u.rad)
                    starsSimbad = custom_simbad.query_region(pointing, radius=self.parameters.radius*1.5 * u.deg)
                    for row in starsSimbad:
                        if len(row[1].split()) == 3:
                            RAstar = (numpy.array(row[1].split()).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum() *360./24.
                        if len(row[1].split()) == 2:
                            RAstar = (numpy.array(row[1].split()).astype(float)*numpy.array([1, 1./60.])).sum()

                        if len(row[2].split()) == 3:
                            DECstar = (numpy.array(row[2].split()).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
                        if len(row[2].split()) == 2:
                            DECstar = (numpy.array(row[2].split()).astype(float)*numpy.array([1, 1./60.])).sum()
                        if numpy.isnan(row[11]) == False:
                            magstar = float(row[11])
                            stars.append([row[0].decode('utf-8'),RAstar,DECstar,magstar])

                except:
                    if self.parameters.verbose:
                        print('Connection to Simbad failed')
                    useBSC=1
                    # use BSC

            if self.parameters.radius >= 2. or useBSC ==1:
                # For large fields of view just use bright stars otherwise plots get cluttered     
                # 
                #   
                # 
                #  
                if self.parameters.verbose:
                    print('Using Bright Star Catalogue')
                with open(os.path.dirname(__file__)+"/data/bsc.dat",'r') as f:     
                    for line in f:
                        # Loop through the Yale Bright Star Catalog, line by line
                        # Ignore blank lines and comment lines
                        if (len(line) < 100) or (line[0] == '#'):
                            continue
                        try:
                            # Read the Henry Draper (i.e. HD) number for this star
                            hd = int(line[25:31])
                            # Read the right ascension of this star (J2000)
                            ra_hrs = float(line[75:77])
                            ra_min = float(line[77:79])
                            ra_sec = float(line[79:82])
                            # Read the declination of this star (J2000)
                            dec_neg = (line[83] == '-')
                            dec_deg = float(line[84:86])
                            dec_min = float(line[86:88])
                            dec_sec = float(line[88:90])
                            # Read the V magnitude of this star
                            mag = float(line[102:107])
                        except ValueError:
                            continue
                        # Turn RA and Dec from sexagesimal units into decimal
                        ra = (numpy.array([ra_hrs,ra_min,ra_sec]).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum() *360./24.
                        dec = (numpy.array([dec_deg,dec_min,dec_sec]).astype(float)*numpy.array([1, 1./60., 1./3600.])).sum()
                        if dec_neg:
                            dec = -dec
                        stars.append([hd,ra,dec,mag])
                
            return stars