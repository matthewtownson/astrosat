# James Osborn
# Plot field of view of telescope with satellite trails overlaid
# Observer and Field parameters defined in params.py

from astrosat import AstroSat
from astrosat import Plot
from astrosat import Parameters
from astrosat import Stars

parameters = Parameters("params.yaml")
AS = AstroSat(parameters)
stars = Stars(parameters)

# timeStep = 5  # 5 second time steps 

# get satellite positions
# use active satellites and visual objects (rocket bodies)
satTLEs = AS.get_TLEs('ACTIVE')
sats = AS.get_satellites(satTLEs)

satDict = AS.find_intercept_sats(Fmodel=None)

# print satellite dictionary
sat_table = AS.print_satellite_dictionary(satDict)

# plot field
plot = Plot(AS, stars)
if len(sat_table)>0:
    plot.plot_satellites(sat_table)
plot.plot_stars()
plot.plot_legend([AS.parameters.Mmin, int(round((AS.parameters.Mmin + AS.parameters.Mmax) / 3.)), 2 * int(round((AS.parameters.Mmin + AS.parameters.Mmax) / 3.))])
plot.make_plot()
plot.save_plot('skyView_%i_%i.png' % (AS.parameters.radius, AS.parameters.date.timestamp()))

