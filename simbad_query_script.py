# coding: utf-8
"""Author: Sasha Safonova
This scrapes the SIMBAD database (http://simbad.u-strasbg.fr) to find any
documented parallaxes for stars in your list.

The input for the script is a csv file named "parallax_query_objects.csv",
placed in the same directory as the script. The csv file should contain three
columns without headers:
    - object name, formatted to your wishes
    - RA in the hour angle format
    - DEC in the degree angle format (plus or minus sign mandatory)

The output of the script is printed to your screen, saved as a numpy array in
a .npy file, and written to a neat csv file.

The code can be modified to look for a quantity other than parallax by
replacing the 'plx' string input, noted inside the code.

The full list of attributes with which SIMBAD works (they may or may not work
with the Simbad Python module) is here:
http://simbad.u-strasbg.fr/simbad/sim-fsam
"""

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import pandas as pd
import astropy.units as u
import numpy as np
import datetime


customSimbad = Simbad()

# if you wish to look for a quantity other than parallax, replace 'plx' with
# that quantity's SIMBAD code below:
customSimbad.add_votable_fields('ra(d)', 'dec(d)', 'plx')
customSimbad.remove_votable_fields('coordinates')
customSimbad.add_votable_fields()

# ----read in the csv and convert to a list of coordinate strings----
star_list = pd.read_csv("parallax_query_objects.csv",
                        header=None,
                        names=['Object_name', 'ra', 'dec', None])
objects = [str(star_list['ra'].as_matrix()[i]) +
           ' ' +
           str(star_list['dec'].as_matrix()[i])
           for i in np.arange(len(star_list['ra'].as_matrix()))]
parallax_objects = pd.DataFrame(columns=['MAIN_ID',
                                         'RA_d',
                                         'DEC_d',
                                         'FOUND_VALUE',
                                         'readable_ra',
                                         'readable_dec'])
# ----loop over the objects provided in "parallax_query_objects.csv"-----
for coords in objects:
    c = SkyCoord(coords, FK5, obstime="J2000", unit=(u.hourangle, u.deg))
    result = customSimbad.query_region(c, radius='0.05 degrees')
    df = result.to_pandas()
    these_have_parallaxes = pd.DataFrame(columns=['MAIN_ID',
                                                  'RA_d',
                                                  'DEC_d',
                                                  'FOUND_VALUE',
                                                  'readable_ra',
                                                  'readable_dec'])
    # this step strips away stars without measured parallaxes
    these_have_parallaxes = pd.concat([these_have_parallaxes,
                                       df[df['FOUND_VALUE'].notnull()]])

    for i in these_have_parallaxes.index:
        addthis = FK5((these_have_parallaxes.loc[i, 'RA_d'] * u.deg),
                      these_have_parallaxes.loc[i, 'DEC_d'] * u.deg)
        these_have_parallaxes['readable_ra'][i]  = addthis.ra.to_string(u.hour)
        these_have_parallaxes['readable_dec'][i] = addthis.dec.to_string(u.deg)
    parallax_objects = pd.concat([parallax_objects, these_have_parallaxes])
    del df
    del these_have_parallaxes

# ----write out the objects with found parallaxes in a variety of formats----
print("parallax objects:")
print(parallax_objects)
np.save("parallax_objects_values.npy", parallax_objects.values)
now = datetime.datetime.now()
parallax_objects.to_csv("parallax_objects{0}.csv".format(now.strftime("%B_%d_%Y_%H:%M:%S")),
                        columns=['MAIN_ID',
                                 'FOUND_VALUE',
                                 'readable_ra',
                                 'readable_dec'],
                        header=['Object name',
                                'Parallax (mas)',
                                'RA',
                                'DEC'],
                        index=False)
