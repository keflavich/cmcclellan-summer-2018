from astropy.table import Table, hstack
import dendrocat
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np




t_adam = Table.read('cat/core_continuum_and_line.ipac', format='ipac')
t_me = Table.read('cat/45-93-226GHz_photometered_adjustedRADEC.dat', format='ascii')
t_me = t_me[t_me['rejected'] == 0]

ra = dendrocat.utils.ucheck(t_adam['RA'], u.deg)
dec = dendrocat.utils.ucheck(t_adam['Dec'], u.deg)

myra = dendrocat.utils.ucheck(t_me['x_cen'], u.deg)
mydec = dendrocat.utils.ucheck(t_me['y_cen'], u.deg)

coords = SkyCoord(ra=ra, dec=dec)
mycoords = SkyCoord(ra=myra, dec=mydec)

tolerance = 0.3*u.arcsec
idx, d2d, d3d = mycoords.match_to_catalog_sky(coords)
matched = d2d < tolerance

t_adam = Table(t_adam[idx], masked=True)
#t_adam.mask[~matched] = True
combinedtable = hstack([t_me, t_adam])

#ind = dendrocat.utils.get_index_masked(combinedtable['SourceID'])
#combinedtable.remove_rows(ind)

combinedtable.write('cat/45-93-226GHz_photometered_matched_to_ALMA.dat', format='ascii')









t_adam = Table.read('cat/EVLA_VLA_PointSourcePhotometry.ipac', format='ipac')
t_me = Table.read('cat/45-93-226GHz_photometered_adjustedRADEC.dat', format='ascii')

ra = dendrocat.utils.ucheck(t_adam['gracen'], u.deg)
dec = dendrocat.utils.ucheck(t_adam['gdeccen'], u.deg)

myra = dendrocat.utils.ucheck(t_me['x_cen'], u.deg)
mydec = dendrocat.utils.ucheck(t_me['y_cen'], u.deg)

coords = SkyCoord(ra=ra, dec=dec, frame='fk5')
mycoords = SkyCoord(ra=myra, dec=mydec, frame='icrs')

tolerance = 0.1*u.arcsec
idx, d2d, d3d = mycoords.match_to_catalog_sky(coords)
matched = d2d < tolerance

t_adam = Table(t_adam[idx], masked=True)
t_adam.mask[~matched] = True
combinedtable = hstack([t_adam, t_me])

ind = dendrocat.utils.get_index_masked(combinedtable['SourceName'])
combinedtable.remove_rows(ind)

combinedtable.write('cat/45-93-226GHz_photometered_matched_to_EVLAVLA.dat', format='ascii')

