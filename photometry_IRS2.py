from astropy.io import fits
import dendrocat
from astropy.table import Table, Column, MaskedColumn, vstack, hstack
import astropy.units as u
n1 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51n_band6_226GHz_adjustedRADEC.fits'))
n2 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51north_sci.spw0_1_2_3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18_19.mfs.I.manual.image.tt0.pbcor.fits'))
n3 = dendrocat.RadioSource(fits.open('/users/bmcclell/Data/irs2/w51n_bandQ_45GHz_adjustedRADEC.fits'))
n3.nu = 45*u.GHz
n3.freq_id = '45.0GHz'
n3.set_metadata()
n2.min_value = 1e-4
n2.min_delta = 1.5 * 1e-4
n1.autoreject()
n2.autoreject()
n2.reject([93001, 93002, 93003, 93004, 93005, 93008, 93009, 93010, 93030, 93031, 93032, 93050, 93052, 93062, 93081, 93092, 93099])
n1.reject([226003, 226004, 226026])
n1.accept([226047])
mc = dendrocat.utils.match(n1, n2)
mc.add_objects(n3)
mc.catalog['_name'][mc.catalog['_name'] == '226007'] = 'w51d2'
dendrocat.utils.saveregions(mc.catalog, '/users/bmcclell/nrao/reg/w51IRS2_photometered.reg')
mc.photometer(dendrocat.aperture.Ellipse)
mc.photometer(dendrocat.aperture.Circle)
mc.photometer(dendrocat.aperture.Annulus)
mc.catalog.write('/users/bmcclell/nrao/cat/w51IRS2_photometered.dat', format='ascii')


