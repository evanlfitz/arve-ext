from barycorrpy import get_BC_vel, JDUTC_to_BJDTDB
from astropy.time import Time
import numpy as np
try:
    from astroquery.simbad import Simbad
except:
    print('Cannot import astroquery.simbad')
import astropy.units as u
from astropy.coordinates import SkyCoord
import os

class get_barycentric_correction:

  def get_barycentric_correction(
    self, 
    OrderJD : float | np.ndarray,
    ) -> float | np.ndarray:
          '''
          # StarName must be Simbad queryable
          # specfile is the spectrum file
          # Order is the specific order you want the correction for, Order 19 is for He10830
          # HPF Orders range from 0 to 27
          # Returns BaryVel in km/s
          '''
          stardatadic, warnings = _get_stellar_data(self.arve.star.target)
          #OrderJD = #float(fits.getval(specfile,'JD_FW{0}'.format(Order)))
          BaryVel = np.zeros(len(OrderJD))
          for i,date in enumerate(OrderJD):
            BaryVel[i] = get_BC_vel(JDUTC=Time(date,format='jd'),lat=30.681,longi=-104.0147,alt=2025,**stardatadic)[0][0]
          #bjd = JDUTC_to_BJDTDB(JDUTC=Time(OrderJD,format='jd'),lat=30.681,longi=-104.0147,alt=2025,**stardatadic)[0][0]
          return BaryVel/1e3  # return in km/s
  
def _get_stellar_data(name=''):
  '''
  Comes from barycorrpy, modified to pull RA and DEC in the new way by querying RA and DEC and getting it back in degrees in the ICRS frame and changing the SkyCoord call accordingly
  Function to query Simbad for following stellar information RA, Dec, PMRA, PMDec, Parallax Epoch
  INPUTS:
      name = Name of source. Example 
  
  
  '''
  warning = []

  #print(name)
  
  customSimbad = Simbad()
  customSimbad.add_votable_fields('ra', 'dec','pmra','pmdec', 'plx_value','rvz_radvel')
  #Simbad.list_votable_fields()
  # customSimbad.remove_votable_fields( 'coordinates')
  #Simbad.get_field_description('orv')
  obj = customSimbad.query_object(name)
  if obj is None:
      raise ValueError('ERROR: {} target not found. Check target name or enter RA,Dec,PMRA,PMDec,Plx,RV,Epoch manually\n\n'.format(name))
  else:        
      warning += ['{} queried from SIMBAD.'.format(name)]

  # Check for masked values
  if all([not x for x in [obj.mask[0][i] for i in obj.colnames]])==False:
      warning += ['Masked values present in queried dataset']


  obj = obj.filled(None)
  
#   pos = SkyCoord(ra=obj['ra'][0]*u.deg,dec=obj['dec'][0]*u.deg,frame='icrs')
  ra = obj['ra'].value[0] #pos.ra.value[0]
  dec = obj['dec'].value[0] #pos.dec.value[0]
  pmra = obj['pmra'][0]
  pmdec = obj['pmdec'][0]
  plx = obj['plx_value'][0]
  rv = obj['rvz_radvel'][0] * 1000 #SIMBAD output is in km/s. Converting to m/s
  epoch = 2451545.0
  
  star = {'ra':ra,'dec':dec,'pmra':pmra,'pmdec':pmdec,'px':plx,'rv':rv,'epoch':epoch}
  
  # Fill Masked values with None. Again. 
  for i in star:
      if star[i] > 1e10:
          star[i] = None           
      
  warning += ['Values queried from SIMBAD are {}'.format(star)]

  
  return star,warning

