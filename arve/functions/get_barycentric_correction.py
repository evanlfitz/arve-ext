from barycorrpy import get_BC_vel, JDUTC_to_BJDTDB, get_stellar_data
from astropy.time import Time

class get_barycentric_correction:
  
  def get_barycentric_correction(
    self, 
    StarName : string,
    OrderJD : float
    ) -> float:
          '''
          # StarName must be Simbad queryable
          # specfile is the spectrum file
          # Order is the specific order you want the correction for, Order 19 is for He10830
          # HPF Orders range from 0 to 27
          # Returns BaryVel in km/s
          '''
          stardatadic, warnings = get_stellar_data(StarName)
          #OrderJD = #float(fits.getval(specfile,'JD_FW{0}'.format(Order)))
          BaryVel = get_BC_vel(JDUTC=Time(OrderJD,format='jd'),lat=30.681,longi=-104.0147,alt=2025,**stardatadic)[0][0]
          #bjd = JDUTC_to_BJDTDB(JDUTC=Time(OrderJD,format='jd'),lat=30.681,longi=-104.0147,alt=2025,**stardatadic)[0][0]
          return BaryVel/1e3  # return in km/s
