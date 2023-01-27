import numpy as np 
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
import pandas as pd
import sys
import matplotlib.pyplot as plt
from astroplan import Observer, FixedTarget, is_observable, is_always_observable, months_observable, observability_table
from astropy.time import Time
from astropy.table import Table
from astroplan import AltitudeConstraint, ObservingBlock
from astropy.io import ascii
from pytz import timezone
from astroplan.plots import plot_sky, plot_schedule_airmass
from astroplan.constraints import TimeConstraint
from astroplan.scheduling import Transitioner, Schedule, SequentialScheduler, PriorityScheduler

mkt_location = {
        "elevation": 1086.599484882955,
        #"elevation": 0,
        "name": "MeerKAT",
        "longitude_unit": "degree",
        "latitude_unit": "degree",
        "latitude": -30.711055553291878,
        "elevation_unit": "meter",
        "longitude": 21.443888889697842,
        "source": "MEERKAT, used in timing mode.\n\n    The origin of this data is unknown but as of 2021 June 8 it agrees exactly with\n    the values used by TEMPO and TEMPO2.\nvia PINT",
	"timezone": "Africa/Johannesburg",
	"aliases": ["MeerKAT"]}

phase_cals = {'J0825-5010':SkyCoord('08:25:26.870000 -50:10:38.50000', 
                                    unit=(u.hourangle, u.deg)),
              'J0906-6829':SkyCoord('09:06:52.440000 -68:29:39.80000', 
                                    unit=(u.hourangle, u.deg)),
              'J1318-4620':SkyCoord('13:18:30.060000 -46:20:34.90000', 
                                    unit=(u.hourangle, u.deg)),
              'J1619-8418':SkyCoord('16:19:33.970000 -84:18:19.10000', 
                                    unit=(u.hourangle, u.deg)),
              'J1726-5529':SkyCoord('17:26:49.700000 -55:29:38.00000', 
                                    unit=(u.hourangle, u.deg)), 
              'J2232+1143':SkyCoord('22:32:36.410000 +11:43:50.90000', 
                                    unit=(u.hourangle, u.deg)),
              'J0828-3731':SkyCoord('08:28:05 -37:31:06',
                                    unit=(u.hourangle, u.deg)),
              'J1154-3505':SkyCoord('11:54:22 -35:05:2',
                                    unit=(u.hourangle, u.deg)),
              'J1323-4452':SkyCoord('13:23:04 -44:52:34',
                                    unit=(u.hourangle, u.deg)),
              'J1424-4913':SkyCoord('14:24:32 -49:13:50',
                                    unit=(u.hourangle, u.deg)),
              'J1733-1304':SkyCoord('17:33:03 -13:04:50',
                                    unit=(u.hourangle, u.deg)),
              'J1744-5144':SkyCoord('17:44:25 -51:44:44',
                                    unit=(u.hourangle, u.deg)),
              'J2158-1501':SkyCoord('22:32:36.410000 +11:43:50.90000',
                                    unit=(u.hourangle, u.deg)),
              }

flux_cal = {'J1939-6342':SkyCoord('19:39:25.0264 -63:42:45.625', 
                                    unit=(u.hourangle, u.deg)),
            'J0408-6545':SkyCoord('04:08:20.3782 -65:45:9.080', 
                                    unit=(u.hourangle, u.deg))}

pol_cal = {'3C286':SkyCoord('13:31:08.2879 +30:30:32.958', 
                                    unit=(u.hourangle, u.deg)),
        '3C138':SkyCoord('05:21:09.886021 +16:38:22.051220', 
                                    unit=(u.hourangle, u.deg))}


# MeerKAT Location
location = EarthLocation.from_geodetic(mkt_location['longitude'] * u.deg, \
    mkt_location['latitude'] * u.deg, mkt_location['elevation'] * u.m)


meerkat = Observer(name='MeerKAT',
               location=location,
               timezone=timezone('Africa/Johannesburg'),
               description="MeerKAT telescope in SA")



# Time range- test
time_range = Time(["2023-02-01 00:00", "2023-02-01 23:59"])



# Read targets, calibrators and a make a list of Skycoord objects
targets_df = pd.read_csv('MMGPS-SBAND_pointing_grid_final.csv')
targets_df = targets_df[['Pointing', 'RA', 'Dec']]
target_table = Table.from_pandas(targets_df)
targets_ra = targets_df['RA'].values
targets_dec = targets_df['Dec'].values


flux_cal_list = []
for key, value in flux_cal.items():
    flux_cal_list.append(FixedTarget(coord = value, name=key))


pol_cal_list = []
for key, value in pol_cal.items():
    pol_cal_list.append(FixedTarget(coord = value, name=key))

phase_cal_list = []
phase_cal_ra = []
phase_cal_dec = []
phase_cal_name = []


for key, value in phase_cals.items():
    phase_cal_list.append(FixedTarget(coord = value, name=key))
   



targets = [FixedTarget(coord=SkyCoord(ra=ra, dec=dec, unit=(u.hourangle, u.deg)), name=name)
            for name, ra, dec in target_table]





# Define exposure times for target and calibrators
target_exposure = 1257 * u.second
pol_calib_exposure = 300 * u.second
flux_calib_exposure = 300 * u.second
phase_cal_exposure = 300 * u.second

#Observing block
blocks = []
fixed_time = Time("2023-02-02 08:00:00")

rise_set_time_df_target = []
rise_set_time_df_phase_cal = []
rise_set_time_df_pol_cal = []
rise_set_time_df_flux_cal = []


rise_set_error_margin = 5*u.minute


for i in range(len(targets)):
    coordinates = targets[i].coord.to_string('hmsdms', sep=':')
    ra = coordinates.split(' ')[0]
    dec = coordinates.split(' ')[1]

    rise_time_utc_20_deg = meerkat.target_rise_time(fixed_time, targets[i], n_grid_points=10, horizon=20 * u.deg) 
    set_time_utc_20_deg = meerkat.target_set_time(fixed_time, targets[i], n_grid_points=10, horizon=20 * u.deg) 

    rise_time_utc_50_deg = meerkat.target_rise_time(fixed_time, targets[i], n_grid_points=10, horizon=50 * u.deg) 
    set_time_utc_50_deg = meerkat.target_set_time(fixed_time, targets[i], n_grid_points=10, horizon=50 * u.deg)
 
    if rise_time_utc_20_deg.mask: 
        rise_time_sidereal_20_deg =  np.ma.masked
        set_time_sidereal_20_deg = np.ma.masked

        if rise_time_utc_50_deg.mask:
            rise_time_sidereal_50_deg = np.ma.masked
            set_time_sidereal_50_deg = np.ma.masked
        rise_set_time_df_target.append([targets[i].name, ra, dec, rise_time_sidereal_20_deg, set_time_sidereal_20_deg, rise_time_sidereal_50_deg, set_time_sidereal_50_deg, 'Target'])

    else:
        rise_time_utc_20_deg += rise_set_error_margin
        set_time_utc_20_deg -= rise_set_error_margin

        # Rise set time Sidereal at 20 deg.
        rise_time_sidereal_20_deg = meerkat.local_sidereal_time(rise_time_utc_20_deg).to_string(sep=':')
        set_time_sidereal_20_deg = meerkat.local_sidereal_time(set_time_utc_20_deg).to_string(sep=':')
       
        # Rise set time Sidereal at 50 deg.
        rise_time_sidereal_50_deg = meerkat.local_sidereal_time(rise_time_utc_50_deg).to_string(sep=':')
        set_time_sidereal_50_deg = meerkat.local_sidereal_time(set_time_utc_50_deg).to_string(sep=':')
        
        rise_set_time_df_target.append([targets[i].name, ra, dec, rise_time_sidereal_20_deg, set_time_sidereal_20_deg, rise_time_sidereal_50_deg, set_time_sidereal_50_deg, 'Target'])

rise_set_time_target = pd.DataFrame(rise_set_time_df_target, columns=['Pointing', 'RA', 'DEC', \
   'Rise_time_sidereal_20_deg', 'Set_time_sidereal_20_deg', 'Rise_time_sidereal_50_deg', 'Set_time_sidereal_50_deg', 'Obs_type'])
rise_set_time_target['Priority'] = 1
rise_set_time_target.to_csv('rise_set_time_target.csv', index=False)

for i in range(len(phase_cal_list)):

    coordinates = phase_cal_list[i].coord.to_string('hmsdms', sep=':')
    ra = coordinates.split(' ')[0]
    dec = coordinates.split(' ')[1]
    
    rise_time_utc = meerkat.target_rise_time(fixed_time, phase_cal_list[i], n_grid_points=10, horizon=20 * u.deg)
    set_time_utc = meerkat.target_set_time(fixed_time, phase_cal_list[i], n_grid_points=10, horizon=20 * u.deg) 
 
    if rise_time_utc.mask: 
        #rise_time_sidereal =  np.ma.masked
        #set_time_sidereal = np.ma.masked
        rise_time_sidereal = '0:00:00.0000000'
        set_time_sidereal = '23:59:59.99999999'
        rise_set_time_df_phase_cal.append([phase_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Phase_Cal'])

    else:
        rise_time_utc += rise_set_error_margin
        set_time_utc -= rise_set_error_margin
        rise_time_sidereal = meerkat.local_sidereal_time(rise_time_utc).to_string(sep=':') 
        set_time_sidereal = meerkat.local_sidereal_time(set_time_utc).to_string(sep=':')
        rise_set_time_df_phase_cal.append([phase_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Phase_Cal'])

rise_set_time_phase_cal = pd.DataFrame(rise_set_time_df_phase_cal, columns=['Pointing', 'RA', 'DEC', 'Rise_time_sidereal_20_deg', 'Set_time_sidereal_20_deg', 'Obs_type'])
rise_set_time_phase_cal['Priority'] = 1
rise_set_time_phase_cal.to_csv('rise_set_time_phase_cal.csv', index=False)

for i in range(len(flux_cal_list)):

    coordinates = flux_cal_list[i].coord.to_string('hmsdms', sep=':')
    ra = coordinates.split(' ')[0]
    dec = coordinates.split(' ')[1]
    
    rise_time_utc = meerkat.target_rise_time(fixed_time, flux_cal_list[i], n_grid_points=10, horizon=20 * u.deg)
    set_time_utc = meerkat.target_set_time(fixed_time, flux_cal_list[i], n_grid_points=10, horizon=20 * u.deg) 
 
    if rise_time_utc.mask: 
        rise_time_sidereal =  np.ma.masked
        set_time_sidereal = np.ma.masked
        rise_set_time_df_flux_cal.append([flux_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Flux_Cal'])


    else:
        rise_time_utc += rise_set_error_margin
        set_time_utc -= rise_set_error_margin
        rise_time_sidereal = meerkat.local_sidereal_time(rise_time_utc).to_string(sep=':') 
        set_time_sidereal = meerkat.local_sidereal_time(set_time_utc).to_string(sep=':')
        rise_set_time_df_flux_cal.append([flux_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Flux_Cal'])

rise_set_time_flux_cal = pd.DataFrame(rise_set_time_df_flux_cal, columns=['Pointing', 'RA', 'DEC', 'Rise_time_sidereal_20_deg', 'Set_time_sidereal_20_deg', 'Obs_type'])
rise_set_time_flux_cal['Priority'] = 1
rise_set_time_flux_cal.to_csv('rise_set_time_flux_cal.csv', index=False)

for i in range(len(pol_cal_list)):

    coordinates = pol_cal_list[i].coord.to_string('hmsdms', sep=':')
    ra = coordinates.split(' ')[0]
    dec = coordinates.split(' ')[1]
    
    rise_time_utc = meerkat.target_rise_time(fixed_time, pol_cal_list[i], n_grid_points=10, horizon=20 * u.deg)
    set_time_utc = meerkat.target_set_time(fixed_time, pol_cal_list[i], n_grid_points=10, horizon=20 * u.deg) 
 
    if rise_time_utc.mask: 
        rise_time_sidereal =  np.ma.masked
        set_time_sidereal = np.ma.masked
        rise_set_time_df_pol_cal.append([pol_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Pol_Cal'])

    else:
        rise_time_utc += rise_set_error_margin
        set_time_utc -= rise_set_error_margin
        rise_time_sidereal = meerkat.local_sidereal_time(rise_time_utc).to_string(sep=':') 
        set_time_sidereal = meerkat.local_sidereal_time(set_time_utc).to_string(sep=':')
        rise_set_time_df_pol_cal.append([pol_cal_list[i].name, ra, dec, rise_time_sidereal, set_time_sidereal, 'Pol_Cal'])       

rise_set_time_pol_cal = pd.DataFrame(rise_set_time_df_pol_cal, columns=['Pointing', 'RA', 'DEC', 'Rise_time_sidereal_20_deg', 'Set_time_sidereal_20_deg', 'Obs_type'])
rise_set_time_pol_cal['Priority'] = 1
rise_set_time_pol_cal.to_csv('rise_set_time_pol_cal.csv', index=False)


# Update target table with nearest phase calibrator

target_ra = rise_set_time_target['RA'].values
target_dec = rise_set_time_target['DEC'].values
phase_cal_ra = rise_set_time_phase_cal['RA'].values
phase_cal_dec = rise_set_time_phase_cal['DEC'].values
phase_cal_name = rise_set_time_phase_cal['Pointing'].values

target_catalog = SkyCoord(ra = target_ra, dec = target_dec, unit=(u.hourangle, u.deg))
phase_cal_catalog = SkyCoord(ra = phase_cal_ra, dec = phase_cal_dec, unit=(u.hourangle, u.deg))

idx, d2d, d3d = target_catalog.match_to_catalog_sky(phase_cal_catalog)

rise_set_time_target['Nearest_Phase_Cal'] = phase_cal_name[idx]
rise_set_time_target['Phase_Cal_Distance_Deg'] = d2d.deg
rise_set_time_target.to_csv('rise_set_time_target.csv', index=False)

print(rise_set_time_target)
