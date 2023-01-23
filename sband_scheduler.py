import numpy as np 
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation
import pandas as pd
import sys, math, os
import matplotlib.pyplot as plt
import datetime
from pytz import timezone
from astroplan import Observer, FixedTarget
from astropy.time import Time
import warnings
import configparser, argparse
from ast import literal_eval




# Define functions for scheduling.
def schedule_flux_cal(flux_cal_df, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule):
    """
    Schedule a flux calibrator
    """
    mask = []
    # Loop through all flux calibrators and check if they are visible at the current time (mask = True)

    for index, row in flux_cal_df.iterrows():
        rise_time = row['Rise_time_sidereal_20_deg'].time()
        # Ensure flux calibrator does not set during the next scan.
        set_time = (row['Set_time_sidereal_20_deg'] - datetime.timedelta(seconds=flux_calib_exposure_time)).time()
        day_end = datetime.time(23,59,59,999999)
        day_start = datetime.time(0,0,0,0)
        no_date_current_time_lst = current_time_lst.time()

        # If source is wrapping around midnight
        if rise_time > set_time:
            # Are we currently before midnight?
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= day_end):
                mask.append(True)
            # Are we after midnight?
            elif (no_date_current_time_lst >= day_start) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)

        # If source is not wrapping around midnight (simple case)

        else:
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)

    flux_cals_available = flux_cal_df[mask]

    if flux_cals_available.empty:
        warnings.warn(f"Warning: No flux calibrator available at {current_time_utc} (UTC), {current_time_lst} (LST)")
        return current_time_utc, current_time_lst, schedule
    
    else:
        # Select the first flux calibrator. Change this later to pick the one with the highest elevation
        flux_cal = flux_cals_available.head(1)
        # Add the flux calibrator to the schedule
        schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), flux_cal['Pointing'].iloc[0], flux_cal['RA'].iloc[0], \
        flux_cal['DEC'].iloc[0], flux_cal['Obs_type'].iloc[0], flux_calib_exposure_time])

         # Update the current time
        current_time_utc = current_time_utc + datetime.timedelta(seconds=flux_calib_exposure_time)
        current_time_lst = current_time_lst + datetime.timedelta(seconds=flux_calib_exposure_time)

        return current_time_utc, current_time_lst, schedule

def schedule_pol_cal(pol_cal_df, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule):
    """
    Schedule a polarisation calibrator
    """
    mask = []
    # Loop through all polarisation calibrators and check if they are visible at the current time (mask = True)

    for index, row in pol_cal_df.iterrows():
        rise_time = row['Rise_time_sidereal_20_deg'].time()
        # Ensure pol calibrator does not set during the next scan.
        set_time = (row['Set_time_sidereal_20_deg'] - datetime.timedelta(seconds=pol_calib_exposure_time)).time()
        day_end = datetime.time(23,59,59,999999)
        day_start = datetime.time(0,0,0,0)
        no_date_current_time_lst = current_time_lst.time()

        # If source is wrapping around midnight

        if rise_time > set_time:
            # Are we currently before midnight?
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= day_end):
                mask.append(True)
            # Are we after midnight?
            elif (no_date_current_time_lst >= day_start) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)

        # If source is not wrapping around midnight (simple case)

        else:
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)
    

    pol_cals_available = pol_cal_df[mask]

    if pol_cals_available.empty:
        warnings.warn(f"Warning: No polarisation calibrator available at {current_time_utc} (UTC), {current_time_lst} (LST)")
        return current_time_utc, current_time_lst, schedule, pol_cal_df['Set_time_sidereal_20_deg'].iloc[0]
    
    else:
        # Select the first polarisation calibrator. Change this later to pick the one with the highest elevation
        pol_cal = pol_cals_available.head(1)
        # Add the polarisation calibrator to the schedule
        schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), pol_cal['Pointing'].iloc[0], pol_cal['RA'].iloc[0], \
        pol_cal['DEC'].iloc[0], pol_cal['Obs_type'].iloc[0], pol_calib_exposure_time])

        # Update the current time
        current_time_utc = current_time_utc + datetime.timedelta(seconds=pol_calib_exposure_time)
        current_time_lst = current_time_lst + datetime.timedelta(seconds=pol_calib_exposure_time)

        return current_time_utc, current_time_lst, schedule, pol_cal['Set_time_sidereal_20_deg'].iloc[0]

def schedule_phase_cal(phase_cal_df, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule, nearest=False, nearest_target=''):
    """
    Schedule a phase calibrator. If nearest=True, the nearest phase calibrator will be scheduled. 
    If nearest=False, the first phase calibrator will be scheduled.
    nearest_target is the last target scheduled. This is used to determine the nearest phase calibrator.
    Should look like a list ['MSGPS_0001', RA, DEC]
    """
    mask = []
    # Loop through all phase calibrators and check if they are visible at the current time (mask = True)

    for index, row in phase_cal_df.iterrows():
        rise_time = row['Rise_time_sidereal_20_deg'].time()
        # Ensure that the source does not set during the phase cal observation
        set_time = (row['Set_time_sidereal_20_deg'] - datetime.timedelta(seconds=phase_calib_exposure_time)).time()
        day_end = datetime.time(23,59,59,999999)
        day_start = datetime.time(0,0,0,0)
        no_date_current_time_lst = current_time_lst.time()
        # If source is wrapping around midnight

        if rise_time > set_time:
            # Are we currently before midnight?
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= day_end):
                mask.append(True)
            # Are we after midnight?
            elif (no_date_current_time_lst >= day_start) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)
    

        # If source is not wrapping around midnight (simple case)

        else:
            if (no_date_current_time_lst >= rise_time) & (no_date_current_time_lst <= set_time):
                mask.append(True)
            else:
                mask.append(False)

    phase_cals_available = phase_cal_df[mask]


    if phase_cals_available.empty:
        warnings.warn(f"Warning: No phase calibrator available at {current_time_utc} (UTC), {current_time_lst} (LST)")

        return current_time_utc, current_time_lst, schedule

    else:

        if nearest:
            phase_cal_available_ra = phase_cals_available['RA'].to_numpy()
            phase_cal_available_dec = phase_cals_available['DEC'].to_numpy()
            phase_cal_available_pointings = phase_cals_available['Pointing'].to_numpy()
            phase_cal_catalog = SkyCoord(ra = phase_cal_available_ra, dec = phase_cal_available_dec, unit=(u.hourangle, u.deg))
            nearest_target = SkyCoord(ra=nearest_target[1], dec=nearest_target[2], unit=(u.hourangle, u.deg))
            nearest_phase_cal = phase_cal_catalog.separation(nearest_target).argmin()
            # Add the phase calibrator to the schedule
            schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), phase_cal_available_pointings[nearest_phase_cal], phase_cal_available_ra[nearest_phase_cal], \
                phase_cal_available_dec[nearest_phase_cal], 'Phase_Cal', phase_calib_exposure_time])

        else:
            # Select the first phase calibrator. Change this later to pick the one with the highest elevation
            phase_cal = phase_cals_available.head(1)
            # Add the phase calibrator to the schedule
            schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], \
            phase_cal['DEC'].iloc[0], phase_cal['Obs_type'].iloc[0], phase_calib_exposure_time])

        # Update the current time
        current_time_utc = current_time_utc + datetime.timedelta(seconds=phase_calib_exposure_time)
        current_time_lst = current_time_lst + datetime.timedelta(seconds=phase_calib_exposure_time)

        return current_time_utc, current_time_lst, schedule


def schedule_target(target_df, target_exposure_time, current_time_utc, current_time_lst, schedule, targets_observed):
    """
    Schedule a target
    targets observed is a list with the targets that have already been observed
    """
    # Ensure that source won't set before the end of a single target scan. If it does, skip it.

    target_df['Rise_time_sidereal_50_deg_minus_target_scan'] = target_df['Rise_time_sidereal_50_deg']
    target_df['Set_time_sidereal_20_deg_minus_target_scan'] = target_df['Set_time_sidereal_20_deg']

    target_df['Rise_time_sidereal_50_deg_minus_target_scan'] = target_df['Rise_time_sidereal_50_deg'] - datetime.timedelta(seconds=target_exposure_time)
    target_df['Set_time_sidereal_20_deg_minus_target_scan'] = target_df['Set_time_sidereal_20_deg'] - datetime.timedelta(seconds=target_exposure_time)


    #source rising from 20-50 degrees or setting down from 50-20 degrees.
    # Targets should be in the 20-50 degree range. If not, skip them.    
    targets_available = target_df.loc[(((target_df['Rise_time_sidereal_20_deg'] < current_time_lst) \
        & (target_df['Rise_time_sidereal_50_deg_minus_target_scan'] > current_time_lst)) | \
            ((target_df['Set_time_sidereal_20_deg_minus_target_scan'] > current_time_lst) \
        & (target_df['Set_time_sidereal_50_deg'] < current_time_lst)))]
    
    # Convert the targets observed list to a list of names
    targets_observed_names = [x[0] for x in targets_observed]
    
    # Remove the targets that have already been observed
    targets_available = targets_available.loc[~targets_available['Pointing'].isin(targets_observed_names)]

    if targets_available.empty:
        raise ValueError(f"Error: No target available at {current_time_utc} (UTC), {current_time_lst} (LST). \
            Change your observing start time or length of this session.")

    else:        
        # Select the first target since they are sorted by RA.
        target = targets_available.head(1)
        # Add the target to the schedule
        schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), target['Pointing'].iloc[0], target['RA'].iloc[0], \
        target['DEC'].iloc[0], target['Obs_type'].iloc[0], target_exposure_time])
        # Add the target to the list of targets observed
        targets_observed.append([target['Pointing'].iloc[0], target['RA'].iloc[0], target['DEC'].iloc[0], current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time()])

        # Update the current time
        current_time_utc = current_time_utc + datetime.timedelta(seconds=target_exposure_time)
        current_time_lst = current_time_lst + datetime.timedelta(seconds=target_exposure_time)
    
        return current_time_utc, current_time_lst, schedule, targets_observed


def schedule_slews(current_time_utc, current_time_lst, slew_time, schedule):
    """
    Schedule slews
    """

    # Update the current time
    current_time_utc = current_time_utc + datetime.timedelta(seconds=slew_time)
    current_time_lst = current_time_lst + datetime.timedelta(seconds=slew_time)

    # Add a slew to the schedule
    schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), 'Slew', 'Slew', 'Slew', 'Slew', slew_time])

    return current_time_utc, current_time_lst, schedule

def run_scheduler(start_time_utc, start_time_lst, end_time_utc, end_time_lst, flux_cals, pol_cals, \
    targets, phase_cals, flux_calib_exposure_time, pol_calib_exposure_time, target_exposure_time, phase_calib_exposure_time, slew_time, target_phase_calib_cycle_length):

    """
    Run the scheduler
    """
    #Build the scheduler
    schedule = []
    targets_observed = []
    # Be careful here. UTC times are given by user hence will have different dates. LST times will always have datetime.today() as the date
    current_time_utc = start_time_utc
    current_time_lst = start_time_lst
    repeat_polarisation_calibrators = True
    # Schedule the flux calibrators
    #while current_time_utc < end_time_utc:

    # Schedule the flux calibrators
    current_time_utc, current_time_lst, schedule = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    # Schedule the polarisation calibrators
    current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    
    # Schedule the phase calibrators
    current_time_utc, current_time_lst, schedule = schedule_phase_cal(phase_cals, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule)
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)


    time_left_for_observation = end_time_utc - current_time_utc
    
    '''
    To figure out how many pairs of phase calibrators and targets can be observed, we need to 
    know how much time is left for the observation and subtract the time for the flux calibrator,
        polarisation calibrator. '''
    
    time_remaining_for_targets_phase_calibrators = time_left_for_observation - datetime.timedelta(seconds=pol_calib_exposure_time) - datetime.timedelta(seconds=flux_calib_exposure_time) - datetime.timedelta(seconds = 2 * slew_time)
    one_group_target_phase_calib_slew_block = datetime.timedelta(seconds = target_exposure_time) * target_phase_calib_cycle_length  + datetime.timedelta(seconds=phase_calib_exposure_time) + datetime.timedelta(seconds = 2 * slew_time)
    number_of_pairs_target_phase_calib_slew_block = time_remaining_for_targets_phase_calibrators // one_group_target_phase_calib_slew_block
    print(time_remaining_for_targets_phase_calibrators, one_group_target_phase_calib_slew_block, number_of_pairs_target_phase_calib_slew_block)
    for i in range(number_of_pairs_target_phase_calib_slew_block):
        if (pol_cal_lst_set_time - current_time_lst < datetime.timedelta(seconds= 6 * pol_calib_exposure_time) and repeat_polarisation_calibrators == True):
            
            # Observe polarisation calibrators again at a different parallactic angle in its last 30-mins in the sky.
            current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_polarisation_calibrators = False

        ''' Schedule targets based on PHASE_CALIBRATOR_CYCLE_LENGTH. If = 1, then we observe target 
        and phase calibrator alternatively. If = 2, then we observe two targets followed by the 
        phase calibrator. '''

        for i in range(target_phase_calib_cycle_length):
            # Schedule the targets. Slew between targets is assumed to be negligible.
            current_time_utc, current_time_lst, schedule, targets_observed = schedule_target(targets, target_exposure_time, current_time_utc, current_time_lst, schedule, targets_observed)

        current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    
        # Schedule the nearest phase calibrators
        current_time_utc, current_time_lst, schedule = schedule_phase_cal(phase_cals, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule, nearest=True, nearest_target = targets_observed[-1])
        current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    #current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    # If less than 6 minutes left, then observe the flux calibrator again.
    #elif (remaining_time < datetime.timedelta(seconds= 1.2 * flux_calib_exposure_time)):

    while current_time_utc < end_time_utc:
        remaining_time = end_time_utc - current_time_utc
        # In case the polarisation calibrator has not been scheduled.
        if repeat_polarisation_calibrators == True:
            current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_polarisation_calibrators = False
        # If there is enough time left, schedule the target/phase calibrator again
        if remaining_time > datetime.timedelta(seconds=target_exposure_time + flux_calib_exposure_time):
            current_time_utc, current_time_lst, schedule, targets_observed = schedule_target(targets, target_exposure_time, current_time_utc, current_time_lst, schedule, targets_observed)
        
        if remaining_time > datetime.timedelta(seconds=phase_calib_exposure_time + flux_calib_exposure_time):
            current_time_utc, current_time_lst, schedule = schedule_phase_cal(phase_cals, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule)
        
        current_time_utc, current_time_lst, schedule = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)

        

    return schedule, targets_observed

def convert_meerkat_lst_to_utc(lst, date, meerkat):

    '''
    Neither astroplan or ephem has a function to convert LST to UTC. Therefore, this is a hack that converts LST to UTC 
    using basic cross multiplication and addition/subtraction of time difference between sidderal day and solar day.
    '''

    # Convert LST to UTC
    year = int(date.split('-')[0])
    month = int(date.split('-')[1])
    day = int(date.split('-')[2])
    hrs = int(lst.split(':')[0])
    mins = int(lst.split(':')[1])
    secs = int(lst.split(':')[2])
    lst = datetime.time(hrs, mins, secs)
    date = datetime.date(year, month, day)
    sidderal_day = datetime.datetime(year, month, day, 23, 56, 4, 90500)
    solar_day = datetime.datetime(year, month, day, 23, 59, 59, 999999)
    factor_diff_solar_sidderal = (solar_day - sidderal_day).total_seconds() + 0.00001
   

    required_lst = datetime.datetime.combine(date, lst)
    temp_lst = meerkat.local_sidereal_time(required_lst).to_string(sep=':')[:-3]
    seconds = int(math.modf(float(temp_lst.split(':')[2]))[1])
    microseconds = int(math.modf(float(temp_lst.split(':')[2]))[0] * 1000000)

    temp_lst = datetime.datetime(year, month, day, int(temp_lst.split(':')[0]), int(temp_lst.split(':')[1]), seconds, microseconds)
   
    if temp_lst.time() > required_lst.time():
        time_diff = (temp_lst - required_lst).total_seconds()
        
        utc_start = required_lst - datetime.timedelta(seconds=time_diff)
        # conversion factor because utc and sidderal time change at a different rate.
        additional_time = factor_diff_solar_sidderal * time_diff/86400
        utc_start = utc_start + datetime.timedelta(seconds=additional_time)
        checking_lst = meerkat.local_sidereal_time(utc_start).to_string(sep=':')
        
        # If the required date is not the same as the utc_start date, then add the sidderal day to utc_start
        if required_lst.date() > utc_start.date():
            date_to_be_added = required_lst.date() - utc_start.date()
            utc_start = utc_start + datetime.timedelta(hours = sidderal_day.hour, minutes = sidderal_day.minute, seconds = sidderal_day.second, microseconds = sidderal_day.microsecond)
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start
        elif required_lst.date() < utc_start.date():
            date_to_be_added = utc_start.date() - required_lst.date()
            utc_start = utc_start - datetime.timedelta(hours = sidderal_day.hour, minutes = sidderal_day.minute, seconds = sidderal_day.second, microseconds = sidderal_day.microsecond)
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start
        else:
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start

    else:
        time_diff = (required_lst - temp_lst).total_seconds()
        utc_start = required_lst + datetime.timedelta(seconds=time_diff)
       
        # conversion factor because utc and sidderal time change at a different rate.
        additional_time = factor_diff_solar_sidderal * time_diff/86400
        utc_start = utc_start - datetime.timedelta(seconds=additional_time)
        checking_lst = meerkat.local_sidereal_time(utc_start).to_string(sep=':')
        
        # If the required date is not the same as the utc_start date, then add the sidderal day to utc_start
        if required_lst.date() > utc_start.date():
            date_to_be_added = required_lst.date() - utc_start.date()
            utc_start = utc_start + datetime.timedelta(hours = sidderal_day.hour, minutes = sidderal_day.minute, seconds = sidderal_day.second, microseconds = sidderal_day.microsecond)
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start
        elif required_lst.date() < utc_start.date():
            date_to_be_added = utc_start.date() - required_lst.date()
            utc_start = utc_start - datetime.timedelta(hours = sidderal_day.hour, minutes = sidderal_day.minute, seconds = sidderal_day.second, microseconds = sidderal_day.microsecond)
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start
        else:
            print('With a requested LST time of ', lst, ' and  UTC start date of ', utc_start.date(), 'the UTC observation start time is ', utc_start.time())
            return utc_start
                    

    

def main(config, output_file):

    sections = config.sections()

    #Define exposure times of targets and calibrators in seconds
    target_exposure_time = literal_eval(config['INTEGRATION_TIME']['TARGET'])
    pol_calib_exposure_time = literal_eval(config['INTEGRATION_TIME']['POLARISATION_CALIBRATOR']) 
    flux_calib_exposure_time = literal_eval(config['INTEGRATION_TIME']['FLUX_CALIBRATOR']) 
    phase_calib_exposure_time = literal_eval(config['INTEGRATION_TIME']['PHASE_CALIBRATOR']) 
    slew_time = literal_eval(config['INTEGRATION_TIME']['SLEW_TIME']) # seconds

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

    # MeerKAT Location
    location = EarthLocation.from_geodetic(mkt_location['longitude'] * u.deg, \
        mkt_location['latitude'] * u.deg, mkt_location['elevation'] * u.m)


    meerkat = Observer(name='MeerKAT',
                location=location,
                timezone=timezone('Africa/Johannesburg'),
                description="MeerKAT telescope in SA")

    # Read the data
    targets = pd.read_csv('rise_set_time_target.csv')
    phase_cals = pd.read_csv('rise_set_time_phase_cal.csv')
    flux_cals = pd.read_csv('rise_set_time_flux_cal.csv')
    pol_cals = pd.read_csv('rise_set_time_pol_cal.csv')
    output_targets_observed ='targets_already_scheduled.csv'
    targets_already_scheduled_list = []
    if os.path.isfile(output_targets_observed):
        targets_already_scheduled = pd.read_csv(output_targets_observed)
        targets_already_scheduled_list = targets_already_scheduled['Pointing'].to_list()
    # If targets have been previously scheduled, then remove them from the list of targets to be scheduled
    if targets_already_scheduled_list:
        targets = targets[~targets['Pointing'].isin(targets_already_scheduled_list)]
        print('Removed %d targets from previous runs that have already been scheduled, we have %d targets left' % (len(targets_already_scheduled_list), len(targets)))
    else:
        print('No targets from previous runs have been scheduled yet. Will use the full list of targets')

    # Convert the rise and set times to datetime64
    targets['Rise_time_sidereal_20_deg'] = targets['Rise_time_sidereal_20_deg'].astype('datetime64[ns]')
    targets['Set_time_sidereal_20_deg'] = targets['Set_time_sidereal_20_deg'].astype('datetime64[ns]')
    targets['Rise_time_sidereal_50_deg'] = targets['Rise_time_sidereal_50_deg'].astype('datetime64[ns]')
    targets['Set_time_sidereal_50_deg'] = targets['Set_time_sidereal_50_deg'].astype('datetime64[ns]')

    flux_cals['Rise_time_sidereal_20_deg'] = flux_cals['Rise_time_sidereal_20_deg'].astype('datetime64[ns]')
    flux_cals['Set_time_sidereal_20_deg'] = flux_cals['Set_time_sidereal_20_deg'].astype('datetime64[ns]')

    pol_cals['Rise_time_sidereal_20_deg'] = pol_cals['Rise_time_sidereal_20_deg'].astype('datetime64[ns]')
    pol_cals['Set_time_sidereal_20_deg'] = pol_cals['Set_time_sidereal_20_deg'].astype('datetime64[ns]')


    # Phase cals have a special case of naT. Source never rises or sets.

    phase_cals['Rise_time_sidereal_20_deg'] = pd.to_datetime(phase_cals['Rise_time_sidereal_20_deg'], errors='coerce')
    phase_cals['Set_time_sidereal_20_deg'] = pd.to_datetime(phase_cals['Set_time_sidereal_20_deg'], errors='coerce')


    session_length = literal_eval(config['OBSERVATION_PLAN']['TOTAL_SESSION_TIME_HRS'])
    requested_lst_time = literal_eval(config['OBSERVATION_PLAN']['LST_START'])
    obs_date = literal_eval(config['OBSERVATION_PLAN']['OBSERVATION_DATE'])
    start_time_utc = convert_meerkat_lst_to_utc(requested_lst_time, obs_date, meerkat)
    target_phase_calib_cycle_length = literal_eval(config['OBSERVATION_PLAN']['PHASE_CALIBRATOR_CYCLE_LENGTH'])
    # Old method of directly giving utc time
    #start_time_utc = '2023-02-01 19:00:00'  
    #start_time_utc = datetime.datetime.strptime(start_time_utc, '%Y-%m-%d %H:%M:%S')

    # Calculate the end time of the session
    end_time_utc = datetime.timedelta(hours = session_length) + start_time_utc

    # Convert utc to local sidereal time

    start_time_lst = meerkat.local_sidereal_time(start_time_utc).to_string(sep=':')[:-3]
    end_time_lst = meerkat.local_sidereal_time(end_time_utc).to_string(sep=':')[:-3]
    
    # Add today's date to the time. We need this later for scheduling. astype('datetime64[ns]') automatically adds today's date
    start_time_lst = datetime.datetime.strptime(str(datetime.date.today()) + ' ' + start_time_lst, '%Y-%m-%d %H:%M:%S.%f')
    end_time_lst = datetime.datetime.strptime(str(datetime.date.today()) + ' ' + end_time_lst, '%Y-%m-%d %H:%M:%S.%f')


    # If Rise time > Set time, then add 1 day to the Set time. 
    # This is because the target is visible after midnight, hence you wrap around.
    flux_cals['Set_time_sidereal_20_deg'] = np.where(flux_cals['Rise_time_sidereal_20_deg'] > flux_cals['Set_time_sidereal_20_deg'], \
        flux_cals['Set_time_sidereal_20_deg'] + datetime.timedelta(days=1), flux_cals['Set_time_sidereal_20_deg']) 

    pol_cals['Set_time_sidereal_20_deg'] = np.where(pol_cals['Rise_time_sidereal_20_deg'] > pol_cals['Set_time_sidereal_20_deg'], \
        pol_cals['Set_time_sidereal_20_deg'] + datetime.timedelta(days=1), pol_cals['Set_time_sidereal_20_deg'])

    phase_cals['Set_time_sidereal_20_deg'] = np.where(phase_cals['Rise_time_sidereal_20_deg'] > phase_cals['Set_time_sidereal_20_deg'], \
        phase_cals['Set_time_sidereal_20_deg'] + datetime.timedelta(days=1), phase_cals['Set_time_sidereal_20_deg'])

    targets['Set_time_sidereal_20_deg'] = np.where(targets['Rise_time_sidereal_20_deg'] > targets['Set_time_sidereal_20_deg'], \
        targets['Set_time_sidereal_20_deg'] + datetime.timedelta(days=1), targets['Set_time_sidereal_20_deg'])

    targets['Set_time_sidereal_50_deg'] = np.where(targets['Rise_time_sidereal_50_deg'] > targets['Set_time_sidereal_50_deg'], \
        targets['Set_time_sidereal_50_deg'] + datetime.timedelta(days=1), targets['Set_time_sidereal_50_deg'])

    # Sort targets by Right Ascension
    targets = targets.sort_values(by=['RA'])
    # Sort flux calibrators by Right Ascension
    flux_cals = flux_cals.sort_values(by=['RA'])
    # Sort polarisation calibrators by Right Ascension
    pol_cals = pol_cals.sort_values(by=['RA'])
    # Sort phase calibrators by Right Ascension
    phase_cals = phase_cals.sort_values(by=['RA'])

    schedule, targets_observed = run_scheduler(start_time_utc, start_time_lst, end_time_utc, end_time_lst, flux_cals, pol_cals, \
        targets, phase_cals, flux_calib_exposure_time, pol_calib_exposure_time, target_exposure_time, phase_calib_exposure_time, slew_time, target_phase_calib_cycle_length)

    final_schedule = pd.DataFrame(schedule, columns=['Start Time UTC', 'Start Time LST', 'Source_Name', 'RA', 'DEC', 'Obs_Type', 'Integration_Time'])
    print_slews = False
    if print_slews:
        final_schedule.to_csv(output_file + '.csv', index=False)
    else:
        final_schedule = final_schedule[final_schedule['Obs_Type'] != 'Slew']
        final_schedule.to_csv(output_file + '.csv', index=False)
    targets_only = final_schedule[final_schedule['Obs_Type'] == 'Target']
    print('Total number of targets observed this session: ', len(targets_only))
    print('Efficiency of this observing session: ', round(len(targets_only)*target_exposure_time * 100/(session_length * 3600), 2), '%')
    targets_observed = pd.DataFrame(targets_observed, columns=['Pointing', 'RA', 'DEC', 'Obs_time_UTC', 'Obs_Time_LST'])
    
    targets_observed.to_csv(output_targets_observed, mode='a', header=not os.path.exists(output_targets_observed), index=False)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Schedule observations for the MMGPS S-BAND telescope')
    parser.add_argument('-c', '--config_file', help='Configuration file for the observation', default='sband_schedule.cfg')
    parser.add_argument('-o', '--output_file', help='Output file for the observation schedule', default='sband_schedule')
    args = parser.parse_args()
    config_file = args.config_file
    output_file = args.output_file
    config = configparser.ConfigParser()
    config.read(config_file)
    main(config, output_file)



