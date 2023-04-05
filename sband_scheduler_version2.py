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
import json



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

        return current_time_utc, current_time_lst, schedule, flux_cal['Set_time_sidereal_20_deg'].iloc[0]

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


def schedule_phase_cal_schedule_block_start(phase_cal_df, phase_calib_exposure_time, target_df, target_exposure_time, targets_observed, current_time_utc, current_time_lst, schedule, cycle_length_duration, last_phase_cal=''):
    """
    Schedule a phase calibrator that has the maximum hits in the available target list 
    and is available for the entire cycle.
    """
    
    mask = []
    # Loop through all phase calibrators and check if they are visible at the current time (mask = True)

    for index, row in phase_cal_df.iterrows():
        rise_time = row['Rise_time_sidereal_20_deg'].time()
        # Ensure that the source does not set during the the cycle length
        set_time = (row['Set_time_sidereal_20_deg'] - datetime.timedelta(seconds=cycle_length_duration)).time()
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
    
    # Start of a new schedule block. Schedule a phase cal based on most hits in the available target list
    # First get the dataframe of targets that are visible at the current time
    phase_cal_max_hits_selected = ''
    targets_available = get_available_target_now(target_df, target_exposure_time, current_time_utc, current_time_lst, targets_observed)
    phase_cal_max_hits = targets_available['Nearest_Phase_Cal'].value_counts().index.tolist()
    for i in range(len(phase_cal_max_hits)):
        if phase_cal_max_hits[i] in phase_cals_available['Pointing'].tolist():
            phase_cal_max_hits_selected = phase_cal_max_hits[i]
            break
    
    if phase_cal_max_hits_selected == last_phase_cal:
        # Don't schedule phase calibrator if it is the same as the last one
        # Will return the last phase calibrator in the schedule
        phase_cal = phase_cals_available[phase_cals_available['Pointing'] == last_phase_cal]
        return current_time_utc, current_time_lst, schedule, phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], phase_cal['DEC'].iloc[0]

    # If no phase calibrator has target hits, then schedule J1619-8418 because it's always visible.
    elif phase_cal_max_hits_selected == '':
        phase_cal_max_hits_selected = 'J1619-8418'

    phase_cal = phase_cals_available[phase_cals_available['Pointing'] == phase_cal_max_hits_selected]
    # Add the phase calibrator to the schedule
    schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], \
    phase_cal['DEC'].iloc[0], phase_cal['Obs_type'].iloc[0], phase_calib_exposure_time])
    # Update the current time
    current_time_utc = current_time_utc + datetime.timedelta(seconds=phase_calib_exposure_time)
    current_time_lst = current_time_lst + datetime.timedelta(seconds=phase_calib_exposure_time)

    return current_time_utc, current_time_lst, schedule, phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], phase_cal['DEC'].iloc[0]

def schedule_phase_cal_schedule_block_end(phase_cal_df, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule, last_phase_cal=''):
    """
    Schedule same phase calibrator to end the schedule block. We already checked for rise and set times.
    """

    phase_cal = phase_cal_df[phase_cal_df['Pointing'] == last_phase_cal]
    # Add the phase calibrator to the schedule
    schedule.append([current_time_utc.strftime(format='%Y-%m-%d %H:%M:%S'), current_time_lst.time(), phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], \
    phase_cal['DEC'].iloc[0], phase_cal['Obs_type'].iloc[0], phase_calib_exposure_time])
    # Update the current time
    current_time_utc = current_time_utc + datetime.timedelta(seconds=phase_calib_exposure_time)
    current_time_lst = current_time_lst + datetime.timedelta(seconds=phase_calib_exposure_time)

    return current_time_utc, current_time_lst, schedule, phase_cal['Pointing'].iloc[0], phase_cal['RA'].iloc[0], phase_cal['DEC'].iloc[0]


def get_available_target_now(target_df, target_exposure_time, current_time_utc, current_time_lst, targets_observed):
    """
    Returns a dataframe of available targets at the current time.
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
        return targets_available


def schedule_target(target_df, target_exposure_time, current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_RA, last_phase_cal_DEC):
    
    """
    Schedule a target
    targets observed is a list with the targets that have already been observed. 
    last_phase_cal is the last phase calibrator scheduled. This is used to determine the nearest target.
    Should look like a list ['J0825-5010', RA, DEC]
    """
    # Ensure that source won't set before the end of a single target scan. If it does, skip it.
    targets_available = get_available_target_now(target_df, target_exposure_time, current_time_utc, current_time_lst, targets_observed)
    
    # Find the nearest target to the last phase calibrator
    targets_nearest_last_phase_calib = targets_available.loc[targets_available['Nearest_Phase_Cal'] == last_phase_cal]      
    if targets_nearest_last_phase_calib.empty:
        print('#######################################################################')
        print(f"The last phase calibrator scheduled was {last_phase_cal}. This is not the nearest phase cal to any available targets. I will re-do the distance calculation to find the next nearest target. We maybe observing the second or third nearest phase calibrator now but we will have the same phase calibrator between any targets")
        print('#######################################################################')
        target_ra = targets_available['RA'].to_numpy()
        target_dec = targets_available['DEC'].to_numpy()
        target_pointings = targets_available['Pointing'].to_numpy()
        target_catalog = SkyCoord(ra = target_ra, dec = target_dec, unit=(u.hourangle, u.deg))
        last_phase_cal_coords = SkyCoord(ra=last_phase_cal_RA, dec=last_phase_cal_DEC, unit=(u.hourangle, u.deg))
        nearest_target = target_catalog.separation(last_phase_cal_coords).argmin()
        nearest_target = target_pointings[nearest_target]
        target = targets_available.loc[targets_available['Pointing'] == nearest_target]
        print(f"The nearest target to the phase calibrator {last_phase_cal} is {nearest_target}")
    else:

        target = targets_available.sort_values('Pointing').head(1)
        #print(f"The nearest target to the phase calibrator {last_phase_cal} is {target['Pointing'].iloc[0]}")

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


def schedule_block_phase_calibrator_and_target(phase_cals, targets, phase_calib_exposure_time, target_exposure_time, slew_time, current_time_utc, current_time_lst, schedule, targets_observed, target_phase_calib_cycle_length, cycle_length_duration, last_phase_cal=''):
    """
    Schedule a block of phase calibrator and target based on the target_phase_calib_cycle_length. 
    Start and stop with the same phase calibrator.
    """

    # Schedule the first phase calibrator only if it is not already scheduled and is available for the entire cycle.
    current_time_utc, current_time_lst, schedule, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec = schedule_phase_cal_schedule_block_start(phase_cals, phase_calib_exposure_time, targets, target_exposure_time, targets_observed, current_time_utc, current_time_lst, schedule, cycle_length_duration, last_phase_cal = last_phase_cal)
    # Slew time
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    # Schedule target(s) based on the cycle length
    for i in range(target_phase_calib_cycle_length):
        current_time_utc, current_time_lst, schedule, targets_observed = schedule_target(targets, target_exposure_time, current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec)
        # Slew between targets is assumed to be negligible
    
    # Slew time
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)

    # Close the schedule with the same phase calibrator
    current_time_utc, current_time_lst, schedule, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec = schedule_phase_cal_schedule_block_end(phase_cals, phase_calib_exposure_time, current_time_utc, current_time_lst, schedule, last_phase_cal = last_phase_cal)
    return current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec


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
    repeat_flux_calibrators = True
    # Schedule the flux calibrators
    #while current_time_utc < end_time_utc:

    # Schedule the flux calibrators
    current_time_utc, current_time_lst, schedule, flux_cal_lst_set_time = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    # Schedule the polarisation calibrators
    current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
    current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
    
    '''
    First calculate the time left for observation assuming one more flux and pol cal and the end
     and the number of schedule blocks that can be done '''

    time_left_for_observation = end_time_utc - current_time_utc
    # Subtract the time for the last flux and pol cal at the end
    time_remaining_for_targets_phase_calibrators = time_left_for_observation - datetime.timedelta(seconds=pol_calib_exposure_time) - datetime.timedelta(seconds=flux_calib_exposure_time) - datetime.timedelta(seconds = 2 * slew_time)
    # Calculate the number of schedule blocks that can be done
    one_group_target_phase_calib_slew_block = datetime.timedelta(seconds = target_exposure_time) * target_phase_calib_cycle_length  + datetime.timedelta(seconds= phase_calib_exposure_time) * 2 + datetime.timedelta(seconds = slew_time) * 2 # Slew from phase cal to target and target to phase cal
   
    number_of_pairs_target_phase_calib_slew_block = time_remaining_for_targets_phase_calibrators // one_group_target_phase_calib_slew_block
    ''' Start the schedule block with phase calibrator, target(s), phase calibrator
         Schedule targets based on PHASE_CALIBRATOR_CYCLE_LENGTH. If = 1, then we observe target 
         and phase calibrator alternatively. If = 2, then we observe two targets followed by the 
         same phase calibrator. '''
    for i in range(number_of_pairs_target_phase_calib_slew_block):

        ''' Depending on the start time of the observation, the polarisation calibrators may set during 
        our target phase calib block or later. If its set during the target phase calib block, we observe   
        it again at a different parallactic angle in its last 30-mins in the sky.
        '''
        
        if i == 0: # No previous phase cal for first block
            current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec  = schedule_block_phase_calibrator_and_target(phase_cals, targets, phase_calib_exposure_time, target_exposure_time, slew_time, current_time_utc, current_time_lst, schedule, targets_observed, target_phase_calib_cycle_length, one_group_target_phase_calib_slew_block.total_seconds())
        

        else:
            current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec = schedule_block_phase_calibrator_and_target(phase_cals, targets, phase_calib_exposure_time, target_exposure_time, slew_time, current_time_utc, current_time_lst, schedule, targets_observed, target_phase_calib_cycle_length, one_group_target_phase_calib_slew_block.total_seconds(), last_phase_cal = last_phase_cal)

        # Schedule the polarisation calibrators if they are in the last 30-mins on the sky
        remaining_pol_cal_time = pol_cal_lst_set_time - current_time_lst
        # taking care of the date wrapping problem
        if remaining_pol_cal_time > datetime.timedelta(days = 1):
            remaining_pol_cal_time = remaining_pol_cal_time - datetime.timedelta(days = int(remaining_pol_cal_time.days))

        if (remaining_pol_cal_time < datetime.timedelta(seconds= 6 * pol_calib_exposure_time) and repeat_polarisation_calibrators == True):
        
            current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_polarisation_calibrators = False  

        remaining_flux_cal_time = flux_cal_lst_set_time - current_time_lst
        # taking care of the date wrapping problem
        if remaining_flux_cal_time > datetime.timedelta(days = 1):
            remaining_flux_cal_time = remaining_flux_cal_time - datetime.timedelta(days = int(remaining_flux_cal_time.days))
                
        # Schedule the flux calibrators if they are in the last 30-mins on the sky
        if (remaining_flux_cal_time < datetime.timedelta(seconds= 6 * flux_calib_exposure_time) and repeat_flux_calibrators == True):
            current_time_utc, current_time_lst, schedule, flux_cal_lst_set_time = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_flux_calibrators = False                                              
    
    
    # The cycle change loop!
    remaining_time = end_time_utc - current_time_utc
    # Remaining time in case flux or pol cal hasn't been scheduled.
    if repeat_flux_calibrators == True:
        remaining_time = remaining_time - datetime.timedelta(seconds=flux_calib_exposure_time) - datetime.timedelta(seconds = slew_time)
    if repeat_polarisation_calibrators == True:
        remaining_time = remaining_time - datetime.timedelta(seconds=pol_calib_exposure_time) - datetime.timedelta(seconds = slew_time)
    # If there is enough time left, schedule the target/phase calibrator with a different cycle length
    while target_phase_calib_cycle_length!= 0:
        remaining_pol_cal_time = pol_cal_lst_set_time - current_time_lst
        # taking care of the date wrapping problem
        if remaining_pol_cal_time > datetime.timedelta(days = 1):
            remaining_pol_cal_time = remaining_pol_cal_time - datetime.timedelta(days = int(remaining_pol_cal_time.days))

        if (remaining_pol_cal_time < datetime.timedelta(seconds= 6 * pol_calib_exposure_time) and repeat_polarisation_calibrators == True):
        
            current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_polarisation_calibrators = False 

        remaining_flux_cal_time = flux_cal_lst_set_time - current_time_lst
        # taking care of the date wrapping problem

        if remaining_flux_cal_time > datetime.timedelta(days = 1):
            remaining_flux_cal_time = remaining_flux_cal_time - datetime.timedelta(days = int(remaining_flux_cal_time.days))
        if (remaining_flux_cal_time < datetime.timedelta(seconds= 6 * flux_calib_exposure_time) and repeat_flux_calibrators == True):
            current_time_utc, current_time_lst, schedule, flux_cal_lst_set_time = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)
            current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
            repeat_flux_calibrators = False            

        one_group_target_phase_calib_slew_block = datetime.timedelta(seconds = target_exposure_time) * target_phase_calib_cycle_length  + datetime.timedelta(seconds= phase_calib_exposure_time) * 2 + datetime.timedelta(seconds = slew_time) * 2 # Slew from phase cal to target and target to phase cal
        if remaining_time > one_group_target_phase_calib_slew_block:
            current_time_utc, current_time_lst, schedule, targets_observed, last_phase_cal, last_phase_cal_ra, last_phase_cal_dec = schedule_block_phase_calibrator_and_target(phase_cals, targets, phase_calib_exposure_time, target_exposure_time, slew_time, current_time_utc, current_time_lst, schedule, targets_observed, target_phase_calib_cycle_length, one_group_target_phase_calib_slew_block.total_seconds(), last_phase_cal = last_phase_cal)
            break
        else:
            target_phase_calib_cycle_length -= 1
    
    # Final scheduling of the polarisation and flux calibrators in case not done earlier.
    if repeat_polarisation_calibrators == True:
        current_time_utc, current_time_lst, schedule, pol_cal_lst_set_time = schedule_pol_cal(pol_cals, pol_calib_exposure_time, current_time_utc, current_time_lst, schedule)
        current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
        repeat_polarisation_calibrators = False
    if repeat_flux_calibrators == True:
        current_time_utc, current_time_lst, schedule, flux_cal_lst_set_time = schedule_flux_cal(flux_cals, flux_calib_exposure_time, current_time_utc, current_time_lst, schedule)
        current_time_utc, current_time_lst, schedule = schedule_slews(current_time_utc, current_time_lst, slew_time, schedule)
        repeat_flux_calibrators = False            


    return schedule, targets_observed, current_time_utc

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
                    

def generate_block(name, ra, dec, tags, ttype, duration):
    """
    Generate a python dict corresponding to an observing block 
    """
    target = {}
    target['name']     = name
    target['ra']       = ra.lstrip(' ').rstrip(' ')
    target['dec']      = dec.lstrip(' ').rstrip(' ')
    target['tags']     = tags
    target['type']     = ttype
    target['duration'] = duration
    #
    return(target)


    

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
    #targets = targets.sort_values(by=['RA'])
    # Sort flux calibrators by Right Ascension
    flux_cals = flux_cals.sort_values(by=['RA'])
    # Sort polarisation calibrators by Right Ascension
    pol_cals = pol_cals.sort_values(by=['RA'])
    # Sort phase calibrators by Right Ascension
    phase_cals = phase_cals.sort_values(by=['RA'])

    schedule, targets_observed, session_end_time = run_scheduler(start_time_utc, start_time_lst, end_time_utc, end_time_lst, flux_cals, pol_cals, \
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
    actual_session_length = (session_end_time - start_time_utc).total_seconds()/3600
    print('Actual session length: ', round(actual_session_length, 2), ' hours')
    print('Efficiency of this observing session: ', round(len(targets_only)*target_exposure_time * 100/(actual_session_length * 3600), 2), '%')
    targets_observed = pd.DataFrame(targets_observed, columns=['Pointing', 'RA', 'DEC', 'Obs_time_UTC', 'Obs_Time_LST'])
    
    targets_observed.to_csv(output_targets_observed, mode='a', header=not os.path.exists(output_targets_observed), index=False)
    print(final_schedule)
    #print(targets_observed)
    full_obs_schedule = []
    for index, row in final_schedule.iterrows():
        if row['Obs_Type'] == 'Flux_Cal':
            tags = ['delaycal', 'fluxcal', 'bpcal']
        elif row['Obs_Type'] == 'Pol_Cal':
            tags = ['polcal']
        elif row['Obs_Type'] == 'Phase_Cal':
            tags = ['gaincal']
        else:
            tags = ['target']
        
        single_block = generate_block(row['Source_Name'], row['RA'], row['DEC'], \
                                      tags, 'track', row['Integration_Time'])
        full_obs_schedule.append(single_block)
    
    # Create the JSON file
    with open('SCI-20200703-MK-02.json', 'r') as json_handle:
        template = json.load(json_handle)
    
    template['activities'] = []
    template['owner'] = 'Vishnu Balakrishnan'
    template['owner_email'] = 'vishnubk93@gmail.com'
    template['id'] = 0
    id_text = datetime.datetime.today().strftime('%Y%m%d') + '_3'

    template['description'] = 'MPIfR Galactic Plane Survey S-band: Setup {}'.format(id_text)
    template['proposal_id'] = 'SCI-20200703-MK-03'
    template['blocks'][0]['name'] = 'Setup {}'.format(id_text)
    template['blocks'][0]['targets'] = full_obs_schedule
    template['instrument']['integration_time']  = '8.0'
    template['instrument']['product']  = 'c875M4k'
    template['instrument']['center_freq']  = '2406.25'


    template['instrument']['pool_resources']    = ['apsuse', 'cbf', 'fbfuse', 'sdp', 'tuse']
    template['instrument']['config_auth_host'] = 'TRAPUM'
    template['horizon'] = 20.
   
    template['lst_start'] = start_time_lst.time().strftime('%H:%M')
    template['lst_start_end'] = end_time_lst.time().strftime('%H:%M')
    template['desired_start_time'] = str(start_time_utc)
    with open('SCI-20200703-MK-03_%s.json'%id_text, 'w') as json_handle:
        json.dump(template, json_handle, indent=4, sort_keys=False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Schedule observations for the MMGPS S-BAND survey')
    parser.add_argument('-c', '--config_file', help='Configuration file for the observation', default='sband_schedule.cfg')
    parser.add_argument('-o', '--output_file', help='Output file for the observation schedule', default='sband_schedule')
    args = parser.parse_args()
    config_file = args.config_file
    output_file = args.output_file
    config = configparser.ConfigParser()
    config.read(config_file)
    main(config, output_file)



