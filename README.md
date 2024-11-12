# Scheduler Software for the MPIfR MeerKAT UHF Survey

The MMGPS-S scheduler software adapted for the UHF survey.

1. Observations start with a noise diode firing, folled by a flux calibrator and a polarisation calibrator.
2. This is followed by target scans, sandwiched by a phase calibrator visit. How many targets observed before a phase calibrator is scheduled is a tunable parameter (PHASE_CALIBRATOR_CYCLE_LENGTH) in the configuration file.
3. We end the observation with another visit of the polarisation and flux calibrator.
4. Targets are scheduled with the constraint that their elevation at MeerKAT needs to be between 15 - 50 deg.
5. Calibrators are scheduled with the constraint that their elevation at MeerKAT needs to be between 15 - 90 deg.
6. For the UHF observations, gain calibrators are tagged as "polcal" to prevent APSUSE updating the gain solutions from unreliable gain targets

# How to use this software?

1. Edit the config file 'uhf_schedule.cfg'.
2. Update the LST time and UTC date when you want to start observing plus the total session length in hours.
3. Then run the code uhf_scheduler.py.

Dependencies:
Pandas, astropy, pytz, numpy
