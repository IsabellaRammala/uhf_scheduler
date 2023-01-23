# Scheduler Software for the MPIfR MeerKAT S-BAND Survey

This is a simple scheduler software for the MPIfR MeerKAT S-BAND Survey.

1. Observations start with a flux calibrator, followed by a polarisation calibrator.
2. This is followed by a cycle of target and phase calibrator visits. How many targets observed before a phase calibrator is scheduled is a tunable parameter (PHASE_CALIBRATOR_CYCLE_LENGTH) in the configuration file.
3. We end the observation with another visit of the polarisation and flux calibrator.
4. Targets are scheduled with the constraint that their elevation at MeerKAT needs to be between 20 - 50 deg.
5. Calibrators are scheduled with the constraint that their elevation at MeerKAT needs to be between 20 - 90 deg.

# How to use this software?

1. Edit the config file 'sband_schedule.cfg'.
2. Update the LST time and UTC date when you want to start observing plus the total session length in hours.
3. Then run the code sband_scheduler.py.

Dependencies:
Pandas, astropy, pytz, numpy
