# Inertial-based algorithms for temporal parameters estimation in running

Written by Myriam Lubrano, Department of Electrical, Electronic and Information Engineering "Guglielmo Marconi", University of Bologna, Bologna, Italy.

This repository accompanies the paper:\
Lubrano M., Rossanigo R., Cereatti A., Cuppini C., Fantozzi S., "_In-field temporal phase analysis of running: performance assessment of 36 inertial-based algorithms across different device positions and running speeds_".
Submitted to Journal of Biomechanics, January 2026.\
The study evaluates the performance of IMU-based algorithms for stance and stride time estimation during running, across different sensor placements (foot, shank, and lower back) and running speeds (15, 20, and 25 km/h), using OptoJump as reference system.

## List of implemented algorithms
•	Aubol 2020\
•	Bailey 2015\
•	Blauberger 2021\
•	Chew 2018\
•	Donahue 2022 (foot)\
•	Donahue 2022 (lower back)\
•	Falbriard 2018\
•	Falbriard 2020\
•	Giandolini 2014\
•	Lee 2010\
•	Maiwald 2015\
•	Mitschke 2018\
•	Mo 2018\
•	Patoz 2022\
•	Reenalda 2019\
•	Schmidt 2016\
•	Strohrmann 2012\
•	Van Werkhoven 2019\
•	Yang 2022

## Purpose
Each function takes accerelometer, gyroscope and/or magnetometer data from the foot, shank or lower back, returning initial contacts (IC) and toe offs (TO) for stride and stance times estimation.

## Dependencies
•	MATLAB's Signal Processing Toolbox\
•	MATLAB's Aerospace Toolbox\
•	MATLAB's Sensor Fusion and Tracking Toolbox\
•	MATLAB's Statistics and Machine Learning Toolbox\


## Support
For questions or technical issues please contact Myriam Lubrano at myriam.lubrano2@unibo.it.



