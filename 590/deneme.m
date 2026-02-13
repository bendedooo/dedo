Re  = 6378.2;
f = 0.00335;
TimeStamp = TimeStamp0;
TimeStampArray = [year(TimeStamp)  month(TimeStamp) day(TimeStamp) hour(TimeStamp) minute(TimeStamp) second(TimeStamp)];
O_ECI2ECEF =  dcmeci2ecef('IAU-2000/2006',TimeStampArray);% ???? Convert Earth-centered inertial (ECI) to Earth-centered Earth-fixed (ECEF) coordinates for January 12, 2000 at 4 hours, 52 minutes, 12.4 seconds and January 12, 2000 at 4 hours, 52 minutes, and 13 seconds. Use the IAU-2000/2006 reduction.
r_ECEF = O_ECI2ECEF * r_ECI_0;
lla = ecef2lla(r_ECEF',f,Re)