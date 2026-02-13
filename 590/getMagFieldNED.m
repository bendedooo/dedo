function out = getMagFieldNED(t,r_ECI,TimeStamp0)
Re  = 6378.2;
f = 0.00335; % flattening constant

TimeStamp = TimeStamp0 + seconds(t);
TimeStampArray = [year(TimeStamp)  month(TimeStamp) day(TimeStamp) hour(TimeStamp) minute(TimeStamp) second(TimeStamp)];
% O_ECI2ECEF =  dcmeci2ecef('IAU-2000/2006',TimeStampArray);% ???? Convert Earth-centered inertial (ECI) to Earth-centered Earth-fixed (ECEF) coordinates for January 12, 2000 at 4 hours, 52 minutes, 12.4 seconds and January 12, 2000 at 4 hours, 52 minutes, and 13 seconds. Use the IAU-2000/2006 reduction.
% r_ECEF = O_ECI2ECEF * r_ECI;
r_ECEF = eci2ecef(TimeStampArray, r_ECI);
lla = ecef2lla(r_ECEF',f,Re); %  [lat(degrees) lon(degrees) meters] r_ECEF in km bu geodetic mi0.00335 flattening constant
% lat = lla(1);
% lon = lla(2);
% alt = lla(3)*1000;
dyear = decyear(TimeStampArray(1),TimeStampArray(2),TimeStampArray(3));
[XYZ_NED,~,~,~,~] = igrfmagm(lla(3)*1000,lla(1),lla(2),dyear,13); % function output row vector, input [meter deg deg]
O_ECEF2NED = Angle2DCM_321(lla(2)*pi/180,-pi/2-lla(1)*pi/180,0);
XYZ_ECEF = O_ECEF2NED' * XYZ_NED';
[XYZ_ECI] = ecef2eci(TimeStampArray,XYZ_ECEF);%O_ECI2ECEF' * XYZ_ECEF;
O_ECI2P = eye(3); % Angle2DCM_313(omega(k),i(k),Omega(k));
XYZ_P = O_ECI2P * XYZ_ECI;

out = [XYZ_NED' XYZ_ECEF XYZ_ECI XYZ_P];




end