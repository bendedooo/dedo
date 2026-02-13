clear all; clc; close all;

h = 6000; %[ft]
mph2fts = 5280/3600; % mph to ft/s
v = (0.003/8.25)^-0.5 %[ft/s]

vz = 8.25-0.3*v+0.003*v^2 %[mph]
t = h/vz %[sec]
x = v*mph2fts*t  %[ft]
x = x/5280 %[miles]
