%% Hw4-1
clear;clc;

% Read file
filename = "Hw4-1.xls"; sheet = "Sheet1"; 
data = xlsread(filename, sheet);
I_data = data(:,1:2); O_data = data(:,3);

R_d = ( (I_data')*I_data)\( (I_data')*O_data);
disp(R_d);