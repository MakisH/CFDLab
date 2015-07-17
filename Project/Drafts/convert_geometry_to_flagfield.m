close all; clear; clc;

%% Define the flags 
% % Flagfield definitions (original)
% FLUID           = 0;
% NO_SLIP         = 1;
% MOVING_WALL     = 2;
% FREE_SLIP       = 3;
% INFLOW          = 4;
% OUTFLOW         = 5;
% PRESSURE_IN     = 6;

% Flagfield definitions (multiple boundaries)
FLUID           = 0;
NO_SLIP         = 1;
FREE_SLIP       = 2;
OUTFLOW         = 3;

MOVING_WALL     = 10;

INFLOW          = 30;
INFLOW_1        = 31;
INFLOW_2        = 32;
INFLOW_3        = 33;
INFLOW_4        = 34;
INFLOW_5        = 35;

PRESSURE_IN     = 50;
PRESSURE_IN_1   = 51;
PRESSURE_IN_2   = 52;
PRESSURE_IN_3   = 53;
PRESSURE_IN_4   = 54;
PRESSURE_IN_5   = 55;

%% Input the 2D geometry (e.g. from ImageJ)
% 2D geometry file from ImageJ (text image)
% fileID          = fopen('geometry_2D.txt','r');

A = importdata('geometry_2D.txt');
res_file_x      = size(A,2);
res_file_y      = size(A,1);

% Scale/map to {0,1}
A = ceil(A / 255);

%% Physical dimensions
px_ratio = 4.464; % px/m

%% Crop to the area of interest
% manual crop
crop_x_min = 127;
crop_x_max = 825;
crop_y_min = 226;
crop_y_max = 350;

% TODO: autocrop

B = A( crop_y_min : crop_y_max , crop_x_min : crop_x_max );


%% Boundary conditions

start_x = 2;
start_y = 3;

% Door from x=10m to x=12m at y=0
B(0 + start_y : 2 + start_y, 10*px_ratio + start_x : 12*px_ratio + start_x) = INFLOW;


%% TODO: Scale
% res_x           =1000;
% scale           = res_x / res_file_x;
% res_y           = floor(scale * res_file_y);

%% Visualize
image(10*B);
axis equal
xlim([0 crop_x_max-crop_x_min])
ylim([0 crop_y_max-crop_y_min])
