clear; clc;

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

PARALLEL_BOUNDARY = 7;


% Set this
x = 8;

% Dimensions for the hall and the wings
% l=length (x), w=width (y), h=height (z), d=distance
hall_l = 32*x;
hall_w = 4*x;
hall_h = 1;  

wing_l = 2*x;
wing_w = 8*x;
wing_h = 1;
wing_d = 4*x;

% Canvas size
Nx = hall_l;
Ny = wing_w - 1 + hall_w - 1 + wing_w;
Nz = hall_h;

% Initialize the canvas (MI building)
MI = NO_SLIP * ones(Nx, Ny, Nz);

% Construct the main hall
Hall = NO_SLIP * ones(hall_l, hall_w, hall_h);
Hall(2:hall_l-1, 2:hall_w-1, :) = FLUID;

% Add the doors
% North door (main entrance)
x_pos = 1 + 2*x;
y_pos = 1;
door_size = x;
Hall(x_pos:x_pos+door_size, y_pos, :) = INFLOW_1;

% East door (LRZ)
x_pos = 1;
y_pos = hall_w - 1 - 2*x;
door_size = x;
Hall(x_pos, y_pos:y_pos+door_size, :) = OUTFLOW;

% Cantine door
x_pos = 1 + wing_l + 2*(wing_l + wing_d) + wing_l + x;
y_pos = hall_w;
door_size = x;
Hall(x_pos:x_pos+door_size, y_pos, :) = INFLOW_2;

% South door (side door, library)
x_pos = hall_l-3*x-1;
y_pos = hall_w;
door_size = x;
Hall(x_pos:x_pos+door_size, y_pos, :) = OUTFLOW;

% West door (back entrance)
x_pos = hall_l;
y_pos = 1 + x + 1;
door_size = x;
Hall(x_pos, y_pos:y_pos+door_size, :) = OUTFLOW;

% Add the tables (at least x=4 is required)
if (x>=4)
    table_l = x/4;
    table_w = x/2;
    table_d = x;
    
    % North-east set (Rechnerhalle)
    start_x = wing_l + wing_d + wing_l + x;
    end_x = start_x + table_l;
    start_y = 1 + x;
    end_y = start_y + table_w;
    Hall(start_x:end_x, start_y:end_y, :) = NO_SLIP;

    for i=1:5
        start_x = start_x + table_d;
        end_x = end_x + table_d;
        Hall(start_x:end_x, start_y:end_y, :) = NO_SLIP;
    end

    % North-west set (Tree)
    start_x = start_x + wing_l + wing_d;
    end_x = start_x + table_l;

    for i=1:5
        start_x = start_x + table_d;
        end_x = end_x + table_d;
        Hall(start_x:end_x, start_y:end_y, :) = NO_SLIP;
    end

    % South set (Cantine)
    start_x = 2*(wing_l + wing_d) + x;
    end_x = start_x + table_l;
    start_y = hall_w -1 -x - table_w + 1;
    end_y = start_y + table_w;
    Hall(start_x:end_x, start_y:end_y, :) = NO_SLIP;

    for i=1:10
        start_x = start_x + table_d;
        end_x = end_x + table_d;
        Hall(start_x:end_x, start_y:end_y, :) = NO_SLIP;
    end
end

% Put the main hall in place
MI(: , wing_w:wing_w+hall_w-1, : ) = Hall;

% Construct a north wing
WingN = NO_SLIP * ones(wing_l, wing_w, wing_h);
WingN(2:wing_l-1, 2:wing_w-1, :) = FLUID;
WingN(2:wing_l-1, wing_w, :) = FLUID;
% WingN(2:wing_l-1, wing_w, :) = PARALLEL_BOUNDARY;

% Put the north wings in place
xpos = wing_l + wing_d + 1; % one wing + distance + east boundary
ypos = 1;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

% Construct a south wing
WingS = NO_SLIP * ones(wing_l, wing_w, wing_h);
WingS(2:wing_l-1, 2:wing_w-1, :) = FLUID;
WingS(2:wing_l-1, 1, :) = FLUID; 
% WingS(2:wing_l-1, 1, :) = PARALLEL_BOUNDARY; 

% Put the south wings in place
xpos = wing_l + 1; % one wing + east boundary
ypos = wing_w -1 + hall_w;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

xpos = xpos + wing_l + wing_d;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

% Visualize the domain
image(10*MI(:,:,1));
axis equal
xlim([0 Ny + 1]) % x for the graph only
ylim([0 Nx + 1]) % y for the graph only

% Write pgm files
% Main Hall
n_cpus = 4;
hall_part_l = hall_l / n_cpus;
if (mod(hall_l, n_cpus) ~= 0)
    printf('Careful! Hall not divisible by number of cpus!')
end
x_start = 1;
x_end = hall_part_l;
y_start = wing_w;
y_end = y_start + hall_w -1;
% ATTENTION: For this, first assign the boundary of each wing to be parallel boundary and take a copy of the wing files.
% I will make it automatic later.
part_0 = [MI(x_start:x_end, y_start:y_end, :); PARALLEL_BOUNDARY*ones(1, hall_w, hall_h)];
dlmwrite('pgm/cpu_0.pgm', part_0, 'delimiter', ' ')

x_start = x_start + hall_part_l;
x_end = x_end + hall_part_l;
part_1 = [PARALLEL_BOUNDARY*ones(1, hall_w, hall_h); MI(x_start:x_end, y_start:y_end, :); PARALLEL_BOUNDARY*ones(1, hall_w, hall_h)];
dlmwrite('pgm/cpu_1.pgm', part_1, 'delimiter', ' ')

x_start = x_start + hall_part_l;
x_end = x_end + hall_part_l;
part_2 = [PARALLEL_BOUNDARY*ones(1, hall_w, hall_h); MI(x_start:x_end, y_start:y_end, :); PARALLEL_BOUNDARY*ones(1, hall_w, hall_h)];
dlmwrite('pgm/cpu_2.pgm', part_2, 'delimiter', ' ')

x_start = x_start + hall_part_l;
x_end = x_end + hall_part_l;
part_3 = [PARALLEL_BOUNDARY*ones(1, hall_w, hall_h); MI(x_start:x_end, y_start:y_end, :)];
dlmwrite('pgm/cpu_3.pgm', part_3, 'delimiter', ' ')

% North wings
WingN_par = [WingN PARALLEL_BOUNDARY * ones(wing_l, 1, wing_h)];
dlmwrite('pgm/cpu_4.pgm', WingN_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_5.pgm', WingN_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_6.pgm', WingN_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_7.pgm', WingN_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_8.pgm', WingN_par, 'delimiter', ' ')

% South wings
WingS_par = [PARALLEL_BOUNDARY * ones(wing_l, 1, wing_h) WingS];
dlmwrite('pgm/cpu_9.pgm',  WingS_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_10.pgm', WingS_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_11.pgm', WingS_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_12.pgm', WingS_par, 'delimiter', ' ')
dlmwrite('pgm/cpu_13.pgm', WingS_par, 'delimiter', ' ')
