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

% Set these
x = 8;
max_flag = 55;

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
x_pos = 1 + 2*x + 1; % East boundary + 2*x + start
y_pos = 1;
door_size = x;
Hall(x_pos:x_pos+door_size-1, y_pos, :) = INFLOW_1;

% East door (LRZ)
x_pos = 1;
y_pos = hall_w - 2*x; % hall width - 2*x
door_size = x;
Hall(x_pos, y_pos:y_pos+door_size-1, :) = OUTFLOW;

% Cantine door
x_pos = 1 + wing_l + 2*(wing_l + wing_d) + wing_l + x;
y_pos = hall_w;
door_size = x;
Hall(x_pos:x_pos+door_size-1, y_pos, :) = INFLOW_2;

% South door (side door, library)
x_pos = hall_l-3*x;
y_pos = hall_w;
door_size = x;
Hall(x_pos:x_pos+door_size-1, y_pos, :) = OUTFLOW;

% West door (back entrance)
x_pos = hall_l;
y_pos = 1 + x + 1;
door_size = x;
Hall(x_pos, y_pos:y_pos+door_size-1, :) = OUTFLOW;

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
% WingN(2:wing_l-1, wing_w, :) = FLUID;
WingN(1:wing_l, wing_w, :) = PARALLEL_BOUNDARY;

% Put the north wings in place
xpos = wing_l + wing_d + 1; % one wing + distance + east boundary
ypos = 1;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 

for i=2:5
    xpos = xpos + wing_l + wing_d;
    MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingN; 
end

% Construct a south wing
WingS = NO_SLIP * ones(wing_l, wing_w, wing_h);
WingS(2:wing_l-1, 2:wing_w-1, :) = FLUID;
% WingS(2:wing_l-1, 1, :) = FLUID; 
WingS(1:wing_l, 1, :) = PARALLEL_BOUNDARY; 

% Put the south wings in place
xpos = wing_l + 1; % one wing + east boundary
ypos = wing_w -1 + hall_w;
MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 

for i=2:5
    xpos = xpos + wing_l + wing_d;
    MI(xpos:xpos+wing_l-1, ypos:ypos+wing_w-1, :) = WingS; 
end

% Visualize the domain
image(rot90(10*MI(:,:,1)));
axis equal
xlim([0 Nx + 1]) % x for the graph only
ylim([0 Ny + 1]) % y for the graph only

% Write pgm files
% Main Hall
n_cpus = 4; % Number of partitions-cpus for the main hall
hall_part_l = hall_l / n_cpus;
if (mod(hall_l, n_cpus) ~= 0)
    printf('Careful! Hall not divisible by number of cpus!')
end
x_start = 1;
x_end = hall_part_l;
y_start = wing_w;
y_end = y_start + hall_w -1;
part_0 = rot90([MI(x_start:x_end, y_start:y_end, :); PARALLEL_BOUNDARY*ones(1, hall_w, hall_h)]);
f_name = 'pgm/cpu_0.pgm';
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(part_0,2)), ' ', num2str(size(part_0,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite('pgm/cpu_0.pgm', part_0, '-append', 'delimiter', ' ')

file_id = 0;

for i=1:n_cpus-2
    file_id = file_id + 1;
    x_start = x_start + hall_part_l;
    x_end = x_end + hall_part_l;
    clear part_i
    part_i = rot90([PARALLEL_BOUNDARY*ones(1, hall_w, hall_h); MI(x_start:x_end, y_start:y_end, :); PARALLEL_BOUNDARY*ones(1, hall_w, hall_h)]);
    f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
    f_id = fopen(f_name,'w');
    fprintf(f_id, '%s\n', ['P3 ', num2str(size(part_i,2)), ' ', num2str(size(part_i,1)), ' ', num2str(max_flag)]);
    fclose(f_id);
    dlmwrite(f_name, part_i, '-append', 'delimiter', ' ')
end

file_id = file_id + 1;
x_start = x_start + hall_part_l;
x_end = x_end + hall_part_l;
part_n = rot90([PARALLEL_BOUNDARY*ones(1, hall_w, hall_h); MI(x_start:x_end, y_start:y_end, :)]);
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(part_n,2)), ' ', num2str(size(part_n,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(['pgm/cpu_',num2str(file_id),'.pgm'], part_n, '-append', 'delimiter', ' ')

% North wings
WingN(1, wing_w, :) = NO_SLIP;
WingN(2:wing_l-1, wing_w, :) = FLUID;
WingN(wing_l, wing_w, :) = NO_SLIP;
WingN_par = rot90([WingN PARALLEL_BOUNDARY * ones(wing_l, 1, wing_h)]);
file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingN_par,2)), ' ', num2str(size(WingN_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingN_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingN_par,2)), ' ', num2str(size(WingN_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingN_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingN_par,2)), ' ', num2str(size(WingN_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingN_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingN_par,2)), ' ', num2str(size(WingN_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingN_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingN_par,2)), ' ', num2str(size(WingN_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingN_par, '-append', 'delimiter', ' ')

% South wings
WingS(1, 1, :) = NO_SLIP;
WingS(2:wing_l-1, 1, :) = FLUID; 
WingS(wing_l, 1, :) = NO_SLIP;
WingS_par = rot90([PARALLEL_BOUNDARY * ones(wing_l, 1, wing_h) WingS]);
file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingS_par,2)), ' ', num2str(size(WingS_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name,  WingS_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingS_par,2)), ' ', num2str(size(WingS_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingS_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingS_par,2)), ' ', num2str(size(WingS_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingS_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingS_par,2)), ' ', num2str(size(WingS_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingS_par, '-append', 'delimiter', ' ')

file_id = file_id + 1;
f_name = ['pgm/cpu_',num2str(file_id),'.pgm'];
f_id = fopen(f_name,'w');
fprintf(f_id, '%s\n', ['P3 ', num2str(size(WingS_par,2)), ' ', num2str(size(WingS_par,1)), ' ', num2str(max_flag)]);
fclose(f_id);
dlmwrite(f_name, WingS_par, '-append', 'delimiter', ' ')

