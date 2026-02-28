%% ALEPH TOWER: Parametric Digital Twin --> MATLAB is the Parametric Modifier
%  Parameter: MeshSize (0-100%), where 100% = Reference Paper Fidelity (0.5m Floor Slab Mesh/0.3m Shear Wall Mesh).
clear; clc; close all;
% --- 1. PARAMETERS ---
MeshSize = 75;    % Mesh Size Parameter -->Parameter 1 (MATLAB), 100 = Paper Fidelity (200m/52 stories =4.4234m height/floor, 4.4234/4=0.9m height/wall). 75 = Optimized (3 subdivs).
VizScale = 5000;  % Visual deformation scale --> Parameter 2, simple magnifier to view modal displacement
% Calculate subdivisions based on percentage (Baseline 100% = 4, atleast 4 ELEMENTS, ie divisions in wall mesh to calculate torision/bending, industry standard)
subdiv = max(1, round(4 * (MeshSize / 100))); 
work_dir = fullfile(pwd, 'ANSYS__Percentage');
if ~exist(work_dir, 'dir'), mkdir(work_dir); end
% --- 2. DEFINE SIMULATION LOGIC (BLOCKS) ---
% Block 1: Physics & Materials
setup_cmds = {
    '/PREP7', 'ANTYPE,MODAL', 'MODOPT,LANB,10', 'MXPAND,10', ... %Setup Phase: Preprocessor set to Modal Analysis Type, with Modal Option for 10 modes using Lanczos Algorithm with Modal Expansion (Displacement) Calculation
    'ET,1,SHELL181', 'SECTYPE,1,SHELL', 'SECDATA,0.4', ...    % Parameters 3, 4 & 5 --> SECDATA - Column Core Wall thickness, Cross section and slab thickness, over next 3 lines. Core Wall Definition: Setting Element Type Shell 181 for in plane stresses and external loads, 0.4m Central Core Shear Wall thickness from reference data
    'ET,2,BEAM188',  'SECTYPE,2,BEAM,RECT', 'SECDATA,0.8,0.8', ... % Column Definition: Setting Element Type Beam 188 for 3D Shearing, with Section Data including L/W/T 3D dimensions, 0.8m Central Load Bearing Column Cross Section from reference data
    'ET,3,SHELL181', 'SECTYPE,3,SHELL', 'SECDATA,0.22', ...   % Slab Definition: Shell 181, used this time for horizontal floor diaphragm-column load transfer.0.22m Floor Slab thickness from reference data
    'MP,EX,1,50e9', 'MP,PRXY,1,0.2', 'MP,DENS,1,2400', ... % Parameters 6 & 7 --> Material Property: Young's Elasticity Modulus for Upper and Base Level Concrete Types and floor over next 3 lines. density and Poisson's Ratio, C80 Concrete More Expensive for Base, think shearing principle at base/viscosity at skin friction level
    'MP,EX,2,44e9', 'MP,PRXY,2,0.2', 'MP,DENS,2,2400', ... % Material Property: Young's Elasticity Modulus, density and Poisson's Ratio, C60 Concrete --> Cheaper and Less stiff for lower walls
    'MP,EX,3,30e9', 'MP,PRXY,3,0.2', 'MP,DENS,3,2500'};  % Material Property: Young's Elasticity Modulus, density and Poisson's Ratio, Low Grade General --> Column takes central load, low cost diagpragm primarily? takes torisional loads
% Block 2: Geometry (Nodes & Elements) - CORRECTED TO FILL CENTER
geo_cmds = {
    'N,1,-6,-6,0', 'N,2,6,-6,0', 'N,3,6,6,0', 'N,4,-6,6,0', ... % Parameters 7 & 8 --> Column and base footprint in coordinates to be extruded later. 12m*12m central column base centered at 0,0, N marks nodal point of structural coordinates
    'N,5,-15,-15,0', 'N,6,15,-15,0', 'N,7,15,15,0', 'N,8,-15,15,0', ... % 30m*30m base centered at 0,0, 6 and 15 indicate Symmetry of Structure
    sprintf('NGEN,%d,1000,1,8,1,0,0,%f', (52*subdiv)+1, 4.423/subdiv), ... % Parameters 9 & 10 --> Geometry Parameter (MATLAB): 52*subdiv is Vertical Element Count, 4.423/subdiv is logically 4.423/subdiv. 4.423 is height/floor. subdiv also meshes the system. sprintf is MATLAB Function, not APDL, sets number string of Nodes to generate from 4 base nodes, %d & %f--> Integer Placeholder and Float Point, i.e. Decimal Placeholder, Relates to Subdiv Mesh Size Parameter to define Meshing for 52 Floors. NGEN, ITIME, INC, NODE1, NODE2, NINC, DX, DY, DZ--> Node Generation Default String, ITIME Set Function in APDL, linked parametrically to %d, and DZ is Z Offset, set parametrically through %f, DX and DY base and width offsets. NINC is node range count, ie, nodal increments from 1 to 8(4 base+4 column).
    'TYPE,1', 'MAT,1', 'SECNUM,1', 'E,1,2,1002,1001', 'E,2,3,1003,1002', 'E,3,4,1004,1003', 'E,4,1,1001,1004', ... % E is element linking nodes. 1-->2-->3--4 Is 2d Square base, 1-->1001 is wall height at one corner from frame of reference, and 1001-->1002 is Ceiling Base/Upper Floor linkage. 1000 arbitrary nodal frame to set nodal offsets for 52 floors.
    'TYPE,2', 'MAT,1', 'SECNUM,2', 'E,5,1005', 'E,6,1006', 'E,7,1007', 'E,8,1008', ... %Type selects element type, recalling what was set in ET, Mat is material set before, SECNUM derives from material property dimensions defined in cecdata in block 1.
    'TYPE,3', 'MAT,3', 'SECNUM,3', ...
    'E,1001,1002,1006,1005', 'E,1002,1003,1007,1006', ... % 5, 6, 7 and 8 are the other side of it.
    'E,1003,1004,1008,1007', 'E,1004,1001,1005,1008', ...
    'E,1001,1002,1003,1004'}; % <--- NEW: CENTRAL CORE CAP (Fills the hollow hole)
% Block 3: Extrusion & BCs
bc_cmds = {
    'ALLSEL,ALL', sprintf('EGEN,%d,1000,ALL', 52*subdiv), ...%Select all elements and extrude by generating elements, 52*Subdiv is primary mesh size scaler for 27 floors, 52*3=156, ALL is ELIST, copies every element
    'NSEL,S,LOC,Z,-0.01,0.01', 'D,ALL,ALL,0', 'ALLSEL,ALL', ... % Parameters 10 & 11--> Defines height (z coord) at base and and material change point in next line. Node select all nodes +-0.01m means only nodes on ground, LOC is APDL 3D Nodal Point selector, here for z axis. After selection, Displacement is set to 0.ALLSEL ALL then selects the whole building again after fixing base.
    'ESEL,S,TYPE,,1', 'ESEL,R,CENT,Z,77.9,200', 'MPCHG,2,ALL', 'ALLSEL,ALL'}; %Select Element with S to select type 1 element set in block 1, and with R to reselect from centroid, ie geometric center onwards to mark a material change of concrete type using MPCHG to 2, ie material 2 defined in block 1. From reference data, first 20 floors need most strength, principally true. 20*4.4239=77.9m. Centroid indicates to select column walls above that., Allsel selects entire figure, priming it for solver.
% Block 4: Solver & Data Export
export_cmds = {
    '/SOLU', 'SOLVE', 'FINISH', '/POST1', 'SET,LIST', ... % Enter Solution Solver engine and initiate solver, Finish cleanup of temp files after results and exit solver. Post is the result entry log viewer, with setlist indicating what result to parse for study. Antype modal in block 1 set the solver for generalized eigenvalue problem at 0 damping assumption for simple hm. Solves for k, m, omega^2-->Eigenvalues, omega-->nf, and associated eigenvectors are the mode shapes
    '*CFOPEN,results,txt', '*DO,m,1,10', '*GET,f,MODE,m,FREQ', '*VWRITE,m,f', '(F5.0,F15.8)', '*ENDDO', '*CFCLOS', ... % Parameter 11--> Mode count loop setting. open command to write results to text, Do APDL Logic loop to record mode results 10 times for 10 modes, set with m, repeat for frequency, and write vector(eigen)into txt.E
    '*CFOPEN,geo,csv', '*GET,n,NODE,0,NUM,MAX', '*DO,i,1,n', '*IF,NSEL(i),EQ,1,THEN', ... %Write Geometry File, Iteratively (i) get node locations from 1 to n (max/total nodeid), If selected node i=1 (binary, ie, does node exist), then,
    '*GET,x,NODE,i,LOC,X', '*GET,y,NODE,i,LOC,Y', '*GET,z,NODE,i,LOC,Z', ... % Get nodal xyz value for that node index number
    '*VWRITE,i,x,y,z', '(F8.0,3F10.3)', '*ENDIF', '*ENDDO', '*CFCLOS', ... % Write Vector Data for a fixed point f floating point number with 8 maximum numbers, and 0 decimals. This is for nodal index vector data., and 3f10.3 is 3 variables x,y and z 10 numbers with max 3 decimals mm resolution. Endif follows from nsel loop, until loop runs out and endo is triggered by closing the commmand file with cfclos to end loop and move to next part.
    '*CFOPEN,modes,csv', '*DO,m,1,10', 'SET,1,m', '*DO,i,1,n', '*IF,NSEL(i),EQ,1,THEN', ... %Same loop but for modal displacement, ie tracking mode shapes whole using xyz tracking.
    '*GET,x,NODE,i,U,X', '*GET,y,NODE,i,U,Y', '*GET,z,NODE,i,U,Z', ...
    '*VWRITE,m,i,x,y,z', '(F5.0,F8.0,3F15.8)', '*ENDIF', '*ENDDO', '*ENDDO', '*CFCLOS', ...
    '*CFOPEN,elems,csv', '*GET,e,ELEM,0,NUM,MAX', '*DO,i,1,e', '*GET,t,ELEM,i,ATTR,TYPE', ... %Same loop to track element shapes constructed nodally in xyz to give initial skeleton to compare against, ELEM is Topology library containing nodal data and element attributes attr for n1,n2...nx. e is the element, attr is attribute type attached to element [4 attribute types - ET,MP, SECNUM, Real Numbers - R,1,x ie. 1 is boolean and x is the dimension/size like 0.5m for example.], basically descriptor of element, ie, describes here only element type et from block 1.
    '*IF,t,GT,0,THEN', '*GET,n1,ELEM,i,NODE,1', '*GET,n2,ELEM,i,NODE,2', '*GET,n3,ELEM,i,NODE,3', '*GET,n4,ELEM,i,NODE,4', ...%If element type t>0, then record element. Element type is the property like wall, floor etc. Ensures structural existence and therefore continuity.
    '*VWRITE,i,t,n1,n2,n3,n4', '(6F10.0)', '*ENDIF', '*ENDDO', '*CFCLOS', 'FINISH'}; %Therefore, if element exists, get node, from nodal list saved in ELEM, for an iterative index i until t=0, then end like before.
% --- 3. WRITE & RUN ---
fid = fopen(fullfile(work_dir, 'model.inp'), 'w'); % Matlab Functions now: Create and open a 'fullfile'-->Platform agnostic readable file directory path, in the working directory to write this fem model's input file.
fprintf(fid, '%s\n', setup_cmds{:}, geo_cmds{:}, bc_cmds{:}, export_cmds{:}); %Serial Write command to take string, move to new line/n. Setup-->Block 1, Geometry--> Block 2, Boundary Conditions --> Block 3, Data Export Parameters -->Block 4. : Opens Cell array, extracts data and comma separates them.
fclose(fid); % Close stream and flush data buffer to disk
fprintf('<strong>[SYSTEM] Solving (Mesh Quality: %d%%)...</strong>\n', MeshSize); %Prints message that solver is running for selected mesh refinement for %d--> Mesh size as integer, extracted from MeshSize, and the /n newline. One % is a command character to insert variable like d. Double %% defines that it should print one percent symbol
exe_list = dir('C:\Program Files\ANSYS Inc\v*\ansys\bin\winx64\ansys*.exe'); % v* indicates highest/latest version to select
if isempty(exe_list), error('ANSYS Executable not found in default path.'); end
ansys_exe = fullfile(exe_list(end).folder, exe_list(end).name); % Pick the last one (highest version)
[~,~] = system(sprintf('taskkill /F /IM "%s" /T', exe_list(end).name)); % Kill processes using ~ ~ output suppression operation to ignore status/result text messages. system specifies string inside to be executed in cmd, /F ignores yes/no permission prompts, and /IM targets the image, ie executable to target for killing. /T instructs entire tree to be killed. These are windows flags
system(sprintf('"%s" -b -dis -np 16 -dir "%s" -j job1 -i model.inp -o out.out', ansys_exe, work_dir)); % Ansys flags: -b: Batch Mode, distributed memory load over ram slots, 16 parallel processors, -dir ("%s" = work_dir) specified at the end, -j creates the job file for processing, split between database and results. -i defines input model as model.inp, which was set to hold all apdl commands from 4 block variables, -o writes output log from ANSYS console

% --- 4. VISUALIZER ---

if exist(fullfile(work_dir, 'results.txt'), 'file') %Checker for file existence
    freqs = load(fullfile(work_dir, 'results.txt')); % Calculated Frequencies stored in results, 10x2 matrix, .csv's were defined in block 4
    nodes = load(fullfile(work_dir, 'geo.csv'));  % Nodes stored in geometry defined in block 2, static xyz coordinates
    modes = load(fullfile(work_dir, 'modes.csv')); %Vibration vectors for nodes Ux,Uy & Uz @ each node, i.e. Nodal Displacement
    elems = load(fullfile(work_dir, 'elems.csv')); %Provides element skeleton of building to compare against
    
    map = containers.Map(nodes(:,1), 1:size(nodes,1)); % Refer to nsel commands in block 4. For 1000 as the base for nodal indexer, not every value is stored -->Sparse Data. These are 0 and MATLAB needs to fill array between the real indices. --> Dense Data,  % Map makes lookup map table to map for empty parts. nodes(:,1) --> Select corresponding row for each nodal column. 1:size(...) --> Size Command measures matrix size in matlab and fills in consequtive integer value
    idx = @(id) arrayfun(@(x) map(x), id); % @ is used to define independant variable for use, id is arbitrary name for variable that will be fed in to the function idx when called. arrayfun is the actual function. Looks up number x (@x) fromvector map of x. x is iterator variable. IMPORTANT: Map(x) gives scalar indexed number by address row to scalar nodal id like 1001. Then arrayfun strings the indexes into single vector column
    
    % Column Core Vizualiser-
    get_faces = @(t) [idx(elems(elems(:,2)==t,3)), idx(elems(elems(:,2)==t,4)), idx(elems(elems(:,2)==t,5)), idx(elems(elems(:,2)==t,6))]; % Function to Parse element types out from t (column 2) defined in block 4 line 49. Indices 3 through 6. Parsing mechanism takes element type 1181 of column IN NEXT LINE and maps to represent the 4-node connectivity for Shel181 elements (n1,n2,n3, n4 of columns 3,4,5,6) as defined in line 49, mapping them into the dense vertex array for patch rendering.
    s_core = get_faces(1); s_slab = get_faces(3); % Calls function above, filtering by type 1 (ET defined line 15) . Filters for 2 parts of building, Type 1 = Vertical Shear High Stress Bearing Core (Red/Orange), Type 3 = Diaphragm Slab (Cyan)
    
    b = elems(elems(:,2)==2, 3:4); b_idx = [idx(b(:,1)), idx(b(:,2))]; %Analogous process for element type 2 - Beam, == Is lhs rhs binary check, Column defined as 1d lines in model. Beam 188 property, thickness is treated as an inherent parameter using secdata in block 1 (0.8*0.8). Defines 2 nodes connecting beams in 1d.
    fprintf('\n<strong>[RESULTS] 10-Mode Scan (Mesh: %d%%):</strong>\n', MeshSize); % Defines command window output
    for k=1:10, fprintf(' Mode %2d: %.4f Hz\n', k, freqs(k,2)); end % Row k is frequency storage. Print mode 1:10 Results %2d is 2 integer number setting, %.4f Hz\n' prints frequency to 4dp, and then newline to print next in iteration. Data parsed from column 2 as defined in line 44.
    for m = 1:10 % Runs loop through u function for modes from 1-10
        U = modes(modes(:,1)==m, 3:5); % Parses modal matrix to find mode for current mode number.
        def = nodes(:, 2:4) + VizScale * U; %Modal displacement x y z real coordinates -->2, 3 & 4 from setup in line 48-49 of modes.csv MATLAB Row organization acts as indexer. Deformation is based on displacement, [Initial Displacement--> nodes(:, 2:4)] + [Scale (Arbitrary Vizscale) * Modal Displacement (U)] 
        is_twist = std(U(nodes(:,4)>190, 1)) > 0.1; % Boolean to measure sway at 10% displacement with standard deviation function std (Twist if yes, sway if no) nodes(:,4) > 190 --> Height z set to 4th column in line 47, 190m is height-40m, ~ Free end of cantilever, least rigidity, max sway/rotation amplitude, U operator selects Displacement matrix for Ux, ie U(..1) to study sway.
        
        u_mag = rssq(U, 2); % Norm to get scalar magnitude of displacement sway for STRESS HEATMAP Computation. For Vector Displacement matrix U, 1 is row, 2 is column. RSSQ (Root Sum Square) takes norm by summing squares of xyz displacement vectors in matrices
        stress_proxy = (max(u_mag) - u_mag) ./ max(u_mag); % Maximum displacement subtracted from U matrix (u_mag). At top, u is max from cantilever principle. Moving to base gives 0 displacement. Therefore subtraction gives inverse. That is, If the displacement is max at the top and min at the bottom, the result is an inversion of the U matrix. Its not inverted as 1/u_mag as that would be undefined, since matrix would get divided by 0 for base. Dividing by u_max normalizes for displacement
        
        f = figure('Name', sprintf('Mode %d', m), 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.9 0.6]); % Mode inserted dynamically by line 85 for loop, white background normalized to screen resolution by screen scaling to 1. Position is arbitrary
        
        subplot(1,3,1); view(3); title(sprintf('Mode %d: %.3f Hz (3D)', m, freqs(m,2))); % Subplot split to put more graphs on one window, view3 is 3D isopmetric View, and title for subplot
        draw_tower(def, s_slab, b_idx, 'c', 0.6); draw_tower(def, s_core, [], [0.85 0.33 0.1], 1); % Draw Cyan Slabs @ 60% Opacity against building core + Orange Core (RGB coordinates), drawn against modal displacement frame defined in def
        
        subplot(1,3,2);
        if is_twist || m >= 4, view(45, 60); title('Iso-Top (Torsion)'); % Parameter 12 --> Torision Detection Binary. Draws from Twisting Mode Boolean Operator in Case 8. Twisting generaly happens after 4th mode as motion causes twisting.. (Higher order or lower? Possibly lower. Study further). Sets isometric view if twisting, and front view at 2d angle for swaying.
        else, view(20, 20); title('Iso-Elev (Sway)'); end
        draw_tower(def, s_slab, b_idx, 'b', 0.4); draw_tower(def, s_core, [], 'r', 0.6); % Draw Blue Slabs + Red Core
        axis equal; grid on; if is_twist || m>=4, zlim([100 220]); end % zlim limits height angle, cutting the rest off from view
        
        subplot(1,3,3); title('Optimal PZT Zones (Stress Hotspots)'); view(3); axis equal; grid on; % Plotting for Piezoelectric Sensor Placement based on Modal Displacement defined in lines 84-85
        patch('Vertices', nodes(:,2:4), 'Faces', [s_slab; s_core], 'FaceVertexCData', stress_proxy, 'FaceColor', 'interp', 'EdgeAlpha', 0.05); % nodes(:2,4) defines vertices xyz data (Static, not def Modal Deformation coords), generates figure for PZT Location Stress Diagram defined in line 85
        colormap('jet'); c = colorbar; c.Label.String = 'Strain Proxy (Inverted Disp)'; clim([0.6 1]); % Square Voltage Cutoff Power [Square Law] proportional to V^2 and Strain^2. So if strain at 50%, Signal Power is 0.5^2, ie 0.25-->25%. 60-70% is a better floor
    end
else
    fprintf('[ERROR] Simulation Failed.\n'); % Part of if or else loop from line 69, where if there is no file, implication is solver failed. Prints failure to window
end

function draw_tower(p, s, b, c, a) % Positional argument, triggered by function draw_tower here(All parts are inputs for positional argument). p = def --> nodes(:,2:4), s = s_core --> s_core , b = b_idx -->  b_idx, c = 'c' --> [0.9 0.9 0.9], a = 0.8 --> 0.2. Note: Function labelled by positional argument = Variable 1 [Data call 1], and --> Variable 2 [Data call 2]
    hold on; grid on; axis equal;
    patch('Vertices', p, 'Faces', s, 'FaceColor', c, 'FaceAlpha', a, 'EdgeColor', 'k', 'EdgeAlpha', 0.1); % Patch indicates to color a section. Vertices are defined from xyz def matrix, faces are drawn over wall/floor slabs
    if ~isempty(b), plot3([p(b(:,1),1) p(b(:,2),1)]', [p(b(:,1),2) p(b(:,2),2)]', [p(b(:,1),3) p(b(:,2),3)]', 'k-'); end % Patch indicates to color a section. Vertices are defined from xyz def matrix. p is def is 3d coordinates and initial static nodal position from nodes.csv. b is beam vector indices of scalar identifiers, in both cases of the argument. [P,(...),1 is x coordinate from def, 2 is y and 3 z]. b(:1) is column 1 of row matrix of nodal id's. Therefore, xyz is defined for lwh 3d plot.
end

% --- 5. APDL SENSOR INJECTION & WIND PHYSICS ---
% This block defines the "Digital Twin" physics:
% 1. Piezoelectric Material (PZT-5H) --> Hard PZT, Low Cost
% 2. Electrical Buffer (10 MOhm Impedance) ... An Approximation to model the effect of JFET as a buffer to allow signal passage.
% 3. Critical Wind Speed---> From Strouhal number 0.15 for square shaped geometry (Vortex Lock-In Calculation)
% Note: ! is APDL Comment Syntax. Changing it will crash solver, so it is maintained as is.
fprintf('\n<strong>[PHYSICS] Injecting Sensor & Environmental Definitions...</strong>\n');

fid = fopen(fullfile(work_dir, 'sensor_physics.inp'), 'w');

% ---> A. Define Piezoelectric Material (PZT-5H) --> Adapted Data
fprintf(fid, '/PREP7 \n');              
fprintf(fid, 'ET, 900, SOLID226, 1001   ! Define Coupled-Field Solid (Volt+Disp)\n');
fprintf(fid, 'MP, PERX, 900, 930        ! Dielectric Permittivity X\n');
fprintf(fid, 'MP, PERY, 900, 930        ! Dielectric Permittivity Y\n');
fprintf(fid, 'MP, PERZ, 900, 857        ! Dielectric Permittivity Z\n');

fprintf(fid, 'TB, PIEZ, 900             ! Initialize Piezo Matrix\n');
fprintf(fid, 'TBDATA, 3, -7.15          ! e13 Piezo Constant (C/m2)\n');
fprintf(fid, 'TBDATA, 6, -7.15          ! e23 Piezo Constant\n');
fprintf(fid, 'TBDATA, 9, 13.70          ! e33 Piezo Constant\n');
fprintf(fid, 'TBDATA, 14, 11.90         ! e15 Piezo Constant\n');
fprintf(fid, 'TBDATA, 16, 11.90         ! e25 Piezo Constant\n');
% fprintf is printing to the inp file in 120
% 123 - Prep7 enters solver back to preprocessor
% Defines element type as 900 (Arbitrary), solid226 is a 3d brick property
% in apdl capable of displacement and voltage physics properties, Keypot 1001 indicates solver to read piezo materials matrix defined below.
% Dielectric Permittivity of x,y and z implies anisotropicity. Polarization potential in electric field.
% TB, PIEZ and 900 set the piezo matrix for data below
% a,b data is coubling coeffs in tensor logic. Remember anisotropicity
% All together encompass direct piezo effect

% ---> B. Nodal PZT Location Placement System
% Base targetted @ z=0, most high strain cantilever zone

fprintf(fid, 'NSEL, S, LOC, Z, 0, 0.5   ! Select nodes at the High-Stress Base\n');
fprintf(fid, 'NSEL, R, LOC, X, 0, 1     ! Refine to Core Wall\n');
fprintf(fid, '*GET, nPOS, NODE, 0, NUM, MIN ! Get ID of the "Positive" Terminal\n');
fprintf(fid, 'nNEG = nPOS               ! In 1-Node simplification, Ground is ref\n');
fprintf(fid, 'ALLSEL                    ! Restore full model\n');
% 147 - NSEL Node Selection for new set 'S', in Z axis loc for height 0-0.5m, ie geometrix base. Can be any arbitrary value < 8m ceiling/floor
% R is same logic for floor  portion on x axis, cuts out diaphragm to isolate core for sensor placement (To ignore floor vibrations)
% npos is set to select pos face on ground. 0 is arbitrary setting to select min numbered node of the nodes on ground. Strain against ground face is
% negative is same as positive. Same node is ground/wired to ground

% ---> C. The "JFET Buffer" (Modelled Resistor - 10MOhm) 
fprintf(fid, '*GET, eid_max, ELEM,, NUM, MAX  ! Find last element ID\n');
fprintf(fid, 'res_id = eid_max + 1            ! Create new ID for Circuit\n');
fprintf(fid, 'ET, res_id, CIRCU94, 0          ! Define Resistor Element\n');
fprintf(fid, 'R, res_id, 10E6                 ! Resistance = 10 MegaOhms\n');
fprintf(fid, 'TYPE, res_id \n');
fprintf(fid, 'REAL, res_id \n');
fprintf(fid, 'E, nPOS, 0                      ! Connect Sensor Node to Ground (0)\n');

fprintf(fid, '/SOLU \n');
fclose(fid);
% 160 - Resistance Element is set as Circuit using CIRCU94 command. 2 lines before take highest element number and add 1 to it for unique id.
% Resistance set to 10Mohm. Prevents instant signal decay. Slow drain. RC time constant. Slow charge accumulation
% APDL Element type pointer. Defines new element from 160 as circuit element for solver
% REAL is a pointer indicating JFET Constant Resistance.
% Node 0 is ground setting to complete circuit

% --- D. Environmental Loading (Vortex Shedding Physics) ---
% Calculates the Critical Wind Speed required to trigger the Fundamental Mode.
if exist('freqs', 'var')
    % 1. Get Fundamental Frequency
    f1 = freqs(1,2);  
    
    % 2. Aerodynamic Constants
    rho = 1.225;      % Air Density (kg/m^3)
    D_width = 30;     % Building Width (m)
    St = 0.14;        % Strouhal Number (0.14 for Square Prism)
    
    % 3. Calculate Critical Wind Speed (Lock-In Velocity) V = (f * D) / St--> f1 is Fundamental Mode, ie, vortex shedding frequency, St is 0.14 for fundamental square shape. v_lock is fundamental point at which sheet detaches (Von-Karman Effect) & @ St, v_lock = v_resonance, fshedding=fresonance, ie, speed at which fundamental frequencty hits resonance
    V_lock = (f1 * D_width) / St; 
    
    % 4. Calculate Vortex Force Magnitude
    % F = 0.5 * rho * V^2 * Area * C_L (Lift Coefficient ~ 0.2 -->Standard Square Shape Reference Value)
    Area = 30 * 230; % Width * Height
    F_vortex = 0.5 * rho * V_lock^2 * Area * 0.2; % Standard Aerodynamic Lift Formula. Used as it incorporates pressure differentials. a=f/m effect on building after.
    
    % Workspace Variables for Simulink Blocks
    assignin('base', 'V_lock_in', V_lock);
    assignin('base', 'F_vortex', F_vortex);
    
    % --- Print Physics Summary ---
    fprintf('   > Material: PZT-5H (SOLID226)\n');
    fprintf('   > Buffer Impedance: 10 Mega-Ohm (Simulated via CIRCU94)\n');
    fprintf('   > Sensor Placement: Auto-snapped to High-Stress Base (Z=0)\n');
    fprintf('   > Critical Wind Speed: %.2f m/s (%.1f km/h)\n', V_lock, V_lock*3.6);
    fprintf('   > Vortex Force: %.2f kN (Oscillating Load)\n', F_vortex/1000);
end

% --- 6. SIMULINK DTIL INTEGRATION (UPDATED WITH PZT PHYSICS) ---> Lumped Model of Approximate Characteristics
if exist('freqs', 'var')
    fprintf('\n<strong>[SYSTEM] Generating Simulink DTIL State-Space Model...</strong>\n');
    
    % --- 1. SYSTEM IDENTIFICATION (Modal Mass & Stiffness) ---
    f1 = freqs(1,2);  
    wn = 2 * pi * f1; % Rad/s conversion for shm
    zeta = 0.015;     % Damping (1.5% for Steel/Concrete -->Industry Standard)
    
    % Mass Calculation (Same as before)
    Vol_Slab = (30 * 30) * 0.22; 
    Vol_Core = (12 * 4) * 0.4 * 4.423;
    Vol_Col  = 4 * (0.8 * 0.8) * 4.423;
    Physical_Mass = (Vol_Slab + Vol_Core + Vol_Col) * 52 * 2400; 
    M_modal = Physical_Mass / 3; % 1/3rd of the top of cantilever vibrates freely with most mass, geometric factor. This fraction is modal mass
    
    % State-Space Matrices (Input: Force -> Output: Strain Proxy)
    % Note: C is scaled to output Strain (approx) rather than just displacement
    % Assumption max base strain is roughly proportional to tip displacement / Height --> Euler-Bernoulli Beam Theory
    A = [0 1; -(wn^2) -2*zeta*wn]; % Stiffness + Damping Variant of Newton's 2nd Law --> Ax+Bu --> Input Eqn --> Transient Characteristics, A is Stiffness+Damping Matrix (1st Derivative), internal vibration physics descriptor, B results in Inertial Acceleration (From mass movement rate derivative) rate of change for both.
    B = [0; 1/M_modal]; % Inertial Acceleration under Eqns of Motion Matrix, 0,1/m is that 2nd derivative has force, first one is velocity eqn, not acc and is thus 0 force. Finally, 1/m as ma/m=a
    C = [1/230 0]; % Approximation: Strain ~ Lateral Displacement / Height (230m) --> Building State Physical Output from Integration, ie total reading, ie, net displacement. [F-->A-->V-->U]
    D = 0; % Cx+Du--> Output Eqn -->Sensor Equation for Reading, Du is instantaneous sensor reading component. -->Strain depends on displacement, not force, hence 1/197, 0 , or anyways at an instant, velocity or position dont change due to instantaneous force before strain hits due to deformation. 
    % State Space --> Ax+Bu=Cx+Du
    % --- 2. VORTEX FORCING FUNCTION ---
    rho = 1.225; D_width = 30; St = 0.14;
    V_lock = f1 * D_width / St; 
    F_vortex = 0.5 * rho * V_lock^2 * (30 * 230) * 0.2; 
    
    % --- 3. PZT SENSOR PHYSICS (DRAWN FROM BLOCK 4.5) ---
    % These values mirror the "sensor_physics.inp" file exactly.
    
    % A. Material Constants (From Block 4.5 TBDATA/MP)
    d31_const = -190e-12;   % derived from e31 = -7.15 (Block 4.5)
    epsilon_r = 857;        % MP,PERZ (Block 4.5)
    epsilon_0 = 8.854e-12;
    
    % B. Buffer Impedance (From Block 4.5 "R, res_id, 10E6")
    R_buffer = 10e6;        % 10 Mega-Ohms (JFET Input Impedance)
    
    % C. Sensor Geometry (Standard Patch for High-Rise)
    L_pzt = 0.05;  % 50mm Length
    W_pzt = 0.03;  % 30mm Width
    t_pzt = 0.0005; % 0.5mm Thickness
    
    % D. Electrical Calculations (The Transfer Function Math)
    % Capacitance of the Sensor: C = (eps * A) / t
    C_pzt = (epsilon_r * epsilon_0 * L_pzt * W_pzt) / t_pzt;
    
    % RC Time Constant --> Used as a Frequency Cutoff Filter
    Tau_c = R_buffer * C_pzt; % R*C Dimensionally leaves seconds. R acts as buffer by limiting drain rate like pipe. Product gives percentage of value still in sensor for unit time. 
    
    % Voltage Sensitivity Gain (Volts per Unit Strain)
    % V = Strain * YoungsModulus * g31 * thickness... simplified to:
    PZT_Gain = 5000; % Approximate Voltage coupling coefficient for PZT-5H
    
    % --- 4. PARAMETER PUSH to Simulink Workspace---
    assignin('base', 'A_tower', A); % Base is matlab as base workspace. Simulink is differnt, variable name and System component. Pushes A-D to simulink
    assignin('base', 'B_tower', B);
    assignin('base', 'C_tower', C);
    assignin('base', 'D_tower', D);
    % The Transfer Function for High-Pass RC Filter: G(s) = Gain * (s / (s + 1/Tau)) -->s is d/dt, ie time derivative. 1/s is Integral by laplace domain transform e^st --> Study Further
    assignin('base', 'Num_PZT', [PZT_Gain 0]); % ---> Gain/[1/tau] = Transfer function. 1/tau to convert sec to hz for freq domain comparison.
    assignin('base', 'Den_PZT', [1 (1/Tau_c)]);%---> s needs to be greater than 1/t for signal to pass, or signal decays before reading, eqn in 266 is given here in matlab readable exponent power syntax. Organised by descending power
    
    % --- 5. SIMULINK MODEL --> Refined from Simulink Extraction
    model_name = 'Aleph_DTIL_Model';
    if bdIsLoaded(model_name), close_system(model_name, 0); end
    if exist([model_name '.slx'], 'file'), delete([model_name '.slx']); end
    
    new_system(model_name);
    open_system(model_name);
    
    % Block A: Source (Wind)
    add_block('simulink/Sources/Sine Wave', [model_name '/Vortex_Force']); %Sine Wave Block
    set_param([model_name '/Vortex_Force'], 'Amplitude', num2str(F_vortex), 'Frequency', num2str(wn)); %Set param injects the actual physics functions defined. Here vortex force and wn rad/s, ie forced vibration
    set_param([model_name '/Vortex_Force'], 'Position', [50, 100, 90, 140]);
    
    % Block B: Digital Twin (Structure)
    add_block('simulink/Continuous/State-Space', [model_name '/Structural_Twin']); %Simple name assigner called from simulink and assigned name for fcn below
    set_param([model_name '/Structural_Twin'], 'A', 'A_tower', 'B', 'B_tower', 'C', 'C_tower', 'D', 'D_tower'); %Sets state space for A--->D Input output system functions. Structural Twin Definer
    set_param([model_name '/Structural_Twin'], 'Position', [150, 95, 250, 145]);
    
    % Block C: PZT Sensor (Physics-Based Transfer Function)
    % This is where the Block 4.5 inputs (R=10M, C, e31) actually live.
    add_block('simulink/Continuous/Transfer Fcn', [model_name '/PZT_Buffer_Circuit']); %Simple name assigner called from simulink and assigned name for fcn below
    set_param([model_name '/PZT_Buffer_Circuit'], 'Numerator', 'Num_PZT', 'Denominator', 'Den_PZT'); % Sets Impedance Piezo Block with JFET modelled as constant Mohm. Func=Vout/strain=[(k*s)/(s+1/t)]
    set_param([model_name '/PZT_Buffer_Circuit'], 'Position', [310, 100, 420, 140]);
    
    % Block D: Output (Voltage)
    add_block('simulink/Sinks/Scope', [model_name '/PZT_Voltage_V']); %Simple name assigner called from simulink and assigned name for fcn below
    set_param([model_name '/PZT_Voltage_V'], 'Position', [480, 100, 520, 140]);
    % All name assigners create simulink blocks
    % Wiring
    add_line(model_name, 'Vortex_Force/1', 'Structural_Twin/1'); %Wire/line connections between blocks
    add_line(model_name, 'Structural_Twin/1', 'PZT_Buffer_Circuit/1');
    add_line(model_name, 'PZT_Buffer_Circuit/1', 'PZT_Voltage_V/1');
    
    save_system(model_name);
    
    fprintf('<strong>[SUCCESS] DTIL Model Generated with Block 4.5 Physics.</strong>\n');
    fprintf('   > Buffer Resistance: %.1f Mega-Ohms\n', R_buffer/1e6);
    fprintf('   > Sensor Capacitance: %.2f nF\n', C_pzt*1e9);
    fprintf('   > RC Time Constant: %.3f s (Cutoff: %.2f Hz)\n', Tau_c, 1/(2*pi*Tau_c));
end
