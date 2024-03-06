%-------------------------------------------------------------------------
%                      Thin_Plate_Solver.m
%                      January 2022
%
% This is the main script where the  boundary conditions for both the 3-
% noded and 6-noded elements are defined, and the functions
% getKUmatrices.m,evalStrain.m and plotDeformation.m are called to
% determine the displacement and strain in the thin plate.
%-------------------------------------------------------------------------
% Variable Name      Definition
%-------------------------Internal----------------------------------------
% E                  Modulus of Elasticity
% nu                 Poisson's Ratio
% thickness          Thickness of thin plate (Z)
% DoF                Total Degrees of Freedom of system
% R                  Radius of cut-out
% P                  Uniform tension applied to plate
% theta              Angle of load from vertical Y direction
%
%----------Three-Noded-Elements-------------------------------------------
%
% total_node3        Total number of nodes in 3-noded system
% n_node3            Number of nodes per element
% Nxy3               [Y,X] matrix of node positions
% H3                 Length of applied uniform tension
% load_nodes3        Array containing indices of nodes with applied load
% elems3             Matrix of node indices of each element
% load_matrix3       [X,Y] matrix of load at each node
% KE3                Global Element Stiffness Matrix
% U3                 Nodal Displacement Matrix
% Be3                Global Strain Matrix
% def_plot3          Axes object for deformation plot
% epsil3             Strain Matrix
%
%----------Six-Noded-Elements---------------------------------------------
%
% total_node6        Total number of nodes in 6-noded system
% n_node6            Number of nodes per element
% Nxy6               [Y,X] matrix of node positions
% H6                 Length of applied uniform tension
% load_nodes6        Array containing indices of nodes with applied load
% elems6             Matrix of node indices of each element
% load_matrix6       [X,Y] matrix of load at each node
% K6                 Global Element Stiffness Matrix
% U6                 Nodal Displacement Matrix
% Be6                Global Strain Matrix
% def_plot6          Axes object for deformation plot
% epsil6             Strain Matrix
% theta_var          Array containing theta test values
% E_var              Array containing E test values
% count              Count number of different combinations of theta and E
% U_var              3D matrix containing displacment values for all test
%                    values
% Be_var             4D matrix containing global strain matrices for all
%                    test values
% epsil_var          3D matrix containing strain matrices for all test
%                    values

% Clear data from previous run and change figure background colour to white
clear all; clc; close all;
set(0,'defaultfigurecolor',[1 1 1])

%% VARIABLES SHARED BETWEEN BOTH THE THREE AND SIX NODED TRIANGULAR ELEMENTS
global E 
E = 100e9;
global nu 
nu = 0.3;
global thickness 
thickness = 2e-3;
global DoF
DoF = 2;
R = 10;
global P
P = 200e6;
global theta
theta = 90;

%% THREE-NODED TRIANGULAR ELEMENTS

total_node3 = 5;
n_node3 = 3;

Nxy3 = [[0 0];               % Node 1
           [30 0];           % Node 2
           [30 (40 - R)];    % Node 3
           [(30 - R) 40];    % Node 4
           [0 40]]*1e-3;     % Node 5
   
H3 = Nxy3(4,1)-Nxy3(5,1);
load_nodes3 = [4 5];
elems3 = [1 3 2;
          1 4 3;
          1 5 4];

% Calculate element stiffness and strain matrices for each element 
load_matrix3 = zeros(total_node3,DoF);
[K3,U3,Be3] = getKUmatrices(Nxy3,elems3,n_node3,total_node3,load_matrix3,load_nodes3,H3);
epsil3 = evalStrain(elems3,Be3,U3);

% Plot overlayed deformation of thin plate
def_plot3 = plotDeformation(elems3,Nxy3,U3);


%% SIX-NODED TRIANGULAR ELEMENTS

total_node6 = 12;
n_node6 = 6;

Nxy6 =     [[0 0];                              % Node 1
           [30 0];                              % Node 2
           [30 (40 - R)];                       % Node 3
           [(30 - R) 40];                       % Node 4
           [0 40];                              % Node 5
           [15 0];                              % Node 6
           [30 ((40-R)/2)];                     % Node 7
           [25 35];                             % Node 8
           [(30-R)/2 40];                       % Node 9
           [0 20];                              % Node 10
           [15 (40-R)/2];                       % Node 11
           [(30-R)/2 20]]*1e-3;                 % Node 12
            
load_nodes6 = [4 5 9];
H6 = Nxy6(4,1)-Nxy6(5,1);
elems6 = [1 11 3 7 2 6;
          1 12 4 8 3 11;
          1 10 5 9 4 12];

% Calculate element stiffness and strain matrices for each element
load_matrix6 = zeros(total_node6,DoF);
[K6,U6,Be6] = getKUmatrices(Nxy6,elems6,n_node6,total_node6,load_matrix6,load_nodes6,H6);
epsil6 = evalStrain(elems6,Be6,U6);

% Plot overlayed deformation of thin plate
def_plot6 = plotDeformation(elems6,Nxy6,U6);


%% CHANGE LOAD AND E CONDITIONS FOR SIX-NODED TRIANGULAR ELEMENTS

% Initialise range of theta and E values
theta_var = [0 30 60 90 120];
E_var = [50 100 150 200 250]*1e9;
count = 1;

% Preallocate matrices for efficiency
U_var = zeros(total_node6,DoF,length(theta_var)*length(E_var));
Be_var = zeros(length(elems6(:,1)),total_node6,length(elems6(:,1)),length(theta_var)*length(E_var));
epsil_var = zeros(length(elems6(:,1)),length(elems6(:,1)),length(theta_var)*length(E_var));

% Loop through to get strain and displacement values for the range of theta
% and E values.
for i = 1:length(theta_var)
    theta = theta_var(i);
    for j = 1:length(E_var)
        E = E_var(j);
        disp([num2str(count),'      ', 'theta = ',num2str(theta),'     ', 'E =',num2str(E)])
        [~,U_var(:,:,count),Be_var(:,:,:,count)] = getKUmatrices(Nxy6,elems6,n_node6,total_node6,load_matrix6,load_nodes6,H6);
        epsil_var(:,:,count) = evalStrain(elems6,Be_var(:,:,:,count),U_var(:,:,count));
        count = count + 1;
    end
end

% Plot displacement and strain against E and theta for Node 4 and Element 3
figure
subplot(2,2,1)
% variable theta, constant E, displacement
plot(theta_var,reshape(U_var(4,1,(2:5:length(theta_var)*length(theta_var)),:),[],5),'k','linewidth',1.5);
hold on
plot(theta_var,reshape(U_var(4,2,(2:5:length(theta_var)*length(theta_var)),:),[],5),'r','linewidth',1.5);
xlabel('\theta (degrees)','Interpreter','Tex')
ylabel('Node 4 Deformation (meters)','Interpreter','Tex')
legend({'X','Y'})
xlim([min(theta_var) max(theta_var)])
subtitle('E = 100 Gpa')
% variable E, constant theta, displacement
subplot(2,2,2)
plot(E_var,reshape(U_var(4,1,(16:20)),[],5),'k','linewidth',1.5);
hold on
plot(E_var,reshape(U_var(4,2,(16:20)),[],5),'r','linewidth',1.5);
xlabel('E (Pa)','Interpreter','Tex')
ylabel('Node 4 Deformation (meters)','Interpreter','Tex')
legend({'X','Y'})
subtitle('\theta = 90^\circ','Interpreter','Tex')
xlim([min(E_var) max(E_var)])
% variable theta, constant E, strain
subplot(2,2,3)
plot(theta_var, reshape(epsil_var(3,1,(2:5:length(theta_var)*length(theta_var)),:),[],5),'k','linewidth',1.5)
hold on
plot(theta_var, reshape(epsil_var(3,2,(2:5:length(theta_var)*length(theta_var)),:),[],5),'r','linewidth',1.5)
hold on
plot(theta_var, reshape(epsil_var(3,3,(2:5:length(theta_var)*length(theta_var)),:),[],5),'color','#4DBEEE','linewidth',1.5)
xlabel('\theta (degrees)','Interpreter','Tex')
ylabel('Element 3 Strain','Interpreter','Tex')
legend({'e_x','e_y','e_{xy}'},'Interpreter','Tex')
xlim([min(theta_var) max(theta_var)])
subtitle('E = 100 Gpa')
% variable E, constant theta, strain
subplot(2,2,4)
plot(E_var,reshape(epsil_var(3,1,(16:20)),[],5),'k','linewidth',1.5)
hold on
plot(E_var,reshape(epsil_var(3,2,(16:20)),[],5),'r','linewidth',1.5)
hold on
plot(E_var,reshape(epsil_var(3,3,(16:20)),[],5),'color','#4DBEEE','linewidth',1.5)
xlabel('E (Pa)','Interpreter','Tex')
ylabel('Element 3 Strain','Interpreter','Tex')
legend({'e_x','e_y','e_{xy}'},'Interpreter','Tex')
xlim([min(E_var) max(E_var)])
subtitle('\theta = 90^\circ','Interpreter','Tex')