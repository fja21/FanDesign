clear
clc

%% Validation Case 1 - S1U

TestFan1(1) = 1.4;      % Volumetric Flow Rate 
TestFan1(2) = 0.495;    % Diameter of Fan 
TestFan1(3) = 0.2475;   % Diameter of Hub 
TestFan1(4) = 0.5;      % Diameter of Duct 
TestFan1(5) = 1500;     % Rotational Speed $Better match for Appendix Data Tables
TestFan1(6) = 9;        % Blade Number 
TestFan1(7) = 0.95;     % Volumetric Efficiency 
TestFan1(8) = 0.9;      % Hydraulic Efficiency
TestFan1(9) = 10;       % Number of Blade Sections 
TestFan1(10) = 0;       % Blade Loading Parameter 
TestFan1(11) = 293.15;  % Ambient Temperature  
TestFan1(12) = 101325;  % Ambient Pressure 
TestFan1(13) = 150;     % Total-to-Static Pressure Difference 

TestFan2 = 'NACA 4510';
%%%% Creates an empty array (alpha,gamma,mu)
TestFan3 = zeros(3,TestFan1(9));

TestFan3(1,:) = [8.0 7.4 6.8 6.3 5.9 5.5 5.1 4.7 4.3 4.0];

% Vortex Parameters
%%%% Creates an empty array of 5 coeff  (FV,CV1,CV2,CV3)
TestFan4 = zeros(1,5);
TestFan4(2) = 1.48;

CaseS1U = FanDesign(TestFan1,TestFan3,TestFan2,TestFan4);

%% %% Validation Case S2B

TestFan1(1) = 1.43;      % Volumetric Flow Rate 
TestFan1(2) = 0.495;    % Diameter of Fan 
TestFan1(3) = 0.2475;   % Diameter of Hub 
TestFan1(4) = 0.5;      % Diameter of Duct 
TestFan1(5) = 1500;     % Rotational Speed
TestFan1(6) = 9;        % Blade Number 
TestFan1(7) = 0.9542;     % Volumetric Efficiency 
TestFan1(8) = 0.9;      % Hydraulic Efficiency
TestFan1(9) = 10;       % Number of Blade Sections 
TestFan1(10) = 1;       % Blade Loading Parameter 
TestFan1(11) = 293.15;  % Ambient Temperature  
TestFan1(12) = 101325;  % Ambient Pressure 
TestFan1(13) = 150;     % Total-to-Static Pressure Difference 
TestFan3 = zeros(3,TestFan1(9));

TestFan3(1,:) = [7.1 6.6 6.2 5.8 5.5 5.1 4.8 4.5 4.3 4.0]; %Angle of Attack
TestFan3(2,:) = [0 -9 -16 -23 -29 -35 -40 -45 -50 -55]; % delta Sweep Angle 
TestFan3(3,:) = [0 -3 -5 -8 -10 -13 -15 -17 -20 -23]; % Dihedral Angle nu



TestFan4 = zeros(1,5);
TestFan4(1) = 3;    % a 
TestFan4(2) = 1.115; %1.115;  % b %1.48;


CaseS2B = FanDesign(TestFan1,TestFan3,TestFan2,TestFan4);

%% %% Validation Case 3 - S2F

TestFan1(1) = 1.4;      % Volumetric Flow Rate 
TestFan1(2) = 0.495;    % Diameter of Fan 
TestFan1(3) = 0.2475;   % Diameter of Hub 
TestFan1(4) = 0.5;      % Diameter of Duct 
TestFan1(5) = 1500;     % Rotational Speed $Better match for Appendix Data Tables
TestFan1(6) = 9;        % Blade Number 
TestFan1(7) = 0.95;     % Volumetric Efficiency 
TestFan1(8) = 0.9;      % Hydraulic Efficiency
TestFan1(9) = 10;       % Number of Blade Sections 
TestFan1(10) = 1;       % Blade Loading Parameter 
TestFan1(11) = 293.15;  % Ambient Temperature  
TestFan1(12) = 101325;  % Ambient Pressure 
TestFan1(13) = 140;     % Total-to-Static Pressure Difference 


TestFan3 = zeros(3,TestFan1(9));

TestFan3(1,:) = [7.1 6.6 6.2 5.8 5.5 5.1 4.8 4.5 4.3 4]; %Angle of Attack
TestFan3(2,:) = [0 8.5 16 22.8 29 34.8 40.3 45.4 50.3 55]; % delta Sweep Angle 
TestFan3(3,:) = [0 3.3 5.4 8.2 9.9 12.5 15 17.4 19.7 23]; % Dihedral Angle nu



TestFan4 = zeros(1,5);
TestFan4(1) = 3;    % a 
TestFan4(2) = 1.02; %1.115;  % b %1.48;

CaseS2F = FanDesign(TestFan1,TestFan3,TestFan2,TestFan4);

return
%%

r = [0.248 0.286 0.320 .350 0.378 0.404 0.429 0.452 0.474 0.495]/2;

Cm2 = [8.727 8.984 9.317 9.670 10.023 10.366 10.697 11.014 11.317 11.606];
Cm2 = [8.727 8.984 9.317 9.670 10.023 10.366 10.697 11.014 11.317 11.606]

V_dpre = trapz(r,Cm2.*r)*2*pi

cu2 = [14.313 12.797 11.763 11.000 10.406 9.928 9.532 9.193 8.908 8.657]
b = (r.*cu2)./(3*r)
cu22 = (3.*r + 1.02)./r