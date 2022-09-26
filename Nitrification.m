%% Initiation
clear       % Clear workspace
close all   % Close open figure windows
clc         % Clear command window

%% Structure of the code
% Define model parameters
% Define initial conditions
% Define ODE functions and equations

%% Define model parameters 
% Parameters are from "Sin2008". Stated below are the parameters for the 
% denitrification steps ammonium-oxidizing bacteria (AOO), 
% nitrite-oxidizing bacteria (NOO) and dissolved oxygen

% Parameter name = "Value" %Short name , range, [unit]

load dataAOONOO 
ts = td;

% AOO
Yield_AOO                   = 0.15; % Y_AOO,      range: 0.11 - 0.21, [ mgCOD/mgN ]            
MaximumGrowthRate_AOO       = 1.5;  % mu_max,AOO, range: 0.5 - 2.1,   [ d-1]   
SubstrateAffinity_AOO       = 0.5;  % ks_AOO,     range: 0.14 - 1,    [ mgNH4-N/L ]
OxygenAffinity_AOO          = 0.7;  % ko_AOO,     range: 0.1 - 1.45,  [ mgO/L ]
DecayRateCoefficient_AOO    = 0.1;  % b_AOO,      range: 0.07 - 0.3,  [ d-1 ]

% NOO
Yield_NOO                   = 0.05; % Y_NOO,      range: 0.03 - 0.09, [ mgCOD/mgN ]       
MaximumGrowthRate_NOO       = 0.5;  % mu_max,NOO, range: 0.4 - 1.05,  [ d-1]
SubstrateAffinity_NOO       = 1.5;  % ks_NOO,     range: 0.1 - 3,     [ mgNO2-N/L ]
OxygenAffinity_NOO          = 1.5; % ko_NOO,     range: 0.3 - 1.5,   [ mgO/L ]
DecayRateCoefficient_NOO    = 0.12; % b_NOO,      range: 0.08 - 0.2,  [ d-1 ]

% Oxygen - assumed known, not estimated here. No range given.
OxygenMassTransfer          = 360;  % kla,        range: 360          [ d-1 ]
OxygenSaturation            = 8;    % o2sat,      range: 8   at 25 C  [ mg O2/L]  

% Parameters collected in vector
Parameters = [Yield_AOO, MaximumGrowthRate_AOO, SubstrateAffinity_AOO, OxygenAffinity_AOO, DecayRateCoefficient_AOO, Yield_NOO, MaximumGrowthRate_NOO, SubstrateAffinity_NOO, OxygenAffinity_NOO, DecayRateCoefficient_NOO, OxygenMassTransfer, OxygenSaturation];

%% Define initial conditions for the experiments
InitialNH4            = 20;  % Initial NH4 concentration [ mgN/L ] 
InitialNO2            = 0;   % Initial NO2 concentration [ mgN/L ]
InitialNO3            = 0;   % Initial NO3 concentration [ mgN/L ] 
InitialOxygen         = 8;   % Initial oxygen concentration [ mgO2/L ]
InitialCOD_AOO        = 75;  % Initial COD concentration in AOO [ mgCOD/L ] 
InitialCOD_NOO        = 25;   % Initial COD concentration in NOO [ mgCOD/L ]

InitialConditions = [...
    InitialNH4;...
    InitialNO2;...
    InitialNO3;...
    InitialOxygen;...
    InitialCOD_AOO;...
    InitialCOD_NOO];

%% Define inputs conditions for the experiments
InputNH4            = 6;  % Input NH4 concentration [ mgN/L ] 
InputNO2            = 0;   % Input NO2 concentration [ mgN/L ]
InputNO3            = 0;   % Input NO3 concentration [ mgN/L ] 
InputOxygen         = 10;   % Input oxygen concentration [ mgO2/L ]
InputCOD_AOO        = 0;  % Input COD concentration in AOO [ mgCOD/L ] 
InputCOD_NOO        = 0;  % Input COD concentration in NOO [ mgCOD/L ]
InputD              = 0; % Input flow rate

Inputs = [...
    InputNH4,...
    InputNO2,...
    InputNO3,...
    InputOxygen,...
    InputCOD_AOO,...
    InputCOD_NOO,...
    InputD];

%% Solve the model
options = odeset('RelTol',1e-7,'AbsTol',1e-8);
[Time1,ModelOutput] = ode45(@NitrificationModel,ts,InitialConditions,options,Parameters,Inputs);

%% Plot the results
% figure(1)
Time=Time1*24; % from d to h.
plot(td*24,nh4,'ko',Time,ModelOutput(:,1),'r','linewidth',1.4)
hold on
plot(td*24,no2,'k*',Time,ModelOutput(:,2),'g','linewidth',1.4)
plot(td*24,no3,'k+',Time,ModelOutput(:,3),'b','linewidth',1.4)
legend('NH4-lab','NH4-sim','NO2-lab','NO2-sim','NO3-lab','NO3-sim')
ylabel('$mg$ N $L^{-1}$','Interpreter','Latex','fontsize',12,'FontWeight','bold')
legend('boxon')
grid on
set(legend,'Interpreter','latex')
xlabel('Time ($h$)','Interpreter','Latex','fontsize',12,'FontWeight','bold')
xlim([0,4])
ylim([-0.1,20.1])
figure(2)
plot(td*24,do,'k*',Time,ModelOutput(:,4),'r','linewidth',1.4)
ylabel('DO','Interpreter','Latex','fontsize',12,'FontWeight','bold') 
legend('DO-lab','DO-sim')
set(legend,'Interpreter','latex')
xlabel('Time ($h$)','Interpreter','Latex','fontsize',12,'FontWeight','bold')
%xlim([0,4])
grid on

%% Define ODE function
function dxdt = NitrificationModel(t,x,Parameters,inputs)

YAOB        = Parameters(1);
mumaxAOB    = Parameters(2);
ks_AOB      = Parameters(3);
ko_AOB      = Parameters(4);
bAOB        = Parameters(5);

YNOB        = Parameters(6);
mumaxNOB    = Parameters(7);
ks_NOB      = Parameters(8);
ko_NOB      = Parameters(9);
bNOB        = Parameters(10);

kla         = Parameters(11);
o2sat       = Parameters(12);

%% ODE model equations
% step 1: NH4 + 3 O_2 --> 2 NO2 + 4 H + 2 H2O    by AOB: Ammonia Oxidizing Bacteria
% step 2: 2 NO2 + O_2  --> 2 NO3                 by NOB: Nitrite Oxidizing Bacteria

% total oxidation reaction: NH4 + 2 O_2 --> NO3 _ 2 H + H2O

% Definition of state variables: x(1): NH4 ; x(2): NO2;  x(3): NO3; x(4):O2; x(5):XAOB; x(6):XNOB.
x1_in = inputs(1);
x2_in = inputs(2);
x3_in = inputs(3);
x4_in = inputs(4);
x5_in = inputs(5);
x6_in = inputs(6);
D     = inputs(7);

% Define growth rate for AOB & NOB respectively
muAOB = mumaxAOB*(x(1)/(x(1)+ks_AOB))*(x(4)/(x(4)+ko_AOB))*x(5); 
muNOB = mumaxNOB*(x(2)/(x(2)+ks_NOB))*(x(4)/(x(4)+ko_NOB))*x(6);

% Define mass balances for batch reactor with pulse addition of substrates
dxdt(1,1)= D*(x1_in - x(1)) + (-1/YAOB)* muAOB ;                           % NH4 balance
dxdt(2,1)= D*(x2_in - x(2)) + (1/YAOB)* muAOB - (1/YNOB)* muNOB ;          % NO2 balance 
dxdt(3,1)= D*(x3_in - x(3)) + (1/YNOB)* muNOB ;                            % NO3 balance 
dxdt(4,1)= D*(x4_in - x(4)) + (1 - 3.43/YAOB) * muAOB - bAOB*x(5) + (1 - 1.14/YNOB) * muNOB - bNOB*x(6)+ kla * (o2sat-x(4)); % oxygen balance 
dxdt(5,1)= D*(x5_in - x(5)) + muAOB - bAOB*x(5);                           % AOB biomass
dxdt(6,1)= D*(x6_in - x(6)) + muNOB - bNOB*x(6);                           % NOB biomass

end