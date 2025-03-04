function [BasicParameters, Balls, ComputationSet, Cleft, SizeCube, Adhesion]=InputParametersSR()

fid = fopen('statisticSR.txt');
s   = textscan(fid,'%s%f', 'CommentStyle','%','Delimiter','=');
fclose(fid);


% Read lines of input parameters from the file statistic.txt
idxnp = strncmpi('np',s{1},2);
idxTrials = strncmpi('Trials',s{1},2);
indMaxProbAdhesive = strncmpi('MaxProbAdhesive',s{1},2);
indTimeINsideAdhesiveZone =strncmpi('TimeINsideAdhesiveZone',s{1},2);
indProbabilityofAstrocytes = strncmpi('ProbabilityofAstrocytes',s{1},2);
indNumbersphere = strncmpi('Numbersphere',s{1},2);
indUnboundProb  = strncmpi('UnboundProb',s{1},2);

% convert the parameters to the numerical value
np = s{2}(idxnp);                  % Number of particles
Trials =  s{2}(idxTrials);
MaxProbAdhesive = s{2}(indMaxProbAdhesive);
TimeINsideAdhesiveZone  = s{2}(indTimeINsideAdhesiveZone);
ProbabilityofAstrocytes  = s{2}(indProbabilityofAstrocytes);
Numbersphere  = s{2}(indNumbersphere);
UnboundProb = s{2}(indUnboundProb);

% ************************************************************************
SizeCube.Size = 4000; %3400
% Size of visualization area***********************************************
SizeCube.AxisSize=SizeCube.Size/2;
%***************************************************************************
% Basic parameters
% Save files 
BasicParameters.SaveFiles = 'y';  % y is save files

BasicParameters.NumberOfParticles = np;   % Number of p

BasicParameters.charge = 0;                % Charge of particles
BasicParameters.TotalTime = 8     ;       % Time of computation in ms 0.21
BasicParameters.Numbersphere = Numbersphere; 
BasicParameters.UnboundProb = UnboundProb;
% 0.2	1400 26 January 2023
% 0.3	1120
% 0.4	850
% 0.5	600
BasicParameters.BeginTime = 0.8;            % Offset  to compute the diffusion coefficient in ms 0.015
BasicParameters.EndTime = 0.9;              % End time to compute the diffusion coefficient in ms 0.02
BasicParameters.Trials = Trials;                 % Number of Runs to compute the statistics
BasicParameters.Control = 2;                % if Control ==1 the algorithm save into the file the final distribution of particles. At the end


% Parameters for Balls
Balls.name = 'Balls';
Balls.MaxBallRadius = 350.00; % Maximum radius 350
Balls.MinBallRadius = 50.0;   % Minimum radius
Balls.BettaControl = 0.8;     % Betta Control - normaly it is switch off
Balls.OverlapingBalls = 'y';  % y is overlapping Balls
Balls.ProbabilityofAstrocytes=ProbabilityofAstrocytes; % 0.1 proporsion of astro

% *************************************************************************
% Computational process parameters
ComputationSet.DifCoefINIT=(0.4*10^6); %diffusion coefficient um^2/sec % 0.08 alpfa 0.2, 0.04 Alpfa 0.1
ComputationSet.GTimeStep=0.0001; % step of Brownian motion in ms 6.25e-05   0.00025*ScalingFactor/4;
ComputationSet.StepRad=sqrt(6*ComputationSet.DifCoefINIT * ComputationSet.GTimeStep); % step of computation in um
ComputationSet.StepLocalMovement=0.1; % zero step-distance  fluctuation in proportion to StepRad


% *************************************************************************
% Geometrical parameters of the computational domain
% Cleft parameters
Cleft.RadiusofCleft= 0.00000100;
Cleft.CleftWidth=0.000000010;
Cleft.DistanceAstrocyte = 10; % Distance location of astrocytes centers from the cleft

%**************************************************************************
% Adhesion properties of astocytes
Adhesion.TimeOfDeAdhesion = 4; % time of detached from ball by ms
%Adhesion.MaxProbAdhesive =0.3; %1 = instant adhesion 0 = non-adhesion
Adhesion.MaxProbAdhesive =MaxProbAdhesive; %1 = instant adhesion 0 = non-adhesion
Adhesion.TimeINsideAdhesiveZone =TimeINsideAdhesiveZone;
%Adhesion.TimeINsideAdhesiveZone =0.001; % Number of steps to reach haft of maximum adgesen is TimeINsideAdhesiveZone/StepRad
                                        % for 0.001 => 160 steps
                                        % ms the constant of the time spent 
                                        % in the zone of the transporters necessary for the glutamate to bind to the transporters. 
                                        % The shorter the time, the better the connection.

% ProbabilityToStuck=0.4;
% DistanceofAdhesion = ProbabilityToStuck*160*ComputationSet.StepRad; % um
% Adhesion.RandomNumber =1 - DistanceofAdhesion;



end