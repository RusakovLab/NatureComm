% close all;
% clear;

TimeGlu = 10; % ms 
BinGlu=4;
timeCalculus = 200; %ms
NPar = 1000;
DistanceReal=(1:1:BinGlu+1)*(1/BinGlu);
InitialDistribution = zeros(NPar, BinGlu);

myFolder = pwd; % or 'C:\Users\yourUserName\Documents\My Pictures' or whatever...
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, 'PD 3 *.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
NumberOfFiles = length(theFiles);
TotalSumMatrix = zeros(NumberOfFiles*NPar, BinGlu);



for K = 1 : NumberOfFiles
%for K = 1 : 3
    baseFileName = theFiles(K).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    thisStructure = load(fullFileName);
    output =thisStructure;
    %InitialDistribution=InitialDistribution + output;
    InitialDistribution=output;


    %***********************************************
    TotalSumMatrix(1+(K-1)*NPar:K*NPar, :) = output;
    %***********************************************
end

%   Free = 1 Bound = 0
Free = 1 ;
Bound = 0;
%***********************
StateParticle = Free;

selected_rowsTotal = TotalSumMatrix(TotalSumMatrix(:, 4) == StateParticle, :);


%****Plot********************************************************
X1 = selected_rowsTotal(:, 1 ); % X coordinates (random values)
Y1 = selected_rowsTotal(:, 2 ); % Y coordinates (random values)
Z1 = selected_rowsTotal(:, 3 ); % Z coordinates (random values)
selected_rowsTotalPrint = TotalSumMatrix(TotalSumMatrix(:, 4) == StateParticle, :);
selected_rowsTotalPrint(:,4) =  sqrt(X1.^2 + Y1.^2 +Z1.^2);
FinalMatrix = selected_rowsTotalPrint(selected_rowsTotalPrint(:, 4) < 3100, :);
% Plot the data as a translucent cloud
figure(1);
scatter3(X1, Y1, Z1, 'filled', 'MarkerFaceAlpha',0.5, 'MarkerEdgeAlpha', 1);
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Data Visualization');
% Adjust plot settings
grid on;
axis equal;

%************************************************************

%InitialDistribution=InitialDistribution/length(theFiles);
% Input data
% Select rows with only 1 in the fourth column
selected_rows = InitialDistribution(InitialDistribution(:, 4) == StateParticle, :);
%selected_rows = InitialDistribution;

X = selected_rows(:, 1 ); % X coordinates (random values)
Y = selected_rows(:, 2 ); % Y coordinates (random values)
Z = selected_rows(:, 3 ); % Z coordinates (random values)

% Plot the data as a translucent cloud
figure(2);
scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);
%scatter3(X, Y, Z, 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Data Visualization');


% Adjust plot settings
grid on;
axis equal;
%  xlim([-2000, 2000])
%  ylim([-2000, 2000])
%  zlim([-2000, 2000])

% Optionally, you can add more visualization options, such as changing the colors or size of markers.
% For example, to change marker size, add the 'SizeData' parameter to the scatter3 function:
% scatter3(X, Y, Z, 'filled', 'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2, 'SizeData', 50);



