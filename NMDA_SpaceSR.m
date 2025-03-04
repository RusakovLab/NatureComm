%Glu = 0.01; % mM

TimeGlu = 8; % ms 
BinGlu=100;
timeCalculus = 200; %ms
DistanceReal=(1:1:BinGlu+1)*(1/BinGlu);
InitialDistribution = zeros(BinGlu, 101);

myFolder = pwd; % or 'C:\Users\yourUserName\Documents\My Pictures' or whatever...
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, 'DistanceFree*.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for K = 1 : length(theFiles)
    baseFileName = theFiles(K).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    thisStructure = load(fullFileName);
    output =thisStructure;
    InitialDistribution=InitialDistribution+output;
end
InitialDistribution=InitialDistribution/length(theFiles);

%NormParameter = 20*4*10^19/(6.023*10^23);
NormParameter = 10^21/(6.023*10^23);
Volume = ((4/3) * pi * ((2 * DistanceReal(2:end)').^3 - (2 * DistanceReal(1:end-1)').^3));
%InitialDistribution(:, 2:end) = NormParameter * InitialDistribution(:, 2:end) ./ Volume;


ConcentrationOneMolecules = 1.66 * 10^-6; % mM 1 molecules creat such concentration in mM 
Con= (1000)*ConcentrationOneMolecules; % uM
InitialDistribution(:, 2:end) = Con*InitialDistribution(:, 2:end) ./ Volume;


figure(1)

contourf(log(InitialDistribution(:, 2:end)),10, 'LineColor','none')
colorbar
clim([-8,4.5]);
%*********************Differencial NMDA ***********************************


BinTotal = round(BinGlu*timeCalculus/TimeGlu); 

NMDATimeBin = linspace(0, timeCalculus, BinTotal); % Generate t for f 

GluTotal=zeros(1,BinTotal);
TwoBount=zeros(100, 100);
TwoBountRaw=1:1:100;
% Maintain original solver but enhance numerical handling
for Jitter = 1:BinGlu
    TestGlu = InitialDistribution(Jitter, 2:end);
    TestGlu(numel(GluTotal)) = 0;
    
    % Customized solver settings
    options = odeset('RelTol', 1e-12, ...       % Balanced precision
                    'AbsTol', 1e-20, ...        % Match your concentration scale
                    'NonNegative', [1 2 3 4 5], ... % Physical constraint
                    'Refine', 10, ...           % Increase output points
                    'MaxStep', 0.1);            % Prevent over-aggressive stepping
    
    % Original-scale computation
    ode = @(t,y) NMDA(t,y,TestGlu',NMDATimeBin, timeCalculus);
    [T,Y] = ode45(ode, [0 timeCalculus], [1 0 0 0 0], options);
    
    % Enhanced interpolation with precision preservation
    TempY = Y(:,4);
    
    % Safe interpolation for vanishing values
    valid_mask = TempY > eps;  % Identify non-zero values
    if any(valid_mask)
        TempY_interp = interp1(T(valid_mask), log10(TempY(valid_mask)), TwoBountRaw, 'pchip', -inf);
        TempY = 10.^TempY_interp;
    else
        TempY = zeros(size(TwoBountRaw));
    end
    
    % Apply physical floor
    TempY(TempY < realmin('double')) = 0;
    
    TwoBount(Jitter, :) = TempY(1:100);
end

axis manual
figure(2);                  % display image
%axis([0 timeCalculus 0 1])

subplot(2, 4, 1); plot(T,Y(:,1))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 2); plot(T,Y(:,2))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 3); plot(T,Y(:,3))
%axis([0 timeCalculus 0 1])

subplot(2, 4, 4); plot(T,Y(:,4))
%axis([0 timeCalculus 0 0.02])

subplot(2, 4, 5); plot(T,Y(:,5))
%axis([0 timeCalculus 0 0.01])
MinPlot=min(min(TwoBount(:,1:end)));
MaxPlot=max(max(TwoBount(:,1:end)));
figure(3)
%TwoBount(1,15)=0.49; SCALING MAX
contourf((TwoBount(:,1:end)), 20, 'LineColor','none')
colorbar; 
clim([0,0.4]);
%caxis([MinPlot MaxPlot]);





