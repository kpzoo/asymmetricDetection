% Reproduce Fig 1 on threshold based action delays
clearvars; clc; close all; 

% Assumptions and notes
% - single true COVID type incidence curve, no R estimation
% - view M curves with either under-reporting, delays or both
% - find time for crossing thresholds for decision making
% - choose under-reporting and delays for COVID-19

% Save data and directories for main code
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(100);

%% Setup single epidemic true simulation

% Choose a scenario and serial interval
epiNo = 4; scenNo = 4;
% Number of replicates and to plot
M = 1000; Mcol = 100;

% Initialise epidemic time and changepoint
tday0 = 1:121; nday0 = length(tday0); chgpt = 50;

% Define possible scenarios for true R and serial interval
scenNam = {'constant', 'cyclic', 'logistic', 'switch', 'boom-bust', 'bottle', '2-step', 'filtered'};
scenChoice = scenNam{scenNo}; disp(['True R scenario: ' scenChoice]);

% Define all SI/generation time distributions
epiNam = {'Marburg', 'MERS', 'Measles', 'COVID-19', 'EVD'};
distChoice = epiNam{epiNo}; disp(['Serial interval: ' distChoice]);

% Simulate epidemic scenarios and truncate initial 0s
Iwarn = 1; % ensure no warnings
while Iwarn
    [Iday, Lam, Rtrue, tday, Iwarn, Pomega0] = epiSimDiseaseChg(scenNo, epiNo, tday0, nday0, 1, chgpt);
end
if Iwarn
    warning('Sequences of zero incidence');
end
% Total number of days and cases
nday = length(tday); totcase = sum(Iday);
% Restrict Pomega
Pomega = Pomega0(1:nday);

% Times of infections of each cases (non-delayed)
tInf = zeros(1, totcase); Icheck = zeros(1, nday);
% Start and end indices for each day
caseStart = Icheck; caseEnd = Icheck;
for i = 1:nday
    % Starting case id for day i
    if i == 1
        caseStart(i) = 1;
    else
        caseStart(i) = caseEnd(i - 1) + 1;
    end
    % Ending case id for day i
    caseEnd(i) = caseStart(i) + Iday(i) - 1;
    
    % All the people infected on day i
    tInf(caseStart(i):caseEnd(i)) = i;
    
    % Check case counts
    Icheck(i) = length(find(tInf == i));
end
% Ensure breakdown is correct by reconstructing Iday
if ~all(Iday == Icheck)
    error('Assignment of case ids incorrect');
end

%% Generate M delayed versions of Iday

% Mean of delay distribution
mtau = 10.8; 
% NegBin delay parameters
r = 10; p = mtau/(r + mtau);

% Draws from this distribution
delay = nbinrnd(r, 1-p, [M, totcase]);
[~, vtau] = nbinstat(r, 1-p);

% PDF of delay distribution 
xdel = 0:50; Pdel = nbinpdf(xdel, r, 1-p);

% Delayed infection ids and days
tDel = zeros(M, totcase); ndaydel = zeros(1, M);
% Delayed and truncated incidence curves
Idel = cell(1, M); Itrunc = Idel;

% For every replicate add delay on infection ids
for i = 1:M
    % Delayed case ids
    tDel(i, :) = tInf + delay(i, :);
    
    % Construct an epi-curve for each delayed version
    ndaydel(i) = max(tDel(i, :));
    % Ensure maximum not < nday (e.g. due to Iday(end) = 0)
    if ndaydel(i) < nday
        ndaydel(i) = nday;
    end
    
    Itemp = zeros(1, ndaydel(i));
    for j = 1:ndaydel(i)
        % Delayed epi-curve
        Itemp(j) = length(find(tDel(i, :) == j));
    end
    % Check total cases conserved
    if sum(Itemp) ~= totcase
        error('Cases not conserved');
    else
        % Store incidence and truncated incidence
        Idel{i} = Itemp; Itrunc{i} = Itemp(1:nday);
    end
end

% Truncated incidence as a matrix
Itrunc = cell2mat(Itrunc');

%% Generate M under-reported versions of Iday

% Mean of sampling distribution
rho = 0.38; b = 20;
% Parameters of beta distribution
fr = rho/(1 - rho); a = fr*b; 

% Under-reporting distribution
xrep = 0:0.01:1; yrep = betapdf(xrep, a, b);

% Under-reported incidence curves
Isamp = zeros(size(Itrunc));
for i = 1:M
    for j = 1:nday
        % Main downsampling
        Isamp(i, j) = binornd(Iday(j), betarnd(a, b));
    end
end

%% Generate M under-reported and delayed versions of Iday

% Apply under-reporting to the already delayed curves
Icomb = zeros(size(Itrunc));
for i = 1:M
    for j = 1:nday
        % Main downsampling
        Icomb(i, j) = binornd(Itrunc(i, j), betarnd(a, b));
    end
end

% Find first time cross a threshold up and down
thresh = [20 40 60 80]; nthresh = length(thresh);
% Crossing times of true infection curve
tup0 = zeros(1, nthresh); tdown0 = tup0;
for i = 1:nthresh
    % Upward change
    tup0(i) = find(Iday >= thresh(i), 1, "first");
    % Downward change
    tdowns = find(Iday <= thresh(i)); tdowns = tdowns - chgpt;
    tdowns = tdowns(tdowns > 0); tdown0(i) = tdowns(1) + chgpt;
end

% Crossings under each type of noise and then both
tupSamp = zeros(M, nthresh); tdownSamp = tupSamp;
tupTrunc = zeros(M, nthresh); tdownTrunc = tupTrunc;
tupComb = zeros(M, nthresh); tdownComb = tupComb;

% For each noise type get upward and downward crossing per trajectory
for j = 1:M
    for i = 1:nthresh
        % Upward changes
        tupSamp(j, i) = find(Isamp(j, :) >= thresh(i), 1, "first");
        tupTrunc(j, i) = find(Itrunc(j, :) >= thresh(i), 1, "first");
        tupComb(j, i) = find(Icomb(j, :) >= thresh(i), 1, "first");
        % Downward changes
        tdowns = find(Isamp(j, :) <= thresh(i)); tdowns = tdowns - chgpt;
        tdowns = tdowns(tdowns > 0); tdownSamp(j, i) = tdowns(1) + chgpt;
        tdowns = find(Itrunc(j, :) <= thresh(i)); tdowns = tdowns - chgpt;
        tdowns = tdowns(tdowns > 0); tdownTrunc(j, i) = tdowns(1) + chgpt;
        tdowns = find(Icomb(j, :) <= thresh(i)); tdowns = tdowns - chgpt;
        tdowns = tdowns(tdowns > 0); tdownComb(j, i) = tdowns(1) + chgpt;
    end
end

% Mean errors in crossing time
eup = tupComb-tup0; edown = tdown0-tdownComb;


%% Visualise all noise and their distributions

% Plot the crossing times for different noise types
figure('Position', [10 10 800 800]);
subplot(1, 2, 1); hold on;
boxchart(tupComb); boxchart(tupTrunc); boxchart(tupSamp);
h = gca; h.XTickLabel = thresh;
stairs(1:4, tup0, 'k', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('threshold, $a$', 'FontSize', fnt); 
ylabel('upward crossing time', 'FontSize', fnt); 
subplot(1, 2, 2); hold on;
boxchart(tdownComb); boxchart(tdownTrunc); boxchart(tdownSamp);
h = gca; h.XTickLabel = thresh;
stairs(1:4, tdown0, 'k', 'LineWidth', 2);
xlabel('threshold, $a$', 'FontSize', fnt); 
ylabel('downward crossing time', 'FontSize', fnt); 
h = legend('combined', 'delayed reports', 'under-reporting', 'perfect');
h.Box = 'off'; hold off; grid off; box off;


figure('Position', [10 10 1000 1000]);
% Under-reported curves
subplot(3, 2, 1); hold on;
for i = 1:Mcol
    stairs(tday, Isamp(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(tday, Iday, 'k', 'LineWidth', 2);
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('time $t$', 'FontSize', fnt);
ylabel('under-reported $C_t$', 'FontSize', fnt);

% Distribution for sampling reporting probabilities
subplot(3, 2, 3);
plot(xrep, yrep/trapz(yrep), 'b', 'LineWidth', 2);
hold on; ax = gca; ax.YGrid = 'on'; 
plot([rho rho], ax.YLim, 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel(['$\rho_t$ $|$ $\bar{\rho}$ = ' num2str(rho)], 'FontSize', fnt);
ylabel('reporting P($\rho_t$)', 'FontSize', fnt);

% Delayed curves
subplot(3, 2, 2); hold on;
for i = 1:Mcol
    stairs(tday, Itrunc(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(tday, Iday, 'k', 'LineWidth', 2);
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('time $t$', 'FontSize', fnt);
ylabel('delayed $C_t$', 'FontSize', fnt);

% Distribution for samplibg delay probabilities
subplot(3, 2, 4);
plot(xdel, Pdel/trapz(Pdel), 'b', 'LineWidth', 2);
hold on; ax = gca; ax.YGrid = 'on'; 
plot([mtau mtau], ax.YLim, 'k--', 'LineWidth', 2);
grid off; box off; hold off;
xlabel(['$\delta_{x}$ $|$ $\bar{\delta}$ = ' num2str(mtau)], 'FontSize', fnt);
ylabel('delay P($\delta_x$)', 'FontSize', fnt);

% Delayed and under-reported curves
subplot(3, 2, [5 6]); hold on;
for i = 1:Mcol
    stairs(tday, Icomb(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(tday, Iday, 'k', 'LineWidth', 2);
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('time $t$', 'FontSize', fnt);
ylabel('combined noisy $C_t$', 'FontSize', fnt);
% Additional axes for R
axes('Position',[0.65 0.18 0.2 0.1]);
plot(tday, Rtrue, 'Color', 'r', 'LineWidth', 2);
hold on; plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2)
grid off; box off; hold off;
xlim([tday(1) tday(end)]);
xlabel('$t$', 'FontSize', fnt);
title('R numbers $R_t$', 'FontSize', fnt);



%% Publishable figure

figure('Position', [10 10 800 800]);
for i = 1:nthresh
    subplot(nthresh, 2, 2*i-1);
    h = histogram(tupComb(:, i)-tup0(i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.3; h.EdgeAlpha = 0;
    hold on;
    h = histogram(tupSamp(:, i)-tup0(i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.2; h.EdgeAlpha = 0;
    h = histogram(tupTrunc(:, i)-tup0(i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.2; h.EdgeAlpha = 0;
    hold off; box off; grid off;
    ylabel(['$a = $ ' num2str(thresh(i))], 'FontSize', fnt);
    if i == nthresh
        xlabel('$\Delta_{grow}$', 'FontSize', fnt);
        legend('combined', 'under-reporting', 'delayed reports');
    end

    subplot(nthresh, 2, 2*i);
    h = histogram(tdown0(i)-tdownComb(:, i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.3; h.EdgeAlpha = 0;
    hold on;
    h = histogram(tdown0(i)-tdownSamp(:, i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.2; h.EdgeAlpha = 0;
    h = histogram(tdown0(i)-tdownTrunc(:, i), 'Normalization','probability', 'BinWidth', 1);
    h.FaceAlpha = 0.2; h.EdgeAlpha = 0;
    hold off; box off; grid off;
    ylabel(['$a = $ ' num2str(thresh(i))], 'FontSize', fnt);
    if i == nthresh
        xlabel('$\Delta_{wane}$', 'FontSize', fnt);
    end
end




