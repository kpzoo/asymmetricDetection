% Reproduce Fig 2 on R credible interval actions
clearvars; clc; close all; 

% Assumptions and notes
% - follows simulations from Fig 2
% - compute FI ratio across time cumulatively

% Save data and directories of code for plotting options
thisDir = cd; cd('main'); mainDir = cd;
cd(thisDir); addpath(genpath(mainDir));

% Default plotting options
[grey1, grey2, cmap, fnt] = defaultSet(100);

%% Setup single epidemic true simulation

% Choose a scenario and serial interval (need constant R)
epiNo = 4; scenNo = 4;
% Number of replicates and to plot
M = 1000; Mcol = 100;

% Initialise epidemic time and changepoint
tday0 = 1:101; nday0 = length(tday0); chgpt = 50;

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
    if max(Iday) < 2000
        Iwarn = 1;
    end
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

%% Proper computation of FI with intervals

% Upward FI goes until chgpt, then downard FI
nup = length(1:chgpt); ndown = length(chgpt+1:nday);
FIup = zeros(1, nup); FIdown = zeros(1, ndown);
% Noisy versions of FI
FIup_noise = zeros(M, nup); FIdown_noise = zeros(M, ndown);

% Compute upward FI treating each point as the present
for i = 1:nup
    % Present time and FI of infection
    T = i; FIup(i) = sum(Lam(1:T)); 
    
    % Delay cumulative probabilities t-s direction
    F_ts = nbincdf(0:T-1, r, 1-p); F_ts = F_ts(end:-1:1);
    % Sampling probabilities
    rho_s = betarnd(a, b, [M T]);

    % FI of cases (noisy)
    for j = 1:M
        FIup_noise(j, i) = sum(rho_s(j, :).*F_ts.*Lam(1:T));
    end
end


% Compute downward FI treating each point as the present
for i = 1:ndown
    % Present time and FI of infection
    T = i; FIdown(i) = sum(Lam(chgpt+1:chgpt+T)); 
    
    % Delay cumulative probabilities t-s direction
    F_ts = nbincdf(0:T-1, r, 1-p); F_ts = F_ts(end:-1:1);
    % Sampling probabilities
    rho_s = betarnd(a, b, [M T]);

    % FI of cases (noisy)
    for j = 1:M
        FIdown_noise(j, i) = sum(rho_s(j, :).*F_ts.*Lam(chgpt+1:chgpt+T));
    end
end

% Take means across the sample probabilities
FIup_noise_m = mean(FIup_noise); FIdown_noise_m = mean(FIdown_noise);
% Find ratio to the perfect surveillance cases
FIup_ratio = FIup_noise./FIup; FIdown_ratio = FIdown_noise./FIdown;

% Quantiles of FI ratios and FI
FIdown_ratio_m = mean(FIdown_ratio); FIup_ratio_m = mean(FIup_ratio);
FIdown_ratio_q = quantile(FIdown_ratio, [0 1]);
FIup_ratio_q = quantile(FIup_ratio, [0 1]);
FIdown_noise_q = quantile(FIdown_noise, [0 1]);
FIup_noise_q = quantile(FIup_noise, [0 1]);

%% Visualise noise and FI ratios

% Perfect and noisy FI in growth
figure('Position', [10 10 600 600]);
subplot(2, 2, 1); hold on;
plot(tday(1:chgpt), FIup_noise, 'LineWidth', 2);
plot(tday(1:chgpt), FIup, 'k', 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$\mathbf{I}(R|C_1^t)$', 'FontSize', fnt);

% Perfect and noisy FI in decline
subplot(2, 2, 2); hold on;
plot(tday(chgpt+1:end), FIdown_noise, 'LineWidth', 2);
plot(tday(chgpt+1:end), FIdown, 'k', 'LineWidth', 2)
hold off; grid off; box off;
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$\mathbf{I}(R|C_1^t)$', 'FontSize', fnt);

% Ratios of the FI streams
subplot(2, 2, 3); hold on;
plot(tday(1:chgpt), FIup_noise./FIup, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$\frac{\mathbf{I}(R|C_1^t)}{\mathbf{I}(R|I_1^t)}$', 'FontSize', fnt);
subplot(2, 2, 4); hold on;
plot(tday(chgpt+1:end), FIdown_noise./FIdown, 'LineWidth', 2);
hold off; grid off; box off;
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$\frac{\mathbf{I}(R|C_1^t)}{\mathbf{I}(R|I_1^t)}$', 'FontSize', fnt);

%% Publishable figure


figure('Position', [10 10 800 800]);
subplot(2, 2, [1 2]);
% FI ratio distributions over time
plotCIRaw(1:chgpt, FIup_ratio_m', FIup_ratio_q(1, :)', FIup_ratio_q(2, :)', 'r');
hold on;
plotCIRaw(chgpt+1:nday, FIdown_ratio_m', FIdown_ratio_q(1, :)', FIdown_ratio_q(2, :)', 'b');
plot(chgpt*ones(1, 2), [0 0.45], 'k--');
hold off; box off; grid off;
ylabel('$\frac{\mathbf{I}(R|C_1^t)}{\mathbf{I}(R|I_1^t)}$', 'FontSize', fnt);
ylim([0 0.45]);

subplot(2, 2, [3 4]);
hold on;
for i = 1:Mcol
    stairs(1:nday, Icomb(i, :), 'Color', cmap(i, :), 'LineWidth', 1);
end
stairs(1:nday, Iday, 'k', 'LineWidth', 2);
plot(chgpt*ones(1, 2), [0 3000], 'k--');
hold off; grid off; box off; xlim([tday(1) tday(end)])
xlabel('$t$ (days)', 'FontSize', fnt);
ylabel('$C_t$', 'FontSize', fnt);

axes('Position',[0.65 0.25 0.2 0.1]);
plot(tday, Rtrue, 'Color', 'r', 'LineWidth', 2);
hold on; plot(tday, ones(size(tday)), 'k--', 'LineWidth', 2)
grid off; box off; hold off;
xlim([tday(1) tday(end)]);
xlabel('$t$', 'FontSize', fnt);
title('R numbers $R_t$', 'FontSize', fnt);




