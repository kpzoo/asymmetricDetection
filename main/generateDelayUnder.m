% Generate delayed and underreported epidemic curves
function [Icomb, Lcomb, ttrunc] = generateDelayUnder(totcase, nday, Iday, M, Pomega)

% Assumptions and notes
% - applies a stochastic delay and a binomial sampling to true incidence
% - also computes total infectiousness of these noisy cases

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

%% Generate M under-reported and delayed versions of Iday

% Mean of delay distribution and NB params
mtau = 10.8; r = 10; p = mtau/(r + mtau);
% Mean of sampling distribution and beta params
rho = 0.38; b = 20; fr = rho/(1 - rho); a = fr*b; 

% Under-reporting distribution
xrep = 0:0.01:1; yrep = betapdf(xrep, a, b);
% Draws from delay distribution
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

% Apply under-reporting to the already delayed curves
Icomb = zeros(size(Itrunc)); Lcomb = Icomb;
for i = 1:M
    for j = 1:nday
        % Main downsampling 
        Icomb(i, j) = binornd(Itrunc(i, j), betarnd(a, b));
        % Compute noisy total infectiousness
        Lcomb(i, j) = sum(Icomb(i, j-1:-1:1).*Pomega(1:j-1));
    end
end

% Find truncation points to protect against initial 0s
ttrunc = zeros(1, M);
for i = 1:M
    ttrunc(i) = find(Lcomb(i, :) > 0, 1, 'first');
end