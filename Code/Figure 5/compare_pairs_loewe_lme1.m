function results = compare_pairs_loewe_lme1(results, varargin)
%COMPARE_PAIRS_LOEWE_LME Pairwise comparisons of Loewe (logCI) between pairs.
%
% Model options (recommended defaults shown):
%   - Use trials directly (no pre-averaging), OR
%   - Use day means but weight by number of trials per (Pair,Date)
%
% LME (effects coding so intercept is grand mean):
%   logCI ~ 1 + Pair + (1|Date)
%
% Adds to results:
%   results.loewe_emm            (EMMs per Pair with SE on log scale)
%   results.loewe_pairwise_lme   (all pairwise contrasts with p & q values)
%
% INPUT
%   results.combo_detail must contain columns:
%     - Pair (string/categorical)
%     - Date (string/categorical)
%     - logCI_Loewe (numeric; per-trial log CI at p_mix)
%     - CI_mix (numeric; used for effect-window filtering)
%
% OPTIONS
%   'MinDays'     (default 2)    : require >= MinDays day means per pair when UseDayMeans=true
%   'DFMethod'    (default 'satterthwaite') : df method for tests
%   'CIBound'     (default 1)    : scalar upper bound, or [low high] effect window on CI_mix
%   'UseDayMeans' (default false): if true, collapse to day means with weights
%
% EXAMPLE
%   results = compare_pairs_loewe_lme(results, 'CIBound',[0.6 1.2], 'UseDayMeans',false);

p = inputParser;
p.addParameter('MinDays', 2, @(x) isscalar(x) && x>=1);
p.addParameter('DFMethod', 'satterthwaite', @(s) ischar(s) || isstring(s));
p.addParameter('CIBound', 1, @(v) isnumeric(v) && (isscalar(v) || numel(v)==2));
p.addParameter('UseDayMeans', false, @(b) islogical(b) || ismember(b,[0 1]));
p.parse(varargin{:});
minDays   = p.Results.MinDays;
dfmethod  = char(p.Results.DFMethod);
CIBound   = p.Results.CIBound;
useMeans  = logical(p.Results.UseDayMeans);

D = results.combo_detail;

% ---------- (0) CI filter (fixed) ----------
if isscalar(CIBound)
    D = D(D.CI_mix <= CIBound, :);
else
    low  = min(CIBound);
    high = max(CIBound);
    D = D(D.CI_mix >= low & D.CI_mix <= high, :);
end

% Keep the columns we need and sanitize
needVars = {'Pair','Date','logCI_Loewe'};
if ~all(ismember(needVars, D.Properties.VariableNames))
    error('combo_detail lacks required columns: %s', strjoin(setdiff(needVars, D.Properties.VariableNames),', '));
end
D = D(isfinite(D.logCI_Loewe), :);
D.Pair = categorical(string(D.Pair));
D.Date = categorical(string(D.Date));

% ---------- (1) Optionally collapse to day means (with weights) ----------
if useMeans
    G = findgroups(string(D.Pair), string(D.AconcNum), string(D.BconcNum));
    daytbl = table();
    daytbl.Pair         = splitapply(@(x) x(1), D.Pair, G);
    daytbl.Date         = splitapply(@(x) x(1), D.Date, G);
    daytbl.logCI_daymean= splitapply(@mean, D.logCI_Loewe, G);
    daytbl.n            = splitapply(@numel, D.logCI_Loewe, G);
    daytbl = daytbl(isfinite(daytbl.logCI_daymean), :);

    if minDays > 1
        [Gc, pairs_c] = findgroups(daytbl.Pair);
        ndays_per_pair = splitapply(@numel, daytbl.logCI_daymean, Gc);
        keep_pairs = pairs_c(ndays_per_pair >= minDays);
        daytbl = daytbl(ismember(daytbl.Pair, keep_pairs), :);
    end

    if height(daytbl) < 2
        error('Not enough day means to fit the model.');
    end

    tbl = daytbl;
    respVar = 'logCI_daymean';
    weights = daytbl.n;
else
    if height(D) < 2
        error('Not enough trials to fit the model.');
    end
    tbl = D;
    tbl.n = ones(height(tbl),1); % for a consistent 'Weights' arg
    respVar = 'logCI_Loewe';
    weights = tbl.n;
end

% ---------- (2) Fit LME with effects coding & consistent DF ----------
lme = fitlme(tbl, sprintf('%s ~ 1 + Pair + (1|Date)', respVar), ...
             'DummyVarCoding','effects', 'FitMethod','REML', 'Weights',weights);

[beta, ~, feStats] = fixedEffects(lme, 'DFMethod', dfmethod); %#ok<ASGLU>
CovB = lme.CoefficientCovariance;
coefNames = string(lme.CoefficientNames);

% Effects coding: for K levels of Pair, there are K-1 coefficients named "Pair_<level>"


idxIntercept = find(coefNames == "(Intercept)");
if isempty(idxIntercept), error('Intercept not found in fixed effects.'); end

% --- after you've built lme, beta, CovB, coefNames, idxIntercept ---

% Use string array for robust equality
pairLevels = string(categories(tbl.Pair));
K = numel(pairLevels);


% Get coefficient index (effects coding) for a given level name.
% Under effects coding with K levels, there are K-1 columns:
%   "Pair_<level1>", ..., "Pair_<level{K-1}>"; the Kth level is the omitted baseline.
idx_for_level = @(lev) ( ...
    find(coefNames == "Pair_" + lev, 1) ...
);
% Fallback relaxed matcher (older MATLAB sometimes formats names differently)
relaxed_idx_for_level = @(lev) ( ...
    find(contains(coefNames, lev) & contains(coefNames, "Pair"), 1) ...
);

% Precompute indices for the K-1 explicit columns (omit the last level)
idxAll = nan(K-1,1);
for k = 1:K-1
    idxk = idx_for_level(pairLevels(k));
    if isempty(idxk)
        idxk = relaxed_idx_for_level(pairLevels(k));
    end
    if isempty(idxk)
        error('Could not find fixed-effect coefficient for Pair level: %s', pairLevels(k));
    end
    idxAll(k) = idxk;
end

% Helper: EMM contrast vector 'a' for a given level (effects coding)
function a = a_vec_for_level(levelName)
    % Ensure string type
    levelName = string(levelName);

    a = zeros(numel(beta),1);
    a(idxIntercept) = 1;

    if levelName == pairLevels(K)  % omitted baseline under effects coding
        % baseline mean = intercept - sum(other level effects)
        a(idxAll) = -1;
    else
        % add the coefficient of this level
        idxk = idx_for_level(levelName);
        if isempty(idxk)
            idxk = relaxed_idx_for_level(levelName);
        end
        if isempty(idxk)
            error('Could not find fixed-effect coefficient for Pair level: %s', levelName);
        end
        a(idxk) = 1;
    end
end

% ---------- (3) EMMs per Pair (log scale) ----------
EMM = nan(K,1);
SE  = nan(K,1);
for k = 1:K
    a = a_vec_for_level(pairLevels(k));
    EMM(k) = a.' * beta;
    SE(k)  = sqrt(max(a.' * CovB * a, 0));
end
emmTbl = table(string(pairLevels), EMM, SE, ...
    'VariableNames', {'Pair','EMM_logCI','SE'});

% ---------- (4) Pairwise contrasts (EMM_i - EMM_j) with consistent DF ----------
rows = {};
for i = 1:K-1
    for j = i+1:K
        Pi = pairLevels(i);
        Pj = pairLevels(j);

        ai = a_vec_for_level(Pi);
        aj = a_vec_for_level(Pj);
        L  = (ai - aj).';                  % row vector over fixed effects

        diff_ij = L * beta;                % on log scale
        se_ij   = sqrt(max(L * CovB * L.', 0));
        tstat   = diff_ij / max(se_ij, eps);

        [p_two, F, df1, df2] = coefTest(lme, L, 0, 'DFMethod', dfmethod); %#ok<ASGLU>

        if df1 == 1 && isfinite(df2)
            p_left  = tcdf(tstat, df2);      % H1: diff < 0  (Pi < Pj)
            p_right = 1 - p_left;            % H1: diff > 0
            df_use  = df2;
        else
            p_left  = NaN; p_right = NaN; df_use = df2;
        end

        rows(end+1,:) = {char(Pi), char(Pj), diff_ij, se_ij, tstat, df_use, p_two, p_left, p_right}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'Pair1','Pair2','Diff_mean_logCI','SE','t','df','p_two_sided','p_left','p_right'});

% ---------- (5) FDR (BH), skipping NaNs safely ----------
T.q_two_sided = nan(height(T),1);
T.q_left      = nan(height(T),1);
T.q_right     = nan(height(T),1);

mask = isfinite(T.p_two_sided);
if any(mask), T.q_two_sided(mask) = mafdr(T.p_two_sided(mask), 'BHFDR', true); end
mask = isfinite(T.p_left);
if any(mask), T.q_left(mask)      = mafdr(T.p_left(mask),      'BHFDR', true); end
mask = isfinite(T.p_right);
if any(mask), T.q_right(mask)     = mafdr(T.p_right(mask),     'BHFDR', true); end

results.loewe_pairwise_lme = T;
results.loewe_emm          = emmTbl;
results.loewe_lme_model    = lme;

results.combo_detail       = D;


end
