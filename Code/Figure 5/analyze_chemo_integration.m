function results = analyze_chemo_integration(T, varargin)
% Analyze integration for A+B combos using single-odor calibrations.
% T is a table with columns (case-insensitive names are okay):
%   Date, Stim, ConcNum, Conc, PlateNum, CI, Component1, Component1Conc, Component2, Component2Conc
%
% Outputs:
%   results.pair_summary   : one row per pair with mean integration metrics
%   results.combo_detail   : one row per combo record with all metrics
%   results.calib_models   : struct with forward & inverse functions per Stim

p = inputParser;
p.addParameter('Eps',1e-4);               % clamp for probabilities/logits
p.addParameter('DoBootstrap',true);
p.addParameter('B',1000);                 % bootstrap draws over days
p.addParameter('Seed',1);
p.parse(varargin{:});
Eps = p.Results.Eps;

%% --- 0) normalize column names & types
% T = standardizeVarNames(T);

% Treat labels as strings
labelVars = ["Stim","Component1","Component2","Date"];
for v = labelVars
    if ismember(v,T.Properties.VariableNames)
        T.(v) = string(T.(v));
    end
end

% collapse technical repeats: same day/plate/condition -> mean CI
T = collapseTechRepeats(T);

% Map CI -> p in [0,1]
T.p = (T.CI);% + 1) / 2;
T.p = min(max(T.p, Eps), 1-Eps); % clamp

% split singles vs combos
isCombo = T.Stim == "Combo";
Ts = T(~isCombo,:);   % singles (calibration)
Tc = T(isCombo,:);    % combos

%% --- 1) Fit single-odor calibrations with a bounded 4-parameter Hill model
stims = unique(Ts.Stim);
calib = struct();

for s = stims'
    rows = Ts.Stim==s;
    cs   = Ts.Conc(rows);
    ps   = Ts.p(rows);

    ok = isfinite(cs) & isfinite(ps);
    cs = cs(ok); ps = ps(ok);
    if isempty(cs)
        warning('Stim %s: no valid calibration rows.', s); continue
    end

    % collapse repeats at same concentration (keep weights)
    [c_uniq,~,g] = unique(cs);
    p_mean = splitapply(@mean, ps, g);
    w_cnt  = splitapply(@numel, ps, g);

    % bounds & inits
    cmin = min(c_uniq);  cmax = max(c_uniq);  cg = geomean(c_uniq);
    L0 = max(0.00, min(p_mean));              % lower asymptote
    U0 = min(1.00, max(p_mean));              % upper asymptote
    p_mid = median(p_mean);
    EC0 = min(max(cg, cmin*1.2), cmax/1.2);   % EC50 near geometric mean
    n0  = 2;                                  % slope

    p0  = [L0, U0, EC0, n0];
    lb  = [0,   0.2, cmin/10,  0.1];
    ub  = [0.9, 1.0, cmax*10,  10];

    % weighted residuals in log-space of concentration
    hill = @(p,c) p(1) + (p(2)-p(1))./(1 + (p(3)./c).^max(p(4),0.1));
    resid = @(p) (sqrt(w_cnt(:)) .* (hill(p, c_uniq(:)) - p_mean(:)));

    % fit with lsqnonlin (Optimization Toolbox) or fminsearch fallback
    haveOPT = exist('lsqnonlin','file')==2;
    if haveOPT
        opts = optimoptions('lsqnonlin','Display','off');
        try
            p_hat = lsqnonlin(@(p) resid(p), p0, lb, ub, opts);
        catch
            % fallback to unconstrained, then clip to bounds
            p_hat = fminsearch(@(p) sum(resid(p).^2), p0, optimset('Display','off'));
            p_hat = max(lb, min(ub, p_hat));
        end
    else
        % simple fallback
        p_hat = fminsearch(@(p) sum(resid(p).^2), p0, optimset('Display','off'));
        p_hat = max(lb, min(ub, p_hat));
    end

    L=p_hat(1); U=p_hat(2); EC=p_hat(3); n=p_hat(4);
    
    % forward & inverse (analytical inverse inside [L,U])
    p_of_c = @(c) L + (U-L)./(1 + (EC./max(c,eps)).^max(n,0.1));
    tol  = 1e-9;                  % small guard to avoid divide-by-zero
    nEff = max(n, 0.1);
    
    c_of_p = @(p) EC .* ( ((max(min(p,U-tol),L+tol) - L) ./ ...
                        max(U - max(min(p,U-tol),L+tol), tol)) ) .^ (1./nEff);


    % store
    calib.(s).model = 'hill4';
    calib.(s).params = struct('L',L,'U',U,'EC50',EC,'n',n);
    calib.(s).p_of_c = @(c) min(max(p_of_c(c),0),1);
    calib.(s).c_of_p = @(p) c_of_p(min(max(p, L+1e-6), U-1e-6)); % clamp to invertible range
    calib.(s).x_min = log10(cmin);  calib.(s).x_max = log10(cmax);
    calib.(s).y_min = logit(max(L,1e-6),Eps);  calib.(s).y_max = logit(min(U,1-1e-6),Eps);
    
end
plot_calibration_grid(Ts, calib, unique(Ts.Stim)');
%% --- 2) Build lookup from (Stim,ConcNum) -> Conc used in singles
keys  = strcat(Ts.Stim,"|",string(Ts.ConcNum));
valsC = Ts.Conc;
concLUT = containers.Map(keys, valsC);

% helper to fetch half-concentration used in combo
getHalfConc = @(stim,concnum) 0.5 * fetchConc(concLUT, stim, concnum);

%% --- 3) Compute predictions & integration metrics per combo row
n = height(Tc);
detail = Tc(:,{'Date','PlateNum','Component1','Component1Conc','Component2','Component2Conc','CI','p'});
detail.Properties.VariableNames = {'Date','Plate','A','AconcNum','B','BconcNum','CI_mix','p_mix'};
detail.A = string(detail.A); detail.B = string(detail.B);

% canonical (unordered) pair label
detail.Pair = arrayfun(@(a,b) sprintf('%s — %s', sort_pair(a,b)), detail.A, detail.B, 'UniformOutput', false);

% numeric concentrations for half-doses
detail.cA_half = arrayfun(@(a,an) getHalfConc(a,an), detail.A, detail.AconcNum);
detail.cB_half = arrayfun(@(b,bn) getHalfConc(b,bn), detail.B, detail.BconcNum);

% predicted p’s from single-odor models
detail.pA_half = arrayfun(@(a,c) forward_p(calib,a,c,Eps), detail.A, detail.cA_half);
detail.pB_half = arrayfun(@(b,c) forward_p(calib,b,c,Eps), detail.B, detail.cB_half);

% Null predictions
detail.p_pred_logit = ilogit( logit(detail.pA_half,Eps) + logit(detail.pB_half,Eps), Eps );                % product rule
detail.p_pred_bliss = detail.pA_half + detail.pB_half - detail.pA_half.*detail.pB_half;                    % Bliss
detail.p_pred_hsa   = max(detail.pA_half, detail.pB_half);                                                 % HSA
detail.p_pred_loewe = arrayfun(@(a,b,x,y) ...
    loewe_effect_at_dose(calib, a, b, x, y), ...
    detail.A, detail.B, detail.cA_half, detail.cB_half);
% Loewe combination index at observed effect p_mix
detail.CI_Loewe = arrayfun(@(a,b,p,cAh,cBh) loeweCI(calib,a,b,p,cAh,cBh,Eps), ...
                           detail.A, detail.B, detail.p_mix, detail.cA_half, detail.cB_half);

% Integration scores
detail.S_logit  = logit(detail.p_mix,Eps) - (logit(detail.pA_half,Eps)+logit(detail.pB_half,Eps));
detail.S_Bliss  = detail.p_mix - detail.p_pred_bliss;
detail.S_HSA    = detail.p_mix - detail.p_pred_hsa;
detail.logCI_Loewe = log(detail.CI_Loewe);

%% --- 4) Aggregate per pair (mean ± bootstrap over days)
pairGroups = findgroups(detail.Pair);
pairs = splitapply(@(x) x(1), detail.Pair, pairGroups);

agg = table(pairs, ...
    splitapply(@mean, detail.S_logit,  pairGroups), ...
    splitapply(@mean, detail.S_Bliss,  pairGroups), ...
    splitapply(@nanmean, detail.logCI_Loewe, pairGroups), ...
    splitapply(@mean, detail.S_HSA,    pairGroups), ...
    splitapply(@numel, detail.S_logit, pairGroups), ...
    'VariableNames',{'Pair','Slogit_mean','Sbliss_mean','logCI_Loewe_mean','SHSA_mean','N'});
detail.Loewe_defined = isfinite(detail.CI_Loewe);

% Pair-level aggregation 
pairGroups = findgroups(detail.Pair);

logCI_vals = detail.logCI_Loewe;
is_def     = detail.Loewe_defined;

logCI_mean_nan = splitapply(@(x) nanmean(x), logCI_vals, pairGroups); % your current
n_defined      = splitapply(@(x) sum(x),       is_def,     pairGroups);
n_total        = splitapply(@numel,            is_def,     pairGroups);
coverage       = n_defined ./ n_total;

agg.logCI_Loewe_mean = logCI_mean_nan;
agg.Loewe_n_defined  = n_defined;
agg.Loewe_n_total    = n_total;
agg.Loewe_coverage   = coverage;

% Optional: mask low coverage so you don’t over-interpret
min_cov = 0.5;  % choose threshold
mask_low = coverage < min_cov;
agg.logCI_Loewe_mean(mask_low) = NaN;   % or keep but mark
% optional bootstrap: resample Combo rows by Date (keeps plate-level structure within a day)
if p.Results.DoBootstrap
    rng(p.Results.Seed);
    days = unique(detail.Date);
    Gd  = findgroups(string(detail.Pair), string(detail.AconcNum), string(detail.BconcNum));

    B     = p.Results.B;
    q025  = @(x) quantile(x,0.025);
    q975  = @(x) quantile(x,0.975);

    [Slog_lo,Slog_hi,Sbl_lo,Sbl_hi,LCI_lo,LCI_hi] = deal(nan(height(agg),1));

    for i=1:height(agg)
        pairMask = strcmp(detail.Pair,agg.Pair(i));
        % which days have this pair?
        days_i = unique(detail.Date(pairMask));
        if numel(days_i)<2
            continue; % not enough days for bootstrap CI
        end
        % indices of rows belonging to those days & this pair
        idx = find(pairMask);

        % bootstrap by days: resample days_i with replacement
        slog = nan(B,1); sbl = nan(B,1); lci = nan(B,1);
        for b=1:B
            boot_rows = [];
            boot_days = days_i(randi(numel(days_i), numel(days_i), 1));
            for d = boot_days'
                boot_rows = [boot_rows; idx(detail.Date(idx)==d)];
            end
            slog(b) = mean(detail.S_logit(boot_rows));
            sbl(b)  = mean(detail.S_Bliss(boot_rows));
            lci(b)  = nanmean(detail.logCI_Loewe(boot_rows));
        end
        Slog_lo(i)=q025(slog);  Slog_hi(i)=q975(slog);
        Sbl_lo(i) =q025(sbl);   Sbl_hi(i) =q975(sbl);
        LCI_lo(i) =q025(lci);   LCI_hi(i) =q975(lci);
    end

    agg.Slogit_CI = [Slog_lo Slog_hi];
    agg.Sbliss_CI = [Sbl_lo  Sbl_hi];
    agg.logCI_Loewe_CI = [LCI_lo LCI_hi];
end

% sort: most synergistic/configural first (positive S_logit / S_Bliss; negative logCI_Loewe)
agg = sortrows(agg, {'Slogit_mean','Sbliss_mean','logCI_Loewe_mean'}, {'descend','descend','ascend'});

%% --- 5) pack outputs
results.pair_summary = agg;
results.combo_detail = detail;
results.calib_models = calib;
end


function T2 = collapseTechRepeats(T)
% average CI within identical condition on same Date × PlateNum
keys = ["Date","PlateNum","Stim","ConcNum","Conc","Component1","Component1Conc","Component2","Component2Conc"];
have = keys(ismember(keys, T.Properties.VariableNames));
[G,~,subs] = unique(T(:,have),'rows');
CImean = splitapply(@mean, T.CI, subs);
T2 = G; T2.CI = CImean;
end

function y = logit(p, eps)
p = min(max(p,eps),1-eps);
y = log(p./(1-p));
end

function p = ilogit(y, eps)
p = 1./(1+exp(-y));
p = min(max(p,eps),1-eps);
end

function c = fetchConc(LUT, stim, concnum)
key = strcat(string(stim),"|",string(concnum));
if isKey(LUT, key)
    c = LUT(key);
else
    error('No single-odor row for key %s', key);
end
end

function p = forward_p(calib, stim, c, eps)
stim = string(stim);
if ~isfield(calib, stim)
    error('No calibration model for stim %s', stim);
end
p = calib.(stim).p_of_c(c);
p = min(max(p,eps),1-eps);
end

function CI = loeweCI(calib, A, B, p_mix, cA_half, cB_half, eps)
% Compute Loewe Combination Index at observed effect p_mix.
A = string(A); B = string(B);
% dose equivalents to reach p_mix alone (invert the single curves)
try
    cA_eq = calib.(A).c_of_p(p_mix);
    cB_eq = calib.(B).c_of_p(p_mix);
catch
    CI = NaN; return
end
CI = (cA_half/max(cA_eq,eps)) + (cB_half/max(cB_eq,eps));
end

function s = sort_pair(a,b)
ab = sort([string(a) string(b)]);
s = strjoin(ab,' — ');
end
function p_pred = loewe_effect_at_dose(calib, A, B, x, y)
    % returns p such that CI_Loewe(p) = 1 for doses x,y (A/2,B/2)
    A=string(A); B=string(B);
    L = max(calib.(A).params.L, calib.(B).params.L) + 1e-9;
    U = min(calib.(A).params.U, calib.(B).params.U) - 1e-9;
    f = @(p) x./calib.(A).c_of_p(p) + y./calib.(B).c_of_p(p) - 1;
    % bracket in [L,U]; fall back to fzero with midpoint if needed
    p_lo=L; p_hi=U;
    if f(p_lo)>0 && f(p_hi)>0, p_pred = NaN; return; end  % outside range
    if f(p_lo)<0 && f(p_hi)<0, p_pred = NaN; return; end
    p_pred = fzero(f, [p_lo p_hi]);  % requires Continuous single-odor fits
end

function plot_calibration_grid(Ts, calib, whichStims)
% whichStims: string array (e.g., ["DA","IAA","NaCl"])

n = numel(whichStims);
ncol = min(3,n); nrow = ceil(n/ncol);
figure('Color','w');

for k = 1:n
    stim = whichStims(k);
    rows = string(Ts.Stim)==stim;
    c = Ts.Conc(rows); p = Ts.p(rows);

    subplot(nrow,ncol,k); hold on;

    scatter(c, p, 14, [0.6 0.6 0.6], 'filled', 'MarkerFaceAlpha', 0.7);

    [c_uniq,~,g] = unique(c);
    p_mean = splitapply(@mean, p, g);
    p_sem  = splitapply(@(x) std(x)/sqrt(numel(x)), p, g);
    errorbar(c_uniq, p_mean, p_sem, 'k.', 'LineWidth', 1.0, 'MarkerSize', 12);

    cmin = min(c_uniq); cmax = max(c_uniq);
    cgrid = logspace(log10(cmin), log10(cmax), 400);
    yfit  = calib.(stim).p_of_c(cgrid);
    plot(cgrid, yfit, 'r-', 'LineWidth', 1.8);

    set(gca,'XScale','log'); xlim([cmin*0.9, cmax*1.1]); ylim([0 1]);
    title(stim, 'Interpreter','none');
    if k> (nrow-1)*ncol, xlabel('conc'); end
    if mod(k-1,ncol)==0, ylabel('p'); end
    grid on; box on;
    h = gca;
    h.XTick = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2];
end
end
