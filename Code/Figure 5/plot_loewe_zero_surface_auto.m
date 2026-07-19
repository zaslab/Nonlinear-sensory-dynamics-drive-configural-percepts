
function surf_dist = plot_loewe_zero_surface_auto(results, overlayTbl, varargin)
% overlayTbl must be a subset of results.combo_detail for exactly ONE pair.
opts = inputParser;
opts.addParameter('NX',120);
opts.addParameter('NP',40);
opts.addParameter('PadP',0.02);
opts.addParameter('UseDayMeans', true);
opts.addParameter('PlotScatter', true);
opts.parse(varargin{:});
NX=opts.Results.NX; NP=opts.Results.NP; pad=opts.Results.PadP;
PlotScatter = opts.Results.PlotScatter;

% Determine A,B orientation from data (so axes match points)
Adata = string(overlayTbl.A(1));
Bdata = string(overlayTbl.B(1));
calib  = results.calib_models;

% Day means if requested
TblUse = overlayTbl;
if opts.Results.UseDayMeans
    TblUse = overlay_to_daymeans(overlayTbl);
end

% Choose p-band covering what we will plot (observed & predicted)
pp = TblUse.p_mix;
if ismember('p_pred_loewe', TblUse.Properties.VariableNames)
    pp = [pp(:); TblUse.p_pred_loewe(:)];
end
pp = pp(isfinite(pp));
LA = calib.(Adata).params.L; UA = calib.(Adata).params.U;
LB = calib.(Bdata).params.L; UB = calib.(Bdata).params.U;

pmin = max([min(pp)-pad, LA, LB]) + 1e-6;
pmax = min([max(pp)+pad, UA, UB]) - 1e-6;
assert(isfinite(pmin) && isfinite(pmax) && pmin < pmax, 'Empty p-band for surface.');

% Draw surface (CI=1) with same A,B orientation as data
plot_loewe_zero_surface(calib, Adata, Bdata, ...
    'PBand',[pmin pmax], 'NX',NX, 'NP',NP, 'Overlay', []); hold on;

% Overlay day-mean points (size encodes #trials that day at that dose)
msz = 28;
if ismember('Nday', TblUse.Properties.VariableNames)
    base = min(TblUse.Nday);
    msz  = 20 + 6*(TblUse.Nday - base);
end
if opts.Results.UseDayMeans
    TblUse = groupsummary(TblUse, ["p_pred_loewe","cA_half", "cB_half"], "mean", ["p_mix"]); 
else
    TblUse = groupsummary(TblUse, ["p_pred_loewe","Date","cA_half", "cB_half"], "mean",  ["p_mix"]); 
end
TblUse.Properties.VariableNames{TblUse.Properties.VariableNames=="mean_p_mix"} = 'p_mix';
if PlotScatter
scatter3(TblUse.cA_half, TblUse.cB_half, TblUse.p_mix, 70, TblUse.p_mix, ...
         'filled','MarkerEdgeColor',[0.2 0.2 0.2],'LineWidth',0.6);

% Optional: day-mean Loewe predictions (should lie on the surface)
if ismember('p_pred_loewe', TblUse.Properties.VariableNames)
    plot3(TblUse.cA_half, TblUse.cB_half, TblUse.p_pred_loewe, ...
          'wo','MarkerSize',5,'LineWidth',1.0);

    plot3([TblUse.cA_half,TblUse.cA_half]', ...
          [TblUse.cB_half,TblUse.cB_half]', ...
          [TblUse.p_pred_loewe,TblUse.p_mix]', ...
          'r','LineWidth',1.5);

end
end
p_mix = TblUse.p_mix;
p_pred = TblUse.p_pred_loewe;
[~, ind] = sort(p_pred);
p_mix = p_mix(ind);
p_pred = p_pred(ind);
surf_dist = (p_mix - p_pred);
surf_dist(:,2) = (p_mix - p_pred)./p_pred;
surf_dist(:,3) = p_pred;
end


% ============================================================
%  CI=1 (Loewe-additive) surface
% ============================================================
function plot_loewe_zero_surface(calib, A, B, varargin)
% Plot the Loewe-additive surface CI(p)=1 for a pair (A,B).
% Axes: x = dose of A (A/2), y = dose of B (B/2), z = effect p.
%
% Usage:
%   plot_loewe_zero_surface(results.calib_models, "DA", "IAA", ...
%       'PBand',[0.6 0.9], 'NX',120, 'NP',40, 'Overlay', overlayTbl)

opts = inputParser;
opts.addParameter('PBand', [], @(v) isempty(v) || (isnumeric(v)&&numel(v)==2));
opts.addParameter('NX', 100, @(x) isscalar(x) && x>=20);
opts.addParameter('NP', 30,  @(x) isscalar(x) && x>=10);
opts.addParameter('Overlay', [], @(t) istable(t) || isempty(t));
opts.addParameter('PlotScatter', true);
opts.parse(varargin{:});
PBand = opts.Results.PBand; NX = opts.Results.NX; NP = opts.Results.NP;
overlayTbl = opts.Results.Overlay;
PlotScatter = opts.Results.PlotScatter;

A = string(A); B = string(B);
CA = calib.(A); CB = calib.(B);

% p-range where BOTH singles invert
LA = CA.params.L; UA = CA.params.U;
LB = CB.params.L; UB = CB.params.U;
pmin = max([LA LB]) + 1e-6;
pmax = min([UA UB]) - 1e-6;
if ~isempty(PBand)
    pmin = max(pmin, PBand(1));
    pmax = min(pmax, PBand(2));
end
assert(isfinite(pmin) && isfinite(pmax) && pmin < pmax, 'Bad/empty p-range.');

% Build surface by sweeping isoboles across p
pgrid = linspace(pmin, pmax, NP);
X = nan(NX, NP); Y = nan(NX, NP); Z = nan(NX, NP);

for j = 1:NP
    p = pgrid(j);
    cAeq = CA.c_of_p(p);
    cBeq = CB.c_of_p(p);
    xLine = linspace(0, cAeq, NX);
    yLine = cBeq * (1 - xLine./cAeq);
    yLine(yLine < 0) = NaN;
    X(:,j) = xLine(:);
    Y(:,j) = yLine(:);
    Z(:,j) = p;
end

figure('Color','w'); hold on;
surf(X, Y, Z, Z, 'EdgeColor','none', 'FaceAlpha', 0.65);
colormap parula; cb = colorbar; ylabel(cb,'p (effect)');
set(gca,'XScale','log','YScale','log');
xlabel(sprintf('%s dose', A),'Interpreter','none');
ylabel(sprintf('%s dose', B),'Interpreter','none');
zlabel('Effect p'); grid on; box on; view(45,25);
title(sprintf('Loewe-additive surface (CI=1): %s — %s', A, B), 'Interpreter','none');

% Optional overlay of points
if ~isempty(overlayTbl) & PlotScatter
    need = {'cA_half','cB_half','p_mix'};
    assert(all(ismember(need, overlayTbl.Properties.VariableNames)), 'Overlay missing cA_half/cB_half/p_mix');
    sc = scatter3(overlayTbl.cA_half, overlayTbl.cB_half, overlayTbl.p_mix, ...
        70, overlayTbl.p_mix, 'filled', 'MarkerEdgeColor',[0.2 0.2 0.2], 'LineWidth',0.6);
    if ismember('p_pred_loewe', overlayTbl.Properties.VariableNames)
        plot3(overlayTbl.cA_half, overlayTbl.cB_half, overlayTbl.p_pred_loewe, ...
              'wo', 'MarkerSize',5, 'LineWidth',1.0);
        legend({'CI=1 surface','Observed mix','Loewe prediction'}, 'Location','northwest');
    else
        legend({'CI=1 surface','Observed mix'}, 'Location','northwest');
    end
end
end

% ============================================================
%  Collapse to day means per (doseA/2, doseB/2)
% ============================================================
function Tday = overlay_to_daymeans(overlayTbl)
% overlayTbl: subset of results.combo_detail for ONE pair
varsNeed = {'Date','A','B','Pair','cA_half','cB_half','p_mix'};
assert(all(ismember(varsNeed, overlayTbl.Properties.VariableNames)), 'Missing vars in overlay.');
hasPred = ismember('p_pred_loewe', overlayTbl.Properties.VariableNames);

G = findgroups(string(overlayTbl.Date), overlayTbl.cA_half, overlayTbl.cB_half);
Tday = table();
Tday.Date    = splitapply(@(x) x(1), string(overlayTbl.Date), G);
Tday.A       = splitapply(@(x) x(1), string(overlayTbl.A),    G);
Tday.B       = splitapply(@(x) x(1), string(overlayTbl.B),    G);
Tday.Pair    = splitapply(@(x) x(1), string(overlayTbl.Pair), G);
Tday.cA_half = splitapply(@(x) x(1), overlayTbl.cA_half,      G);
Tday.cB_half = splitapply(@(x) x(1), overlayTbl.cB_half,      G);
Tday.p_mix   = splitapply(@mean,      overlayTbl.p_mix,       G);
if hasPred
    Tday.p_pred_loewe = splitapply(@mean, overlayTbl.p_pred_loewe, G);
end
Tday.Nday    = splitapply(@numel,     overlayTbl.p_mix,       G);
end

