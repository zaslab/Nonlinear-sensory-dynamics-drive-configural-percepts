% ============================================================
%  OPTIONAL: Loewe-null predicted effect at (x,y) doses
% ============================================================
function p_pred = loewe_effect_at_dose(calib, A, B, x, y)
% Solve for p such that x/cA*(p) + y/cB*(p) = 1 (the Loewe-null effect).
A=string(A); B=string(B);
LA = calib.(A).params.L; UA = calib.(A).params.U;
LB = calib.(B).params.L; UB = calib.(B).params.U;
p_lo = max([LA LB]) + 1e-9;
p_hi = min([UA UB]) - 1e-9;
if ~(isfinite(p_lo) && isfinite(p_hi) && p_lo<p_hi)
    p_pred = NaN; return
end
f = @(p) x./calib.(A).c_of_p(p) + y./calib.(B).c_of_p(p) - 1;
% If not bracketing due to numerical quirks, try mid-point seed
if (f(p_lo)>0 && f(p_hi)>0) || (f(p_lo)<0 && f(p_hi)<0)
    pm = 0.5*(p_lo+p_hi);
    try
        p_pred = fzero(f, pm);
    catch
        p_pred = NaN;
    end
else
    p_pred = fzero(f, [p_lo p_hi]);
end
end
