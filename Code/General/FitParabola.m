function [p1,p2,e1,e2,mu,alpha,beta,mse1,mse2,ind_s,pts] = FitParabola(pts,shuffle)

if shuffle
    for i = 1:length(shuffle)
    pts(:,shuffle(i)) = pts(randperm(size(pts,1), size(pts,1)),shuffle(i));
    end
    % pts(:,2) = pts(randperm(size(pts,1), size(pts,1)),2);
    % pts(:,3) = pts(randperm(size(pts,1), size(pts,1)),3);
end
mu = mean(pts,1)';  
[~,~,V] = svd(pts - mu','econ');
e1    = V(:,1);    % first plane axis
e2    = V(:,2);    % second plane axis
% [coef, score] = pca((pts));
% e1 = coef(:,1);
% e2 = coef(:,2);
% (V(:,3) is the normal to that plane)
% 1) center
Xc = pts - mu';      % N×3

% 2) project onto the plane basis
alpha =  Xc * e1;   % N×1 coordinate along e1
beta  =  Xc * e2;   % N×1 coordinate along e2

data2D = [alpha, beta];                % from your PCA projection step
[~, dist] = knnsearch(pts, pts, 'K', 6);
avgDist = mean(dist(:,2:end),2);
weights = avgDist;                     % N×1 weight per point
w_sqrt = sqrt(weights); 
ind_s = 1:length(weights);
% [alpha, ind_s] = sort(alpha);
% beta = beta(ind_s);
% w_sqrt  = diff(alpha);
% w_sqrt (end+1) = mean(w_sqrt);

% 3) rebuild in 3D but only in the e1–e2 directions
pts_flat = mu' ...
         + alpha*e1' ...     % moves you along e1
         + beta *e2';        % plus moves along e2


% p1 = polyfit(alpha, beta, 1);   % [A, B, C]
a1 = fit(alpha,beta,"poly1",'Weights', w_sqrt);
p1 = coeffvalues(a1);




y1 =  p1(1).*alpha + p1(2);

mse1 = rmse(beta, y1);


% p2 = polyfit(alpha, beta, 2);   % [A, B, C]
a2 = fit(alpha,beta,"poly2",'Weights', w_sqrt);
p2 = coeffvalues(a2);

y2 = p2(1).*alpha.^2 + p2(2).*alpha+ p2(3);

mse2 = rmse(beta, y2);

end