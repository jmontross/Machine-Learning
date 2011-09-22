%Tarin Ziyaee
%Farthest-First Traversal Implementation of K-center Clustering Algorithm

clear all; hold off;
color = 'y';

N.ptsInS = 100;
N.k = 7;
N.dim = 2;

%% Three easy to see gaussian clusters
N.subPtsInS = round(N.ptsInS/3);
S = [randn(N.dim, N.subPtsInS) (4 + .5.*randn(N.dim, N.subPtsInS)) ([-4;6]*ones(1, N.subPtsInS) +  .2.*randn(N.dim, N.subPtsInS)) ];
%     S = [randn(N.dim, N.subPtsInS) (4 + .5.*randn(N.dim, N.subPtsInS)) ([-4;6;-1]*ones(1, N.subPtsInS) +  .2.*randn(N.dim, N.subPtsInS))];

%% Just one simple random gaussian cluster
% S = randn(N.dim, N.ptsInS);

%% Another easy to see gaussian cluster but with it being rotated, skewed, etc.
% S = [1.*randn(1, N.ptsInS); 2.*randn(1, N.ptsInS)] + repmat([1;1],1, N.ptsInS);
% theta = randi(180)*pi/180;
% rotationMatrix = [cos(theta) -sin(theta); sin(theta) cos(theta)];
% for ii = 1:N.ptsInS
%     S(:,ii) = rotationMatrix * S(:,ii);
% end



N.ptsInS = size(S,2)
clear T
%T = S(:,randi(N.k));
T = S(:,N.ptsInS*ceil(rand(1,1)));
tCounter = 1;

axisSettings = [-10 10 -10 10];
plot([0 0], [axisSettings(1:2)], '-.k', 'linewidth', 1); hold on;
plot([axisSettings(1:2)], [0 0], '-.k', 'linewidth', 1);
plot(S(1,:),S(2,:),'*k', 'linewidth', 3); grid on; 
plot(T(1,:), T(2,:),'*', 'color', color,'linewidth', 4);
set(gcf,'color', [1 1 1]);
set(gca,'color', [0 .5 .9]);
axis(axisSettings);

while ( size(T,2) < N.k)    
    
    for pp = 1:N.ptsInS       
        minDistances(pp) = rhoDistance(2,S(:,pp), T);
    end
    tCounter = tCounter + 1;
    [theMax maxInd] = max(minDistances);
    T(:,tCounter) = S(:,maxInd);
       
    pause(1)
    plot(T(1,:), T(2,:), '*', 'color', color,'linewidth', 4);
    
end
eVal = max(minDistances) %The maximum e-cover value. 



%WEKA