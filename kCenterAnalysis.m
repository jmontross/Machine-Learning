%Tarin Ziyaee
%Analysis of Farthest-First Traversal Implementation of K-center Clustering Algorithm

clear all;
color = 'y';

N.ptsInS = 100;
N.k = 10;
N.dim = 2;
N.subPtsInS = round(N.ptsInS/3);
N.monteCarlo = 1000;

for kk = 1:N.monteCarlo
    clear S
   S = [randn(N.dim, N.subPtsInS) (4 + .5.*randn(N.dim, N.subPtsInS)) ([-4;6]*ones(1, N.subPtsInS) +  .2.*randn(N.dim, N.subPtsInS))];
%     S = [randn(N.dim, N.subPtsInS) (4 + .5.*randn(N.dim, N.subPtsInS)) ([-4;6;-1]*ones(1, N.subPtsInS) +  .2.*randn(N.dim, N.subPtsInS))];
   
N.ptsInS = size(S,2);
    %S = randn(N.dim, N.ptsInS);

    clear T
    T = S(:,randi(N.k));
    tCounter = 1;   

    while ( size(T,2) < N.k)       
        for pp = 1:N.ptsInS       
            minDistances(pp) = rhoDistance(2,S(:,pp), T);
        end
        tCounter = tCounter + 1;
        [theMax maxInd] = max(minDistances);
        T(:,tCounter) = S(:,maxInd);       
    end
    eVals(kk) = max(minDistances);
    hold off;
end

%Last instantiation of S and T
figure(1);
plot(S(1,:),S(2,:),'*k', 'linewidth', 3); grid on; hold on;
plot(T(1,:), T(2,:),'*', 'color', color,'linewidth', 4);
set(gca,'color', [0 .5 .9]);
set(gcf,'color', [1 1 1]);
axis([-10 10 -10 10]);

figure(2);
subplot(2,1,1);
plot(eVals,'.'); grid on;
set(gcf,'color', [1 1 1]);
set(gca,'color', [1 1 1]);
axis([0 N.monteCarlo 0 max(eVals) + 1]);
subplot(2,1,2);
[aa bb] = hist(eVals, 10);
hist(eVals, 10); grid on;
set(gcf,'color', [1 1 1]);
set(gca,'color', [1 1 1]);
axis([0 2*max(bb) 0 max(aa) + 10]);



