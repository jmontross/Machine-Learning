
%% Distance of a point (x) to a set (T), defined as the minimal distance from it to T. (Meaning, to the point 
% in T that has the shortest distance to x.
% Distance can be L1 norm, L2 norm, or Linf norm

% x is a d-dimensional column vector, representing a point in the space
% T is a d-dimensional matrix, (# of rows = d, number of columns = number
% of points in the T-space), where each column represents a point in the
% space.
% L: Norm being used. '1' for L1 norm, '2' for L2 norm, 'inf' for Linf norm

function [D] = rhoDistance(L, x, T)    
    
    switch L
        case 1
            D = min((sum( abs(repmat(x, 1, size(T,2)) - T)) )); %L1
        case 2
            D = min(sqrt(sum( (repmat(x, 1, size(T,2)) - T).^2) )); %L2        
        case inf
            D = min(max( abs(repmat(x, 1, size(T,2)) - T)) ); %L2                         
    end
end