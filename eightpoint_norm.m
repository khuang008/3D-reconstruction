function F = eightpoint_norm(pts1, pts2, normalization_constant)
%%estimate the fundamental matrix by using eight-point algorithm
% input:  
%   pts1: 2 by N matrice corresponding to the (x, y) coordinates (one point
%         per column) of the points in the first image
%   pts2: 2 by N matrice corresponding to the (x, y) coordinates (one point
%         per column) of the points in the second image
%   normalization constant: to avoid numerical issues
% output: 
%   F: the estimation of fundamental matrix 

    N = size(pts1,2);% Number of point correspondence
    
    %Compute  normalize transformation
    T=[1./normalization_constant  0         0;
        0       1./normalization_constant   0;
        0       0                           1];
    %normalize to avoid numerical issues
    pts1_norm=[pts1./normalization_constant;ones(1,N)];
    pts2_norm=[pts2./normalization_constant;ones(1,N)];
    x1=pts1_norm(1,:)';
    y1=pts1_norm(2,:)';
    x2=pts2_norm(1,:)';
    y2=pts2_norm(2,:)';
    
    A=[x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(N,1)];
    [ ~,~,V] = svd(A);
    % get the singular vector corresponding to smallest singular value to 
    % form F
    F=vec2mat(V(:,9),3,3);
    %Compute the singular value decomposition Udiag(r, s, t)V^T of F, and
    %set F= Udiag(r, s, 0)V^T
    [U,S,V]=svd(F);
    S(3,3)=0;
    F=U*S*V';
    % get the final estimation of the fundamental matrix
    F=T'*F*T;
end
