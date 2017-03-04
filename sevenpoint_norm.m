function F = sevenpoint_norm(pts1, pts2, normalization_constant)
%%estimate the fundamental matrix by using seven-point algorithm
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
T=[1./ normalization_constant 0 0;
        0  1./normalization_constant 0;
        0 0 1];
    
%normalize to avoid numerical issues
pts1=[pts1./normalization_constant;ones(1,N)];
pts2=[pts2./normalization_constant;ones(1,N)];
x1=pts1(1,:)';
y1=pts1(2,:)';
x2=pts2(1,:)';
y2=pts2(2,:)';

U=[x2.*x1 x2.*y1 x2 y2.*x1 y2.*y1 y2 x1 y1 ones(N,1)];
[~, ~, V] = svd(U);

% get the singular vectors corresponding to smallest and second smallest singular value 
% and compute the F1 and F2
F1 = vec2mat( V(:,9),3,3);
F2 = vec2mat(V(:,8),3,3);

%cumpute the coefficient of det|x*F1+(1-x)*F2|

a=F1(1,1);b=F1(1,2);c=F1(1,3);
d=F1(2,1);e=F1(2,2);f=F1(2,3);
g=F1(3,1);h=F1(3,2);i=F1(3,3);

j=F2(1,1);k=F2(1,2);l=F2(1,3);
m=F2(2,1);n=F2(2,2);o=F2(2,3);
p=F2(3,1);q=F2(3,2);r=F2(3,3);

% I use this kind of 'hard code' way to get the compute faster so I can
% use more iterations.
f_coeffs =[ a*e*i - a*f*h - b*d*i + b*f*g + c*d*h - c*e*g - a*e*r + a*f*q + a*h*o - a*i*n + b*d*r - b*f*p - b*g*o + b*i*m - c*d*q + c*e*p + c*g*n - c*h*m - d*h*l + d*i*k + e*g*l - e*i*j - f*g*k + f*h*j + a*n*r - a*o*q - b*m*r + b*o*p + c*m*q - c*n*p - d*k*r + d*l*q + e*j*r - e*l*p - f*j*q + f*k*p + g*k*o - g*l*n - h*j*o + h*l*m + i*j*n - i*k*m - j*n*r + j*o*q + k*m*r - k*o*p - l*m*q + l*n*p, 
  a*e*r - a*f*q - a*h*o + a*i*n - b*d*r + b*f*p + b*g*o - b*i*m + c*d*q - c*e*p - c*g*n + c*h*m + d*h*l - d*i*k - e*g*l + e*i*j + f*g*k - f*h*j - 2*a*n*r + 2*a*o*q + 2*b*m*r - 2*b*o*p - 2*c*m*q + 2*c*n*p + 2*d*k*r - 2*d*l*q - 2*e*j*r + 2*e*l*p + 2*f*j*q - 2*f*k*p - 2*g*k*o + 2*g*l*n + 2*h*j*o - 2*h*l*m - 2*i*j*n + 2*i*k*m + 3*j*n*r - 3*j*o*q - 3*k*m*r + 3*k*o*p + 3*l*m*q - 3*l*n*p, 
  a*n*r - a*o*q - b*m*r + b*o*p + c*m*q - c*n*p - d*k*r + d*l*q + e*j*r - e*l*p - f*j*q + f*k*p + g*k*o - g*l*n - h*j*o + h*l*m + i*j*n - i*k*m - 3*j*n*r + 3*j*o*q + 3*k*m*r - 3*k*o*p - 3*l*m*q + 3*l*n*p, 
  j*n*r - j*o*q - k*m*r + k*o*p + l*m*q - l*n*p];

% get the real solutions
R=roots(f_coeffs);
real_roots=[];
for i = 1:length(R)
    if isreal(R(i))
       real_roots=[real_roots R(i)]; 
    end
end
% compute fundamental matrix estimations
F = cell(1,length(real_roots));
    for i = 1: length(F)
        F{i} = real_roots(i).*F1 + (1-real_roots(i)).*F2;
        F{i} = T'*F{i}*T;
    end
    
end
