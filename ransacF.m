function [F, inliers] = ransacF(pts1, pts2, normalization_constant)
%%estimate the fundamental matrix by using  RANSAC
% input:  
%   pts1: 2 by N matrice corresponding to the (x, y) coordinates (one point
%         per column) of the points in the first image
%   pts2: 2  by N matrice corresponding to the (x, y) coordinates (one point
%         per column) of the points in the second image
%   normalization constant: to avoid numerical issues
% output: 
%   F: the estimation of fundamental matrix 
%   inliers: indices of the inliers of F .


% threshold for inliers 
threshold=0.001;
% number of iterations
num_of_iter=10000;

inliers_temp=[];
inliers=[];
most_inliers_count=0;
N=size(pts1,2);
x1=[pts1; ones(1,N)];
x2=[pts2; ones(1,N)];
   for i = 1:num_of_iter 
       r=randperm(N,7);
       F_temp=sevenpoint_norm(pts1(:,r), pts2(:,r), normalization_constant);
        
        for j=1:length(F_temp) 
            inliers_count=0;
            inliers_temp=[];
            for k = 1:N
                dist=abs(x2(:,k)'* F_temp{j}*x1(:,k));               
                if(dist<threshold)
                    inliers_temp=[inliers_temp k];
                    inliers_count=inliers_count+1;
                end
            end
 
            if inliers_count>most_inliers_count
                
                inliers=inliers_temp;
                F=F_temp{j};
                most_inliers_count=inliers_count;
            end
        end
   end
end