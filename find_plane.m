function [plane inliers] = find_plane(pts)
%%find a plane from the given 3D points by using RANSAC
% input:  
%   pts: 3 by N 3D points coordinates
%   
% output: 
%   plane: the estimation of fundamental matrix 
%   inliers: indices of the inliers of F .


% threshold for inliers 
threshold=0.01; 
num_of_iter=20000;          
inliers_temp=[];
inliers=[];
most_inliers_count=0;
N=size(pts,2);

for i = 1:num_of_iter        
        r=randperm(N,3);
        points=pts(:,r);
        plane_formular=get_plane_equation(points);
        
        inliers_count=0;
        inliers_temp=[];
        for j=1:N           
            dist=abs(distance_point_to_plane(pts(:,j),plane_formular));
            if(dist<threshold)
                    inliers_temp=[inliers_temp j];
                    inliers_count=inliers_count+1;
            end
        end      
        if inliers_count>most_inliers_count
                inliers=inliers_temp;
                plane=plane_formular;         
                most_inliers_count=inliers_count;
        end
end