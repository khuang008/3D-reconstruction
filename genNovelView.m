function genNovelView
	addpath(genpath('.'));
	load('data/K.mat'); %intrinsic parameters K
	i1 = imread('data/i1.jpg');
	i2 = imread('data/i2.jpg');
    load('data/noisy_correspondences.mat');
    normalization_constant = max(size(i1));
    % compute the fundamental matrix
    [F,inliers] = ransacF(pts1, pts2, normalization_constant);
    R=eye(3);
    t=zeros(3,1);
    M1=K*[R t];
    % computer M2, assume intrinsic parameters is equal to camera 1
    M2=camera2(F,K,K,pts1,pts2);
    P=triangulate(M1, pts1, M2, pts2);
    [plane1 inliers1]=find_plane(P);
    
    remaining_points = P;
    remaining_points(:,inliers1)=[];
    
    [plane2, inliers2] = find_plane(remaining_points);    
    h=figure;
    set(h,'name','View 1 of Smith Hall');
    frame= drawNovelView(plane1, plane2, M1);
    imshow(frame);
    h=figure;
    set(h,'name','View 2 of Smith Hall');
    frame2= drawNovelView(plane1, plane2, M2);
    imshow(frame2);
    
    h=figure;
    set(h,'name','Novel view 1 of Smith Hall');
    %get the rotation matrix about x axis
    R1=[  1 0 0 ;
        0 cos(pi/6) -sin(pi/6);
        0 sin(pi/6)  cos(pi/6)];
    % get the translation 
    t1=[0;2;2];
    M3=K*[R1 t1];
    frame3= drawNovelView(plane1, plane2, M3);
    imshow(frame3);
    h=figure;
    set(h,'name','Novel view 2 of Smith Hall');
    %get the rotation matrix about y axis
    R2=[cos(-pi/12) 0 sin(-pi/12) ;
        0  1 0;
        -sin(-pi/12) 0 cos(-pi/12)];
    % get the translation
    t2=[2;0;2];
    
    M4=K*[R2 t2];
    frame= drawNovelView(plane1, plane2, M4);
    imshow(frame);

end

