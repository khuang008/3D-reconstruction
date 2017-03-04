clear;
close all;
addpath('lib');

load('data/noisy_correspondences.mat');

i1=imread('data/i1.jpg');
i2=imread('data/i2.jpg');
[F_7_RANSCA, inliers] = ransacF(pts1, pts2, max(size(i1)));
disp('F_7_RANSCA=');
disp(F_7_RANSCA);
displayEpipolarF(i1, i2, F_7_RANSCA)