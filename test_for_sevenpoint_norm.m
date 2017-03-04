clear;
close all;
addpath('lib');
load('data/clean_correspondences.mat');


i1=imread('data/i1.jpg');
i2=imread('data/i2.jpg');

F_7_clean = sevenpoint_norm(pts1(:,1:7), pts2(:,1:7), max(size(i1)));
disp('F_7_clean=');
disp(F_7_clean{1});
displayEpipolarF(i1, i2, F_7_clean{1});