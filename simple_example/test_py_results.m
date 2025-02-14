clear 
clc 

% test 

addpath(genpath('../src'))
addpath('../../subspace_MoM/test/')
addpath('../../finufft/matlab/')


% N = 1e4; % sample size 
% SNR = 1; % signal-to-noise ratil 
n = 64; % size of downsampled images 
% use_precomputed_moments = false;

load('vol.mat')
load('vol_vec.mat')
vol = double(vol);


% rng(1)
% V = double(ReadMRC('emd_34948.map'));
% V = V/norm(V(:));
% % D = size(V,1);
% V_ds = vol_downsample(V, n);
% load('vol.mat')
% V_ds = double(vol);

L = 3; 
c = 0.5;

S = zeros(1,4);
S(1) = 31;
S(2) = 31;
S(3) = 31;
S(4) = 30;



vol_coef = vol_t_vol_coeffs(vol, L, S, c, 0);
load('vol_coef1.mat')
norm(vol_coef(:)-vol_coef1(:))









 
