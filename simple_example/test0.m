clear 
clc 

addpath(genpath('../src'))
addpath('../../subspace_MoM/test/')


load('vmf_params.mat')

mus = double(centers'); 
w = double(w_vmf);
k = double(kappa); 
b = b(:);
P = 2;
view_coef = vMF_t_rot_coef(mus,w,k,P);

get_vMF_density(mus,w,[sin(1)*cos(2);sin(1)*sin(2);cos(1)],k)*4*pi 