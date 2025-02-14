clear 
clc 

addpath(genpath('../src'))
% addpath('../../subspace_MoM/test/')
addpath ../../subspace_MoM/test/


load('vmf_params.mat')

mus = double(centers'); 
w = double(w_vmf);
k = double(kappa); 
% b = b(:);
P = 2;
view_coef = vMF_t_rot_coef(mus,w,k,P);

func = @(alpha, beta, gamma) 4*pi*get_vMF_density(mus,w,[sin(beta)*cos(alpha);sin(beta)*sin(alpha);cos(beta)],k);
func(1,2,3)

