clear 
clc 



addpath(genpath('../src'))
addpath ~/finufft/matlab

N = 1e4;
SNR = 1; 
out = generate_mom(N,SNR);
save('out_N_1e4_SNR_1.mat','out')

N = 1e4;
SNR = 0.1; 
out = generate_mom(N,SNR);
save('out_N_1e4_SNR_0.1.mat','out')

N = 1e4;
SNR = 0.01; 
out = generate_mom(N,SNR);
save('out_N_1e4_SNR_0.01.mat','out')


function out = generate_mom(N,SNR)
rng(1)
%% generate ground truth volume
V = double(ReadMRC('emd_34948.map'));
V = V/norm(V(:));
n = 64;

n_vMF = 12;
k = 5;
M = 5;
mus = randn(3,n_vMF);
for i=1:n_vMF
    mus(:,i) = mus(:,i)/norm(mus(:,i));
end
w = ones(1,n_vMF)/n_vMF;
Rots = reject_sampling_vMF(N,mus,w,k,M);


% form moments 
params.s2 = 250;
params.s3 = 120;
params.r2_max = 250 ;
params.r3_max = 120 ;
params.tol2 = 1e-12;
params.tol3 = 1e-12;
params.sigma = 0;
params.nstd = get_estimated_std(V, SNR);
params.augsize = 5;
out = form_estimated_moments(V, n, Rots, params);

end