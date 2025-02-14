clear 
clc 

% test 

addpath(genpath('../src'))
addpath('../../subspace_MoM/test/')
addpath('../../finufft/matlab/')



n = 64;
L = 3;
c = 0.5; 
S = zeros(1,4);
S(1) = 31;
S(2) = 31;
S(3) = 31;
S(4) = 30;

load('vol_coef.mat')
load('rot_coef.mat')

vol_coef = vol_coef(:);
rot_coef = rot_coef(:);
a_lms = a(:);
b_pu = b(:);




d = n^2;


% prepare the volume coefficients
real_t_complex = get_vol_real_t_complex(L,S);
a_lms = real_t_complex*a_lms;
A_lms_tns = vol_coeffs_vec_t_tns(a_lms, L, S);


sphbes_zeros = zeros([L+1,max(S)]);
for l=0:L
    sphbes_zeros(l+1,:) = besselzero(l+.5,max(S),1);
end
if mod(n,2)==0
    k = ((-n/2):(n/2-1))/n;  
else
    k = ((-(n-1)/2):((n-1)/2))/n; 
end

[kx,ky] = meshgrid(k);
kx = kx(:); ky = ky(:); 


% form basis matrix
grids = [kx';ky';zeros(size(ky'))]; 
grids = grids';
[r,th,phi] = Cart2Sph(grids(:,1), grids(:,2), grids(:,3));

load('grid_info.mat')
% norm(r(:)-rs(:))
% norm(th(:)-ths(:))
% norm(phi(:)-phs(:))


n_vol = get_num_basis(L,S);
base_mat = zeros(n^2,n_vol);
rad_part = zeros(n^2,n_vol);
ang_part = zeros(n^2,n_vol);
idx = 1;

load('basis_funcs.mat')

for l=0:L

    Yl = sph_harmonics(l, th, phi);
    zl = sphbes_zeros(l+1,:);
    cl = sqrt(2/c^3)./abs(sphbes(l, zl, true));

    fl = sphbes(l, r*zl/c, false)*diag(cl); 
    Fl = fl;
    fl(r>c) = 0;
    Fl(r>c,:) = 0;


    for s=1:S(l+1)
        fls = fl(:,s);
        Fls = Fl(:,s);
        for m=-l:l
            Ylm = Yl(:,m+l+1);
            base_mat(:,idx) = (fls.*Ylm);
            rad_part(:,idx) = fls(:);
            ang_part(:,idx) = Ylm(:);
            idx = idx + 1;
            
        end
    end
end


figure
imagesc(abs(fl))
colorbar


% load('precomp_vol_basis.mat')
% norm(precomp_vol_basis(:)-base_mat(:))/norm(base_mat(:))
% norm(radial_part(:)-rad_part(:))/norm(radial_part(:))
% norm(angular_part(:)-ang_part(:))/norm(ang_part(:))
% err = precomp_vol_basis - base_mat;
% figure 
% imagesc(log10(abs(err)))
% colorbar 










 
