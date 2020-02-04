clear all;close all;

% Construct a phantom of size NxN
N = 128;
f =phantom(N);

% Choose angles
nTheta = 180;
theta = [0:5:nTheta-1];

g = radon(f,theta);

% Add noise to the data
noiseLev=0.01;
g_noise = g + noiseLev*randn(size(g));
%% Projected Grad descent for Tikhonov

fk=zeros(size(f));
alpha=0.01;

lambda = 0.01;
for iii = 1: 100000
 f0= fk;   
 
fk =  max(fk - lambda *(gradFid(fk,g_noise) + alpha* 2*fk),0);

%  norm_diff= norm(f0-fk)
 if norm(f0-fk) < 1e-5
              break;
        end

end

% % Compute relative error
err_squ = norm(f(:)-fk(:))/norm(f(:));

figure(1)
clf
imagesc([f,fk],[0,1])
colormap gray
axis equal
axis off
title(['Tikhonov regularization error: ', num2str(round(err_squ*100)), '%, Number of steps: ',num2str((iii)) ])

    



