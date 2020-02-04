clear all;close all;

% Construct a phantom of size NxN
N = 128;
f =phantom(N);

% Choose angles
nTheta = 180;
theta = [0:5:nTheta-1];

g = radon(f,theta);
[srow1,scol1] = size(g);

% Add noise to the data
 g_noise = g(:)+ 0.001*max(abs(g(:)))*randn(size(g(:)));

f = f(:);
f = reshape(f,N,N);
Q_f= Q(f);

x = zeros(N^2,1);
b = At(reshape(g_noise,185,36),theta);

b = b(:);

x_opt = cg(f,b,x);
recn = reshape(x_opt,N,N);

% % Compute relative error
err_squ = norm(f(:)-recn(:))/norm(f(:));

figure(1)
clf
imagesc([f,recn],[0,1])
colormap gray
axis equal
axis off
title(['Tikhonov regularization: error ', num2str(round(err_squ*100)), '%'])
