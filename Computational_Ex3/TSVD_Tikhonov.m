clear all;close all;
% Construct a phantom of size NxN
N = 32;
ph =phantom(N);

% Choose angles
nTheta1 = 180;
nTheta2 = 90;
theta1 = [0:10:nTheta1-1];
theta2 = [0:5:nTheta2-1];

sinogram1 = radon(ph,theta1);
[srow1,scol1] = size(sinogram1);

sinogram2 = radon(ph,theta2);
[srow2,scol2] = size(sinogram2);

% Inversion by filtered back-projection. 
% FBP = iradon(sinogram1,theta1);

% Add noise to the data
g_noise1 = sinogram1(:)+0.01*max(abs(sinogram1(:)))*randn(size(sinogram1(:)));
g_noise2 = sinogram2(:)+0.01*max(abs(sinogram2(:)))*randn(size(sinogram2(:)));

%%  construct a matrix model 

% Initialize system matrix
A1 = zeros(srow1*scol1,N^2);

% Construct columns of A 
for iii = 1:N^2
    % Construct current "unit phantom"
    tmpim = zeros(N);
    tmpim(iii) = 1;
    
    % Compute sinogram
    tmpsino = radon(tmpim,theta1);
    
    % Insert the result into matrix
    A1(:,iii) = tmpsino(:);
end


% Initialize system matrix
A2 = zeros(srow1*scol1,N^2);

% Construct columns of A one by one
for iii = 1:N^2
    % Construct  "unit phantom"
    tmpim = zeros(N);
    tmpim(iii) = 1;
    
    % Compute sinogram
    tmpsino = radon(tmpim,theta2);
    
    % Insert the result into matrix
    A2(:,iii) = tmpsino(:);
end

%% Singular Value Decomposition
% A = U*D*V.'

[U1,D1,V1] = svd(A1);
svals1 = diag(D1);

[U2,D2,V2] = svd(A2);
svals2 = diag(D2);

%% Truncated SVD reconstruction

% Choose how many singular values to use
ind_alpha = 100;

% Reconstruct from data 
D_alpha1 = zeros(size(D1.'));
for jjj = 1:ind_alpha
    D_alpha1(jjj,jjj) = 1/svals1(jjj);
end
rec1_1 = V1*D_alpha1*U1.'*g_noise1;
rec1_1 = reshape(rec1_1,[N,N]);

% Reconstruct from data 
D_alpha2 = zeros(size(D2.'));
for jjj = 1:ind_alpha
    D_alpha2(jjj,jjj) = 1/svals1(jjj);
end
rec1_2 = V2*D_alpha2*U2.'*g_noise2;
rec1_2 = reshape(rec1_2,[N,N]);


figure(2)
% clf
 subplot(121)
imagesc([ph,rec1_1])
title(['TSVD : 0-180'],'fontsize',20)
axis equal
axis off
subplot(122)
imagesc([ph,rec1_2])
title(['TSVD : 0-90 '],'fontsize',20)
axis equal
axis off

%% Tikhonov regularization
% Goal is to solve the equation
%  ((A.')*A + alpha*I) f = (A.')*m

% Regularization parameter
alpha = 1;

% Compute regularized reconstruction using
% matrix inversion
rec2_1 = inv((A1'*A1 + alpha*eye(N^2)))*(A1')*g_noise1;
rec2_1 = reshape(rec2_1,[N,N]);

rec2_2 = inv((A2'*A2 + alpha*eye(N^2)))*(A2')*g_noise2;
rec2_2 = reshape(rec2_2,[N,N]);

err2_1=norm(rec2_1-ph)/norm(rec2_1);
% original&reconstruction
figure(3)
clf
imagesc([ph,rec2_1])
title(['Tikhonov  : 0-180'],'fontsize',20)
axis equal
axis off

figure(4)
clf
imagesc([ph,rec2_2])
title(['Tikhonov  : 0-90'],'fontsize',20)
axis equal
axis off






