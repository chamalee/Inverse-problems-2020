%% Radon transform
% Set up phantom
close all, clear

n=128;
%Standard test phantom (Shepp-Logan)
f_true = phantom(n);

figure(1),
imagesc(f_true),
axis equal, axis off
colormap gray
title('f true')

%% Compute transform

% Choose angles
nTheta = 180;
theta = [0:nTheta-1];

sinogram = radon(f_true,theta);
[row_sin,col_sin] = size(sinogram); % rows:185 (no of angles) % col:180

figure(2),
imagesc(sinogram),
colormap gray
title('Sinogram')

%% Backprojection
At = @(g,theta) iradon(g,theta,'linear','none',1,n);

%Let's backproject angle by angle
for iii = 2:nTheta
    fBack = At(sinogram(:,1:iii),theta(1:iii));
    figure(3),
    imagesc(fBack)
    axis equal, axis off, colormap gray
    title(['Back Projection--Angle: ' num2str(iii)])
%     pause(0.1)
end
[row_back,col_back] = size(fBack) % rows:128  % col:128

%%
At = @(g,theta) iradon(g,theta,'linear','none',1,n);
FBP = @(g,theta) iradon(g,theta,'linear','ram-lak',1,n);

for iii = 2:nTheta
    ffBack = FBP(sinogram(:,1:iii),theta(1:iii));
    figure(4),
    imagesc(ffBack)
    axis equal, axis off, colormap gray
    title(['Filtered Back Projection--Angle: ' num2str(iii)])
   % pause(0.01)
end
[row_fback,col_fback] = size(ffBack); % rows:128  % col:128

%%
gNoisy=sinogram+0.01*randn(size(sinogram))*max(abs(sinogram(:)));

FBP = @(g,theta) iradon(g,theta,'linear','ram-lak',1,n);
figure(5),
imagesc(gNoisy)
title('With noise')
figure(6),
thetaS = theta(1:5:180);
imagesc(FBP(gNoisy(:,1:5:180), thetaS))
axis equal, axis off
colormap gray
title('Filtered Back Projection: With noise')

