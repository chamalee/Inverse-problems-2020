
close all, clear

n=32;
%Standard test phantom (Shepp-Logan)
f = phantom(n);
 xx=linspace(0,1,n);

%% Compute transform

% Choose angles
nTheta1 = 180;
nTheta2 = 90;
theta1 = [0:10:nTheta1-1];
theta2 = [0:nTheta2-1];

g1 = radon(f,theta1);
g2 = radon(f,theta2);

[srow1,scol1] = size(g1);
[srow2,scol2] = size(g2);

%% FBP

FBP1 = @(g,theta) iradon(g1,theta1,'linear','ram-lak',1,n);
FBP2 = @(g,theta) iradon(g2,theta2,'linear','ram-lak',1,n);

gNoisy1=g1+0.01*randn(size(g1))*max(abs(g1(:)));
gNoisy2=g2+0.01*randn(size(g2))*max(abs(g2(:)));

for iii = 2:18
    fBack1 = FBP1(gNoisy1(:,1:iii),theta1(10:iii)); 
    figure(5),
    imagesc(fBack1)
    axis equal, axis off, colormap gray
    title(['Angle: ' num2str(iii)])
%   pause(0.01)
end
% 
for iii = 2:nTheta2
    fBack2 = FBP2(gNoisy2(:,1:iii),theta2(1:iii));
    figure(6),
    imagesc(fBack2)
    axis equal, axis off, colormap gray
    title(['Angle: ' num2str(iii)])
%   pause(0.01)
end




