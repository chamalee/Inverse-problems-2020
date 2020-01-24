close all, clear,

%Number of evaluation points and grid
n=100; 
xx = linspace(-1,1,n);
%Init function as zeros
f = zeros(n,1);

%Create function
f(xx>-0.95 & xx<=-0.6) = 1;
f(xx>-0.6 & xx<=-0.2) = 0.2;
f(xx>-0.2 & xx<=0.2) = -0.5;
f(xx>0.4 & xx<=0.6) = 0.7;
f(xx>0.6 & xx<=1) = -0.7;

%plot function
figure(1), clf
plot(xx,f,'linewidth',3), hold on
set(gca,'ylim',[-1,1.1])
title('Signal f')
%% Define the convolution kernel

sigma_arry = [0.05,0.1,0.2];  % Sigma values


deltaX = xx(2)-xx(1); %grid size
plotIndex = 1;

G = zeros(length(sigma_arry),length(xx))

for sigmaIndex =1: length(sigma_arry)
%sigmaIndex = sigmaIndex + 1;
sigma = sigma_arry(sigmaIndex);
convKer=@(x,sigma) deltaX/(sqrt(2*pi)*sigma)*exp(-(x).^2 / (2*sigma^2)  );

for j = 1:length(xx)
G(sigmaIndex,j) = convKer( xx(j),sigma); 
end

A = zeros(length(sigma_arry),n,n);


for p = 1: n
    for q = 1 : n 
       A(sigmaIndex,p,q) =  deltaX/(sqrt(2*pi)*sigma)*exp(-(xx(p)-xx(q)).^2 / (2*sigma^2)  );
    end
end

A_sq= squeeze(A(sigmaIndex,:,:));
[U,W,V] = svd(A_sq);

rank_A = rank(A_sq)
diag_vec = diag(W);
non_zero = diag_vec(1:rank_A);
plotIndex = plotIndex+1;
figure(plotIndex), clf
plot((1:rank_A),log(non_zero));
title(['Singular values ( \sigma = ' num2str(sigma) ' )'])


g_conv = f'*A_sq;
plotIndex = plotIndex+1;
figure(plotIndex), clf
plot(xx,f,'linewidth',3), hold on
set(gca,'ylim',[-1,1.1])
plot(xx,g_conv,'linewidth',3)
title(['Convolution by matrix multiplication ( \sigma = ' num2str(sigma) ' )'])

%plot function
plotIndex = plotIndex+1;
figure(plotIndex), clf
plot(xx,convKer(xx,sigma),'linewidth',3)
title(['Kernel k ( \sigma = ' num2str(sigma) ' )'])

%% Convolution by multiplication in Fourier space
%close all;
g = conv(f,convKer(xx,sigma),'same');
%Fourier transformations of f and k
Ffx=fftshift(fft(fftshift(f)));
Fker=fftshift(fft(fftshift(convKer(xx,sigma))));

g_byF=(fftshift(ifft(fftshift(Ffx.*Fker.'))));
 
plotIndex = plotIndex+1;
figure(plotIndex), clf
plot(xx,f,'linewidth',3), hold on
plot(xx,g,'linewidth',3)
plot(xx,g_byF,'g--','linewidth',3);
set(gca,'ylim',[-1 1.1])
title(['Convolution by multiplication Fourier Space ( \sigma = ' num2str(sigma) ' )'])

end
 
