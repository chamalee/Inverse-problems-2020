
close all, clear,

%Number of evaluation points and grid
n=100;
xx = linspace(-1,1,n);
%% Define the convolution kernel

%we use a function handle (with @) to define a Gaussian kernel
sigma = 0.2; 
deltaX = xx(2)-xx(1); %grid size
convKer=@(x,sigma) deltaX/(sqrt(2*pi)*sigma)*exp(-(x).^2 / (2*sigma^2)  );

%plot function
figure(1), clf
plot(xx,convKer(xx,sigma),'linewidth',3)
% set(gca,'ylim',[0,1.1])
title('Gaussian vector')

%% Create convolution matrix
A = zeros(n,n);

for i = 1: n
    for j = 1 : n 
       A(i,j) =  deltaX/(sqrt(2*pi)*sigma)*exp(-(xx(i)-xx(j)).^2 / (2*sigma^2)  );
    end
end

%% Draw convolution matrix as an image
imagesc(A);
title('Convolution Matrix as an image')

%% SVD
[U,W,V] = svd(A);

norm(U*W*V'-A)

%% Pseudo inverse
Pseudo_inv = V*pinv(W)*U';

norm((pinv(W)*W)-eye(n))

norm((Pseudo_inv*A)-eye(n))

%%  First 9 columns of U , Last 9 columns of V , Singular values
%close all;
figure(3),clf
plot((1:n), U(:, 1:9));
title('First 9 columns of U')

figure(4),clf
plot((1:n), V(:,91:100));
title(' Last 9 columns of V')

rank_A = rank(A);
diag_vec = diag(W);
non_zero = diag_vec(1:rank_A);
figure(5),clf
plot((1:rank_A),log(non_zero));
title('Singular values')

