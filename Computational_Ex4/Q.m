function [Qf] = Q(f)
alpha = 10;
nTheta = 180;
theta = [0:5:nTheta-1];
Qf = At(A(f,theta),theta) +  alpha*f;
end

