function [grad] = gradFid(fk,g)
%UNTITLED2 Summary of this function goes here
%  A'*(A*f-g)
nTheta = 180;
theta = [0:5:nTheta-1];
grad = At(A(fk,theta)-g,theta);
end

