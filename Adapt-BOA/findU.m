function [ u ] = findU(X, a, k, m, psize )
%FÝNDU Summary of this function goes here
%   Detailed explanation goes here

u = zeros(psize,1);
for j = 1:psize
    for i = 1:30
        if X(j,i) > a
            u(j) = u(j) + k.*(X(j,i) - a).^m;
        elseif (-a <= X(j,i)) && (X(j,i) <= a)
            u(j) = 0;
        elseif X(j,i) < -a
            u(j) = u(j) + k.*(-X(j,i) - a).^m;
        else
            ...
        end
    end
end

