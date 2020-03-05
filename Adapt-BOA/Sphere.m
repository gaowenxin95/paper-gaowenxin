%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA102
% Project Title: Implementation of Particle Swarm Optimization in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

% function z=Sphere(x)
% 
%     z=sum(x.^2);
% 
% end


function z=Sphere(x)

     z=0;
     for i=1:dim
         z1=-20*exp(-0.2*sqrt(sum(x(i).^2)/dim))-exp(sum(cos(2*pi*x(i)))/dim)+20+exp(1);
         z=z+z1;
end



