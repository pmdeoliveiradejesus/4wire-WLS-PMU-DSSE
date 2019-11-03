%  pol2rec(r,theta) creates complex number r angle theta
% Written by S Williams 9/8/14
%University of Technology, Sydney
%Good for quick conversion of polar co-ordinates to complex number
% In complex form, 1 angle 45 is:
%a=pol2rec(1,45) % answer is 0.7071+.7071i (double)
%vpa(a) will then show more decimal places if required (symbolic)
    function [ans] = pol2rec(r,theta)
     ans=r*(exp(j*theta*pi/180));% this notation is in radians, facilitates conversion from degrees
    

    end
     