function [G, symbol] = gaborfilter(scale, orientation,type);
% generate Gabor filter at fixed scale and orientation
% "G" is the Gabor pair
% "symbol" is the bar for display

expand = 12; b = floor(scale * expand+.5);  % b is the half size or radius
theta = (pi * orientation)/180;  % theta is orientation
Gauss = zeros(b+b+1,type);  % Gaussian function
Gcos = zeros(b+b+1,type); Gsin = zeros(b+b+1,type); % Gabor cos and sine pair
symbol = zeros(b+b+1,type); 
for x0 = -b : b
    for y0 = -b : b
        if (x0^2+y0^2>b^2) 
            fac = 0;  % zeros for pixels outside the circle
        else
            fac = 1; 
        end
        x = (x0 * cos(theta) + y0 * sin(theta))/scale; 
        y = (y0 * cos(theta) - x0 * sin(theta))/scale;
        g = exp(-(4*x^2+y^2)/100)/50/pi/scale^2; 
        Gauss(b+x0+1,b+y0+1) = g*fac; 
        Gcos(b+x0+1, b+y0+1) = g*cos(x)*fac;
        Gsin(b+x0+1, b+y0+1) = g*sin(x)*fac;
        symbol(b+x0+1, b+y0+1) = (abs(x)<b/2.5)*fac; % make a bar
    end
end
s = sum(Gauss(:)); sc = sum(Gcos(:)); r = sc/s; 
Gcos = Gcos - Gauss*r;   % mean is 0 by substracting DC component
Scos = sqrt(sum(sum(Gcos.^2))); Ssin = sqrt(sum(sum(Gsin.^2)));
Gcos = Gcos/Scos; Gsin = Gsin/Ssin; % l_2 norm is 1. 
G = Gcos + sqrt(-1)*Gsin;  % sine and cosine pair







