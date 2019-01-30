function [D,mss,mssSlope] = getMss(points, lengthDivisor, micronPerPixel,secPerFrame);
% calculates the diffusion constant(s), the mss and the mssSlope.
% D(3) is the normal 2D-diffusion constant
% mssSlope is the slope of the BEST LINEAR FIT to the mss-values
 
% points: trajectory data in the format (npoints x 2)

% scales with O(n^2)
% THIS IS VERY SLOW, FIRST THINK IF YOU REALLY NEED TO CALL THIS FUNCTION!!

nMoments = 7;

l = length(points);
l3 = floor(l/lengthDivisor);
frameshifts = 1:1:l3;

% Dies ist ungef�hr O(n^2), also relativ teuer
meanMoments = zeros(nMoments,length(frameshifts));

for iF = 1:length(frameshifts) % verschiedene Fensterbreiten
    frameshift = frameshifts(iF);
    % neu anlegen f�r jeden Frameshift
    moments = zeros(nMoments,l-frameshift);

    dx = points(1:l-frameshift,1) - points(1+frameshift:l,1);
    dy = points(1:l-frameshift,2) - points(1+frameshift:l,2);    
    d_sq = (dx.*dx + dy.*dy)*micronPerPixel*micronPerPixel;
    
    for iMoment = 1:nMoments
       moments(iMoment,:) = d_sq(:).^((iMoment-1)/2);
    end        
    % mitteln �ber ganzen Track
    meanMoments(:,iF) = mean(moments,2);
end

% least square fit to log(meanMoments) vs. log(delta_t) for all orders
delta_t = frameshifts*secPerFrame;
for iMoment = 1:nMoments
    x = log10(delta_t);
    y = log10(meanMoments(iMoment,:));
    
    % find polynomial y = c(1)*x + c(2)
    c = polyfit(x,y,1);
    scalingCoefficients(iMoment) = c(1);
    intercepts(iMoment) = c(2);
end
mss = scalingCoefficients;
D = (10.^(intercepts))/4;

% find slope of moment scaling spectrum 
x = 0:(nMoments-1);
y = mss;

% CASE 1: (straight line through origin)

% options = optimset('TolFun',1.0e-12);
% a = fminsearch(@(a) myfun(a,x,y),0.5,options);
% mssSlope = a;


% CASE 2: LS fit to data
c = polyfit(x,y,1);    
mssSlope = c(1);



function f = myfun(a,x,y)
% squared difference of data to a straight line through the origin
f = sum((y-a*x).^2);