function [gf] = make_gabor(x,y,x0,y0, sigma, theta, lambda,varargin)

    if ~isempty(varargin)
        gamma = varargin{1};
    else
        gamma=1;
    end
    [xn,yn] = meshgrid(x-x0,y-y0);
    xnew = xn*cos(theta)+yn*sin(theta);
    ynew = -xn*sin(theta) + yn*cos(theta);
    
    gf = exp( -(xnew.^2+gamma^2*ynew.^2)/(2*sigma^2)).*exp(1i*2*pi*xnew/lambda);

end
