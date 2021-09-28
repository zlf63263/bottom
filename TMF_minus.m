function [a_minus, r_c] = TMF_minus(u_s, Polarization)

global TotalLayers Thickness Lambda n_cplx EMLLocation z0

%%   Defining Fresnel coefficients, last indices are layer numbers.
%   ex: t(:,5,6) is t from 5th layer to 6th layer for all wavelengths.
%   There are (TotalLayers-EMLLocation) number of interfaces, 
%   from EML to the bottom air layer.
%   I, L matrices initialized following the same convention
%   I - 2 x 2 x Lambda x (i-th) x (i+1-th) matrix
%   L - 2 x 2 x Lambda x (i-th) matrix
t = zeros(size(Lambda,1),TotalLayers-1,TotalLayers);
r = zeros(size(Lambda,1),TotalLayers-1,TotalLayers);
Theta = zeros(size(Lambda,1),TotalLayers);
I = zeros(2,2,size(Lambda,1),TotalLayers-1,TotalLayers);
L = zeros(2,2,size(Lambda,1),TotalLayers-1);
L_total = zeros(2,2,size(Lambda,1),TotalLayers-1);
%M_downward = zeros(2,2,size(Lambda,1));
%Total_downward = zeros(2,2,size(Lambda,1));
a_minus = zeros(size(Lambda,1),1);
r_c = zeros(size(Lambda,1),1);

%%   Angle in the b-th layer, starting from top(sub) to bottom(air) side
for b = 1:TotalLayers
    Theta(:,b) = asin(u_s.*n_cplx(:,1)./n_cplx(:,b));
end

%   For loop for t_p,r_p,t_s,r_s, starting from top (1st) to the bottom (8th)
for b = 1:TotalLayers-1
    if Polarization == 'TM' % p-polarization
        % t(1,2) = 2n1cos01/(n1cos02+n2cos01)
        % r(1,2) = (n2cos01-n1cos02)/(n1cos02+n2cos01)
        t(:,b,b+1) = (2*n_cplx(:,b).*cos(Theta(:,b)))./(n_cplx(:,b).*cos(Theta(:,b+1))+n_cplx(:,b+1).*cos(Theta(:,b)));
        r(:,b,b+1) = (n_cplx(:,b+1).*cos(Theta(:,b))-n_cplx(:,b).*cos(Theta(:,b+1)))./(n_cplx(:,b).*cos(Theta(:,b+1))+n_cplx(:,b+1).*cos(Theta(:,b)));
    elseif Polarization == 'TE' % s-polarization
        % t(1,2) = 2n1cos01/(n1cos01+n2cos02)
        % r(1,2) = (n1cos01-n2cos02)/(n1cos01+n2cos02)
        t(:,b,b+1) = (2*n_cplx(:,b).*cos(Theta(:,b)))./(n_cplx(:,b).*cos(Theta(:,b))+n_cplx(:,b+1).*cos(Theta(:,b+1)));
        r(:,b,b+1) = (n_cplx(:,b).*cos(Theta(:,b))-n_cplx(:,b+1).*cos(Theta(:,b+1)))./(n_cplx(:,b).*cos(Theta(:,b))+n_cplx(:,b+1).*cos(Theta(:,b+1)));
    end
    
    % Calculate interface/propagation matrices
    I(1,1,:,b,b+1) = 1./t(:,b,b+1);
    I(1,2,:,b,b+1) = r(:,b,b+1)./t(:,b,b+1);
    I(2,1,:,b,b+1) = r(:,b,b+1)./t(:,b,b+1);
    I(2,2,:,b,b+1) = 1./t(:,b,b+1);
    
    L_total(1,1,:,b) = exp(-1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
    L_total(1,2,:,b) = 0;
    L_total(2,1,:,b) = 0;
    L_total(2,2,:,b) = exp(1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
    
    if b == EMLLocation
        L(2,2,:,b) = exp(1i*2*pi*n_cplx(:,b)./Lambda(:)*(Thickness(b)-z0).*cos(Theta(:,b)));
    end
    
    if b > EMLLocation
        L(1,1,:,b) = exp(-1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
        L(1,2,:,b) = 0;
        L(2,1,:,b) = 0;
        L(2,2,:,b) = exp(1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
    end
    
end

M_downward = I(:,:,:,EMLLocation,EMLLocation+1);
Total_downward = I(:,:,:,1,2);

%for b = EMLLocation+1:TotalLayers-1
for b = 2:TotalLayers-1
    for c = 1:size(Lambda,1)
        if b >= EMLLocation+1
            M_downward(:,:,c) = M_downward(:,:,c)*L(:,:,c,b)*I(:,:,c,b,b+1);
        end
        Total_downward(:,:,c) = Total_downward(:,:,c)*L_total(:,:,c,b)*I(:,:,c,b,b+1);
    end
end

a_minus(:,1) = M_downward(2,1,:)./M_downward(1,1,:).*(L(2,2,:,EMLLocation)).^2;
r_c(:,1) = Total_downward(2,1,:)./Total_downward(1,1,:);
