function [a_plus, T_plus] = TMF_plus2(u_e, u_s, Polarization)

global Thickness Lambda n_cplx EMLLocation z0

%%   Defining Fresnel coefficients, last indices are layer numbers.
%   ex: t(:,5,4) is t from 5th layer to 4th layer for all wavelengths.
%   There are (EMLLocation-1) number of interfaces, 
%   from EML to the top substrate layer.
%   I, L matrices initialized following the same convention
%   I - 2 x 2 x Lambda x (i-th) x ('i-1'-th) matrix
%   L - 2 x 2 x Lambda x (i-th) matrix
t = zeros(size(Lambda,1),EMLLocation,EMLLocation-1);
r = zeros(size(Lambda,1),EMLLocation,EMLLocation-1);
Theta = zeros(size(Lambda,1),EMLLocation);
I = zeros(2,2,size(Lambda,1),EMLLocation,EMLLocation-1);
L = zeros(2,2,size(Lambda,1),EMLLocation);
a_plus = zeros(size(Lambda,1),1);
T_temp = zeros(size(Lambda,1),1);
T_plus = zeros(size(Lambda,1),1);

%%   Angle in the b-th layer, starting from EML to the top(subs) side
for b = 1:EMLLocation
    Theta(:,b) = asin(u_e.*n_cplx(:,EMLLocation)./n_cplx(:,b));
end

%   For loop for t_p,r_p,t_s,r_s, starting from EML (5) to the top (1)
for b = EMLLocation:-1:2
    if strcmp(Polarization,'TM') == 1 % p-polarization
        % t(2,1) = 2n2cos02/(n2cos01+n1cos02)
        % r(2,1) = (n1cos02-n2cos01)/(n2cos01+n1cos02)
        t(:,b,b-1) = (2*n_cplx(:,b).*cos(Theta(:,b)))./(n_cplx(:,b).*cos(Theta(:,b-1))+n_cplx(:,b-1).*cos(Theta(:,b)));
        r(:,b,b-1) = (n_cplx(:,b-1).*cos(Theta(:,b))-n_cplx(:,b).*cos(Theta(:,b-1)))./(n_cplx(:,b).*cos(Theta(:,b-1))+n_cplx(:,b-1).*cos(Theta(:,b)));
    elseif strcmp(Polarization,'TE') == 1 % s-polarization
        % t(2,1) = 2n2cos02/(n2cos02+n1cos01)
        % r(2,1) = (n2cos02-n1cos01)/(n2cos02+n1cos01)
        t(:,b,b-1) = (2*n_cplx(:,b).*cos(Theta(:,b)))./(n_cplx(:,b).*cos(Theta(:,b))+n_cplx(:,b-1).*cos(Theta(:,b-1)));
        r(:,b,b-1) = (n_cplx(:,b).*cos(Theta(:,b))-n_cplx(:,b-1).*cos(Theta(:,b-1)))./(n_cplx(:,b).*cos(Theta(:,b))+n_cplx(:,b-1).*cos(Theta(:,b-1)));
    end
    
    % Calculate interface/propagation matrices
    I(1,1,:,b,b-1) = 1./t(:,b,b-1);
    I(1,2,:,b,b-1) = r(:,b,b-1)./t(:,b,b-1);
    I(2,1,:,b,b-1) = r(:,b,b-1)./t(:,b,b-1);
    I(2,2,:,b,b-1) = 1./t(:,b,b-1);
    
    if b == EMLLocation
        L(2,2,:,b) = exp(1i*2*pi*n_cplx(:,b)./Lambda(:)*z0.*cos(Theta(:,b)));
    end
    
    if b < EMLLocation
        L(1,1,:,b) = exp(-1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
        L(1,2,:,b) = 0;
        L(2,1,:,b) = 0;
        L(2,2,:,b) = exp(1i*2*pi*n_cplx(:,b)./Lambda(:)*Thickness(b).*cos(Theta(:,b)));
    end
    
end

M_upward = I(:,:,:,EMLLocation,EMLLocation-1);

for b = EMLLocation-1:-1:2
    for c = 1:size(Lambda,1)
        M_upward(:,:,c) = M_upward(:,:,c)*L(:,:,c,b)*I(:,:,c,b,b-1);
    end
end

a_plus(:,1) = M_upward(2,1,:)./M_upward(1,1,:).*(L(2,2,:,EMLLocation)).^2;

if strcmp(Polarization,'TM') == 1 % p-polarization
    T_temp(:,1) = (abs(1./M_upward(1,1,:))).^2;
    %   PRB notation
    %T_plus(:,1) = T_temp(:,1).*((n_cplx(:,EMLLocation)./n_cplx(:,1))).^2.*(n_cplx(:,1).*sqrt(1-u_s(:,1).^2))./(n_cplx(:,EMLLocation)*sqrt(1-u_e^2));
    %   Hyunsu's notation
    %T_plus(:,1) = T_temp(:,1).*n_cplx(:,1).*sqrt(1-u_s(:,1).^2)./n_cplx(:,EMLLocation)./sqrt(1-(u_e)^2);
    
    %if u_s < 1
        T_plus(:,1) = T_temp(:,1).*n_cplx(:,1).*sqrt(1-u_s(:,1).^2)./n_cplx(:,EMLLocation)./sqrt(1-(u_e)^2);
    %else
    %    T_plus(:,1) = 0;
    %end
    
elseif strcmp(Polarization,'TE') == 1 % s-polarization
    T_temp(:,1) = (abs(1./M_upward(1,1,:))).^2;
    %   PRB notation
    %T_plus(:,1) = T_temp(:,1).*(n_cplx(:,EMLLocation)*sqrt(1-u_e^2))./((n_cplx(:,1)).*sqrt(1-u_s(:,1).^2));
    %   Hyunsu's notation
    %T_plus(:,1) = T_temp(:,1).*n_cplx(:,1).*sqrt(1-u_s(:,1).^2)./n_cplx(:,EMLLocation)./sqrt(1-(u_e)^2);
    
    %if u_s < 1
        T_plus(:,1) = T_temp(:,1).*n_cplx(:,1).*sqrt(1-u_s(:,1).^2)./n_cplx(:,EMLLocation)./sqrt(1-(u_e)^2);
    %else
    %    T_plus(:,1) = 0;
    %end
end
