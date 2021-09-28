function calculate(aniso, etlThick, emlLoc, ...
    dipoleLoc, angleMax, angleRes, quantumY, ...
    layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9, layer10, ...
    thick1, thick2, thick3, thick4, thick5, thick6, thick7, thick8, thick9, thick10, ...
    aInit, aFin, aStep, bInit, bFin, bStep, cInit, cFin, cStep, ...
    fig1On, fig2On, fig3On, wav1, fig4On, wav2 )
% To save angular profile simulation data, uncomment line near 342
% To save EQE data, uncomment line near 436
% To save spectrum vs. emission angle data, uncomment line near 425
% To save Power density vs. unit wavevector *To be Added*
% To save Power density polarizations vs. unit wavevector *To be Added*
% Check the spectrum file input near line 120 to match emitter. PL_perov is
% specifically for MAPbI3 emitter (peak near 760nm)

% Load library of materials
library = load('library-Diane.mat');

global TotalLayers Thickness Lambda n_cplx EMLLocation z0 n

% Set variable values from UI inputs
anisotropy = aniso; 
ETL_thickness = etlThick;

% Emissive layer location from the topmost layer (substrate) in stack
EMLLocation = emlLoc;
    
% Maximum in-plane wavevector value
aMax = angleMax;

% Resolution of in-plane wavevector sweep
AngleResolution = angleRes;

% Store EQE for each index2 iteration
EQEdrift = []; 

for i = 1:length(ETL_thickness)
    
    % Load real component of layer materials
    re1 = realComponent(layer1);
    re2 = realComponent(layer2);
    re3 = realComponent(layer3);
    re4 = realComponent(layer4);
    re5 = realComponent(layer5);
    re6 = realComponent(layer6);
    re7 = realComponent(layer7);
    re8 = realComponent(layer8);
    re9 = realComponent(layer9);
    re10 = realComponent(layer10);
    
    % Load imaginary component of layer materials
    im1 = imaginaryComponent(layer1);
    im2 = imaginaryComponent(layer2);
    im3 = imaginaryComponent(layer3);
    im4 = imaginaryComponent(layer4);
    im5 = imaginaryComponent(layer5);
    im6 = imaginaryComponent(layer6);
    im7 = imaginaryComponent(layer7);
    im8 = imaginaryComponent(layer8);
    im9 = imaginaryComponent(layer9);
    im10 = imaginaryComponent(layer10);
    
    % Layers from top to bottom
    n = [re1 re2 re3 re4 re5 re6 re7 re8 re9 re10];
    k = [im1 im2 im3 im4 im5 im6 im7 im8 im9 im10];
    
    % Set extinction coefficient of EML to zero
    k(:,EMLLocation) = zeros(451, 1);
    
    % Follow n+ik convention
    n_cplx = n+1i*k; 

    TotalLayers = size(n,2);
     
    %n_scan = n(1,1:EMLLocation);
    %n_smallest = min(n_scan);
    %n_smallest_location = find(n_scan == n_smallest);

    Lambda = transpose(colon(400,850));
    Spectrum = library.('PL_perov');
    decayingSpectra = zeros(6,length(Lambda)); %store the spectra generated for each angle from the norm, you MUST have the first index equal the number of angles you sample plus one (for the wavelengths), and the second index must be the length of the Lambda array.
    decayingSpectra(1,:) = Lambda;
    
    A_initial = aInit; A_final = aFin; A_step = aStep;
    B_initial = bInit; B_final = bFin; B_step = bStep;
    C_initial = cInit; C_final = cFin; C_step = cStep;
    
    F_prime = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    U_prime = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    n_out = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Portion_outcoupled = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Portion_substrate_trapped = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Portion_waveguided_new = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Portion_SPP = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Portion_absorption_new = zeros((A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);
    Spectrum_data = zeros(size(Lambda,1),(A_final-A_initial)/A_step+1,(B_final-B_initial)/B_step+1,(C_final-C_initial)/C_step+1);

    Radiance_angle = zeros(1,10);
    out = zeros(1,5);
    
    for param_C = C_initial:C_step:C_final
        for param_B = B_initial:B_step:B_final 
            for param_A = A_initial:A_step:A_final
                
                % Set layer thickness if value is a parameter
                thickness1 = updateParam(thick1, param_A, param_B, param_C);
                thickness2 = updateParam(thick2, param_A, param_B, param_C);
                thickness3 = updateParam(thick3, param_A, param_B, param_C);
                thickness4 = updateParam(thick4, param_A, param_B, param_C);
                thickness5 = updateParam(thick5, param_A, param_B, param_C);
                thickness6 = updateParam(thick6, param_A, param_B, param_C);
                thickness7 = updateParam(thick7, param_A, param_B, param_C);
                thickness8 = updateParam(thick8, param_A, param_B, param_C);
                thickness9 = updateParam(thick9, param_A, param_B, param_C);
                thickness10 = updateParam(thick10, param_A, param_B, param_C);
                
                % Create thickness array
                Thickness = [thickness1 thickness2 thickness3 thickness4...
                    thickness5 thickness6 thickness7 thickness8 thickness9 thickness10];
              
                % Set QY value 
                QY = updateParam(quantumY, param_A, param_B, param_C);
                
                % Set z0 value
                % Location of dipole sheet (in nm), from the top edge of the EML
                z0 = updateParam(dipoleLoc, param_A, param_B, param_C);

                a_TM_plus = zeros(size(Lambda,1),1); a_TE_plus = zeros(size(Lambda,1),1);
                a_TM_minus = zeros(size(Lambda,1),1); a_TE_minus = zeros(size(Lambda,1),1);
                T_TM_plus = zeros(size(Lambda,1),1); T_TE_plus = zeros(size(Lambda,1),1);
                K_TMv = zeros(size(Lambda,1),aMax); K_TMv_prime = zeros(size(Lambda,1),aMax);
                K_TMh = zeros(size(Lambda,1),aMax); K_TMh_prime = zeros(size(Lambda,1),aMax);
                K_TEh = zeros(size(Lambda,1),aMax); K_TEh_prime = zeros(size(Lambda,1),aMax);
                K = zeros(size(Lambda,1),aMax); K_temp = zeros(size(Lambda,1),aMax);
                K_prime = zeros(size(Lambda,1),aMax); K_prime_temp = zeros(size(Lambda,1),aMax);
                K_out = zeros(size(Lambda,1),aMax); K_out_temp = zeros(size(Lambda,1),aMax);
                U_temp = zeros(size(Lambda,1),aMax);
                K_sub = zeros(size(Lambda,1),aMax); K_SPP = zeros(size(Lambda,1),aMax); K_absorption = zeros(size(Lambda,1),aMax);
                K_waveguided = zeros(size(Lambda,1),aMax);
                t_TM_so = zeros(size(Lambda,1),aMax); r_TM_so = zeros(size(Lambda,1),aMax); r_TM_c = zeros(size(Lambda,1),aMax);
                t_TE_so = zeros(size(Lambda,1),aMax); r_TE_so = zeros(size(Lambda,1),aMax); r_TE_c = zeros(size(Lambda,1),aMax);
                T_TM_so = zeros(size(Lambda,1),aMax); R_TM_so = zeros(size(Lambda,1),aMax); R_TM_c = zeros(size(Lambda,1),aMax);
                T_TE_so = zeros(size(Lambda,1),aMax); R_TE_so = zeros(size(Lambda,1),aMax); R_TE_c = zeros(size(Lambda,1),aMax);
                TM_ThickSubstrate = zeros(size(Lambda,1),aMax); TE_ThickSubstrate = zeros(size(Lambda,1),aMax);
                Power_sub = zeros(size(Lambda,1),1); Power_SPP = zeros(size(Lambda,1),1); Power_absorption = zeros(size(Lambda,1),1);
                
                %  Main calculation

                % Create progress bar in output window
                f = uifigure;
                d = uiprogressdlg(f,'Title','Progress Bar',...
                'Message','Calculating','Cancelable','on');
                
                for a = 1:aMax
                    u_e = (a-1)/AngleResolution;   % u within EML, u=sin(theta_EML)
                    u_s = n_cplx(:,EMLLocation)./n_cplx(:,1).*u_e ;  % u within substrate
                    u_a = n_cplx(:,EMLLocation)./n_cplx(:,TotalLayers).*u_e;
                    [a_TM_plus, T_TM_plus] = TMF_plus2(u_e, u_s, 'TM');    % Calculates A4 and A12
                    [a_TE_plus, T_TE_plus] = TMF_plus2(u_e, u_s, 'TE');    % Calculates A4 and A13
                    [a_TM_minus, r_TM_c] = TMF_minus(u_s, 'TM');    % Calculates A5
                    [a_TE_minus, r_TE_c] = TMF_minus(u_s, 'TE');    % Calculates A5

                    % A1, A2, A3, A6
                    if a ~= AngleResolution+1   % to avoid NaN
                        K_TMv(:,a) = 3/4*real(u_e^2/sqrt(1-u_e^2).*(1+a_TM_plus(:,1)).*(1+a_TM_minus(:,1))./(1-a_TM_plus(:,1).*a_TM_minus(:,1)));
                        K_TMh(:,a) = 3/8*real(sqrt(1-u_e^2).*(1-a_TM_plus(:,1)).*(1-a_TM_minus(:,1))./(1-a_TM_plus(:,1).*a_TM_minus(:,1)));
                        K_TEh(:,a) = 3/8*real(1./sqrt(1-u_e^2).*(1+a_TE_plus(:,1)).*(1+a_TE_minus(:,1))./(1-a_TE_plus(:,1).*a_TE_minus(:,1)));
                        K(:,a) = anisotropy.*K_TMv(:,a)+(1-anisotropy)*K_TMh(:,a)+ (1-anisotropy).*K_TEh(:,a);
                        K_temp(:,a) = K(:,a)*u_e;
                    end

                    % A14 : calculating Fresnel coefficients for glass/air interface:
                    % (t_TM_so, r_TM_so, t_TE_so, r_TE_so)
                    % ## TM polarized mode: t(1,2) = 2n1cos01/(n1cos02+n2cos01), r(1,2) = (n2cos01-n1cos02)/(n1cos02+n2cos01)
                    %t_TM_so(:,a) = (2*n_cplx(:,1).*cos(asin(u_s(:,1))))./(n_cplx(:,1).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1)))+n_cplx(:,TotalLayers).*cos(asin(u_s(:,1))));
                    t_TM_so(:,a) = (2*n_cplx(:,1).*cos(asin(u_s(:,1))))./(n_cplx(:,1).*cos(asin(u_a(:,1)))+n_cplx(:,TotalLayers).*cos(asin(u_s(:,1))));
                    r_TM_so(:,a) = (n_cplx(:,TotalLayers).*cos(asin(u_s(:,1)))-n_cplx(:,1).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))))./(n_cplx(:,1).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1)))+n_cplx(:,TotalLayers).*cos(asin(u_s(:,1))));
                    % ## TE polarized mode: t(1,2) = 2n1cos01/(n1cos01+n2cos02), r(1,2) = (n1cos01-n2cos02)/(n1cos01+n2cos02)
                    t_TE_so(:,a) = (2*n_cplx(:,1).*cos(asin(u_s(:,1))))./(n_cplx(:,1).*cos(asin(u_s(:,1)))+n_cplx(:,TotalLayers).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))));
                    r_TE_so(:,a) = (n_cplx(:,1).*cos(asin(u_s(:,1)))-n_cplx(:,TotalLayers).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))))./(n_cplx(:,1).*cos(asin(u_s(:,1)))+n_cplx(:,TotalLayers).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))));
                    % Calculating transmittance/reflectance for glass/air interface:
                    T_TM_so(:,a) = real((n_cplx(:,TotalLayers).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))))./(n_cplx(:,1).*cos(asin(u_s(:,1)))).*(abs(t_TM_so(:,a)).^2));
                    R_TM_so(:,a) = abs(r_TM_so(:,a)).^2;
                    T_TE_so(:,a) = real((n_cplx(:,TotalLayers).*cos(asin(n_cplx(:,1)./n_cplx(:,TotalLayers).*u_s(:,1))))./(n_cplx(:,1).*cos(asin(u_s(:,1)))).*(abs(t_TE_so(:,a)).^2));
                    R_TE_so(:,a) = abs(r_TE_so(:,a)).^2;
                    % Calculating total reflection for OLED stacks from sub
                    % (R_TM_c, R_TE_c)
                    R_TM_c(:,a) = abs(r_TM_c(:,1)).^2;
                    R_TE_c(:,a) = abs(r_TE_c(:,1)).^2;
                    TM_ThickSubstrate(:,a) = T_TM_so(:,a)./(1-R_TM_so(:,a).*R_TM_c(:,a));
                    TE_ThickSubstrate(:,a) = T_TE_so(:,a)./(1-R_TE_so(:,a).*R_TE_c(:,a));

                    % A8, A9, A10, A11

                    if a ~= AngleResolution+1   % to avoid NaN
                        K_TMv_prime(:,a) = 3/8*u_e^2/sqrt(1-u_e^2).*(abs(1+a_TM_minus(:,1)).^2)./(abs(1-a_TM_plus(:,1).*a_TM_minus(:,1)).^2).*T_TM_plus(:,1);
                        K_TMh_prime(:,a) = 3/16*sqrt(1-u_e^2).*(abs(1-a_TM_minus(:,1)).^2)./(abs(1-a_TM_plus(:,1).*a_TM_minus(:,1)).^2).*T_TM_plus(:,1);
                        K_TEh_prime(:,a) = 3/16/sqrt(1-u_e^2).*(abs(1+a_TE_minus(:,1)).^2)./(abs(1-a_TE_plus(:,1).*a_TE_minus(:,1)).^2).*T_TE_plus(:,1);
                        K_prime(:,a) = anisotropy*K_TMv_prime(:,a)+ (1-anisotropy).*K_TMh_prime(:,a)+ (1-anisotropy).*K_TEh_prime(:,a);
                        K_prime_temp(:,a) = K_prime(:,a)*u_e;
                        K_out(:,a) = anisotropy*K_TMv_prime(:,a).*TM_ThickSubstrate(:,a)+ (1-anisotropy).*K_TMh_prime(:,a).*TM_ThickSubstrate(:,a)+ (1-anisotropy).*K_TEh_prime(:,a).*TE_ThickSubstrate(:,a);
                        K_out_temp(:,a) = ( anisotropy*K_TMv_prime(:,a).*TM_ThickSubstrate(:,a)+ (1-anisotropy).*K_TMh_prime(:,a).*TM_ThickSubstrate(:,a)+ (1-anisotropy).*K_TEh_prime(:,a).*TE_ThickSubstrate(:,a))*u_e;
                    end

                    % Check for 'cancel' button press
                    if d.CancelRequested
                        break
                    end

                    % Update progress bar
                    d.Value = a/aMax;
                end

                % Close progress bar and corresponding window
                close(d);
                delete(f);

                % Fill out AngleResolution+1-th term by averaging nearby terms
                K_temp(:,AngleResolution+1) = (K_temp(:,AngleResolution)+K_temp(:,AngleResolution+2))/2;
                K_prime_temp(:,AngleResolution+1) = (K_prime_temp(:,AngleResolution)+K_prime_temp(:,AngleResolution+2))/2;
                K_out(:,AngleResolution+1) = (K_out(:,AngleResolution)+K_out(:,AngleResolution+2))/2;
                K_out_temp(:,AngleResolution+1) = (K_out_temp(:,AngleResolution)+K_out_temp(:,AngleResolution+2))/2;
 
                % A7
                F = sum(K_temp,2)*((aMax-1)/AngleResolution)/(aMax-1)*2; 

                for b = 1:size(Lambda,1)
                    EML_to_air_crit = floor(n_cplx(b,TotalLayers)/n_cplx(b,EMLLocation)*AngleResolution);
                    U_temp(b,1:EML_to_air_crit) = K_out_temp(b,1:EML_to_air_crit);
                    EML_to_sub_crit = floor(real(n_cplx(b,1)/n_cplx(b,EMLLocation)*AngleResolution));
                    K_sub(b,1:EML_to_sub_crit) = K_prime_temp(b,1:EML_to_sub_crit);
                    K_waveguided(b,EML_to_sub_crit:AngleResolution) = K_temp(b,EML_to_sub_crit:AngleResolution);
                    K_SPP(b,1:AngleResolution+1) = K_temp(b,AngleResolution+1:aMax);
                    K_absorption(b,1:EML_to_sub_crit) = K_temp(b,1:EML_to_sub_crit)-K_prime_temp(b,1:EML_to_sub_crit);
                end
                

                U = sum(U_temp,2)*((aMax-1)/AngleResolution)/(aMax-1)*2;
                Power_sub(:,1) = sum(K_sub,2)*((aMax-1)/AngleResolution)/(aMax-1)*2;
                Power_SPP(:,1) = sum(K_SPP,2)*((aMax-1)/AngleResolution)/(aMax-1)*2;
                Power_absorption(:,1) = sum(K_absorption,2)*((aMax-1)/AngleResolution)/(aMax-1)*2;
                Power_waveguided(:,1) = sum(K_waveguided,2)*((aMax-1)/AngleResolution)/(aMax-1)*2;

                Portion_outcoupled((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*U./F);
                Portion_substrate_trapped((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*(Power_sub-U)./F);
                Portion_waveguided_new((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*Power_waveguided./F);
                Portion_SPP((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*Power_SPP./F);
                Portion_absorption_new((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*(F-Power_sub-Power_SPP-Power_waveguided)./F);

                F_prime((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*F);
                U_prime((param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = sum(Spectrum.*U);
                n_out_lambda = Spectrum.*U./F;

                n_EQE = Spectrum.*(QY*F)./(1-QY+QY*F).*U./F;
                EQE = sum(n_EQE,1);

                u_e_ext0 = floor(n_cplx(:,TotalLayers)*sin(0/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext10 = floor(n_cplx(:,TotalLayers)*sin(10/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext20 = floor(n_cplx(:,TotalLayers)*sin(20/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext30 = floor(n_cplx(:,TotalLayers)*sin(30/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext40 = floor(n_cplx(:,TotalLayers)*sin(40/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext50 = floor(n_cplx(:,TotalLayers)*sin(50/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext60 = floor(n_cplx(:,TotalLayers)*sin(60/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext70 = floor(n_cplx(:,TotalLayers)*sin(70/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext80 = floor(n_cplx(:,TotalLayers)*sin(80/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);
                u_e_ext90 = floor(n_cplx(:,TotalLayers)*sin(90/180*pi)./n_cplx(:,EMLLocation)*AngleResolution);

                I_0d = zeros(size(Lambda,1),1);
                I_10d = zeros(size(Lambda,1),1);
                I_20d = zeros(size(Lambda,1),1);
                I_30d = zeros(size(Lambda,1),1);
                I_40d = zeros(size(Lambda,1),1);
                I_50d = zeros(size(Lambda,1),1);
                I_60d = zeros(size(Lambda,1),1);
                I_70d = zeros(size(Lambda,1),1);
                I_80d = zeros(size(Lambda,1),1);

                for b = 1:size(Lambda,1)
                    I_0d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(0/180*pi)/pi.*K_out(b,u_e_ext0(b)+1)/F(b);
                    I_10d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(10/180*pi)/pi.*K_out(b,u_e_ext10(b)+1)/F(b);
                    I_20d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(20/180*pi)/pi.*K_out(b,u_e_ext20(b)+1)/F(b);
                    I_30d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(30/180*pi)/pi.*K_out(b,u_e_ext30(b)+1)/F(b);
                    I_40d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(40/180*pi)/pi.*K_out(b,u_e_ext40(b)+1)/F(b);
                    I_50d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(50/180*pi)/pi.*K_out(b,u_e_ext50(b)+1)/F(b);
                    I_60d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(60/180*pi)/pi.*K_out(b,u_e_ext60(b)+1)/F(b);
                    I_70d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(70/180*pi)/pi.*K_out(b,u_e_ext70(b)+1)/F(b);
                    I_80d(b,1) = 1240/Lambda(b)*Spectrum(b)*(QY*F(b))/(1-QY+QY*F(b))*(n_cplx(b, TotalLayers)./n_cplx(b,EMLLocation)).^2*cos(80/180*pi)/pi.*K_out(b,u_e_ext80(b)+1)/F(b);
                end

                % Plot spectrum v. viewing angle if selected on UI
                if( fig1On )
                    figure('Name', 'Spectrum v. Viewing Angle');
                    subplot(3,2,1)
                    plot(Lambda,I_0d(:,1)/sum(I_0d(:,1),1),Lambda,I_0d(:,1)/sum(I_0d(:,1),1));
                    legend('measured', 'simulated')
                    title('at 0d')
                    subplot(3,2,2)
                    decayingSpectra(2,:) = I_0d(:,1)/sum(I_0d(:,1));
                    
                    plot(Lambda,I_0d(:,1)/sum(I_0d(:,1),1),Lambda,I_20d(:,1)/sum(I_20d(:,1),1));
                    legend('measured', 'simulated')
                    title('at 20d')
                    subplot(3,2,3)
                    decayingSpectra(3,:) = I_20d(:,1)/sum(I_20d(:,1));
                    
                    plot(Lambda,I_0d(:,1)/sum(I_0d(:,1),1),Lambda,I_40d(:,1)/sum(I_40d(:,1),1));
                    legend('measured', 'simulated')
                    title('at 40d')
                    subplot(3,2,4)
                    decayingSpectra(4,:) = I_40d(:,1)/sum(I_40d(:,1));
                    
                    plot(Lambda,I_0d(:,1)/sum(I_0d(:,1),1),Lambda,I_60d(:,1)/sum(I_60d(:,1),1));
                    legend('measured', 'simulated')
                    title('at 60d')
                    subplot(3,2,5)
                    decayingSpectra(5,:) = I_60d(:,1)/sum(I_60d(:,1));
                    
                    plot(Lambda,I_0d(:,1)/sum(I_0d(:,1),1),Lambda,I_80d(:,1)/sum(I_80d(:,1)));
                    legend('measured', 'simulated')
                    title('at 80d')
                    subplot(3,2,6)
                    decayingSpectra(6,:) = I_80d(:,1)/sum(I_80d(:,1));
                    
                    plot(Lambda,I_0d(:,1), Lambda, I_20d(:,1), Lambda, I_40d(:,1), Lambda, I_60d(:,1), Lambda, I_80d(:,1));
                    legend('0d', '20d', '40d', '60d', '80d')
                    title('calculated I')
                end

                % Plot intensity v. viewing angle if selected on UI
                Radiance_angle(i,:) = [max(I_0d),max(I_10d),max(I_20d),max(I_30d),max(I_40d),max(I_50d),max(I_60d),max(I_70d),max(I_80d),0];
                Radiance_angle(i,:) = Radiance_angle(i,:)/max(Radiance_angle(i,:))
                if( fig2On )
                    figure('Name', 'Intensity v. Viewing Angle');
                    hold on;
                    theta_rad = linspace(0,pi/2,100);
                    theta = theta_rad*180/pi;
                    plot(theta,cos(theta_rad));
                    AngularIntensity = Radiance_angle(i,:)
                    %Uncomment the line below to save Angular Profile data,
                    %change "param_A" to the correct iterative variable. A
                    %list of 10 numbers, each point the normalized intensity at
                    %the angles 0-90 in steps of 10 degrees.
                    %save(strcat(int2str(param_A),'_ETLvsA.mat'),'AngularIntensity')
                    plot(0:10:90,AngularIntensity,'-o');
                    legend('Lambertian','Simulated');
                    hold off;
                end
                
                % Plot spectral power density if selected on UI
                if( fig3On )
                    figure('Name', 'Spectral Power Density');
                    hold on; 
                    
                    % Create normalized in-plane wavevector
                    u = 1:aMax;
                   
                    % Plot data on log axis
                    semilogy(u/AngleResolution, K_temp(wav1 - 400, u));
                    set(gca,'XLim',[0 aMax/AngleResolution]);
                    set(gca, 'YScale', 'log');
                    
                    % Set axes limits (optional)
                    xlim( [0 1.4] );
                    ylim([1e-2 1e2]);
                    
                    % Label axes
                    xlabel( 'Normalized in-plane wavevector' );
                    ylabel( 'Spectral power density' );
                    
                    % Mark regions of plot
                    xline( EML_to_air_crit/AngleResolution, '--' );
                    xline( EML_to_sub_crit/AngleResolution, '--' );
                    xline( AngleResolution/AngleResolution, '--' );
                    
                    hold off;
                end
                
                % Plot power dissipation spectrum if selected on UI
                if( fig4On )
                    figure('Name', 'Power Dissipation Spectrum');
                    
                    % Mark regions of plot
                    xline( EML_to_air_crit/AngleResolution, '--' );
                    xline( EML_to_sub_crit/AngleResolution, '--' );
                    xline( AngleResolution/AngleResolution, '--' );
                    
                    hold on;
                    
                    % Create normalized in-plane wavevector
                    u = 1:aMax;
                    
                    % Plot data on log axis
                    s1 = semilogy(u/AngleResolution, K_TMv(wav2 - 400, u));
                    s2 = semilogy(u/AngleResolution, K_TMh(wav2 - 400, u));
                    s3 = semilogy(u/AngleResolution, K_TEh(wav2 - 400, u));
                    set(gca,'XLim',[0 aMax/AngleResolution]);
                    set(gca, 'YScale', 'log')
                    
                    % Label axes
                    xlabel( 'Normalized in-plane wavevector' );
                    ylabel( 'Power Dissipation Spectrum' );
                    
                    % Set axes limits (optional)
                    % xlim( [0 1.4] );
                    ylim([1e-3 1e2]);
                    
                    hold off;

                    % Add legend
                    legend([s1 s2 s3], {'TMv', 'TMh', 'TEh'});
                    
                end    

                Spectrum_data(:,(param_A-A_initial+A_step)/A_step,(param_B-B_initial+B_step)/B_step,(param_C-C_initial+C_step)/C_step) = I_0d;

                param_A
                param_B
                param_C

                out(i,:) = [Portion_outcoupled, Portion_substrate_trapped, Portion_waveguided_new, Portion_SPP, Portion_absorption_new]
                EQEdrift = [EQEdrift;real(EQE)]
            end
        end
    end
    % Uncomment the line below to save spectrum vs. emission angle data
    % save(strcat(int2str(index2(i)),'decayingSpectrum.mat'),'decayingSpectra')
end

% POWER (Radiance/cossine)
angle = [(0:10:80)*pi/180 0];
angle = repmat(angle,i,1);

Power_angle = Radiance_angle./cos(angle);

% Uncomment the line below to save EQE data; a list of EQE values for
% each iteration of simulation is stored
% save(strcat(int2str(index2(i)),"_index_",int2str(10*param_C),'_QY_',int2str(param_B),'_Z0_','EQEvsIndex.mat'),"EQEdrift")