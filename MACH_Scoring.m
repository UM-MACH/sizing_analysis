% ===========================================================
%                                         .:,:`           
%                                      `:,:               
%     ...       ...     `;;;          :,,     ;;;     ,;; 
%     ,,,:     ,,,,     ####        `,,:      ###     ;## 
%     ,,,,,   .,,,,     #####      `,,:       ###     ;## 
%     ,,,,,. `,,,,,    +##+##      ,,,        ###     ;## 
%     ,,,,,, ,,,,,,    ### ###    ,,,.        ###     ;## 
%     ,,, ,,,,, ,,,   .##` ###   ,,,,  #+     ########### 
%     ,,,  ,,,  ,,,   ###  `##;  *######+.,,,,########### 
%     ,,,  `,.  ,,,   ###   ###   ::,       .,,:.     :## 
%     ,,,   `   ,,,  ##########.  :,,         '`,,,,, ;## 
%     ,,,       ,,,  ###########   ,,,        ##@.,,,,, # 
%     ,,,       ,,, .##;     ###    ,,        ###  ,,,,,` 
%     ,,,       ,,: ###       ###    :,       ###  :,,,,, 
%                                   `:                 
% ===========================================================
%
% SUMMARY: Scoring analysis for the 2017-2018 DBF competition        

%% Initialization
% close all
clear all
% clc
%% =========================== Varibles ===================================== %%

span_wing = linspace(20,50,20);
passengers = linspace(0,40,20);
[passengers, span_wing] = meshgrid(passengers,span_wing);

g = 9.81;
num_wings = 2;

weight_fuselage_intial = 10;  % fuse empty  Newtons
% weight_fuselage = weight_fuselage_intial + passengers*0.0283495 * g; % Newtons
weight_fuselage = weight_fuselage_intial + (passengers*0.0283495 * g).^(1/2); % Newtons

% dens_lin_wing = 0.2875* 0.0283495 * g;  % density N/in

dens_lin_wing = 0.2875* 0.0283495 * g;  % density N/in

thrust_to_weight = .7;
Takeoff_velocity = 12; %m/s
CD_0 = 0.06;

e = 0.80;
lap_length = 5000; % ft
airfoil_Cl_max = 1.4;

%historical data
rho = 1.17; %air (kg/m^3)
% rho = 1.225; %air (kg/m^3)

%symbolic variables
v = sym('v','real');





% AR = ;  % inital guess for AR
Sref = 0.8; % inital guess for Sref
weight_propulsion = 12; % inital guess for Propulsion System weight in Newtons



%% ================ Data for Propulsion System Weight ======================= %%

thrust_data = [ 0.5, 1, 2, 3, 4, 5, 15]*4.44822; % Thrust needed in Newtons
weight_data =[9, 12, 25.9, 37.5, 42, 55, 113]*0.28; % Weight of prop. sys. in N

poly_coeff = polyfit(thrust_data,weight_data, 1 );


%% ========================= Size Aircraft ================================= %%
% Iterative method to solve for Sref, AR, weight_propulsion, MTOW, weight_empty
% Note: this method is not very sophisticated and prone to error 

for i = 0:1000

    % calculate the aspect ratio based on teh Sref
    AR = (span_wing.*0.0254).^2./Sref*num_wings;

    % wing weight in  Newtons 
    % weight_wings = num_wings*span_wing*dens_lin_wing.*(span_wing./AR /12).^0.585; 
    %weight_wings = num_wings*(span_wing*dens_lin_wing).*1.2072.*(span_wing./AR /12).^0.6076;
    weight_wings = num_wings*(span_wing*dens_lin_wing).*(span_wing./AR /12).^0.585;

    
    for k = 0:100

        weight_empty = weight_fuselage + weight_wings + weight_propulsion;
        MTOW = weight_empty + passengers*(1.3*0.0283495*g);
        
        thrust = thrust_to_weight*MTOW;

        err = sum(sum(abs(polyval(poly_coeff, thrust) - weight_propulsion)));
        % fprintf('err %f\n',err);

        if err < 1e-8
            % fprintf('thrust converged after %d\n',k);
            break

        end
        
       weight_propulsion = weight_propulsion + ...
        ( polyval(poly_coeff, thrust) - weight_propulsion) ;

    end
    
   
    Cl_stall = airfoil_Cl_max * AR ./ (AR + 2); % finite wing correction 
    Cl_takeoff = Cl_stall/(1.1^2); % equation from 481
    S_req = 2*MTOW./(Cl_takeoff.*rho*(Takeoff_velocity^2));


    err_sref = sum(sum(abs(Sref - S_req)));
    % fprintf('err  %d  %f\n ', i, err_sref);

    if err_sref < 1e-8
        fprintf('Sref converged after %d\n',i);
        break

    end

    Sref = Sref + 0.1*(S_req - Sref);  %m^2


    % % AR = (span_wing.*0.0254).^2./Sref*2;
    % % EW = (Weight_Fuselage + 1 * 0.0283495 * passengers) + (dens_lin_wing* span_wing).*(span_wing./AR / 12 * 0.5 + 0.5) + P;
    % % MTOW = (EW+passengers*2*0.0283495)*g


    % if sum(sum( abs( (span_wing.*0.0254).^2./Sref*2 - AR) )) < 1e-2
    %     disp('break')
    %     break
    % end
    % i
end
% return
AR( AR<= 0) = 0.1
AR(isnan(AR)==1) = 0.1
MTOW( MTOW<= 0) = 1e6
MTOW(isnan(MTOW)==1) = 0.1
Sref( Sref<= 0) = inf
Sref( isnan(Sref)== 1) = 0.01


% %% ======================= Calculations ================================ %%


% calculate Sref based on throwing speed 

% span_wingpan = (AR .* Sref).^0.5;
% Chord_root = Sref./span_wingpan/(1+lambda)*2;

% Chord = sqrt(Sref./AR); %meters
% span_wingpan = Sref./Chord; %wingspan in meters
%     Aspect_R = span_wingpan./Chord


T_0 = thrust_to_weight*MTOW;
T_1 = -0.060;
T_2 = -0.015;


for j = [1:size(passengers,1)]
    j
    for i = [1:size(span_wing,2)]
        v_cruise(j,i) = double(max(vpasolve( v.^2*T_2 + v*T_1 + T_0(j,i) ==...
        0.5*rho*v.^2*Sref(j,i).*...
        (CD_0 + 1./(pi*e*AR(j,i))*(2*MTOW(j,i)./(rho*v^2*Sref(j,i))).^2)...
        , v)));

        if imag(v_cruise(j,i)) || real(v_cruise(j,i)) < 0
            v_cruise(j,i) = 0;
%             disp('yes')
        end
        
    end
end

laps = (v_cruise)./(lap_length/3.28) * 60 * 10;
t_3laps = 3*(lap_length/3.28)./v_cruise*0.8;
RAC = weight_empty*0.224809.*span_wing;



%only account for well define aircraft
for j = [1:size(passengers,1)]
    for i = [1:size(span_wing,2)]
        if(j*3-45 > i)
            t_3laps(j,i) = 100000000;
            laps(j,i) = 0;
            
        end
    end
end        

v_cruise(v_cruise>25) = 25;


M1 = zeros(size(passengers));
M1(t_3laps/60 < 5) = 1;

M2 = 2*passengers./(t_3laps)./max(max(passengers./(t_3laps))); %mission 2 score
M3 =4*(passengers/2.*(passengers*2/2).*laps)./max(max((passengers/2.*(passengers*2/2).*laps))) + 2;

score = (M1 + M2 + M3)./RAC;
% score(M2 == 0) = 0;
% maximum = max(max(score));
% score = score./maximum;
% %% ========================== Plotting ================================= %%

figure
shading interp;
[ANALY2 , ANALY2] = contourf( passengers, span_wing, score);
set(ANALY2,'edgecolor','none');

title('Scoring Analysis','FontSize',23);


ylabel('b_s','FontSize',20);
xlabel('passengers','FontSize',20);
hold on 

max_score = max(max(score));
scatter(passengers(score== max_score), span_wing(score== max_score),250, 'ko','filled')

c = colorbar;
c.Label.String = 'Normalized Score';

fprintf('Max Score %.2f\n',max_score);
fprintf('RAC %.2f\n', RAC(score== max_score));
fprintf('passengers %.2f\n', passengers(score== max_score));
fprintf('span_wing [in]%.2f\n', span_wing(score== max_score));
fprintf('MTOW [N]%.2f\n', MTOW(score== max_score));
fprintf('Sref [m^2]%.2f\n', Sref(score== max_score));
fprintf('span_wingpan [in] %.2f\n', span_wing(score== max_score));
fprintf('Chord [in] %.2f\n', span_wing(score== max_score)/AR(score== max_score));
fprintf('AR %.2f\n', AR(score== max_score));
fprintf('V %.2f\n', v_cruise(score== max_score));
fprintf('Cl %.2f\n', Cl_takeoff(score== max_score));




figure
contourf( passengers, span_wing, RAC);
title('RAC','FontSize',23);
ylabel('b_s','FontSize',20);
xlabel('passengers','FontSize',20);
hold on 
c = colorbar;


% figure
% contourf( passengers, span_wing, (M1 + M2 + M3));
% title('(M1 + M2 + M3)','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, MTOW);
% title('MTOW','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, weight_empty);
% title('weight_empty','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, AR);
% title('AR','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, Sref);
% title('s','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;


% figure
% contourf( passengers, span_wing, v_cruise);
% title('v','FontSize',23);
% ylabel('b_s','FontSize',20);
% xlabel('passengers','FontSize',20);
% hold on 
% c = colorbar;
