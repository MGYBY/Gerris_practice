clc;
clear;

%%%%% setup in Barker & Gray 2017 %%%%%
mu_d = 0.557;
mu_s = 0.342;
% d = 0.5*10^(-3);
i_node = 0.069;
phi_param = 0.60;
normal_depth = 4.0*10^(-3);
zeta = 0.4974188;%0.418; % 0.49742;
mu = tan(zeta);
mu_inf = 0.05;
gravity = 9.81;
grain_dia = 0.5*10^(-3);% 0.5*10^(-3); % 1.43*10^(-4); 
%i_alpha = i_node*((mu-mu_s)/(mu_s+(mu_d-mu_s)-mu));
%Reg i_alpha
%i_alpha = 0.4046
%% things could not match, don't know why

%%%%% another setup in Barker et al. 2020 %%%%%
mu_d = 0.557;
mu_s = 0.342;
i_node = 0.069;
phi_param = 0.60;
normal_depth = 5.0*10^(-3);
zeta = 0.41888;
mu = tan(zeta);
mu_inf = 0.05;
gravity = 9.81;
grain_dia = 1.43*10^(-4); % 0.5*10^(-3);

i_alpha = (mu-mu_d+sqrt(((mu_d-mu)^2)+4.0*i_node*mu_inf*(mu-mu_s)))/(2.0*mu_inf);
% bagnold_vel_max = (2.0)/(3.0)*i_alpha*sqrt(gravity*phi_param*grain_dia*cos(zeta)*(normal_depth/grain_dia)^3)*1.0;
bagnold_vel_max = (2.0)/(3.0)*i_alpha/grain_dia*sqrt(gravity*phi_param*cos(zeta))*((normal_depth)^(3.0/2.0));
froude = bagnold_vel_max/sqrt(gravity*cos(zeta)*normal_depth);

i_alpha_noReg = i_node*(mu-mu_s)/(mu_s+(mu_d-mu_s)-mu);
bagnold_vel_max_noReg = (2.0)/(3.0)*i_alpha_noReg/grain_dia*sqrt(gravity*phi_param*cos(zeta))*((normal_depth)^(3.0/2.0));
froude_noReg = bagnold_vel_max_noReg/sqrt(gravity*cos(zeta)*normal_depth);

