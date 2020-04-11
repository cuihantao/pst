% svm_a1.m -- Joe Chow 8/91
% m file to generate a state space model of area 1 of
%   NPCC 16 machine system
clear
jay = sqrt(-1);

pst_var                                % set up global variables

data16m                                % load input data from m.file
basrad  = 2*pi*60;                     % system frequency is 60 Hz
basmva  = 100;                         % 100 MVA base
%syn_ref = 1;                          % absolute machine angle (zero reference)
syn_ref = 1;                           % machine 1 is reference

lfw = input('Do you need to recompute a loadflow solution? [y/n]   ','s');

if lfw == 'y',                         % solve for loadflow
                                       % loadflow parameter
  tol      = 1e-9;                     % tolerance for convergence
  iter_max = 30;                       % maximum number of iterations
  vmin     = 0.5;                      % voltage minimum
  vmax     = 1.5;                      % voltage maximum
  acc      = 1.0;                      % acceleration factor
  [bus_sol,line_flw] = ...
     loadflow(bus,line,tol,iter_max,vmin,vmax,acc,'n',2); 
  bus = bus_sol;                       % loadflow solution for initialization
  save data16m bus line mac_con asyn_damp exc_con bus_int
else
  load data16m
end
                                       % list of area 1 buses
bus_list = [ [1:1:29] [53:1:61] ];
[bus1,line1,bbus1,mac_con1,exc_con1] = ...
                            island(bus,line,bus_list);
bus1 = binject(bus1,line1,bbus1);
[nbus1 dum] = size(bus1);

while(0)
                                       % check loadflow
  tol      = 1e-9;                     % tolerance for convergence
  iter_max = 30;                       % maximum number of iterations
  vmin     = 0.5;                      % voltage minimum
  vmax     = 1.5;                      % voltage maximum
  acc      = 1.0;                      % acceleration factor

  bus1(nbus1,10) = 1;                  % change last bus to a swing bus
  [bus1_sol,line1_flw] = ...
     loadflow(bus1,line1,tol,iter_max,vmin,vmax,acc,'n',2); 
  bus1 = bus1_sol;                     % loadflow solution for initialization

nbbus = length(bbus1);                 % construct nc_load data
load_con1 = [ bbus1 zeros(nbbus,2) ones(nbbus,2)];
save d_16ma1 bus1 line1 mac_con1 exc_con1 load_con1 ...
              bus_list
end                                    % while(0)
load d_16ma1

bus = bus1; line = line1;              % rename variables
mac_con = mac_con1; exc_con = exc_con1; 
% exc_con = exc_con2;
load_con = load_con1;
[num_mach dummy] = size(mac_con);

                                       % simulation
stepsize = 0.01;                       % integration stepsize

% step 1: construct reduced Y matrix - load at bus 3 
%         modeled as constant PQ load
[Y11,Y12,Y21,Y22,rec_V1,rec_V2,bus_ord] = ...
                  red_Ybus(bus,line);  % pre-fault 
                                       % admittance matrix

                                       % step 2: initialization
flag = 0;
f = mac_tra(0,1,bus,flag);             % machine 1
f = exc_dc12(0,1,bus,flag,stepsize);   % machine 1 exciter
mach_ref(1) = 0;
V_nc = nc_load(bus,flag,Y22,Y21);
V_o = V_nc;
  y1 = abs  ( V_nc(1,1) );
  y2 = angle( V_nc(1,1) );
  y3 = abs  ( V_nc(2,1) );
  y4 = angle( V_nc(2,1) );
  u1 = real ( load_pot(1,2) );
  u2 = imag ( load_pot(1,2) );
  u3 = real ( load_pot(2,2) );
  u4 = imag ( load_pot(2,2) );
                                       % check for steady state initialization
                                       % network solution
  k = 1;
  flag = 1;
  f = mac_tra(0,k,bus,flag);           % machine 1
  f = exc_dc12(0,k,bus,flag,stepsize); % machine 1 exciter
  psi = psi_re(:,k) + jay*psi_im(:,k);
  V_nc = V_o;
  V_nc = nc_load(bus,flag,Y22,Y21,psi,V_nc,1e-10);
  cur = Y11*psi + Y12*V_nc;
  V_b = rec_V1*psi + rec_V2*V_nc;
  cur_re(:,k) = real(cur); cur_im(:,k) = imag(cur);
% v(bus_ord,1) = [V_nc; V_b];
                                       % compute dynamics
  pmech(:,k) = pmech(:,1);
  exc_sig(:,k) = zeros(9,1);
%  vex(1,k) = vex(1,1);
  flag = 2;
  f = mac_tra(0,k,bus,flag);           % machine 1
  f = exc_dc12(0,k,bus,flag,stepsize); % machine 1 exciter

 keyboard
                                       % end steady state initialization check

                                       % step 3: perform perturbation
k = 1;                                 % k is time step
pert_val = 0.001;                       % 1% perturbation
                                       % machine angle
for i = 1:num_mach                     % i is index on machine number
  j = i;
  V_nc = V_o;
  pert = pert_val*mac_ang(i,1);         % one percent perturbation
  if abs(pert) < .00001
     pert = 0.00001
  end
  nominal = mac_ang(i,1);
  mac_ang(i,1) = mac_ang(i,1) + pert;
  p_16ma1                              % m file of perturbations
  mac_ang(i,1) = nominal; 
end

                                       % machine speed
for i = 1:num_mach                     % i is index on machine number
  j = num_mach + i;;
  V_nc = V_o;
  pert = pert_val*mac_spd(i,1);         % one percent perturbation
  nominal = mac_spd(i,1);
  mac_spd(i,1) = mac_spd(i,1) + pert;
  p_16ma1                              % m file of perturbations
  mac_spd(i,1) = nominal; 
end

                                       % eqprime
for i = 1:num_mach                     % i is index on machine number
  j = 2*num_mach + i;
  V_nc = V_o;
  pert = pert_val*eqprime(i,1);        % one percent perturbation
  nominal = eqprime(i,1);
  eqprime(i,1) = eqprime(i,1) + pert;
  p_16ma1                              % m file of perturbations
  eqprime(i,1) = nominal; 
end

                                       % edprime
for i = 1:num_mach                     % i is index on machine number
  j = 3*num_mach + i;
  V_nc = V_o;
  pert = pert_val*edprime(i,1);        % one percent perturbation
  nominal = edprime(i,1);
  edprime(i,1) = edprime(i,1) + pert;
  p_16ma1                              % m file of perturbations
  edprime(i,1) = nominal; 
end

[num_exc dum] = size(exc_con);
                                       % Efd
for i = 1:num_exc                      % i is index on exciter number
  j = 4*num_mach + i;
  V_nc = V_o;
  pert = pert_val*Efd(i,1);            % one percent perturbation
  nominal = Efd(i,1);
  Efd(i,1) = Efd(i,1) + pert;
  p_16ma1                              % m file of perturbations
  Efd(i,1) = nominal; 
end

                                       % V_R
for i = 1:num_exc                      % i is index on exciter number
  j = 4*num_mach + num_exc + i;
  V_nc = V_o;
  pert = pert_val*V_R(i,1);            % one percent perturbation
  nominal = V_R(i,1);
  V_R(i,1) = V_R(i,1) + pert;
  p_16ma1                              % m file of perturbations
  V_R(i,1) = nominal; 
end

                                       % R_f
for i = 1:num_exc                      % i is index on exciter number
  j = 4*num_mach + 2*num_exc + i;
  V_nc = V_o;
  pert = pert_val*R_f(i,1);             % one percent perturbation
  nominal = R_f(i,1);
  R_f(i,1) = R_f(i,1) + pert;
  p_16ma1                              % m file of perturbations
  R_f(i,1) = nominal; 
end

if j == 63,
      j = 64
end

ii = 1;
      nominal = load_pot(1,2);         % perturbing active part of bus 1 current
      pert = pert_val * real(load_pot(1,2));
      load_pot(1,2) = (real(load_pot(1,2)) + pert) + imag(nominal) * jay;
      p_16ma1

ii = 2;
                                       % perturbing reactive part of bus 1 current
      pert = pert_val * imag(load_pot(1,2));
      load_pot(1,2) = real(nominal) + (imag(load_pot(1,2)) + pert) * jay;
      p_16ma1
      load_pot(1,2) = nominal;

ii = 3;
      nominal = load_pot(2,2);         % perturbing active part of bus 9 current
      pert = pert_val * real(load_pot(2,2));
      load_pot(2,2) = (real(load_pot(2,2)) + pert) + imag(nominal) * jay;
      p_16ma1

ii = 4;
                                       % perturbing reactive part of bus 9 current
      pert = pert_val * imag(load_pot(2,2));
      load_pot(2,2) = real(nominal) + (imag(load_pot(2,2)) + pert) * jay;
      p_16ma1
      load_pot(2,2) = nominal;
end
