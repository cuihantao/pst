function [n_bus,n_line,nmac_con] = ..
                       s_coh1(bus,line,area,nmach_a,basemva)
% Syntax   : [n_bus,n_line,nmac_con] = ..
%                      s_coh1(bus,line,area,nmach_a,basemva)
%
% Purpose  : To aggregate coherent machines using the
%            slow coherency method. In this scheme, only 
%            xdprime is corrected. In s_coh2, connection 
%            between the generator terminal buses are also
%            included for higher accuracy. A solved loadflow
%            input is required.
%
% Input    : bus      - bus data
%            line     - line data
%            area     - matrix of coherent machines
%            nmach_a  - vector of number of machines in each
%                         coherent area
%            basemva  - base mva (optional)
% Output   : n_bus  - new system bus data
%            n_line - new line data
%            nmac_con - aggregate generator data
%
% See also : podmore, i_agg, s_coh2
%
% Calls    :
%
% Call by  : 
%

% (c) Copyright 1991-2 Joe H. Chow - All Rights Reserved
%
% History ( in reverse chronologocal order )
%
% Version  : 1.0
% Author   : Joe H. Chow
% Date     : March 11, 1992

jay=sqrt(-1);
if nargin == 4
  basemva = basmva;         % from global variable
end

n_bus  = bus;              % create new bus, line, 
n_line = line;             %   and machine data
n_mac  = mac_con;
[nrow,ncol] = size(mac_con);
if ncol <= 21
  n_mac = [n_mac ones(nrow,23-ncol)];
end

% set up internal bus index vector
nbus = length(n_bus(:,1));
bus_int = zeros(round(max(n_bus(:,1))),1);
for i = 1:nbus
  bus_int(n_bus(i,1)) = i;
end
% set up internal generator index vector
tot_mac  = length(mac_con(:,1)); % total number of machines
mac_int = zeros(round(max(mac_con(:,1))),1);
for i = 1:tot_mac
  mac_int(mac_con(i,1)) = i;
end

bus_vol  = n_bus(:,2);        % system bus voltages
bus_ang  = n_bus(:,3);        % system bus angles
num_area = length(area(:,1)); % number of coherent areas

nmac_con = [];
for r=1:num_area,             % cycle thru all areas
  if nmach_a(r) == 1          % single machine area
    nmac_con = [nmac_con; n_mac(mac_int(area(r,1)),:)];
    n_mac(mac_int(area(r,1)),1) = 0;  % set machine # to 0
   else                  % areas with more than 1 machine
    num_mach = nmach_a(r);   % # of coherent machines
    com_bus  = max(n_bus(:,1)) + 1;    % common bus 
       % number, make to one higher than largest bus number
    mac_list  = mac_int(area(r,1:nmach_a(r)))';   
                                 % coherent machine numbers
    nlines   = length(n_line(:,1));    % # of lines
    bus_list = bus_int(mac_con(mac_list,2)); 
               % coherent machines bus numbers. It is okay
               %   to have identical buses in bus_list
    bus_type = 2;                     % generator bus
    for ii=1:num_mach,
      if n_bus(bus_list(ii,1),10) == 1,
        bus_type = 1;
      end
    end
    % inertia weighted aggregate machine internal voltage 
    %   and angle
    m = mac_con(mac_list,3).*mac_con(mac_list,16)/basemva;
    ma = sum(m);
    Pg = n_bus(bus_list,4).*mac_con(mac_list,22);
    Qg = n_bus(bus_list,5).*mac_con(mac_list,23);
    vb = n_bus(bus_list,2);
    angb = n_bus(bus_list,3)*pi/180;
    vbx = vb.*exp(jay*angb);
    xdp = basemva*mac_con(mac_list,7)./mac_con(mac_list,3);
    %  current and voltage computations include
    %  possibilities of multiple generator on the same bus
    int_cur = conj((Pg+jay*Qg)./vbx);
    int_vol = vbx + jay*xdp.*int_cur;
    % common bus (aggregate machine bus) voltage magnitude
    %   and angle
    mag_cbus = abs(int_vol)'*m/ma;   
    ang_cbus = angle(int_vol)'*m/ma;  
    vg = abs(int_vol);
    delta = angle(int_vol);
    % linearization
    K1 = -diag(vg.*vb.*cos(delta-angb)./xdp./m);
    K2 = [-diag(Pg./vb./m) -K1];

    % slow coherency transformation matrix
    T = [m'/ma; -ones(num_mach-1,1) eye(num_mach-1)]; 
                           % T = [C; G]

    % transform K1 and K2
    Tinv = inv(T);
    K1n = T*K1*Tinv;
    K2n = T*K2;

    % slow aggregate model
    invK1n22 = inv(K1n(2:num_mach,2:num_mach));
    K1a = ma*(K1n(1,1) - ..
              K1n(1,2:num_mach)*invK1n22*K1n(2:num_mach,1));
    K2a = ma*(K2n(1,:) - ..
              K1n(1,2:num_mach)*invK1n22*K2n(2:num_mach,:));
%keyboard
    % recover line parameters
    % aggregate generator
    beta1 = -K2a(1,1:num_mach);
    beta2 = K2a(1,num_mach+1:2*num_mach);
    phia = -atan2(beta1',beta2'./vb) + ..
                  ang_cbus*ones(angb) - angb;
    ra = sqrt( beta1'.^2 + (beta2'./vb).^2 );
%    tapa = mag_cbus*ones(vg)./vg;
%    xdpa = vg./ra;
    tapa = ones(vg);
    xdpa = mag_cbus./ra;
    v1 = mag_cbus*exp(jay*ang_cbus)*ones(vg);
    v2 = vb.*exp(jay*angb);
    [s1,s2] = line_pq(v1,v2,zeros(vg),xdpa,zeros(vg), ..
                      tapa,phia*180/pi);

    n_bus(com_bus,1)   = com_bus;     % add common bus to
                                      % bus data list
    n_bus(com_bus,2)   = mag_cbus;    % common bus voltage
    n_bus(com_bus,3)   = ang_cbus*180/pi; % angle
    n_bus(com_bus,4:9) = zeros(1,6);
    n_bus(com_bus,10)  = 3;           % bus type

    n_bus(bus_list,4:5)= zeros(num_mach,2); 
                                      % remove P & Q 
    for j = 1:num_mach
      n_bus(bus_list(j),6:7)= n_bus(bus_list(j),6:7) - ..
        [Pg(j)+real(s2(j)) Qg(j)+imag(s2(j))]; % adjust load
    end
    n_bus(bus_list,10) = ones(bus_list)*3;   
                                      % change load type

    rl = nlines + num_mach;
    n_line(nlines+1:rl,1) = ones(bus_list)*com_bus; 
                                          % from bus
    n_line(nlines+1:rl,2) = n_bus(bus_list,1); % to bus;
    n_line(nlines+1:rl,4) = xdpa;
    n_line(nlines+1:rl,6) = tapa; 
    n_line(nlines+1:rl,7) = phia*180/pi; 
    
    xdeq = 1/sum(ones(num_mach,1)./xdp);
    
    term_bus = max(n_bus(:,1))+1;     % new terminal bus # 
    com_vol = mag_cbus*exp(jay*ang_cbus);
    new_cur = conj(sum(s1)/com_vol);
    term_vol = com_vol -jay * xdeq * new_cur;

    n_bus(term_bus,1)   = term_bus;
    n_bus(term_bus,2)   = abs(term_vol);
    n_bus(term_bus,3)   = angle(term_vol)*180/pi;
    n_bus(term_bus,4)   = real(term_vol * conj(new_cur));
    n_bus(term_bus,5)   = imag(term_vol * conj(new_cur)); 
    n_bus(term_bus,6:9) = zeros(1,4); 
    n_bus(term_bus,10)  = bus_type;

    n_line(rl+1,1)      = com_bus;
    n_line(rl+1,2)      = term_bus;
    n_line(rl+1,3)      = 0.0;
    n_line(rl+1,4)      = -xdeq;
    n_line(rl+1,5:7)    = zeros(1,3);

    %  perform machine aggregation
    agg_mac = eqgen(n_mac,mac_list,basemva, ..
                                        term_bus,area(r,1));
    nmac_con = [nmac_con; agg_mac];
    n_mac(mac_list,1) = zeros(num_mach,1);  % set machine # 
                                            % to zero
  end
end

%  organize machine data
for i = 1:tot_mac
  if n_mac(i,1) ~= 0
    nmac_con = [nmac_con; n_mac(i,:)]
   end
end
