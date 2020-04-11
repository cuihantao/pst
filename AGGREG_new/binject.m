function [bus_bal] = binject(bus,line,bbus)
% Syntax:  [bus_bal] = binject(bus,line,bbus)
%
% Purpose: balance the boundary bus injections to 
%          obtain a solved loadflow.
%
% Input:   bus       - bus data
%          line      - line data
%          bbus      - boundary buses (optional). If not
%                      input, all buses are boundary buses. 
%
% Output:  bus_bal   - balanced loadflow bus data
%
% See also: island
%
% Calls:
%
% Call By: 

% (c) Copyright 1991 Joe H. Chow - All Rights Reserved
%
% History (in reverse chronological order)
%
% Version:   1.0
% Author:    Joe H. Chow
% Date:      August 1991
%
% ************************************************************
pst_var
jay = sqrt(-1);
if isempty(line) == 1
  bus_bal = bus;
  bus_bal(:,6) = bus(:,4) - bus(:,8).*bus(:,2).^2;
  bus_bal(:,7) = bus(:,5) + bus(:,9).*bus(:,2).^2; 
  return  
end
Y = ybus(bus,line);
V = bus(:,2); ang = bus(:,3)*pi/180;
% voltage in rectangular coordinate
V_rect = V.*(cos(ang)+jay*sin(ang));  
% bus current injection
cur_inj = Y*V_rect;
% power output 
S = V_rect.*conj(cur_inj);
P = real(S); Q = imag(S);

Pg = bus(:,4); Qg = bus(:,5);
if exist('bbus') == 1
    for i = 1:length(bbus)
      j = bus_int(bbus(i));
      bus(j,6) = Pg(j) - P(j);  % bus active power load
      bus(j,7) = Qg(j) - Q(j);  % bus reactive power load
    end
  else
    bus(:,6) = Pg - P;          % bus active power load
    bus(:,7) = Qg - Q;          % bus reactive power load
end
bus_bal = bus;
