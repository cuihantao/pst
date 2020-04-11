function[agg_exc,exc_list,exca_n]=eqexc(n_mac,n_exc,mac_list,agg_ind)
% purpose:to construct an equivalent exciter for a coherent area
% to be used with i_agg.m
%
% syntax :[agg_exc,exc_list,exca_n]=eqexc(n_mac,n_exc,mac_list,agg_ind) 
% input :
%        n_mac  :  machine  data 
%        n_exc  :  exciter data
%        mac_list: coherent machines
%        agg_ind  : aggregation indicator:
%                 : agg_ind=1 BMI
%                 : agg_ind=2 Xd'_w
%                 : default: BMI + K_aw (no input required)
% output: 
%        agg_exc :equivalent exciter constants
%        exc_list: coherent area exciter list
%        exca_n :number of exciters in coherent area
%<><><><><><><><><><><><><><><><><><><><><><><><><><><

exca_c=[];
macha_ce=[];
macha_n=length(mac_list);      %# of machines in coherent area
exc_bus=n_exc(:,2)';
flag=0;
% 
for i=1:macha_n,
  k=element(mac_list(i),exc_bus);
  if k ==1,
    flag=1;   % exciter present in area
    macha_ce=[macha_ce;n_mac(mac_list(i),:)]; % coherent area machine constant
    exca_c=[exca_c;n_exc(mac_list(i),:)];     % coherent area exciter constants
  end
end
[exc_list]=exca_c(:,2)
%keyboard
if flag ==1,      %  aggregate exciters
  % find machine with biggest H 
  bigm_h=0;
  [machae_n,dum]=size(macha_ce);
  for i=1:machae_n,   
    if bigm_h <= macha_ce(i,16),
      bigm_h=macha_ce(i,16);
      bigm_n=i;
    end
  end 
  agg_exc=exca_c(bigm_n,:);  % make eq. exc_con equal to that of
                            % machine with biggest H
  [exca_n,dum]=size(exca_c);
%keyboard
  if nargin==3,
    num_ka=0;num_ta=0;num_ke=0;num_te=0;num_kf=0;num_tf=0;
    den=0;
    for i= 1:exca_n
      num_ka=num_ka+exca_c(i,4)./macha_ce(i,7);
      num_ta=num_ka+exca_c(i,5)./macha_ce(i,7);
      num_ke=num_ka+exca_c(i,10)./macha_ce(i,7);
      num_te=num_ka+exca_c(i,11)./macha_ce(i,7);
      num_kf=num_ka+exca_c(i,16)./macha_ce(i,7);
      num_tf=num_ka+exca_c(i,17)./macha_ce(i,7);
      den=den+(1/macha_ce(i,7));
    end 
    agg_exc(:,4)=num_ka/den;
    agg_exc(:,5)=num_ta/den;
    agg_exc(:,10)=num_ke/den;
    agg_exc(:,11)=num_te/den;
    agg_exc(:,16)=num_kf/den;
    agg_exc(:,17)=num_tf/den;
    else
    if agg_ind == 2,          % BMI+k_aw aggregation
      num=0;
      den=0;
      for i= 1:exca_n
        num=num+exca_c(i,4)./macha_ce(i,7);
        den=den+(1/macha_ce(i,7));
      end 
      agg_exc(:,4)=num/den;   % weights K_a
    end
  end
  else
    agg_exc=[];
end
disp('end of eqexc')
keyboard
return
