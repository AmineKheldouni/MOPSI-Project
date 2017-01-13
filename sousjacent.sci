
exec('/home/amine/Documents/2A/Projet MOPSI/MOPSI-Project/code_ameliore.sci', -1)

function [Y]=W_t(t)
  Y = grand(1,1,"nor",0,t);
endfunction


function [Y]= S_t(x,r,t,sigma)
  Y=x*exp((r-0.5*(sigma^2))*t+sigma*W_t(t));
endfunction

T_tot = 360;
temps = zeros(1,T_tot);
S = zeros(1,T_tot);

for i=1:T_tot do
  temps(1,i) = i*1/T_tot;
  S(1,i) = S_t(20,0,temps(1,i),0.5);
end

//plot(temps,S);
