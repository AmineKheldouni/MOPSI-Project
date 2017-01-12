funcprot(0);

function [Y]=d1(t,T,K,r,sigma,x)
  if t==T
   if x>K then Y=9999999;else Y=-9999999;end
  else
   Y=(log(x/K)+(r+(sigma^2)/2)*(T-t))/(sigma*sqrt(T-t));
  end
endfunction

function [Y]=Dd1(t,T,K,r,sigma,x,d)
   if d>2*p then Y = 0
   elseif d>p then
       //Y = -(log(x/K)+(r+(sigma^2)/2)*t)/(sqrt(T-t)*sigma^2)
       //+((log(x/K)+sigma*t)/(sigma*sqrt(T-t)));
       Y = -d1(t,T,K,r,sigma,x)/sigma+sqrt(T-t);
   else
       Y = 1/(x*sigma*sqrt(T-t));
   end
endfunction


function [y]=N(x)
// Fonction de répartition de la loi normale
  [y,Q]=cdfnor("PQ",x,0,1);
endfunction

// Fonction dfct_normale :
function [Y]=dN(x)
    Y = exp(-0.5*x.^2)/sqrt(2*%pi)
endfunction

function [Y]= prix_call(t,T,K,r,sigma,x)
// Prix du call dans la formule de Black et Scholes
    d_1=d1(t,T,K,r,sigma,x);
    d_2=d_1 - sigma*sqrt(T-t);
    Y=x*N(d_1)-K*exp(-r*(T-t))*N(d_2);
endfunction

function [Y]= prix_call_der(t, T, K, r, sigma, x, d)
  d_1 = d1(t,T,K,r,sigma,x)
  if d > 2*p then Y = 0
  elseif d > p then
    Y = x*Dd1(t,T,K,r,sigma,x,d)*dN(d_1)-K*exp(-r*(T-t))*(Dd1(t,T,K,r,sigma,x,d)-sqrt(T-t))*dN(d_1 - sigma*sqrt(T-t));
  else
    Y = N(d_1)+x*Dd1(t,T,K,r,sigma,x,d)*dN(d_1) - K*exp(-r*(T-t))*Dd1(t,T,K,r,sigma,x,d)*dN(d_1 - sigma*sqrt(T-t));
  end
endfunction

function [Y]= prix_put(t,T,K,r,sigma,x)
// Prix du put dans la formule de Black et Scholes
    d_1=d1(t,T,K,r,sigma,x);
    d_2=d_1 - sigma*sqrt(T-t);
    Y=-x*N(-d_1)+K*exp(-r*(T-t))*N(-d_2);
endfunction

function [Y]= prix_call_melange(t,T,K,r,proba,sigma,x)
// Prix du call dans la formule de Black et Scholes
  Y=0;
  for i=[1:prod(size(sigma))] do
    Y =  Y + proba(i) * prix_call(t,T,K,r,sigma(i),x(i))
  end
endfunction

function [Y]= prix_put_melange(t,T,K,r,proba,sigma,x)
// Prix du call dans la formule de Black et Scholes
  Y=0;
  for i=[1:prod(size(sigma))] do
    Y =  Y + proba(i) * prix_put(t,T,K,r,sigma(i),x(i));
  end
endfunction

function [Y]= prix_melange(type_operation, t,T,K,r,proba,sigma,x)
// Prix du call dans la formule de Black et Scholes
  if type_operation == "C" then
    Y=prix_call_melange(t,T,K,r,proba,sigma,x);
  else // probablement un put
    Y=prix_put_melange(t,T,K,r,proba,sigma,x);
  end
endfunction
