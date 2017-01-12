
exec('/home/amine/Documents/2A/Projet MOPSI/MOPSI-Project/sousjacent.sci', -1)


Call_market = 32.4;

function [Y]=f(t,T,K,r,sigma,x)
  Y = prix_call(t,T,K,r,sigma,x)-Call_market;
endfunction

function [Y]= vega(t,T,r,K,sigma,x)
  Y = x*sqrt(T)*dN(d1(t,T,K,r,sigma,x))
endfunction


function [Y]=volatilite_implicite(t,T,K,r,sigma,x)
  Y = zeros(1,M+1);
  eps = 10^(-5);
  for k=1:M+1 do
    x0 = sigma;
    n = 0;
    while abs(f(t,T,K(k),r,x0,x))> eps
      n = n+1;
      disp(k);
      disp(x0);
      x0 = x0 -(f(t,T,K(k),r,x0,x))/vega(t,T,r,K(k),x0,x);
    end
    Y(1,k) = x0;
  end
endfunction


n = prod(size(K));


function [sigmaC,sigmaP]=bsimpvol(option,S,K,r,T,sigma0);
  function [Y]=difference(s);
    d1=-((log(K/S)-(r+1/2*s^2)*T)/(s*sqrt(T)));
    d2=-((log(K/S)-(r-1/2*s^2)*T)/(s*sqrt(T)));
    Y=segno*S*cdfnor('PQ',segno*d1,0,1)-segno*K*exp(-r*T)*cdfnor('PQ',segno*d2,0,1)-option;
  endfunction
  segno=1;
  [sigmaC,d,inf]=fsolve(sigma0,difference);
  [sigmaP,d,inf]=fsolve(sigma0,difference);
endfunction

pas = 100;
m = zeros(1,M+1);
tt = zeros(1, pas);

for i=1:M+1
  m(1,i) = K(i)/S0;
end

for i=1:pas
  tt(1,i) = i*T/pas;
end

sg_imp = zeros(pas,M+1);

for i=1:pas
  for j=1:M+1
    sg_imp(i,j) = bsimpvol(Call_market,S0,K(j),r,tt(1,i),0.2)
  end
end

cmap=hotcolormap(40);
f=gcf();//figure courante
f.color_map=cmap;

plot3d(tt,m,sg_imp);
//plot(m,sg_imp(1,:));
