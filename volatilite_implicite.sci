
exec('/home/amine/Documents/2A/Projet MOPSI/MOPSI-Project/sousjacent.sci', -1)


Call_market = 20;

function [Y]=f(t,T,K,r,sigma,x)
  Y = prix_call(t,T,K,r,sigma,x)-Call_market;
endfunction

function [Y]= vega(t,T,r,K,sigma,x)
  Y = x*sqrt(T)*dN(d1(t,T,K,r,sigma,x))
endfunction


function [Y]=volimpl(t,T,K,r,sgm,x)
  a1 = 0.001;
  a2 = sgm;
  tmp = (a1+a2)/2;
  i = 0;
  while i<300 then
    if f(t,T,K,r,a1,x)*f(t,T,K,r,tmp,x) > 0 then
      a1 = tmp;
    else
      a2 = tmp;
    end
    tmp = (a1+a2)/2;
    i = i+1;
  end
  Y = tmp;
endfunction

vi =volimpl(0,T,K(100),r,1,10);
disp(vi);

n = prod(size(K));


function [sigmaC]=bsimpvol(option,S,K,r,T,sigma0);
  function [Y]=difference(s);
    d1=-((log(K/S)-(r+1/2*s^2)*T)/(s*sqrt(T)));
    d2=-((log(K/S)-(r-1/2*s^2)*T)/(s*sqrt(T)));
    Y=segno*S*cdfnor('PQ',segno*d1,0,1)-segno*K*exp(-r*T)*cdfnor('PQ',segno*d2,0,1)-option;
  endfunction
  segno=1;
  [sigmaC,d,inf]=fsolve(sigma0,difference);
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

//cmap=hotcolormap(40);
//f=gcf();//figure courante
//f.color_map=cmap;

plot3d(tt,m,sg_imp);
xtitle('Nappe de volatilité implicite','Maturité T', 'Moneyness K/S0','volatilité implicite')
//plot(m,sg_imp(1,:));
