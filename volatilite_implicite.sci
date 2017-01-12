
exec('/home/amine/Documents/2A/Projet MOPSI/MOPSI-Project/code_lapeyre.sci', -1)


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
    while abs(f(t,T,K(k),r,x0,x))> eps,
      n = n+1;
      x0 = x0 -(f(t,T,K(k),r,x0,x))/vega(t,T,r,K(k),x0,x);
    end
    Y(1,k) = x0;
  end
endfunction

function [Y]=volatilite_implicite2(K)
  x0 = 0;
  y0 = 1;
  for i=1:50 do
    z = 0.5*(x0+y0);
    if f(0,1,K,0,z/(1-z),S0)<0 then
      x0 = z;
    else
      y0 = z;
    end
  end
  z = 0.5*(x0+y0);
  Y = z/(1-z);
endfunction

n = prod(size(K));
sgimpl2 = zeros(1,n);
for i=1:n
  sgimpl2(1,i) = volatilite_implicite2(K(1,i));
end

//plot(K,sgimpl2);

C_bs = zeros(1,n);
for i=1:n
  C_bs(i) = f(0,1,K(i),0,sgimpl2(i),S0);
  ;
end;
//plot(sgimpl2,C_bs);


Moneyness = zeros(1,M+1);
I = 4;
Maturity = zeros(1,I);
for k=1:M+1
  Moneyness(1,k) = K(k)/S0;
end

for k=1:I
  Maturity(1,k) = k;
end

//sg_impl = zeros(I,M+1);
//for j=[1:I]
//  Y = volatilite_implicite(0,j*I/5,K,0,0.3,S0);
//    sg_impl(j,1:M+1) = Y;
//  end


//plot3d(Maturity,Moneyness,sg_impl)

S = zeros(I,M+1)
//for j=[1:I] do
//  S(j,1:M+1) = fsolve(0.3,f)
//end

//plot(K,sg_impl);



function [sigmaC,sigmaP]=bsimpvol(option,S,K,r,T,sigma0);

// PURPOSE: Compute the implicit volatility for the Black and Scholes
//          model when the price of an option is known, for both
//          call option and put option
//------------------------------------------------------------------
// INPUT:
// * option  = option price
// * S       = price of the underlying asset
// * K       = strike price
// * r       = riskless interest rate
// * T       = time to maturity
// * sigma0  = starting guess value for the iterations
// -----------------------------------------------------------------
// OUTPUT:
// * sigmaC  = implied volatility of a call option
// * sigmaP  = implied volatility of a put option
// -----------------------------------------------------------------
// Francesco Menoncin (2010)

  function [Y]=difference(s);
    d1=-((log(K/S)-(r+1/2*s^2)*T)/(s*sqrt(T)));
    d2=-((log(K/S)-(r-1/2*s^2)*T)/(s*sqrt(T)));
    Y=segno*S*cdfnor('PQ',segno*d1,0,1)-segno*K*exp(-r*T)*cdfnor('PQ',segno*d2,0,1)-option;
  endfunction
  segno=1;
  [sigmaC,d,inf]=fsolve(sigma0,difference);
  //if inf==1 then disp('Call: good convergence'); else disp('Call: bad convergence'); end
  //segno=-1;
  [sigmaP,d,inf]=fsolve(sigma0,difference);
  //if inf==1 then disp('Put: good convergence'); else disp('Put: bad convergence'); end
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
