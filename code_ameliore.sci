exec('/home/amine/Documents/2A/Projet MOPSI/MOPSI-Project/fct_BS.sci', -1)


function valeur=residus(Beta_c,C,K)
//A Beta connu, renvoit Somme sur les strikes de (E(phi_k)-alpha_k)**2
  valeur=0;
  x0=Beta_c(1:p);
  sigma=Beta_c(p+1:2*p);
  proba=Beta_c(2*p+1:3*p);
  type_operation="C";
  C_fct = zeros(1,M+1)
  for i=1:p
    if Beta_c(2*p+i)<0 then
      valeur = 10000000000000000000
    elseif Beta_c(2*p+i)>1 then
      valeur = 10000000000000000000
    end
  end
  if valeur < 10000000000000000000 then
    for i=1:M+1
      C_fct(1,i) = prix_melange(type_operation,0,T,K(i),r,proba,sigma,x0)
    end
    valeur = sum((C_fct-C).^2);
  end
endfunction

//r = 0;
//T = 0.1;
r = 0.033;
T = 2;
//S0 = 100;// Valeurs du sous jacent (l'actif) :
//S0 = 630.15;
S0 = 505.15;

p = 2;// Nombre de modèles mélangés :
sigma = [0.2, 0.4];// volatilités sur un an dans le modèle de Black Scholes
proba = ones(1,p)/p;
Beta_opt = [S0, S0, sigma(1), sigma(2), proba(1), proba(2)];
// Prix d'un call : (alpha_k d'un Call)
M= 200;
C = zeros(1, M+1);
K = zeros(1, M+1);
for i=[1:M+1] do
  K(i)= 500+(i-1);
  C(i)=prix_melange("C",0,T,K(i),r,proba,sigma,S0*ones(1,p));
end

P = zeros(1, M+1);
for i=[1:M+1] do
  P(i)=prix_melange("P",0,T,K(i),r,proba,sigma,S0*ones(1,p));
end

clf();
subplot(211);
plot(K, C);
xtitle('Prix Call contre le strike K','K','Prix Call');

subplot(212);
plot(K, P);
xtitle('Prix Put contre le strike K','K','Prix Put');

//Beta0 = [100, 100, 0.2, 0.4, 0.5, 0.5];
Beta0 = [90, 90, 0.1,  0.5, 0.5, 0.5];
Beta_c = Beta0(1:prod(size(Beta0))-1);


function [Y]=grad_residus(Beta,C,K)
  Y = zeros(1,3*p-1);
  G = zeros(p+1,M+1);
  C_fct = zeros(1,M+1);
  for i=1:M+1 do
    C_fct(1,i) = prix_melange("C",0,T,K(i),r,Beta(2*p+1:3*p),Beta(p+1:2*p),Beta(1:p));
  end
  sigma = Beta(p+1:2*p);
  x = Beta(1:p);
  proba = Beta(2*p+1:3*p);
  for d=1:3*p-1 do
    if d>2*p then
        for j=1:M+1
          G(p+1,j) = 2*(C_fct(j)-C(j))*prix_call(0,T,K(j),r,sigma(modulo(d,p)),x(modulo(d,p)));
        end
        Y(d) = sum(G(p+1,:));
    elseif d>p then
      for i=1:p
        for j=1:M+1
          G(i,j) = 2*(C_fct(j)-C(j))*proba(i)*prix_call_der(0,T,K(j),r,sigma(i),x(i),d);
        end
      end
      Y(d) = sum(G);
    else
      for i=1:p
        for j=1:M+1
          G(i,j) = 2*(C_fct(j)-C(j))*proba(i)*prix_call_der(0,T,K(j),r,sigma(i),x(i),d);
        end
      end
      Y(d) = sum(G);
    end
    G = zeros(p+1,M+1);
  end
endfunction

function [f,g,ind]=fct_objective(Beta_c,ind)
  // Gestion des contraintes sur les probabilités p_i
  Beta=[Beta_c,1-sum(Beta_c(2*p+1:3*p-1))];
  f=residus(Beta,C,K);
  g=grad_residus(Beta,C,K);
endfunction

Beta_inf = [90, 90, 0.1, 0.1, 0.5]
Beta_sup = [110, 110, 0.4, 0.4, 0.5]
//[fopt,xopt] = optim(fct_objective,"b", Beta_inf, Beta_sup, Beta_c)

//disp(xopt)
//disp(fopt)
