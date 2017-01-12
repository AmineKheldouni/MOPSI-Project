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

r = 0;
T = 1;
S0 = 100;// Valeurs du sous jacent (l'actif) :

p = 2;// Nombre de modèles mélangés :
sigma = [0.2, 0.4];// volatilités sur un an dans le modèle de Black Scholes
proba = ones(1,p)/p;
Beta_opt = [100, 100, 0.2, 0.4, 1/2, 1/2];
// Prix d'un call : (alpha_k d'un Call)
M=20;
C = zeros(1, M+1);
K = zeros(1, M+1);
for i=[1:M+1] do
  K(i)=50+(i-1);
  C(i)=prix_melange("C",0,T,K(i),r,proba,sigma,S0*ones(1,p));
end

//plot(K,C);

//Beta0 = [100, 100, 0.2, 0.4, 0.5, 0.5];
Beta0 = [100, 100, 0.3, 0.5, 0.3, 0.7];
Beta_c = Beta0(1:prod(size(Beta0))-1);


function [Y]=der_prix_melange(type_operation, t, T, r, sigma, x, d)
  Y = zeros(1,M+1);
  if d>2*p then
    for k=1:M+1 do
        Y(1,k) = prix_call(t,T,K(k),r,sigma(modulo(d,p)+1),x(modulo(d,p)+1));
    end
  elseif d > p then
    for j=1:M+1 do
      for k=1:p do
          Y(1,j) =  Y(1,j) + proba(k) * prix_call_der(t,T,K(j),r,sigma(k),x(k),d);
      end
    end
  else
    for j=1:M+1 do
      for i=1:p do
        Y(1,j) =  Y(1,j) + proba(i) * prix_call_der(t,T,K(j),r,sigma(i),x(i),d);
      end
    end
  end
endfunction

function [Y]=grad_residus(Beta,C,K)
  Y = zeros(1,3*p);
  C_fct = zeros(1,M+1);
  for i=1:M+1 do
    C_fct(1,i) = prix_melange("C",0,T,K(i),r,Beta(2*p+1:3*p),Beta(p+1:2*p),Beta(1:p));
  end
  for i=1:3*p do
    Y(1,i) = 2*sum((C_fct-C).*der_prix_melange("C",0,T,r,Beta(p+1:2*p),Beta(1:p),i));
  end
endfunction

function [f,g,ind]=fct_objective(Beta_c,ind)
  // Gestion des contraintes sur les probabilités p_i
  Beta=[Beta_c,1-sum(Beta_c(2*p+1:3*p-1))];
  f=residus(Beta,C,K);
  g=grad_residus(Beta,C,K);
endfunction

Beta_inf = [30; 30; 0.1; 0.1; 0]
Beta_sup = [200; 200; 5; 5; 2]
//[fopt,xopt] = optim(fct_objective,"b", Beta_inf, Beta_sup, Beta_c)