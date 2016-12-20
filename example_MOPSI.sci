exec("utils.sci");

//https://www.rocq.inria.fr/mathfi/Premia/

d=30;rho=0.0;
[mu,Gamma] = caracs_sans_actif_sans_risque(d,rho);

function [f,g,ind]=cost(x,ind)
   p=prod(size(x))+1;// dimension de lambda = 1+dim x
   lambda=[1-sum(x);x];// construction de lambda a partir de x
                       // en tenant compte de la contrainte sum(lambda)=1
   ps=mu*lambda;
   var=lambda'*Gamma*lambda;
   f=ps^2 / var;// fonction (lambda.mu)^2 / lambda'.Gamma.lambda
   k=(2*ps/var)*mu - (2*ps^2/var^2)*lambda'*Gamma;// derivee en fonction
de lambda
   g=k(2)-k(1);// derivee par rapport a x (en fonction de la derivee
en lambda)
   f=-f;g=-g;// On maximise mais Scicos suppose que l'on minimise ...
endfunction

x0=ones(d-1,1)/d;
[f,xopt]=optim(cost,x0);
Xopt=[1-sum(xopt);xopt];
Fopt=sqrt(-f)

// Est on bien entre 0 et 1 ?
// C'est toujours le cas pour Rho diagonale
// mais pas toujours dans le cas non diagonal
ok = and(0<=Xopt) & and(Xopt<=1)

// Lorsque rho=0, il y a une solution explicite (exercice)
// lambda_i = alpha * mu_i/sigma_i^2, renormalisé
// On verifie ...
sigma=sqrt(diag(Gamma))';
x=(mu ./ sigma^2);
x=x'/sum(x);// noramisation
norm(x - Xopt) // lorsque la matrice Rho est diagonale
                // ca devrait etre petit
