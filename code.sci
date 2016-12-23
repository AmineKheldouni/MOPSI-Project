///////////////////////////// Définiion des variables ///////////////////////
funcprot(0);
r = 0;
t = 1;
// Nombre de modèles mélangés :
p = 2;
M = 10;
// Valeurs du sous jacent (l'actif) :
lbd = [100, 100];


// Valeurs des écarts types sur un an dans le modèle de Black Scholes :
// (Les volatilités)
sg = [0.2, 0.4];


// Probabilité de la variable aléatoire N :
p1 = 1/2;
p_n = [p1, 1-p1];

// Simulation d'une loi Gaussienne centrée réduite :
G = grand(1, p, 'nor', 0, 1);

// Put ou Call 
type_operation = "C"

// Tirage aléatoire d'une variable uniforme avec probabilité p_n :
function n_valeur=N(liste_p)
    alpha = grand(1,1, 'def');
    l = cumsum(liste_p)
    n_valeur=min(find(alpha<l))*max(liste_p)
endfunction

// Définition du paramètre beta à retrouver :
Beta_opt = zeros(1, 3*p);

for i=1:p
    Beta_opt(1,i)= lbd(i);
    Beta_opt(1,i+p)= sg(i);
    Beta_opt(1,i+2*p)= p_n(i);
end

disp("Beta* = ")
disp(Beta_opt)

S0 = 0
liste_alea = zeros(1, p)
liste_gauss = zeros(1,p)
for i=1:p
   liste_alea(1,i) = N(p_n)
   liste_gauss(1,i) = G(i)
end
for i=1:p
    S0 = S0 + liste_alea(1,i)*Beta_opt(i)*exp(Beta_opt(i+p)*liste_gauss(1,i));
end
//S0 = 100
disp("S0 = ")
disp(S0)
// Création de la liste des strikes K_k:
K = zeros(1,M+1);
K(1,1)=50;
for i=2:M+1
    K(1,i)= K(1,i-1)+1;
end



// Création de la fonction pour le calcul des Payoffs :
function phi_k=payoff(x,k,s)
    if s=="P" then phi_k = max(0, -x+K(1,1+k))
    else  
        phi_k = max(0, x-K(1,1+k))
    end
endfunction



// Densité de la loi normale :
mu = 0;
sigma = 1;
function resultat=fct_normale(x)
    resultat = cdfnor("PQ", x, mu, sigma);
endfunction



// Variables pour la formule de Black-Scholes des Payoffs :
d1 = zeros(M+1, p);
d2 = zeros(M+1, p);

for i=1:p
    for j=1:M+1   
        d1(j, i) = (log(S0/K(1, j))+(r+0.5*sg(i)**2)*t)/(sg(i)*sqrt(t));
        d2(j, i) = d1(j, i) - sg(i)*sqrt(t);
    end
end


function valeur=dS0(Beta,i)
    if i>2*p then valeur = 0
    elseif i<=p then
        valeur = liste_alea(1,i)*exp(Beta(i+p)*liste_gauss(1,i))        
    else 
        valeur = S0 - liste_alea(1,i-p)*Beta(i)*(1-liste_gauss(1,i-p))*exp(Beta(i+p)*liste_gauss(1,i-p))
    end
endfunction

// Fonction d1, d2 :
function valeur=fct_d1(Beta)
    valeur = zeros(M+1,p)
    for i=1:p
        for j=1:M+1
            valeur(j,i) = (log(S0/K(1,j))+(r+0.5*Beta(1,i+p).^2)*t)/(Beta(1,i+p)*sqrt(t));
    end
    end
endfunction


function valeur=fct_d2(Beta)
     valeur = zeros(M+1,p)
    for i=1:p
        for j=1:M+1
            valeur(j,i) = (log(S0/K(1,j))+(r-0.5*Beta(1,i+p).^2)*t)/(Beta(1,i+p)*sqrt(t));
    end
    end
endfunction

//Fonction Dd1, Dd2 :
function valeur=Dd1(Beta,i)
    valeur = zeros(M+1,1)
    if i>2*p then valeur = zeros(M+1,1)
    elseif i>p then 
        for j=1:M+1
            valeur(j,1) = ((dS0(Beta,i)/S0+Beta(i)*t)*Beta(i)*sqrt(t)-sqrt(t)*(log(S0/K(1,j))+(r+0.5*Beta(i).^2)*t))/(Beta(i).^2*t)
            end
    else
        for j=1:M+1
            valeur(j,1) = (dS0(Beta,i)/S0)/(Beta(i+p)*sqrt(t))
            end
    end
endfunction

function valeur=Dd2(Beta,i)
    valeur = zeros(M+1,1)
    if i<=p then valeur = Dd1(Beta,i)
    elseif i<=2*p then
        valeur = Dd1(Beta,i)-sqrt(t)
    end
endfunction
// Fonction dfct_normale :
function valeur=df_normale(x)
    valeur = exp(-0.5*x.^2)/sqrt(2*%pi)
endfunction

// Prix d'un call : (alpha_k d'un Call)
Ci = zeros(p,M+1)
for i=1:M+1
    for j=1:p
        Ci(j, i) = p_n(j)*(S0*fct_normale(d1(i, j)) - K(1,i)*exp(-r*t)*fct_normale(d2(i, j)));
    end
end
C = zeros(1,M+1)
for i=1:M+1
    for j=1:p
        C(1,i) = C(1,i)+Ci(j,i) 
    end
end

clf();
subplot(121);
plot(K, C);
xtitle("Call Price versus Strike K")

// Prix d'un Put : (alpha_k d'un Put)
Pi = zeros(p,M+1)
for i=1:M+1
    for j=1:p
    Pi(j, i) = p_n(j)*(-S0*fct_normale(-d1(i, j)) + K(1,i)*exp(-r*t)*fct_normale(-d2(i, j)));
    end
end
P = zeros(1,M+1)
for i=1:M+1
    for j=1:p
        P(1,i) = P(1,i)+ Pi(j,i) 
    end
end


subplot(122);
plot(K, P);
xtitle("Put Price versus Strike K")

// Définition des fonctions de payoff à présent inconnues :

// Maitenant que l'on a les valeurs des alpha_k, on va essayer de retrouver nos
// conditions initiales lambda, sigma, p


// Optimisation en minimisant les résidus des moindres carrés :
// Fonction dS0 :

function val_payoff=fct_payoff(Beta)
    d1 = fct_d1(Beta)
    d2 = fct_d2(Beta)
    if type_operation == "C" then
            Ci_fct = zeros(p,M+1)
            for i=1:M+1
                for j=1:p
                Ci_fct(j, i) = Beta(2*p+j)*(S0*fct_normale(d1(i, j)) - K(1,i)*exp(-r*t)*fct_normale(d2(i, j)));
                end
            end
            val_payoff = zeros(1,M+1)
            for i=1:M+1
                for j=1:p
                    val_payoff(1,i) = val_payoff(1,i)+Ci_fct(j,i) 
                end
            end
        else
            Pi_fct = zeros(p,M+1)
            for i=1:M+1
                for j=1:p
                Pi_fct(j, i) = Beta(2*p+j)*(-S0*fct_normale(-d1(i, j)) + K(1,i)*exp(-r*t)*fct_normale(-d2(i, j)));
                end
            end
            val_payoff = zeros(1,M+1)
            for i=1:M+1
                for j=1:p
                    val_payoff(1,i) = val_payoff(1,i)+ Pi_fct(j,i) 
                end
            end
        end
endfunction

function valeur=residus(Beta) 
    //A Beta connu, renvoit Somme sur les strikes de (E(phi_k)-alpha_k)**2
    // Calcul de S0 :
    valeur = 100
    // Gestion des contraintes sur les probabilités p_i
    Beta(3*p) <> 1-sum(Beta(2*p+1:3*p-1))
    for i=1:p
        if Beta(2*p+i)<0 then valeur = 1000000000000000000000000000000000000
        elseif Beta(2*p+i)>1 then valeur = 100000000000000000000000000000000
        elseif Beta(i+p) == 0 then valeur = 1000000000000000000000000000000000000
        end
     end
     if valeur < 1000000000000 then
    // Calcul de d1 & d2 :
            
        if type_operation == "C" then
            C_fct = fct_payoff(Beta)
            valeur = sum((C_fct-C).^2);
        else
            P_fct = fct_payoff(Beta)
            valeur = sum((P_fct-P).^2);
        end
    end
endfunction

disp(residus(Beta_opt))


function gradient=grad_residus(Beta)
    liste_partielle = zeros(M+1,3*p)
    gradient = zeros(1,3*p)
    for i=1:p
        for j=1:M+1
             if type_operation == "C" then
                    d1 = fct_d1(Beta)
                    d2 = fct_d2(Beta)
                    der_d1 = Dd1(Beta,i)
                    der_d2 = Dd2(Beta,i)
                    der_d1p = Dd1(Beta,i+p)
                    der_d2p = Dd2(Beta,i+p)
                    pay_off = (fct_payoff(Beta)-C)
                    liste_partielle(j,i) = 2*Beta(2*p+i)*pay_off(1,j)*(dS0(Beta,i)*fct_normale(d1(j,i))+S0*df_normale(d1(j,i))*der_d1(j,1)-K(1,j)*exp(-r*t)*der_d2(j,1)*df_normale(d2(j,i)))
                    liste_partielle(j,i+p) = 2*Beta(2*p+i)*pay_off(1,j)*(dS0(Beta,i+p)*fct_normale(d1(j,i))+S0*df_normale(d1(j,i))*der_d1p(j,1)-K(1,j)*exp(-r*t)*der_d2p(j,1)*df_normale(d2(j,i)))
                    liste_partielle(j,i+2*p) = 2*pay_off(1,j)*(S0*fct_normale(d1(j,i))-K(1,j)*exp(-r*t)*fct_normale(d2(j,i)))
                end
               
    end
end
for i=1:p
    for j=1:M+1
        gradient(1,i) = gradient(1,i) + liste_partielle(j,i)
        gradient(1,i+p) = gradient(1,i+p) + liste_partielle(j,i+p)
        gradient(1,i+2*p) = gradient(1,i+2*p) + liste_partielle(j,i+2*p)
    end
end
endfunction

function [f,g,ind]=fct_objective(x, ind)
    f = residus(x);
    g = grad_residus(x)
endfunction

Beta0 = [70, 70, 0.01, 0.01, 0, 0]

Beta_inf = [60, 60, 0, 0, 0, 0]
Beta_sup = [100, 100, 1, 1, 1, 1]

[fopt,xopt] = optim(fct_objective, "b", Beta_inf, Beta_sup, Beta0)

disp("xopt = ")
disp(xopt)

disp(residus(xopt))
/////////////////////////////////////////////////////////////////////////////
// Mouvement Brownien
mu = 0;
T = 54;
vect_temps = zeros(1, T);
for i=1:T
	vect_temps(i)=i
end

sigma_b = t;
function resultat=fct_normale(x)
    resultat = cdfnor("PQ", x, mu, sigma_b);
endfunction

1
