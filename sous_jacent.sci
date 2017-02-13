

function [M] = marcheAleatoire(N)
 Y = 2*grand(N,1,'bin',1,1/2) - 1;
 M = [0;cumsum(Y)];
endfunction

N = 100000;
abscisse = 0 :1 :N;
M = marcheAleatoire(N);

function [B] = brownien(N,couleur)
   B = M(1 :N+1)/sqrt(N) ;
   abscisse = 0 :1/N :1 ;
   plot2d(abscisse,B,style=couleur)
endfunction


function S = BS(mu, sigma)
 S = exp(mu*t-((sigma^2)/2)*t+sigma*B);
endfunction


t=(0:1/N:1)';
B=marcheAleatoire(N)/sqrt(N);
// xset('window',0)
clf();

liste_mu = [0, 0.02, 0.05, 0.1, 0.15, 0.2]
liste_sigma = [0.01, 0.05, 0.1, 0.15, 0.2,0.4]
subplot(211);

for i=1:6 do
   plot2d(t, BS(liste_mu(i),liste_sigma(6)),i+1);
end
legends(['mu=0%';'mu=2%';'mu=5%';'mu=10%';'mu=15%';'mu=20%'],[2,3,4,5,6,7],'ur');

xtitle('Valeur du sous-jacent dans le temps pour différentes dérives');

subplot(212);

for i=1:6 do
   plot2d(t, BS(liste_mu(4),liste_sigma(i)),i+1);
end

legends(['sigma=1%';'sigma=5%';'sigma=10%';'sigma=15%';'sigma=20%';'sigma=30%'],[2,3,4,5,6,7],'ur');

xtitle('Valeur du sous-jacent dans le temps pour différentes volatilités');
