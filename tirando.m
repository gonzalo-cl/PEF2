%nota para el medonasso. Me ha sorpresido ver que la reconstruccion de la
%imagen sale en el plot 4 (wtf) que es el ifft de cross_aux. Pero mira
%funsiona asi que no me voy a quejar. Lo unico que sale reescalado en
%anchura, pero eso es tipico al estar haciendo transformaciones de fourier
%discreta (estamos tomando como uno lo que realmente es la distancia de
%separacion entre emisores o receptores)


clear all; close all;
% Parametros del setting. 
%Posicion emisor y receptor. Ampltiud cada emisor
Nemi = 5000;
Demi = 0.001; %%Distancia entre emisores. Muy pequenia
xemi = Demi*(1:Nemi);
%yemi = linspace(-1000, 1000, Nemi);
yemi = Demi*(1:Nemi);
%ampemi = [100 * ones(1, Nemi/8), 10000*ones(1,Nemi/4), 1000*ones(1,Nemi/8), 100000 * ones(1, Nemi/2)];
%ampemi = (1+sin(yemi)).*yemi*10000;
ampemi = 1+4*exp(-(yemi/Demi-70).^2/4096) + exp(-(yemi/Demi-200).^2/196);
%ampemi = [10000*ones(1,Nemi/4), 100000 * ones(1, Nemi*3/4)];
%ampemi = 1+sin(yemi/50);
figure(10),plot(yemi, ampemi),title('amplitud teorica');
Nrec = 256;
L = 25600;
lambda = 0.00001;
xrec = L*ones(1, Nrec);
yrec = (1:Nrec);

for i = 1:Nemi
    for j = 1:Nrec
        d(i,j) = sqrt((xemi(i) - xrec(j)).^2 + (yemi(i) - yrec(j)).^2);
        % observa que, para L>>distanciaentreemisores y L>>distanciaentrereceptores
        % se satisface d(i,j) ~ L + (yemi(i)-yrec(j))^2/(2*L);
    end
end

% Emision a cada t siguiendo una gaussiana
t_sample = 5000;
emis = zeros(Nemi, t_sample);
for i = 1:Nemi
    emis(i, :) = normrnd(0, ampemi(i), [1, t_sample]);
end

% Recepcion a cada t sumando contribuciones
k = 2*pi/lambda;
recep = zeros(Nrec, t_sample);
for j = 1:Nrec
  %  recep(j,:) = real(( exp(k*sqrt(-1)*d(:,j)) ./ d(:,j) )' * emis );
    recep(j,:) = ( exp(k*sqrt(-1)*d(:,j)) ./ d(:,j) )' * emis ;
end

% Calculamos las correlaciones
% aqui, en corr2, he dividido las correlaciones por un termino corrector,
% para que solo dependan de la diferencia entre antenas. no se, si se hacen
% los numeritos sale, se me hace dificil explicarlo bien por aqui
% (o sea, usando la aproximacion que hay en el calculo de las distancias
% sale bien)
corr = zeros(Nrec,Nrec);
corr2 = corr;
for j1 = 1:Nrec
    for j2 = 1:Nrec
        corr(j1,j2) = sum( recep(j1,:) .* conj(recep(j2,:)) )/t_sample;
        corr2(j1,j2) = corr(j1,j2) * exp( -sqrt(-1)*k*(yrec(j1)^2-yrec(j2)^2)/2/xrec(1) );
    end
end

%Brute-force-approach para reconstrur. Para cada posible diferencia
%receptor i - receptor j nos la guardamos. mult tiene cuantas veces aparece
%cada diferencia (entre -10 y 10). Para poder meterla como indice le
%sumamos 11. Es decir mult(5) = 3 significa que hay tres pares de
%receptores con diferencia en posicion 5-11 = -6. Luego, para cada
%diferencia promediamos ampltidues (al no haber error y suponiendo
%rayos paralelos, seran todas parecidas). Finalmente reconstruimos Fourier

% Receta de encode: Cogete frec(j1,j2) y sumale el menor valor posible para
% que el nuevo minimo sea 0. Ahora multiplica por la "resolucion
% frecuencial" es decir, haz que el segundo valor sea 1. Por ultimo suma 1
res_frec = (max(yrec)-min(yrec))/(Nrec-1);

mult = zeros(1, 2*Nrec - 1);
amplitud = zeros(1, 2*Nrec - 1);
amplitud2 = amplitud; % calculada con corr2, en vez de con corr
for j1 = 1:Nrec
    for j2 = 1:Nrec
        frec = (yrec(j2) - yrec(j1) -min(yrec)+max(yrec)) / res_frec + 1;
        mult(frec) = mult(frec) + 1;
        amplitud(frec) = amplitud(frec) + corr(j1,j2);
        amplitud2(frec) = amplitud2(frec) + corr2(j1,j2);
    end
end
hold on;
amplitud = amplitud ./ mult;
amplitud2 = amplitud2 ./ mult;

figure(20), plot(abs(amplitud2)) % fig20 & fig21 no son nada bueno
figure(21), plot(abs(ifft(amplitud2)))

amplitud2 = amplitud2(1:Nrec); % aqui cojo solo esta mitad porque asi iba bien

% aqui me habia calculado una cosa auxiliar para ir viendo si las cosas
% iban bien. serian las correlaciones (amplitud2) asumiendo que el caracter
% estocastico de los ampemi salen perfectos
corr_aux = zeros(1,Nrec);
for t=1:Nrec
    for x=1:Nemi
        corr_aux(t) = corr_aux(t) + abs(ampemi(x))^2/L^2*exp(-sqrt(-1)*k*yemi(x)*t/L);
    end
end

% Drec = yrec(2)-yrec(1);
% res = lambda*L/Drec/Nrec
% visib = lambda*L/Drec


figure(1), plot(abs(amplitud2));
figure(2), plot(abs(corr_aux));
figure(3), plot(L*sqrt(abs(ifft(amplitud2))));
figure(4), plot(L*sqrt(abs(ifft(corr_aux))));

% % % 
% % % % COSA QUE DEBERIAMOS ARREGLAR:
% % % 
% % % dadas las correlaciones entre distancias 0, Dy, 2Dy, ..., (N-1)Dy
% % % la ifft nos devuelve como son las emisiones como si tuvieramos N emisores en
% % % 0, x, 2x, ..., (N-1)x, con x = lambda*L / Dy*N
% % %     
% % % entonces si situamos Nrec emisores (o menos) y estan en las buenas posiciones
% % % (ie, Dx*Dy*N/lambda*L = 1) sale bastante bien. sino no.