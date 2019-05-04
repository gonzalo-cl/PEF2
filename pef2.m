% Parametros del setting. 
%Posicion emisor y receptor. Ampltiud cada emisor
Nemi = 100;
xemi = zeros(1,Nemi);
yemi = linspace(-1, 1, Nemi);
ampemi = [10000 * ones(1, Nemi/2), 1000000 * ones(1, Nemi/2)];
Nrec = 161;
xrec = 100*ones(1, Nrec);
yrec = linspace(-20, 20, 161);
for i = 1:Nemi
    for j = 1:Nrec
        d(i,j) = sqrt((xemi(i) - xrec(j)).^2 + (yemi(i) - yrec(j)).^2);
    end
end

% Emision a cada t siguiendo una gaussiana
t_sample = 50;
emis = zeros(Nemi, t_sample);
for i = 1:Nemi
    emis(i, :) = normrnd(0, ampemi(i), [1, t_sample]);
end

% Recepcion a cada t sumando contribuciones
recep = zeros(Nrec, t_sample);
for j = 1:Nrec
    recep(j,:) = ( exp(2*pi*0.5*sqrt(-1)*d(:,j)) ./ d(:,j) )' * emis;
end

% Calculamos las correlaciones
for j1 = 1:Nrec
    for j2 = 1:Nrec
        corr(j1,j2) = sum( recep(j1,:) .* conj(recep(j2,:)) );
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

mult = zeros(1, 2*Nrec - 1);
amplitud = zeros(1, 2*Nrec - 1);
for j1 = 1:Nrec
    for j2 = 1:Nrec
        frec(j1,j2) = (yrec(j2) - yrec(j1) + 40) * 4 + 1;
        mult(frec(j1,j2))  = mult(frec(j1,j2)) + 1;
        amplitud(frec(j1,j2)) = amplitud(frec(j1,j2)) + corr(j1,j2);
    end
end
hold on;
amplitud = amplitud ./ mult;
plot(abs(amplitud));
hold on;
figure(2);
plot(abs(ifft(amplitud)));
        
    