clear all; close all;

% Parametos
% emisores
emi = 0:1:54264; % emisores distribuidos de forma continua
Ampemi = @(x) x .* (54264 - x) .* (sin(x/5000)+3); % amplitud de la emision
ampemi = Ampemi(emi);
Nemi = length(emi);

L = 763000; % distancia emisores y receptores
lambda = 0.2128;

% Receptores. Comentar y descomentar para el caso N = 8
Nrec = 64;
Drec = 100; % distancia maxima entre receptores
drec = Drec/(Nrec-1); % distancia minima entre receptores
rec = 0:drec:Drec;
%rec = (100/54)*[0, 19, 29, 33, 34, 36, 42, 54];
%Nrec = 8;
%Drec = 100;
%drec = 100/54;



% resolucion que tendremos
dx = L*lambda/Drec
% la cosa maxima que podremos ver bien
cosa_max = (Nrec-1)*dx


% distancias emisores-receptores
d = @(i,j) sqrt( L^2 + (emi(i)-rec(j))^2 );
dist = zeros(Nemi, Nrec);
for i=1:Nemi
    for j=1:Nrec
        dist(i,j) = d(i,j);
    end
end

% emision a cada t siguiendo una gaussiana
t_sample = 1250;
emis = zeros(Nemi, t_sample);
for i=1:Nemi
    emis(i,:) = normrnd(0, ampemi(i), [1,t_sample]);
end

% recepcion a cada t
k = 2*pi/lambda;
recep = zeros(Nrec, t_sample);
for j=1:Nrec
    recep(j,:) = ( exp(1i*k*dist(:,j)) ./ dist(:,j) )' * emis ;
end

% correlaciones
corr = zeros(Nrec,Nrec);
for j1 = 1:Nrec
    for j2 = 1:Nrec
        corr(j1,j2) = sum( recep(j1,:) .* conj(recep(j2,:)) )/t_sample;
        corr(j1,j2) = corr(j1,j2) * exp( -1i*k*(rec(j1)^2-rec(j2)^2)/(2*L) );
    end
end


% reconstruccion
nDiferencias = round(2*Drec/drec + 1); % numero de diferencias tanto en positivo como en negativo
nDif = round(Drec/drec +1);
mult = zeros(1, nDiferencias);
amplitud = zeros(1, nDiferencias);
for j1 = 1:Nrec
    for j2 = 1:Nrec
        frec = (rec(j2)-rec(j1)+Drec) / drec +1;
        frec = round(frec);
        mult(frec) = mult(frec)+1;
        amplitud(frec) = amplitud(frec) + corr(j1,j2);
    end
end
mult(find(mult == 0)) = 1;
amplitud = amplitud./mult;

%%
rec_xi = 0:1:54264;
rec_emis = zeros(1,length(rec_xi));
for i=1:length(rec_xi)
    for j=-(nDif-1):(nDif-1)
        rec_emis(i) = rec_emis(i) + amplitud(j+nDif)*exp( 1i*k/L*rec_xi(i)*j*drec);
    end
end
salida = L*sqrt(abs(rec_emis));

salida_aver = sqrt(mean(salida.^2));
entrada_aver = L * sqrt(mean(diag(corr))) / sqrt(length(rec_emis));
salida = (entrada_aver / salida_aver) * salida;

figure(6);
plot(rec_xi, salida, 'b');
hold on;
plot(rec_xi, ampemi, 'r');
grid on;
title('Noisy reconstruction (N = 64, 1250 samples)');
ylabel('Amplitude');
xlabel('Distance (m)');
legend('Reconstructed pattern', 'Emitted pattern');


nucleo = ceil(dx);
salida_correg = salida;
for i = 1:length(salida)
    ext_left = max([i-nucleo, 1]);
    ext_right = min([length(salida), i+nucleo]);
    salida_correg(i) = sum(salida(ext_left:ext_right)) / (2*nucleo+1);
end

resto = abs(salida_correg - ampemi);
figure(5);
plot(emi, ampemi, 'r');
hold on;
plot(rec_xi, salida_correg, 'b');
plot(rec_xi, resto, 'g');
title('Emission pattern');
ylabel('Amplitude');
xlabel('Distance (m)');
legend('Emitted pattern', 'Reconstructed pattern', 'Error');
grid on;
sqrt(mean(resto(200:5000).^2))
