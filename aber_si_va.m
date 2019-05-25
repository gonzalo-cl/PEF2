clear all; close all;

% Parametos
% emisores
emi = 0:1:54264; % muchos emisores pa que la cosa sea mas continua
Ampemi = @(x) x.*(sin(x/10)+4); % amplitud de la emision
ampemi = Ampemi(emi);
Nemi = length(emi);
L = 100000; % distancia emisores receptores
lambda = 0.2128;
% receptores: empezamos con una distribucion uniforme
Nrec = 256;
Drec = 100; % distancia maxima entre receptores
drec = Drec/(Nrec-1); % distancia minima entre receptores
rec = 0:drec:Drec;

% resolucion que tendremos
dx = L*lambda/Drec
% la cosa maxima que podremos ver bien
cosa_max = (Nrec-1)*dx
if (Nrec-1)*dx < emi(end)
    printf('cuidao');
    pause;
end

% distancias emisores-receptores
d = @(i,j) sqrt( L^2 + (emi(i)-rec(j))^2 );
dist = zeros(Nemi, Nrec);
for i=1:Nemi
    for j=1:Nrec
        dist(i,j) = d(i,j);
    end
end

% emision a cada t siguiendo una gaussiana
t_sample = 1000;
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
nDiferencias = 2*Drec/drec + 1; % numero de diferencias tanto en positivo como en negativo
nDif = Drec/drec +1;
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
amplitud = amplitud./mult;

%%
%%%%%%%%%%%%%%
figure(1), plot(emi, ampemi), title('amplitud teorica');
%figure(2), plot(
rec_xi = 0:1:54264;
rec_emis = zeros(1,length(rec_xi));
for i=1:length(rec_xi)
    for j=-(nDif-1):(nDif-1)
        rec_emis(i) = rec_emis(i) + amplitud(j+nDif)*exp( 1i*k/L*rec_xi(i)*j*drec);
    end
end
figure(3), plot(rec_xi, L*sqrt(abs(rec_emis)));

%% a partir de aqui es spam

coso = amplitud(1:nDif);
%figure(3), plot(abs(coso));
figure(4), plot(L*sqrt(abs(ifft(coso))));


ifft_ampl = ifft(amplitud);
%ifft1 = L*sqrt(abs(ifft(amplitud(1:nDif))));
ifft2 = L*sqrt(abs(ifft(amplitud(nDif:end)))); % crec que aquests son els bons




puntos_ifft = 0:dx:(Nrec-1)*dx;


figure(2), plot(puntos_ifft, ifft1)
figure(3), plot(puntos_ifft, ifft2);
    
