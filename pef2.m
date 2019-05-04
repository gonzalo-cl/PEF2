% Parametros del setting. 
%Posicion emisor y receptor. Ampltiud cada emisor
Nemi = 5;
xemi = [0, 0, 0, 0, 0];
yemi = [-2, -1, 0, 1, 2];
ampemi = [1, 5, 10, 10, 5];
Nrec = 9;
xrec = [30, 30, 30, 30, 30, 30, 30, 30, 30];
yrec = [-5, -4, -3, -1, 0, 1, 3, 4, 5];
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
    recep(j,:) = (1./d(:,j))' * emis;
end

    