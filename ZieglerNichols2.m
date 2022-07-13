%Ziegler - Nichols 2 
clc
num = 100;
den = [1 7 12 0]; % Nuestro ejercicio Mp 45% y ts < 1.5s
G=tf(num, den);
% Respuesta a una entrada escalón
syms s t w K;
[nums, dens] = tfdata(G);
Gsym = poly2sym(cell2mat(nums),s)/simplify(poly2sym(cell2mat(dens),s));
Gsym = simplify(Gsym);
fprintf('Planta : \n');
pretty(Gsym);
Ga = feedback(G,1);
[nums, dens] = tfdata(Ga);
numK = coeffs(poly2sym(cell2mat(nums),s));
numK(length(numK)) = numK(length(numK))*K;
denK = coeffs(poly2sym(cell2mat(dens),s));
denK = fliplr(denK);
denK(length(denK)) = denK(length(denK))*K;
Gasym = poly2sym(numK,s)/poly2sym(denK,s);
Gasym = expand(Gasym);
fprintf('Planta para criterio de Routh - Horwitz: \n');
pretty(Gasym);
M = sym(zeros(length(denK),length(denK)-1));

if (mod(length(denK),2)==0)
    limC = length(denK)/2 + 1;
else
    limC = ceil(length(denK)/2);
end

for i=1:length(denK)
    M(i,1) = s^(length(denK)-i);
end
cont = 1;
for i=2:limC
    M(1,i) = denK(cont);
    M(2,i) = denK(cont+1);
    cont=cont+2;
end

for i=3:length(denK)
    cont = 2;
    while(cont < limC)
        M(i,cont) = (M(i-1,2)*M(i-2, cont+1)-M(i-2,2)*M(i-1,cont+1))/M(i-1,2);
        cont = cont+1;
    end 
end
filling = 1;
Kcr = zeros(length(denK));
for i=1:length(denK)
    for j=2:limC
        if (symvar(M(i,j))==K)
            if (solve(M(i,j)==0, K) ~= 0)
                Kcr(filling) = solve(M(i,j)==0, K);
                filling = filling + 1;
            end
        end
    end
end

fprintf('Resultados del criterio de Routh - Horwitz:\n');
disp(M);
if (length(Kcr)>1)
     fprintf('Los valores de Kcr (=/= de 0) que generan polos complejos dominantes, son:\n');
    for i=1:length(Kcr)
        fprintf('   %d. Kcr = %.6f \n', i, Kcr(i))
    end
    for i=1:length(Kcr)
        if (Kcr(i)>0)
            Kcr = Kcr(i);
            break
        end
    end
    fprintf('Eligiendo el valor positivo: Kcr = ')
else
    fprintf('Solo se encontró un valor positivo de Kcr que genera polos complejos dominantes, Kcr = ');
end
fprintf('%.6f.\n\n', Kcr);
fprintf('Dada la expresion: \n');
denK = subs(denK,K,Kcr);
pol = poly2sym(denK,s);
pretty(pol)
fprintf('Haciendo la conversion s = jw, se obtiene: \n');
poljw = subs(pol,s,1j*w);
pretty(expand(poljw))
fprintf('Hallando sus raices. \n');
fprintf('Parte real: 0 = \n');
repoljw = subs(poljw,1i,0);
pretty(repoljw);
fprintf('Parte imaginaria: 0 = \n');
impoljw = poljw - subs(poljw,1i,0);
pretty(impoljw)
fprintf('\n');
Wcrre = solve(repoljw==0);
Wcrim = solve(impoljw==0);
Wcr = [Wcrre; Wcrim];
if (length(Wcr)>1)
     fprintf('Los valores de Wcr (=/= de 0) raíces del polinomio real, son:\n');
    for i=1:length(Wcr)
        fprintf('   %d. Wcr = %.6f \n', i, Wcr(i))
    end
    for i=1:length(Wcr)
        if (Wcr(i)>0)
            Wcr = Wcr(i);
            break
        end
    end
    fprintf('\nEligiendo el valor positivo: Wcr = ')
else
    fprintf('\nSolo se encontró un valor positivo de Wcr raíz, Wcr = ');
end
fprintf('%.6f.\n\n', Wcr);
Pcr = 2*pi/Wcr;
fprintf('P crítico calculado Pcr = %.6f.\n\n', Pcr);

%% Diseño del controlador PID
fprintf('--- Dis. del controlador PID ---\n\n');

Kp = 0.6*Kcr;
Ti = 0.5*Pcr;
Td = 0.125*Pcr;
fprintf('Valor de Kp = 0.6*Kcr = %.6f\n', Kp);
fprintf('Valor de Ti = 0.5*Pcr = %.6f\n', Ti);
fprintf('Valor de Td = 0.125*Pcr = %.6f\n', Td);

Kd = Kp*Td;
Ki = Kp/Ti;
fprintf('Valor de Kd = Kp*Td = %.6f\n', Kd);
fprintf('Valor de Ki = Kp/Ti = %.6f\n', Ki);

numc = double([Kd Kp Ki]);
denc = [1 0];
GcPID = tf(numc,denc);
[nums, dens] = tfdata(GcPID);
Gcsym = poly2sym(cell2mat(nums),s)/poly2sym(cell2mat(dens),s);
fprintf('\nFuncion del controlador obtenida :\n');
pretty(vpa(Gcsym));
Glc=feedback(GcPID*G,1);
[nums, dens] = tfdata(Glc);
Glcsym = poly2sym(cell2mat(nums),s)/poly2sym(cell2mat(dens),s);
fprintf('\nLa funcion en lazo cerrado corresponde a: Glc(Gc(s)*G(s)) = \n');
pretty(vpa(Glcsym));
% Obtención de los coeficientes del denominador para factorizar 
denv=double(coeffs(poly2sym(cell2mat(dens),s)));
raicesden = roots(denv);
denfac=1;
cont=1;
while cont-1<numel(raicesden)
    if (imag(raicesden(cont))==0)
        denfac=denfac*(s-(raicesden(cont)));
    else
        denfac=denfac*(s^2-(raicesden(cont)+raicesden(cont+1))*s+raicesden(cont)*raicesden(cont+1));
        cont=cont+1;
    end
    cont=cont+1;
end
% Con esto, se arregla el denominador y queda de forma factorizada.
Glcsym = poly2sym(cell2mat(nums),s)/denfac; %=X(s)Y(s)
fprintf('\nFuncion en lazo cerrado simplificada: Glc = \n');
pretty(vpa(Glcsym));

%% Sintonizacion de parametros para obtener los resultados deseados
Kpm = 0.6*Kcr;
Tim = 0.5*Pcr;
Tdm = 0.125*Pcr;
fprintf('Valor de Kp modificado = %.6f\n', Kpm);
fprintf('Valor de Ti modificado = %.6f\n', Tim);
fprintf('Valor de Td modificado = %.6f\n', Tdm);
%% Termina modificación de parámetros

Kd = Kpm*Tdm;
Ki = Kpm/Tim;
fprintf('Valor de Kd modificado = Kp*Td = %.6f\n', Kd);
fprintf('Valor de Ki modificado = Kp/Ti = %.6f\n', Ki);
numc = double([Kd Kpm Ki]);
denc = [1 0];
GcmPID = tf(numc,denc); %Controlador "sintonizado"
[nums, dens] = tfdata(GcmPID);
Gcmsym = poly2sym(cell2mat(nums),s)/poly2sym(cell2mat(dens),s);
fprintf('\nFuncion del controlador modificado obtenida :\n');
pretty(vpa(Gcsym));
Glcm = feedback(GcmPID*G,1);
[nums, dens] = tfdata(Glcm);
Glcmsym = poly2sym(cell2mat(nums),s)/poly2sym(cell2mat(dens),s);
fprintf('\nLa funcion en lazo cerrado modificada corresponde a: Glc(Gcm(s)*G(s)) = \n');
pretty(vpa(Glcmsym));
% Obtención de los coeficientes del denominador para factorizar 
denv=double(coeffs(poly2sym(cell2mat(dens),s)));
raicesden = roots(denv);
denfac=1;
cont=1;
while cont-1<numel(raicesden)
    if (imag(raicesden(cont))==0)
        denfac=denfac*(s-(raicesden(cont)));
    else
        denfac=denfac*(s^2-(raicesden(cont)+raicesden(cont+1))*s+raicesden(cont)*raicesden(cont+1));
        cont=cont+1;
    end
    cont=cont+1;
end
% Con esto, se arregla el denominador y queda de forma factorizada.
Glcmsym = poly2sym(cell2mat(nums),s)/denfac; %=X(s)Y(s)
fprintf('\nFuncion en lazo cerrado modificada simplificada: Glc = \n');
pretty(vpa(Glcmsym));

% Grafica de la planta en lazo cerrado step(feedback(G,1))
% Grafica del primer controlador en lazo cerrado step(Glc)
% Grafica del controlador sintonizado step(Glcm)
fprintf('\nLa respuesta de este sistema a una entrada escalón, es: y(t) = \n');
yt = ilaplace(Glcsym);
pretty(vpa(yt));

opts = timeoptions('cstprefs');
opts.TimeUnits = 'auto';
opts.Grid='on';
opts.YLabel.String = 'Amplitud';
opts.XLabel.String = 'Tiempo';
opts.Title.FontSize = 12;

figure (2);
suptitle('Respuestas a entradas escalón - ZN2.');
subplot(3,1,1);
opts.Title.String = 'Planta en lazo cerrado';
step(feedback(G,1),'y-',opts);
info = stepinfo(feedback(G,1));
texto = strcat('\uparrow Tiempo de asentamiento t_{s}:', {'  '}, num2str(info.SettlingTime), 's.');
text(info.SettlingTime, 0.5, texto);
texto = strcat('  \leftarrow Máximo sobreimpulso M_{p}:', {'  '}, num2str(info.Overshoot),'%.');
text(info.PeakTime, info.Peak, texto);
subplot(3,1,2);
opts.Title.String = 'Primer controlador';
step(Glc,'c-',opts);
info = stepinfo(Glc);
texto = strcat('\uparrow Tiempo de asentamiento t_{s}:', {'  '}, num2str(info.SettlingTime), 's.');
text(info.SettlingTime, 0.5, texto);
texto = strcat('  \leftarrow Máximo sobreimpulso M_{p}:', {'  '}, num2str(info.Overshoot),'%.');
text(info.PeakTime, info.Peak, texto);
subplot(3,1,3);
opts.Title.String = 'Controlador Sintonizado';
step(Glcm,'g-',opts);
info = stepinfo(Glcm);
texto = strcat('\uparrow Tiempo de asentamiento t_{s}:', {'  '}, num2str(info.SettlingTime), 's.');
text(info.SettlingTime, 0.5, texto);
texto = strcat('\uparrow Máximo sobreimpulso M_{p}:', {'  '}, num2str(info.Overshoot),'%.');
text(info.PeakTime, info.Peak-0.2, texto);