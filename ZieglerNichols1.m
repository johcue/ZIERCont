%Ziegler - Nichols 1 
clc
num = [1 5];
den = [1 3 3 1]; % Nuestro ejercicio Mp 0% y ts < 5s
G=tf(num, den);
% Respuesta a una entrada escalón
syms s t;
[nums, dens] = tfdata(G);
Gsym = poly2sym(cell2mat(nums),s)/poly2sym(cell2mat(dens),s);
Gsym = simplify(Gsym);
fprintf('Planta : \n');
pretty(Gsym);
fprintf('Respuesta a una entrada escalón: C(t) = Transformada inversa de \n');
pretty(Gsym/s);
Ct = ilaplace(Gsym/s);
fprintf('La función de optimización es: C(t) = \n');
pretty(Ct);
fprintf('Hallando el punto de inflexión: \n\n');
fprintf('Condición de inflexión: diff2(C(t)) = 0, t = ?...\n')
fprintf('Derivada segunda de C(t): \n')
C1t = diff(Ct);
C2t = diff(C1t);
pretty(C2t)
fprintf('Igualando a cero, se obtiene que el posible (o los posibles valores) valor de t puede ser: \n')
ti = solve(C2t==0);
if (length(ti)>1)
    for i=1:length(ti)
        fprintf('   %d. ti = %.6f \n', i, ti(i))
    end
    for i=1:length(ti)
        if (ti(i)>0)
            ti = ti(i);
        end
    end
    fprintf('Eligiendo el valor positivo: ')
end
fprintf('ti = %.6f \n\n', ti)
m = double(subs(C1t,t,ti));
fprintf('Con este valor de ti, se obtiene que\n el valor de la pendiente de\n la recta de opt.,\n diff(C(t)) en t = ti,\n es: %.6f.\n\n', m);
b = double(subs(Ct,t,ti));
fprintf('Tambien, con este valor de ti, el\n valor del punto de corte en y,\n C(t) en t = ti, es: %.6f. \n', b);
gammat = vpa(m*t - m*ti + b);
fprintf('Hallando los valores de L y T, secantes a la recta de infl. : ');
T = solve(gammat == 1);
fprintf('T: t cuando y(t) = 1, T = %.6f\n\n', T);
L = double(solve(subs(gammat,t,'L')==0));
fprintf('L: t cuando y(t) = 0, L = %.6f\n\n', L);

%% Diseño del controlador PID
fprintf('--- Dis. del controlador PID ---\n\n');

gammat = m*t - m*ti + b;
fprintf('la recta de inflexión queda de la forma: %.6f*t%.6f\n\n', m, -m*ti+b);
Kp = 1.2*T/L;
Ti = 2*L;
Td = 0.5*L;
fprintf('Valor de Kp = 1.2*T/L = %.6f\n', Kp);
fprintf('Valor de Ti = 2*L = %.6f\n', Ti);
fprintf('Valor de Td = 0.5*L = %.6f\n', Td);

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
Kpm = 12.7*T/L;
Tim = 11.3*L;
Tdm = 10.5*L;
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

figure (1);
suptitle('Respuestas a entradas escalón - ZN1.');
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