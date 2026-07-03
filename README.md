close all
clear all
clc

% Trabalho - sistema Bola-Viga 
%% Equações do modelo não-linear do sistema

% Variáveis

% i: corrente; 
% xb: distancia entre o centro da bola e a extremidade da viga]
% alpha: angulo da viga com a horizontal
% theta: posição angular do motor
% uN: tensão na entrada do motor

mb  = 0.01;          % Massa da bola
b   = 2e-3;          % Coef. de atrito viscoso do ar com a bola
r   = 0.02;          % Raio da bola
ml  = 0.1;           % Massa da viga
l   = 1;             % Comprimento da viga
Jm  = 4e-4;          % Momento de inercia do motor
Jb  = (2/5)*mb*r^2;  % Momento de inercia da bola
Jl  = ml*l^2/3;      % Momneto de inercia da viga
Bm  = 1.1e-2;        % Coef. de atrito do motor
Rm  = 1;             % Resistencia da armadura
Km  = 0.03;          % Constante do motor
d   = 0.04;          % Distancia do eixo do motor e a haste
g   = 9.8;           % Aceleracao da gravidade


% Definindo as variáveis de estado e suas derivadas
syms xb xb_dot xb_ddot alpha alpha_ddot theta theta_dot theta_ddot f i Lm di dt uN


eq1_NL = (mb+Jb/r^2)*xb_ddot + b*xb_dot - mb*g*sin(alpha) == 0;                                         % Equação da bola
eq2_NL = Jl*alpha_ddot + ml*g*cos(alpha)*(l/2) + mb*g*cos(alpha)*(l-xb) + f*l*cos(alpha) == 0;          % Equação da viga
eq3_NL = l*sin(alpha) - d*sin(theta) == 0;                                                              % Equação da relação entre theta e alpha
eq4_NL = Jm*theta_ddot + Bm*theta_dot - f*cos(theta)*d - Km*i == 0;                                     % Equação do motor de corrente contínua
eq5_NL = uN - Rm*i - Lm*(di/dt) - Km*theta_dot == 0;                                                    % Equação do motor de corrente contínua


%% 1) Obtenha o modelo linearizado para pequenas perturbações, 
% considerando as aproximações sen(alpha)~alpha, cos(alpha)~1 e 
% sen(theta)~theta, cos(theta)~1. Além disso, considere desprezível a indutância de 
% armadura, isto é, Lm ~ 0. 
% Nessas condições, reescreva o modelo matemático do sistema (1)-(5)
% somente em função das variáveis xb, theta e uN. 


% Equações Linearizadas para pequenos deslocamentos:
eq1_L = (mb+Jb/r^2)*xb_ddot + b*xb_dot - mb*g*alpha == 0;                            % Equação da bola
eq2_L = Jl*alpha_ddot + ml*g*(l/2) + mb*g*1*(l-xb) + f*l*1 == 0;                                  % Equação da viga
eq3_L = l*alpha - d*theta == 0;                                                      % Equação da relação entre theta e alpha
eq4_L = Jm*theta_ddot + Bm*theta_dot - f*1*d - Km*i == 0;                            % Equação do motor de corrente contínua
eq5_L = uN - Rm*i - Km*theta_dot == 0;                                               % Equação do motor de corrente contínua

% Reescrevendo as equações em função de xb, theta e uN:
% -> Da eq3 temos que alpha = d*theta/l. 
% Substituindo alpha nas equações
eq1_L = subs(eq1_L, alpha, d*theta/l);
eq2_L = subs(eq2_L, alpha_ddot, d*theta_ddot/l);

% -> Da eq2 temos que: f = (ml*g*1*(xb-l) - Jl*alpha_ddot)/l;
% -> Da eq5 temos que: i = (uN - Km*theta_dot)/Rm
eq4_L = subs(eq4_L, f, (mb*g*1*(xb-l) - ml*g*(l/2) - Jl*alpha_ddot)/l);
eq4_L = subs(eq4_L, [alpha_ddot, i], [d*theta_ddot/l, (uN - Km*theta_dot)/Rm]);


%% 2) Utilizando o modelo obtido no item anterior, determine a tensão uo 
% para que o sistema permaneça em equilíbrio no ponto de operação 
% (xb_o; theta_o) = (0,5, 0°)

% Definindo os pontos de operação
xb_o = 0.5;  % Posição de equilíbrio da bola
theta_o = 0; % Posição de equilíbrio do motor


% Calculando a tensão uo para o ponto de operação. No ponto de equilíbrio
% TODAS as DERIVADAS são nulas, isto é, iguais a zero, logo:

alpha_o = d*theta_o/l;                                  % Posição de equilíbrio da viga
f_o = (mb*g*(xb_o - l) - ml*g*(l/2))/l;               % Força no ponto de operação
i_o = -(f_o*1*d)/Km;                                    % Corrente no ponto de operação
uo = Rm*i_o;                                            % Tensão no ponto de operação

%% 3) Determine a representação em espaço de estado do modelo linearizado
% x_dot = Ax + Bu, x(0) = 0; 
% y = Cx + Du
% considerando o seguinte vetor de estado x1 = xb-xbo, x2 = xb_dot, 
% x3 = theta, x4 = theta_dot, a entrada de controle u = uN-uo e a saída 
% dada por y = x1.

% eq1_nova = (mb+Jb/r^2)*xb_ddot + b*xb_dot - mb*g*d*theta/l == 0;                                                                                      % Equação da bola
% eq2_nova = Jl*d*theta_ddot/l - Jl*(d*theta_ddot/l)*l*1 == 0;                                                                                          % Equação da viga                                                                                     % Equação da relação entre theta e alpha
% eq3_nova = Jm*theta_ddot + Bm*theta_dot - ((ml*g*1*(xb-l) - Jl*(d*theta_ddot/l))/l)*1*d - Km*(uN - Km*theta_dot)/Rm == 0;                             % Equação do motor de corrente contínua 

% Variáveis auxiliares para simplificar as matrizes:
M_eq = mb + Jb/r^2;            % Inércia equivalente da bola
J_eq = Jm + (Jl * d^2) / l^2;  % Inércia equivalente do motor/viga
B_eq = (Bm + Km^2/Rm); 

A = [0, 1, 0, 0;
    0, -b/M_eq, (mb*g*d)/(l*M_eq), 0;
    0, 0, 0, 1;
    (mb*g*d)/(l*J_eq), 0, 0, -B_eq/J_eq];

B = [0;
    0;
    0;
    Km/(Rm*J_eq)];

C = [1, 0, 0, 0];
D = 0;

%% 4) Representação em Espaço de Estado para o sistema não-linear utilizando 
% as seguintes variaveis de estado z1 = xb, z2 = xb_dot, z3 = theta, 
% z4 = theta_dot a entrada uN e a saída incremental dada por y = z1 - xbo

k = d/l;
delta  = 1 - k^2*(sin(theta))^2;
cos(alpha) = sqrt(delta);
alpha_ddot = (k*sin(theta)*theta_ddot^2*(k^2-1) + k*cos(theta)*theta_ddot*delta)/(delta^(3/2));
sin(alpha) = k*sin(theta);

% eq1_NL = (mb+Jb/r^2)*xb_ddot + b*xb_dot - mb*g*sin(alpha) == 0;                                       % Equação da bola
% eq2_NL = Jl*alpha_ddot + ml*g*cos(alpha)*(l/2) + mb*g*cos(alpha)*(l-xb) + f*l*cos(alpha) == 0;        % Equação da viga
% eq3_NL = l*sin(alpha) - d*sin(theta) == 0;                                                            % Equação da relação entre theta e alpha
% eq4_NL = Jm*theta_ddot + Bm*theta_dot - f*cos(theta)*d - Km*i == 0;                                   % Equação do motor de corrente contínua
% eq5_NL = uN - Rm*i - Lm*(di/dt) - Km*theta_dot == 0;                                                  % Equação do motor de corrente contínua

Lm = 0;                                                                          % Aproximação dada pelo enunciado
i  = (uN - Km*theta_dot)/Rm;
f  = (1/d)*(Jm*theta_ddot + B*theta_dot - Km*(uN - Km*theta_dot)/Rm)/cos(theta); 


eq1_NL = subs(eq1_NL, sin(alpha), k*sin(theta));
eq2_NL = subs(eq2_NL, f,(1/d)*(Jm*theta_ddot + B*theta_dot - Km*(uN - Km*theta_dot)/Rm)/cos(theta)); 
eq4_NL = subs(eq1_NL, i, (uN - Km*theta_dot)/Rm);

% Isolando xb_ddot na eq1_NL
xb_ddot = (mb*g*k*sin(theta) - b*xb_dot)/(mb + Jb/r^2);

% Isolando theta_ddot na eq2_NL
% theta_ddot = Jl*k*theta_dot^2*(k^2-1)*sin(theta) + (delta^2)*(l/2*ml*g + (l - xb)*mb*g + theta_dot((B+(Km^2/Rm))/k*cos(theta)) - (Km*uN)/(k*Rm*cos(theta)) / (Jl*k*delta*cos(theta) + (delta^2*Jm/(k*cos(theta)))); % para simplificar multiplica-se por kcos(theta)/kcos(theta)
theta_ddot = Jl*k^2*theta_dot^2*(k^2-1)*sin(theta)*cos(theta) + (delta^2)*(k*cos(theta)*g*(l/2*ml + (l - xb)*mb) + theta_dot*(B+(Km^2/Rm)) - (Km*uN)/Rm) / (Jl*k^2*delta*cos(theta)^2 + (delta^2*Jm));

% Representação em Espaços de Estados
syms z1 z2 z3 z4 z1_dot z2_dot z3_dot z4_dot

z1 = xb;
z2 = xb_dot;
z3 = theta;
z4 = theta_dot;

z1_dot  = xb_dot;
z2_dot  = (mb*g*k*sin(z3) - b*z2)/(mb + Jb/r^2); % = xb_ddot 
z3_dot  = z4; % = theta_dot  
z4_dot  = theta_ddot;

Z_dot = [z1_dot
         z2_dot
         z3_dot
         z4_dot];

%% 5) Função de Transferência G(s) = X1(s) / U(s) e Estabilidade

syms x1 x3 x1_dot x3_dot x1_ddot x3_ddot u X1 X3 s U G


eq1_sistema = M_eq*x1_ddot + b*x1_dot == mb*g*k*x3;
eq2_sistema = J_eq*x3_ddot + B_eq*x3_dot - (k*mb*g)*x1 == (Km/Rm)*u;

% Aplicando a Transformada de Laplace (Condições Iniciais Nulas)
eq1_s = subs(eq1_sistema, [x1_ddot, x1_dot, x3], [s^2*X1, s*X1, X3]);
eq2_s = subs(eq2_sistema, [x3_ddot, x3_dot, x1, u], [s^2*X3, s*X3, X1, U]);

X3_expr = solve(eq1_s, X3);
eq3_s = subs(eq2_s, X3, X3_expr);
X1_expr = solve(eq3_s, X1);
G = X1_expr / U;

[num, den] = numden(G);          % Extrai numerador e denominador

% Encontra o coeficiente do termo de maior grau (s^4)
coef_lider = coeffs(den, s, 'All');
fator_divisao = coef_lider(1); 

% Normaliza e agrupa o numerador e o denominador nas potências de 's'
num_norm = collect(num / fator_divisao, s);
den_norm = collect(den / fator_divisao, s);

% Monta a função de transferência e aplica o VPA na expressão COMPLETA
G_normalizado = vpa(num_norm / den_norm, 4)

%% Análise de Estabilidade 

polos = double(solve(den_norm == 0, s));

disp('Os polos da Função de Transferência são:');
disp(polos);

%------------------------------------------------
% 6) Matrizes de Controlabilidade e Observabilidade
%------------------------------------------------

% Cria o objeto do sistema LTI no MATLAB
sys = ss(A, B, C, D);

%% 4. Análise de Controlabilidade (Questão 6)
% Método 1: Usando a função nativa do MATLAB
Mc = ctrb(A, B);

% Cálculo do Posto e Determinante
posto_Mc = rank(Mc);
det_Mc = det(Mc);

%% 5. Análise de Observabilidade (Questão 6)
% Método 1: Usando a função nativa do MATLAB
Mo = obsv(A, C);

% Cálculo do Posto e Determinante
posto_Mo = rank(Mo);
det_Mo = det(Mo);

%% 6. Exibição dos Resultados na Tela

fprintf('\n----------------------------------------------------\n\n');

disp('Matriz de Controlabilidade (Mc):');
disp(Mc);
fprintf('Posto de Mc = %d (Esperado: 4)\n', posto_Mc);
fprintf('Determinante de Mc = %e\n\n', det_Mc);

if posto_Mc == size(A,1)
    disp('-> CONCLUSÃO: O sistema é COMPLETAMENTE CONTROLÁVEL.');
else
    disp('-> CONCLUSÃO: O sistema NÃO é totalmente controlável.');
end

fprintf('\n----------------------------------------------------\n\n');

disp('Matriz de Observabilidade (Mo):');
disp(Mo);
fprintf('Posto de Mo = %d (Esperado: 4)\n', posto_Mo);
fprintf('Determinante de Mo = %e\n\n', det_Mo);

if posto_Mo == size(A,1)
    disp('-> CONCLUSÃO: O sistema é COMPLETAMENTE OBSERVÁVEL.');
else
    disp('-> CONCLUSÃO: O sistema NÃO é totalmente observável.');
end

%% PARTE 2 - CONTROLE CLÁSSICO

% Lugar das Raízes de G(s)

% Extrai o numerador e denominador simbólicos de G_normalizado
[num_sym, den_sym] = numden(G_normalizado);

% Converter os polinômios simbólicos em vetores de números
num_num = sym2poly(num_sym);
den_num = sym2poly(den_sym);

LR_G = tf(num_num, den_num);

figure;
rlocus(LR_G);
grid on;
title('Lugar das Raízes - Sistema Bola-Viga');

% Alocação de Polos (Equação Diofantina)

% Polinômio P(s) escolhido para 1 + CG
P = (s+1)^2*(s+2)*(s+4)*(s+5)*(s+6)*(s+8)*(s+10)^2;

% Expande o polinômio P e extrai os coeficientes
P = expand(P);
pv = coeffs(P, s, 'All')';

% Expande o polinômio do denominador da planta D_G e obtém-se os
% coeficientes
D_G = expand(den_sym*s); % Multiplica-se por s para garantir um sistema do tipo 1
dgv = coeffs(D_G, s, 'All');

% Expande o polinômio do numerador da planta N_G e obtém-se os coeficientes
N_G = expand(num_sym);
ngv = coeffs(N_G, s, 'All');

% Definindo a matriz de Silvester para os polinômios numerador e
% denominador
n_C = 4; % Ordem do denominador do controlador
m_C = 4; % Ordem do numerador do controlador

S = matriz_sylvester_controle(dgv, ngv, n_C, m_C)

% Cálculo do vetor de coeficientes do controlador com a equação Diofantina
cv = (S'*S)\(S'*pv(end:-1:1));
%cv = S\pv;

% Conversão dos coeficientes no controlador C_0
C_0 = poly2sym(cv(end:-1:n_C + 2), s) / poly2sym(cv(n_C + 1:-1:1), s);

% Cálculo do controlador C com base em C_0 / s
C = vpa(C_0 / s, 4)




