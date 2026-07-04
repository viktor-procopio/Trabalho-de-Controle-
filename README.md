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
syms xb xb_dot xb_ddot alpha alpha_ddot theta theta_dot theta_ddot f i Lm t uN


eq1_NL = (mb+Jb/r^2)*xb_ddot + b*xb_dot - mb*g*sin(alpha) == 0;                                         % Equação da bola
eq2_NL = Jl*alpha_ddot + ml*g*cos(alpha)*(l/2) + mb*g*cos(alpha)*(l-xb) + f*l*cos(alpha) == 0;          % Equação da viga
eq3_NL = l*sin(alpha) - d*sin(theta) == 0;                                                              % Equação da relação entre theta e alpha
eq4_NL = Jm*theta_ddot + Bm*theta_dot - f*cos(theta)*d - Km*i == 0;                                     % Equação do motor de corrente contínua
eq5_NL = uN - Rm*i - Lm*diff(i,t) - Km*theta_dot == 0;                                                    % Equação do motor de corrente contínua


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
%cos(alpha) = sqrt(delta);
%alpha_ddot = (k*sin(theta)*theta_ddot^2*(k^2-1) + k*cos(theta)*theta_ddot*delta)/(delta^(3/2));
%sin(alpha) = k*sin(theta);

%Lm = 0;                                                                          % Aproximação dada pelo enunciado
%i  = (uN - Km*theta_dot)/Rm;
% Isolamento de i a partir da equação 5
i_eq_5 = solve(subs(eq5_NL, Lm, 0), i);

% Aplicação da simplificação nas equações iniciais
eq1_NL_simp = subs(eq1_NL, sin(alpha), k*sin(theta));
eq2_NL_simp = subs(eq2_NL, [alpha_ddot, cos(alpha)], [(k*sin(theta)*theta_dot^2*(k^2-1) + k*cos(theta)*theta_ddot*delta)/(delta^(3/2)), sqrt(delta)]);

% Isolando f de eq2_NL_simp
f_eq_2 = solve(eq2_NL_simp, f);

% Isolando xb_ddot na eq1_NL_simp
xb_ddot_NL = solve(eq1_NL_simp, xb_ddot);

% Isolando theta_ddot na eq4_NL, substituindo f por f_eq_2 e i por i_eq_5
theta_ddot_eq_2 = solve(subs(eq4_NL, [f, i], [f_eq_2, i_eq_5]), theta_ddot);

% Representação em Espaços de Estados

z1_dot  = xb_dot;
z2_dot  = xb_ddot_NL; % = xb_ddot 
z3_dot  = theta_dot; % = theta_dot  
z4_dot  = theta_ddot_eq_2;

Z_dot = [z1_dot;
         z2_dot;
         z3_dot;
         z4_dot];


% Convertendo Z_dot para uma função do MATLAB
Z_dot_func = matlabFunction(Z_dot, 'Vars', {xb, xb_dot, theta, theta_dot, uN});

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

[num, den] = numden(G);          % Pegar o numerador e denominador

% coeficiente do termo de maior grau
coef_lider = coeffs(den, s, 'All');
fator_divisao = coef_lider(1); 
num_norm = collect(num / fator_divisao, s);
den_norm = collect(den / fator_divisao, s);

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

% Parâmetros dos Polos Dominantes
wn = 0.45;          % Frequência natural (velocidade do sistema)
zeta = 0.75;      % Fator de amortecimento (mínimo de projeto é 0.3 para Margem de Fase >=30°)
P_dominante = s^2 + 2*zeta*wn*s + wn^2;

% Polos Não-Dominantes (devem ser mais rápidos/distantes da origem)
P_rapidos = (s+9)^2 * (s+12)^2 * (s+15)^2 * (s+18); % Inicial P_rapidos = (s+15)^2 * (s+18)^2 * (s+20)^2 * (s+25); u_N = 97.8196 V

% Polinômio de Malha Fechada Final (Grau 9)
P = P_dominante * P_rapidos;

% Expande o polinômio P e extrai os coeficientes
P = expand(P);
pv = coeffs(P, s, 'All')';

% Expande o polinômio do denominador da planta D_G e obtém-se os coeficientes
D_G = expand(den_sym*s); % Multiplica-se por s para garantir um sistema do tipo 1
dgv = coeffs(D_G, s, 'All');

% Expande o polinômio do numerador da planta N_G e obtém-se os coeficientes
N_G = expand(num_sym);
ngv = coeffs(N_G, s, 'All');

% Definindo a matriz de Silvester para os polinômios numerador e denominador
n_C = 4; % Ordem do denominador do controlador
m_C = 4; % Ordem do numerador do controlador
S = matriz_sylvester_controle(dgv, ngv, n_C, m_C);

cv = S \ pv(end:-1:1);

% Conversão dos coeficientes no controlador C_0
C_0 = poly2sym(cv(end:-1:n_C + 2), s) / poly2sym(cv(n_C + 1:-1:1), s);

% Cálculo do controlador C com base em C_0 / s
C = vpa(C_0 / s, 4);

%% Verificação se o controlador está de acordo com as especificações
% CÁLCULO DA MARGEM DE FASE E FREQUÊNCIA DE CRUZAMENTO
[num_C, den_C] = numden(C);
num_C_num = sym2poly(num_C);
den_C_num = sym2poly(den_C);

C_tf = tf(num_C_num, den_C_num);
Malha_Aberta = C_tf * LR_G;

[MG, MF, Wcg, Om_g] = margin(Malha_Aberta);

fprintf('\n--- RESULTADOS DA MALHA ABERTA ---\n');
fprintf('Frequência de cruzamento de ganho (Om_g): %.4f rad/s\n', Om_g);
fprintf('Margem de Fase (MF): %.4f graus\n', MF);

figure;
margin(Malha_Aberta);
grid on;
title('Diagrama de Bode - Sistema em Malha Aberta (C * G)');

figure;
nyquist(Malha_Aberta);
grid on;
title('Diagrama de Nyquist - Sistema em Malha Aberta (C * G)');

Malha_Fechada = feedback(Malha_Aberta, 1);

polos_MF = pole(Malha_Fechada);

fprintf('\n--- POLOS DE MALHA FECHADA ---\n');
disp(polos_MF);

fprintf('\n--- ANÁLISE DE AMORTECIMENTO (DAMP) ---\n');
damp(Malha_Fechada);

% Análise de Esforço de Controlo 
% A função de transferência de R(s) para U(s) é C(s) / (1 + C(s)G(s))
% O comando feedback(C, G) faz exatamente a malha: C / (1 + C*G)
Tr_u = feedback(C_tf, LR_G);

% Definindo a amplitude da referência r_chapeu = 0.02
opt = stepDataOptions('StepAmplitude', 0.02);

% Simulando o esforço de controlo no tempo (vamos olhar os primeiros 20 segundos)
[u_t, t_u] = step(Tr_u, 10, opt);

% Encontrando o pico máximo absoluto do esforço de controlo
pico_u = max(abs(u_t));

% Exibindo o resultado na Command Window
fprintf('\n--- VERIFICAÇÃO DO ESFORÇO DE CONTROLO ---\n');
fprintf('Pico de tensão no motor (max |u_N|): %.4f V\n', pico_u);

if pico_u <= 50
    fprintf('STATUS: APROVADO! O pico é menor que 50 V.\n');
else
    fprintf('STATUS: REPROVADO! O pico ultrapassou 50 V.\n');
end

% Plotando o gráfico do Esforço de Controlo para o Relatório
figure;
plot(t_u, u_t, 'b', 'LineWidth', 1.5);
hold on;
yline(50, 'r--', 'Limite de +50V', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
yline(-50, 'r--', 'Limite de -50V', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
grid on;
title('Esforço de Controlo u_N(t) para Degrau de Referência de 0.02m');
xlabel('Tempo (s)');
ylabel('Tensão no Motor - u_N (V)');

% Plotando a resposta da malha fechada em relação à entrada de referência
% de 0,02
% Simulando a resposta da malha fechada a um degrau de referência
[y_t, t] = step(Malha_Fechada, 20, opt);

% Plotando a resposta da malha fechada
figure;
plot(t, y_t, 'r', 'LineWidth', 1.5);
grid on;
title('Resposta da Malha Fechada a um Degrau de Referência de 0.02m');
xlabel('Tempo (s)');
ylabel('Saída do Sistema - y(t)');

%% Simulação da aplicação de C(s) para a representação em espaço de estados não linear

% 1. Transformação para Espaço de Estados do Controlador
[A_C, B_C, C_C, D_C] = tf2ss(num_C_num, den_C_num);

% 2. Definição das Condições Iniciais
z_0 = [0.5; 0; 0; 0]; % Condição inicial da planta (já como vetor coluna)
x_C_0 = zeros(size(A_C, 1), 1); % Condição inicial do controlador
zeta_0 = [z_0; x_C_0]; % Estado estendido

% 3. Integração Numérica
t_span = [0 20]; % Intervalo de 20 segundos

% Passamos A_C, B_C, C_C, e D_C como argumentos fixos para a função
[t_out, zeta_out] = ode45(@(t, zeta) rep_nao_linear(t, zeta, Z_dot_func, A_C, B_C, C_C, D_C), t_span, zeta_0);

% 4. Plotagem do resultado
figure(5);
plot(t_out, zeta_out(:, 1), 'g', 'LineWidth', 1.5);
title('Resposta do Sistema Não Linear em Malha Fechada');
xlabel('Tempo (s)');
ylabel('Posição zeta_1 (m)');
grid on;

%% -------------------------------------------------------------------------
% Função de Equações Diferenciais (Otimizada para Velocidade e Precisão)
function zeta_dot = rep_nao_linear(t, zeta, Z_dot_func, A_C, B_C, C_C, D_C)
    
    % Inicializa o vetor coluna de derivadas
    zeta_dot = zeros(length(zeta), 1);
    
    % Separação dos estados para clareza
    x_planta = zeta(1:4);
    x_controlador = zeta(5:end);
    
    % Cálculo do Erro (Referência - Posição Atual)
    ref = 0.02;
    erro = ref - x_planta(1);
    
    % Cálculo da Ação de Controle: u = C_c * x_c + D_c * erro
    u = C_C * x_controlador + D_C * erro;
    
    % Derivadas da Planta Não Linear
    % Passamos os 4 estados e a entrada de controle calculada
    zeta_dot(1:4, 1) = Z_dot_func(x_planta(1), x_planta(2), x_planta(3), x_planta(4), u);
    
    % Derivadas do Controlador Linear: x_c_dot = A_c * x_c + B_c * erro
    zeta_dot(5:end, 1) = A_C * x_controlador + B_C * erro;
    
end
