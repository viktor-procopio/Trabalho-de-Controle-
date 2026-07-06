% ========================================================================
% PROJETO BOLA-VIGA - CORREÇÃO DE DIMENSÃO DO LSIM (ITENS 10 A 13)
% ========================================================================
clear; clc; close all;

% ========================================================================
% DEFINIÇÃO PARÂMETROS FÍSICOS (Conforme Tabela 1 do roteiro)
% ========================================================================
mb = 0.01; b = 2e-3; r = 0.02; ml = 0.1; l = 1; Jm = 4e-4; 
Jb = (2/5)*mb*r^2; Jl = ml*l^2/3; Bm = 1.1e-2; Rm = 1; Km = 0.03; 
d = 0.04; g = 9.8;

% Ponto de Operação Absoluto
xb_o = 0.5; 
theta_o = 0;

% ========================================================================
% MODELAGEM SIMBÓLICA NÃO-LINEAR (Item 4)
% ========================================================================
syms xb xb_dot theta theta_dot uN

k_param = d/l; 
delta = 1 - k_param^2 * (sin(theta))^2;

cos_alpha = sqrt(delta);
sin_alpha = k_param * sin(theta);

alpha_ddot_v = (k_param * sin(theta) * (k_param^2 - 1)) / (delta^(3/2));
alpha_ddot_a = (k_param * cos(theta) * delta) / (delta^(3/2)); 

i_arm = (uN - Km*theta_dot) / Rm;

M_bola = mb + Jb/r^2;
xb_ddot_NL = (mb*g*sin_alpha - b*xb_dot) / M_bola;

Termo_J_viga = (Jl * alpha_ddot_a * cos(theta) * d) / (l * cos_alpha);
J_total_NL = Jm + Termo_J_viga;

f_resto = -(Jl*alpha_ddot_v*theta_dot^2 + ml*g*cos_alpha*(l/2) + mb*g*cos_alpha*(l-xb)) / (l*cos_alpha);
Torques_resto = - Bm*theta_dot + f_resto*cos(theta)*d + Km*i_arm;

theta_ddot_NL = Torques_resto / J_total_NL;

Z_dot = [xb_dot; xb_ddot_NL; theta_dot; theta_ddot_NL];
Z_dot_func = matlabFunction(Z_dot, 'Vars', {xb, xb_dot, theta, theta_dot, uN});

% ========================================================================
% 10. CONTROLADOR (PARTE 3 - ITEM 10)
% ========================================================================
M_eq = mb + Jb/r^2; 
J_eq = Jm + (Jl * d^2) / l^2; 
B_eq = (Bm + Km^2/Rm);

A = [0, 1, 0, 0; 
     0, -b/M_eq, (mb*g*d)/(l*M_eq), 0; 
     0, 0, 0, 1; 
     (mb*g*d)/(l*J_eq), 0, 0, -B_eq/J_eq];
B = [0; 0; 0; Km/(Rm*J_eq)];
C = [1, 0, 0, 0]; 
D = 0;

% Requisitos de Polos Dominantes (TS < 4s, xi = 0.6, wn = 2.5)

% Calculando os polos com LQR
R = 1/50; % Ro definido para garantir entrada menor que 50 V
Q1 = 2*C'*C; % Definição de Q
[K,~,polos_planta] = lqr(A,B,Q1,R); % Cálculo do ganho matricial K e dos polos da planta
M = zeros(size(K, 2), 1);
M(1) = 1 / (C * inv(B*K - A) * B * K(1));

% Determinação de uo de equilíbrio estático
f_o = (mb*g*(xb_o - l) - ml*g*(l/2))/l;
i_o = -(f_o*d)/Km;
uo = Rm*i_o;

% Simulação Linear Incremental Original do Item 10 via lsim
t_base = 0:0.01:6; t = t_base(:); r_mag = 0.02;
sys_mf = ss(A - B*K, B*K*M, C, 0); % CORRIGIDO: B*M garante entrada escalar simples (1 canal)
[y_inc, ~, x_inc] = lsim(sys_mf, r_mag*ones(size(t)), t);

x_b_absoluto_10 = xb_o + y_inc;
u_N_absoluto_10 = (K * (M * r_mag - x_inc'))' + uo;

% --- Passo 5: Análise Automática de Requisitos ---
fprintf('\n--- Passo 5: Análise Automática de Requisitos ---\n');
min_xb = min(x_b_absoluto_10); % Ajustado para _10
max_xb = max(x_b_absoluto_10); % Ajustado para _10

if min_xb >= 0 && max_xb <= 1
    fprintf('[OK] Restrição da Barra Atendida! Posição ficou entre %0.3fm e %0.3fm.\n', min_xb, max_xb);
else
    fprintf('[FALHA] A bola caiu da barra! Posição extrapolou limites (Min: %0.3fm, Max: %0.3fm).\n', min_xb, max_xb);
end

max_tensao = max(abs(u_N_absoluto_10)); % Ajustado para _10
if max_tensao <= 50
    fprintf('[OK] Restrição de Tensão Atendida! O pico foi de %0.2f Volts.\n', max_tensao);
else
    fprintf('[FALHA] O motor saturou! Pico de tensão calculado de %0.2f Volts excedeu os 50V.\n', max_tensao);
end
disp(max_tensao)

% Plotagem Gráfica Exclusiva do Item 10 (Apenas Linear)
figure('Name', 'Item 10 - Apenas Linear', 'NumberTitle', 'off');
subplot(2,1,1); plot(t, x_b_absoluto_10, 'b', 'LineWidth', 2); hold on;
yline(xb_o + r_mag, 'k--', 'Referência Alvo (0.52m)', 'LineWidth', 1.2);
title('Resposta Temporal da Posição Absoluta da Bola $x_b(t)$ (Item 10)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Posição $x_b$ (m)', 'Interpreter', 'latex'); grid on;

subplot(2,1,2); plot(t, u_N_absoluto_10, 'g', 'LineWidth', 2); hold on;
yline(50, 'r--', 'Limite Superior (+50V)', 'LineWidth', 1.2);
yline(-50, 'r--', 'Limite Inferior (-50V)', 'LineWidth', 1.2);
title('Esforço de Controle Absoluto - Tensão Aplicada $u_N(t)$ (Item 10)', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Tensão $u_N$ (V)', 'Interpreter', 'latex'); grid on;
ylim([-60, 60]);

% Configuração das tolerâncias do ODE45
opcoes_ode45 = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep', 0.005);
% ========================================================================
% CÁLCULO DA MARGEM DE FASE (Análise de Frequência do Item 10)
% ========================================================================
% Criação do modelo de malha aberta: L(s) = K*(sI - A)^-1*B
sys_ol = ss(A, B, K, 0); 

% A função 'margin' calcula as margens de ganho e fase e as frequências de cruzamento
[Gm, Pm, Wcg, Wcp] = margin(sys_ol);

% Exibição dos resultados no Command Window
fprintf('\n--- Análise de Estabilidade no Domínio da Frequência ---\n');
fprintf('Margem de Ganho (Gm): %.2f (Absoluto) -> %.2f dB\n', Gm, 20*log10(Gm));
fprintf('Margem de Fase (Pm): %.2f graus\n', Pm);
fprintf('Frequência de Cruzamento de Fase (Wcg): %.2f rad/s\n', Wcg);
fprintf('Frequência de Cruzamento de Ganho (Wcp): %.2f rad/s\n', Wcp);

% Plotagem do Diagrama de Bode para visualização das Margens
figure('Name', 'Margens de Estabilidade (Item 10)', 'NumberTitle', 'off');
margin(sys_ol);
% ========================================================================
% 11. ANÁLISE DE DESEMPENHO (ITEM 11)
% ========================================================================
fprintf('\n--- Executando Integração Numérica ODE45 (Item 11) ---\n');
xi0_abs = [xb_o; 0; 0; 0];

[t_nl_11, X_nl_11] = ode45(@(t, xi) dinamica_planta_simbolica(t, xi, Z_dot_func, K, M, r_mag, uo, xb_o), t_base, xi0_abs, opcoes_ode45);

x_b_nao_linear_11 = X_nl_11(:, 1);
X_inc_nl_11 = [X_nl_11(:,1)-xb_o, X_nl_11(:,2), X_nl_11(:,3), X_nl_11(:,4)];
u_N_nao_linear_11 = K * (M * r_mag - X_inc_nl_11') + uo;

u_N_nao_linear_11(u_N_nao_linear_11 > 50) = 50;
u_N_nao_linear_11(u_N_nao_linear_11 < -50) = -50;

figure('Name', 'Item 11(b) - Posição da Bola', 'NumberTitle', 'off');
plot(t, x_b_absoluto_10, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_11, x_b_nao_linear_11, 'b-', 'LineWidth', 2);
yline(xb_o + r_mag, 'k--', 'Referência Desejada (0.52m)', 'LineWidth', 1.2);
title('Item 11: Posição Absoluta da Bola $x_b(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Posição (m)');
legend('Modelo Linearizado', 'Modelo Não-Linear Real', 'Location', 'Best'); grid on;

figure('Name', 'Item 11(b) - Esforço de Controle', 'NumberTitle', 'off');
plot(t, u_N_absoluto_10, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_11, u_N_nao_linear_11, 'g-', 'LineWidth', 1.5);
yline(50, 'k--', 'Limite Máximo (+50V)');
yline(-50, 'k--', 'Limite Mínimo (-50V)');
title('Item 11: Esforço de Controle $u_N(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Tensão (V)');
legend('Modelo Linearizado', 'Modelo Não-Linear Real', 'Location', 'Best'); grid on;

% ========================================================================
% 12. PROJETO DO OBSERVADOR (ITEM 12)
% ========================================================================
fator_rapidez = 6; 
polo_base_obs = fator_rapidez * max(polos_planta); 
polos_observador = [polo_base_obs, polo_base_obs-1, polo_base_obs-2, polo_base_obs-3];
L = place(A', C', polos_observador)';

% ========================================================================
% 13. SIMULAÇÃO COM OBSERVADOR EM MALHA FECHADA (ITEM 13)
% ========================================================================
fprintf('\n--- Executando Integração Numérica ODE45 (Item 13) ---\n');
xi0_comp = [xb_o - 0.01; 0; 0; 0]; 
xo0_comp = [0; 0; 0; 0];           
X0_composto = [xi0_comp; xo0_comp];     

% --- SIMULAÇÃO LINEAR COMPREENSIVA VIA LSIM ---
x0_lin_13 = [-0.01; 0; 0; 0; 0; 0; 0; 0]; 
A_comp_lin = [A, -B*K; L*C, A - B*K - L*C];
B_comp_lin = [B*K*M; B*K*M]; % CORRIGIDO: Entrada escalar única coerente com 1 canal
sys_lin_13 = ss(A_comp_lin, B_comp_lin, [C, zeros(1,4)], 0);
[Y_lin_13, t_lin_13, X_lin_13] = lsim(sys_lin_13, r_mag*ones(size(t_base)), t_base, x0_lin_13);
x_b_linear_13 = xb_o + Y_lin_13;
u_N_linear_13 = (K * (M * r_mag - X_lin_13(:, 5:8)'))' + uo;

% --- SIMULAÇÃO NÃO-LINEAR REAL VIA ODE45 ---
f_composto = @(t, X) dinamica_composta_simbolica(t, X, Z_dot_func, A, B, C, K, L, M, r_mag, uo, xb_o);
[t_nl_13, X_nl_13] = ode45(f_composto, t_base, X0_composto, opcoes_ode45);

x_b_nao_linear_13 = X_nl_13(:, 1);                       
X_obs_13 = X_nl_13(:, 5:8);
u_N_nao_linear_13 = (K * (M * r_mag - X_obs_13'))' + uo;     

u_N_nao_linear_13(u_N_nao_linear_13 > 50) = 50;
u_N_nao_linear_13(u_N_nao_linear_13 < -50) = -50;

figure('Name', 'Item 13(b) - Posição y(t) com Observador Corrigido', 'NumberTitle', 'off');
plot(t_lin_13, x_b_linear_13, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_13, x_b_nao_linear_13, 'b-', 'LineWidth', 2);
yline(xb_o + r_mag, 'k--', 'Referência Desejada (0.52m)', 'LineWidth', 1.2);
title('Item 13: Comparação da Posição Absoluta da Bola $x_b(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Posição (m)');
legend('Modelo Linearizado', 'Modelo Não-Linear Real', 'Location', 'Best'); grid on;

figure('Name', 'Item 13(b) - Esforço u_N(t) com Observador Corrigido', 'NumberTitle', 'off');
plot(t_lin_13, u_N_linear_13, 'r--', 'LineWidth', 1.5); hold on;
plot(t_nl_13, u_N_nao_linear_13, 'g-', 'LineWidth', 2);
yline(50, 'k--', 'Limite Máximo (+50V)');
yline(-50, 'k--', 'Limite Mínimo (-50V)');
title('Item 13: Comparação do Esforço de Controle $u_N(t)$', 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Tempo (s)'); ylabel('Tensão (V)');
legend('Modelo Linearizado', 'Modelo Não-Linear Real', 'Location', 'Best'); grid on;

% ========================================================================
% FUNÇÕES AUXILIARES
% ========================================================================
function dxi = dinamica_planta_simbolica(~, xi, Z_dot_func, K, M, r_mag, uo, xb_o)
    x_inc = [xi(1) - xb_o; xi(2); xi(3); xi(4)];
    u_N = K * (M * r_mag - x_inc) + uo;
    if u_N > 50,  u_N = 50;  end
    if u_N < -50, u_N = -50; end
    dxi = Z_dot_func(xi(1), xi(2), xi(3), xi(4), u_N);
end

% Nota explicativa sobre a modelagem incremental do observador:
% O observador estima puramente os desvios incrementais em relação ao ponto de equilíbrio.
function dX = dinamica_composta_simbolica(~, X, Z_dot_func, A, B, C, K, L, M, r_mag, uo, xb_o)
    xi = X(1:4);  
    xo = X(5:8);  
    
    u_N = K * (M * r_mag - xo) + uo;
    if u_N > 50,  u_N = 50;  end
    if u_N < -50, u_N = -50; end
    
    dxi = Z_dot_func(xi(1), xi(2), xi(3), xi(4), u_N);
    
    y_incremental_real = xi(1) - xb_o; 
    u_incremental = u_N - uo;
    dxo = A * xo + B * u_incremental + L * (y_incremental_real - C * xo);
    
    dX = [dxi; dxo];
end
