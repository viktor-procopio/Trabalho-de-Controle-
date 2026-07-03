function S = matriz_sylvester_controle(D_G, N_G, n_C, m_C)
    % matriz_sylvester_controle: Formula a Matriz de Sylvester generalizada
    % para a Equação Diofantina Polinomial D_G*X + N_G*Y = C_d
    %
    % Entradas:
    % D_G - Vetor com coeficientes do denominador da planta (maior p/ menor grau)
    % N_G - Vetor com coeficientes do numerador da planta (maior p/ menor grau)
    % n_C - Ordem desejada para o denominador do controlador, X(s)
    % m_C - Ordem desejada para o numerador do controlador, Y(s)
    %
    % Saída:
    % S   - Matriz de Sylvester

    % Extraindo os graus dos polinômios da planta (tamanho do vetor - 1)
    n_G = length(D_G) - 1;
    m_G = length(N_G) - 1;

    % O grau do polinômio de malha fechada C_d(s) será o maior grau resultante
    % da soma das parcelas D_G(s)X(s) e N_G(s)Y(s)
    grau_malha_fechada = max(n_G + n_C, m_G + m_C);

    % O número de linhas (equações) é o grau da malha fechada + 1
    num_linhas = grau_malha_fechada + 1;

    % O número de colunas (incógnitas) são os coeficientes de X(s) e Y(s)
    num_colunas = (n_C + 1) + (m_C + 1);

    % Inicializa a matriz de Sylvester com zeros
    S = zeros(num_linhas, num_colunas);

    % 1. Preenchendo o bloco esquerdo (relativo a D_G multiplicado por X)
    % Haverá (n_C + 1) colunas deslocadas
    for i = 1:(n_C + 1)
        % Insere os coeficientes de D_G na coluna 'i', deslocando 'i-1' zeros acima
        S(i : i + n_G, i) = D_G(end:-1:1); 
    end

    % 2. Preenchendo o bloco direito (relativo a N_G multiplicado por Y)
    % Haverá (m_C + 1) colunas deslocadas
    offset_col = n_C + 1; % Para não sobrescrever o bloco de D_G
    for j = 1:(m_C + 1)
        % Insere os coeficientes de N_G na coluna 'offset_col + j'
        S(j : j + m_G, offset_col + j) = N_G(end:-1:1);
    end
end
