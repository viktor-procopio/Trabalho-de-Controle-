function S = matriz_sylvester(p, q, var)
    % matriz_sylvester Cria a matriz de Sylvester para dois polinômios.
    %
    % Entradas:
    % p, q - Polinômios simbólicos
    % var  - Variável simbólica em relação à qual os polinômios são definidos
    %
    % Saída:
    % S    - Matriz de Sylvester (simbólica)

    % Extrai os coeficientes dos polinômios, incluindo os termos nulos (graus ausentes).
    % O parâmetro 'All' garante que todos os coeficientes sejam retornados,
    % ordenados do termo de maior grau para o de menor grau.
    cp = coeffs(p, var, 'All');
    cq = coeffs(q, var, 'All');

    % O grau do polinômio é o número de coeficientes menos 1
    n = length(cp) - 1; % Grau de p
    m = length(cq) - 1; % Grau de q

    % Inicializa a matriz de Sylvester com zeros simbólicos
    % A dimensão da matriz será (n + m) x (n + m)
    S = sym(zeros(n + m));

    % Preenche as primeiras 'm' linhas com os coeficientes de 'p'
    for i = 1:m
        S(i, i:(i + n)) = cp(end:-1:1);
    end

    % Preenche as 'm' linhas seguintes com os coeficientes de 'q'
    for j = 1:n
        S(m + j, j:(j + m)) = cq(end:-1:1);
    end
end
