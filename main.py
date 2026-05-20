import math
import matplotlib.pyplot as plt


# ============================================================
# 1. PROPRIEDADES FÍSICAS
# ============================================================

k = 3.0
rho = 10500.0
Cp = 300.0

alpha = k / (rho * Cp)

# ============================================================
# 2. GEOMETRIA
# ============================================================

L = 0.005

# ============================================================
# 3. CONVECÇÃO
# ============================================================

h = 10000.0
T_inf = 300.0

# ============================================================
# 4. CONDIÇÃO INICIAL
# ============================================================

T_i = 300.0

# ============================================================
# 5. MALHA NUMÉRICA
# ============================================================

N_nodes = 21

dx = L / (N_nodes - 1)

# ============================================================
# 6. TEMPO
# ============================================================

dt = 0.1

t_final = 20.0

# ============================================================
# 7. PARÂMETROS ADIMENSIONAIS
# ============================================================

Fo = (alpha * dt) / (dx**2)

Bi_mesh = (h * dx) / k

Bi_global = (h * L) / k

# ============================================================
# 8. ALGORITMO DE THOMAS
# ============================================================

def thomas_algorithm(e, f, g, r):

    n = len(f)

    x = [0.0] * n

    e_c = list(e)
    f_c = list(f)
    r_c = list(r)

    # --------------------------------------------
    # DECOMPOSIÇÃO LU
    # --------------------------------------------

    for i in range(1, n):

        e_c[i] = e_c[i] / f_c[i - 1]

        f_c[i] = f_c[i] - e_c[i] * g[i - 1]

    # --------------------------------------------
    # SUBSTITUIÇÃO PROGRESSIVA
    # --------------------------------------------

    for i in range(1, n):

        r_c[i] = r_c[i] - e_c[i] * r_c[i - 1]

    # --------------------------------------------
    # SUBSTITUIÇÃO REGRESSIVA
    # --------------------------------------------

    x[n - 1] = r_c[n - 1] / f_c[n - 1]

    for i in range(n - 2, -1, -1):

        x[i] = (r_c[i] - g[i] * x[i + 1]) / f_c[i]

    return x

# ============================================================
# 9. SOLUÇÃO ANALÍTICA
# ============================================================

def obter_raizes(Bi, quantidade=6):

    raizes = []

    passo = 0.005

    for n in range(quantidade):

        inicio = max(0.001, n * math.pi)

        fim = inicio + math.pi / 2 - 0.01

        valor = inicio

        while valor < fim:

            f1 = valor * math.tan(valor) - Bi

            f2 = (valor + passo) * math.tan(valor + passo) - Bi

            if f1 * f2 < 0:

                a = valor
                b = valor + passo

                for _ in range(40):

                    c = (a + b) / 2

                    if (c * math.tan(c) - Bi) * (a * math.tan(a) - Bi) < 0:
                        b = c
                    else:
                        a = c

                raizes.append((a + b) / 2)

                break

            valor += passo

    return raizes

lambdas = obter_raizes(Bi_global)

# ============================================================
# 10. SOLUÇÃO ANALÍTICA
# ============================================================

def solucao_analitica(x, t):

    if t == 0:
        return T_i

    Fo_global = (alpha * t) / (L**2)

    soma = 0.0

    x_star = x / L

    for lam in lambdas:

        numerador = 4.0 * math.sin(lam)

        denominador = 2.0 * lam + math.sin(2.0 * lam)

        coeficiente = numerador / denominador

        soma += (
            coeficiente
            * math.cos(lam * x_star)
            * math.exp(-(lam**2) * Fo_global)
        )

    return T_inf + (T_i - T_inf) * soma

# ============================================================
# 11. SOLVER NUMÉRICO
# ============================================================

def resolver_sistema(q_dot_val):

    T = [T_i] * N_nodes

    historico = [list(T)]

    A_termo = (q_dot_val * dt) / (rho * Cp)

    e = [0.0] * N_nodes
    f = [0.0] * N_nodes
    g = [0.0] * N_nodes

    # --------------------------------------------
    # NÓ 1
    # --------------------------------------------

    f[0] = 1.0 + 2.0 * Fo

    g[0] = -2.0 * Fo

    # --------------------------------------------
    # NÓS INTERNOS
    # --------------------------------------------

    for i in range(1, N_nodes - 1):

        e[i] = -Fo

        f[i] = 1.0 + 2.0 * Fo

        g[i] = -Fo

    # --------------------------------------------
    # NÓ FINAL
    # --------------------------------------------

    e[-1] = -2.0 * Fo

    f[-1] = 1.0 + 2.0 * Fo + 2.0 * Fo * Bi_mesh

    passos_tempo = int(t_final / dt)

    for p in range(passos_tempo):

        b_vetor = [0.0] * N_nodes

        b_vetor[0] = T[0] + A_termo

        for i in range(1, N_nodes - 1):

            b_vetor[i] = T[i] + A_termo

        b_vetor[-1] = (
            T[-1]
            + A_termo
            + 2.0 * Fo * Bi_mesh * T_inf
        )

        T_novo = thomas_algorithm(e, f, g, b_vetor)

        T = list(T_novo)

        historico.append(T)

    return historico

# ============================================================
# 12. EXECUÇÃO DAS SIMULAÇÕES
# ============================================================

X_posicoes = [i * dx for i in range(N_nodes)]

X_mm = [x * 1000 for x in X_posicoes]

# --------------------------------------------
# VALIDAÇÃO
# --------------------------------------------

print('\nExecutando validação numérica...\n')

resultado_validacao = resolver_sistema(q_dot_val=0.0)

# --------------------------------------------
# CASO COM GERAÇÃO INTERNA
# --------------------------------------------

q_reator = 2.0e8

print('Executando caso com geração interna...\n')

resultado_operacional = resolver_sistema(q_dot_val=q_reator)

# ============================================================
# 13. GRÁFICOS
# ============================================================

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# ============================================================
# VALIDAÇÃO ANALÍTICA
# ============================================================

T_num = resultado_validacao[-1]

T_exato = [
    solucao_analitica(x, t_final)
    for x in X_posicoes
]

ax1.plot(
    X_mm,
    T_num,
    'go',
    label='Numérico',
    markersize=6
)

ax1.plot(
    X_mm,
    T_exato,
    'k--',
    label='Analítico',
    linewidth=1.5
)

ax1.set_title('Validação Numérica')

ax1.set_xlabel('Posição x [mm]')

ax1.set_ylabel('Temperatura [°C]')

ax1.legend()

ax1.grid(True)

# ============================================================
# EVOLUÇÃO TRANSIENTE
# ============================================================

instantes = [0, 1, 3, 5, 10, 20]

for t_alvo in instantes:

    indice = int(t_alvo / dt)

    ax2.plot(
        X_mm,
        resultado_operacional[indice],
        label=f't = {t_alvo} s'
    )

ax2.set_title('Perfis de Temperatura')

ax2.set_xlabel('Posição x [mm]')

ax2.set_ylabel('Temperatura [°C]')

ax2.legend()

ax2.grid(True)

plt.tight_layout()

plt.savefig('resultados_ppc3.png', dpi=300)

print('Gráfico salvo com sucesso.')

plt.show()

# ============================================================
# 14. TEMPERATURAS FINAIS
# ============================================================

print('\nTemperaturas finais:\n')

T_final = resultado_operacional[-1]

for i in range(N_nodes):

    print(f'Nó {i+1}: {T_final[i]:.2f} °C')
