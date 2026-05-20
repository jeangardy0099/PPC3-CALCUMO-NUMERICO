# Descrição

Este trabalho apresenta a solução numérica do problema de distribuição de temperatura transiente em uma parede plana unidimensional utilizando o Método das Diferenças Finitas.

O problema foi resolvido utilizando:

- esquema implícito;
- algoritmo de Thomas (TDMA);
- discretização espacial e temporal;
- solução analítica para validação numérica.

---

# Equação Governante

\[
\frac{\partial T}{\partial t}
=
\alpha
\frac{\partial^2 T}{\partial x^2}
+
\frac{\dot q}{\rho C_p}
\]

onde:

- \(T\) = temperatura;
- \(\alpha\) = difusividade térmica;
- \(\dot q\) = geração interna de calor;
- \(\rho\) = massa específica;
- \(C_p\) = calor específico.

---

# Métodos Utilizados

O programa foi desenvolvido utilizando:

- Método das Diferenças Finitas;
- Esquema Implícito;
- Algoritmo de Thomas para sistemas tridiagonais;
- Série de Fourier para validação analítica.

O algoritmo de Thomas foi implementado manualmente, sem utilização de bibliotecas externas de integração numérica, conforme solicitado na atividade.

---

# Bibliotecas Utilizadas

O programa utiliza apenas bibliotecas básicas:

- math
- matplotlib

---

# Estrutura do Código

O código foi dividido nas seguintes etapas:

1. Definição das propriedades físicas;
2. Definição da malha numérica;
3. Implementação do algoritmo de Thomas;
4. Implementação da solução analítica;
5. Solver numérico implícito;
6. Pós-processamento dos resultados;
7. Geração de gráficos.

---

# Como Executar

## 1. Instalar dependências

```bash
pip install matplotlib
