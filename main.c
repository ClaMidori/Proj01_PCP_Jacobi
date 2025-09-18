#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <omp.h>

typedef struct {
    const char* name;
    int threads;
    const char* schedule;
    int chunk;
} Config;

#define LIMITE_ERRO 1e-5
#define LIMITE_ITERACOES 1000000

double** aloca_matriz(int n) {
    double** mat = (double**) malloc(n * sizeof(double*));
    if (!mat) { perror("malloc"); exit(1); }
    for (int i = 0; i < n; i++) {
        mat[i] = (double*) malloc(n * sizeof(double));
        if (!mat[i]) { perror("malloc row"); exit(1); }
    }
    return mat;
}

void libera_matriz(double** mat, int n) {
    if (!mat) return;
    for (int i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

double** copiar_matriz(double** original, int n) {
    double** copia = aloca_matriz(n);
    for (int i = 0; i < n; i++)
        memcpy(copia[i], original[i], n * sizeof(double));
    return copia;
}

double* copiar_vetor(double* original, int n) {
    double* copia = (double*) malloc(n * sizeof(double));
    memcpy(copia, original, n * sizeof(double));
    return copia;
}

// pivo para evitar diagonal zero
int ensure_nonzero_diagonal(double** A, double* b, int n) {
    for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) > 1e-14) continue;
        int found = -1;
        for (int r = i+1; r < n; r++) {
            if (fabs(A[r][i]) > 1e-14) { found = r; break; }
        }
        if (found == -1) {
            for (int r = 0; r < n; r++) {
                if (r!=i && fabs(A[r][i]) > 1e-14) { found = r; break; }
            }
        }
        if (found == -1) return 0; 
        // troca linha i com found
        double* tmp = A[i];
        A[i] = A[found];
        A[found] = tmp;
        double tb = b[i]; b[i] = b[found]; b[found] = tb;
    }
    return 1;
}

//pré-condicionamento diagonal
// retorna 1 se OK, 0 se diagonal for muito pequena
int precondiciona_diagonal(double** A, double* b, int n, int verbose) {
    // primeiro garanta diagonais não-zero (tenta trocar linhas)
    if (!ensure_nonzero_diagonal(A, b, n)) {
        if (verbose) printf("Precond: não foi possível garantir diagonal não-zero.\n");
        return 0;
    }
    for (int i = 0; i < n; i++) {
        double d = A[i][i];
        if (fabs(d) < 1e-14) return 0;
        double inv = 1.0 / d;
        for (int j = 0; j < n; j++) A[i][j] *= inv;
        b[i] *= inv;
    }
    return 1;
}

//verificaça se é diagonalmente dominante
int eh_diagonalmente_dominante(double** A, int n) {
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        for (int j = 0; j < n; j++)
            if (i != j) soma += fabs(A[i][j]);
        if (fabs(A[i][i]) <= soma) return 0;
    }
    return 1;
}

// reordena linhas das equações
void reordenar_para_dominancia(double** A, double* b, int n, int verbose) {
    if (verbose) printf("Tentando reordenar matriz para dominância diagonal\n");
    int* linha_usada = (int*) calloc(n, sizeof(int));
    int* nova_ordem = (int*) malloc(n * sizeof(int));
    for (int col = 0; col < n; col++) {
        int melhor_linha = -1;
        double melhor_valor = -1.0;
        for (int lin = 0; lin < n; lin++) {
            if (!linha_usada[lin]) {
                double valor_abs = fabs(A[lin][col]);
                if (valor_abs > melhor_valor) {
                    melhor_valor = valor_abs;
                    melhor_linha = lin;
                }
            }
        }
        if (melhor_linha == -1) break;
        linha_usada[melhor_linha] = 1;
        nova_ordem[col] = melhor_linha;
    }
    double** A_temp = aloca_matriz(n);
    double* b_temp = (double*) malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        int linha_original = nova_ordem[i];
        for (int j = 0; j < n; j++) {
            A_temp[i][j] = A[linha_original][j];
        }
        b_temp[i] = b[linha_original];
    }
    for (int i = 0; i < n; i++) {
        memcpy(A[i], A_temp[i], n * sizeof(double));
        b[i] = b_temp[i];
    }
    libera_matriz(A_temp, n);
    free(b_temp);
    free(linha_usada);
    free(nova_ordem);
    if (verbose) printf("Reordenamento concluído.\n");
}

//função de verificação
void verificar_solucao(double** A, double* b, double* x, int n) {
    double* Ax = (double*) calloc(n, sizeof(double));
    double max_erro = 0.0;
    double norma_b = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Ax[i] += A[i][j] * x[j];
        double erro = fabs(Ax[i] - b[i]);
        if (erro > max_erro) max_erro = erro;
        if (fabs(b[i]) > norma_b) norma_b = fabs(b[i]);
    }
    double erro_relativo = max_erro / (norma_b + 1e-12);
    printf("Erro relativo = %.3e\n", erro_relativo);
    free(Ax);
}

int jacobi_principal(double** A, double* b, double* x, int n, int verbose, int* out_iters) {
    double* novo_x = (double*) calloc(n, sizeof(double));
    double* residuo = (double*) calloc(n, sizeof(double));
    for (int iter = 0; iter < LIMITE_ITERACOES; iter++) {
        for (int i = 0; i < n; i++) {
            double soma = 0.0;
            for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
            if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo);}
            novo_x[i] = (b[i] - soma) / A[i][i];
        }
        double max_residuo = 0.0;
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            residuo[i] = b[i];
            for (int j = 0; j < n; j++) residuo[i] -= A[i][j] * novo_x[j];
            double diff = fabs(novo_x[i] - x[i]);
            if (diff > max_diff) max_diff = diff;
            if (fabs(residuo[i]) > max_residuo) max_residuo = fabs(residuo[i]);
        }
        memcpy(x, novo_x, n * sizeof(double));
        if (max_diff < LIMITE_ERRO && max_residuo < LIMITE_ERRO) {
            if (out_iters) *out_iters = iter + 1;
            free(novo_x); free(residuo); return 1;
        }
        if (verbose && (iter < 10 || iter % 1000 == 0)) {
            printf("Iteracao %d: Erro=%.3e  Residuo=%.3e\n", iter, max_diff, max_residuo);
        }
    }
    if (out_iters) *out_iters = LIMITE_ITERACOES;
    free(novo_x); free(residuo);
    return 0;
}

//aplica precondicionamento (se necessário) e chama jacobi_principal
int jacobi(double** A_in, double* b_in, double* x, int n, int verbose, int* out_iters) {
    double** A = copiar_matriz(A_in, n);
    double* b = copiar_vetor(b_in, n);
    int convergiu = 0;

    //tenta reordenação para dominância (heurística) - não obrigatória, só tentativa
    if (!eh_diagonalmente_dominante(A, n)) {
        if (verbose) printf("Matriz não é diagonalmente dominante: tentando reordenar e precondicionar.\n");
        reordenar_para_dominancia(A, b, n, verbose);
    }

    //aplicar pré-condicionamento diagonal (divide linha i por A[i][i])
    if (!precondiciona_diagonal(A, b, n, verbose)) {
        if (verbose) printf("Precondicionamento falhou, tentando usar Jacobi sem precondicionamento.\n");
        libera_matriz(A, n);
        free(b);
        return jacobi_principal(A_in, b_in, x, n, verbose, out_iters);
    }

    convergiu = jacobi_principal(A, b, x, n, verbose, out_iters);

    libera_matriz(A, n);
    free(b);
    return convergiu;
}

//Jacobi paralelo
int jacobi_omp(double** A_in, double* b_in, double* x, int n, int verbose, int num_threads, const char* sched_type, int chunk, int* out_iters) {
    double** A = copiar_matriz(A_in, n);
    double* b = copiar_vetor(b_in, n);
    int convergiu = 0;

    if (!eh_diagonalmente_dominante(A, n)) {
        if (verbose) printf("Matriz não diagonalmente dominante: reordenando e precondicionando...\n");
        reordenar_para_dominancia(A, b, n, verbose);
    }

    if (!precondiciona_diagonal(A, b, n, verbose)) {
        if (verbose) printf("Precondicionamento falhou; tentando sem precond.\n");
        libera_matriz(A, n);
        free(b);
        // usa versão paralela sem precondicionamento mas ainda tenta pôr diagonal não-nula
        A = copiar_matriz(A_in, n);
        b = copiar_vetor(b_in, n);
        if (!ensure_nonzero_diagonal(A, b, n)) {
            if (verbose) printf("Nao foi possivel garantir diagonal nao-zero.\n");
            libera_matriz(A, n); free(b);
            return 0;
        }
    }

    double* novo_x = (double*) calloc(n, sizeof(double));
    double* residuo = (double*) calloc(n, sizeof(double));

    omp_set_num_threads(num_threads);

    // Predefinir schedule no openmp
    for (int iter = 0; iter < LIMITE_ITERACOES; iter++) {
        // atualização do novo x
        if (sched_type == NULL || strcmp(sched_type, "static_simple") == 0) {
            #pragma omp parallel for
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else if (strcmp(sched_type, "static_chunk") == 0) {
            #pragma omp parallel for schedule(static, 1)
            for (int i = 0; i < n; i++) { 
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else if (strcmp(sched_type, "static") == 0) {
            #pragma omp parallel for schedule(static,1)
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else if (strcmp(sched_type, "staticN") == 0) {
            #pragma omp parallel for schedule(static, 25)
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else if (strcmp(sched_type, "dynamic25") == 0) {
            #pragma omp parallel for schedule(dynamic,25)
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else if (strcmp(sched_type, "guided25") == 0) {
            #pragma omp parallel for schedule(guided,25)
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        } else {
            #pragma omp parallel for
            for (int i = 0; i < n; i++) {
                double soma = 0.0;
                for (int j = 0; j < n; j++) if (j != i) soma += A[i][j] * x[j];
                if (fabs(A[i][i]) < 1e-14) { free(novo_x); free(residuo); libera_matriz(A, n); free(b);}
                novo_x[i] = (b[i] - soma) / A[i][i];
            }
        }

        double max_residuo = 0.0;
        double max_diff = 0.0;
        // calcular diffs (paralelo)
        #pragma omp parallel for reduction(max:max_diff,max_residuo)
        for (int i = 0; i < n; i++) {
            residuo[i] = b[i];
            double rowdot = 0.0;
            for (int j = 0; j < n; j++) rowdot += A[i][j] * novo_x[j];
            residuo[i] -= rowdot;
            double diff = fabs(novo_x[i] - x[i]);
            if (diff > max_diff) max_diff = diff;
            if (fabs(residuo[i]) > max_residuo) max_residuo = fabs(residuo[i]);
        }

        //copia novo x
        #pragma omp parallel for
        for (int i = 0; i < n; i++) x[i] = novo_x[i];

        if (max_diff < LIMITE_ERRO && max_residuo < LIMITE_ERRO) {
            if (out_iters) *out_iters = iter + 1;
            convergiu = 1;
            break;
        }
    }

    if (!convergiu && out_iters) *out_iters = LIMITE_ITERACOES;

    free(novo_x); free(residuo);
    libera_matriz(A, n);
    free(b);
    return convergiu;
}

int executa_testes(double** A_orig, double* b_orig, int n, Config cfg, int run_idx,FILE* csv, int verbose) {
    double** A = copiar_matriz(A_orig, n);
    double* b = copiar_vetor(b_orig, n);
    double* x = (double*) calloc(n, sizeof(double));

    struct timespec t0, t1;
    int iters = 0;
    int convergiu = 0;
    double elapsed = 0.0;

    if (strcmp(cfg.name, "sequencial") == 0) {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        convergiu = jacobi(A, b, x, n, verbose ? 1 : 0, &iters);
        clock_gettime(CLOCK_MONOTONIC, &t1);
    } else {
        const char* sched_code = NULL;
        if (strcmp(cfg.schedule, "static_simple") == 0) sched_code = "static_simple";
        else if (strcmp(cfg.schedule, "static") == 0 && cfg.chunk > 0) sched_code = "staticN";
        else if (strcmp(cfg.schedule, "static") == 0) sched_code = "static";
        else if (strcmp(cfg.schedule, "dynamic") == 0) sched_code = "dynamic25";
        else if (strcmp(cfg.schedule, "guided") == 0) sched_code = "guided25";
        else sched_code = "static_simple";

        clock_gettime(CLOCK_MONOTONIC, &t0);
        convergiu = jacobi_omp(A, b, x, n, verbose ? 1 : 0,
                               cfg.threads, sched_code, cfg.chunk, &iters);
        clock_gettime(CLOCK_MONOTONIC, &t1);
    }

    elapsed = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) * 1e-9;

    //salva resultado e valida solução
    if (convergiu) {
        char output[256];
        snprintf(output, sizeof(output), "testArchives/resultado%d.dat", n);
        FILE* out = fopen(output, "w");
        if (out) {
            for (int i = 0; i < n; i++) fprintf(out, "%.12lf\n", x[i]);
            fclose(out);
        }

        //valida solução Ax=b
        if (verbose) {
            printf("Validacao da solucao (%s, threads=%d, schedule=%s, rodada=%d):\n",
                   cfg.name, cfg.threads, cfg.schedule ? cfg.schedule : "none", run_idx+1);
        }
        verificar_solucao(A_orig, b_orig, x, n);
    }

    // escrevee linha CSV
    fprintf(csv, "%s,%d,%s,%d,%d,%.9f,%d,%d\n",
            cfg.name, cfg.threads, cfg.schedule ? cfg.schedule : "none",
            cfg.chunk, run_idx + 1, elapsed, iters, convergiu);

    fflush(csv);

    libera_matriz(A, n);
    free(b);
    free(x);
    return convergiu;
}


int main(int argc, char* argv[]) {
    int verbose = 0;
    if (argc >= 3) verbose = atoi(argv[2]);

    FILE* f = fopen(argv[1], "r");
    if (!f) { perror("Erro ao abrir arquivo"); return 1; }

    // encontra N pela primeira linha
    char first_line[65536];
    if (!fgets(first_line, sizeof(first_line), f)) { fclose(f); printf("Erro leitura\n"); return 1; }
    int n = 0;
    char* token = strtok(first_line, " \t\n");
    while (token != NULL) {
        char* endptr; strtod(token, &endptr);
        if (*endptr != '\0') { printf("Erro: formato inválido\n"); fclose(f); return 1; }
        n++; token = strtok(NULL, " \t\n");
    }
    rewind(f);
    if (verbose) printf("Sistema %dx%d encontrado\n", n, n);

    double** A = aloca_matriz(n);
    double* b = (double*) malloc(n * sizeof(double));
    double* x_init = (double*) calloc(n, sizeof(double));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(f, "%lf", &A[i][j]) != 1) { printf("Erro na leitura da matriz\n"); fclose(f); libera_matriz(A,n); free(b); free(x_init); return 1; }
        }
    }
    for (int j = 0; j < n; j++) {
        if (fscanf(f, "%lf", &b[j]) != 1) { printf("Erro na leitura do vetor\n"); fclose(f); libera_matriz(A,n); free(b); free(x_init); return 1; }
    }
    fclose(f);

    // matriz A e vetor b originais
    double** A_orig = copiar_matriz(A, n);
    double* b_orig = copiar_vetor(b, n);

    // combinações para teste
    Config configs[13] = {
        {"sequencial", 1, NULL, 0},
        {"openmp", 2, "static_simple", 0},
        {"openmp", 4, "static_simple", 0},
        {"openmp", 8, "static_simple", 0},
        {"openmp", 2, "static", 25},
        {"openmp", 4, "static", 25},
        {"openmp", 8, "static", 25},
        {"openmp", 2, "dynamic", 25},
        {"openmp", 4, "dynamic", 25},
        {"openmp", 8, "dynamic", 25},
        {"openmp", 2, "guided", 25},
        {"openmp", 4, "guided", 25},
        {"openmp", 8, "guided", 25}
    };

    char csvname[128];
    snprintf(csvname, sizeof(csvname), "benchmarks/resultados_%d.csv", n);
    FILE* csv = fopen(csvname, "w");
    if (!csv) { perror("Erro criando CSV"); return 1; }
    fprintf(csv, "config,threads,schedule,chunk,run,time_seconds,iterations,converged\n");

    // executa 3 rodadas de medicao
    for (int i = 0; i < 13; i++) {
        for (int r = 0; r < 3; r++) {
            if (verbose) printf("Combinacao %d/%d (threads=%d schedule=%s chunk=%d) rodada %d/3\n",
                               i+1, 13, configs[i].threads,
                               configs[i].schedule ? configs[i].schedule : "none",
                               configs[i].chunk, r+1);
            executa_testes(A_orig, b_orig, n, configs[i], r, csv, verbose);
        }
    }

    fclose(csv);
    printf("Benchmark concluido!\n");

    libera_matriz(A_orig, n);
    free(b_orig);
    libera_matriz(A, n);
    free(b);
    free(x_init);
    return 0;
}
