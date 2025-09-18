/* CONDIÇÕES DE TESTE 

Os tempos de execução devem ser medidos para cada condiçãoao de teste a seguir, sendo que para cada uma delas deve ser apresentada a média de três execuções:

    1. Versão sequencial pura;
    2. Versão openMP com 2 threads, schedule estático simples
    3. Versão openMP com 4 threads, schedule estático simples
    4. Versão openMP com 8 threads, schedule estático simples
    5. Versão openMP com 2 threads, schedule estático, chunk de 25
    6. Versão openMP com 4 threads, schedule estático, chunk de 25
    7. Versão openMP com 8 threads, schedule estático, chunk de 25
    8. Versão openMP com 2 threads, schedule dinâmico, chunk de 25
    9. Versão openMP com 4 threads, schedule dinâmico, chunk de 25
    10. Versão openMP com 8 threads, schedule dinâmico, chunk de 25
    11. Versão openMP com 2 threads, schedule guided, chunk de 25
    12. Versão openMP com 4 threads, schedule guided, chunk de 25
    13. Versão openMP com 8 threads, schedule guided, chunk de 25
pra rodar precisa colocar o -fopenmp:
    gcc -o paralelo paralelo.c -fopenmp */

// troca OMP_NUM_THREADS e OMP_SCHEDULE antes de rodar.
// --> OMP_NUM_THREADS=4 OMP_SCHEDULE="static,25" ./jacobi_omp_to_file entrada.txt

#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define MAXN 2000
#define TOL 1e-5
#define MAX_ITERS 100000
#define IDX(i,j,n) ((size_t)(i)*(size_t)(n) + (size_t)(j))

static void die(const char *m){ fprintf(stderr,"%s\n",m); exit(EXIT_FAILURE); }

static void read_input(const char *path, int *out_n, double **out_A, double **out_b){
    FILE *fin = fopen(path,"r");
    if (!fin) die("Nâo foi possível abrir o arquivo de entrada");
    size_t bufsize = 1<<20;
    char *line = (char*)malloc(bufsize);
    if (!line) die("Falha ao alocar buffer");
    if (!fgets(line,(int)bufsize,fin)) die("Arquivo vazio");

    int n=0;
    { char *tmp=strdup(line); if(!tmp) die("Falha strdup");
      for (char *tok=strtok(tmp," \t\r\n"); tok; tok=strtok(NULL," \t\r\n")) n++;
      free(tmp); }
    if (n<=0 || n>MAXN) die("Numero de colunas invalido (1..2000)");

    double *A = malloc((size_t)n*(size_t)n*sizeof(double));
    double *b = malloc((size_t)n*sizeof(double));
    if (!A||!b) die("Falha ao alocar A/b");

    // linha 0 de A
    { int j=0;
      for (char *tok=strtok(line," \t\r\n"); tok && j<n; tok=strtok(NULL," \t\r\n"), j++)
        A[IDX(0,j,n)] = strtod(tok,NULL);
      if (j!=n) die("Primeira linha com numero de colunas != de n"); }

    // linhas 1..n-1 de A
    for (int i=1;i<n;i++){
        if (!fgets(line,(int)bufsize,fin)) die("Linhas insuficientes para A");
        int j=0;
        for (char *tok=strtok(line," \t\r\n"); tok && j<n; tok=strtok(NULL," \t\r\n"), j++)
            A[IDX(i,j,n)] = strtod(tok,NULL);
        if (j!=n) die("Número de colunas != n");
    }

    // última linha: b
    if (!fgets(line,(int)bufsize,fin)) die("Faltou a linha com coeficientes do vetor b");
    { int j=0;
      for (char *tok=strtok(line," \t\r\n"); tok && j<n; tok=strtok(NULL," \t\r\n"), j++)
        b[j] = strtod(tok,NULL);
      if (j!=n) die("Vetor b != n valores"); }

    fclose(fin);
    free(line);
    *out_n = n; *out_A = A; *out_b = b;
}

int main(int argc, char **argv){
    if (argc < 2) die("Chamada errada --> ./paralelo <arquivo_entrada>");

    int n; double *A,*b;
    read_input(argv[1], &n, &A, &b);

    double *x    = calloc((size_t)n,sizeof(double));
    double *xnew = malloc((size_t)n*sizeof(double));
    double *invD = malloc((size_t)n*sizeof(double));
    if (!x||!xnew||!invD) die("Falha ao alocar matrizes");

    for (int i=0;i<n;i++){
        double dii = A[IDX(i,i,n)];
        if (fabs(dii) < 1e-15) die("Diagonal zero (Jacobi indefinido)");
        invD[i] = 1.0/dii;
    }

    double t0 = omp_get_wtime();
    int converged = 0;
    for (int iter=0; iter<MAX_ITERS; iter++){
        double maxdiff = 0.0;
        #pragma omp parallel for default(none) shared(A,b,x,xnew,invD,n) reduction(max:maxdiff) schedule(runtime)
        for (int i=0;i<n;i++){
            const size_t row=(size_t)i*(size_t)n;
            double sum=0.0;
            for (int j=0;j<n;j++){ if (j==i) continue; sum += A[row+(size_t)j]*x[j]; }
            double xi_new = (b[i]-sum)*invD[i];
            xnew[i]=xi_new;
            double d=fabs(xi_new - x[i]); if (d>maxdiff) maxdiff=d;
        }
        if (maxdiff < TOL){ converged=1; break; }
        double *tmp=x; x=xnew; xnew=tmp;
    }
    double t1 = omp_get_wtime();
    if (converged){ for (int i=0;i<n;i++) x[i]=xnew[i]; }

    // ---- escreve em arquivo .txt ----
    FILE *fout = fopen("saida.txt","w");
    if (!fout) die("Não foi possível abrir arquivo de saída");

    fprintf(fout,"Tempo: %.6f segundos\n",(t1-t0));
    for (int i=0;i<n;i++){
        if (i) fputc(' ', fout);
        double v=x[i]; if (fabs(v)<0.00005) v=0.0;
        fprintf(fout,"%.4f", v);
    }
    fputc('\n', fout);
    fclose(fout);

    free(A); free(b); free(x); free(xnew); free(invD);
    return 0;
}