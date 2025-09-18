# Projeto 01 de Programação Concorrente e Paralela
Implementação do método de Jacobi com openMP paralela e sequencialmente em arquivos fonte distintos, que resolvam sistemas lineares com solução factível, cuja saída deve retornar numa única linha a solução do sistema com quatro casas decimais, posteriormente será efetuada a comparação entre as versões.

## Funcionalidades
- Geração de sistemas lineares com **diagonal dominante** (garantia de convergência do Jacobi).  
- Implementações do Método de Jacobi: **sequencial** e **paralela**.    
- Comparação de desempenho entre as versões.  
- Verbosidade ajustável para depuração.

## Exemplo

### 1. Gera um sistema linear aleatório com diagonal dominante
```bash
python gera_sistema.py N
```
&emsp;Onde N é o tamanho do sistema (ex.: python gera_sistema.py 1000).

### 2. Compilar a versão paralela com OpenMP
```bash
gcc -fopenmp main.c -o main -lm
```
&emsp;Certifique-se de que o compilador gcc está instalado e com suporte a OpenMP.

### 3. Executar o programa
```bash
./main.exe testArchives/linear1000.dat 2
```
&emsp;O **primeiro argumento** é o caminho para o arquivo de entrada (`linearN.dat`).  

&emsp;O **segundo argumento** define o nível de **verbosidade**:  

&emsp;&emsp; `0` → execução silenciosa (apenas resultado final).  
&emsp;&emsp; `1` → modo detalhado para depuração.