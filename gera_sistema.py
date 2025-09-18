import numpy as np
import sys

def gerar_sistema_diagonalmente_dominante(N):
    A = np.random.rand(N, N)
    for i in range(N):
        soma_linha = np.sum(np.abs(A[i, :])) - np.abs(A[i, i])
        A[i, i] = soma_linha + 1.0
    x_true = np.arange(1, N + 1, dtype=float)
    
    b = A.dot(x_true)
    
    return A, b, x_true

def salvar(A, b, filename):
    with open(filename, "w") as f:
        for i in range(A.shape[0]):
            linha = " ".join([f"{val:.12f}" for val in A[i]])
            f.write(linha + "\n")
        linha_b = " ".join([f"{val:.12f}" for val in b])
        f.write(linha_b + "\n")

if __name__ == "__main__":
    N = int(sys.argv[1])
    
    print(f"Gerando sistema {N}x{N}")
    A, b, x_true = gerar_sistema_diagonalmente_dominante(N)
    
    filename = f"testArchives/linear{N}.dat"
    salvar(A, b, filename)
    
    print(f"salvo em {filename}")