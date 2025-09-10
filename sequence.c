#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define size 10
#define size2 size + 1

void initVet(float vet[size])
{
    for (int i = 0; i < size; i++)
    {
        vet[i] = 0;
    }
}

void initMatrix(float matrix[size2][size])
{
    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = 0;
        }
    }
}

void copyMatrix(float matrix[size2][size])
{

    FILE *f = fopen("/home/allan/Documentos/Programação/C/projetoPCP/linear10.dat", "r");
    if (f == NULL)
    {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < size2; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fscanf(f, "%f", &matrix[i][j]);
        }
    }
    fclose(f);
    return;
}

void calculateResultX(float matrix[size2][size], float result[size])
{
    int verify = 0;
    float buffer[size];
    printf("oláaaa");
    initVet(buffer);

    do
    {

        for (int i = 0; i < size; i++)
        {
            float sum = 0;
            for (int j = 0; j < size; j++)
                if (i != j)
                    sum += (matrix[i][j] * result[j]);
            buffer[i] = (matrix[size][i] - sum) / (matrix[i][i]);
        }

        for (int i = 0; i < size; i++)
            if (fabs((result[i] - buffer[i])) > pow(10, -5))
            {
                verify = 1;
                break;
            }
        for (int i = 0; i < size; i++)
            result[i] = buffer[i];

    } while (verify == 1);
    printf("oláaaa2");

    return;
}

int main()
{
    printf("oi4");

    float matrix[size2][size];
    float result[size];

    initMatrix(matrix);
    initVet(result);
    printf("oi1");
    copyMatrix(matrix);
    printf("oi2");
    calculateResultX(matrix, result);
    printf("oi3");

    for (int i = 0; i < size; i++)
        printf("X%d = %.4f\n", i + 1, result[i]);

    return 0;
}