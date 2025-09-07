#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define size 10

void initVet(float *vet[size])
{
    for (int i = 0; i < size; i++)
    {
        *vet[i] = 0;
    }
}

void copyMatrix(float *matrix[size + 1][size])
{

    FILE *f = fopen("testArchives/linear10.dat", "r");
    if (f == NULL)
    {
        perror("Error opening file");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < size + 1; i++)
    {
        for (int j = 0; j < size; j++)
        {
            fscanf(f, "%f", &matrix[i][j]);
        }
    }
    fclose(f);
}

void calculateResultX(float *matrix[size + 1][size], float *result[size])
{
    float sum = 0;
    int count = 0;
    float buffer[size];
    int verify = 1;

    initVet(buffer);

    while (count < size)
    {
        for (int i = 0; i < size; i++)
            if (count != i)
                sum += *matrix[count][i] * (*result[i]);

        buffer[count] = (*matrix[size + 1][count] - sum) / (*matrix[count][count]);
        count++;
    }

    for (int i = 0; i < size; i++)
        if (fabs(*result[i] - buffer[i] < pow(10, -5)))
            verify = 0;

    for (int i = 0; i < size; i++)
        *result[i] = buffer[i];

    return verify ? calculateResultX(matrix, result) : 0;
}

int main()
{
    float matrix0[size + 1][size];
    float result[size];

    copyMatrix(matrix0);
    initVet(result);
    calculateResultX(matrix0, result);

    return EXIT_SUCCESS;
}