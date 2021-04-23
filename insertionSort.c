#include <stdio.h>

#define DIM 15

/* Stampa Array */
void stampa(int A[], int dim);

/* // Insertion Sort
void InsertionSort(int A[], int low, int high);
*/

/*void insertionSort(int arr[], int n)  
{  
    int i, key, j;  
    for (i = 1; i < n; i++) 
    {  
        key = arr[i];  
        j = i - 1;  
  
         Move elements of arr[0..i-1], that are  greater than key, to one position ahead  of their current position 
        while (j >= 0 && arr[j] > key) 
        {  
            arr[j + 1] = arr[j];  
            j = j - 1;  
        }  
        arr[j + 1] = key;  
    }  
}  
*/

void insertionSort(int *A, const int low, const int high) // implementato
{
    int i, j, key = 0;
    for (j = low; j <= high; j++)
    {
        key = A[j];
        i = j - 1;
        while ((i >= low) && (A[i] > key))
        {
            A[i + 1] = A[i];
            i--;
        }
        A[i + 1] = key;
    }
}

/* Main */
int main(void)
{
    // Provo l'esecuzione di Insertion Sort
    // int i = 0;
    int A[DIM] = {3, 1, 9, 6, 0, 15, 11, 18, 10, 8, 25, 3, 2, 20, -1};

    fprintf(stdout, "In tutto ci sono %d elementi\n", DIM);
    fprintf(stdout, "Vettore non ordinato:\n");
    // Stampa vettore non ordinato
    stampa(A, DIM);
    // Chiamata ad IS
    insertionSort(A, 0, DIM);
    // InsertionSort(A, 0, DIM);
    fprintf(stdout, "Vettore ordinato:\n");
    // Stampa vettore ordinato
    stampa(A, DIM);

    return 0;
}

void stampa(int A[], int dim)
{
    int i;
    for (i = 0; i < dim; i++)
    {
        printf("%d ", A[i]);
    }
    printf("\n");
}

/*
void InsertionSort(int A[], int low, int high)
{
    int i, j, key = 0;
    for (j = 1; j < (high - low + 1); j++)
    {
        key = A[j];
        i = j - 1;
        while ((i > -1) && (A[i] > key))
        {
            A[i + 1] = A[i];
            i--;
        }
        A[i + 1] = key;
    }
}
*/