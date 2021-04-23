#include <stdio.h>
#include <stdlib.h>

#define DIM 15

/** 
 * @brief Copia dell'elemento di indice indOrigine dell'array Orgine, nell'elemento di indice indDestinazione dell'array Destinazione
 * @param Origine Array dal quale si vuole copiare
 * @param Destinazione Array sul quale si vuole copiare
 * @param indOrigine indice dell'elemento da copiare dell'array d'origine
 * @param indDestinazione posizione sulla quale copiare nell'array di destinazione
 * 
*/
void CopyFrom(int *Origine, int *Destinazione, int indOrigine, int indDestinazione)
{

    Destinazione[indDestinazione] = Origine[indOrigine];
}

void merge(int *A, int low, int mid, int high)
{
    int i;
    int j;
    int k;                  // serve come indice di controllo
    int n1 = mid - low + 1; // trovo la dimensione n1 (primo array)
    int n2 = high - mid;    // trovo la dimensione n2 (secondo array)

    // Let L[1,...,n1], R[1,...,n2] be new array
    int L[n1], R[n2];

    /* Effettuo la copia degli elementi di A nei due nuovi array esterni */
    for (i = 0; i < n1; i++)
    {
        L[i] = A[low + i]; // il -1 dello pseudocodice si toglie...
    }
    for (j = 0; j < n2; j++)
    {
        R[j] = A[mid + 1 + j];
    }
    i = 0; // indice iniziale del primo sottoarray di sinistra
    j = 0; // indice iniziale del secondo sottoarray di destra
    // k è l'indice iniziale del sottoarray unito
    for (k = low; k <= high; k++) // inizializzo k che va da low ad high (k<=high)
    {
        if (i < n1) // Sei arrivato a fine corsa con l'indice i?
        {
            if (j < n2) // Sei arrivato a fine corsa con l'indice j?
            {
                // Se non sono arrivato a fine corsa con nessuno dei due, allora chi decide che cosa ci va
                // nella posizione k-esima è il confronto tra L[i] ed R[j]
                if (L[i] <= R[j])
                {
                    CopyFrom(L, A, i, k); //A[k] = L[i]; // Copio L[i] in A[k]
                    i++;
                }
                else // Se (L[i] >= R[j] )
                {
                    CopyFrom(R, A, j, k); //A[k] = R[j];
                    j++;
                }
            }
            /* Altrimenti copio gli elementi rimanenti */
            else
            {
                CopyFrom(L, A, i, k); //A[k] = L[i];
                i++;
            }
        }
        else
        {
            CopyFrom(R, A, j, k); // A[k] = R[j];
            j++;
        }
    } // fine for k
} // fine Merge

void mergeSort(int *A, int low, int high)
{
    int mid;
    if (low < high)
    {
        mid = (low + high) / 2;
        mergeSort(A, low, mid);
        mergeSort(A, mid + 1, high);
        merge(A, low, mid, high);
    }
}

int main(void)
{
    srand(11);
    /* int A[DIM] = {
        52,
        39,
        35,
        68,
        32,
        21,
        10,
        4,
        1,
        24,
        12,
        17,
        98,
        7,
        19}; */
    int i;
    int A[DIM];
    for (i = 0; i < DIM; i++)
        A[i] = +rand() % 100; // genero numeri pseudo-casuali da 1 a 100
    printf("Array DISORDINATO:\n");
    for (i = 0; i < DIM; i++)
        printf("%d ", A[i]);
    printf("\n");
    /* Chiamata a mergeSort() */
    mergeSort(A, 0, DIM - 1);
    printf("Array ORDINATO:\n");
    for (i = 0; i < DIM; i++)
        printf("%d ", A[i]);
    printf("\n");
    return 0;
}
