/**
 * @brief Problem 1, Laboratory of Algorithms and Data Structures.
 * @author SCIAVICCO Guido (guido.sciavicco@unife.it)
 * @author STAN Ionel Eduard (ioneleduard.stan@unife.it)
 * @version Student
 */

// ----- INCLUSIONE LIBRERIE ----- //

// Standard input-output library (e.g., fprintf).
#include <stdio.h>
// Time library (e.g., time, clock()).
#include <time.h>
// Standard library (e.g., rand, srand).
#include <stdlib.h>
// Boolean library (e.g., bool).
#include <stdbool.h>
// String library (e.g., memcpy, strcmp)
#include <string.h>

// ----- fine INCLUSIONE LIBRERIE ----- //

// ----- STRUTTURE DATI AUSILIARIE ----- //

/**
 * @brief Enumeration data type for the output.
 */
// Struttura dati che decide se voglio il risultato sullo schermo oppure su file
typedef enum
{
    ONCONSOLE, // On console. = 0
    ONFILE     // On file. = 1
} outputEnumType;

/**
 * @brief Pair data structure.
 */
// Introduco un nuovo tipo: "pairType(tipoCoppia) per chiedere che ogni alg. di ordinamento restituisca questa coppia"
typedef struct
{
    clock_t time;  // Time needed to sort. Il tempo dell'algoritmo
    bool isSorted; // Flag representing if the algorithm performed correctly its task. (Risultato della funzione antagonista)
} pairType;

// ----- fine STRUTTURE DATI AUSILIARIE ----- //

// ----- VARIABILI GLOBALI ----- //

// IMPORTANT: queste costanti sono inserite solo per propositi didattici, devo trovare quelle giuste!
/* Ricordo che: il mio esperimento va da minSize a maxSize, saltando di granularità */
// Seed (important for reproducibility).
time_t SEED = 20;
// Minimum size of the array. Il mio esperimento va da minSize
const int minSize = 5; // 10
// Maximum size of the array. Ed arriva a maxSize
const int maxSize = 1200; //  500 oppure 2000
// Number of experiments. Fissata una dimensione, questa è testata numExperiments volte
const int numExperiments = 100; //  100
// Granularity of the experiment. Di quanto avanzo ogni volta
const int granularity = 5; // 10
// Maximum random integer allowed when generating random numbers.
const int maxRandInt = 1000000; // 1000000
// Thereshold parameter for the base case of HybridSort.
// IMPORTANT: this is the result of the first part of the experiment!
// threshold (scalino, livello) è la costante k che devo trovare dal primo esperimento
// questo punto di incrocio lo vedo dal confronto tra IS e MS e lo devo mettere qui per poter fare la 2° parte dell'esercizio
const int threshold = 235; // K=235
// Output type.
const outputEnumType outputType = ONCONSOLE;
// Output pointer (for printing).
FILE *outputPointer;

// ----- fine VARIABILI GLOBALI ----- //

// ----- FUNZIONI AUSILIARIE ----- //

/**
 * @brief Generate a collection of random numbers into an array A of size n.
 * @param A Array of random numbers.
 * @param n Size of the array.
 */
void generateRandomArray(int *A, const int n)
{
    // For each i in 0..n-1, generate a random number in the interval [1,maxRandInt].
    for (int i = 0; i < n; i++)
        A[i] = rand() % maxRandInt + 1;
}

/**
 * @brief Display the array A of size n.
 * @param A Array to be displayed.
 * @param n Size of the array.
 */
void displayArray(const int *A, const int n)
{
    // For each i in 0..n-1, print A[i].
    for (int i = 0; i < n; i++)
        printf("%d ", A[i]);
    printf("\n");
}

/** 
 * @brief Copia di UN elemento da un array Origine ad un array Destinazione
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

/**
 * @brief Copia GLI elementi del vettore B nel vettore A.
 * @param A Vettore di origine
 * @param B Vettore di destinazione
 * @param dim Dimensione uguale di entrambi i vettori
 */
void Copy(int *A, int *B, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
    {
        B[i] = A[i]; // copio elemento dell'origine nell'elemento destinazione
    }
    // Per correttezza, controllo se gli elementi di posto uguale sono di valore uguale
    // commento questa funzione, perchè incide negativamente sul tempo di esecuzione dell'esperimento
    /*
    for (i = 0; i < dim; i++)
    {
        if (A[i] == B[i])
            ;
        else
            break;
    }
    */
}

// ----- fine FUNZIONI AUSILIARIE ----- //

// ----- FUNZIONI ANTAGONISTE ----- //

/**
 * @brief Unit test: check if the input array A of size n is sorted.
 * @param A Array to be checked if sorted.
 * @param n Size of the array.
 * @return true if it is sorted; otherwise, false
 */
bool isSorted(const int *A, const int n)
{
    // For each i in 0..n-2, if the current element is greater than the next one,
    // then it is unsorted.
    for (int i = 0; i < n - 1; i++)
        if (A[i] > A[i + 1])
            return false;
    // Otherwise it is.
    return true;
}

// ----- fine FUNZIONI ANTAGONISTE ----- //

// ----- CORE FUNCTIONS ----- //

/**
 * @brief InsertionSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes time O(n^2), where n=high-low+1 is the size of the input array.
 */
void insertionSort(int *A, const int low, const int high) // implementato
{
    int i, j, key = 0;
    for (j = low; j <= high; j++) // for (j = low; j <high+1; j++)
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

/**
 * @brief Merge algorithm.
 * @param A Array to be merged.
 * @param low Left-end index of the array.
 * @param mid Mid index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n), where n=high-low+1 is the size of the input array.
 */
void merge(int *A, const int low, const int mid, const int high) //implementato
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
            /* Altrimenti copio gli elementi rimanenti da L*/
            else
            {
                CopyFrom(L, A, i, k); //A[k] = L[i];
                i++;
            }
        }
        /* Altrimenti copio gli elementi rimanenti da R*/
        else
        {
            CopyFrom(R, A, j, k); // A[k] = R[j];
            j++;
        }
    } // fine for k
} // fine Merge

/**
 * @brief MergeSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n*logn), where n=high-low+1 is the size of the input array.
 */
void mergeSort(int *A, const int low, const int high) // implementato
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

/**
 * @brief HybridSort algorithm.
 * @param A Array of random numbers to be sorted.
 * @param low Left-end index of the array.
 * @param high Right-end index of the array.
 * @property It takes O(n*logn), where n=high-low+1 is the size of the input array.
 */

void hybridSort(int *A, const int low, const int high) // implementato
{
    int mid;
    if ((high - low + 1) > threshold) // Per dimensioni maggiori di threshold chiamo MS
    {
        mid = (low + high) / 2;
        hybridSort(A, low, mid);
        hybridSort(A, mid + 1, high);
        merge(A, low, mid, high);
    }
    else // Per dimensioni minori di threshold chiamo IS
        insertionSort(A, low, high);
}

/**
 * @brief Polymorphic function that calls different sorting algorithms.
 * @param randomArray Array of random numbers to be sorted.
 * @param dim Dimension of the input array.
 * @param algo Sorting algorithm to be called. The possible values are: insertionSort,
 *             mergeSort and hybridSort.
 * @return Pair containing the total time needed to sort and the isSorted flag.
 */
/* Funzione "polimorfica": quando chiamerò un algoritmo di ordinamento, in realtà chiamerò la funzione sortArray()
dicendo quale algoritmo voglio, usando la stringa algo.
L'ordinamento è effettuato una volta sola, in base al tipo di parametro scelto (algo) */
pairType sortArray(const int *randomArray, const int dim, const char *algo)
{
    // Initiliazation of a pairType with values time = 0 and isSorted = true.
    pairType pair = {0, true};

    // Start and end time.
    clock_t startTime, endTime = 0;

    // Allocate memory for dim integers. Array che mi serve per copiarci la dimensione giusta di tutto l'array allocato
    int *sliceRandomArray = malloc(dim * sizeof(int));

    // Put every i-th element of randomArray into sliceRandomArray for each i in [0, ..., dim-1].
    for (int i = 0; i < dim; i++)
        sliceRandomArray[i] = randomArray[i];

    // Start the clock. Prendo il tempo per l'algoritmo scelto
    startTime = clock();

    // Use InsertionSort. Se la stringa inserita è insertionSort allora chiamo la funzione insertionSort
    if (strcmp(algo, "insertionSort") == 0)
        insertionSort(sliceRandomArray, 0, dim - 1);
    // Use MergeSort. Analogo discorso per mergeSort
    else if (strcmp(algo, "mergeSort") == 0)
        mergeSort(sliceRandomArray, 0, dim - 1);
    // Use HybridSort. Analogo discorso per hybridSort
    else if (strcmp(algo, "hybridSort") == 0)
        hybridSort(sliceRandomArray, 0, dim - 1);
    // Error
    else
    {
        fprintf(stderr, "ERROR: There is no such sorting algorithm called %s\n", algo);
        exit(1);
    }
    // Stop the clock. Prendo il tempo finale per l'algoritmo scelto
    endTime = clock();

    // Total time needed to sort.
    pair.time = endTime - startTime;
    // Have we sorted the instance? If not, then the flag is set to false; otherwise, it remains true.
    if (!isSorted(sliceRandomArray, dim))
        pair.isSorted = false;

    // Free sliceRandomArray.
    free(sliceRandomArray);

    /* Ritorno la coppia che viene costruita con:
       la differenza di tempi
       il booleano: è o no ordinato
    */
    return pair;
}

// ----- End CORE FUNCTIONS ----- //

// ----- MAIN FUNCTION ----- //

/**
 * @brief Main function.
 * @return Exit code 0.
 */
int main()
{
    // Initialize the random seed only once.
    srand(SEED);

    // Accumulated times.
    clock_t timeIS = 0; // InsertionSort
    clock_t timeMS = 0; // MergeSort
    clock_t timeHS = 0; // HybridSort

    // Flags saying either that the output of the execution is correct or not.
    bool isSortedIS = true; // InsertionSort
    bool isSortedMS = true; // MergeSort
    bool isSortedHS = true; // HybridSort

    // Creo la coppia per ognuno dei 3 algoritmi.
    pairType pairIS; // InsertionSort
    pairType pairMS; // MergeSort
    pairType pairHS; // HybridSort

    // Allocate an array of maxSize*sizeof(int) cells on the heap.
    // We use this array as a container.
    int *randomArray = malloc(maxSize * sizeof(int));
    // Array necessario per la funzione Copy()
    int *A = malloc(maxSize * sizeof(int));
    int *B = malloc(maxSize * sizeof(int));
    // Vedo se ho scelto di vedere su console o su file
    if (outputType == ONCONSOLE || outputType == ONFILE)
    {
        // On console.
        if (outputType == ONCONSOLE)
            outputPointer = stdout;
        // On file.
        else
        {
            // Open file.
            outputPointer = fopen("results.txt", "w");
            // Have we opened the file?
            if (outputPointer == NULL)
            {
                fprintf(stderr, "Error: The outputPointer has not been created\n");
                exit(1);
            }
        }
    }
    // Error
    else
    {
        fprintf(stderr, "Error: The outputType can be only ONCONSOLE or ONFILE\n");
        exit(1);
    }

    // // Print the header, only if it is on console.
    if (outputType == ONCONSOLE)
    {
        fprintf(outputPointer, "+-----------+-------------------------------+-------------------------------+-------------------------------+\n");
        fprintf(outputPointer, "| ######### | InsertionSort                 | MergeSort                     | HybridSort                    |\n");
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
        fprintf(outputPointer, "| Dimension | Time              | isSorted? | Time              | isSorted? | Time              | isSorted? |\n");
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");
    }

    // Going from minSize to maxSize with step equal to granularity.
    for (int dim = minSize; dim <= maxSize; dim += granularity)
    {
        // Reset the accumulated times from one experiment to another.
        timeIS = 0; // InsertionSort
        timeMS = 0; // MergeSort
        timeHS = 0; // HybridSort
        // Reset the isSorted flag for InsertionSort from one experiment to another.
        // We set this flag to true at the beginning, then if at least one time
        // the InsertionSort algorithm fails to sort the input such flag will be
        // set to false; otherwise, it remains true. Similarly for the other algorithms.
        isSortedIS = true; // InsertionSort
        isSortedMS = true; // MergeSort
        isSortedHS = true; // HybridSort

        // Repeat the experiment a numExperiments times for the fixed size (dim).
        for (int exper = 0; exper < numExperiments; exper++)
        {
            // Fill the array with (pseudo-) random numbers. That is, initialize only the
            // prefix of size dim (<= maxSize) with random numbers.
            generateRandomArray(randomArray, dim);
         /* Copy(randomArray, A, dim); // per mergeSort
            Copy(randomArray, B, dim); // per hybridSort */
            // InsertionSort.
            pairIS = sortArray(randomArray, dim, "insertionSort"); // chiamo IS
            timeIS += pairIS.time;                                 // calcolo il tempo di IS
            isSortedIS = pairIS.isSorted;                          // vedo se è o no ordinato

            // Seconda parte: devo chiamare sortArray con mergeSort e con hybridSort

            // MergeSort.
            // generateRandomArray(randomArray, dim);             // rigenero lo stesso array disordinato
            pairMS = sortArray(randomArray, dim, "mergeSort"); // chiamo MS
            timeMS += pairMS.time;                   // calcolo il tempo di MS
            isSortedMS = pairMS.isSorted;

            // HybridSort.
            //generateRandomArray(randomArray, dim);              // rigenero lo stesso array disordinato
            pairHS = sortArray(randomArray, dim, "hybridSort"); // chiamo HS
            timeHS += pairHS.time;                    // calcolo il tempo di HS
            isSortedHS = pairHS.isSorted;
        }
        // Printing the (sample mean as) result. Use TAB (\t) on file.
        // Stampo le medie e le scrivo o sulla console o sul file a seconda dei casi
        if (outputType == ONCONSOLE)
            fprintf(outputPointer, "| %9d | %17f | %9s | %17f | %9s | %17f | %9s |\n",
                    dim,
                    (double)timeIS / numExperiments, isSortedIS ? "true" : "false",  // InsertionSort
                    (double)timeMS / numExperiments, isSortedMS ? "true" : "false",  // MergeSort
                    (double)timeHS / numExperiments, isSortedHS ? "true" : "false"); // HybridSort
        else
            fprintf(outputPointer, "%9d\t %17f\t %9s\t %17f\t %9s\t %17f\t %9s\n",
                    dim,
                    (double)timeIS / numExperiments, isSortedIS ? "true" : "false",  // InsertionSort
                    (double)timeMS / numExperiments, isSortedMS ? "true" : "false",  // MergeSort
                    (double)timeHS / numExperiments, isSortedHS ? "true" : "false"); // HybridSort
    }

    // Print the ending part, only if it is on console.
    if (outputType == ONCONSOLE)
        fprintf(outputPointer, "+-----------+-------------------+-----------+-------------------+-----------+-------------------+-----------+\n");

    // Free the allocated memory.
    free(randomArray);
    free(A);
    free(B);

    // If the output is on file, we need to close the file.
    if (outputType == ONFILE)
        fclose(outputPointer);

    // We are done.
    return 0;
}

// ----- End MAIN FUNCTION ----- //

/* PER DEBUGGARE USARE COSTANTI BASSE:
time_t SEED = 20;
const int minSize = 1; // 10
const int maxSize = 6; //  500 oppure 2000
const int numExperiments = 5; //  100
const int granularity = 1; // 10
const int maxRandInt = 1000000;
const int threshold = 0; // K=235
const outputEnumType outputType = ONCONSOLE;
 */