#include <iostream>
#include <cstdlib>  // for rand() and srand()
#include <cmath>   // for sqrt, log, log10 etc

using namespace std;


/* --------------------------------------------------------------------- */
double ran01()
{
       // [0,1) random number using rand()
       // initialize by calls to srand()

       double x, xmax;

       xmax=double(RAND_MAX);
       x=double(rand());

       return x/xmax;

}

/* --------------------------------------------------------------------- */
double randg(const double sigma)
{
         // This function generates Gaussian random numbers in pairs
         // (e.g., y1, y2). A function that returns a single random number
         // would only calculate y1 or y2.
         // Algorithm by Dr. Everett (Skip) Carter, Jr.
         // here ran01() is the routine to obtain a random number uniformly
         // distributed in [0,1).

         double x1, x2, w, y1, y2;
 
         do {
                 x1 = 2.0 * ran01() - 1.0;
                 x2 = 2.0 * ran01() - 1.0;
                 w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = std::sqrt( (-2.0 * std::log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;

         return sigma*y1;
}


/* --------------------------------------------------------------------- */
void testRandNum()
{
  cout << "[0,1) random number ...\n";
  for (int nCount=0; nCount < 100; ++nCount)
    {
        cout << ran01() << "\t";
 
        // If we've printed 5 numbers, start a new column
        if ((nCount+1) % 5 == 0)
            cout << endl;
    }

  cout << "Gaussian random number ...\n";
  for (int nCount=0; nCount < 100; ++nCount)
    {
        cout << randg(0.5) << "\t";
 
        // If we've printed 5 numbers, start a new column
        if ((nCount+1) % 5 == 0)
            cout << endl;
    }
}


/* --------------------------------------------------------------------- */
void sort(int n, double arr[])
{
   /*Sorts an array arr[1..n] into ascending numerical order using the
     Quicksort algorithm. n is input; arr is replaced on output by its
     sorted rearrangement. Quicksort sorts by employing a divide and conquer
     strategy to divide a list into two sub-lists. The steps are:
    
     1.) Pick an element, called a pivot, from the list.
     2.) Reorder the list so that all elements with values less than the
         pivot come before the pivot, while all elements with values greater
         than the pivot come after it (equal values can go either way).
         After this partitioning, the pivot is in its final position.
         This is called the partition operation.
     3.) Recursively sort the sub-list of lesser elements and the sub-list
         of greater elements.
    
     n is input; arr is replaced on output by its sorted rearrangement.
     Parameters: M is the size of subarrays sorted by straight insertion;
                 NSTACK is maximum size of the required auxiliary storage;
                 istack stores the index ranges of the unprocessed subarrays. */

   #define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
   #define M 7
   #define NSTACK 50

   int i,ir=n-1,j,k,l=0;    // l & ir are the left and right indexes
   int jstack=0,istack[NSTACK];       // of the sub-array
   double a,temp;

   for (;;) {                   // Insertion sort when subarray small enough.
       if (ir-l < M) {
          for (j=l+1;j<=ir;j++) {
              a=arr[j];
              for (i=j-1;i>=l;i--) {
                  if (arr[i] <= a) break;
                  arr[i+1]=arr[i];
              }
              arr[i+1]=a;
          }
          if (jstack == 0) break;
          ir=istack[jstack--]; // Pop stack and begin a new round of partitioning.
          l=istack[jstack--];
       } else {                // Choose median of left, center, and right elements
                               // as partitioning element a. Also rearrange so that
                               // a[l] < a[l+1] < a[ir].
          k=int((l+ir)/2);
          SWAP(arr[k],arr[l+1]);
          if (arr[l] > arr[ir]) {
             SWAP(arr[l],arr[ir]);
          }
          if (arr[l+1] > arr[ir]) {
             SWAP(arr[l+1],arr[ir]);
          }
          if (arr[l] > arr[l+1]) {
             SWAP(arr[l],arr[l+1]);
          }
          i=l+1;                  // Initialize pointers for partitioning.
          j=ir;
          a=arr[l+1];             // Partitioning element.
          for (;;) {              // Beginning of innermost loop.
              do i++; while (arr[i] < a);  // Scan up to find element > a.
              do j--; while (arr[j] > a);  // Scan down to find element < a.
              if (j < i) break;            // Pointers crossed. Partitioning complete.
              SWAP(arr[i],arr[j]);         // Exchange elements.
           }                      // End of innermost loop.
           arr[l+1]=arr[j];       // Insert partitioning element.
           arr[j]=a;
           jstack += 2;
                                  // Push pointers to larger subarray on stack,
                                  // process smaller subarray immediately.
           if (jstack > NSTACK) {
              cout << "Error: NSTACK too small in sort.\n";
              exit(1);
           }
           if (ir-i+1 >= j-l) {
              istack[jstack]=ir;
              istack[jstack-1]=i;
              ir=j-1;
           } else {
              istack[jstack]=j-1;
              istack[jstack-1]=l;
              l=i;
           }
       }
    }
}
