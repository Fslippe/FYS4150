#include <time.h>

int main ()
{

  // ...

  // Start measuring time
  clock_t t1 = clock();

  //
  // The code you want to perform timing on
  //

  // Stop measuring time
  clock_t t2 = clock();

  // Calculate the elapsed time.
  double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

  // ...
}
