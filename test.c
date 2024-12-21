
#include <stddef.h>
#include <stdio.h>
#include <math.h>

#define N 1000
#define PI 3.14159265358979323844

int main() {
  for (ptrdiff_t i = 0; i < N; i++) {
    double f = -2.0 * cos(2.0 * PI * (double) i / (double) N);
    printf("%f\n", f);
  }
  return 0;
}
