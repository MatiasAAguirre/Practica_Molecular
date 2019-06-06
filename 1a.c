#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float set_box(float *r, int N, float L);

int main(int argc, char const *argv[]) {
  int N=8, i;
  float L=4;
  float *r;

  r = malloc(3*N*sizeof(float));

  for(i=0;i<3*N;i++) {
    *(r+i) = 0;
  }

  set_box(r,N,L);

  for(i=1; i<3*N+1;i++) {
    printf("%f ", *(r+i-1));
    if(i%3==0) {
      printf("\n");
    }
  }

  free(r);

  return 0;
}

float set_box(float *r, int N, float L) {
  int n = cbrt(N), i=0, x, y, z;
  float dL=L/n;

  for(x=0; x<n; x++) {
    for(y=0; y<n; y++) {
      for(z=0; z<n; z++) {

          *(r+3*i) = dL*(x+0.5);
          *(r+3*i+1) = dL*(y+0.5);
          *(r+3*i+2) = dL*(z+0.5);

          i++;

      }
    }
  }
  return 0;
}
