#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

float set_box(float *r, int N, float L);
float set_v(float *v, int N, float T);
float Gaussiana(float nu, float sigma);
float fuerza(int N, int L, float *r, float *rij, float *Fij, float *Fr);

int main(int argc, char const *argv[]) {
  int N=8, i;
  float L=4, T=6;
  float *r, *v, *rij, *rij2, *Fij, *Vij, *Fr;

  srand(time(NULL));

  r = malloc(3*N*sizeof(float));
  v = malloc((3*N+3)*sizeof(float));
  Fr = malloc(3*N*sizeof(float));

  for(i=0; i<3*N; i++) {
    *(Fr + i) = 0.0;
  }

//############Tabla para interpolar###################//

  rij = malloc(100000*sizeof(float));
  rij2 = malloc(100000*sizeof(float));
  Fij = malloc(100000*sizeof(float));
  Vij = malloc(100000*sizeof(float));

  for(i=0; i<100000; i++) {
    *(rij + i) = 0.01 + (2.5/99999) * i;
  }

  for(i=0; i<100000; i++) {
    *(rij2 + i) = *(rij + i) * *(rij + i);
  }

  for(i=0; i<100000; i++) {
    *(Fij + i) = 24 * (2/ (*(rij+i) * pow(*(rij2 + i),6))  - 1/ (*(rij+i) * pow(*(rij2 + i),3)));
  }

  for(i=0; i<100000; i++) {
    *(Vij + i) = 4 * (1/(pow(*(rij2 + i),6))  - 1/ (pow(*(rij2 + i),3)));
  }

  for(i=0;i<3*N;i++) {
    *(r+i) = 0;
  }

//####################################################//


  set_box(r,N,L);   //Coloca las particulas en su posición inicial.
  set_v(v,N,T);     //Crea las velocidades iniciales.

  printf("Imprimo los rp's y luego las veces que sumo F's\n");
  fuerza(N, L, r, rij2, Fij, Fr);    //Calcula las fuerzas.

  printf("\n");



  //x(t+h) = x(t) + v(t)*h + 1/2 * F(t)/m * h*h


                    //Calcula las nuevas posiciones         <-----------
                    //Calculo de fuerzas.                              |
                    //Calcula las nuevas velocidades                  --
  printf("Fuerza final\n");
  for(i=0; i<3*N;i++) {
    printf("%f\n", *(Fr+i));//, *(Fr+i));
  }

  printf("\n");
  printf("Posiciones partículas\n");
  for(i=0; i<N;i++) {
    printf("%f %f %f\n", *(r+3*i), *(r+3*i+1), *(r+3*i+2));
  }


  free(r);
  free(rij);
  free(v);
  free(rij2);
  free(Fij);
  free(Vij);
  free(Fr);

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

float Gaussiana(float nu, float sigma) {
  int n = 10, i;
  float z = 0;

  for(i=0; i<n; i++) {
    z += (float)rand()/(float)RAND_MAX;
  }

  z = sqrt(12*n)*(z/n - 0.5);

  return z*sigma+nu;
}

float set_v(float *v, int N, float T) {
  int i,j;
  float sigma = sqrt(T);
  float *VCM;

  VCM = malloc(3*sizeof(float));

  *(VCM) = 0;
  *(VCM + 1) = 0;
  *(VCM + 2) = 0;

  for(i=0; i<3; i++) {
    *(v+i) = Gaussiana(0.0,sigma);
  }

  for(i=0; i<N; i++) {
    for(j=0; j<3; j++) {
      *(VCM + j) += *(v + 3*i+j)/N;
    }
  }

  for(i=0; i<N; i++) {
    for(j=0; j<3; j++) {
      *(v + 3*i+j) -= *(VCM + j);
    }
  }

  return 0;

}

float fuerza(int N, int L, float *r, float *rij2, float *Fij, float *Fr) {
  int i, j, k, *cont;
  float rp2, rpx2, rpy2, rpz2, ra, rb, Fa, Fb;

  cont = malloc((N-1)*sizeof(int));

  for(i=0; i<N-1; i++) {
    *(cont+i)=0;
    for(j=i+1; j<N; j++) {
      if(i!=j) {
          rpx2 = pow(*(r + 3*i) - *(r + 3*j),2);
          rpy2 = pow(*(r + 3*i+1) - *(r + 3*j+1),2);
          rpz2 = pow(*(r + 3*i+2) - *(r + 3*j+2),2);

          rp2 = rpx2 + rpy2 + rpz2;

          printf("%f %f %f\n", rpx2, rpy2, rpz2);
          printf("%f\n", rp2);
          printf("\n");


          if(rp2 < L*L/2) {
            *(cont+i)+=1;
            for(k=0; k<100000; k++) {
              if(rpx2>=*(rij2+k)) {
                ra = *(rij2+k-1);
                rb = *(rij2+k);
                Fa = *(Fij + k-1);
                Fb = *(Fij + k);

                *(Fr + 3*i) += Fa + (rpx2 - ra) * (Fb - Fa)/(rb - ra);

                break;
              }
            }
          }

          if(rpy2 < L*L/4) {
            for(k=0; k<100000; k++) {
              if(rpy2>=*(rij2+k)) {
                ra = *(rij2+k-1);
                rb = *(rij2+k);
                Fa = *(Fij + k-1);
                Fb = *(Fij + k);

                *(Fr + 3*i+1) += Fa + (rpy2 - ra) * (Fb - Fa)/(rb - ra);

                break;
              }
            }
          }

          if(rpz2 < L*L/4) {
            for(k=0; k<100000; k++) {
              if(rpz2>=*(rij2+k)) {
                ra = *(rij2+k-1);
                rb = *(rij2+k);
                Fa = *(Fij + k-1);
                Fb = *(Fij + k);

                *(Fr + 3*i+2) += Fa + (rpz2 - ra) * (Fb - Fa)/(rb - ra);

                break;
              }
            }
          }
      }
    }
  }

  printf("\n");

  for(i=0; i<N-1; i++) {
    printf("%i\n", *(cont+i));
  }

  return 0;
}
