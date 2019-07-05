#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double set_box(double *r, int N, double L);
double set_v(double *v, int N, double T);
double Gaussiana(double nu, double sigma);
double fuerza(int N, int L, double *r, double *rij, double *rij2, double *Fij, double *Vij, double *Fr);
int posiciones(int N, double *r, double *v, double *Fr);
double velocity_verlet(int N, double *v, double *Fr, double *Frv);

int main(int argc, char const *argv[]) {
  int N=8, i, j, it=3;
  double L=4.0, T=6.0, Vs=0.0, E_c=0.0;
  double *r, *v, *rij, *rij2, *Fij, *Vij, *Fr, *Frv;
  FILE *video, *energias;

  srand(time(NULL));

  r = malloc(3*N*sizeof(double));
  v = malloc((3*N+3)*sizeof(double));
  Fr = malloc(3*N*sizeof(double));
  Frv = malloc(3*N*sizeof(double));

  video = fopen("./video.lammpstrj", "a");
  energias = fopen("./energias", "a");

  for(i=0; i<3*N; i++) {
    *(Fr + i) = 0.0;
  }

//############Tabla para interpolar###################//

  rij = malloc(100000*sizeof(double));
  rij2 = malloc(100000*sizeof(double));
  Fij = malloc(100000*sizeof(double));
  Vij = malloc(100000*sizeof(double));

  for(i=0; i<99999; i++) {
    *(rij + i) = (2.5/99999.0) + (2.5/99999.0) * i;
  }

  for(i=0; i<99999; i++) {
    *(rij2 + i) = *(rij + i) * *(rij + i);
  }

  for(i=0; i<99999; i++) {
    *(Fij + i) = 24.0 * (2.0/ (*(rij+i) * pow(*(rij2 + i),6))  - 1.0/ (*(rij+i) * pow(*(rij2 + i),3)));
  }

  for(i=0; i<99999; i++) {
    *(Vij + i) = 4.0 * (1.0/(pow(*(rij2 + i),6))  - 1.0/ (pow(*(rij2 + i),3)));
  }

  for(i=0;i<3*N;i++) {
    *(r+i) = 0;
  }

//####################################################//


  set_box(r,N,L);   //Coloca las particulas en su posición inicial.
  Vs = fuerza(N, L, r, rij, rij2, Fij, Vij, Fr);  //Calculo de Vs.
  E_c = set_v(v,N,T);     //Crea las velocidades iniciales.

  fprintf(video, "ITEM: TIMESTEP\n");
  fprintf(video, "0\n");
  fprintf(video, "ITEM: NUMBER OF ATOMS\n");
  fprintf(video, "%i\n", N);
  fprintf(video, "ITEM: BOX BOUNDS pp pp pp\n");
  fprintf(video, "0 10.000000\n");
  fprintf(video, "0 10.000000\n");
  fprintf(video, "0 10.000000\n");
  fprintf(video, "ITEM: ATOMS id x y z vx vy vz \n");

  fprintf(energias, "%i %f %f %f\n", 0, Vs, E_c, Vs+E_c);

  for(i=0; i<N; i++) {
    fprintf(video, "%i %f %f %f %f %f %f\n", i, *(r+3*i), *(r+3*i+1), *(r+3*i+2), *(v+3*i), *(v+3*i+1), *(v+3*i+2));
  }

  for(i=0; i<it; i++) {
    fprintf(video, "ITEM: TIMESTEP\n");
    fprintf(video, "%i\n", i+1);
    fprintf(video, "ITEM: NUMBER OF ATOMS\n");
    fprintf(video, "%i\n", N);
    fprintf(video, "ITEM: BOX BOUNDS pp pp pp\n");
    fprintf(video, "0 10.000000\n");
    fprintf(video, "0 10.000000\n");
    fprintf(video, "0 10.000000\n");
    fprintf(video, "ITEM: ATOMS id x y z vx vy vz\n");


    posiciones(N, r, v, Fr);         //Calcula las nuevas posiciones  <----
    for(j=0; j<3*N; j++) {                                              //|
      *(Frv + j) = *(Fr + j);        //Fuerza anterior.                 //|
    }                                                                   //|
    Vs = fuerza(N, L, r, rij, rij2, Fij, Vij, Fr);  //Calculo de fuerzas.                |
    E_c = velocity_verlet(N, v, Fr, Frv);      //Calcula las nuevas velocidades   ---


    for(j=0; j<N; j++) {
      fprintf(video, "%i %f %f %f %f %f %f\n", j, *(r+3*j), *(r+3*j+1), *(r+3*j+2), *(v+3*j), *(v+3*j+1), *(v+3*j+2));
    }

    fprintf(energias, "%i %f %f %f\n", i+1, Vs, E_c, Vs+E_c );

  }

  fclose(video);
  fclose(energias);

  free(r);
  free(rij);
  free(v);
  free(rij2);
  free(Fij);
  free(Vij);
  free(Fr);
  free(Frv);

  return 0;
}


//Ecuaciones Secundarias:
//##########################################################//
double set_box(double *r, int N, double L) {
  int n = cbrt(N), i=0, x, y, z;
  double dL=L/n;

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
//##########################################################//

//##########################################################//
double Gaussiana(double nu, double sigma) {
  int n = 10, i;
  double z = 0;

  for(i=0; i<n; i++) {
    z += (double)rand()/(double)RAND_MAX;
  }

  z = sqrt(12*n)*(z/n - 0.5);

  return z*sigma+nu;
}
//##########################################################//

//##########################################################//
double set_v(double *v, int N, double T) {
  int i,j;
  double sigma = sqrt(T), E_c=0.0;
  double *VCM;

  VCM = malloc(3*sizeof(double));

  *(VCM) = 0;
  *(VCM + 1) = 0;
  *(VCM + 2) = 0;

  for(i=0; i<3*N; i++) {
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

  for(i=0; i<N; i++) {
    E_c += 0.5 * (*(v+3*i) * *(v+3*i) + *(v+3*i+1) * *(v+3*i+1) + *(v+3*i+2) * *(v+3*i+2));
  }

  return E_c;
}
//##########################################################//

//##########################################################//
double fuerza(int N, int L, double *r, double *rij, double *rij2, double *Fij, double *Vij, double *Fr) {
  int i, j, k;
  double rp, rp2, rpx, rpy, rpz, ra, rb, Fa, Fb, Fs, Va, Vb, Vs=0.0, dr2=pow(2.5,2)*pow(99999,2);

  for(i=0; i<N-1; i++) {
    for(j=i+1; j<N; j++) {

        rpx = *(r + 3*i) - *(r + 3*j);
        rpy = *(r + 3*i+1) - *(r + 3*j+1);
        rpz = *(r + 3*i+2) - *(r + 3*j+2);

        rp = sqrt(rpx*rpx + rpy*rpy + rpz*rpz);
        rp2 = sqrt(rpx*rpx + rpy*rpy + rpz*rpz);

        if(rp < L) {
          if(rp < -L/2) {
            rp2 = (rpx+L)*(rpx+L) + (rpy+L)*(rpy+L) + (rpz+L)*(rpz+L);

            k = (int)((rp2-dr2)/dr2);
            ra = *(rij2 + k-1);
            rb = *(rij2 + k);
            Fa = *(Fij + k-1);
            Fb = *(Fij + k);
            Va = *(Vij + k-1);
            Vb = *(Vij + k);

            Fs = Fa + (sqrt(rp2) - ra) * (Fb - Fa)/(rb - ra);
            Vs += Va + (sqrt(rp2) - ra) * (Vb - Va)/(rb - ra);

            *(Fr + 3*i) += Fs * (rpx+L)/rp2;
            *(Fr + 3*i+1) += Fs * (rpy+L)/rp2;
            *(Fr + 3*i+2) += Fs * (rpz+L)/rp2;

            //printf("rpx: %f rij: %fn", rpx, *(rij + k));

            *(Fr + 3*j) -= Fs * rpx/rp2;
            *(Fr + 3*j+1) -= Fs * rpy/rp2;
            *(Fr + 3*j+2) -= Fs * rpz/rp2;
            }

          if ((rp > L/2)) {
            rp2 = (rpx-L)*(rpx-L) + (rpy-L)*(rpy-L) + (rpz-L)*(rpz-L);

            k = (int)((rp2-dr2)/dr2);
            ra = *(rij2 + k-1);
            rb = *(rij2 + k);
            Fa = *(Fij + k-1);
            Fb = *(Fij + k);
            Va = *(Vij + k-1);
            Vb = *(Vij + k);


            Fs = Fa + (sqrt(rp2) - ra) * (Fb - Fa)/(rb - ra);
            Vs += Va + (sqrt(rp2) - ra) * (Vb - Va)/(rb - ra);

            *(Fr + 3*i) += Fs * (rpx - L)/rp2;
            *(Fr + 3*i+1) += Fs * (rpy - L)/rp2;
            *(Fr + 3*i+2) += Fs * (rpz - L)/rp2;

            *(Fr + 3*j) -= Fs * (rpx - L)/rp2;
            *(Fr + 3*j+1) -= Fs * (rpy - L)/rp2;
            *(Fr + 3*j+2) -= Fs * (rpz - L)/rp2;
            }
          }
      }
    }

  return Vs;
}
//##########################################################//

//##########################################################//
int posiciones(int N, double *r, double *v, double *Fr) {
  int i;
  double h=0.01;

  for(i=0; i<N; i++) {
    *(r+3*i) = *(r+3*i) + *(v+3*i) * h + (1/2) * *(Fr+3*i) * h*h;
    *(r+3*i+1) = *(r+3*i+1) + *(v+3*i+1) * h + (1/2) * *(Fr+3*i+1) * h*h;
    *(r+3*i+2) = *(r+3*i+2) + *(v+3*i+2) * h + (1/2) * *(Fr+3*i+2) * h*h;
  }

  return 0;
}
//##########################################################//

//##########################################################//
double velocity_verlet(int N, double *v, double *Fr, double *Frv) {
  int i;
  double h=0.01, E_c=0.0;

  for(i=0; i<N; i++) {
    //printf("%f %fn", *(Frv+3*i), *(Fr+3*i));
    *(v+3*i) = *(v+3*i) + (*(Frv+3*i) + *(Fr+3*i)) * h;
    *(v+3*i+1) = *(v+3*i+1) + (*(Frv+3*i+1) + *(Fr+3*i+1)) * h;
    *(v+3*i+2) = *(v+3*i+2) + (*(Frv+3*i+2) + *(Fr+3*i+2)) * h;

    E_c += 0.5 * (*(v+3*i) * *(v+3*i) + *(v+3*i+1) * *(v+3*i+1) + *(v+3*i+2) * *(v+3*i+2));
  }

  return E_c;
}
//##########################################################//
