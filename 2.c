#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double set_box(double *r, int N, double L);
double set_v(double *v, int N, double T);
double Gaussiana(double nu, double sigma);
double fuerza(int N, int L, double *r, double *Fr, double *Frv);
int posiciones(int N, int L, double *r, double *v, double *Fr);
double velocity_verlet(int N, double T, double *v, double *Fr, double *Frv);
double Coef_V(int N, int L,double *r);
double presion(int N, double L, double T, double *r, double *Fr);

int main(int argc, char const *argv[]) {
  int N=125, i, it=4000;
  double L, T, Vs=0.0, E_c=0.0, p=0.0;
  double *r, *v, *Fr, *Frv;
  FILE *presiones;

  srand(time(NULL));

  r = malloc(3*N*sizeof(double));
  v = malloc(3*N*sizeof(double));
  Fr = malloc(3*N*sizeof(double));
  Frv = malloc(3*N*sizeof(double));

  presiones = fopen("./presion", "a");

  for(i=0; i<3*N; i++) {
    *(Fr + i) = 0.0;
  }

for(L=5; L<7; L+=0.1) {
    fprintf(presiones, "%f\n", (double)N/(L*L*L));
    for(T=0.4; T<2; T+=0.1) {
      Vs = 0.0;
      E_c = 0.0;

      fprintf(presiones, "Temperatura %f\n", T);
      set_box(r,N,L);   //Coloca las particulas en su posición inicial.
      Vs = fuerza(N, L, r, Fr, Frv);  //Calculo de Vs.
      p = presion(N, L, T, r, Fr);
      E_c = set_v(v,N,T);     //Crea las velocidades iniciales.

      fprintf(presiones, "Datos %f %f\n", Vs+E_c, p);

      for(i=0; i<it; i++) {
        posiciones(N, L, r, v, Fr);            //Calcula las nuevas posiciones  <----
        E_c = 0.0;                                                                //|
        Vs = fuerza(N, L, r, Fr, Frv) ;         //Calculo de fuerzas.                |
        p = presion(N, L, T, r, Fr);
        E_c = velocity_verlet(N, T, v, Fr, Frv);  //Calcula las nuevas velocidades   ---

        fprintf(presiones, "Datos %f %f\n", Vs+E_c, p);

      }
    }
  }

  fclose(presiones);

  free(r);
  free(v);
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

  E_c /= N;

  return E_c;
}
//##########################################################//

//##########################################################//
double fuerza(int N, int L, double *r, double *Fr, double *Frv) {
  int i, j;
  double rp, rp3, rp6, rpx, rpy, rpz, Fs=0.0, Vs=0.0, corr_V=4.0 * (1.0/59605.0 - 1.0/244.0);

  for(j=0; j<3*N; j++) {                                                        //|
    *(Frv + j) = *(Fr + j);                  //Fuerza anterior.                   |
  }

  for(j=0; j<3*N; j++) {
    *(Fr + j) = 0.0;
  }

  for(i=0; i<N-1; i++) {
    for(j=i+1; j<N; j++) {

        rpx = *(r + 3*i) - *(r + 3*j);
        rpy = *(r + 3*i+1) - *(r + 3*j+1);
        rpz = *(r + 3*i+2) - *(r + 3*j+2);

        if(rpx > 0.5 * L) rpx -= L;
        if(rpy > 0.5 * L) rpy -= L;
        if(rpz > 0.5 * L) rpz -= L;

        if(rpx < -0.5 * L) rpx += L;
        if(rpy < -0.5 * L) rpy += L;
        if(rpz < -0.5 * L) rpz += L;

        rp = sqrt(rpx*rpx + rpy*rpy + rpz*rpz);

        rp3 = rp * rp * rp;
        rp6 = rp3 * rp3;

        if(rp < 2.5) {
            Fs = 24.0 * (2.0/ (rp * rp6*rp6) - 1.0/ (rp * rp6));
            Vs += 4.0 * (1.0/(rp6*rp6) - 1.0/rp6);
            Vs -= corr_V;

            *(Fr + 3*i) += Fs * (rpx/rp);
            *(Fr + 3*i+1) += Fs * (rpy/rp);
            *(Fr + 3*i+2) += Fs * (rpz/rp);

            *(Fr + 3*j) -= Fs * (rpx/rp);
            *(Fr + 3*j+1) -= Fs * (rpy/rp);
            *(Fr + 3*j+2) -= Fs * (rpz/rp);
            }
      }
    }


  Vs /= (double)N;

  return Vs;
}
//##########################################################//

//##########################################################//
int posiciones(int N, int L, double *r, double *v, double *Fr) {
  int i;
  double h=0.001;

  for(i=0; i<N; i++) {
    *(r+3*i) = *(r+3*i) + *(v+3*i) * h + 0.5 * (*(Fr+3*i)) * h*h;
    if(*(r+3*i) > L) *(r+3*i) -= L;
    if(*(r+3*i) < 0.0) *(r+3*i) += L;
    *(r+3*i+1) = *(r+3*i+1) + *(v+3*i+1) * h + 0.5 * (*(Fr+3*i+1)) * h*h;
    if(*(r+3*i+1) > L) *(r+3*i+1) -= L;
    if(*(r+3*i+1) < 0.0) *(r+3*i+1) += L;
    *(r+3*i+2) = *(r+3*i+2) + *(v+3*i+2) * h + 0.5 * (*(Fr+3*i+2)) * h*h;
    if(*(r+3*i+2) > L) *(r+3*i+2) -=L;
    if(*(r+3*i+2) < 0.0) *(r+3*i+2) += L;
  }


  return 0;
}
//##########################################################//

//##########################################################//
double velocity_verlet(int N, double T, double *v, double *Fr, double *Frv) {
  int i;
  double h=0.001, E_c=0.0, Ta, vx, vy, vz;

  for(i=0; i<N; i++) {
    *(v+3*i) = *(v+3*i) + (*(Frv+3*i) + *(Fr+3*i)) * h * 0.5;
    *(v+3*i+1) = *(v+3*i+1) + (*(Frv+3*i+1) + *(Fr+3*i+1)) * h * 0.5;
    *(v+3*i+2) = *(v+3*i+2) + (*(Frv+3*i+2) + *(Fr+3*i+2)) * h * 0.5;

    vx = *(v+3*i);
    vy = *(v+3*i+1);
    vz = *(v+3*i+2);

    Ta = 0.33 * (vx*vx + vy*vy + vz*vz);

    *(v+3*i) = (*(v+3*i)) * sqrt(T/Ta);
    *(v+3*i+1) = (*(v+3*i+1)) * sqrt(T/Ta);
    *(v+3*i+2) = (*(v+3*i)) * sqrt(T/Ta);

    E_c += 0.5 * ((*(v+3*i)) * (*(v+3*i)) + (*(v+3*i+1)) * (*(v+3*i+1)) + (*(v+3*i+2)) * (*(v+3*i+2)));
  }

  E_c /= (double)N;

  return E_c;
}
//##########################################################//

//##########################################################//
double Coef_V(int N, int L, double *r) {
  int i;
  double lambda_x=0.0, lambda_y=0.0, lambda_z=0.0, lambda_t=0.0, a=L/cbrt(N);

  for(i=0; i<N-1; i++) {
      lambda_x += cos((2*3.1415/a) * ((*(r+3*i)) - (a * 0.5)));
      lambda_y += cos((2*3.1415/a) * ((*(r+3*i+1)) - (a * 0.5)));
      lambda_z += cos((2*3.1415/a) * ((*(r+3*i+2)) - (a * 0.5)));
      }

    lambda_x /= (double)N;
    lambda_y /= (double)N;
    lambda_z /= (double)N;

    lambda_t = 0.33 * (lambda_x + lambda_y + lambda_z);

  return lambda_t;
}
//##########################################################//

//##########################################################//
double presion(int N, double L, double T, double *r, double *Fr) {
  int i,j;
  double rho=(double)N/(L*L*L), V=L*L*L, p, w=0, rpx, rpy, rpz, rp, Frx, Fry, Frz, Frs;

  for(i=0; i<N-1; i++) {
    for(j=i+1; j<N; j++) {
      rpx = *(r + 3*i) - *(r + 3*j);
      rpy = *(r + 3*i+1) - *(r + 3*j+1);
      rpz = *(r + 3*i+2) - *(r + 3*j+2);
      rp = sqrt(rpx*rpx + rpy*rpy + rpz*rpz);

      Frx = *(Fr + 3*i) - *(Fr + 3*j);
      Fry = *(Fr + 3*i+1) - *(Fr + 3*j+1);
      Frz = *(Fr + 3*i+2) - *(Fr + 3*j+2);
      Frs = sqrt(Frx*Frx + Fry*Fry + Frz*Frz);

      w += (rp*Frs)/N;
    }
  }

  p = rho*T + w/(3*V);

  return p;
}
