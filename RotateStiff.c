#include "radiats.h"

/*******************************************************************************
* Compute Stiffness tensor components in a rotated coordinate system 0=x,1=y,2=z.
*	alpha, phi	tilt angle and azimuth (in degrees.) of the old z axis.
*******************************************************************************/
void RotateStiff(float alpha, float phi, float stiff[6][6]) {
  int i,j,k,l,ij,kl,p,q,s,t;
  float rot[3][3],c[3][3][3][3];
  for(i=0;i<3;i++) {
  for(j=0;j<3;j++) {
     ij = indxij[i][j];
     for(k=0;k<3;k++) {
     for(l=0;l<3;l++) {
        kl = indxij[k][l];
        c[i][j][k][l] = stiff[ij][kl];
     }}
  }}
  // rotation transformation
  rotationMatrix(alpha, phi, rot);
  for(ij=0;ij<6;ij++) {
     i = indx_i[ij];
     j = indx_j[ij];
     for(kl=ij;kl<6;kl++) {
	k = indx_i[kl];
	l = indx_j[kl];
        stiff[ij][kl] = 0.;
        for(p=0; p<3; p++) {
        for(q=0; q<3; q++) {
        for(s=0; s<3; s++) {
        for(t=0; t<3; t++) {
           stiff[ij][kl] += rot[i][p]*rot[j][q]*rot[k][s]*rot[l][t]*c[p][q][s][t];
        }}}}
	stiff[kl][ij] = stiff[ij][kl];
     }
  }
}


/*******************************************************************************
* Compute the rotation matrix rot[3][3] from (x',y',z') to (x,y,z) where
* y'=z'Xz (so y' is in the x-y plane). 0<=alpha<=180 is the angle between z'
* and z and 0<=phi<=360 is the angle of y' clockwise from y.
* First rotate (x',z') clockwise around y' by alpha:
*  |x"|   | cos(alpha) 0 sin(alpha) | |x'|
*  |y'| = |     0      1      0     |*|y'| 
*  |z |   |-sin(alpha) 0 cos(alpha) | |z'|,
* then rotate (x",y') counter-clockwise around z by phi:
*  |x|   | cos(phi) -sin(phi) 0 | |x"|
*  |y| = | sin(phi)  cos(phi) 0 |*|y'| 
*  |z|   |     0       0      1 | |z |.
* So,
*  |x|   | cos(phi)*cos(alpha) -sin(phi) cos(phi)*sin(alpha) | |x'|
*  |y| = | sin(phi)*cos(alpha)  cos(phi) sin(phi)*sin(alpha) |*|y'| 
*  |z|   |    -sin(alpha)         0          cos(alpha)      | |z'|.
* Vector transformation:
*	x[i] = rot[i][j] * x'[j],
* and tensor transformation:
*    T[i][j] = rot[i][p]*rot[j][q]*T'[p][q].
*
* In an x=N,y=E,z=D geographic coordinates, phi is the azimuth of z'-axis from
* North and alpha is its tilt angle from the vertical down.
*******************************************************************************/
void rotationMatrix(float alpha, float phi, float rot[3][3]) {
  alpha *= DEG2RAD;
  phi *= DEG2RAD;
  rot[2][2] = cos(alpha);
  rot[0][2] = sin(alpha)*cos(phi);
  rot[1][2] = sin(alpha)*sin(phi);
  rot[2][0] =-sin(alpha);
  rot[0][0] = cos(alpha)*cos(phi);
  rot[1][0] = cos(alpha)*sin(phi);
  rot[2][1] = 0.;
  rot[0][1] =-sin(phi);;
  rot[1][1] = cos(phi);
}
