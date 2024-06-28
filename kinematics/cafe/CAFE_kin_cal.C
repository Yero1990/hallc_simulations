// Write a code to calculate the kinematic for CAFE experiment
// This can be used later
#include <iostream>
using namespace std;

const double Mp = 0.938; // Proton mass
void CAFE_kin_cal(double E0=10.6, double xb=0.976, double Q2=2.1) {

  
  double Ee = 0;
  double theta_e_rad = 0;
  double nu = 0;
  double nu_cyero = 0;
  double theta_e_deg = 0;

  
  nu = Q2 / 2 / Mp / xb;
  nu_cyero = Q2 / (xb * 2. * Mp);
  
  std::cout << "nu = " << nu < endl;
  //cout << "nu_cyero = " << nu_cyero << endl;

  /*  
  Ee = E0 - nu;

  theta_e_rad = 2 * asin(sqrt(Q2 / 4 / E0 / Ee));
  theta_e_deg = theta_e_rad * 180. / 3.14;

  cout << "Electron kine: Ee, theta_e_deg: " << Ee << "  " << theta_e_deg
       << endl;
  */
}

/*
void Cal_x_Q2(double E0, double theta_e_deg, double Ee) {
  double nu = E0 - Ee;

  double theta_e_rad = theta_e_deg * 3.14 / 180.;

  double sin2 = sin(theta_e_rad / 2) * sin(theta_e_rad / 2);

  double Q2 = 4 * E0 * Ee * sin2;

  double xb = Q2 / 2 / Mp / nu;

  cout << "xb, Q2: " << xb << "  " << Q2 << endl;
}
*/
