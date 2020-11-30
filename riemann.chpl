prototype module Riemann

  proc Roe(const ref uL, uR : [1..3] real(64)) : [1..3] real(64)
  {
    use Input;

    var gm1 : real(64) = fGamma - 1.0;
    var gp1 : real(64) = fGamma + 1.0;

    //Primitive and other variables.
    //  Left state
    var rhoL : real(64) = uL[1];
    var vL   : real(64) = uL[2]/uL[1];
    var pL   : real(64) = pressure(uL);
    var aL   : real(64) = sqrt(fGamma*pL/rhoL);
    var HL   : real(64) = ( uL[3] + pL ) / rhoL;
    //  Right state
    var rhoR : real(64) = uR[1];
    var vR   : real(64) = uR[2]/uR[1];
    var pR   : real(64) = pressure(uR);
    var aR   : real(64) = sqrt(fGamma*pR/rhoR);
    var HR   : real(64) = ( uR[3] + pR ) / rhoR;

    //First compute the Roe Averages
    var RT  : real(64) = sqrt(rhoR/rhoL);
    var rho : real(64) = RT*rhoL;
    var v   : real(64) = (vL+RT*vR)/(1.0+RT);
    var H   : real(64) = (HL+RT*HR)/(1.0+RT);
    var a   : real(64) = sqrt( (fGamma-1.0)*(H-0.5*v*v) );

    //Differences in primitive variables.
    var drho = rhoR - rhoL;
    var du   =   vR - vL;
    var dP   =   pR - pL;

    //Wave strength (Characteristic Variables).
    var dV : [1..3] real(64);
    dV[1] =  0.5*(dP-rho*a*du)/(a*a);
    dV[2] = -( dP/(a*a) - drho );
    dV[3] =  0.5*(dP+rho*a*du)/(a*a);

    //Absolute values of the wave speeds (Eigenvalues)
    var ws : [1..3] real(64);
    ws[1] = abs(v-a);
    ws[2] = abs(v  );
    ws[3] = abs(v+a);

    //Modified wave speeds for nonlinear fields (to remove expansion shocks).
    //There are various ways to implement an entropy fix. This is just one
    //example.
    var Da : real(64) = max(0.0, 4.0*((vR-aR)-(vL-aL)) );
    if (ws(1) < 0.5*Da) then ws(1) = ws(1)*ws(1)/Da + 0.25*Da;
    Da = max(0.0, 4.0*((vR+aR)-(vL+aL)) );
    if (ws(3) < 0.5*Da) then ws(3) = ws(3)*ws(3)/Da + 0.25*Da;

    //Right eigenvectors
    var R : [1..3,1..3] real(64);
    R[1,1] = 1.0;
    R[2,1] = v - a;
    R[3,1] = H - v*a;

    R[1,2] = 1.0;
    R[2,2] = v;
    R[3,2] = 0.5*v*v;

    R[1,3] = 1.0;
    R[2,3] = v + a;
    R[3,3] = H + v*a;

    //Compute the average flux.
    Roe = 0.5*( physical_flux(uL) + physical_flux(uR) );

    //!Add the matrix dissipation term to complete the Roe flux.
    for j in 1..3 do {
      for k in 1..3 do {
        Roe[j] = Roe[j] - 0.5*ws[k]*dV[k]*R[j,k];
      }
    }
  }
}
