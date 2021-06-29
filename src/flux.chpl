prototype module Flux
{
  proc pressure_cv(cons : [] real) : real
  {
    import Input.fGamma;
    import LinearAlgebra.dot;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var pressure : real = (fGamma-1.0)*(cons[idxEne] - 0.5*dot(cons[idxMom],cons[idxMom])/cons[idxRho]);

    return pressure;
  }

  proc temperature_cv(cons : [] real) : real
  {
    import Input.fR;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density

    var p : real = pressure_cv(cons);

    var temperature : real = p/(cons[idxRho]*fR);

    return temperature;
  }

  proc enthalpy_cv(cons : [] real) : real
  {
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var p : real = pressure_cv(cons);

    var enthalpy : real = (cons[idxEne] + p);

    return enthalpy;
  }

  proc internal_energy_cv(cons : [] real) : real
  {
    import Input.fCv;
    import Input.fR;

    var p : real = pressure_cv(cons);

    var internalEnergy : real = fCv*p/fR;

    return internalEnergy;
  }

  proc sound_speed_cv(cons : [] real) : real
  {
    import Input.fGamma;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var p : real = pressure_cv(cons);

    return sqrt(fGamma*p/cons[idxRho]);
  }

  proc mach_cv(cons : [] real) : real
  {
    import LinearAlgebra.norm;

    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var vel : [idxMom] real = cons[idxMom]/cons[idxRho];

    var mach : real = norm(vel) / sound_speed_cv(cons);

    return mach;
  }

  proc ener_pv(prim : [] real) : real
  {
    import Input.fGamma;
    import LinearAlgebra.dot;

    var idxDens : int   = prim.domain.dim(0).low;           // First element is density
    var idxVel  : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxPres : int   = prim.domain.dim(0).high;          // Last element is energy

    var ener : real = prim[idxPres]/(fGamma-1) + 0.5*prim[idxDens]*dot(prim[idxVel], prim[idxVel]);

    return ener;
  }

  proc convection_flux_cv(cons : real) : [] real
  {
    import Input.convectionSpeed;

    return cons*convectionSpeed;
  }

  proc euler_flux_cv_1d(cons : [1..3] real) : [1..3] real
  {
    var euler_flux_cv : [cons.domain] real;
    var p : real = pressure_cv(cons);

    euler_flux_cv[1] = cons[2];
    euler_flux_cv[2] = cons[2]*cons[2]/cons[1] + p;
    euler_flux_cv[3] = cons[2]/cons[1]*(cons[3] + p);

    return euler_flux_cv;
  }

  proc euler_flux_cv(cons : [] real) : [] real
  {
    var idxRho : int   = cons.domain.dim(0).low;           // First element is density
    var idxMom : range = cons.domain.dim(0).expand(-1);    // Intermediary elements are the velocities
    var idxEne : int   = cons.domain.dim(0).high;          // Last element is energy

    var euler_flux_cv : [idxMom-1, cons.domain.dim(0)] real;

    var vel : [idxMom] real = cons[idxMom] / cons[idxRho];
    var p   : real = pressure_cv(cons);

    for i in idxMom-1
    {
      euler_flux_cv[i, idxRho] = cons[i+1];
      euler_flux_cv[i, idxMom] = cons[i+1]*vel[idxMom];
      euler_flux_cv[i, idxEne] = vel[i+1]*(cons[idxEne] + p);

      euler_flux_cv[i, i+1] += p;
    }

    return euler_flux_cv;
  }

  proc euler_flux_pv(prim : [] real) : [] real
  {
    use Input;
    import LinearAlgebra.dot;

    var idxRho : int   = prim.domain.dim(0).low;           // First element is density
    var idxVel : range = prim.domain.dim(0).expand(-1);    // Intermediary elements are the velocity components
    var idxPre : int   = prim.domain.dim(0).high;          // Last element is pressure

    var euler_flux_pv : [idxVel-1, prim.domain.dim(0)] real;

    for i in idxVel-1
    {
      euler_flux_pv[i, idxRho] = prim[idxRho+i]*prim[idxRho];
      euler_flux_pv[i, idxVel] = prim[idxRho+i]*prim[idxRho]*prim[idxVel];
      euler_flux_pv[i, idxPre] = prim[idxRho+i]*(prim[idxPre]*fGamma/(fGamma-1.0)+0.5*prim[idxRho]*dot(prim[idxVel],prim[idxVel]));

      euler_flux_pv[i, i+1] += prim[idxPre];
    }

    return euler_flux_pv;
  }

  proc visc_flux_cv(cons : [] real, consGrad : [] real) : [] real
  {
  }

  proc visc_flux_pv(prim : [] real, promGrad : [] real) : [] real
  {
  }

  proc main()
  {
    var cons1d : [1..3] real = [1.225, 250.1160830494510, 278846.40];
    var cons2d : [1..4] real = [1.225, 249.9637190895680, 8.728925415626780, 278846.40];
    var cons3d : [1..5] real = [1.225, 249.9637190895680, 8.728925415626780, 0.0, 278846.40];

    var prim1d : [1..3] real = [1.225, 204.1763943260830, 101325.0];
    var prim2d : [1..4] real = [1.225, 204.0520155833210, 7.125653400511660, 101325.0];
    var prim3d : [1..5] real = [1.225, 204.0520155833210, 7.125653400511660, 0.0, 101325.0];

    writeln("Conserverd variables, Density (kg/m³), XYZ Momentum components (kg*m/s*m³), Energy (J/m³):");
    writeln("1D: ", cons1d);
    writeln("2D: ", cons2d);
    writeln("3D: ", cons3d);
    writeln();
    writeln("Primitive Variables, Density (kg/m³), XYZ Velocity components (m/s), Pressure (Pa):");
    writeln("1D: ", prim1d);
    writeln("2D: ", prim2d);
    writeln("3D: ", prim3d);
    writeln();
    writeln("Pressure and Temperature functions:");
    writeln("1D - Pressure (Pa): ", pressure_cv(cons1d), ", Temperature (K): ", temperature_cv(cons1d));
    writeln("2D - Pressure (Pa): ", pressure_cv(cons2d), ", Temperature (K): ", temperature_cv(cons2d));
    writeln("3D - Pressure (Pa): ", pressure_cv(cons3d), ", Temperature (K): ", temperature_cv(cons3d));
    writeln();
    writeln("Enthalpy and Internal Energy functions:");
    writeln("1D - Enthalpy (J/m³): ", enthalpy_cv(cons1d), ", Internal Energy (J/m³): ", internal_energy_cv(cons1d));
    writeln("2D - Enthalpy (J/m³): ", enthalpy_cv(cons2d), ", Internal Energy (J/m³): ", internal_energy_cv(cons2d));
    writeln("3D - Enthalpy (J/m³): ", enthalpy_cv(cons3d), ", Internal Energy (J/m³): ", internal_energy_cv(cons3d));
    writeln();
    writeln("Mach and Speed of Sound functions:");
    writeln("1D - Mach (non-dimensional): ", mach_cv(cons1d), ", Speed of Sound (m/s): ", sound_speed_cv(cons1d));
    writeln("2D - Mach (non-dimensional): ", mach_cv(cons2d), ", Speed of Sound (m/s): ", sound_speed_cv(cons2d));
    writeln("3D - Mach (non-dimensional): ", mach_cv(cons3d), ", Speed of Sound (m/s): ", sound_speed_cv(cons3d));
    writeln();
    writeln("Euler 1D Flux:\n", euler_flux_cv_1d(cons1d));
    writeln();
    writeln("Generic Euler flux function (conserver vars):");
    writeln("  Euler 1D Flux:\n", euler_flux_cv(cons1d));
    writeln("  Euler 2D Flux:\n", euler_flux_cv(cons2d));
    writeln("  Euler 3D Flux:\n", euler_flux_cv(cons3d));
    writeln();
    writeln("Generic Euler flux function (primitive vars):");
    writeln("  Euler 1D Flux:\n", euler_flux_pv(prim1d));
    writeln("  Euler 2D Flux:\n", euler_flux_pv(prim2d));
    writeln("  Euler 3D Flux:\n", euler_flux_pv(prim3d));
  }
}
