prototype module Init
{
  use Parameters.ParamInput;
  use Input;
  use Config;

  proc initial_condition(IC : int, const ref xyz : [] real) : [] real
  {
    var sol : [xyz.domain.dim(0), 1..3] real;
    select IC
    {
      when IC_SHOCKTUBE
      {
        var xMin : real = 0.0;
        var xMax : real = 1.0;

        var rhoHi : real = 2.0;
        var eHi   : real = 2.0;
        var rhoLo : real = 1.0;
        var eLo   : real = 1.0;

        for i in xyz.domain.dim(0)
        {
          if xyz[i,1] < 0.5*(xMin+xMax) then
          {
            sol[i,1] = rhoHi;
            sol[i,2] = 0.0;
            sol[i,3] = eHi;
          }
          else
          {
            sol[i,1] = rhoLo;
            sol[i,2] = 0.0;
            sol[i,3] = eLo;
          }
        }
      }
      when IC_1D_NOZZLE
      {
        const xT : real = 0.4;
        const aE : real = 0.4;
        const aT : real = 0.2;
        var areaCrit : real = nozzle_ratio(xT); // Smooth transonic flow = Throat Area

        for i in xyz.domain.dim(0)
        {
          const rho0 = 1;
          const p0 = 1;
          const t0 = 1;

          var area = nozzle_ratio(xyz[i,1]);
          var mach = nozzle_mach_sub(area);

          var aux = 1 + (fGamma-1.0)*mach**2/2;

          var rho = rho0 * aux**(-1/(fGamma-1));
          var p = p0 * aux**(-fGamma/(fGamma-1));
          var t = t0 * aux**(-1);

          var a = sqrt(fGamma*p/rho);

          sol[i,1] = rho;  // Density
          sol[i,2] = rho*mach*a;  // Momentum
          sol[i,3] = p*(1/(fGamma-1) + fGamma*mach**2);  // Energy
        }
      }
      otherwise {}
    }

    return sol;
  }

  proc nozzle_ratio(in x : real) : real
  {
    // Nozzle shape (from "I do Like CFD", vol 1 by Katate Masatsuka (Hiroaki Nishikawa))
    // A(x) = 25/9(Ae-At)(x-0.4)^2 + At
    // aE = Exit area
    // aT = Throat area
    const xT : real = 0.4;
    const aE : real = 0.4;
    const aT : real = 0.2;

    return (25.0/9.0*(aE-aT)*(x-xT)**2 + aT)/aT;
  }

  proc nozzle_mach_sub(in areaRatio : real) : real
  {
    use IO;

    var machLo : real = 0.0;
    var machHi : real = 1.0;
    var mach : real = 0.5*(machLo + machHi);

    while (mach != machLo) && (mach != machHi)
    {
      var res = (1.0/mach)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res > areaRatio**2 then
        machLo = mach;
      else
        machHi = mach;

      mach = 0.5*(machLo + machHi);
    }

    return mach;
  }

  proc nozzle_mach_sup(in areaRatio : real) : real
  {
    var machLo : real =   1.0;
    var machHi : real = 100.0;
    var mach : real = 0.5*(machLo + machHi);

    while (mach != machLo) && (mach != machHi)
    {
      var res = (1.0/mach)*((2.0+(fGamma-1.0)*mach**2)/(fGamma+1.0))**((fGamma+1.0)/(fGamma-1.0));

      if res < areaRatio**2 then
        machLo = mach;
      else
        machHi = mach;

      mach = 0.5*(machLo + machHi);
    }

    return mach;
  }

  proc entropy_wave_1d(x : real) : real
  {
    
  }

  proc main()
  {
    use IO;
    use Parameters.ParamInput;
    use Flux;

    const nNodes = 21;
    var xyz : [1..nNodes, 1..3] real;
    var sol : [1..nNodes, 1..3] real;

    writeln("1D Meshes");

    for i in xyz.domain.dim(0)
    {
      xyz[i, 1] = (i-1.0)/(nNodes-1.0);
      xyz[i, 2] = 0.0;
      xyz[i, 3] = 0.0;
    }

    writeln();
    writeln("1D Shock Tube");
    sol = initial_condition(IC_SHOCKTUBE, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3]");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1], xyz[i, 2], xyz[i, 3],
                                                                                    sol[i, 1], sol[i, 2], sol[i, 3]));

    writeln();
    writeln("1D Nozzle Flow");
    sol = initial_condition(IC_1D_NOZZLE, xyz);
    writeln("Point #,    X-Coord,    Y-Coord,    Z-Coord,     Sol[1],     Sol[2],     Sol[3],   Pressure,       Mach");
    for i in xyz.domain.dim(0) do
      writeln("%7i, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er, %10.3er".format(i, xyz[i, 1],
            xyz[i, 2], xyz[i, 3], sol[i, 1], sol[i, 2], sol[i, 3], pressure_cv(sol[i,..]), mach_cv(sol[i,..])));
  }
}
