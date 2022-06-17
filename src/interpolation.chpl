module Interpolation
{
  use Random;
  use UnitTest;
  use Set;

  class interpolation_coefficients_c
  {
    var coefs_d : domain(2); // {nFPs, nSPs}
    var coefs: [coefs_d] real;
  }

  class derivation_coefficients_c
  {
    var coefs_d : domain(3); // {nFPs, nDims, nSPs}
    var coefs: [coefs_d] real;
  }

  // Define type for the interpolation structure. Add "?" to allow default initialization to nil.
  type interpolation_coefficients_t = unmanaged interpolation_coefficients_c?;
  type derivation_coefficients_t    = unmanaged derivation_coefficients_c?;

  // Domains
  var sp2fpInterp_d   : domain(2*int);
  var sp2spDeriv_d    : domain(2*int);
  var sp2nodeInterp_d : domain(2*int);

  // Coefficient structures
  var sp2fpInterp   : [sp2fpInterp_d] interpolation_coefficients_t; // Cell to face interpolation
  var sp2spDeriv    : [sp2spDeriv_d]  derivation_coefficients_t;    // Cell to face derivative
  var sp2nodeInterp : [sp2spDeriv_d]  interpolation_coefficients_t; // Cell to node interpolation

  // Alternative coefficient structures.
  //    1- It's assumed interpolation is always from the cellTopo SPs to somewhere else.
  //    2- Interpolating to a lower order is possible but obviously does not preserve order of accuracy.
  //    3- fromPtIdx is typically the cell SP index
  //    4-   toPtIdx can be either the cell vertices, the face FPs, or an arbitrary distribution of points in the edge.
  //
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toCell[          solOrder ].coefs[toPtIdx, fromPtIdx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toFace[(faceIdx, solOrder)].coefs[toPtIdx, fromPtIdx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toEdge[(edgeIdx, solOrder)].coefs[toPtIdx, fromPtidx]
  //interpolate[(cellTopo, spDistribution, spSolOrder)].toNode[(nodeIdx, solOrder)].coefs[toPtIdx, fromPtidx]
  //
  // Ex: interpolate[(TOPO_PRIS, 4)].toFace[(1, 4)].

  //////////////////////////////////
  //   Initialization Procedure   //
  //////////////////////////////////

  proc init_sp2fpInterp(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;
    use Mesh;

    writeln();
    writeln("Initializing SP -> FP Interpolation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2fpInterp_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2fpInterp.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;
          var fpCnt : int = 2;

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spDistLine : [1..spCnt] real = nodes_legendre_gauss(spCnt);
          var fpDistLine : [1..fpCnt] real = [-1.0, 1.0];

          for fp in 1..fpCnt do
            sp2fpInterp[(cellTopo, interpOrder)]!.coefs[{fp..#1, 1..spCnt}] = reshape(
                eval_LagrangePoly1D_array(fpDistLine[fp], spDistLine), {fp..#1, 1..spCnt});
        }
        //when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;
          var fpCnt : int = (interpOrder+1)*4;

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var fpDistLine : [1..2] real = [-1.0, 1.0];

          for cellFace in 1..elem_faces(TOPO_QUAD) do
            for faceFP in 1..interpOrder+1
            {
              var  xiFP : real;
              var etaFP : real;

              select cellFace
              {
                when 1
                {
                  xiFP  = spDistLine[faceFP];
                  etaFP = fpDistLine[1];
                }
                when 2
                {
                  xiFP  = fpDistLine[2];
                  etaFP = spDistLine[faceFP];
                }
                when 3
                {
                  xiFP  = spDistLine[interpOrder+2-faceFP];
                  etaFP = fpDistLine[2];
                }
                when 4
                {
                  xiFP  = fpDistLine[1];
                  etaFP = spDistLine[interpOrder+2-faceFP];
                }
              }

              for spIdx in 1..spCnt
              {
                var i : int = 1 + ( (spIdx-1)/(interpOrder+1)**0 ) %(interpOrder+1);
                var j : int = 1 + ( (spIdx-1)/(interpOrder+1)**1 ) %(interpOrder+1);

                var fpIdx : int = faceFP + (cellFace-1)*(interpOrder+1);

                sp2fpInterp[(cellTopo, interpOrder)]!.coefs[fpIdx, spIdx] = eval_LagrangePoly1D( xiFP, i, spDistLine)
                                                                           *eval_LagrangePoly1D(etaFP, j, spDistLine);
              }
            }
        }
        //when TOPO_TETR {}
        //when TOPO_PYRA {}
        //when TOPO_PRIS {}
        when TOPO_HEXA
        {
          var spCnt : int = (interpOrder+1)**3;
          var fpCnt : int = ((interpOrder+1)**2)*6;

          sp2fpInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..fpCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var fpDistLine : [1..2] real = [-1.0, 1.0];

          for cellFace in 1..elem_faces(TOPO_HEXA) do
            for faceFP in 1..(interpOrder+1)**2
            {
              var   xiFP : real;
              var  etaFP : real;
              var zetaFP : real;

              select cellFace
              {
                when 1
                {
                  xiFP   = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                  etaFP  = spDistLine[(faceFP-1)/(interpOrder+1)+1];
                  zetaFP = fpDistLine[1];
                }
                when 2
                {
                  xiFP   = spDistLine[(faceFP-1)/(interpOrder+1)+1];
                  etaFP  = fpDistLine[1];
                  zetaFP = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                }
                when 3
                {
                  xiFP   = fpDistLine[1];
                  etaFP  = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                  zetaFP = spDistLine[(faceFP-1)/(interpOrder+1)+1];
                }
                when 4
                {
                  xiFP   = spDistLine[(faceFP-1)/(interpOrder+1)+1];
                  etaFP  = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                  zetaFP = fpDistLine[2];
                }
                when 5
                {
                  xiFP   = spDistLine[((interpOrder+1)**2-faceFP)/(interpOrder+1)+1]; // Reverse order
                  etaFP  = fpDistLine[2];
                  zetaFP = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                }
                when 6
                {
                  xiFP   = fpDistLine[2];
                  etaFP  = spDistLine[(faceFP-1)/(interpOrder+1)+1];
                  zetaFP = spDistLine[(faceFP-1)%(interpOrder+1)+1];
                }
              }

              for spIdx in 1..spCnt
              {
                var i : int = 1 + ( (spIdx-1)/(interpOrder+1)**0 ) %(interpOrder+1);
                var j : int = 1 + ( (spIdx-1)/(interpOrder+1)**1 ) %(interpOrder+1);
                var k : int = 1 + ( (spIdx-1)/(interpOrder+1)**2 ) %(interpOrder+1);

                var fpIdx : int = faceFP + (cellFace-1)*(interpOrder+1)**2;

                sp2fpInterp[(cellTopo, interpOrder)]!.coefs[fpIdx, spIdx] = eval_LagrangePoly1D(  xiFP, i, spDistLine)
                                                                           *eval_LagrangePoly1D( etaFP, j, spDistLine)
                                                                           *eval_LagrangePoly1D(zetaFP, k, spDistLine);
              }
            }
        }
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  proc init_sp2spDeriv(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing SP -> SP' Differentiation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2spDeriv_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2spDeriv.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;

          sp2spDeriv[(cellTopo, interpOrder)] = new derivation_coefficients_t({1..spCnt, 1..1, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spLoc : [1..spCnt] real = nodes_legendre_gauss(spCnt);

          for sp in 1..spCnt do
            sp2spDeriv[(cellTopo, interpOrder)]!.coefs[{sp..#1, 1..#1, 1..spCnt}] =
                  reshape(eval_DLagrangeDx_array(spLoc[sp], spLoc), {sp..#1, 1..#1, 1..spCnt});
        }
        //when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;

          sp2spDeriv[(cellTopo, interpOrder)] = new derivation_coefficients_t({1..spCnt, 1..2, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

          for toSP in 1..spCnt
            {
              var  xi : real = spDistLine[1 + (toSP-1)%(interpOrder+1)];
              var eta : real = spDistLine[1 + (toSP-1)/(interpOrder+1)];

              for fromSP in 1..spCnt
              {
                var i : int = 1 + (fromSP-1)%(interpOrder+1);
                var j : int = 1 + (fromSP-1)/(interpOrder+1);

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 1, fromSP] = eval_DLagrangeDx(xi, i, spDistLine)
                                                                             *eval_LagrangePoly1D(eta, j, spDistLine);

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 2, fromSP] = eval_LagrangePoly1D(xi, i, spDistLine)
                                                                             *eval_DLagrangeDx(eta, j, spDistLine);
              }
            }
        }
        //when TOPO_TETR {}
        //when TOPO_PYRA {}
        //when TOPO_PRIS {}
        when TOPO_HEXA
        {
          var spCnt : int = (interpOrder+1)**3;

          sp2spDeriv[(cellTopo, interpOrder)] = new derivation_coefficients_t({1..spCnt, 1..3, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

          for toSP in 1..spCnt
            {
              var   xi : real = spDistLine[1 + ( (toSP-1)/(interpOrder+1)**0 ) %(interpOrder+1)];
              var  eta : real = spDistLine[1 + ( (toSP-1)/(interpOrder+1)**1 ) %(interpOrder+1)];
              var zeta : real = spDistLine[1 + ( (toSP-1)/(interpOrder+1)**2 ) %(interpOrder+1)];

              for fromSP in 1..spCnt
              {
                var i : int = 1 + ( (fromSP-1)/(interpOrder+1)**0 ) %(interpOrder+1);
                var j : int = 1 + ( (fromSP-1)/(interpOrder+1)**1 ) %(interpOrder+1);
                var k : int = 1 + ( (fromSP-1)/(interpOrder+1)**2 ) %(interpOrder+1);


                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 1, fromSP] = eval_DLagrangeDx(     xi, i, spDistLine)
                                                                             *eval_LagrangePoly1D( eta, j, spDistLine)
                                                                             *eval_LagrangePoly1D(zeta, k, spDistLine);

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 2, fromSP] = eval_LagrangePoly1D(  xi, i, spDistLine)
                                                                             *eval_DLagrangeDx(    eta, j, spDistLine)
                                                                             *eval_LagrangePoly1D(zeta, k, spDistLine);

                sp2spDeriv[(cellTopo, interpOrder)]!.coefs[toSP, 3, fromSP] = eval_LagrangePoly1D(  xi, i, spDistLine)
                                                                             *eval_LagrangePoly1D( eta, j, spDistLine)
                                                                             *eval_DLagrangeDx(   zeta, k, spDistLine);
              }
            }
        }
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  proc init_sp2nodeInterp(minOrder : int, maxOrder : int, cellTopos : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Polynomials;

    writeln();
    writeln("Initializing SP -> Node Differentiation matrices");
    writeln("    Cell Topologies: ", cellTopos);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for cellTopo in cellTopos do
      for interpOrder in minOrder..maxOrder do
        sp2nodeInterp_d.add((cellTopo, interpOrder));

    // Calculate all relevant coefficients
    for (cellTopo, interpOrder) in sp2nodeInterp.domain
    {
      select cellTopo
      {
        when TOPO_LINE
        {
          var spCnt : int = interpOrder+1;
          var nodeCnt : int = 2;//elem_vertices(cellTopo);

          sp2nodeInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..nodeCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots. xD
          var spDistLine : [1..spCnt] real = nodes_legendre_gauss(spCnt);
          var nodeDistLine : [1..nodeCnt] real = [-1.0, 1.0];

          for node in 1..nodeCnt do
            sp2nodeInterp[(cellTopo, interpOrder)]!.coefs[{node..#1, 1..spCnt}] = reshape(
                eval_LagrangePoly1D_array(nodeDistLine[node], spDistLine), {node..#1, 1..spCnt});
        }
        //when TOPO_TRIA {}
        when TOPO_QUAD
        {
          var spCnt : int = (interpOrder+1)**2;
          var nodeCnt : int = 4;//elem_vertices(cellTopo);

          sp2nodeInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..nodeCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var nodeDistLine : [1..2] real = [-1.0, 1.0];

          for cellNode in 1..nodeCnt
          {
            var  xiNode : real;
            var etaNode : real;

            select cellNode
            {
              when 1
              {
                xiNode  = nodeDistLine[1];
                etaNode = nodeDistLine[1];
              }
              when 2
              {
                xiNode  = nodeDistLine[2];
                etaNode = nodeDistLine[1];
              }
              when 3
              {
                xiNode  = nodeDistLine[2];
                etaNode = nodeDistLine[2];
              }
              when 4
              {
                xiNode  = nodeDistLine[1];
                etaNode = nodeDistLine[2];
              }
            }

            for spIdx in 1..spCnt
            {
              var i : int = (spIdx-1)%(interpOrder+1)+1;
              var j : int = (spIdx-1)/(interpOrder+1)+1;

              sp2nodeInterp[(cellTopo, interpOrder)]!.coefs[cellNode, spIdx] = eval_LagrangePoly1D( xiNode, i, spDistLine)
                                                                              *eval_LagrangePoly1D(etaNode, j, spDistLine);
            }
          }
        }
        //when TOPO_TETR {}
        //when TOPO_PYRA {}
        //when TOPO_PRIS {}
        when TOPO_HEXA
        {
          var spCnt : int = (interpOrder+1)**3;
          var nodeCnt : int = 8;//elem_vertices(cellTopo);

          sp2nodeInterp[(cellTopo, interpOrder)] = new interpolation_coefficients_t({1..nodeCnt, 1..spCnt})!;

          // Need to build an appropriate way to query the point location for each element.
          // Initially assume the whole mesh uses the same base distribution specified in input file.
          // Even more initially assume the whole mesh has SPs on Legendre roots.
          var spDistLine : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
          var nodeDistLine : [1..2] real = [-1.0, 1.0];

          for cellNode in 1..nodeCnt
          {
            var   xiNode : real;
            var  etaNode : real;
            var zetaNode : real;

            select cellNode
            {
              when 1
              {
                xiNode   = nodeDistLine[1];
                etaNode  = nodeDistLine[1];
                zetaNode = nodeDistLine[1];
              }
              when 2
              {
                xiNode   = nodeDistLine[2];
                etaNode  = nodeDistLine[1];
                zetaNode = nodeDistLine[1];
              }
              when 3
              {
                xiNode   = nodeDistLine[2];
                etaNode  = nodeDistLine[2];
                zetaNode = nodeDistLine[1];
              }
              when 4
              {
                xiNode   = nodeDistLine[1];
                etaNode  = nodeDistLine[2];
                zetaNode = nodeDistLine[1];
              }
              when 5
              {
                xiNode   = nodeDistLine[1];
                etaNode  = nodeDistLine[1];
                zetaNode = nodeDistLine[2];
              }
              when 6
              {
                xiNode   = nodeDistLine[2];
                etaNode  = nodeDistLine[1];
                zetaNode = nodeDistLine[2];
              }
              when 7
              {
                xiNode   = nodeDistLine[2];
                etaNode  = nodeDistLine[2];
                zetaNode = nodeDistLine[2];
              }
              when 9
              {
                xiNode   = nodeDistLine[1];
                etaNode  = nodeDistLine[2];
                zetaNode = nodeDistLine[2];
              }
            }

            for spIdx in 1..spCnt
            {
              var i : int = 1 + ( (spIdx-1)/(interpOrder+1)**0 ) %(interpOrder+1);
              var j : int = 1 + ( (spIdx-1)/(interpOrder+1)**1 ) %(interpOrder+1);
              var k : int = 1 + ( (spIdx-1)/(interpOrder+1)**2 ) %(interpOrder+1);

              sp2nodeInterp[(cellTopo, interpOrder)]!.coefs[cellNode, spIdx] = eval_LagrangePoly1D(  xiNode, i, spDistLine)
                                                                              *eval_LagrangePoly1D( etaNode, j, spDistLine)
                                                                              *eval_LagrangePoly1D(zetaNode, k, spDistLine);
            }
          }
        }
        otherwise do writeln("Unsupported mesh element found at interpolation initialization.");
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  //////////////////////////////////////////
  //   Lagrange Interpolation functions   //
  //////////////////////////////////////////

  proc eval_LagrangePoly1D(x : real, k : int, xi : [] real) : real
  {
    // Evaluate the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the point x.
    //
    //               nfp
    //              _____
    //               | |     x  - xi(i)
    //    L_k (x) =  | |  -------------
    //               i=1  xi(k) - xi(i)
    //              i/=k
    //
    //  x  : Coordinate at which to evaluate the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates
    //
    // The k-th basis vector of a Lagrange basis is a polynomial function that evaluates to 0 at all nodes xi[] except
    // at xi[k] where it evaluates to 1.

    var return_value : real = 1.0;

    for i in xi.domain do
      if i != k then
        return_value *= (x-xi[i]) / (xi[k]-xi[i]);

    return return_value;
  }

  proc eval_LagrangePoly1D_array(x : real, xi : [] real) : [] real
  {
    // Get the values of the 1D Lagrange polynomials at point x
    //
    //  x  : Coordinate at which to evaluate the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates
    //
    // This evaluates the Lagrange polynomial at the location x. The Lagrange polynomial corresponds to a function that
    // is 0.0 at all iterpolation points xi(i) for i=1,nfp except at xi(k) where the function is equal to 1.

    var return_value : [xi.domain] real;

    for k in xi.domain do
      return_value[k] = eval_LagrangePoly1D(x, k, xi);

    return return_value;
  }

  proc eval_DLagrangeDx(x : real, k : int, xi : [] real) : real
  {
    // Evaluate the derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the
    // point x.
    //
    //                    nfp                  nfp
    //                    ___                 _____
    //     d              \          1         | |     x  - xi(j)
    //    -- [L_k (x)] =  /__  -------------   | |  -------------
    //    dx                   xi(k) - xi(i)   j=1  xi(k) - xi(j)
    //                    i=1                 j/=i
    //                   i/=k                 j/=k
    //
    //  x  : Coordinate at which to evaluate the derivative of the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates

    var evalDLagrangeDx : real = 0.0;
    var xiSliced : [xi.domain.low..xi.domain.high-1] real;

    // Populate xiSliced with all nodes except the k-th
    for i in xiSliced.domain do
      if i<k then
        xiSliced[i] = xi[i];
      else
        xiSliced[i] = xi[i+1];

    for i in xiSliced.domain
    {
      var aux = xiSliced[i];
      xiSliced[i] = xi[k];
      evalDLagrangeDx += eval_LagrangePoly1D(x, i, xiSliced) / (xiSliced(i)-aux);
      xiSliced[i] = aux;
    }

    return evalDLagrangeDx;
  }

  proc eval_DLagrangeDx_array(x : real, xi : [] real) : [] real
  {
    // Evaluate the derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[] at the
    // point x.
    //
    //  x  : Coordinate at which to evaluate the derivative of the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates

    var evalDLagrangeDxArray : [xi.domain] real;

    for k in xi.domain do
      evalDLagrangeDxArray[k] = eval_DLagrangeDx(x, k, xi);

    return evalDLagrangeDxArray;
  }

  proc eval_D2LagrangeDx2(x : real, k : int, xi : [] real) : real
  {
    // Evaluate the second derivative of the k-th basis vector of the 1D Lagrange basis defined by the set of nodes xi[]
    // at the point x.
    //
    //                      nfp                nfp                  nfp
    //                      ___                ___                 _____
    //     d²               \         1        \          1         | |     x  - xi(n)
    //    --- [L_k (x)]  =  /__ -------------  /__  -------------   | |  -------------
    //    dx²                   xi(k) - xi(l)       xi(k) - xi(m)   n=1  xi(k) - xi(n)
    //                      l=1                m=1                 n/=l
    //                     l/=k               m/=l                 n/=m
    //                                        m/=k                 n/=k
    //
    //  x  : Coordinate at which to evaluate the second derivative of the k-th Lagrange basis polynomial
    //  k  : Index of the basis vector being evaluated
    //  xi : Array of all the interpolation nodes coordinates

    var evalD2LagrangeDx2 : real = 0.0;
    var xiSliced : [xi.domain.low..xi.domain.high-1] real;

    // Populate xiSliced with all nodes except the k-th
    for i in xiSliced.domain do
      if i<k then
        xiSliced[i] = xi[i];
      else
        xiSliced[i] = xi[i+1];

    for i in xiSliced.domain do
    {
      var aux = xiSliced[i];
      xiSliced[i] = xi[k];
      evalD2LagrangeDx2 += eval_DLagrangeDx(x, i, xiSliced) / (xiSliced(i)-aux);
      xiSliced[i] = aux;
    }

    return evalD2LagrangeDx2;
  }

  proc eval_D2LagrangeDx2_array(x : real, k : int, xi : real) : [xi.domain] real
  {
    // Get value of second derivative of k-th Lagrange polynomial at point x.
    //
    //  x  : Coordinate at which to evaluate the second derivative of the k-th Lagrange basis polynomial
    //  xi : Array of all the interpolation nodes coordinates

    var evalD2LagrangeDx2Array : [xi.domain] real = 0.0;

    for k in xi.domain do
      evalD2LagrangeDx2Array[k] = eval_D2LagrangeDx2(x, k, xi[xi.domain]);

    return evalD2LagrangeDx2Array;
  }

  proc eval_LagrangePoly2D(x : [1..2] real, k : [1,.2] int, xi : [] real, eta : [] real) : real
  {
    return  eval_LagrangePoly1D(x[1], k[1],  xi[])
           *eval_LagrangePoly1D(x[2], k[2], eta[]);
  }

  proc eval_LagrangePoly3D(x : [1..3] real, k : [1..2] int, xi : [] real, eta : [] real, zeta : [] real) : real
  {
    return  eval_LagrangePoly1D(x[1], k[1],   xi[])
           *eval_LagrangePoly1D(x[2], k[2],  eta[])
           *eval_LagrangePoly1D(x[3], k[3], zeta[]);
  }

  ///////////////////////////////
  //   Module Test Procedure   //
  ///////////////////////////////

  proc main()
  {
    use IO;
    use Testing;
    use Parameters.ParamTest;
    use Parameters.ParamMesh;
    use Polynomials;

    var randStream = new RandomStream(real);
    var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

    // Test the primitive Lagrange functions
    {
      // Define test parameters
      var minDegree : int = 0;
      var maxDegree : int = 9;

      for interpDegree in minDegree..maxDegree
      {
        writef("Degree %2i test:\n", interpDegree);

        // Generate test polynomials by randomizing the coefficients if the canonical form
        var coef : [0..interpDegree] real;
        randStreamSeeded.fillRandom(coef);
        write("Test Poly:  ");
        for idx in coef.domain do
          writef(" %{+.4dr}*x^%i", coef[idx], idx);
        writeln();

        // Generate a random set of interpolation nodes
        var gaussNodes : [0..interpDegree] real = nodes_legendre_gauss(interpDegree+1);
        var randNodes : [0..interpDegree] real;
        randStreamSeeded.fillRandom(randNodes);
        write("Gauss Nodes: ");
        for idx in gaussNodes.domain do
          writef(" %{ .4dr}", gaussNodes[idx]);
        writeln();
        write("Rand Nodes:  ");
        for idx in randNodes.domain do
          writef(" %{ .4dr}", randNodes[idx]);
        writeln();

        // Evaluate the test polynomial at the interpolation nodes
        var interpPolyGauss : [0..interpDegree] real = 0;
        var interpPolyRand  : [0..interpDegree] real = 0;
        for nodeIdx in gaussNodes.domain do
          for deg in 0..interpDegree do
            interpPolyGauss[nodeIdx] += coef[deg]*gaussNodes[nodeIdx]**deg;
        for nodeIdx in randNodes.domain do
          for deg in 0..interpDegree do
            interpPolyRand[nodeIdx] += coef[deg]*randNodes[nodeIdx]**deg;

        // Get a random coordinate and evaluate the test polynomial at it
        var xTest : real = randStreamSeeded.getNext();
        var yDirect : real = 0;
        var dydxDirect : real = 0;
        var d2ydx2Direct : real = 0;
        for deg in 0..interpDegree do
          yDirect += coef[deg]*xTest**deg;
        for deg in 1..interpDegree do
          dydxDirect += deg*coef[deg]*xTest**(deg-1);
        for deg in 2..interpDegree do
          d2ydx2Direct += deg*(deg-1)*coef[deg]*xTest**(deg-2);

        // Interpolate the polynomial to the test coordinate
        var yInterpGauss : real = 0;
        var yInterpRand  : real = 0;
        for nodeIdx in gaussNodes.domain do
          yInterpGauss += interpPolyGauss[nodeIdx]*eval_LagrangePoly1D(xTest, nodeIdx, gaussNodes);
        for nodeIdx in randNodes.domain do
          yInterpRand  += interpPolyRand[nodeIdx]*eval_LagrangePoly1D(xTest, nodeIdx, randNodes);

        // Interpolate the first derivative to the test coordinate
        var dydxInterpGauss : real = 0;
        var dydxInterpRand  : real = 0;
        for nodeIdx in gaussNodes.domain do
          dydxInterpGauss += interpPolyGauss[nodeIdx]*eval_DLagrangeDx(xTest, nodeIdx, gaussNodes);
        for nodeIdx in randNodes.domain do
          dydxInterpRand  += interpPolyRand[nodeIdx]*eval_DLagrangeDx(xTest, nodeIdx, randNodes);

        // Interpolate the second derivative to the test coordinate
        var d2ydx2InterpGauss : real = 0;
        var d2ydx2InterpRand  : real = 0;
        for nodeIdx in gaussNodes.domain do
          d2ydx2InterpGauss += interpPolyGauss[nodeIdx]*eval_D2LagrangeDx2(xTest, nodeIdx, gaussNodes);
        for nodeIdx in randNodes.domain do
          d2ydx2InterpRand  += interpPolyRand[nodeIdx]*eval_D2LagrangeDx2(xTest, nodeIdx, randNodes);

        // Core Functions
        // Evaluate the 'k'-th component of the basis defined by the 'xi' interpolation nodes at the point 'x'
        //eval_LagrangePoly1D(x : real, k : int, xi : [] real)
        //eval_DLagrangeDx(x : real, k : int, xi : [] real)
        //eval_D2LagrangeDx2(x : real, k : int, xi : real)

        // Compare interpolated polynomial to direct evaluation of the test polynomial
        writef("Gauss Nodes: x=%.4dr, yDirect=%.8dr, yInterp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, yDirect, yInterpGauss, abs(yDirect-yInterpGauss), abs((yDirect-yInterpGauss)/yDirect));
        writef("Rand  Nodes: x=%.4dr, yDirect=%.8dr, yInterp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, yDirect, yInterpRand, abs(yDirect-yInterpRand), abs((yDirect-yInterpRand)/yDirect));
        writef("Gauss Nodes: x=%.4dr, dydxDirect=%.8dr, dydxInterp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, dydxDirect, dydxInterpGauss, abs(dydxDirect-dydxInterpGauss), abs((dydxDirect-dydxInterpGauss)/dydxDirect));
        writef("Rand  Nodes: x=%.4dr, dydxDirect=%.8dr, dydxInterp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, dydxDirect, dydxInterpRand, abs(dydxDirect-dydxInterpRand), abs((dydxDirect-dydxInterpRand)/dydxDirect));
        writef("Gauss Nodes: x=%.4dr, d2ydx2Direct=%.8dr, d2ydx2Interp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, d2ydx2Direct, d2ydx2InterpGauss, abs(d2ydx2Direct-d2ydx2InterpGauss), abs((d2ydx2Direct-d2ydx2InterpGauss)/d2ydx2Direct));
        writef("Rand  Nodes: x=%.4dr, d2ydx2Direct=%.8dr, d2ydx2Interp=%.8dr, Abs Error=%.4er, Rel Error=%.4er\n",
            xTest, d2ydx2Direct, d2ydx2InterpRand, abs(d2ydx2Direct-d2ydx2InterpRand), abs((d2ydx2Direct-d2ydx2InterpRand)/d2ydx2Direct));
        writeln();

        // Multi Dimensional versions
        //eval_LagrangePoly2D(x : [1..2] real, k : [1,.2] int, xi : [] real, eta : [] real)
        //eval_LagrangePoly3D(x : [1..3] real, k : [1..2] int, xi : [] real, eta : [] real, zeta : [] real)

        // Array versions
        //eval_LagrangePoly1D_array(x : real, xi : [] real)
        //eval_DLagrangeDx_array(x : real, xi : [] real)
        //eval_D2LagrangeDx2_array(x : real, k : int, xi : real)
      }
    }

    // Define FR parameters
    var minDegree : int = 0;
    var maxDegree : int = 9;
    var cellTopos : set(int);
    cellTopos.add(TOPO_LINE);
    cellTopos.add(TOPO_QUAD);
    cellTopos.add(TOPO_HEXA);

    // Initialize the FR structures
    init_sp2fpInterp(  minDegree, maxDegree, cellTopos);
    init_sp2spDeriv(   minDegree, maxDegree, cellTopos);
    init_sp2nodeInterp(minDegree, maxDegree, cellTopos);

    // Test the FR structures
    {
      //// Generate test polynomials by randomizing the coefficients if the canonical form
      //var coef : [0..maxDegree, 0..maxDegree, 0..maxDegree] real;
      //randStreamSeeded.fillRandom(coef);

      //writeln();
      //writeln("Testing sp2fpInterp");
      //writeln();
      //for cellTopo in cellTopos do
      //  for interpDegree in minDegree..maxDegree
      //  {
      //    // Generate a set of SP coordinates
      //    var xyzSP : [1..3, 1..n_cell_sps(cellTopo, interpDegree)] real;
      //    for i in 1..

      //    // Evaluate the test polynomial at the interpolation nodes
      //    var solSP : [1..n_cell_sps(cellTopo, interpDegree)] real;
      //    for sol in solSP
      //    {
      //      sol = xyzSP * coef;
      //    }

      //    // Get a random coordinate and evaluate the test polynomial at it
      //    // Interpolate the polynomial to the test coordinate
      //    // Interpolate the first derivative to the test coordinate
      //  }

      //writeln();
      //writeln("Testing sp2spDeriv");
      //writeln();
      //for cellTopo in cellTopos do
      //  for interpOrder in minOrder..maxOrder
      //  {
      //  }

      //writeln();
      //writeln("Testing sp2nodeInterp");
      //writeln();
      //for cellTopo in cellTopos do
      //  for interpOrder in minOrder..maxOrder
      //  {
      //  }
    }
  }
}
