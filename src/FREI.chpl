/* Documentation for FREI */
prototype module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";
  config const inputMesh : string = "mesh.mesh";

  proc main() {
    use IO;
    use Time;
    use Parameters.ParamInput;
    use Config;
    use Input;
    use Flux;
    use Riemann;
    use Interpolation;
    use Output;
    use Gmesh;
    use Mesh;
    use FRMesh;
    use Boundary;
    use Correction;
    use Init;
    use Quadrature;
    use Projection;
    use Limiter;
    use FR;
    use LinearAlgebra;
    use SourceTerm;
    use Temporal_Methods;

    // Timing variables
    var initTime : real = 0.0;
    var iterTime : real = 0.0;
    var srcTermTime : real = 0.0;
    var dscFluxTime : real = 0.0;
    var cntFluxTime : real = 0.0;
    var timeStepTime : real = 0.0;
    var stabilizeTime : real = 0.0;
    var stopwatch : Timer;
    var iterTimer : Timer;
    stopwatch.start();

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();
    var timeStepStages : int;
    select timeScheme
    {
      when TIME_TVDRK_O2S2 do
       timeStepStages = 2;
      when TIME_TVDRK_O2S3 do
       timeStepStages = 3;
      when TIME_TVDRK_O2S4 do
       timeStepStages = 4;
      when TIME_TVDRK_O2SN do
       timeStepStages = 5;
      when TIME_TVDRK_O3S3 do
       timeStepStages = 3;
      when TIME_TVDRK_O3S4 do
       timeStepStages = 4;
      when TIME_TVDRK_O3S5 do
       timeStepStages = 5;
      when TIME_TVDRK_O4S5 do
       timeStepStages = 5;
      otherwise do
       timeStepStages = 1;
    }

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c();
    select Input.meshFormat
    {
      when MESH_GENERATE
      {
        select Input.meshingScheme
        {
          when MESH_GEN_UNIFORM do
            gmesh2.uniform1D(Input.nCells, Input.xMin, Input.xMax);
          when MESH_GEN_RANDOM do
            gmesh2.random1D(Input.nCells, Input.xMin, Input.xMax);
        }
      }
      when MESH_GMESH do
        gmesh2.read_gmesh_file(Input.meshFileName);
      when MESH_CGNS {}
    }

    // 4. Convert input mesh to solver mesh
    var frMesh = new unmanaged fr_mesh_c(nDims=gmesh2.mesh_dimension(), nVars=Input.nEqs, solOrder=Input.iOrder);
    frMesh.import_gmesh2(gmesh2);   // Convert mesh to native format
    frMesh.set_families(famlList);  // Get families data from input file and write to mesh

    // 5. Initialize FR mesh
    frMesh.allocate_fr_vars();      // Allocate SP and FP solution/flux/residue arrays
    frMesh.set_points_locations();  // Calculate coordinate transformations and point coordinates

    // 6. Save mesh file in internal format

    // 7. Initialize the FR solver, pre calculate coefficients and stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_quadratureWeights(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_polyProj(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction(Input.minOrder+1, Input.maxOrder+1, frMesh.cellTopos);

    // Initialize solution
    for cellIdx in frMesh.cellList.domain
    {
      ref familyIdx = frMesh.cellList[cellIdx].family;

      ref familyType = frMesh.famlList[familyIdx].bocoType;
      ref familySubType = frMesh.famlList[familyIdx].bocoSubType;
      ref familyParameters = frMesh.famlList[familyIdx].bocoProperties;

      ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
      ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

      frMesh.solSP[cellSPini.. #cellSPcnt, ..] = flow_condition(familySubType,
                                                                familyParameters,
                                                                frMesh.xyzSP[cellSPini.. #cellSPcnt, ..]);
    }

    // Output initial state
    iterOutput(iteration, frMesh);

    // Stabilize Solution
    {
      // Loop through cells
      for cellIdx in frMesh.cellList.domain
      {
        var cellSPini = frMesh.cellSPidx[cellIdx, 1];
        var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

        for varIdx in 1..frMesh.nVars
        {
          var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                        cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                        solDegree = iOrder);

          if stableDegree < iOrder then
            frMesh.solSP[cellSPini.. #cellSPcnt, varIdx] = projection_limiter(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                                              solDegree = iOrder,
                                                                              projDegree = stableDegree);
        }
      }
    }

    // Save restart file

    // Initialize convergence monitoring variables
    var l2DeltaIni           : [1..frMesh.nVars] real;
    var l2RelativeDeltaIni   : [1..frMesh.nVars] real;
    var convergenceLog : file;
    try {
      convergenceLog = open("convengence.dat" , iomode.cw);
    } catch {
      stdout.writeln("Unknown Error opening convergence log file.");
      stderr.writeln("Unknown Error opening convergence log file.");
    }
    var convergenceLogChan = convergenceLog.writer();

    initTime = stopwatch.elapsed(TimeUnits.milliseconds);
    writef("Stopwatch - Init    : %10.2dr ms\n", initTime);
    writef("Start Iterating\n");

    // Solve flow
    iterTimer.start();
    for iteration in 1..Input.maxIter
    {
      iterTimer.clear();

      // Zero out residue
      frMesh.resSP = 0.0;

      // Save initial solution
      frMesh.oldSolSP = frMesh.solSP;

      for stage in 1..timeStepStages
      {
        // Calculate residue for this iteration
        {
          // The residue has 3 components:
          //   1. Continuous Flux
          //   2. Discontinuous Flux
          //   3. Source terms
          //
          // The residual array is reset in the time stepping procedure

          // Component 1: Source Term
          {
            stopwatch.clear();

            for spIdx in frMesh.resSP.domain.dim(0) do
              frMesh.resSP[spIdx..#1, ..] = -source_term(frMesh.xyzSP[spIdx..#1, ..],
                                                         frMesh.solSP[spIdx..#1, ..],
                                                         Input.eqSet                )
                                            * frMesh.jacSP[spIdx];

            srcTermTime += stopwatch.elapsed(TimeUnits.milliseconds);
          }

          // Component 2: Discontinuous Flux
          {
            stopwatch.clear();

            // Calculate flux at SPs and it´s divergence
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              // Allocate temporary flux array
              var flxSP : [cellSPini.. #cellSPcnt, 1..frMesh.nVars] real;

              // Calculate fluxes
              for meshSP in cellSPini.. #cellSPcnt do
                select Input.eqSet
                {
                  when EQ_CONVECTION do
                    flxSP[meshSP, ..] = convection_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_INVBURGERS do
                    flxSP[meshSP, ..] = burgers_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_EULER do
                    flxSP[meshSP, ..] = euler_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                  when EQ_QUASI_1D_EULER do
                    flxSP[meshSP, ..] = euler_flux_cv_1d(frMesh.solSP[meshSP, ..]);
                }

              // Interpolate fluxes to FPs
              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP = cellFace;
                  frMesh.flxFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), iOrder)]!.coefs(cellFP, ..),
                                                           flxSP[cellSPini..#cellSPcnt,..]                              );
                }
              }

              // Convert fluxes from physical to computational domain.
              // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determiant
              for meshSP in cellSPini.. #cellSPcnt do
                flxSP[meshSP, ..] = dot(flxSP[meshSP, ..], frMesh.metSP[meshSP, 1, 1]**(-1))*frMesh.jacSP[meshSP];

              // Calculate flux divergence
              for cellSP in 1..cellSPcnt
              {
                var meshSP = cellSPini + cellSP - 1;
                frMesh.resSP[meshSP,..] += dot(sp2spDeriv[(thisCell.elemTopo(), iOrder)]!.coefs(cellSP, ..),
                                               flxSP[cellSPini..#cellSPcnt,..]                             );
              }
            }
            dscFluxTime += stopwatch.elapsed(TimeUnits.milliseconds);
          }

          // Component 3: Continuous Flux
          {
            stopwatch.clear();

            // Interpolate solution to FPs
            for cellIdx in frMesh.cellList.domain
            {
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP = cellFace;
                  frMesh.solFP[meshFP, faceSide, ..] = dot(sp2fpInterp[(thisCell.elemTopo(), iOrder)]!.coefs(cellFP, ..),
                                                           frMesh.solSP[cellSPini..#cellSPcnt,..]                       );
                }
              }
            }

            // Apply boundary conditions
            for faceIdx in frMesh.faceList.domain
            {
              // Get loop variables
              var thisFace = frMesh.faceList[faceIdx];
              var faceFPini = frMesh.faceFPidx[faceIdx, 1];
              var faceFPcnt = frMesh.faceFPidx[faceIdx, 2];

              // Check if the face´s right neighbor is a Boundary Condition
              if frMesh.faceList[faceIdx].cells[2] < 0
              {
                // Yep, it is, lets get some local iteration variables
                var bocoIdx = -frMesh.faceList[faceIdx].cells[2];
                var thisBoco = frMesh.bocoList[bocoIdx];

                var famlIdx = thisBoco.family;
                var thisFaml = frMesh.famlList[famlIdx];

                // Iterate through the FPs on this face
                for meshFP in faceFPini.. #faceFPcnt
                {
                  // Calculate the Boundary Condition using the solution at the left neighbor´s corresponding FP
                  frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml);
                }
              }
            }

            // Calculate interface correction
            for cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              var thisCell = frMesh.cellList[cellIdx];
              var cellSPini = frMesh.cellSPidx[cellIdx, 1];
              var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              for cellFace in thisCell.faces.domain
              {
                var faceIdx  = thisCell.faces[cellFace];
                var thisFace = frMesh.faceList[faceIdx];
                var faceSide = thisCell.sides[cellFace];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  // For 1D each face has 1 FP therefore the FP and the Face have the same index Relative to it's
                  // position in the cell
                  var cellFP = cellFace;

                  // Calculate the flux jump = -1*(local_flux) + numerical_flux
                  var jump : [1..frMesh.nVars] real = -frMesh.flxFP[meshFP, faceSide, ..];
                  select Input.eqSet
                  {
                    when EQ_CONVECTION do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_INVBURGERS do
                      jump += upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_QUASI_1D_EULER do
                      jump += roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_EULER do
                      jump += roe(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                  }

                  // Convert fluxes from physical to computational domain.
                  // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determiant
                  jump[..] = dot(jump[..], frMesh.metFP[meshFP, faceSide, 1, 1]**(-1))*frMesh.jacFP[meshFP, faceSide];

                  if faceSide == 1 && cellFace == 1 then
                    jump = -jump;
                  if faceSide == 2 && cellFace == 2 then
                    jump = -jump;

                  // The correction function was calculated in the computational domain already, therefore no
                  // transformation is required.
                  frMesh.resSP[cellSPini.. #cellSPcnt, ..] += outer(
                      flux_correction[(thisCell.elemTopo(), iOrder+1)]!.correction[cellFP, 1..cellSPcnt],
                      jump[..]);
                }
              }
            }
            cntFluxTime += stopwatch.elapsed(TimeUnits.milliseconds);
          }
        }

        // Advance RK Stage
        {
          stopwatch.clear();

          // Loop through cells
          for cellIdx in frMesh.cellList.domain
          {
            var cellSPini = frMesh.cellSPidx[cellIdx, 1];
            var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

            // Calculate dt for this cell
            var dt : real = Input.timeStep;
            //if variableTimeStep then
            //  dt = time_step();

            // Convert the residual from the computational to the physical domain
            for meshSP in cellSPini.. #cellSPcnt do
              frMesh.resSP[meshSP, ..] /= frMesh.jacSP[meshSP];


            // Update solution
            frMesh.solSP[cellSPini.. #cellSPcnt, ..] = time_advance(frMesh.oldSolSP[cellSPini.. #cellSPcnt, ..],
                                                                    frMesh.solSP[cellSPini.. #cellSPcnt, ..],
                                                                    frMesh.resSP[cellSPini.. #cellSPcnt, ..],
                                                                    dt, stage, Input.timeScheme);
          }
          timeStepTime += stopwatch.elapsed(TimeUnits.milliseconds);
        }

        // Stabilize Solution
        {
          stopwatch.clear();

          // Loop through cells
          for cellIdx in frMesh.cellList.domain
          {
            var cellSPini = frMesh.cellSPidx[cellIdx, 1];
            var cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

            for varIdx in 1..frMesh.nVars
            {
              var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                            cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                            solDegree = iOrder);

              if stableDegree < iOrder then
                frMesh.solSP[cellSPini.. #cellSPcnt, varIdx] = projection_limiter(solPoly = frMesh.solSP[cellSPini.. #cellSPcnt, varIdx],
                                                                                  cellTopo = frMesh.cellList[cellIdx].elemTopo(),
                                                                                  solDegree = iOrder,
                                                                                  projDegree = stableDegree);
            }
          }

          stabilizeTime += stopwatch.elapsed(TimeUnits.milliseconds);
        }
      }

      // Print solver status / log

      // Save restart file

      // Calculate and print convergence metrics
      {
        var l1Delta              : [1..frMesh.nVars] real;
        var l2Delta              : [1..frMesh.nVars] real;
        var lInfDelta            : [1..frMesh.nVars] real;
        var l1RelativeDelta      : [1..frMesh.nVars] real;
        var l2RelativeDelta      : [1..frMesh.nVars] real;
        var lInfRelativeDelta    : [1..frMesh.nVars] real;

        // Calculate solution delta from previous iteration
        for varIdx in 1..frMesh.nVars
        {
          l1Delta[varIdx]           = + reduce (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l2Delta[varIdx]           = sqrt(+ reduce (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])**2);
          lInfDelta[varIdx]         = max reduce abs(frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l1RelativeDelta[varIdx]   = + reduce (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx]);
          l2RelativeDelta[varIdx]   = sqrt(+ reduce (frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])**2);
          lInfRelativeDelta[varIdx] = max reduce abs((frMesh.oldSolSP[.., varIdx] - frMesh.solSP[.., varIdx])/frMesh.oldSolSP[.., varIdx]);
        }

        // Save values from first iterations as reference
        if iteration == 1
        {
          l2DeltaIni         = l2Delta;
          l2RelativeDeltaIni = l2RelativeDelta;
        }

        // Output summarized convergence metrics to stdOut
        writef("Iteration %9i | Time %{ 10.2dr}ms | Log10(L2(ΔSol)/L2(ΔSol0)) = %{ 7.4dr}", iteration,
            iterTimer.elapsed(TimeUnits.milliseconds), log10(norm(l2Delta)/norm(l2DeltaIni)));

        // Output full state to log file
        log_convergence(convergenceLogChan, iteration, l1Delta, l2Delta, lInfDelta, l1RelativeDelta, l2RelativeDelta, lInfRelativeDelta);

        if iteration % ioIter == 0 then
          writef(" | Saving solution file\n");
        else
          writef("\n");
      }

      // Check if we should write the solution this iteration
      if iteration % ioIter == 0 then
        iterOutput(iteration, frMesh);

      // Check if input file changed
    }

    // Output the final solution
    //iterOutput(iteration, frMesh);

    writeln();
    writef("Time splits:\n");
    writef("  Stopwatch - Init    : %11.2dr ms\n", initTime);
    writef("  Stopwatch - Src Term: %11.2dr ms\n", srcTermTime);
    writef("  Stopwatch - Dsc Flux: %11.2dr ms\n", dscFluxTime);
    writef("  Stopwatch - Cnt Flux: %11.2dr ms\n", cntFluxTime);
    writef("  Stopwatch - Stabiliz: %11.2dr ms\n", stabilizeTime);
    writef("  Stopwatch - Timestep: %11.2dr ms\n", timeStepTime);

    writeln();
    writeln("Fin");
  }
}
