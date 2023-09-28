/* Documentation for FREI */
module FREI
{
  // Run-time constants
  config const inputFile : string = "input.toml";
  config const inputMesh : string = "mesh.mesh";

  proc main() {
    use IO;
    use Time;
    use Parameters.ParamInput;
    use Parameters.ParamMesh;
    use Config;
    use Input;
    use Dimensional;
    use Flux;
    use Riemann;
    use Interpolation;
    use Output;
    use ErrorCalc;
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
    use Temporal;
    import Math.log10;

    /////////////////////////////////
    // Declare profiling variables //
    /////////////////////////////////

    // Main program stowatch
    var programWatch  : stopwatch;
    var programTime   : real = 0.0;

    // Program steps stopwatch
    var totalWatch    : stopwatch;
    var initTime      : real = 0.0;
    var solveTime     : real = 0.0;
    var outputTime    : real = 0.0;

    // Iteration time stopwatch
    var solveWatch    : stopwatch;
    var oldSolTime    : real = 0.0;
    var residueTime   : real = 0.0;
    var timeStepTime  : real = 0.0;
    var stabilizeTime : real = 0.0;
    var ioIterTime    : real = 0.0;

    // Stopwatch and variables for residue's components
    var residueWatch : stopwatch;
    var srcTermTime : real = 0.0;
    var dscFluxTime : real = 0.0;
    var cntFluxTime : real = 0.0;

    // Stopwatch and variables for discontinuous flux component's internal steps
    //var dscFluxWatch : stopwatch;
    //var dscFluxTime1 : real = 0.0;
    //var dscFluxTime2 : real = 0.0;
    //var dscFluxTime3 : real = 0.0;
    //var dscFluxTime4 : real = 0.0;

    // Stopwatch and variables for continuous flux component's internal steps
    var cntFluxWatch : stopwatch;
    var cntFluxTime1 : real = 0.0;
    var cntFluxTime2 : real = 0.0;
    var cntFluxTime3 : real = 0.0;

    // Stopwatch and variables for the components of the 3rd step of the continuous flux calculation
    //var correctionWatch : stopwatch;
    //var riemTime : real = 0.0;
    //var jumpTime : real = 0.0;
    //var corrTime : real = 0.0;

    // Start stopwatches
    programWatch.start();
    totalWatch.start();

    ///////////////////////////
    // Solver Initialization //
    ///////////////////////////

    // 0. Initialize iteration count and stopwatch
    var iterWatch : stopwatch;
    var lastIter  : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();
    init_scales(lengRef=Input.lengRef, velMRef=Input.velMRef, tempRef=Input.tempRef, presRef=Input.presRef);
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
    const gmesh2 = new unmanaged gmesh2_c();
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
    //const frMesh = new unmanaged fr_mesh_c(nDims=gmesh2.mesh_dimension(), nVars=Input.nEqs, solOrder=Input.iOrder);
    var frMesh = new unmanaged fr_mesh_c(mesh=gmesh2, nVars=Input.nEqs, solOrder=Input.iOrder);
    frMesh.import_gmesh2(gmesh2);   // Convert mesh to native format
    frMesh.set_families(famlList);  // Get families data from input file and write to mesh

    // 5. Initialize FR mesh
    frMesh.allocate_fr_vars();      // Allocate SP and FP solution/flux/residue arrays
    frMesh.set_points_locations();  // Calculate coordinate transformations and point coordinates
    frMesh.build_cell_char_leng();  // Calculate the characteristic length of the mesh cells

    // 6. Save mesh file in internal format

    // 7. Initialize the FR solver, pre calculate coefficients and stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2nodeInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_quadratureWeights(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_polyProj(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction(Input.minOrder+1, Input.maxOrder+1, frMesh.cellTopos);

    // 8.a Initialize solution
    forall cellIdx in frMesh.cellList.domain
    {
      const ref familyIdx = frMesh.cellList[cellIdx].family;

      const ref familyType = frMesh.famlList[familyIdx].bocoType;
      const ref familySubType = frMesh.famlList[familyIdx].bocoSubType;
      const ref familyParameters = frMesh.famlList[familyIdx].bocoProperties;

      const ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
      const ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];
      const ref thisCell = frMesh.cellList[cellIdx];

      frMesh.solSP[.., cellSPini.. #cellSPcnt] = flow_condition(familySubType                           ,
                                                                familyParameters                        ,
                                                                frMesh.xyzSP[cellSPini.. #cellSPcnt, ..]).T;
    }
    // 8.b Interpolate solution to FPs
    forall cellIdx in frMesh.cellList.domain
    {
      const ref thisCell = frMesh.cellList[cellIdx];
      const ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
      const ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

      for cellFace in thisCell.faces.domain
      {
        const ref faceIdx  : int = thisCell.faces[cellFace];
        const ref faceSide : int = thisCell.sides[cellFace];
        const ref thisFace = frMesh.faceList[faceIdx];

        for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
        {
          var cellFP : int;
          if faceSide == 1 then
            cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
          else
            cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

          frMesh.solFP[meshFP, faceSide, ..] = 0;
          for varIdx in 1..frMesh.nVars do
            for cellSP in 1.. #cellSPcnt do
              frMesh.solFP[meshFP, faceSide, varIdx] += frMesh.solSP[varIdx, cellSPini+cellSP-1]
                                                       *sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, cellSP];
        }
      }
    }
    // 8.c Calculate ghost solutions
    forall faceIdx in frMesh.faceList.domain
    {
      // Check if the face's right neighbor is a Boundary Condition
      if frMesh.faceList[faceIdx].cells[2] < 0
      {
        // Yep, it is, lets get some local iteration variables
        ref faceFPini : int = frMesh.faceFPidx[faceIdx, 1];
        ref faceFPcnt : int = frMesh.faceFPidx[faceIdx, 2];

        ref thisBoco = frMesh.bocoList[-frMesh.faceList[faceIdx].cells[2]];
        ref thisFaml = frMesh.famlList[thisBoco.family];

        // Iterate through the FPs on this face
        for meshFP in faceFPini.. #faceFPcnt
        {
          // Calculate the boundary condition using the solution at the left neighbor´s corresponding FP
          frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml             ,
                                                          frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
        }
      }
    }

    // 9. Output initial solution
    iterOutput(lastIter, frMesh, flagNormals = true);

    // 10. Stabilize initial solution
    if Input.limiterScheme != LIMITER_NONE
    {
      // Loop through cells
      forall cellIdx in frMesh.cellList.domain
      {
        ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
        ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

        for varIdx in 1..frMesh.nVars
        {
          var stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                        jacobian = frMesh.jacSP[cellSPini.. #cellSPcnt]       ,
                                                        cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                        solDegree = frMesh.solOrder                           );

          if stableDegree < frMesh.solOrder then
            frMesh.solSP[varIdx, cellSPini.. #cellSPcnt] = projection_limiter(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                                              solDegree = frMesh.solOrder                           ,
                                                                              projDegree = stableDegree                             );
        }
      }
    }

    // 11. Save first restart file

    // 12. Initialize convergence monitoring variables
    var l1ResAbs : [1..frMesh.nVars] real;
    var l2ResAbs : [1..frMesh.nVars] real;
    var lfResAbs : [1..frMesh.nVars] real;
    var l1ResRel : [1..frMesh.nVars] real;
    var l2ResRel : [1..frMesh.nVars] real;
    var lfResRel : [1..frMesh.nVars] real;

    var convergenceLogFile : file;
    var     residueLogFile : file;
    var       errorLogFile : file;
    try! {
      convergenceLogFile = open("convergence.dat", ioMode.cw);
          residueLogFile = open(    "residue.dat", ioMode.cw);
            errorLogFile = open(      "error.dat", ioMode.cw);
    } catch {
      try! stdout.writeln("Unknown Error opening convergence log file.");
      try! stderr.writeln("Unknown Error opening convergence log file.");
    }
    var convergenceLogWriter = try! convergenceLogFile.writer();
    var     residueLogWriter = try!     residueLogFile.writer();
    var       errorLogWriter = try!       errorLogFile.writer();

    writeln();
    initTime = totalWatch.elapsed();
    writef("Initialization Time: %10.2dr ms\n", initTime*1000);

    ///////////////////////////
    // Main Solver Iteration //
    ///////////////////////////

    writeln();
    writef("Start Iterating\n");
    totalWatch.restart();

    label TIME_ITER for iteration in 1..Input.maxIter
    {
      iterWatch.restart();
      solveWatch.restart();

      // Save initial solution
      frMesh.oldSolSP = frMesh.solSP;
      oldSolTime += solveWatch.elapsed();

      // Iterate RK stages
      label RK_STAGE for rkStage in 1..timeStepStages
      {
        // Calculate residue for this iteration
        {
          solveWatch.restart();
          // The residue has 3 components:
          //   1. Continuous Flux
          //   2. Discontinuous Flux
          //   3. Source terms
          //
          // The residual array is reset in the time stepping procedure

          // Component 1: Source Term
          {
            residueWatch.restart();

            if Input.eqSet == EQ_QUASI_1D_EULER then
              forall spIdx in frMesh.resSP.domain.dim(1) do
                frMesh.resSP[.., spIdx.. #1] += -source_term(frMesh.xyzSP[spIdx.. #1, ..],
                                                             frMesh.solSP[.., spIdx.. #1],
                                                             Input.eqSet                 )
                                               * frMesh.jacSP[spIdx];

            srcTermTime += residueWatch.elapsed();
          }

          // Component 2: Discontinuous Flux
          {
            residueWatch.restart();

            // Calculate flux at SPs and it's divergence
            forall cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              const ref thisCell = frMesh.cellList[cellIdx];
              const     cellTopo = thisCell.elemTopo();
              const ref cellSPini = frMesh.cellSPidx[cellIdx, 1];
              const ref cellSPcnt = frMesh.cellSPidx[cellIdx, 2];

              // Allocate temporary flux array
              var flxSP : [1..frMesh.nDims, 1..frMesh.nVars, 1..cellSPcnt] real;

              // Step 1: Calculate fluxes at SPs
              //dscFluxWatch.restart();
              for meshSP in cellSPini.. #cellSPcnt do
                select Input.eqSet
                {
                  when EQ_CONVECTION do
                    flxSP[ 1, .., meshSP-cellSPini+1] = convection_flux_cv_1d(frMesh.solSP[.., meshSP], Input.convectionSpeed);
                  when EQ_INVBURGERS do
                    flxSP[ 1, .., meshSP-cellSPini+1] = burgers_flux_cv_1d(frMesh.solSP[.., meshSP]);
                  when EQ_QUASI_1D_EULER do
                    flxSP[ 1, .., meshSP-cellSPini+1] = euler_flux_cv_1d(frMesh.solSP[.., meshSP], Input.fGamma);
                  when EQ_EULER do
                    flxSP[.., .., meshSP-cellSPini+1] = euler_flux_cv(frMesh.solSP[.., meshSP], Input.fGamma);
                }
              //dscFluxTime1 += dscFluxWatch.elapsed();

              // Step 2: Interpolate fluxes to FPs and save the FP normal flux
              //dscFluxWatch.restart();
              for cellFace in thisCell.faces.domain
              {
                // Get loop variables
                const ref faceIdx  = thisCell.faces[cellFace];
                const ref thisFace = frMesh.faceList[faceIdx];
                const ref faceSide = thisCell.sides[cellFace];

                // Iterate though all FPs on this face
                for meshFP in frMesh.faceFPidx[faceIdx, 1].. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var cellFP : int;
                  if faceSide == 1 then
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                  else
                    cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  const uniNrm = frMesh.nrmFP[meshFP, ..]/norm(frMesh.nrmFP[meshFP, ..]);

                  frMesh.flxFP[meshFP, faceSide, ..] = 0;
                  for dimIdx in 1..frMesh.nDims do
                    for varIdx in 1..frMesh.nVars do
                      for cellSPidx in 1..cellSPcnt do
                        frMesh.flxFP[meshFP, faceSide, varIdx] += flxSP[ dimIdx, varIdx, cellSPidx]
                                                                 *sp2fpInterp[(cellTopo, frMesh.solOrder)]!.coefs[cellFP, cellSPidx]
                                                                 *uniNrm[dimIdx];
                }
              }
              //dscFluxTime2 += dscFluxWatch.elapsed();

              // Step 3: Convert fluxes from physical to computational domain
              //dscFluxWatch.restart();
              for cellSPidx in 1..cellSPcnt
              {
                const meshSPidx = cellSPini + cellSPidx - 1;

                // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                //var flxsp = flxSP[.., .., meshSP-cellSPini+1];
                //flxSP[.., .., cellSPidx] = dot( frMesh.metSP[meshSPidx, .., ..], flxsp)*frMesh.jacSP[meshSPidx];

                for varIdx in 1..frMesh.nVars
                {
                  var compFlxSP : [1..frMesh.nDims] real = 0;

                  for compDimIdx in 1..frMesh.nDims
                  {
                    // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                    for physDimIdx in 1..frMesh.nDims do
                      compFlxSP[compDimIdx] += flxSP[physDimIdx, varIdx, cellSPidx]
                                              *frMesh.metSP[meshSPidx, compDimIdx, physDimIdx];
                  }

                  flxSP[.., varIdx, cellSPidx] = compFlxSP * frMesh.jacSP[meshSPidx];
                }
              }
              //dscFluxTime3 += dscFluxWatch.elapsed();

              // Step 4: Calculate flux divergence
              //dscFluxWatch.restart();
              for cellSPidx in 1..cellSPcnt
              {
                const meshSPidx = cellSPini + cellSPidx - 1;

                for varIdx in 1..frMesh.nVars do
                  for dimIdx in 1..frMesh.nDims do
                    for spIdx in 1.. #cellSPcnt do
                      frMesh.resSP[varIdx, meshSPidx] += flxSP[dimIdx, varIdx, spIdx]
                                                        *sp2spDeriv[(cellTopo, frMesh.solOrder)]!.coefs[cellSPidx, dimIdx, spIdx];
              }
              //dscFluxTime4 += dscFluxWatch.elapsed();
            }

            dscFluxTime += residueWatch.elapsed();
          }

          // Component 3: Continuous Flux
          {
            residueWatch.restart();

            // Step 1: Interpolate solution to FPs
            cntFluxWatch.restart();
            forall cellIdx in frMesh.cellList.domain
            {
              const ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              const ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              const ref thisCell = frMesh.cellList[cellIdx];

              forall cellFace in thisCell.faces.domain
              {
                const ref faceIdx  : int = thisCell.faces[cellFace];
                const ref faceSide : int = thisCell.sides[cellFace];
                const ref thisFace = frMesh.faceList[faceIdx];

                forall meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  const cellFP : int = if faceSide == 1
                    then (cellFace-1)*(frMesh.solOrder+1) +  meshFP - frMesh.faceFPidx[faceIdx, 1] + 1
                    else (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] - (meshFP - frMesh.faceFPidx[faceIdx, 1]));

                  frMesh.solFP[meshFP, faceSide, ..] = 0;
                  for varIdx in 1..frMesh.nVars do
                    for spIdx in 1.. #cellSPcnt do
                      frMesh.solFP[meshFP, faceSide, varIdx] += frMesh.solSP[varIdx, cellSPini + spIdx - 1]
                                                               *sp2fpInterp[(thisCell.elemTopo(), frMesh.solOrder)]!.coefs[cellFP, spIdx];
                }
              }
            }
            cntFluxTime1 += cntFluxWatch.elapsed();

            // Step 2: Apply boundary conditions
            cntFluxWatch.restart();
            forall faceIdx in frMesh.faceList.domain
            {
              // Get loop variables
              ref faceFPini : int = frMesh.faceFPidx[faceIdx, 1];
              ref faceFPcnt : int = frMesh.faceFPidx[faceIdx, 2];
              ref thisFace = frMesh.faceList[faceIdx];

              // Check if the face´s right neighbor is a Boundary Condition
              if frMesh.faceList[faceIdx].cells[2] < 0
              {
                // Yep, it is, lets get some local iteration variables
                ref thisBoco = frMesh.bocoList[-frMesh.faceList[faceIdx].cells[2]];
                ref thisFaml = frMesh.famlList[thisBoco.family];

                // Iterate through the FPs on this face
                forall meshFP in faceFPini.. #faceFPcnt
                {
                  // Calculate the boundary condition using the solution at the left neighbor´s corresponding FP
                  frMesh.solFP[meshFP, 2, ..] = Boundary.boundary(frMesh.solFP[meshFP, 1, ..], thisFaml             ,
                                                                  frMesh.xyzFP[meshFP, ..], frMesh.nrmFP[meshFP, ..]);
                }
              }
            }
            cntFluxTime2 += cntFluxWatch.elapsed();

            // Step 3: Calculate interface correction
            cntFluxWatch.restart();
            forall cellIdx in frMesh.cellList.domain
            {
              // Get loop variables
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];
              ref thisCell = frMesh.cellList[cellIdx];
              const cellTopo : int = thisCell.elemTopo();

              for cellFace in thisCell.faces.domain
              {
                ref faceIdx  : int = thisCell.faces[cellFace];
                ref faceSide : int = thisCell.sides[cellFace];
                ref thisFace = frMesh.faceList[faceIdx];

                for meshFP in frMesh.faceFPidx[faceIdx, 1] .. #frMesh.faceFPidx[faceIdx, 2]
                {
                  var jump : [1..frMesh.nVars] real;

                  // Operation 1: Calculate Riemann flux at the FP
                  //correctionWatch.restart();
                  select Input.eqSet
                  {
                    when EQ_CONVECTION do
                      jump = upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_INVBURGERS do
                      jump = upwind_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_QUASI_1D_EULER do
                      jump = roe_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                    when EQ_EULER do
                      jump = roe(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..], frMesh.nrmFP[meshFP, ..]);
                  }
                  //riemTime += correctionWatch.elapsed();

                  // Operation 2: Calculate jump at a FP and convert it to the physical domain
                  //correctionWatch.restart();
                  {
                    // Calculate the flux jump = -1*(local_flux) + numerical_flux
                    jump -= frMesh.flxFP[meshFP, faceSide, ..];

                    // Convert fluxes from physical to computational domain.
                    // Multiply the flux vector by the inverse Jacobian matrix and by the Jacobian determinant
                    jump[..] = jump[..] * norm(frMesh.nrmFP[meshFP, ..], normType.norm2);

                    select cellTopo
                    {
                      when TOPO_LINE
                      {
                        if      faceSide == 1 && cellFace == 1 then
                          jump *= -1;
                        else if faceSide == 2 && cellFace == 2 then
                          jump *= -1;
                      }
                      when TOPO_QUAD
                      {
                        if      faceSide == 1 && (cellFace == 1 || cellFace == 4) then
                          jump *= -1;
                        else if faceSide == 2 && (cellFace == 2 || cellFace == 3) then
                          jump *= -1;
                      }
                    }
                  }
                  //jumpTime += correctionWatch.elapsed();

                  // Operation 3: Apply the correction to the residue matrix
                  //correctionWatch.restart();
                  {
                    // For 1D each face has 1 FP therefore the FP and the Face have the same index Relative to it's
                    // position in the cell
                    var cellFP : int;
                    const faceFP : int = meshFP - frMesh.faceFPidx[faceIdx, 1] + 1;
                    if faceSide == 1 then
                      cellFP = (cellFace-1)*(frMesh.solOrder+1) +  faceFP;
                    else
                      cellFP = (cellFace-1)*(frMesh.solOrder+1) + (frMesh.faceFPidx[faceIdx, 2] + 1 - faceFP);

                    // The correction function was calculated in the computational domain already, therefore no
                    // transformation is required.
                    frMesh.resSP[.., cellSPini.. #cellSPcnt] += outer(jump[..]                          ,
                        flux_correction[(cellTopo, frMesh.solOrder+1)]!.correction[cellFP, 1..cellSPcnt]);
                  }
                  //corrTime += correctionWatch.elapsed();
                }
              }
            }
            cntFluxTime3 += cntFluxWatch.elapsed();

            cntFluxTime += residueWatch.elapsed();
          }

          residueTime += solveWatch.elapsed();
        }

        // Advance RK Stage
        {
          solveWatch.restart();

          if rkStage == 1 then frMesh.calc_time_step();

          // Loop through cells
          forall cellIdx in frMesh.cellList.domain
          {
            ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
            ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

            // Convert the residual from the computational to the physical domain
            forall meshSP in cellSPini.. #cellSPcnt do
              frMesh.resSP[.., meshSP] /= frMesh.jacSP[meshSP];

            // Update solution
            frMesh.solSP[.., cellSPini.. #cellSPcnt] = time_advance(frMesh.oldSolSP[.., cellSPini.. #cellSPcnt],
                                                                    frMesh.solSP[   .., cellSPini.. #cellSPcnt],
                                                                    frMesh.resSP[   .., cellSPini.. #cellSPcnt],
                                                                    frMesh.cellTimeStep[cellIdx],
                                                                    rkStage, Input.timeScheme                   );
          }

          // Get the residue norms at the first stage of the RK Iteration
          if rkStage == 1
          {
            l1ResAbs = [varIdx in 1..frMesh.nVars]
              + reduce abs( frMesh.resSP[varIdx, ..])
              /frMesh.resSP.domain.dim(1).size;
            l1ResRel = [varIdx in 1..frMesh.nVars]
              + reduce abs( frMesh.resSP[varIdx, ..] / frMesh.oldSolSP[varIdx, ..])
              /frMesh.resSP.domain.dim(1).size;

            l2ResAbs = [varIdx in 1..frMesh.nVars]
              sqrt(   + reduce    ( frMesh.resSP[varIdx, ..]**2))
              /frMesh.resSP.domain.dim(1).size;
            l2ResRel = [varIdx in 1..frMesh.nVars]
              sqrt(   + reduce    ((frMesh.resSP[varIdx, ..] / frMesh.oldSolSP[varIdx, ..])**2))
              /frMesh.resSP.domain.dim(1).size;

            lfResAbs = [varIdx in 1..frMesh.nVars] max reduce abs( frMesh.resSP[varIdx, ..]);
            lfResRel = [varIdx in 1..frMesh.nVars] max reduce abs( frMesh.resSP[varIdx, ..] / frMesh.oldSolSP[varIdx, ..]);

            // Output all residues to log file
            log_convergence(residueLogWriter, resToggle=true, iteration, l1ResAbs, l2ResAbs, lfResAbs, l1ResRel, l2ResRel, lfResRel);
          }

          // Zero out residue
          frMesh.resSP = 0.0;

          timeStepTime += solveWatch.elapsed();
        }

        // Stabilize Solution
        {
          solveWatch.restart();

          if Input.limiterScheme != LIMITER_NONE
          {
            // Loop through cells
            forall cellIdx in frMesh.cellList.domain
            {
              ref cellSPini : int = frMesh.cellSPidx[cellIdx, 1];
              ref cellSPcnt : int = frMesh.cellSPidx[cellIdx, 2];

              for varIdx in 1..frMesh.nVars
              {
                const stableDegree : int = troubled_cell_marker(solPoly = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                              jacobian = frMesh.jacSP[cellSPini.. #cellSPcnt]       ,
                                                              cellTopo = frMesh.cellList[cellIdx].elemTopo()        ,
                                                              solDegree = frMesh.solOrder                           );

                if stableDegree < frMesh.solOrder then
                  frMesh.solSP[varIdx, cellSPini.. #cellSPcnt] = projection_limiter(solPoly    = frMesh.solSP[varIdx, cellSPini.. #cellSPcnt],
                                                                                    cellTopo   = frMesh.cellList[cellIdx].elemTopo()         ,
                                                                                    solDegree  = frMesh.solOrder                             ,
                                                                                    projDegree = stableDegree                                );
              }
            }
          }

          stabilizeTime += solveWatch.elapsed();
        }
      }

      // IO
      {
        solveWatch.restart();

        // Save restart file

        // Check if we should write the solution this iteration
        if (ioIter > 0 && iteration % ioIter == 0) then
            iterOutput(iteration, frMesh);

        // Calculate solution error and write to error log
        if Input.outError > 0
        {
          const errors = error_calc(Input.outError, frMesh);
          print_log(errorLogWriter, iteration, errors);
        }

        // Calculate and print convergence metrics
        {
          // Calculate solution delta from previous iteration
          const l1SolDeltaAbs : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars]
              + reduce abs(frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..])
            /frMesh.solSP.domain.dim(1).size;
          const l1SolDeltaRel : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars]
              + reduce abs((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..])
            /frMesh.solSP.domain.dim(1).size;

          const l2SolDeltaAbs : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars] sqrt(
              + reduce    (frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..])**2)
            /frMesh.solSP.domain.dim(1).size;
          const l2SolDeltaRel : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars] sqrt(
              + reduce    ((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..] )**2)
            /frMesh.solSP.domain.dim(1).size;

          const lfSolDeltaAbs : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars]
              max reduce abs(frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]);
          const lfSolDeltaRel : [1..frMesh.nVars] real =
            [varIdx in 1..frMesh.nVars]
              max reduce abs((frMesh.solSP[varIdx, ..] - frMesh.oldSolSP[varIdx, ..]) / frMesh.oldSolSP[varIdx, ..]);

          // Output full state to log file
          log_convergence(convergenceLogWriter, resToggle=false, iteration,
                          l1SolDeltaAbs, l2SolDeltaAbs, lfSolDeltaAbs,
                          l1SolDeltaRel, l2SolDeltaRel, lfSolDeltaRel);

          // Output summarized convergence metrics to stdOut
          writef("Iteration %9i | Time %{ 10.2dr}ms | Log10(L2(Res)) = %{ 8.4dr} | Log10(L2(ΔSol)) = %{ 8.4dr}",
              iteration, iterWatch.elapsed()*1000, log10(norm(l2ResAbs)), log10(norm(l2SolDeltaAbs)));

          if (ioIter > 0 && iteration % ioIter == 0) then
            writef(" | Solution file saved\n");
          else
            writef("\n");

          // Iteration stop criteria
          if !( && reduce isFinite(frMesh.solSP) ) // Solution Diverged
          {
            writef("\nSolver diverged, restoring previous solution\n");

            // Recover last stable solution and let the finalization output print the pre-divergence solution
            frMesh.solSP = frMesh.oldSolSP;

            lastIter = iteration-1;
            break TIME_ITER;
          }
          else if (    norm(l2ResAbs     ) <= Input.l2ResStop
                    || norm(l2SolDeltaAbs) <= Input.l2SolStop ) // Solution converged
          {
            writef("\nSolver reached convergence criteria\n");

            // Export iteration count to outside the loop
            lastIter = iteration;
            break TIME_ITER;
          }
          else if (iteration == Input.maxIter) // Maximum iteration count reached
          {
            writef("\nSolver reached maxIter\n");

            // If this is the last iteration then export the count to outside the loop
            lastIter = iteration;
          }
        }

        // Check if input file changed
        {}

        ioIterTime += solveWatch.elapsed();
      }
    }

    solveTime = totalWatch.elapsed();

    /////////////////////////////
    // Output and Finalization //
    /////////////////////////////

    // Output the final solution
    totalWatch.restart();
    if (lastIter > 0) then iterOutput(lastIter, frMesh);
    outputTime = totalWatch.elapsed();

    // Finalize program stopwatch and calculate agregate times for major program steps
    programTime = programWatch.elapsed();

    writeln();
    writef("Time splits:\n");
    writef("- Init        : %11.2dr ms - %4.1dr%% of Run-Time\n",       initTime*1000,      initTime/programTime  *100);
    writef("- Iterations  : %11.2dr ms - %4.1dr%% of Run-Time\n",      solveTime*1000,     solveTime/programTime  *100);
    writef("  - Residue   : %11.2dr ms - %4.1dr%% of Iteration\n",   residueTime*1000,   residueTime/solveTime    *100);
    writef("    - Src Term: %11.2dr ms - %4.1dr%% of Residue\n",     srcTermTime*1000,   srcTermTime/residueTime  *100);
    writef("    - Dsc Flux: %11.2dr ms - %4.1dr%% of Residue\n",     dscFluxTime*1000,   dscFluxTime/residueTime  *100);
  //writef("      - Step 1: %11.2dr ms - %4.1dr%% of Dsc Flux\n",   dscFluxTime1*1000,  dscFluxTime1/dscFluxTime  *100);
  //writef("      - Step 2: %11.2dr ms - %4.1dr%% of Dsc Flux\n",   dscFluxTime2*1000,  dscFluxTime2/dscFluxTime  *100);
  //writef("      - Step 3: %11.2dr ms - %4.1dr%% of Dsc Flux\n",   dscFluxTime3*1000,  dscFluxTime3/dscFluxTime  *100);
  //writef("      - Step 4: %11.2dr ms - %4.1dr%% of Dsc Flux\n",   dscFluxTime4*1000,  dscFluxTime4/dscFluxTime  *100);
    writef("    - Cnt Flux: %11.2dr ms - %4.1dr%% of Residue\n",     cntFluxTime*1000,   cntFluxTime/residueTime  *100);
    writef("      - Step 1: %11.2dr ms - %4.1dr%% of Cnt Flux\n",   cntFluxTime1*1000,  cntFluxTime1/cntFluxTime  *100);
    writef("      - Step 2: %11.2dr ms - %4.1dr%% of Cnt Flux\n",   cntFluxTime2*1000,  cntFluxTime2/cntFluxTime  *100);
    writef("      - Step 3: %11.2dr ms - %4.1dr%% of Cnt Flux\n",   cntFluxTime3*1000,  cntFluxTime3/cntFluxTime  *100);
  //writef("        - Op 1: %11.2dr ms - %4.1dr%% of St3 Flux\n",       riemTime*1000,      riemTime/cntFluxTime3 *100);
  //writef("        - Op 2: %11.2dr ms - %4.1dr%% of St3 Flux\n",       jumpTime*1000,      jumpTime/cntFluxTime3 *100);
  //writef("        - Op 3: %11.2dr ms - %4.1dr%% of St3 Flux\n",       corrTime*1000,      corrTime/cntFluxTime3 *100);
    writef("  - Time-Step : %11.2dr ms - %4.1dr%% of Iteration\n",  timeStepTime*1000,  timeStepTime/solveTime    *100);
    writef("  - Stabilize : %11.2dr ms - %4.1dr%% of Iteration\n", stabilizeTime*1000, stabilizeTime/solveTime    *100);
    writef("  - Iter IO   : %11.2dr ms - %4.1dr%% of Iteration\n",    ioIterTime*1000,    ioIterTime/solveTime    *100);
    writef("- Output      : %11.2dr ms - %4.1dr%% of Run-Time\n",     outputTime*1000,    outputTime/programTime  *100);
    writef("---------------------------------------------------------\n");
    writef("  Run-time    : %11.2dr ms\n", programTime*1000);

  //var  dscFluxSumTime : real = dscFluxTime1 + dscFluxTime2 + dscFluxTime3 + dscFluxTime4;
    var  cntFluxSumTime : real = cntFluxTime1 + cntFluxTime2 + cntFluxTime3 ;
  //var cntStep3SumTime : real = riemTime + jumpTime + corrTime;
    var  residueSumTime : real = srcTermTime + dscFluxTime + cntFluxTime;
    var     iterSumTime : real = residueTime + stabilizeTime + timeStepTime + ioIterTime;
    var  programSumTime : real = initTime + solveTime + outputTime;

    writeln();
    writef("Split verifications:\n");
  //writef("  Dsc Flux  Sum: %11.2dr ms - %4.1dr%% of Dsc Flux\n",   dscFluxSumTime*1000,  dscFluxSumTime/dscFluxTime *100);
    writef("  Cnt Flux  Sum: %11.2dr ms - %4.1dr%% of Cnt Flux\n",   cntFluxSumTime*1000,  cntFluxSumTime/cntFluxTime *100);
  //writef("  Cnt Step3 Sum: %11.2dr ms - %4.1dr%% of Cnt Step3\n", cntStep3SumTime*1000, cntStep3SumTime/cntFluxTime3*100);
    writef("  Residue   Sum: %11.2dr ms - %4.1dr%% of Residue\n",    residueSumTime*1000,  residueSumTime/residueTime *100);
    writef("  Iter      Sum: %11.2dr ms - %4.1dr%% of Solve\n",         iterSumTime*1000,     iterSumTime/solveTime   *100);
    writef("  Program   Sum: %11.2dr ms - %4.1dr%% of Program\n",    programSumTime*1000,  programSumTime/programTime *100);

    writeln();
    writeln("Fin");
  }
}
