/* Documentation for FREI */
module FREI
{
  //Runtime constants
  config const inputFile : string = "input.toml";

  proc main() {
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
    use FR;
    use LinearAlgebra;
    use SourceTerm;

    var iteration : int = 0;

    // 1. Read input data
    indat(inputFile);

    // 2. Process input data and configure program
    //configure();

    // 3. Read / define mesh
    var gmesh2 = new unmanaged gmesh2_c(nNodes=Input.nCells+1, nElements=Input.nCells+2, nFamilies=Input.nFaml);
    gmesh2.random1D(Input.xMin, Input.xMax);

    // 5. Convert input mesh to solver mesh
    var frMesh = new unmanaged fr_mesh_c(nDims=1, nVars=3, solOrder=Input.iOrder-1);
    frMesh.import_gmesh2(gmesh2);   // Convert mesh to native format
    frMesh.set_families(famlList);  // Get families data from input file and write to mesh
    frMesh.allocate_fr_vars();      // Allocate SP and FP solution/flux/residue arrays
    frMesh.set_points_locations();  // Calculate coordinate trasnformations and point coordinates

    // 4. Initialize the solver, pre calculate stuff
    init_sp2fpInterp(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_sp2spDeriv(Input.minOrder, Input.maxOrder, frMesh.cellTopos);
    init_correction();

    // Save mesh file in internal format

    // Initialize solution
    frMesh.solSP = initial_condition(Input.initialCondition, frMesh.xyzSP);

    // Save restart file

    // Output initial state
    iterOutput(iteration, frMesh);

    writeln("Start Iterating");

    // Solve flow
    for iteration in 1..Input.maxIter
    {
      // Save initial solution
      frMesh.oldSolSP = frMesh.solSP;

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
          frMesh.resSP = -source_term(frMesh.xyzSP, frMesh.solSP, Input.eqSet) * frMesh.jacSP;
        }

        // Component 2: Discontinuous Flux
        {
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
              flxSP[meshSP, ..] = invs_flux_cv_1d(frMesh.solSP[meshSP, ..]);

            // Convert fluxes from physical to computational domain
            flxSP[cellSPini.. #cellSPcnt, ..] *= frMesh.metSP[cellSPini.. #cellSPcnt, 1, 1];

            // Calculate flux divergence
            for cellSP in 1..cellSPcnt
            {
              var meshSP = cellSPini + cellSP - 1;
              frMesh.resSP[meshSP,..] += dot(sp2spDeriv[(thisCell.elemTopo(), iOrder)]!.coefs(cellSP, ..),
                                             flxSP[cellSPini..#cellSPcnt,..]                             );
            }
          }
        }

        // Component 3: Continuous Flux
        {
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

          // Calculate numerical flux at interfaces
          for meshFP in frMesh.flxFP.domain.dim(0)
          {
            // Riemann flux
            frMesh.flxFP[meshFP, ..] = rusanov_1d(frMesh.solFP[meshFP, 1, ..], frMesh.solFP[meshFP, 2, ..]);
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
                // For 1D each face has 1 FP therefore the FP and the Face have the same index relative to it's
                // position in the cell
                var cellFP = cellFace;

                // Calculate the flux jump = numerical_flux - local_flux
                var jump : [1..frMesh.nVars] real = frMesh.flxFP[meshFP, ..] - invs_flux_cv_1d(frMesh.solFP[meshFP, faceSide, ..]);

                // Convert flux jump from physical to computational domains
                jump *= frMesh.metFP[meshFP, faceSide, 1, 1];

                frMesh.resSP[cellSPini.. #cellSPcnt, ..] += outer(jump[..],
                    flux_correction[thisCell.elemTopo(), iOrder]!.correction[1..cellSPcnt, cellFP]);
              }
            }
          }
        }
      }

      // Print solver status / log

      // Save restart file

      // Output solution state
      if iteration % ioIter == 0 then
        iterOutput(iteration, frMesh);

      // Check input changes
    }

    // Output the final solution
    //iterOutput(iteration, frMesh);

    writeln();
    writeln("Fin.");
  }
}
