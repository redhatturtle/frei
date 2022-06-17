module Gmesh
{
  use Time;
  use Random;
  use UnitTest;

  config const testMesh : string = "test-cases/gmesh-1d-test.msh";

  // Planned supported mesh formats
  //
  // - 1D Internal generated mesh  (Frei 1.0)
  // - 1D/2D Gmesh version 2/4     (Frei 2.0)
  // - 2D Internal generated Gmesh (Frei 2.x)
  // - 3D Gmesh version 2/4        (Frei 3.0)
  // - 3D CGNS                     (Frei 4.0)
  // - 2D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated CGNS  (Frei 4.x)
  // - 3D Internal generated Gmesh (Frei 4.x)
  //
  // It would probably be better to have a separate module for each IO standard (CGNS, Gmesh, Plot3D, etc). At the
  //    moment shit's just pilled up here for compactness since most of it is just general concepts and sketches.
  //
  // The idea behind having a separate class for each IO format it to isolate the conversion from the solver mesh format
  //    to the IO format from the reading and writing of the IO files. It might waste to much memory for just code
  //    organization's sake and should be reviewed in the future.

  record gmesh_element_r
  {
    var elemType : int;
    var tags_d   : domain(rank=1, idxType=int);
    var tags     : [tags_d] int;
    var nodes_d  : domain(rank=1, idxType=int);
    var nodes    : [nodes_d] int;

    proc elemDim() : int
    {
      return gmsh_elem_dimension(this.elemType);
    }

    proc setNodes()
    {
      use Parameters.ParamGmesh;

      select this.elemType {
        // Point
        when GMESH_PNT_1    do this.nodes_d = {1..1   };
        // Line
        when GMESH_LIN_2    do this.nodes_d = {1..2   };
        when GMESH_LIN_3    do this.nodes_d = {1..3   };
        when GMESH_LIN_4    do this.nodes_d = {1..4   };
        when GMESH_LIN_5    do this.nodes_d = {1..5   };
        when GMESH_LIN_6    do this.nodes_d = {1..6   };
        when GMESH_LIN_7    do this.nodes_d = {1..7   };
        when GMESH_LIN_8    do this.nodes_d = {1..8   };
        when GMESH_LIN_9    do this.nodes_d = {1..9   };
        when GMESH_LIN_10   do this.nodes_d = {1..10  };
        // Lagrange Triangles
        when GMESH_TRI_3    do this.nodes_d = {1..3   };
        when GMESH_TRI_6    do this.nodes_d = {1..6   };
        when GMESH_TRI_10   do this.nodes_d = {1..10  };
        when GMESH_TRI_15   do this.nodes_d = {1..15  };
        when GMESH_TRI_21   do this.nodes_d = {1..21  };
        when GMESH_TRI_28   do this.nodes_d = {1..28  };
        when GMESH_TRI_36   do this.nodes_d = {1..36  };
        when GMESH_TRI_45   do this.nodes_d = {1..45  };
        when GMESH_TRI_55   do this.nodes_d = {1..55  };
        // Lagrange Quadrilaterals
        when GMESH_QUA_4    do this.nodes_d = {1..4   };
        when GMESH_QUA_9    do this.nodes_d = {1..9   };
        when GMESH_QUA_16   do this.nodes_d = {1..16  };
        when GMESH_QUA_25   do this.nodes_d = {1..25  };
        when GMESH_QUA_36   do this.nodes_d = {1..36  };
        when GMESH_QUA_49   do this.nodes_d = {1..49  };
        when GMESH_QUA_64   do this.nodes_d = {1..64  };
        when GMESH_QUA_81   do this.nodes_d = {1..81  };
        when GMESH_QUA_100  do this.nodes_d = {1..100 };
        // Lagrange Tetrahedra
        when GMESH_TET_4    do this.nodes_d = {1..4   };
        when GMESH_TET_10   do this.nodes_d = {1..10  };
        when GMESH_TET_20   do this.nodes_d = {1..20  };
        when GMESH_TET_35   do this.nodes_d = {1..35  };
        when GMESH_TET_56   do this.nodes_d = {1..56  };
        when GMESH_TET_84   do this.nodes_d = {1..84  };
        when GMESH_TET_120  do this.nodes_d = {1..120 };
        when GMESH_TET_165  do this.nodes_d = {1..165 };
        when GMESH_TET_220  do this.nodes_d = {1..220 };
        // Lagrange Pyramids
        when GMESH_PYR_5    do this.nodes_d = {1..5   };
        when GMESH_PYR_14   do this.nodes_d = {1..14  };
        when GMESH_PYR_30   do this.nodes_d = {1..30  };
        when GMESH_PYR_55   do this.nodes_d = {1..55  };
        when GMESH_PYR_91   do this.nodes_d = {1..91  };
        when GMESH_PYR_140  do this.nodes_d = {1..140 };
        when GMESH_PYR_204  do this.nodes_d = {1..204 };
        when GMESH_PYR_285  do this.nodes_d = {1..285 };
        when GMESH_PYR_385  do this.nodes_d = {1..385 };
        // Lagrange Prisms
        when GMESH_PRI_6    do this.nodes_d = {1..6   };
        when GMESH_PRI_18   do this.nodes_d = {1..18  };
        when GMESH_PRI_40   do this.nodes_d = {1..40  };
        when GMESH_PRI_75   do this.nodes_d = {1..75  };
        when GMESH_PRI_126  do this.nodes_d = {1..126 };
        when GMESH_PRI_196  do this.nodes_d = {1..196 };
        when GMESH_PRI_288  do this.nodes_d = {1..288 };
        when GMESH_PRI_405  do this.nodes_d = {1..405 };
        when GMESH_PRI_550  do this.nodes_d = {1..550 };
        // Lagrange Hexahedra
        when GMESH_HEX_8    do this.nodes_d = {1..8   };
        when GMESH_HEX_27   do this.nodes_d = {1..27  };
        when GMESH_HEX_64   do this.nodes_d = {1..64  };
        when GMESH_HEX_125  do this.nodes_d = {1..125 };
        when GMESH_HEX_216  do this.nodes_d = {1..216 };
        when GMESH_HEX_343  do this.nodes_d = {1..343 };
        when GMESH_HEX_512  do this.nodes_d = {1..512 };
        when GMESH_HEX_729  do this.nodes_d = {1..729 };
        when GMESH_HEX_1000 do this.nodes_d = {1..1000};
      }
    }
  }

  record gmesh_family_r
  {
    var tag  : int;
    var nDim : int;
    var name : string;
  }

  class gmesh2_c
  {
    var nodes_d    : domain(rank=2, idxType=int);
    var elements_d : domain(rank=1, idxType=int);
    var families_d : domain(rank=1, idxType=int);

    var nodes    : [nodes_d] real;
    var elements : [elements_d] gmesh_element_r;
    var families : [families_d] gmesh_family_r;

    proc random1D(nCells : int, xMin : real = -1.0, xMax : real= 1.0)
    {
      use Random;
      use Parameters.ParamTest;
      use Parameters.ParamGmesh;
      import Sort.sort;

      var randStreamSeeded = new RandomStream(real, RANDOM_SEED);

      // Allocate mesh elements
      this.nodes_d = {1..nCells+1, 1..3};
      this.elements_d = {1..nCells+2};
      this.families_d = {1..3};

      var nodePermutation : [this.nodes.domain.dim(0)] int;
      var elemPermutation : [this.elements.domain] int;
      var x = this.nodes[..,1];
      var cells : [1..nCells, 1..2] int;

      // Get random values from xMin to xMax
      x[1] = 0;
      x[2] = 1;
      randStreamSeeded.fillRandom(x[3..]);
      x = x*(xMax-xMin) + xMin;

      // Sort nodes to build elements and generate a random permutation
      sort(x);
      permutation(nodePermutation, RANDOM_SEED);
      permutation(elemPermutation, RANDOM_SEED);

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..nCells
      {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].tag  = 1;
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].tag  = 2;
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].tag  = 3;
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list in randomized order
      for i in this.nodes.domain.dim(0) do
        this.nodes[nodePermutation[i],1] = x[i];

      // Fill element list with the internal elements in random order
      for i in this.elements.domain.dim(0).expand(-1)
      {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes();
        this.elements[elemPermutation[i]].tags_d = {1..1};
        this.elements[elemPermutation[i]].tags[1] = 1; // Family
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[cells[i-1,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[cells[i-1,2]];
      }

      // Add left boundary point to elements list
      this.elements[elemPermutation[1]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[1]].setNodes();
      this.elements[elemPermutation[1]].tags_d = {1..1};
      this.elements[elemPermutation[1]].tags[1] = 2; // Family
      this.elements[elemPermutation[1]].nodes[1] = nodePermutation[1];

      // Add right boundary point to elements list
      this.elements[elemPermutation[nCells+2]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[nCells+2]].setNodes();
      this.elements[elemPermutation[nCells+2]].tags_d = {1..1};
      this.elements[elemPermutation[nCells+2]].tags[1] = 3; // Family
      this.elements[elemPermutation[nCells+2]].nodes[1] = nodePermutation[this.nodes.domain.dim(0).high];
    }

    proc uniform1D(nCells : int, xMin: real = -1.0, xMax: real = 1.0)
    {
      use Parameters.ParamGmesh;
      import Sort.sort;

      // Allocate mesh elements
      this.nodes_d = {1..nCells+1, 1..3};
      this.elements_d = {1..nCells+2};
      this.families_d = {1..3};

      var x = this.nodes[..,1];
      var cells : [1..nCells, 1..2] int;

      // Get uniform nodes from xMin to xMax
      var step : real = (xMax - xMin)/nCells;
      for node in this.nodes.domain.dim(0) do
        x[node] = xMin + (node-1)*step;

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..nCells
      {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].tag  = 1;
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].tag  = 2;
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].tag  = 3;
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list
      for i in this.nodes.domain.dim(0) do
        this.nodes[i,1] = x[i];

      // Fill element list with the internal elements
      for i in this.elements.domain.dim(0).expand(-1)
      {
        this.elements[i].elemType = GMESH_LIN_2;
        this.elements[i].setNodes();
        this.elements[i].tags_d = {1..1};
        this.elements[i].tags[1] = 1; // Family
        this.elements[i].nodes[1] = cells[i-1,1];
        this.elements[i].nodes[2] = cells[i-1,2];
      }

      // Add left boundary point to elements list
      this.elements[1].elemType = GMESH_PNT_1;
      this.elements[1].setNodes();
      this.elements[1].tags_d = {1..1};
      this.elements[1].tags[1] = 2; // Family
      this.elements[1].nodes[1] = 1;

      // Add right boundary point to elements list
      this.elements[nCells+2].elemType = GMESH_PNT_1;
      this.elements[nCells+2].setNodes();
      this.elements[nCells+2].tags_d = {1..1};
      this.elements[nCells+2].tags[1] = 3; // Family
      this.elements[nCells+2].nodes[1] = this.nodes.domain.dim(0).high;
    }

    proc read_gmesh_file(meshFileName : string)
    {
      use IO;
      use SysError;

      enum section {Main, MeshFormat, PhysicalNames, Nodes, Elements, Periodic, NodeData, ElementData,
                    ElementNodeData, InterpolationScheme};

      var meshFile : file;
      var procWatch    : Timer;
      var sectionWatch : Timer;

      // Open the mesh file for reading only
      try {
        writeln();
        procWatch.start();
        sectionWatch.start();
        writeln("Opening Gmesh2 formated mesh file");
        meshFile = open(meshFileName, iomode.r);
      } catch e : FileNotFoundError {
        writeln("Critical Error: Mesh file not found.");
      } catch {
        writeln("Unknown Error opening mesh file.");
      }

      // Open reader channel to gmesh file
      var meshReaderChannel = try! meshFile.reader();

      // Set initial state
      var state = section.Main;
      var nLines : int = 0;
      var lineIdx : int = 1;

      // Start reading and parsing file
      for line in meshReaderChannel.lines()
      {
        select state
        {
          when section.Main
          {
            // Update state with the section starting at this line
            select line
            {
              when "$MeshFormat\n" do
              {
                sectionWatch.clear();
                write("    Verifying Mesh Format  ...");
                state = section.MeshFormat;
              }
              when "$PhysicalNames\n"
              {
                sectionWatch.clear();
                write("    Reading Physical Names ...");
                state = section.PhysicalNames;
              }
              when "$Nodes\n"
              {
                sectionWatch.clear();
                write("    Reading Mesh Nodes     ...");
                state = section.Nodes;
              }
              when "$Elements\n"
              {
                sectionWatch.clear();
                write("    Reading Mesh Elements  ...");
                state = section.Elements;
              }
              when "$Periodic\n"
              {
                sectionWatch.clear();
                write("    Reading Periodic BCs   ...");
                state = section.Periodic;
              }
              when "$NodeData\n"
              {
                sectionWatch.clear();
                write("    Reading Node Data      ...");
                state = section.NodeData;
              }
              when "$ElementData\n"
              {
                sectionWatch.clear();
                write("    Reading Elem Data      ...");
                state = section.ElementData;
              }
              when "$ElementNodeData\n"
              {
                sectionWatch.clear();
                write("    Reading Elem Node Data ...");
                state = section.ElementNodeData;
              }
              when "$InterpolationScheme\n"
              {
                sectionWatch.clear();
                write("    Reading Interp Scheme  ...");
                state = section.InterpolationScheme;
              }
              otherwise do
                writeln("Unexpected line on main section");
            }
          }
          when section.MeshFormat
          {
            if line == "$EndMeshFormat\n"
            {
                // Reset state to main section
                writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
                state = section.Main;
                nLines = 0;
            }
            else if line != "2.2 0 8\n"
            {
              // Validate if gmesh format is compatible with this class
              writeln("Unsuported gmesh format version");
              writeln("Expected '2.2 0 8'");
            }
          }
          when section.PhysicalNames
         {
            if line == "$EndPhysicalNames\n"
            {
                // Reset state to main section
                writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Physical Names Count` to int");
              }
              this.families_d = {1..nLines};
            }
            else
            {
              var physicalDim  : int;
              var physicalTag  : int;
              var physicalName : string;
              var valueIdx : int = 1;

              for value in line.split()
              {
                 if valueIdx == 1 then
                    try! this.families[lineIdx].nDim = value:int;
                 if valueIdx == 2 then
                    try! this.families[lineIdx].tag = value:int;
                 if valueIdx == 3 then
                    try! this.families[lineIdx].name = value.strip("\"");

                 valueIdx += 1;
              }

              lineIdx += 1;
            }
          }
          when section.Nodes
          {
            if line == "$EndNodes\n"
            {
                // Reset state to main section
                writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Node Count` to int");
              }
              this.nodes_d = {1..nLines, 1..3};
            }
            else
            {
              var nodeIdx   : int;
              var nodeCoord : [1..3] real;

              var valueIdx  : int = 1;
              for value in line.split()
              {
                 if valueIdx == 1 then
                    try! nodeIdx = value : int;
                 if valueIdx == 2 then
                    try! nodeCoord[1] = value : real;
                 if valueIdx == 3 then
                    try! nodeCoord[2] = value : real;
                 if valueIdx == 4 then
                    try! nodeCoord[3] = value : real;

                 valueIdx += 1;
              }

              this.nodes[nodeIdx, 1..3] = nodeCoord[1..3];

              lineIdx += 1;
            }
          }
          when section.Elements
          {
            if line == "$EndElements\n"
            {
                // Reset state to main section
                writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
                state = section.Main;
                nLines = 0;
                lineIdx = 1;
            }
            else if nLines == 0
            {
              // If it's the first line of the section get the number of Physical Namas and allocate families
              try! {
                nLines = line : int;
              }
              catch {
                writeln("Failed to cast `Number of Elements` to int");
              }
              this.elements_d = {1..nLines};
            }
            else
            {
              var elemIdx   : int;

              var valueIdx  : int = 1;
              for value in line.split()
              {
                 // Element Number
                 if valueIdx == 1 then
                    try! elemIdx = value : int;
                 // Element Type
                 if valueIdx == 2
                 {
                    try! this.elements[elemIdx].elemType = value : int;
                    this.elements[elemIdx].setNodes();
                 }
                 // Number of Tags
                 if valueIdx == 3 then
                    try! this.elements[elemIdx].tags_d = {1..value:int};
                 // Tags
                 if ((valueIdx > 3) && (valueIdx < 4+this.elements[elemIdx].tags.domain.dim(0).high)) then
                    try! this.elements[elemIdx].tags[valueIdx-3] = value : int;
                 // Node List
                 if valueIdx >= 4+this.elements[elemIdx].tags.domain.dim(0).high then
                   try! this.elements[elemIdx].nodes[valueIdx-(3+this.elements[elemIdx].tags_d.dim(0).high)] = value : int;

                 valueIdx += 1;
              }

              lineIdx += 1;
            }
          }
          when section.Periodic
          {
            if line == "$EndPeriodic\n"
            {
              writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
              state = section.Main;
            }
          }
          when section.NodeData
          {
            if line == "$EndNodeData\n"
            {
              writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
              state = section.Main;
            }
          }
          when section.ElementData
          {
            if line == "$EndElementData\n"
            {
              writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
              state = section.Main;
            }
          }
          when section.ElementNodeData
          {
            if line == "$EndElementNodeData\n"
            {
              writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
              state = section.Main;
            }
          }
          when section.InterpolationScheme
          {
            if line == "$EndInterpolationScheme\n"
            {
              writef(" done in %7.1dr ms\n", sectionWatch.elapsed(TimeUnits.milliseconds));
              state = section.Main;
            }
          }
        }
      }

      writeln("    Families found:");
      for family in this.families do
        writef("        %1i-D Elements, Tag: %2i, Name: %s\n", family.nDim, family.tag, family.name);
      writef("    Reading Gmsh2 file done in %7.1dr ms\n", procWatch.elapsed(TimeUnits.milliseconds));
      writeln();
    }

    proc write_gmesh_file() {}

    proc mesh_dimension()
    {
      var meshDim : int = 0;

      for family in this.families do
        if family.nDim > meshDim then
          meshDim = family.nDim;

      return meshDim;
    }
  }

  class gmesh4_c
  {
    var nNodes : int;
    var nElements : int;
    var nFamilies : int;

    var nodes    : [1..nNodes, 1..3] real;
    var elements : [1..nElements] gmesh_element_r;
    var families : [1..nFamilies] gmesh_family_r;

    proc random1D(xMin: real, xMax: real)
    {
      use Random;
      use Parameters.ParamTest;
      use Parameters.ParamGmesh;
      import Sort.sort;

      var randStreamSeeded = new RandomStream(real, RANDOM_SEED);
      var x = this.nodes[..,1];
      var cells : [1..this.nElements-2, 1..2] int;
      var nodePermutation : [1..this.nNodes] int;
      var elemPermutation : [1..this.nElements] int;

      // Get random values from xMin to xMax
      x[1] = 0;
      x[2] = 1;
      randStreamSeeded.fillRandom(x[3..]);
      x = x*(xMax-xMin) + xMin;

      // Sort nodes to build elements and generate a random permutation
      sort(x);
      permutation(nodePermutation, RANDOM_SEED);
      permutation(elemPermutation, RANDOM_SEED);

      // Fill element list with non overlapping elements oriented from left to right
      for i in 1..this.nElements-2 {
        cells[i,1] = i;
        cells[i,2] = i+1;
      }

      // Commit values to object

      // Set the boundary and internal families
      this.families[1].nDim = 1;
      this.families[1].name = "flow";
      this.families[2].nDim = 0;
      this.families[2].name = "left";
      this.families[3].nDim = 0;
      this.families[3].name = "right";

      // Fill node list in randomized order
      for i in 1..this.nNodes do
        this.nodes[nodePermutation[i],1] = x[i];

      // Fill element list with the internal elements in random order
      for i in 2..this.nElements-1 {
        this.elements[elemPermutation[i]].elemType = GMESH_LIN_2;
        this.elements[elemPermutation[i]].setNodes();
        this.elements[elemPermutation[i]].tags[1] = 1; // Family
        this.elements[elemPermutation[i]].nodes[1] = nodePermutation[cells[i-1,1]];
        this.elements[elemPermutation[i]].nodes[2] = nodePermutation[cells[i-1,2]];
      }

      // Add left boundary point to elements list
      this.elements[elemPermutation[1]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[1]].setNodes();
      this.elements[elemPermutation[1]].tags[1] = 2; // Family
      this.elements[elemPermutation[1]].nodes[1] = nodePermutation[1];

      // Add right boundary point to elements list
      this.elements[elemPermutation[this.nElements]].elemType = GMESH_PNT_1;
      this.elements[elemPermutation[this.nElements]].setNodes();
      this.elements[elemPermutation[this.nElements]].tags[1] = 3; // Family
      this.elements[elemPermutation[this.nElements]].nodes[1] = nodePermutation[this.nNodes];
    }

    proc read_gmesh_file() {}
    proc write_gmesh_file() {}
  }

  //////////////////////////////////////////////////////////
  //  Gmsh Element node ordering                          //
  //////////////////////////////////////////////////////////

  proc elem_node_order_gmsh(elemType : int) : [] int
  {
    // Convert element node indexing from Gmesh style recursive numbering to Cartesian position
    use Parameters.ParamMesh;
    use Mesh;

    // Get basic element properties
    var elemDim    : int = gmsh_elem_dimension(elemType);
    var elemDegree : int = gmsh_elem_degree(elemType);
    var elemTopo   : int = gmsh_elem_topology(elemType);
    var nodeCnt    : int = gmsh_elem_nodes(elemType);

    // Allocate return vector with the
    var nodeOrder : [1..nodeCnt] int;

    // Select function based of element type
    select elemType
    {
      // Line
      when TYPE_LINE_2 do nodeOrder = node_order_line_gmsh(elemDegree);
      when TYPE_LINE_3 do nodeOrder = node_order_line_gmsh(elemDegree);
      when TYPE_LINE_4 do nodeOrder = node_order_line_gmsh(elemDegree);
      when TYPE_LINE_5 do nodeOrder = node_order_line_gmsh(elemDegree);
      // Tria Full
      // Tria Edge, high-order triangle with nodes only at edges
      // Quad Full
      when TYPE_QUAD_4  do nodeOrder = node_order_quad_gmsh(elemDegree);
      when TYPE_QUAD_9  do nodeOrder = node_order_quad_gmsh(elemDegree);
      when TYPE_QUAD_16 do nodeOrder = node_order_quad_gmsh(elemDegree);
      when TYPE_QUAD_25 do nodeOrder = node_order_quad_gmsh(elemDegree);
      // Quad Edge, high-order quadrilateral with nodes only at edges
      // Tetr Full
      // Tetr Face
      // Tetr Edge, high-order tetrahedral with nodes only at edges
      // Pyra Full
      // Pyra Face
      // Pyra Edge, high-order pyramidal with nodes only at edges
      // Pris Full
      // Pris Face
      // Pris Edge, high-order prismatic with nodes only at edges
      // Hexa Full
      when TYPE_HEXA_8   do nodeOrder = node_order_hexa_gmsh(elemDegree);
      when TYPE_HEXA_27  do nodeOrder = node_order_hexa_gmsh(elemDegree);
      when TYPE_HEXA_64  do nodeOrder = node_order_hexa_gmsh(elemDegree);
      when TYPE_HEXA_125 do nodeOrder = node_order_hexa_gmsh(elemDegree);
      // Hexa Face
      // Hexa Edge, high-order hexahedral with nodes only at edges
      otherwise do writeln("Well, screw node order then...");
    }

    return nodeOrder;
  }

  proc node_order_line_gmsh(elemDegree : int) : [] int
  {
    // Provide an array mapping the cell node index
    //   - From the Gmsh recursive order
    //   - To a Cartesian order. IE nodeIdx[i,j,k] = k*n² + j*n¹ + i
    var nodeCnt   : int = elemDegree+1;
    var nodeOrder : [1..nodeCnt] int;

    // Reorder distribution to node order
    nodeOrder[1] = 1;
    nodeOrder[2] = nodeCnt;
    forall idx in 3..nodeCnt do
      nodeOrder[idx] = idx-1;

    return nodeOrder;
  }

  proc node_order_tria_gmsh(elemDegree : int) : [] int {}

  proc node_order_quad_gmsh(elemDegree : int) : [] int
  {
    // Provide an array mapping the cell node index
    //   - From the Gmsh recursive order
    //   - To a Cartesian order. IE nodeIdx[i,j,k] = k*n² + j*n¹ + i
    use Mesh;

    var nodeCnt   : int = (elemDegree+1)**2;
    var nodeOrder : [1..nodeCnt] int;
    var nodeIdx   : int = 0;

    // Main recursion by order of the element
    for n in 0..elemDegree by -2
    {
      var lo : int = 1            + (elemDegree-n)/2;
      var hi : int = elemDegree+1 - (elemDegree-n)/2;

      if n == 0
      { // Recursion ends with 1 central node
        nodeIdx += 1;
        nodeOrder[nodeIdx] = (lo-1)*(elemDegree+1)+lo;
      }
      else
      { // Fill in corner nodes
        nodeIdx += 1;
        nodeOrder[nodeIdx] = (elemDegree+1)**1 * (lo-1)
                            +(elemDegree+1)**0 *  lo;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (elemDegree+1)**1 * (lo-1)
                            +(elemDegree+1)**0 *  hi;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (elemDegree+1)**1 * (hi-1)
                            +(elemDegree+1)**0 *  hi;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (elemDegree+1)**1 * (hi-1)
                            +(elemDegree+1)**0 *  lo;
      }

      if n > 1
      { // Then we have mid edge nodes at this level
        for k in lo+1..hi-1
        {
          nodeIdx += 1;

          // Bottom edge
          nodeOrder[nodeIdx]         = (elemDegree+1)**1 * (lo-1)
                                      +(elemDegree+1)**0 *  k;
          // Right edge
          nodeOrder[nodeIdx+1*(n-1)] = (elemDegree+1)**1 * ( k-1)
                                      +(elemDegree+1)**0 *  hi;
          // Top edge
          nodeOrder[nodeIdx+2*(n-1)] = (elemDegree+1)**1 * (hi-1)
                                      +(elemDegree+1)**0 * ((elemDegree+1)-(k-1));
          // Left edge
          nodeOrder[nodeIdx+3*(n-1)] = (elemDegree+1)**1 * (elemDegree+1-k)
                                      +(elemDegree+1)**0 *  lo;
        }
        nodeIdx += 3*(n-1);
      }
    }

    return nodeOrder;
  }

  proc node_order_tetr_gmsh(elemDegree : int) : [] int {}
  proc node_order_pyra_gmsh(elemDegree : int) : [] int {}
  proc node_order_pris_gmsh(elemDegree : int) : [] int {}

  proc node_order_hexa_gmsh(elemDegree : int) : [] int
  {
    // Provide an array mapping the cell node index
    //   - From the Gmsh recursive ordering
    //     - Corner Nodes. Back face first following the same ordering as a 2D face
    //     - Mid Edge Nodes
    //       - An edge is oriented from the node with the lowest to the highest index
    //     - Mid Face Nodes
    //       - A face is oriented such that the computed normal points outward
    //       - The starting point is the node with the lowest index
    //     - Mid Volume Nodes
    //       - Recursivelly apply the indexing method described
    //   - To a Cartesian order. IE nodeIdx[n] = k*n² + j*n¹ + i.
    //     Where [i,j,k] are the relative positions of each node in an [x,y,z] coordinate system
    use Mesh;

    var nodeCnt   : int = (elemDegree+1)**3;
    var nodeOrder : [1..nodeCnt] int;
    var nodeIdx   : int = 0;

    // Main recursion by order of the element
    for n in 0..elemDegree by -2
    {
      var lo : int = 1            + (elemDegree-n)/2;
      var hi : int = elemDegree+1 - (elemDegree-n)/2;

      if n == 0
      { // Recursion ends with 1 central node
        nodeIdx += 1;
        nodeOrder[nodeIdx] = (elemDegree+1)**2 * (lo-1)
                            +(elemDegree+1)**1 * (lo-1)
                            +(elemDegree+1)**0 *  lo;
      }
      else
      {
        // Fill in corner nodes
        {
          // Left, lower, back
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (lo-1)
                              +(elemDegree+1)**1 * (lo-1)
                              +(elemDegree+1)**0 *  lo;

          // Right, lower, back
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (lo-1)
                              +(elemDegree+1)**1 * (lo-1)
                              +(elemDegree+1)**0 *  hi;

          // Right, upper, back
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (lo-1)
                              +(elemDegree+1)**1 * (hi-1)
                              +(elemDegree+1)**0 *  hi;

          // Left, upper, back
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (lo-1)
                              +(elemDegree+1)**1 * (hi-1)
                              +(elemDegree+1)**0 *  lo;

          // Left, lower, front
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (hi-1)
                              +(elemDegree+1)**1 * (lo-1)
                              +(elemDegree+1)**0 *  lo;

          // Right, lower, front
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (hi-1)
                              +(elemDegree+1)**1 * (lo-1)
                              +(elemDegree+1)**0 *  hi;

          // Right, upper, front
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (hi-1)
                              +(elemDegree+1)**1 * (hi-1)
                              +(elemDegree+1)**0 *  hi;

          // Left, upper, front
          nodeIdx += 1;
          nodeOrder[nodeIdx] = (elemDegree+1)**2 * (hi-1)
                              +(elemDegree+1)**1 * (hi-1)
                              +(elemDegree+1)**0 *  lo;
        }

        // Is this a high-order element?
        if n > 1
        { // Then we have a few more nodes to add

          // Mid edge nodes
          {
            for k in lo+1..hi-1
            {
              nodeIdx += 1;

              // Edge  1: Bottom, Back
              //x=+k, y=lo, z=lo
              nodeOrder[nodeIdx]         = (elemDegree+1)**2 * (lo-1)
                                          +(elemDegree+1)**1 * (lo-1)
                                          +(elemDegree+1)**0 *   k;
              // Edge  2: Left, Back
              //x=lo, y=+k, z=lo
              nodeOrder[nodeIdx+1*(n-1)] = (elemDegree+1)**2 * (lo-1)
                                          +(elemDegree+1)**1 * ( k-1)
                                          +(elemDegree+1)**0 *  lo;
              // Edge  3: Left, Bottom
              //x=lo, y=lo, z=+k
              nodeOrder[nodeIdx+2*(n-1)] = (elemDegree+1)**2 * ( k-1)
                                          +(elemDegree+1)**1 * (lo-1)
                                          +(elemDegree+1)**0 *  lo;
              // Edge  4: Right, Back
              //x=hi, y=+k, z=lo
              nodeOrder[nodeIdx+3*(n-1)] = (elemDegree+1)**2 * (lo-1)
                                          +(elemDegree+1)**1 * ( k-1)
                                          +(elemDegree+1)**0 *  hi;
              // Edge  5: Right, Bottom
              //x=hi, y=lo, z=+k
              nodeOrder[nodeIdx+4*(n-1)] = (elemDegree+1)**2 * ( k-1)
                                          +(elemDegree+1)**1 * (lo-1)
                                          +(elemDegree+1)**0 *  hi;
              // Edge  6: Top, Back
              //x=-k, y=hi, z=lo
              nodeOrder[nodeIdx+5*(n-1)] = (elemDegree+1)**2 * (lo-1)
                                          +(elemDegree+1)**1 * (hi-1)
                                          +(elemDegree+1)**0 * ((elemDegree+1) - (k-1));
              // Edge  7: Right, top
              //x=hi, y=hi, z=+k
              nodeOrder[nodeIdx+6*(n-1)] = (elemDegree+1)**2 * ( k-1)
                                          +(elemDegree+1)**1 * (hi-1)
                                          +(elemDegree+1)**0 *  hi;
              // Edge  8: Left, Top
              //x=lo, y=hi, z=+k
              nodeOrder[nodeIdx+7*(n-1)] = (elemDegree+1)**2 * ( k-1)
                                          +(elemDegree+1)**1 * (hi-1)
                                          +(elemDegree+1)**0 *  lo;
              // Edge  9: Bottom, Front
              //x=+k, y=lo, z=hi
              nodeOrder[nodeIdx+8*(n-1)] = (elemDegree+1)**2 * (hi-1)
                                          +(elemDegree+1)**1 * (lo-1)
                                          +(elemDegree+1)**0 *   k;
              // Edge 10: Left, Front
              //x=lo, y=+k, z=hi
              nodeOrder[nodeIdx+9*(n-1)] = (elemDegree+1)**2 * (hi-1)
                                          +(elemDegree+1)**1 * ( k-1)
                                          +(elemDegree+1)**0 *  lo;
              // Edge 11: Right, Front
              //x=hi, y=+k, z=hi
              nodeOrder[nodeIdx+10*(n-1)] = (elemDegree+1)**2 * (hi-1)
                                           +(elemDegree+1)**1 * ( k-1)
                                           +(elemDegree+1)**0 *  hi;
              // Edge 12: Top, Front
              //x=-k, y=hi, z=hi
              nodeOrder[nodeIdx+11*(n-1)] = (elemDegree+1)**2 * (hi-1)
                                           +(elemDegree+1)**1 * (hi-1)
                                           +(elemDegree+1)**0 * ((elemDegree+1) - (k-1));
            }

            // Update idx counter to end of edge nodes
            nodeIdx += 11*(n-1);
          }

          // Mid Face nodes
          {
            // Get the node ordering fora n-th degree element
            var quadNodeOrder  : [1..(n-1)**2] int = node_order_quad_gmsh(n-2);
            var quadNodeCoord1 : [1..(n-1)**2] int = lo + (quadNodeOrder-1)%(n-1)+1;
            var quadNodeCoord2 : [1..(n-1)**2] int = lo + (quadNodeOrder-1)/(n-1)+1;

            nodeIdx += 1;

            // Back   face (face 1)
            {
              nodeOrder[nodeIdx+0*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (lo-1)
                                                         +(elemDegree+1)**1 * (quadNodeCoord1-1)
                                                         +(elemDegree+1)**0 * quadNodeCoord2;
            }
            // Bottom face (face 2)
            {
              nodeOrder[nodeIdx+1*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (quadNodeCoord2-1)
                                                         +(elemDegree+1)**1 * (lo-1)
                                                         +(elemDegree+1)**0 * quadNodeCoord1;
            }
            // Left   face (face 3)
            {
              nodeOrder[nodeIdx+2*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (quadNodeCoord1-1)
                                                         +(elemDegree+1)**1 * (quadNodeCoord2-1)
                                                         +(elemDegree+1)**0 * lo;
            }
            // Right  face (face 4)
            {
              nodeOrder[nodeIdx+3*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (quadNodeCoord2-1)
                                                         +(elemDegree+1)**1 * (quadNodeCoord1-1)
                                                         +(elemDegree+1)**0 * hi;
            }
            // Top    face (face 5)
            {
              nodeOrder[nodeIdx+4*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (quadNodeCoord2-1)
                                                         +(elemDegree+1)**1 * (hi-1)
                                                         +(elemDegree+1)**0 * ((elemDegree+1) - (quadNodeCoord1-1));
            }
            // Front  face (face 6)
            {
              nodeOrder[nodeIdx+5*(n-1)**2.. #(n-1)**2] = (elemDegree+1)**2 * (hi-1)
                                                         +(elemDegree+1)**1 * (quadNodeCoord2-1)
                                                         +(elemDegree+1)**0 * quadNodeCoord1;
            }

            nodeIdx += 6*(n-1)**2-1;
          }
        }
      }
    }

    return nodeOrder;
  }

  //////////////////////////////////////////////////////////
  //  Gmsh mesh element inquiry functions                 //
  //////////////////////////////////////////////////////////

  proc gmsh_elem_dimension(in elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1    do return 0;
      when GMESH_LIN_2    do return 1;
      when GMESH_LIN_3    do return 1;
      when GMESH_LIN_4    do return 1;
      when GMESH_LIN_5    do return 1;
      when GMESH_LIN_6    do return 1;
      when GMESH_LIN_7    do return 1;
      when GMESH_LIN_8    do return 1;
      when GMESH_LIN_9    do return 1;
      when GMESH_LIN_10   do return 1;
      when GMESH_LIN_11   do return 1;
      when GMESH_TRI_3    do return 2;
      when GMESH_TRI_6    do return 2;
      when GMESH_TRI_10   do return 2;
      when GMESH_TRI_15   do return 2;
      when GMESH_TRI_21   do return 2;
      when GMESH_TRI_28   do return 2;
      when GMESH_TRI_36   do return 2;
      when GMESH_TRI_45   do return 2;
      when GMESH_TRI_55   do return 2;
      when GMESH_TRI_66   do return 2;
      when GMESH_QUA_4    do return 2;
      when GMESH_QUA_9    do return 2;
      when GMESH_QUA_16   do return 2;
      when GMESH_QUA_25   do return 2;
      when GMESH_QUA_36   do return 2;
      when GMESH_QUA_49   do return 2;
      when GMESH_QUA_64   do return 2;
      when GMESH_QUA_81   do return 2;
      when GMESH_QUA_100  do return 2;
      when GMESH_QUA_121  do return 2;
      when GMESH_TET_4    do return 3;
      when GMESH_TET_10   do return 3;
      when GMESH_TET_20   do return 3;
      when GMESH_TET_35   do return 3;
      when GMESH_TET_56   do return 3;
      when GMESH_TET_84   do return 3;
      when GMESH_TET_120  do return 3;
      when GMESH_TET_165  do return 3;
      when GMESH_TET_220  do return 3;
      when GMESH_TET_286  do return 3;
      when GMESH_PYR_5    do return 3;
      when GMESH_PYR_14   do return 3;
      when GMESH_PYR_30   do return 3;
      when GMESH_PYR_55   do return 3;
      when GMESH_PYR_91   do return 3;
      when GMESH_PYR_140  do return 3;
      when GMESH_PYR_204  do return 3;
      when GMESH_PYR_285  do return 3;
      when GMESH_PYR_385  do return 3;
      when GMESH_PRI_6    do return 3;
      when GMESH_PRI_18   do return 3;
      when GMESH_PRI_40   do return 3;
      when GMESH_PRI_75   do return 3;
      when GMESH_PRI_126  do return 3;
      when GMESH_PRI_196  do return 3;
      when GMESH_PRI_288  do return 3;
      when GMESH_PRI_405  do return 3;
      when GMESH_PRI_550  do return 3;
      when GMESH_HEX_8    do return 3;
      when GMESH_HEX_27   do return 3;
      when GMESH_HEX_64   do return 3;
      when GMESH_HEX_125  do return 3;
      when GMESH_HEX_216  do return 3;
      when GMESH_HEX_343  do return 3;
      when GMESH_HEX_512  do return 3;
      when GMESH_HEX_729  do return 3;
      when GMESH_HEX_1000 do return 3;
      otherwise return -1;
    }
  }

  proc gmsh_elem_topology(in elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1    do return GMESH_PNT;
      when GMESH_LIN_2    do return GMESH_LIN;
      when GMESH_LIN_3    do return GMESH_LIN;
      when GMESH_LIN_4    do return GMESH_LIN;
      when GMESH_LIN_5    do return GMESH_LIN;
      when GMESH_LIN_6    do return GMESH_LIN;
      when GMESH_LIN_7    do return GMESH_LIN;
      when GMESH_LIN_8    do return GMESH_LIN;
      when GMESH_LIN_9    do return GMESH_LIN;
      when GMESH_LIN_10   do return GMESH_LIN;
      when GMESH_LIN_11   do return GMESH_LIN;
      when GMESH_TRI_3    do return GMESH_TRI;
      when GMESH_TRI_6    do return GMESH_TRI;
      when GMESH_TRI_10   do return GMESH_TRI;
      when GMESH_TRI_15   do return GMESH_TRI;
      when GMESH_TRI_21   do return GMESH_TRI;
      when GMESH_TRI_28   do return GMESH_TRI;
      when GMESH_TRI_36   do return GMESH_TRI;
      when GMESH_TRI_45   do return GMESH_TRI;
      when GMESH_TRI_55   do return GMESH_TRI;
      when GMESH_TRI_66   do return GMESH_TRI;
      when GMESH_QUA_4    do return GMESH_QUA;
      when GMESH_QUA_9    do return GMESH_QUA;
      when GMESH_QUA_16   do return GMESH_QUA;
      when GMESH_QUA_25   do return GMESH_QUA;
      when GMESH_QUA_36   do return GMESH_QUA;
      when GMESH_QUA_49   do return GMESH_QUA;
      when GMESH_QUA_64   do return GMESH_QUA;
      when GMESH_QUA_81   do return GMESH_QUA;
      when GMESH_QUA_100  do return GMESH_QUA;
      when GMESH_QUA_121  do return GMESH_QUA;
      when GMESH_TET_4    do return GMESH_TET;
      when GMESH_TET_10   do return GMESH_TET;
      when GMESH_TET_20   do return GMESH_TET;
      when GMESH_TET_35   do return GMESH_TET;
      when GMESH_TET_56   do return GMESH_TET;
      when GMESH_TET_84   do return GMESH_TET;
      when GMESH_TET_120  do return GMESH_TET;
      when GMESH_TET_165  do return GMESH_TET;
      when GMESH_TET_220  do return GMESH_TET;
      when GMESH_TET_286  do return GMESH_TET;
      when GMESH_PYR_5    do return GMESH_PYR;
      when GMESH_PYR_14   do return GMESH_PYR;
      when GMESH_PYR_30   do return GMESH_PYR;
      when GMESH_PYR_55   do return GMESH_PYR;
      when GMESH_PYR_91   do return GMESH_PYR;
      when GMESH_PYR_140  do return GMESH_PYR;
      when GMESH_PYR_204  do return GMESH_PYR;
      when GMESH_PYR_285  do return GMESH_PYR;
      when GMESH_PYR_385  do return GMESH_PYR;
      when GMESH_PRI_6    do return GMESH_PRI;
      when GMESH_PRI_18   do return GMESH_PRI;
      when GMESH_PRI_40   do return GMESH_PRI;
      when GMESH_PRI_75   do return GMESH_PRI;
      when GMESH_PRI_126  do return GMESH_PRI;
      when GMESH_PRI_196  do return GMESH_PRI;
      when GMESH_PRI_288  do return GMESH_PRI;
      when GMESH_PRI_405  do return GMESH_PRI;
      when GMESH_PRI_550  do return GMESH_PRI;
      when GMESH_HEX_8    do return GMESH_HEX;
      when GMESH_HEX_27   do return GMESH_HEX;
      when GMESH_HEX_64   do return GMESH_HEX;
      when GMESH_HEX_125  do return GMESH_HEX;
      when GMESH_HEX_216  do return GMESH_HEX;
      when GMESH_HEX_343  do return GMESH_HEX;
      when GMESH_HEX_512  do return GMESH_HEX;
      when GMESH_HEX_729  do return GMESH_HEX;
      when GMESH_HEX_1000 do return GMESH_HEX;
      otherwise return -1;
    }
  }

  proc gmsh_elem_degree(in elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1    do return  0; // Nodes have no degree
      when GMESH_LIN_2    do return  1; //  1st degree Edge
      when GMESH_LIN_3    do return  2; //  2nd degree Full Edge
      when GMESH_LIN_4    do return  3; //  3rd degree Full Edge
      when GMESH_LIN_5    do return  4; //  4th degree Full Edge
      when GMESH_LIN_6    do return  5; //  5th degree Full Edge
      when GMESH_LIN_7    do return  6; //  6th degree Full Edge
      when GMESH_LIN_8    do return  7; //  7th degree Full Edge
      when GMESH_LIN_9    do return  8; //  8th degree Full Edge
      when GMESH_LIN_10   do return  9; //  9th degree Full Edge
      when GMESH_LIN_11   do return 10; // 10th degree Full Edge
      when GMESH_TRI_3    do return  1; //  1st degree Triangle
      when GMESH_TRI_6    do return  2; //  2nd degree Full Triangle
      when GMESH_TRI_10   do return  3; //  3rd degree Full Triangle
      when GMESH_TRI_15   do return  4; //  4th degree Full Triangle
      when GMESH_TRI_21   do return  5; //  5th degree Full Triangle
      when GMESH_TRI_28   do return  6; //  6th degree Full Triangle
      when GMESH_TRI_36   do return  7; //  7th degree Full Triangle
      when GMESH_TRI_45   do return  8; //  8th degree Full Triangle
      when GMESH_TRI_55   do return  9; //  9th degree Full Triangle
      when GMESH_TRI_66   do return 10; // 10th degree Full Triangle
      when GMESH_QUA_4    do return  1; //  1st degree Quadrilateral
      when GMESH_QUA_9    do return  2; //  2nd degree Full Quadrilateral
      when GMESH_QUA_16   do return  3; //  3rd degree Full Quadrilateral
      when GMESH_QUA_25   do return  4; //  4th degree Full Quadrilateral
      when GMESH_QUA_36   do return  5; //  5th degree Full Quadrilateral
      when GMESH_QUA_49   do return  6; //  6th degree Full Quadrilateral
      when GMESH_QUA_64   do return  7; //  7th degree Full Quadrilateral
      when GMESH_QUA_81   do return  8; //  8th degree Full Quadrilateral
      when GMESH_QUA_100  do return  9; //  9th degree Full Quadrilateral
      when GMESH_QUA_121  do return 10; // 10th degree Full Quadrilateral
      when GMESH_TET_4    do return  1; //  1st degree Tetrahedron
      when GMESH_TET_10   do return  2; //  2nd degree Full Tetrahedron
      when GMESH_TET_20   do return  3; //  3rd degree Full Tetrahedron
      when GMESH_TET_35   do return  4; //  4th degree Full Tetrahedron
      when GMESH_TET_56   do return  5; //  5th degree Full Tetrahedron
      when GMESH_TET_84   do return  6; //  6th degree Full Tetrahedron
      when GMESH_TET_120  do return  7; //  7th degree Full Tetrahedron
      when GMESH_TET_165  do return  8; //  8th degree Full Tetrahedron
      when GMESH_TET_220  do return  9; //  9th degree Full Tetrahedron
      when GMESH_TET_286  do return 10; // 10th degree Full Tetrahedron
      when GMESH_PYR_5    do return  1; //  1st degree Pyramid
      when GMESH_PYR_14   do return  2; //  2nd degree Full Pyramid
      when GMESH_PYR_30   do return  3; //  3rd degree Full Pyramid
      when GMESH_PYR_55   do return  4; //  4th degree Full Pyramid
      when GMESH_PYR_91   do return  5; //  5th degree Full Pyramid
      when GMESH_PYR_140  do return  6; //  6th degree Full Pyramid
      when GMESH_PYR_204  do return  7; //  7th degree Full Pyramid
      when GMESH_PYR_285  do return  8; //  8th degree Full Pyramid
      when GMESH_PYR_385  do return  9; //  9th degree Full Pyramid
      when GMESH_PRI_6    do return  1; //  1st degree Prism
      when GMESH_PRI_18   do return  2; //  2nd degree Full Prism
      when GMESH_PRI_40   do return  3; //  3rd degree Full Prism
      when GMESH_PRI_75   do return  4; //  4th degree Full Prism
      when GMESH_PRI_126  do return  5; //  5th degree Full Prism
      when GMESH_PRI_196  do return  6; //  6th degree Full Prism
      when GMESH_PRI_288  do return  7; //  7th degree Full Prism
      when GMESH_PRI_405  do return  8; //  8th degree Full Prism
      when GMESH_PRI_550  do return  9; //  9th degree Full Prism
      when GMESH_HEX_8    do return  1; //  1st degree Hexahedron
      when GMESH_HEX_27   do return  2; //  2nd degree Full Hexahedron
      when GMESH_HEX_64   do return  3; //  3rd degree Full Hexahedron
      when GMESH_HEX_125  do return  4; //  4th degree Full Hexahedron
      when GMESH_HEX_216  do return  5; //  5th degree Full Hexahedron
      when GMESH_HEX_343  do return  6; //  6th degree Full Hexahedron
      when GMESH_HEX_512  do return  7; //  7th degree Full Hexahedron
      when GMESH_HEX_729  do return  8; //  8th degree Full Hexahedron
      when GMESH_HEX_1000 do return  9; //  9th degree Full Hexahedron
      otherwise return -1;
    }
  }

  proc gmsh_elem_nodes(in elemType : int) : int
  {
    use Parameters.ParamGmesh;

    select elemType {
      when GMESH_PNT_1    do return    1;
      when GMESH_LIN_2    do return    2;
      when GMESH_LIN_3    do return    3;
      when GMESH_LIN_4    do return    4;
      when GMESH_LIN_5    do return    5;
      when GMESH_LIN_6    do return    6;
      when GMESH_LIN_7    do return    7;
      when GMESH_LIN_8    do return    8;
      when GMESH_LIN_9    do return    9;
      when GMESH_LIN_10   do return   10;
      when GMESH_LIN_11   do return   11;
      when GMESH_TRI_3    do return    3;
      when GMESH_TRI_6    do return    6;
      when GMESH_TRI_10   do return   10;
      when GMESH_TRI_15   do return   15;
      when GMESH_TRI_21   do return   21;
      when GMESH_TRI_28   do return   28;
      when GMESH_TRI_36   do return   36;
      when GMESH_TRI_45   do return   45;
      when GMESH_TRI_55   do return   55;
      when GMESH_TRI_66   do return   66;
      when GMESH_QUA_4    do return    4;
      when GMESH_QUA_9    do return    9;
      when GMESH_QUA_16   do return   16;
      when GMESH_QUA_25   do return   25;
      when GMESH_QUA_36   do return   36;
      when GMESH_QUA_49   do return   49;
      when GMESH_QUA_64   do return   64;
      when GMESH_QUA_81   do return   81;
      when GMESH_QUA_100  do return  100;
      when GMESH_QUA_121  do return  121;
      when GMESH_TET_4    do return    4;
      when GMESH_TET_10   do return   10;
      when GMESH_TET_20   do return   20;
      when GMESH_TET_35   do return   35;
      when GMESH_TET_56   do return   56;
      when GMESH_TET_84   do return   84;
      when GMESH_TET_120  do return  120;
      when GMESH_TET_165  do return  165;
      when GMESH_TET_220  do return  220;
      when GMESH_TET_286  do return  286;
      when GMESH_PYR_5    do return    5;
      when GMESH_PYR_14   do return   14;
      when GMESH_PYR_30   do return   30;
      when GMESH_PYR_55   do return   55;
      when GMESH_PYR_91   do return   91;
      when GMESH_PYR_140  do return  140;
      when GMESH_PYR_204  do return  204;
      when GMESH_PYR_285  do return  285;
      when GMESH_PYR_385  do return  385;
      when GMESH_PRI_6    do return    6;
      when GMESH_PRI_18   do return   18;
      when GMESH_PRI_40   do return   40;
      when GMESH_PRI_75   do return   75;
      when GMESH_PRI_126  do return  126;
      when GMESH_PRI_196  do return  196;
      when GMESH_PRI_288  do return  288;
      when GMESH_PRI_405  do return  405;
      when GMESH_PRI_550  do return  550;
      when GMESH_HEX_8    do return    8;
      when GMESH_HEX_27   do return   27;
      when GMESH_HEX_64   do return   64;
      when GMESH_HEX_125  do return  125;
      when GMESH_HEX_216  do return  216;
      when GMESH_HEX_343  do return  343;
      when GMESH_HEX_512  do return  512;
      when GMESH_HEX_729  do return  729;
      when GMESH_HEX_1000 do return 1000;
      otherwise return -1;
    }
  }

  //////////////////////////////////////////////////////////
  //  Module Test                                         //
  //////////////////////////////////////////////////////////

  proc main()
  {
    {
      writeln("Test 1: Random 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.random1D(nCells=6, xMin=-1, xMax=2);
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 2: Uniform 1D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.uniform1D(nCells=6, xMin=-1, xMax=2);
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 3: Read 2D mesh - Gmsh2:");
      var test_gmesh2 = new gmesh2_c();
      test_gmesh2.read_gmesh_file(testMesh);
      writeln(test_gmesh2);
      writeln();
    }

    {
      writeln("Test 4: Node Ordering on Gmsh Elements:");
      writeln("Line q1: ", node_order_line_gmsh(1));
      writeln("Line q2: ", node_order_line_gmsh(2));
      writeln("Line q3: ", node_order_line_gmsh(3));
      writeln("Line q4: ", node_order_line_gmsh(4));
      writeln();
      writeln("Quad q1: ", node_order_quad_gmsh(1));
      writeln("Quad q2: ", node_order_quad_gmsh(2));
      writeln("Quad q3: ", node_order_quad_gmsh(3));
      writeln("Quad q4: ", node_order_quad_gmsh(4));
      writeln();
      writeln("Hexa q1: ", node_order_hexa_gmsh(1));
      writeln("Hexa q2: ", node_order_hexa_gmsh(2));
      writeln("Hexa q3: ", node_order_hexa_gmsh(3));
      writeln("Hexa q4: ", node_order_hexa_gmsh(4));
    }
  }
}
