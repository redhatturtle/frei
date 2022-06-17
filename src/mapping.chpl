module Mapping
{
  use Random;
  use UnitTest;
  use Set;

  class element_mapping_c
  {
    var coefs_d : domain(2); // {sp_id, elem_node}
    var coefs   : [coefs_d] real;
  }

  class element_mapping_metrics_c
  {
    var coefs_d : domain(3); // {derivative_direction, sp_id, elem_node}
    var coefs   : [coefs_d] real;
  }

  type element_mapping_t = unmanaged element_mapping_c?;
  type element_mapping_metrics_t = unmanaged element_mapping_metrics_c?;

  // Domains
  var mapping_d : domain(2*int); // {element_type, interpolation_order}
  var mappingMetrics_d : domain(2*int); // {element_type, interpolation_order}

  // Mapping coefficients structure
  var mapping : [mapping_d] element_mapping_t;
  var mappingMetrics : [mappingMetrics_d] element_mapping_metrics_t;

  proc init_mapping(minOrder : int, maxOrder : int, elemTypes : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    writeln();
    writeln("Initializing Mapping matrices");
    writeln("    Element types: ", elemTypes);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for elemType in elemTypes do
      for interpOrder in minOrder..maxOrder
      {
        // Get mapping matrix dimensions
        var interpPtCnt : int = n_cell_sps(elem_topology(elemType), interpOrder);
        var elemNodeCnt : int = elem_nodes(elemType);

        // Allocate mapping matrix
        mapping_d.add((elemType, interpOrder));
        mapping[(elemType, interpOrder)] = new element_mapping_t({1..interpPtCnt, 1..elemNodeCnt})!;
      }

    for (elemType, interpOrder) in mapping.domain
    {
      select elem_topology(elemType)
      {
        when TOPO_NODE
        {
          mapping[(elemType, interpOrder)]!.coefs = 1.0;
        }
        when TOPO_LINE
        {
          // Check if the element type has a complete set of nodes, is face defined or edge defined
          // if (elemComplete == true)
          {
            // Assuming complete set of nodes
            mapping[(elemType, interpOrder)]!.coefs = mapping_line(elemType, interpOrder);
          }
        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          // Check if the element type has a complete set of nodes, is face defined or edge defined
          // if (elemComplete == true)
          {
            // Assuming a Quad mesh element with a full set of nodes
            // Assuming a tensor product SP set
            mapping[(elemType, interpOrder)]!.coefs = mapping_quad(elemType, interpOrder);
          }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA
        {
          // Check if the element type has a complete set of nodes, is face defined or edge defined
          // if (elemComplete == true)
          {
            // Assuming a Quad mesh element with a full set of nodes
            // Assuming a tensor product SP set
            mapping[(elemType, interpOrder)]!.coefs = mapping_hexa(elemType, interpOrder);
          }
        }
        otherwise
        {
          writeln("Unsupported mesh element found at element mapping initialization.");
          writeln("Element Type: ", elemType, ", Element Topology: ", elem_topology(elemType));
          writeln();
        }
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  proc init_mapping_metrics(minOrder : int, maxOrder : int, elemTypes : set(int))
  {
    use Time;
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    writeln();
    writeln("Initializing Mapping Metrics matrices");
    writeln("    Element types: ", elemTypes);
    writeln("    Minimum Polynomial Degree: ", minOrder);
    writeln("    Maximum Polynomial Degree: ", maxOrder);
    var stopwatch : Timer;
    stopwatch.start();

    // Add all combination of cell topology and interpolation order to the domain
    for elemType in elemTypes do
      for interpOrder in minOrder..maxOrder
      {
        // Get mapping matrix dimensions
        var interpPtCnt : int = n_cell_sps(elem_topology(elemType), interpOrder);
        var elemNodeCnt : int = elem_nodes(elemType);
        var elemDims    : int = elem_dimension(elem_topology(elemType));

        // Allocate mapping metrics matrix
        mappingMetrics_d.add((elemType, interpOrder));
        mappingMetrics[(elemType, interpOrder)] = new element_mapping_metrics_t({1..elemDims, 1..interpPtCnt, 1..elemNodeCnt})!;
      }

    for (elemType, interpOrder) in mappingMetrics.domain
    {
      select elem_topology(elemType)
      {
        when TOPO_NODE {} // There are no metrics for a single node mapping
        when TOPO_LINE
        {
          // Check if the element type has a complete set of nodes, is face defined or edge defined
          // if (elemComplete == true)
          {
            // Assuming complete set of nodes
            mappingMetrics[(elemType, interpOrder)]!.coefs = mapping_metrics_line(elemType, interpOrder);
          }
        }
        when TOPO_TRIA {}
        when TOPO_QUAD
        {
          // Check if the element type has a complete set of nodes, is face defined or edge defined
          // if (elemComplete == true)
          {
            // Assuming a Quad mesh element with a full set of nodes
            // Assuming a tensor product SP set
            mappingMetrics[(elemType, interpOrder)]!.coefs = mapping_metrics_quad(elemType, interpOrder);
          }
        }
        when TOPO_TETR {}
        when TOPO_PYRA {}
        when TOPO_PRIS {}
        when TOPO_HEXA {}
        otherwise
        {
          writeln("Unsupported mesh element found at mapping metrics initialization.");
          writeln("Element Type: ", elemType, ", Element Topology: ", elem_topology(elemType));
          writeln();
        }
      }
    }

    writef("    Initialized in  %6.1dr ms\n", stopwatch.elapsed(TimeUnits.milliseconds));
  }

  //////////////////////////////////////////////////////////
  //  Build SP map for each element                       //
  //////////////////////////////////////////////////////////

  /*
   Calculate mapping matrix for the interpolation points is a line element
   */
  proc mapping_line(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions
    var interpPtCnt : int = n_cell_sps(TOPO_LINE, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);

    // Allocate mapping matrix
    var mappingCoefs : [1..interpPtCnt, 1..elemNodeCnt] real;

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemNodeCnt]   real = nodes_uniform_lobatto(elemNodeCnt);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
    //    Need to build an appropriate way to query the point location for each element.
    //    Initially assume the whole mesh uses the same base distribution specified in input file.
    //    Even more initially assume the whole mesh has SPs on Legendre roots.

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi  : real = spLoc[ptIdx];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = nodeOrder[nodeIdx];

        mappingCoefs[ptIdx, nodeIdx] = eval_LagrangePoly1D(xi, i, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_tria(elemType : int, interpOrder : int) : [] real {}

  proc mapping_quad(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions and allocate it
    var interpPtCnt : int = n_cell_sps(TOPO_QUAD, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);

    // Allocate mapping matrix
    var mappingCoefs : [1..interpPtCnt, 1..elemNodeCnt] real;

    var elemDegree   : int = elem_degree(elemType);

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemDegree+1]  real = nodes_uniform_lobatto(elemDegree+1);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi  : real = spLoc[1 + (ptIdx-1)%(interpOrder+1)];
      var eta : real = spLoc[1 + (ptIdx-1)/(interpOrder+1)];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = 1 + (nodeOrder[nodeIdx]-1)%(elemDegree+1);
        var j : int = 1 + (nodeOrder[nodeIdx]-1)/(elemDegree+1);

        mappingCoefs[ptIdx, nodeIdx] = eval_LagrangePoly1D(xi , i, nodeLoc)
                                      *eval_LagrangePoly1D(eta, j, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_tetr(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pyra(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pris(elemType : int, interpOrder : int) : [] real {}
  proc mapping_hexa(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions and allocate it
    var interpPtCnt : int = n_cell_sps(TOPO_HEXA, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);

    // Allocate mapping matrix
    var mappingCoefs : [1..interpPtCnt, 1..elemNodeCnt] real;

    var elemDegree   : int = elem_degree(elemType);

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemDegree+1]  real = nodes_uniform_lobatto(elemDegree+1);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi   : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**0)%(interpOrder+1)];
      var eta  : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**1)%(interpOrder+1)];
      var zeta : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**2)%(interpOrder+1)];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**0)%(elemDegree+1);
        var j : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**1)%(elemDegree+1);
        var k : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**2)%(elemDegree+1);

        mappingCoefs[ptIdx, nodeIdx] = eval_LagrangePoly1D(xi  , i, nodeLoc)
                                      *eval_LagrangePoly1D(eta , j, nodeLoc)
                                      *eval_LagrangePoly1D(zeta, k, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_tetr_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pyra_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pris_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_hexa_face(elemType : int, interpOrder : int) : [] real {}

  proc mapping_tria_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_quad_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_tetr_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pyra_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_pris_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_hexa_edge(elemType : int, interpOrder : int) : [] real {}

  //////////////////////////////////////////////////////////
  //  Build SP map for elem                               //
  //////////////////////////////////////////////////////////

  /*
   Calculate mapping matrix for the interpolation points is a line element
   */
  proc mapping_metrics_line(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions
    var interpPtCnt : int = n_cell_sps(TOPO_LINE, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);
    var elemDims    : int = elem_dimension(TOPO_LINE);

    // Allocate mapping matrix
    var mappingCoefs : [1..elemDims, 1..interpPtCnt, 1..elemNodeCnt] real;

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemNodeCnt] real = nodes_uniform_lobatto(elemNodeCnt);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);
    //    Need to build an appropriate way to query the point location for each element.
    //    Initially assume the whole mesh uses the same base distribution specified in input file.
    //    Even more initially assume the whole mesh has SPs on Legendre roots.

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi  : real = spLoc[ptIdx];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = nodeOrder[nodeIdx];

        mappingCoefs[1, ptIdx, nodeIdx] = eval_DLagrangeDx(xi, i, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_metrics_tria(elemType : int, interpOrder : int) : [] real {}

  proc mapping_metrics_quad(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions and allocate it
    var interpPtCnt : int = n_cell_sps(TOPO_QUAD, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);
    var elemDims    : int = elem_dimension(TOPO_QUAD);

    // Allocate mapping matrix
    var mappingCoefs : [1..elemDims, 1..interpPtCnt, 1..elemNodeCnt] real;

    var elemDegree   : int = elem_degree(elemType);

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemDegree+1]  real = nodes_uniform_lobatto(elemDegree+1);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi  : real = spLoc[1 + (ptIdx-1)%(interpOrder+1)];
      var eta : real = spLoc[1 + (ptIdx-1)/(interpOrder+1)];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = 1 + (nodeOrder[nodeIdx]-1)%(elemDegree+1);
        var j : int = 1 + (nodeOrder[nodeIdx]-1)/(elemDegree+1);

        mappingCoefs[1, ptIdx, nodeIdx] = eval_DLagrangeDx(   xi , i, nodeLoc)
                                         *eval_LagrangePoly1D(eta, j, nodeLoc);

        mappingCoefs[2, ptIdx, nodeIdx] = eval_LagrangePoly1D(xi , i, nodeLoc)
                                         *eval_DLagrangeDx(   eta, j, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_metrics_tetr(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pyra(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pris(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_hexa(elemType : int, interpOrder : int) : [] real
  {
    use Parameters.ParamMesh;
    use Mesh;
    use FRMesh;
    use Polynomials;
    use Interpolation;

    // Get mapping matrix dimensions and allocate it
    var interpPtCnt : int = n_cell_sps(TOPO_HEXA, interpOrder);
    var elemNodeCnt : int = elem_nodes(elemType);
    var elemDims    : int = elem_dimension(TOPO_HEXA);

    // Allocate mapping matrix
    var mappingCoefs : [1..elemDims, 1..interpPtCnt, 1..elemNodeCnt] real;

    var elemDegree   : int = elem_degree(elemType);

    // Get the mesh element and interpolation nodes coordinates/distribution
    var nodeLoc : [1..elemDegree+1]  real = nodes_uniform_lobatto(elemDegree+1);
    var spLoc   : [1..interpOrder+1] real = nodes_legendre_gauss(interpOrder+1);

    // Get the node ordering for this element type
    var nodeOrder : [1..elemNodeCnt] int = node_order_elem(elemType);

    // Build mapping matrix
    for ptIdx in 1..interpPtCnt
    {
      // Computational coordinate of the SP spIdx
      var xi   : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**0)%(interpOrder+1)];
      var eta  : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**1)%(interpOrder+1)];
      var zeta : real = spLoc[1 + ((ptIdx-1)/(interpOrder+1)**2)%(interpOrder+1)];

      for nodeIdx in 1..elemNodeCnt
      {
        // Relative position of the mesh element's node
        var i : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**0)%(elemDegree+1);
        var j : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**1)%(elemDegree+1);
        var k : int = 1 + ((nodeOrder[nodeIdx]-1)/(interpOrder+1)**2)%(elemDegree+1);

        mappingCoefs[1, ptIdx, nodeIdx] = eval_DLagrangeDx(   xi  , i, nodeLoc)
                                         *eval_LagrangePoly1D(eta , j, nodeLoc)
                                         *eval_LagrangePoly1D(zeta, k, nodeLoc);

        mappingCoefs[2, ptIdx, nodeIdx] = eval_LagrangePoly1D(xi  , i, nodeLoc)
                                         *eval_DLagrangeDx(   eta , j, nodeLoc)
                                         *eval_LagrangePoly1D(zeta, k, nodeLoc);

        mappingCoefs[3, ptIdx, nodeIdx] = eval_LagrangePoly1D(xi  , i, nodeLoc)
                                         *eval_LagrangePoly1D(eta , j, nodeLoc)
                                         *eval_DLagrangeDx(   zeta, k, nodeLoc);
      }
    }

    return mappingCoefs;
  }

  proc mapping_metrics_tetr_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pyra_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pris_face(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_hexa_face(elemType : int, interpOrder : int) : [] real {}

  proc mapping_metrics_tria_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_quad_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_tetr_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pyra_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_pris_edge(elemType : int, interpOrder : int) : [] real {}
  proc mapping_metrics_hexa_edge(elemType : int, interpOrder : int) : [] real {}

  //////////////////////////////////////////////////////////
  //  Element nodes ordering                              //
  //////////////////////////////////////////////////////////

  proc node_order_elem(elemType : int) : [] int
  {
    // Convert element node indexing from Gmesh style recursive numbering to Cartesian position
    use Parameters.ParamMesh;
    use Mesh;

    // Get basic element properties
    var elemDim   : int = elem_dimension(elemType);
    var elemDegree : int = elem_degree(elemType);
    var elemTopo  : int = elem_topology(elemType);
    var nodeCnt   : int = elem_nodes(elemType);

    // Allocate return vector with the
    var nodeOrder : [1..nodeCnt] int;

    // Select function based of element type
    select elemType
    {
      // Line
      when TYPE_LINE_2 do nodeOrder = node_order_line(elemDegree);
      when TYPE_LINE_3 do nodeOrder = node_order_line(elemDegree);
      when TYPE_LINE_4 do nodeOrder = node_order_line(elemDegree);
      when TYPE_LINE_5 do nodeOrder = node_order_line(elemDegree);
      // Tria Full
      // Tria Edge, high-order triangle with nodes only at edges
      // Quad Full
      when TYPE_QUAD_4  do nodeOrder = node_order_quad(elemDegree);
      when TYPE_QUAD_9  do nodeOrder = node_order_quad(elemDegree);
      when TYPE_QUAD_16 do nodeOrder = node_order_quad(elemDegree);
      when TYPE_QUAD_25 do nodeOrder = node_order_quad(elemDegree);
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
      when TYPE_HEXA_8   do nodeOrder = node_order_hexa(elemDegree);
      when TYPE_HEXA_27  do nodeOrder = node_order_hexa(elemDegree);
      when TYPE_HEXA_64  do nodeOrder = node_order_hexa(elemDegree);
      when TYPE_HEXA_125 do nodeOrder = node_order_hexa(elemDegree);
      // Hexa Face
      // Hexa Edge, high-order hexahedral with nodes only at edges
      otherwise do writeln("Well, screw node order then...");
    }

    return nodeOrder;
  }

  proc node_order_line(elemDegree : int) : [] int
  {
    var nodeCnt   : int = elemDegree+1;
    var nodeOrder : [1..nodeCnt] int;

    // Reorder distribution to node order
    nodeOrder[1] = 1;
    nodeOrder[2] = nodeCnt;
    forall idx in 3..nodeCnt do
      nodeOrder[idx] = idx-1;

    return nodeOrder;
  }

  proc node_order_tria(elemDegree : int) : [] int {}

  proc node_order_quad(elemDegree : int) : [] int
  {
    // Provide an array mapping the cell node index
    //   - From the mesh recursive order
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
        nodeOrder[nodeIdx] = (lo-1)*(elemDegree+1)+lo;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (lo-1)*(elemDegree+1)+hi;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (hi-1)*(elemDegree+1)+hi;

        nodeIdx += 1;
        nodeOrder[nodeIdx] = (hi-1)*(elemDegree+1)+lo;
      }

      if n > 1
      { // Then we have mid edge nodes at this level
        for k in lo+1..hi-1
        {
          nodeIdx += 1;

          // Bottom edge
          nodeOrder[nodeIdx]         = (lo-1)*(elemDegree+1)+k;
          // Right edge
          nodeOrder[nodeIdx+1*(n-1)] = ( k-1)*(elemDegree+1)+hi;
          // Top edge
          nodeOrder[nodeIdx+2*(n-1)] = (hi-1)*(elemDegree+1)+((elemDegree+1)-(k-1));
          // Left edge
          nodeOrder[nodeIdx+3*(n-1)] = (elemDegree+1-k)*(elemDegree+1)+lo;
        }
        nodeIdx += 3*(n-1);
      }
    }

    return nodeOrder;
  }

  proc node_order_tetr(elemDegree : int) : [] int {}
  proc node_order_pyra(elemDegree : int) : [] int {}
  proc node_order_pris(elemDegree : int) : [] int {}

  proc node_order_hexa(elemDegree : int) : [] int
  {
    // Provide an array mapping the cell node index
    //   - From the native mesh ordering
    //       - Corner Nodes. Back face first following the same ordering as a 2D face
    //       - Mid Edge Nodes (Following edge orientation)
    //       - Mid Face Nodes (Cartesian and following face orientation)
    //       - Mid Volume Nodes (Cartesian and following volume orientation)
    //   - To a Cartesian order. IE nodeIdx[n] = k*n² + j*n¹ + i.
    //     Where [i,j,k] are the relative positions of each node in an [x,y,z] coordinate system
    use Mesh;

    var nodeCnt   : int = (elemDegree+1)**3;
    var nodeOrder : [1..nodeCnt] int;
    var nodeIdx   : int = 0;

    var lo = 1;
    var hi = elemDegree+1;

    // Fill in corner nodes
    {
      // Left, lower, back
      nodeOrder[1] = (elemDegree+1)**2 * (lo-1)*
                    +(elemDegree+1)**1 * (lo-1)*
                    +(elemDegree+1)**0 *  lo;

      // Right, lower, back
      nodeOrder[2] = (elemDegree+1)**2 * (lo-1)
                    +(elemDegree+1)**1 * (lo-1)
                    +(elemDegree+1)**0 *  hi;

      // Right, upper, back
      nodeOrder[3] = (elemDegree+1)**2 * (lo-1)
                    +(elemDegree+1)**1 * (hi-1)
                    +(elemDegree+1)**0 *  hi;

      // Left, upper, back
      nodeOrder[4] = (elemDegree+1)**2 * (lo-1)
                    +(elemDegree+1)**1 * (hi-1)
                    +(elemDegree+1)**0 *  lo;

      // Left, lower, front
      nodeOrder[5] = (elemDegree+1)**2 * (hi-1)
                    +(elemDegree+1)**1 * (lo-1)
                    +(elemDegree+1)**0 *  lo;

      // Right, lower, front
      nodeOrder[6] = (elemDegree+1)**2 * (hi-1)
                    +(elemDegree+1)**1 * (lo-1)
                    +(elemDegree+1)**0 *  hi;

      // Right, upper, front
      nodeOrder[7] = (elemDegree+1)**2 * (hi-1)
                    +(elemDegree+1)**1 * (hi-1)
                    +(elemDegree+1)**0 *  hi;

      // Left, upper, front
      nodeOrder[8] = (elemDegree+1)**2 * (hi-1)
                    +(elemDegree+1)**1 * (hi-1)
                    +(elemDegree+1)**0 *  lo;

      // Update idx counter to end of the corner nodes
      nodeIdx = 8;
    }

    // Is this a high-order element?
    if elemDegree > 1
    {
      // Then we have a few more nodes to add

      // Mid edge nodes
      {
        for k in lo+1..hi-1
        {
          nodeIdx += 1;

          // Back edges (z=zMin). Oriented counter clockwise
          {
            // Bottom edge
            nodeOrder[nodeIdx+0*(elemDegree-1)] = (elemDegree+1)**2 * (lo-1)
                                                 +(elemDegree+1)**1 * (lo-1)
                                                 +(elemDegree+1)**0 *   k;
            // Right edge
            nodeOrder[nodeIdx+1*(elemDegree-1)] = (elemDegree+1)**2 * (lo-1)
                                                 +(elemDegree+1)**1 * ( k-1)
                                                 +(elemDegree+1)**0 *  hi;
            // Top edge
            nodeOrder[nodeIdx+2*(elemDegree-1)] = (elemDegree+1)**2 * (lo-1)
                                                 +(elemDegree+1)**1 * (hi-1)
                                                 +(elemDegree+1)**0 * ((elemDegree+1)-(k-1));
            // Left edge
            nodeOrder[nodeIdx+3*(elemDegree-1)] = (elemDegree+1)**2 * (lo-1)
                                                 +(elemDegree+1)**1 * ((elemDegree+1)-(k-1)-1)
                                                 +(elemDegree+1)**0 *  lo;
          }

          // Z edges. Oriented back to front
          {
            // Bottom-Left edge
            nodeOrder[nodeIdx+4*(elemDegree-1)] = (elemDegree+1)**2 * ( k-1)
                                                 +(elemDegree+1)**1 * (lo-1)
                                                 +(elemDegree+1)**0 *  lo;
            // Bottom-Right edge
            nodeOrder[nodeIdx+5*(elemDegree-1)] = (elemDegree+1)**2 * ( k-1)
                                                 +(elemDegree+1)**1 * (lo-1)
                                                 +(elemDegree+1)**0 *  hi;
            // Top-Right edge
            nodeOrder[nodeIdx+6*(elemDegree-1)] = (elemDegree+1)**2 * ( k-1)
                                                 +(elemDegree+1)**1 * (hi-1)
                                                 +(elemDegree+1)**0 *  hi;
            // Top-Left edge
            nodeOrder[nodeIdx+7*(elemDegree-1)] = (elemDegree+1)**2 * ( k-1)
                                                 +(elemDegree+1)**1 * (hi-1)
                                                 +(elemDegree+1)**0 *  lo;
          }

          // Front edges (z=zMax). Oriented counter clockwise
          {
            // Bottom edge
            nodeOrder[nodeIdx+8*(elemDegree-1)] = (elemDegree+1)**2 * (hi-1)
                                                 +(elemDegree+1)**1 * (lo-1)
                                                 +(elemDegree+1)**0 *   k;
            // Right edge
            nodeOrder[nodeIdx+9*(elemDegree-1)] = (elemDegree+1)**2 * (hi-1)
                                                 +(elemDegree+1)**1 * ( k-1)
                                                 +(elemDegree+1)**0 *  hi;
            // Top edge
            nodeOrder[nodeIdx+10*(elemDegree-1)] = (elemDegree+1)**2 * (hi-1)
                                                  +(elemDegree+1)**1 * (hi-1)
                                                  +(elemDegree+1)**0 * ((elemDegree+1)-(k-1));
            // Left edge
            nodeOrder[nodeIdx+11*(elemDegree-1)] = (elemDegree+1)**2 * (hi-1)
                                                  +(elemDegree+1)**1 * ((elemDegree+1)-(k-1)-1)
                                                  +(elemDegree+1)**0 *  lo;
          }
        }

        // Update idx counter to end of edge nodes
        nodeIdx += 11*(elemDegree-1);
      }

      // Mid Face nodes
      {
        // Back   face (face 1)
        {
          var z = lo;
          for x in lo+1..hi-1 do
            for y in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
        // Bottom face (face 2)
        {
          var y = lo;
          for z in lo+1..hi-1 do
            for x in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
        // Left   face (face 3)
        {
          var x = lo;
          for y in lo+1..hi-1 do
            for z in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
        // Front  face (face 4)
        {
          var z = hi;
          for y in lo+1..hi-1 do
            for x in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
        // Top    face (face 5)
        {
          var y = hi;
          for x in lo+1..hi-1 do
            for z in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
        // Right  face (face 6)
        {
          var x = hi;
          for z in lo+1..hi-1 do
            for y in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
        }
      }

      // Mid Volume nodes
      {
        for z in lo+1..hi-1 do
          for y in lo+1..hi-1 do
            for x in lo+1..hi-1
            {
              nodeIdx += 1;
              nodeOrder[nodeIdx] = (elemDegree+1)**2 * z
                                  +(elemDegree+1)**1 * y
                                  +(elemDegree+1)**0 * x;
            }
      }
    }

    return nodeOrder;
  }

  proc node_order_tetr_face(elemDegree : int) : [] int {}
  proc node_order_pyra_face(elemDegree : int) : [] int {}
  proc node_order_pris_face(elemDegree : int) : [] int {}
  proc node_order_hexa_face(elemDegree : int) : [] int {}

  proc node_order_tria_edge(elemDegree : int) : [] int {}
  proc node_order_quad_edge(elemDegree : int) : [] int {}
  proc node_order_tetr_edge(elemDegree : int) : [] int {}
  proc node_order_pyra_edge(elemDegree : int) : [] int {}
  proc node_order_pris_edge(elemDegree : int) : [] int {}
  proc node_order_hexa_edge(elemDegree : int) : [] int {}

  //////////////////////////////////////////////////////////
  //  Module Test                                         //
  //////////////////////////////////////////////////////////

  proc main()
  {
    use IO;
    use Set;
    use LinearAlgebra;
    use Testing;
    use Parameters.ParamTest;
    use Parameters.ParamMesh;
    use Polynomials;
    import FRMesh.determinant;

    // Create a set with the element types contained in a hypothetical test mesh
    var cellTypes : set(int);
    cellTypes.add(TYPE_LINE_2);
    cellTypes.add(TYPE_LINE_3);
    cellTypes.add(TYPE_LINE_4);
    cellTypes.add(TYPE_LINE_5);
    cellTypes.add(TYPE_QUAD_4);
    cellTypes.add(TYPE_QUAD_9);
    cellTypes.add(TYPE_QUAD_16);
    cellTypes.add(TYPE_QUAD_25);

    // Calculate mapping matrices for all these elements
    init_mapping(minOrder=0, maxOrder=3, cellTypes);
    init_mapping_metrics(minOrder=0, maxOrder=3, cellTypes);

    writeln();
    writeln("Element mapping:");
    writeln(mapping);
    writeln();
    writeln("Element mapping metrics:");
    writeln(mappingMetrics);
    writeln();
    writeln("------------------------------------------------------------");

    // Test mapping a 2x3 rectangle
    {
      // Get the coordinates of these nodes
      var xyzMshNodes : [1..2, 1..4] real = 0.0;
      xyzMshNodes[1, 1..4] = [0.0, 3.0, 3.0, 0.0];
      xyzMshNodes[2, 1..4] = [0.0, 0.0, 2.0, 2.0];

      writeln();
      writeln("Node Coordinates:");
      writeln(xyzMshNodes);

      for solOrder in 0..2
      {
        var cellType : 2*int = (TYPE_QUAD_4, solOrder);
        var spCnt : int = (solOrder+1)**2;

        var xyzSP : [1..2, 1..spCnt] real;
        var metSP : [1..2, 1..2, 1..spCnt] real;
        var jacSP : [1..spCnt] real;

        writeln();
        writeln("------------------------------------------------------------");

        //////////////////////////
        ///   SP coordinates   ///
        //////////////////////////
        for spIdx in 1..spCnt do
          xyzSP[.., spIdx] = dot( mapping[cellType]!.coefs[spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Coordinates:");
        writef("%.4t\n", xyzSP);

        /////////////////////////////////////
        ///   Calculate mapping metrics   ///
        /////////////////////////////////////
        for rstDim in 1..2 do
          for spIdx in 1..spCnt do
            metSP[rstDim, .., spIdx] = dot( mappingMetrics[cellType]!.coefs[rstDim, spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Metric Terms:");
        for xyzDim in 1..2 do
          for rstDim in 1..2 do
            writef("%.4t\n", metSP[rstDim, xyzDim, ..]);

        //////////////////////////////
        ///   Calculate Jacobian   ///
        //////////////////////////////
        // Calculate the Jacobian at SPs
        for spIdx in 1..spCnt do
          jacSP[spIdx] = determinant(metSP[.., .., spIdx]);

        writeln();
        writef("Jacobians: %.4t", jacSP);
      }
    }

    // Test mapping a 1x2 parallelogram
    {
      // Get the coordinates of these nodes
      var xyzMshNodes : [1..2, 1..4] real = 0.0;
      xyzMshNodes[1, 1..4] = [0.0, 2.0, 3.0, 1.0];
      xyzMshNodes[2, 1..4] = [0.0, 0.0, 1.0, 1.0];

      writeln();
      writeln("Node Coordinates:");
      writeln(xyzMshNodes);

      for solOrder in 0..2
      {
        var cellType : 2*int = (TYPE_QUAD_4, solOrder);
        var spCnt : int = (solOrder+1)**2;

        var xyzSP : [1..2, 1..spCnt] real;
        var metSP : [1..2, 1..2, 1..spCnt] real;
        var jacSP : [1..spCnt] real;

        writeln();
        writeln("------------------------------------------------------------");

        //////////////////////////
        ///   SP coordinates   ///
        //////////////////////////
        for spIdx in 1..spCnt do
          xyzSP[.., spIdx] = dot( mapping[cellType]!.coefs[spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Coordinates:");
        writef("%.4t\n", xyzSP);

        /////////////////////////////////////
        ///   Calculate mapping metrics   ///
        /////////////////////////////////////
        for rstDim in 1..2 do
          for spIdx in 1..spCnt do
            metSP[rstDim, .., spIdx] = dot( mappingMetrics[cellType]!.coefs[rstDim, spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Metric Terms:");
        for xyzDim in 1..2 do
          for rstDim in 1..2 do
            writef("%.4t\n", metSP[rstDim, xyzDim, ..]);

        //////////////////////////////
        ///   Calculate Jacobian   ///
        //////////////////////////////
        // Calculate the Jacobian at SPs
        for spIdx in 1..spCnt do
          jacSP[spIdx] = determinant(metSP[.., .., spIdx]);

        writeln();
        writef("Jacobians: %.4t", jacSP);
      }
    }

    // Test mapping a trapezoid quadrilateral
    {
      // Get the coordinates of these nodes
      var xyzMshNodes : [1..2, 1..4] real = 0.0;
      xyzMshNodes[1, 1..4] = [0.0, 3.0, 2.0, 1.0];
      xyzMshNodes[2, 1..4] = [0.0, 0.0, 1.0, 1.0];

      writeln();
      writeln("Node Coordinates:");
      writeln(xyzMshNodes);

      for solOrder in 0..2
      {
        var cellType : 2*int = (TYPE_QUAD_4, solOrder);
        var spCnt : int = (solOrder+1)**2;

        var xyzSP : [1..2, 1..spCnt] real;
        var metSP : [1..2, 1..2, 1..spCnt] real;
        var jacSP : [1..spCnt] real;

        writeln();
        writeln("------------------------------------------------------------");

        //////////////////////////
        ///   SP coordinates   ///
        //////////////////////////
        for spIdx in 1..spCnt do
          xyzSP[.., spIdx] = dot( mapping[cellType]!.coefs[spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Coordinates:");
        writef("%.4t\n", xyzSP);

        /////////////////////////////////////
        ///   Calculate mapping metrics   ///
        /////////////////////////////////////
        for rstDim in 1..2 do
          for spIdx in 1..spCnt do
            metSP[rstDim, .., spIdx] = dot( mappingMetrics[cellType]!.coefs[rstDim, spIdx, ..], xyzMshNodes[..,..].T );

        writeln();
        writeln("SP Metric Terms:");
        for xyzDim in 1..2 do
          for rstDim in 1..2 do
            writef("%.4t\n", metSP[rstDim, xyzDim, ..]);

        //////////////////////////////
        ///   Calculate Jacobian   ///
        //////////////////////////////
        // Calculate the Jacobian at SPs
        for spIdx in 1..spCnt do
          jacSP[spIdx] = determinant(metSP[.., .., spIdx]);

        writeln();
        writef("Jacobians: %.4t", jacSP);
      }
    }
  }
}
