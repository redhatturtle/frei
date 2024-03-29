module Parameters
{
  module ParamTest
  {
    param RANDOM_SEED : int = 47;
  }

  module ParamInput
  {
    // Equation Sets
    param EQ_CONVECTION     : int = 1;    // Convection Eq            du/dt + c*(du/dx) = 0
    param EQ_INVBURGERS     : int = 2;    // Inviscid Burgers Eq      du/dt + u*(du/dx) = 0
    param EQ_DIFFUSION      : int = 3;    // Diffusion Eq             du/dt             + k*(ddu/dxx) = 0
    param EQ_LINBURGERS     : int = 4;    // Linear Burgers Eq        du/dt + c*(du/dx) + k*(ddu/dxx) = 0
    param EQ_VISBURGERS     : int = 5;    // Viscous Burgers Eq       du/dt + u*(du/dx) + k*(ddu/dxx) = 0
    param EQ_EULER          : int = 6;    // Euler Eq
    param EQ_NAVIERSTOKES   : int = 7;    // Navier-Stokes Eq
    param EQ_QUASI_1D_EULER : int = 8;    // Euler equations for the Quasi 1D Nozzle flow

    // Parametric internal 1D meshing methods
    param MESH_GENERATE    : int = 0;
    param MESH_GMESH       : int = 1;
    param MESH_CGNS        : int = 2;

    param MESH_GEN_UNIFORM : int = 1;
    param MESH_GEN_RANDOM  : int = 2;

    // Spatial Schemes
    param SPATIAL_BEAMWARMING       : int =  1;
    param SPATIAL_LAXWENDROFF       : int =  2;
    param SPATIAL_MACCORMAK         : int =  3;

    param SPATIAL_STEGERWARMING_O1  : int =  4;
    param SPATIAL_STEGERWARMING_O2  : int =  5;

    param SPATIAL_STEGERWARMINGO1   : int =  4;
    param SPATIAL_STEGERWARMINGO2   : int =  5;

    param SPATIAL_STEGER_WARMING_O1 : int =  4;
    param SPATIAL_STEGER_WARMING_O2 : int =  5;

    param SPATIAL_VANLEER_O1        : int =  6;
    param SPATIAL_VANLEER_O2        : int =  7;
    param SPATIAL_AUSM_O1           : int =  8;
    param SPATIAL_AUSMPLUS_O1       : int =  9;
    param SPATIAL_ROE               : int = 10;
    param SPATIAL_FR                : int = 11;

    // Dissipation Scheme
    param DISS_NONE    : int = 0;
    param DISS_SECOND  : int = 1;
    param DISS_FOURTH  : int = 2;
    param DISS_JAMESON : int = 3;

    // Limiter Implementation
    param LIMITER_NONE              : int = 0;
    param LIMITER_PROJECTION        : int = 1;

    // FR correction functions
    param FR_DG  : int = 1;
    param FR_GA  : int = 2; // Lumping at Gauss Points. Similar to SD with internal FPs at Gauss Points.
    param FR_G2  : int = 3;
    param FR_SD  : int = 4; // Lumping at SomeOther Points. Similar to SD with internal FPs at Gauss Points.
    param FR_G3  : int = 5;

    // Inviscid Numerical Flux Schemes
    param FLUX_RUSANOV : int = 1;
    param FLUX_ROE     : int = 2;
    param FLUX_HLL     : int = 3;
    param FLUX_HLLC    : int = 4;
    param FLUX_RHLL    : int = 5;

    // Viscou Numerical Flux Scheme
    param VISC_BR1 : int = 1;
    param VISC_BR2 : int = 2;
    param VISC_LDG : int = 3;

    // Solution Point distributions
    param PTS_UNIFORM          : int = 1;    // Uniform
    param PTS_LEGENDRE         : int = 2;    // Gauss-Legendre
    param PTS_LEGENDRELOBATTO  : int = 3;    // Gauss-Legendre-Lobatto
    param PTS_CHEBYSHEV        : int = 4;    // Gauss-Chebyshev
    param PTS_CHEBYSHEVLOBATTO : int = 5;    // Gauss-Chebyshev-Lobatto

    // Time Schemes
    param TIME_EULER      : int = 0;
    param TIME_RK_CLASSIC : int = 1;
    param TIME_TVDRK_O2S2 : int = 2;
    param TIME_TVDRK_O2S3 : int = 3;
    param TIME_TVDRK_O2S4 : int = 4;
    param TIME_TVDRK_O2SN : int = 5;
    param TIME_TVDRK_O3S3 : int = 6;
    param TIME_TVDRK_O3S4 : int = 7;
    param TIME_TVDRK_O3S5 : int = 8;
    param TIME_TVDRK_O4S5 : int = 9;

    // Time step methods
    param DT_GLOBAL_CONST: int = 0;    // Constant time step for all cells on all interations
    param DT_GLOBAL_CFL  : int = 1;    // Time step calculated at each cell on all iterations, minimum used for all mesh
    param DT_LOCAL_CFL   : int = 2;    // Time step calculated at each cell on all iterations

    //////////////////////
    //   Family Types   //
    //////////////////////

    param BC_TYPE_FLOW                : int =  0;
    param BC_TYPE_INFLOW              : int =  1;
    param BC_TYPE_OUTFLOW             : int =  2;
    param BC_TYPE_OPENING             : int =  3;
    param BC_TYPE_WALL                : int =  4;
    param BC_TYPE_SPECIAL             : int =  9;

    //////////////////////////
    //   Family Sub-Types   //
    //////////////////////////

    // Convection Waves
    param IC_SINUSOIDAL      : int = 11;    // Sinusoidal Wave
    param IC_GAUSSPULSE      : int = 12;    // Gaussian Pulse Wave
    param IC_ELLIPTICALPULSE : int = 13;    // Elliptic Pulse Wave
    param IC_SQUARE          : int = 14;    // Square Wave
    param IC_MIXEDWAVE       : int = 15;    // Mixed Wave

    // 1D Euler
    param IC_SHOCKTUBE       : int = 61;    // Shock Tube

    // 2D Euler
    param IC_RINGLEB         : int = 62;

    // Quasi 1D Nozzle
    param IC_1D_NOZZLE_SUBSONIC          : int = 81;    // Quasi 1D Nozzle
    param IC_1D_NOZZLE_SMOOTH_TRANSONIC  : int = 82;    // Quasi 1D Nozzle
    param IC_1D_NOZZLE_SHOCKED_TRANSONIC : int = 83;    // Quasi 1D Nozzle

    // High-Order Workshop test cases
    param IC_GENERIC_MEANFLOW        : int = 60;
    param IC_CHANNEL_FLOW            : int;
    param IC_INVISCID_GAUSSIAN_BUMP  : int;
    param IC_LAMINAR_BNDLYR          : int;
    param IC_SHU_VORTEX              : int;
    param IC_VORTEX_TRANSPORT        : int;
    param IC_DENSITY_TRANSPORT       : int;
    param IC_ENTROPY_TRANSPORT       : int;
    param IC_TAYLOR_GREEN_VORTEX     : int;
    param IC_NOZZLE_JET              : int;
    param IC_DOUBLE_MACH_REFLECTION  : int;
    param IC_INFINITE_CYLINDER       : int;
    param IC_FREESTREAM_PRESERVATION : int;

    // Inflow BCs
    param BC_SUBTYPE_GENERIC_INFLOW   : int = 11;
    param BC_SUBTYPE_SUB_INFLOW       : int = 12;
    param BC_SUBTYPE_SUP_INFLOW       : int = 13;
    param BC_SUBTYPE_MDOT_INFLOW      : int = 14;
    // Outflow BCs
    param BC_SUBTYPE_GENERIC_OUTFLOW  : int = 21;
    param BC_SUBTYPE_SUB_OUTFLOW      : int = 22;
    param BC_SUBTYPE_SUP_OUTFLOW      : int = 23;
    param BC_SUBTYPE_MDOT_OUTFLOW     : int = 24;
    // Generic Inflow/Outflow BCs
    param BC_SUBTYPE_RIEMANN          : int = 31;
    param BC_SUBTYPE_GENERIC_FREEFLOW : int = 32;
    param BC_SUBTYPE_FREESTREAM       : int = 33;
    param BC_SUBTYPE_FIXED            : int = 34;
    // Wall BCs
    param BC_SUBTYPE_SLIP_WALL        : int = 41;
    param BC_SUBTYPE_EULER_WALL       : int = 42;
    param BC_SUBTYPE_ADIABATIC_WALL   : int = 43;
    param BC_SUBTYPE_ISOTHERMAL_WALL  : int = 44;
    param BC_SUBTYPE_DEFAULT_WALL     : int = 45;
    // Special BCs
    param BC_SUBTYPE_SYMMETRY         : int = 91;
    param BC_SUBTYPE_PERIODIC         : int = 92;
    param BC_SUBTYPE_DIRICHLET        : int = 93;
    param BC_SUBTYPE_1D_NOZZLE_CRIT_INFLOW     : int = 94;
    param BC_SUBTYPE_1D_NOZZLE_SUBSONIC_INFLOW : int = 95;
    param BC_SUBTYPE_1D_NOZZLE_SHOCKED_INFLOW  : int = 96;
    param BC_SUBTYPE_RINGLEB_DIRICHLET: int = 97;
    param BC_SUBTYPE_MMS_DIRICHLET    : int = 99;
  }

  module ParamConstants
  {
    param PI : real = 3.1415926535897932384626433832795028841971;

    param EPS4 :  real = 1.0e-4;
    param EPS6  : real = 1.0e-6;
    param EPS8  : real = 1.0e-8;
    param EPS10 : real = 1.0e-10;
    param EPS12 : real = 1.0e-12;
    param EPS13 : real = 1.0e-13;
    param EPS14 : real = 1.0e-14;
    param EPS15 : real = 1.0e-15;
    param EPS16 : real = 1.0e-16;
    param EPS17 : real = 1.0e-17;
  }

  module ParamMesh
  {
    // Mesh element topologies
    param TOPO_NODE : int = 1; // Point
    param TOPO_LINE : int = 2; // Line
    param TOPO_TRIA : int = 3; // Triangle
    param TOPO_QUAD : int = 4; // Quadrangle
    param TOPO_TETR : int = 5; // Tetrahedron
    param TOPO_PYRA : int = 6; // Pyramid
    param TOPO_PRIS : int = 7; // Prism
    param TOPO_HEXA : int = 8; // Hexahedron

    // Mesh element types
    param TYPE_NODE     : int = 10;
    param TYPE_LINE_2   : int = 21; // 1st Order Line
    param TYPE_LINE_3   : int = 22; // 2nd Order Line
    param TYPE_LINE_4   : int = 23; // 3rd Order Line
    param TYPE_LINE_5   : int = 24; // 4th Order Line
    param TYPE_TRIA_3   : int = 31; // 1st Order Triangle
    param TYPE_TRIA_6   : int = 32; // 2nd Order Triangle
    param TYPE_TRIA_10  : int = 33; // 3rd Order Triangle
    param TYPE_TRIA_15  : int = 34; // 4th Order Triangle
    param TYPE_QUAD_4   : int = 41; // 1st Order Quadrilateral
    param TYPE_QUAD_9   : int = 42; // 2nd Order Quadrilateral
    param TYPE_QUAD_16  : int = 43; // 3rd Order Quadrilateral
    param TYPE_QUAD_25  : int = 44; // 4th Order Quadrilateral
    param TYPE_TETR_4   : int = 51; // 1st Order Tetrahedra
    param TYPE_TETR_10  : int = 52; // 2nd Order Tetrahedra
    param TYPE_TETR_20  : int = 53; // 3rd Order Tetrahedra
    param TYPE_TETR_35  : int = 54; // 4th Order Tetrahedra
    param TYPE_PYRA_5   : int = 61; // 1st Order Pyramid
    param TYPE_PYRA_14  : int = 62; // 2nd Order Pyramid
    param TYPE_PYRA_30  : int = 63; // 3rd Order Pyramid
    param TYPE_PYRA_55  : int = 64; // 4th Order Pyramid
    param TYPE_PRIS_6   : int = 71; // 1st Order Prism
    param TYPE_PRIS_18  : int = 72; // 2nd Order Prism
    param TYPE_PRIS_40  : int = 73; // 3rd Order Prism
    param TYPE_PRIS_75  : int = 74; // 4th Order Prism
    param TYPE_HEXA_8   : int = 81; // 1st Order Hexahedron
    param TYPE_HEXA_27  : int = 82; // 2nd Order Hexahedron
    param TYPE_HEXA_64  : int = 83; // 3rd Order Hexahedron
    param TYPE_HEXA_125 : int = 84; // 4th Order Hexahedron

    {
      const NODE_ORDER_QUAD_4   : [1..4] int = [1,2,
                                                4,3];

      const NODE_ORDER_QUAD_9   : [1..9] int = [1,5,2,
                                                8,9,6,
                                                4,7,3];

      const NODE_ORDER_QUAD_16  : [1..16] int = [ 1, 5, 6, 2,
                                                 12,13,14, 7,
                                                 11,16,15, 8,
                                                  4,10, 9, 3];

      const NODE_ORDER_QUAD_25  : [1..25] int = [ 1, 5, 6, 7, 2,
                                                 16,17,21,18, 8,
                                                 15,24,25,22, 9,
                                                 14,20,23,19,10,
                                                  4,13,12,11, 3];
    }
  }

  module ParamCGNS
  {
    // CGNS Element Types    - https://cgns.github.io/CGNS_docs_current/sids/conv.html
    // Relevant CPEX Proposals:
    //  - P4 Extension Proposal - https://cgns.github.io/ProposedExtensions/GNS_P4_elem_defn2.pdf
    //  - Generalized curved mesh and cell-wise polynomial - https://cgns.github.io/ProposedExtensions/CPEX0045_HighOrder_v2.pdf
    //
    //Dimension | Shape         | Linear    | Quadratic           | Cubic                        | Quartic
    //---------------------------------------------------------------------------------------------------------------------------
    //  0-D     |   Point       |   NODE    |   NODE              |   NODE                       |   NODE
    //  1-D     |   Line        |   BAR_2   |   BAR_3             |   BAR_4                      |   BAR_5
    //  2-D     |   Triangle    |   TRI_3   |   TRI_6             |   TRI_9    TRI_10            |   TRI_12     TRI_15
    //          |   Quadrangle  |   QUAD_4  |   QUAD_8   QUAD_9   |   QUAD_12  QUAD_16           |   QUAD_P4_16 QUAD_25
    //  3-D     |   Tetrahedron |   TETRA_4 |   TETRA_10          |   TETRA_16 TETRA_20          |   TETRA_22   TETRA_34 TETRA_35
    //          |   Pyramid     |   PYRA_5  |   PYRA_13  PYRA_14  |   PYRA_21  PYRA_29  PYRA_30  |   PYRA_P4_29 PYRA_50  PYRA_55
    //          |   Pentahedron |   PENTA_6 |   PENTA_15 PENTA_18 |   PENTA_24 PENTA_38 PENTA_40 |   PENTA_33   PENTA_66 PENTA_75
    //          |   Hexahedron  |   HEXA_8  |   HEXA_20  HEXA_27  |   HEXA_32  HEXA_56  HEXA_64  |   HEXA_44    HEXA_98  HEXA_125
    param CGNS_NODE       : int;

    param CGNS_BAR_2      : int;
    param CGNS_BAR_3      : int;
    param CGNS_BAR_4      : int;
    param CGNS_BAR_5      : int;

    param CGNS_TRI_3      : int;
    param CGNS_TRI_6      : int;
    param CGNS_TRI_9      : int;
    param CGNS_TRI_10     : int;
    param CGNS_TRI_12     : int;
    param CGNS_TRI_15     : int;

    param CGNS_QUAD_4     : int;
    param CGNS_QUAD_8     : int;
    param CGNS_QUAD_9     : int;
    param CGNS_QUAD_12    : int;
    param CGNS_QUAD_16    : int;
    param CGNS_QUAD_P4_16 : int;
    param CGNS_QUAD_25    : int;

    param CGNS_TETRA_4    : int;
    param CGNS_TETRA_10   : int;
    param CGNS_TETRA_16   : int;
    param CGNS_TETRA_20   : int;
    param CGNS_TETRA_22   : int;
    param CGNS_TETRA_34   : int;
    param CGNS_TETRA_35   : int;

    param CGNS_PYRA_5     : int;
    param CGNS_PYRA_13    : int;
    param CGNS_PYRA_14    : int;
    param CGNS_PYRA_21    : int;
    param CGNS_PYRA_29    : int;
    param CGNS_PYRA_30    : int;
    param CGNS_PYRA_P4_29 : int;
    param CGNS_PYRA_50    : int;
    param CGNS_PYRA_55    : int;

    param CGNS_PENTA_6    : int;
    param CGNS_PENTA_15   : int;
    param CGNS_PENTA_18   : int;
    param CGNS_PENTA_24   : int;
    param CGNS_PENTA_38   : int;
    param CGNS_PENTA_40   : int;
    param CGNS_PENTA_33   : int;
    param CGNS_PENTA_66   : int;
    param CGNS_PENTA_75   : int;

    param CGNS_HEXA_8     : int;
    param CGNS_HEXA_20    : int;
    param CGNS_HEXA_27    : int;
    param CGNS_HEXA_32    : int;
    param CGNS_HEXA_56    : int;
    param CGNS_HEXA_64    : int;
    param CGNS_HEXA_44    : int;
    param CGNS_HEXA_98    : int;
    param CGNS_HEXA_125   : int;
  }

  module ParamGmesh
  {
    // Gmesh element topologies
    param GMESH_PNT     : int =  1;
    param GMESH_LIN     : int =  2;
    param GMESH_TRI     : int =  3;
    param GMESH_QUA     : int =  4;
    param GMESH_TET     : int =  5;
    param GMESH_PYR     : int =  6;
    param GMESH_PRI     : int =  7;
    param GMESH_HEX     : int =  8;
    param GMESH_POLYG   : int =  9;
    param GMESH_POLYH   : int = 10;
    param GMESH_XFEM    : int = 11;
    param GMESH_MINI    : int = 12;
    param GMESH_TRIH    : int = 13;

    // Gmesh Element Types
    //
    //Dimension | Shape         | Linear    | Quadratic           | Cubic                        | Quartic
    //---------------------------------------------------------------------------------------------------------------------------
    //  1-D     |   Line        |   LIN_2   |   LIN_3             |   LIN_4                      |   LIN_5
    //  2-D     |   Triangle    |   TRI_3   |   TRI_6             |   TRI_9    TRI_10            |   TRI_12     TRI_15
    //          |   Quadrangle  |   QUA_4   |   QUA_8    QUA_9    |   QUA_12   QUA_16            |   QUA_16I    QUAD_25
    //  3-D     |   Tetrahedron |   TET_4   |   TET_10            |   TET_16   TET_20            |   TET_22     TET_34   TET_35
    //          |   Pyramid     |   PYR_5   |   PYR_13   PYR_14   |   PYR_21   PYR_29   PYR_30   |   PYR_29I    PYR_50   PYR_55
    //          |   Prism       |   PRI_6   |   PRI_15   PRI_18   |   PRI_24   PRI_38   PRI_40   |   PRI_33     PRI_66   PRI_75
    //          |   Hexahedron  |   HEX_8   |   HEX_20   HEX_27   |   HEX_32   HEX_56   HEX_64   |   HEX_44     HEX_98   HEX_125

    param GMESH_PNT_1    : int = 15;
    param GMESH_TRIH_4   : int = 140; // Trihedron

    //     Line:                 Line3:          Line4:
    //
    //       v
    //       ^
    //       |
    //       |
    // 0-----+-----1 --> u   0----2----1     0---2---3---1

    // First order element
    param GMESH_LIN_2    : int =   1;
    // Complete high-order elements
    param GMESH_LIN_3    : int =   8;
    param GMESH_LIN_4    : int =  26;
    param GMESH_LIN_5    : int =  27;
    param GMESH_LIN_6    : int =  28;
    param GMESH_LIN_7    : int =  62;
    param GMESH_LIN_8    : int =  63;
    param GMESH_LIN_9    : int =  64;
    param GMESH_LIN_10   : int =  65;
    param GMESH_LIN_11   : int =  66;
    // 0th order elements
    param GMESH_LIN_1    : int =  84;
    // Exotic elements
    param GMESH_LIN_B    : int =  67;
    param GMESH_LIN_C    : int =  70;
    param GMESH_LIN_SUB  : int = 134;

    // Triangle:               Triangle6:          Triangle9/10:          Triangle12/15:
    //
    // v
    // ^                                                                   2
    // |                                                                   | \
    // 2                       2                    2                      9   8
    // |`\                     |`\                  | \                    |     \
    // |  `\                   |  `\                7   6                 10 (14)  7
    // |    `\                 5    `4              |     \                |         \
    // |      `\               |      `\            8  (9)  5             11 (12) (13) 6
    // |        `\             |        `\          |         \            |             \
    // 0----------1 --> u      0-----3----1         0---3---4---1          0---3---4---5---1

    // First order element
    param GMESH_TRI_3    : int =   2;
    // Complete high-order elements
    param GMESH_TRI_6    : int =   9;
    param GMESH_TRI_10   : int =  21;
    param GMESH_TRI_15   : int =  23;
    param GMESH_TRI_21   : int =  25;
    param GMESH_TRI_28   : int =  42;
    param GMESH_TRI_36   : int =  43;
    param GMESH_TRI_45   : int =  44;
    param GMESH_TRI_55   : int =  45;
    param GMESH_TRI_66   : int =  46;
    // Edge based high-order elements
    param GMESH_TRI_9    : int =  20;
    param GMESH_TRI_12   : int =  22;
    param GMESH_TRI_15I  : int =  24;
    param GMESH_TRI_18   : int =  52;
    param GMESH_TRI_21I  : int =  53;
    param GMESH_TRI_24   : int =  54;
    param GMESH_TRI_27   : int =  55;
    param GMESH_TRI_30   : int =  56;
    // 0th order elements
    param GMESH_TRI_1    : int =  85;
    // Exotic elements
    param GMESH_TRI_B    : int =  68;
    param GMESH_TRI_SUB  : int = 135;
    param GMESH_TRI_MINI : int = 138;

    // Quadrangle:            Quadrangle8:            Quadrangle9:
    //
    //       v
    //       ^
    //       |
    // 3-----------2          3-----6-----2           3-----6-----2
    // |     |     |          |           |           |           |
    // |     |     |          |           |           |           |
    // |     +---- | --> u    7           5           7     8     5
    // |           |          |           |           |           |
    // |           |          |           |           |           |
    // 0-----------1          0-----4-----1           0-----4-----1

    // First order element
    param GMESH_QUA_4    : int =  3;
    // Complete high-order elements
    param GMESH_QUA_9    : int = 10;
    param GMESH_QUA_16   : int = 36;
    param GMESH_QUA_25   : int = 37;
    param GMESH_QUA_36   : int = 38;
    param GMESH_QUA_49   : int = 47;
    param GMESH_QUA_64   : int = 48;
    param GMESH_QUA_81   : int = 49;
    param GMESH_QUA_100  : int = 50;
    param GMESH_QUA_121  : int = 51;
    // Edge based high-order elements
    param GMESH_QUA_8    : int = 16;
    param GMESH_QUA_12   : int = 39;
    param GMESH_QUA_16I  : int = 40;
    param GMESH_QUA_20   : int = 41;
    param GMESH_QUA_24   : int = 57;
    param GMESH_QUA_28   : int = 58;
    param GMESH_QUA_32   : int = 59;
    param GMESH_QUA_36I  : int = 60;
    param GMESH_QUA_40   : int = 61;
    // 0th order elements
    param GMESH_QUA_1    : int = 86;

    // Tetrahedron:                          Tetrahedron10:
    //
    //                    v
    //                  .
    //                ,/
    //               /
    //            2                                     2
    //          ,/|`\                                 ,/|`\
    //        ,/  |  `\                             ,/  |  `\
    //      ,/    '.   `\                         ,6    '.   `5
    //    ,/       |     `\                     ,/       8     `\
    //  ,/         |       `\                 ,/         |       `\
    // 0-----------'.--------1 --> u         0--------4--'.--------1
    //  `\.         |      ,/                 `\.         |      ,/
    //     `\.      |    ,/                      `\.      |    ,9
    //        `\.   '. ,/                           `7.   '. ,/
    //           `\. |/                                `\. |/
    //              `3                                    `3
    //                 `\.
    //                    ` w

    // First order element
    param GMESH_TET_4    : int =   4;
    // Complete high-order elements
    param GMESH_TET_10   : int =  11;
    param GMESH_TET_20   : int =  29;
    param GMESH_TET_35   : int =  30;
    param GMESH_TET_56   : int =  31;
    param GMESH_TET_84   : int =  71;
    param GMESH_TET_120  : int =  72;
    param GMESH_TET_165  : int =  73;
    param GMESH_TET_220  : int =  74;
    param GMESH_TET_286  : int =  75;
    // Edge based high-order elements
    // 2nd order edge based element is identical to complete element
    param GMESH_TET_16   : int = 137;
    param GMESH_TET_22   : int =  32;
    param GMESH_TET_28   : int =  33;
    param GMESH_TET_34   : int =  79;
    param GMESH_TET_40   : int =  80;
    param GMESH_TET_46   : int =  81;
    param GMESH_TET_52   : int =  82;
    param GMESH_TET_58   : int =  83;
    // 0th order elements
    param GMESH_TET_1    : int =  87;
    // Exotic elements
    param GMESH_TET_SUB  : int = 136;
    param GMESH_TET_MINI : int = 139;

    // Pyramid:                     Pyramid13:                   Pyramid14:
    //
    //                4                            4                            4
    //              ,/|\                         ,/|\                         ,/|\
    //            ,/ .'|\                      ,/ .'|\                      ,/ .'|\
    //          ,/   | | \                   ,/   | | \                   ,/   | | \
    //        ,/    .' | `.                ,/    .' | `.                ,/    .' | `.
    //      ,/      |  '.  \             ,7      |  12  \             ,7      |  12  \
    //    ,/       .' w |   \          ,/       .'   |   \          ,/       .'   |   \
    //  ,/         |  ^ |    \       ,/         9    |    11      ,/         9    |    11
    // 0----------.'--|-3    `.     0--------6-.'----3    `.     0--------6-.'----3    `.
    //  `\        |   |  `\    \      `\        |      `\    \     `\        |      `\    \
    //    `\     .'   +----`\ - \ -> v  `5     .'        10   \      `5     .' 13     10   \
    //      `\   |    `\     `\  \        `\   |           `\  \       `\   |           `\  \
    //        `\.'      `\     `\`          `\.'             `\`         `\.'             `\`
    //           1----------------2            1--------8-------2           1--------8-------2
    //                     `\
    //                       u

    // First order element
    param GMESH_PYR_5    : int =   7;
    // Complete high-order elements
    param GMESH_PYR_14   : int =  14;
    param GMESH_PYR_30   : int = 118;
    param GMESH_PYR_55   : int = 119;
    param GMESH_PYR_91   : int = 120;
    param GMESH_PYR_140  : int = 121;
    param GMESH_PYR_204  : int = 122;
    param GMESH_PYR_285  : int = 123;
    param GMESH_PYR_385  : int = 124;
    // Edge based high-order elements
    param GMESH_PYR_13   : int =  19;
    param GMESH_PYR_21   : int = 125;
    param GMESH_PYR_29   : int = 126;
    param GMESH_PYR_37   : int = 127;
    param GMESH_PYR_45   : int = 128;
    param GMESH_PYR_53   : int = 129;
    param GMESH_PYR_61   : int = 130;
    param GMESH_PYR_69   : int = 131;
    // 0th order elements
    param GMESH_PYR_1    : int = 132;

    // Prism:                      Prism15:               Prism18:
    //
    //            w
    //            ^
    //            |
    //            3                       3                      3
    //          ,/|`\                   ,/|`\                  ,/|`\
    //        ,/  |  `\               12  |  13              12  |  13
    //      ,/    |    `\           ,/    |    `\          ,/    |    `\
    //     4------+------5         4------14-----5        4------14-----5
    //     |      |      |         |      8      |        |      8      |
    //     |    ,/|`\    |         |      |      |        |    ,/|`\    |
    //     |  ,/  |  `\  |         |      |      |        |  15  |  16  |
    //     |,/    |    `\|         |      |      |        |,/    |    `\|
    //    ,|      |      |\        10     |      11       10-----17-----11
    //  ,/ |      0      | `\      |      0      |        |      0      |
    // u   |    ,/ `\    |    v    |    ,/ `\    |        |    ,/ `\    |
    //     |  ,/     `\  |         |  ,6     `7  |        |  ,6     `7  |
    //     |,/         `\|         |,/         `\|        |,/         `\|
    //     1-------------2         1------9------2        1------9------2

    // First order element
    param GMESH_PRI_6    : int =   6;
    // Complete high-order elements
    param GMESH_PRI_18   : int =  13;
    param GMESH_PRI_40   : int =  90;
    param GMESH_PRI_75   : int =  91;
    param GMESH_PRI_126  : int = 106;
    param GMESH_PRI_196  : int = 107;
    param GMESH_PRI_288  : int = 108;
    param GMESH_PRI_405  : int = 109;
    param GMESH_PRI_550  : int = 110;
    // Edge based high-order elements
    param GMESH_PRI_15   : int =  18;
    param GMESH_PRI_24   : int = 111;
    param GMESH_PRI_33   : int = 112;
    param GMESH_PRI_42   : int = 113;
    param GMESH_PRI_51   : int = 114;
    param GMESH_PRI_60   : int = 115;
    param GMESH_PRI_69   : int = 116;
    param GMESH_PRI_78   : int = 117;
    // 0th order elements
    param GMESH_PRI_1    : int =  89;

    // Hexahedron:             Hexahedron20:          Hexahedron27:
    //
    //        v
    // 3----------2            3----13----2           3----13----2
    // |\     ^   |\           |\         |\          |\         |\
    // | \    |   | \          | 15       | 14        |15    24  | 14
    // |  \   |   |  \         9  \       11 \        9  \ 20    11 \
    // |   7------+---6        |   7----19+---6       |   7----19+---6
    // |   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
    // 0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
    //  \  |    \  \  |         \  17      \  18       \ 17    25 \  18
    //   \ |     \  \ |         10 |        12|        10 |  21    12|
    //    \|      w  \|           \|         \|          \|         \|
    //     4----------5            4----16----5           4----16----5

    // First order element
    param GMESH_HEX_8    : int =   5;
    // Complete high-order elements
    param GMESH_HEX_27   : int =  12;
    param GMESH_HEX_64   : int =  92;
    param GMESH_HEX_125  : int =  93;
    param GMESH_HEX_216  : int =  94;
    param GMESH_HEX_343  : int =  95;
    param GMESH_HEX_512  : int =  96;
    param GMESH_HEX_729  : int =  97;
    param GMESH_HEX_1000 : int =  98;
    // Edge based high-order elements
    param GMESH_HEX_20   : int =  17;
    param GMESH_HEX_32   : int =  99;
    param GMESH_HEX_44   : int = 100;
    param GMESH_HEX_56   : int = 101;
    param GMESH_HEX_68   : int = 102;
    param GMESH_HEX_80   : int = 103;
    param GMESH_HEX_92   : int = 104;
    param GMESH_HEX_104  : int = 105;
    // 0th order elements
    param GMESH_HEX_1    : int =  88;
  }
}
