Graphite plugin that implements P1 FEM for Poisson equation on tetrahedral meshes.

The interest of this code is mainly to show how to use the Geogram and Graphite
APIs, functions and classes. The FEM part is very basic and simple.

See the [tutorial](docs/tutorial.md) for more information.

The FEM implementation is in the file `algo/femb.cpp`

To run simulations from command line:

   graphite  data/cube_s12.meshb data/sinbump_test.lua 
   graphite  data/cube_s12.meshb data/neumann_test.lua 


Validation: same results than with *MFEM* on `sinbump_test.lua` and `neumann_test.lua`




