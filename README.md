Graphite plugin that implements P1 FEM for Poisson equation on tetrahedral meshes.

The interest of this code is mainly to show how to use the Geogram and Graphite
APIs, functions and classes. The FEM part is very basic and simple.

If you want to learn how to develop a plugin for *Graphite* and *Geogram*, please use the
[tutorial](docs/tutorial.md) associated to this plugin.

If you are only interested in the FEM implementation (~200 lines), you can directly look 
at the file `algo/femb.cpp`

Validation: same results than with *MFEM* on `sinbump_test.lua` and `neumann_test.lua`

### Installation

Assuming that you have a working *Graphite* setup, the steps are:

    // get the plugin source code
    mkdir /path/to/this/plugin
    git clone https://github.com/mxncr/femb.git /path/to/this/plugin

    // tell graphite to build the new plugin
    ln -s /path/to/this/plugin /path/to/graphite/plugins/OGF/femb
    echo "add_subdirectory(femb)" >> /path/to/graphite/plugins/OGF/Plugins.txt

Now you should:
- configure and build Graphite (which will build this plugin)
- launch Graphite, go to *Files > Preferences > Plugins*, enter *femb* and click *Add*, then *Save configuration file*
- restart Graphite, from now on you should have the *femb* plugin working (right click on a mesh in Scene and check if the FEM command menu is present)

### Run FEM simulations

Two possibilities:

**a)** From command line, you can load a mesh and run a Lua script. See the Lua scripts in `data/` for examples. This will open Graphite with the simulation output.

    graphite  data/cube_s12.meshb data/sinbump_test.lua 
    graphite  data/cube_s12.meshb data/neumann_test.lua 

**b)** From the Graphite GUI: open Graphite, load a mesh (right click on Scene), right click on the new mesh
and use the *FEM* menu.

After the simulation, a new mesh is created in the scene (with the suffix \_fem) with the solution as an vertex attribute (named `u` by default). Right click
on the new mesh and use the *Properties* menu to visualize the attributes.


