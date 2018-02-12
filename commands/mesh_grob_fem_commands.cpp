
/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Maxence Reberol
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */
 


#include <OGF/femb/commands/mesh_grob_fem_commands.h>

#include <OGF/femb/algo/femb.h>

namespace OGF {

    MeshGrobFemCommands::MeshGrobFemCommands() { 
    }
        
    MeshGrobFemCommands::~MeshGrobFemCommands() { 
    }        

    void MeshGrobFemCommands::fem(
            const std::string& dirichlet_region,
            const std::string& dirichlet_value,
            const std::string& neumann_region,
            const std::string& neumann_value,
            const std::string& sourceterm_value,
            const std::string& diffusion_value,
            const std::string& solution_name) {

        /* Copy current mesh to a new one for the FEM simulation */
        std::string Mo_name = mesh_grob()->name() + "_fem";
        MeshGrob* Mo = MeshGrob::find_or_create(scene_graph(), Mo_name);
        Mo->copy(*mesh_grob());

        femb::fem_simulation(*Mo, dirichlet_region, dirichlet_value,
                neumann_region, neumann_value, sourceterm_value,
                diffusion_value, solution_name);

        Mo->update();
    }

}
