scene_graph.current().query_interface("OGF::MeshGrobFemCommands").fem({sourceterm_value="0", dirichlet_region="if(y > 0.95, 1, -1)", neumann_value="1.", solution_name="u", diffusion_value="1", dirichlet_value="0", neumann_region="if(y < 0.10, 1, -1)"})

-- Visualisation
scene_graph.current_object = scene_graph.current_object .. "_fem"
scene_graph.current().shader.mesh_style='true; 0 0 0 1; 2'
sg_shader_mgr.current().painting="ATTRIBUTE"
sg_shader_mgr.current().attribute="vertices.u"
main.camera().set_property("draw_selected_only","true")
