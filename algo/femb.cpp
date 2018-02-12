#include "OGF/femb/algo/femb.h"

#include "OGF/femb/algo/functions.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_remesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_tetrahedralize.h>
#include <geogram/points/co3ne.h>
#include <geogram/third_party/PoissonRecon/poisson_geogram.h>
#include <geogram/basic/matrix.h>
#include <geogram/basic/logger.h>
#include <geogram/NL/nl.h>

namespace femb {
    using namespace GEO;

    std::vector<double> triangle_o2 = { 
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    }; 

    std::vector<double> triangle_o4 = { 
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    }; 

    std::vector<double> tetra_o2 = { 
        0.0416666666666667, 0.138196601125011, 0.138196601125011, 0.138196601125011,
        0.0416666666666667, 0.138196601125011, 0.138196601125011, 0.585410196624968,
        0.0416666666666667, 0.138196601125011, 0.585410196624968, 0.138196601125011,
        0.0416666666666667, 0.585410196624968, 0.138196601125011, 0.138196601125011
    }; 

    std::vector<double> tetra_o4 = { 
        0.00762222222222222, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
        0.00762222222222222, 0.0714285714285714, 0.0714285714285714, 0.785714285714286,
        0.00762222222222222, 0.0714285714285714, 0.785714285714286, 0.0714285714285714,
        0.00762222222222222, 0.785714285714286, 0.0714285714285714, 0.0714285714285714,
        -0.0131555555555556, 0.25, 0.25, 0.25,
        0.0248888888888889, 0.100596423833201, 0.100596423833201, 0.399403576166799,
        0.0248888888888889, 0.100596423833201, 0.399403576166799, 0.100596423833201,
        0.0248888888888889, 0.399403576166799, 0.100596423833201, 0.100596423833201,
        0.0248888888888889, 0.399403576166799, 0.399403576166799, 0.100596423833201,
        0.0248888888888889, 0.399403576166799, 0.100596423833201, 0.399403576166799,
        0.0248888888888889, 0.100596423833201, 0.399403576166799, 0.399403576166799
    }; 

    typedef GEO::Matrix<3,double> Mat;

    double determinant(Mat M) {
        return det3x3(M.data()[0], M.data()[1], M.data()[2], M.data()[3],
                M.data()[4], M.data()[5], M.data()[6], M.data()[7],
                M.data()[8]);
    }

    vec3 mapping_tet_transform(const Mesh& M, index_t c, vec3 x_ref) {
        const vec3 p0 = M.vertices.point(M.cells.vertex(c,0));
        const vec3 p1 = M.vertices.point(M.cells.vertex(c,1));
        const vec3 p2 = M.vertices.point(M.cells.vertex(c,2));
        const vec3 p3 = M.vertices.point(M.cells.vertex(c,3));
        return (1.-x_ref.x-x_ref.y-x_ref.z)*p0 + x_ref.x*p1 + x_ref.y * p2 + x_ref.z * p3;
    }

    /* does not depend on x_ref because P1 mappings */
    Mat mapping_tet_jacobian(const Mesh& M, index_t c) {
        /* Mapping: p = (1-u-v-w)*p_0 + u*p_1 + v*p_2 + w*p_3
         * jacobian matrix: J(i,j) = dF_i/dx_j */
        const vec3 p0 = M.vertices.point(M.cells.vertex(c,0));
        const vec3 p1 = M.vertices.point(M.cells.vertex(c,1));
        const vec3 p2 = M.vertices.point(M.cells.vertex(c,2));
        const vec3 p3 = M.vertices.point(M.cells.vertex(c,3));
        Mat J; 
        J.load_zero();
        J(0,0) = - p0.x + p1.x; J(0,1) = - p0.x + p2.x; J(0,2) = - p0.x + p3.x; 
        J(1,0) = - p0.y + p1.y; J(1,1) = - p0.y + p2.y; J(1,2) = - p0.y + p3.y; 
        J(2,0) = - p0.z + p1.z; J(2,1) = - p0.z + p2.z; J(2,2) = - p0.z + p3.z; 
        return J;
    }

    double phi_tet_eval(index_t i, vec3 x_ref) {
        geo_debug_assert(i >= 0 && i < 4);
        switch( i ) {
            case 0:
                return 1. - x_ref.x - x_ref.y - x_ref.z ;
            case 1:
                return x_ref.x ;
            case 2:
                return x_ref.y ;
            case 3:
                return x_ref.z ;
        }
        return 0. ; // should not be reached
    }

    /* does not depend on x_ref because P1 tets */
    vec3 grad_phi_tet_eval(index_t i) {
        geo_debug_assert(i >= 0 && i < 4);
        switch( i ) {
            case 0:
                return vec3(-1.,-1.,-.1);
            case 1:
                return vec3(1.,0.,0.);
            case 2:
                return vec3(0.,1.,0.);
            case 3:
                return vec3(0.,0.,1.);
        }
        return vec3(0.,0.,0.) ; // should not be reached
    }

    vec3 matvec(Mat M, vec3 v) {
        vec3 r(0.,0.,0.);
        mult(M,v.data(),r.data()) ;
        return r;
    }

    bool fem_simulation(Mesh& M,
            const std::string& dirichlet_region,
            const std::string& dirichlet_value,
            const std::string& neumann_region,
            const std::string& neumann_value,
            const std::string& sourceterm_value,
            const std::string& diffusion_value,
            const std::string& solution_name
            ) {

        /* ------ MESH PREPROCESSING ------ */
        /* Closed surface from point cloud if no cells/facets */
        if (M.cells.nb() == 0 && M.facets.nb() == 0 && M.vertices.nb() > 0) {
            Logger::out("fem") << "building a closed surface from point cloud" << std::endl;
            double radius = 5.;    /* reconstruction parameter */
            double R = bbox_diagonal(M);
            mesh_repair(M, GEO::MESH_REPAIR_COLOCATE, 1e-6*R);
            radius *= 0.01 * R;
            Co3Ne_reconstruct(M, radius);

            /* Smooth remeshing of the reconstructed surface */
            index_t nb_points = M.vertices.nb() * 2;
            index_t Lloyd_iter = 5;
            index_t Newton_iter = 30;
            index_t Newton_m = 7;
            Mesh M_tmp(3); /* because output of remeshing is a different mesh*/
            remesh_smooth(M, M_tmp, nb_points, 0, Lloyd_iter, Newton_iter, Newton_m);
            M.copy(M_tmp); /* result of remeshing is copied in initial mesh */
        }

        /* Try to mesh the volume with TetGen */
        if (M.cells.nb() == 0 && M.facets.nb() > 0) {
            /* Triangulate surface mesh (useful if input is quad or polygonal mesh) */
            if (!M.facets.are_simplices()) {
                M.facets.triangulate();
            }

            /* Colocate vertices (useful if input is STL mesh) */
            Logger::out("fem") << "colocating vertices of triangulated mesh" << std::endl;
            double R = bbox_diagonal(M);
            mesh_repair(M, GEO::MESH_REPAIR_COLOCATE, 1e-6*R);

            /* TetGen call */
            Logger::out("fem") << "tetrahedralize interior of triangulated mesh" << std::endl;
            mesh_tetrahedralize(M, false, true, 0.7);
        }

        /* Verify input */
        if (!(M.cells.nb() > 0 && M.cells.are_simplices())) {
            Logger::out("fem") << "invalid input mesh" << std::endl;
            return false;
        }

        /* Ensure right adjacencies and matching between volume (cells) 
         * and boundary (facets).
         * Useful if input is tet mesh without or with bad boundaries */
        M.cells.connect();
        M.facets.clear();
        M.cells.compute_borders(); /* add facets to M on the boundary */
        M.vertices.remove_isolated(); /* avoid hanging vertices in the mesh */
        /* ------ end of mesh pre-processing ------ */

        /* ------ FEM ASSEMBLY ------ */
        const index_t n = M.vertices.nb(); /* # degrees of freedom */
        ScalarFunction is_dirichlet(dirichlet_region);
        ScalarFunction is_neumann(neumann_region);
        ScalarFunction f(sourceterm_value);
        ScalarFunction coef(diffusion_value);
        ScalarFunction dirichlet(dirichlet_value);
        ScalarFunction neumann(neumann_value);

        /* the solution is an attribute of the Geogram mesh attached to vertices */
        Attribute<double> x(M.vertices.attributes(), solution_name);
        x.fill(0.);

        /* Initialize the OpenNL linear system */
        nlNewContext();
        nlSolverParameteri(NL_NB_VARIABLES, (NLint) n);

        nlBegin(NL_SYSTEM);
        nlBegin(NL_MATRIX);

        /* Loop on cells */
        Logger::out("fem") << "assembly, loop on " << M.cells.nb() << " cells .. (" << M.vertices.nb() << " dofs)" << std::endl;
        for (index_t c = 0; c < M.cells.nb(); ++c) {
            const index_t ln = M.cells.nb_vertices(c); /* 4 local dofs per tet */
            Mat J = mapping_tet_jacobian(M, c);
            double detJ = std::abs(determinant(J));
            geo_debug_assert(detJ > 0.);
            Mat iJt = J.inverse().transpose();

            /* Loop over grad-grad quadrature points */
            std::vector<double> Ke(ln*ln, 0.); /* Elementary matrix Ke */
            for (index_t k = 0; k < tetra_o2.size() / 4; ++k ) {
                double w = tetra_o2[4*k];
                vec3 x_r = vec3(tetra_o2[4*k+1], tetra_o2[4*k+2], tetra_o2[4*k+3]);
                vec3 x = mapping_tet_transform(M, c, x_r);
                double pre = coef(x.data()) * w * detJ;
                /* Loop over local dofs */
                for (index_t li = 0; li < ln; ++li) {
                    vec3 grad_phi_i = grad_phi_tet_eval(li);
                    for (index_t lj = 0; lj < ln; ++lj) {
                        vec3 grad_phi_j = grad_phi_tet_eval(lj);
                        Ke[ln*li+lj] += pre * dot(matvec(iJt, grad_phi_i), matvec(iJt, grad_phi_j));
                    }
                }
            }

            /* Loop over rhs quadrature points */
            std::vector<double> Fe(ln, 0.); /* Elementary vector */
            for (index_t k = 0; k < tetra_o4.size() / 4; ++k ) {
                double w = tetra_o4[4*k];
                vec3 x_r = vec3(tetra_o4[4*k+1], tetra_o4[4*k+2], tetra_o4[4*k+3]);
                vec3 x = mapping_tet_transform(M, c, x_r);
                for (index_t li = 0; li < ln; ++li) {
                    Fe[li] += w * f(x.data()) * phi_tet_eval(li, x_r) * detJ;
                    // printf("%.2e, %.2e, %.2e, %.2e\n", w, f(x.data()), phi_tet_eval(li, x_r), detJ);
                }
            }

            /* Contribution from local to global matrix */
            for (index_t li = 0; li < ln; ++li) {
                for (index_t lj = 0; lj < ln; ++lj) {
                    nlAddIJCoefficient(M.cells.vertex(c,li), M.cells.vertex(c,lj),
                            Ke[ln*li+lj]);
                    // printf("mat: %i,%i <- %i,%i | +%.2e\n", M.cells.vertex(c,li),
                    //         M.cells.vertex(c,lj), li, lj, Fe[li]);
                }
                nlAddIRightHandSide(M.cells.vertex(c,li), Fe[li]);
                // printf("rhs: %i <- %i | +%.2e\n", M.cells.vertex(c,li), li, Fe[li]);
            }

        }

        /* Before loop on facets,
         * We flag facets with a constant per-facet attribute
         *  0: zero neumann (zero flux)
         *  1: dirichlet (imposed value)
         *  2: non-zero neumann (imposed normal gradient)
         * This attribute is useful for visualization / debugging
         * */
        Attribute<int> bc(M.facets.attributes(), "bc_type");
        bc.fill(0.);
        for (index_t f = 0; f < M.facets.nb(); ++f) {
            vec3 center = 1./3. * (
                      M.vertices.point(M.facets.vertex(f,0))
                    + M.vertices.point(M.facets.vertex(f,1))
                    + M.vertices.point(M.facets.vertex(f,2)));
            if (is_dirichlet(center.data()) > 0.) {
                bc[f] = 1;
            } else if (is_neumann(center.data()) > 0.) {
                bc[f] = 2;
            }
        }

        /* Apply Dirichlet boundary conditions (with the penaly method) */
        Logger::out("fem") << "assembly, enforcing dirichlet BCs ..";
        const double penalty = 100000. ;
        std::vector<bool> dirichlet_locked(M.vertices.nb(), false);
        index_t nb_dirichlet = 0;
        for (index_t f = 0; f < M.facets.nb(); ++f) {
            if (bc[f] == 1) {
                for (index_t lv = 0; lv < M.facets.nb_vertices(f); ++lv) {
                    index_t v = M.facets.vertex(f,lv);
                    if (!dirichlet_locked[v]) {
                        double value = dirichlet(M.vertices.point_ptr(v));
                        nlAddIJCoefficient(v, v, penalty);
                        nlAddIRightHandSide(v, penalty*value);
                        dirichlet_locked[v] = true;
                        nb_dirichlet += 1;
                    }
                }
            }
        }
        Logger::out("fem") << " (" << nb_dirichlet << " dirichlet dofs)" << std::endl;

        nlEnd(NL_MATRIX);
        nlEnd(NL_SYSTEM);

        /* Call the OpenNL solver */
        Logger::out("fem") << "OpenNL solver .." << std::endl;
        nlSolve();
        
        /* Copy the result from the solver to mesh vertex attribute */
        double min = DBL_MAX;
        double max = -DBL_MAX;
        for (index_t i = 0; i < M.vertices.nb(); ++i) {
            x[i] = nlGetVariable(i);
            min = geo_min(min, x[i]);
            max = geo_max(max, x[i]);
        }
        Logger::out("fem") << "solution: min=" << min << ", max=" << max << std::endl;

        nlDeleteContext(nlGetCurrent());  /* OpenNL cleanup */

        return true;
    }
}
