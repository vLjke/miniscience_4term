#include <set>
#include <cmath>
#include <gmsh.h>

namespace model = gmsh::model;

int main(int argc, char **argv) {
    gmsh::initialize(argc, argv);
    // Function to create torus
    auto createTorus = [](){
        // Clear all models and create a new one
        gmsh::clear();
        model::add("Torus");
        // Set centimeter and torus radius
        double cm = 1e-2;
        double R = 100 * cm;
        double r1 = 30 * cm;
        double r2 = 40 * cm;
        // Set center of torus
        model::occ::addPoint(0, 0, 0, 2 * cm, 1);
        // Set center of сut
        model::occ::addPoint(R, 0, 0, 2 * cm, 2);
        // Get section's center coordinates
        model::occ::synchronize();
        std::vector<double> xyz;
        model::getValue(0, 2, {}, xyz);
        // Add point to create circle arcs
        // Smaller radius circle points
        model::occ::addPoint(xyz[0] - r1, xyz[1], xyz[2], 2 * cm, 3);
        model::occ::addPoint(xyz[0], xyz[1] - r1, xyz[2], 2 * cm, 4);
        model::occ::addPoint(xyz[0] + r1, xyz[1], xyz[2], 2 * cm, 5);
        model::occ::addPoint(xyz[0], xyz[1] + r1, xyz[2], 2 * cm, 6);
        // Bigger radius circle points
        model::occ::addPoint(xyz[0] - r2, xyz[1], xyz[2], 2 * cm, 7);
        model::occ::addPoint(xyz[0], xyz[1] - r2, xyz[2], 2 * cm, 8);
        model::occ::addPoint(xyz[0] + r2, xyz[1], xyz[2], 2 * cm, 9);
        model::occ::addPoint(xyz[0], xyz[1] + r2, xyz[2], 2 * cm, 10);
        // Create smaller circle and bigger circle arcs
        for (int i = 3; i < 6; ++i) {
            model::occ::addCircleArc(i, 2, i + 1, i - 2);
            model::occ::addCircleArc(i + 4, 2, i + 5, i + 2);
        }
        model::occ::addCircleArc(6, 2, 3, 4);
        model::occ::addCircleArc(10, 2, 7, 8);
        // Create smaller circle and bigger circle
        model::occ::addCurveLoop({1, 2, 3, 4}, 1);
        model::occ::addCurveLoop({5, 6, 7, 8}, 2);
        // Create surface between circles
        model::occ::addPlaneSurface({-1, 2}, 1);
        // Revolve surface
        double angle = M_PI;
        model::getValue(0, 1, {}, xyz);
        std::vector<std::pair<int, int>> rot;
        gmsh::model::occ::revolve({{2, 1}}, xyz[0], xyz[1], xyz[2], 0, 1, 0, angle, rot);
        // Generate torus
        gmsh::model::occ::synchronize();
        gmsh::model::mesh::generate(3);
        gmsh::write("Torus.msh");
    };

    // Function to create geometry from stl file
    auto createFromStl = []() {
        // Clear all models and create a new one
        gmsh::clear();
        model::add("Leg");
        try {
            gmsh::merge("../stl_model/leg.stl");
        } catch(...) {
            gmsh::logger::write("Could not load STL mesh!");
            gmsh::finalize();
            return 0;
        }
        double angle = 40;
        bool forceParametrizablePatches = false;
        bool includeBoundary = true;
        double curveAngle = 180;
        model::mesh::classifySurfaces(angle * M_PI / 180, includeBoundary,
                                      forceParametrizablePatches, curveAngle * M_PI / 180);
        model::mesh::createGeometry();
        // Create a volume from all the surfaces
        std::vector<std::pair<int, int> > s;
        gmsh::model::getEntities(s, 2);
        std::vector<int> sl;
        for(auto surf : s) sl.push_back(surf.second);
        int l = gmsh::model::geo::addSurfaceLoop(sl);
        gmsh::model::geo::addVolume({l});
        gmsh::model::geo::synchronize();

        // We specify element sizes imposed by a size field, just because we can :-)
        int f = gmsh::model::mesh::field::add("MathEval");
        gmsh::model::mesh::field::setString(f, "F", "10");
        gmsh::model::mesh::field::setAsBackgroundMesh(f);

        gmsh::model::geo::synchronize();
        gmsh::model::mesh::generate(3);
        gmsh::write("Leg.msh");
    };

    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) {
        gmsh::fltk::run();
    }
    gmsh::finalize();
    return 0;
}