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
    createTorus();

    // Launch the GUI to see the results:
    std::set<std::string> args(argv, argv + argc);
    if(!args.count("-nopopup")) {
        gmsh::fltk::run();
    }
    gmsh::finalize();
    return 0;
}