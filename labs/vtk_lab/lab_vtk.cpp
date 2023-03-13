#include <set>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

#include <gmsh.h>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

namespace model = gmsh::model;
// Node class
class CalcNode {
    friend class CalcMesh;
protected:
    // Node coordinates
    double x, y, z;
    // Node velocity
    double vx, vy, vz;
    // Node temperature
    double temp;
    // Time of evolution of system
    double evolutionTime = 0;
public:
    // Constructor
    CalcNode() = default;
    CalcNode(double x, double y, double z, double vx, double vy, double vz, double temp):
        x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), temp(temp) {}

    // Move method
    void move(double dt) {
        this->evolutionTime += dt / 2;
        // Center of rotation of the foot
        double centerY = 400;
        double centerZ = 400;
        // Rotate only ankle
        if (this->y < centerY) {
            // Angular velocity
            double omega = M_PI / 70;
            // Velocity projections on Y and Z
            this->vy = omega * (this->z - centerZ);
            this->vz = -omega * (this->y - centerY);
        }
        // Change node coordinates
        this->x += vx * dt;
        this->y += vy * dt;
        this->z += vz * dt;
    }
    // Temperature change method
    void changeTemp(double dt) {
        this->evolutionTime += dt / 2;
        // Some temperature distribution throughout the time
        this->temp = 5 * cos(this->evolutionTime * M_PI / 5 / dt + this->x * 2 * M_PI / 0.1) +
                sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2)) / 1e3 + 270;
    }
};
// Mesh element class
class Element {
    friend class CalcMesh;
protected:
    unsigned long nodesIds[4];
};
// Mesh class
class CalcMesh {
protected:
    std::vector<CalcNode> nodes;
    std::vector<Element> elements;
public:
    CalcMesh(const std::vector<double>& nodesCoords, const std::vector<std::size_t>& tetrsPoints) {
        // Create nodes
        this->nodes.resize(nodesCoords.size() / 3);
        for (int i = 0; i < this->nodes.size(); ++i) {
            double pX = nodesCoords[i * 3];
            double pY = nodesCoords[i * 3 + 1];
            double pZ = nodesCoords[i * 3 + 2];
            // Some initial temperature distribution
            double temperature = 5 * cos(pX * M_PI / 0.2) + sqrt(pow(pX, 2) + pow(pY, 2) + pow(pZ, 2)) / 1e3 + 270;
            this->nodes[i] = CalcNode(pX, pY, pZ, 0., 0., 0., temperature);
        }
        // Create mesh elements(tetrs)
        this->elements.resize(tetrsPoints.size() / 4);
        for (int i = 0; i < this->elements.size(); ++i) {
            this->elements[i].nodesIds[0] = tetrsPoints[i * 4] - 1;
            this->elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
            this->elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
            this->elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
        }
    }

    // Do time step method
    void doTimeStep(double dt) {
        for (auto& node : this->nodes) {
            node.move(dt);
            node.changeTemp(dt);
        }
    }
    // Do VTK snapshot method
    void doSnapshot(unsigned int snap_number) {
        // VTK grid
        vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        // Points in VTK grid
        vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

        // Temperature field
        auto T = vtkSmartPointer<vtkDoubleArray>::New();
        T->SetName("Temperature");

        // Velocity field
        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("Velocity");
        vel->SetNumberOfComponents(3);

        // Set grid points and points' vel and temp
        for (auto& node : this->nodes) {
            dumpPoints->InsertNextPoint(node.x, node.y, node.z);

            double _vel[3] = {node.vx, node.vy, node.vz};
            vel->InsertNextTuple(_vel);

            T->InsertNextValue(node.temp);
        }
        unstructuredGrid->SetPoints(dumpPoints);

        // Connect nodes with their temp and vel
        unstructuredGrid->GetPointData()->AddArray(vel);
        unstructuredGrid->GetPointData()->AddArray(T);

        // Create grid from points
        for(auto& element : this->elements) {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId( 0, element.nodesIds[0] );
            tetra->GetPointIds()->SetId( 1, element.nodesIds[1] );
            tetra->GetPointIds()->SetId( 2, element.nodesIds[2] );
            tetra->GetPointIds()->SetId( 3, element.nodesIds[3] );
            unstructuredGrid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        // Write snapshot in the file
        std::string fileName = "../frames_vtk/leg3d-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(unstructuredGrid);
        writer->Write();
    }
};


int main(int argc, char **argv) {
    // Mesh step
    double h = 7.;
    // Time step
    double tau = 0.1;

    const unsigned int GMSH_TETR_CODE = 4;

    gmsh::initialize();
    gmsh::model::add("leg");
    // Read file .stl
    try {
        gmsh::merge("../../gmsh_lab/stl_model/leg.stl");
    } catch(...) {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }
    // Create geometry
    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();
    // Create volume
    std::vector<std::pair<int, int> > s;
    gmsh::model::getEntities(s, 2);
    std::vector<int> sl;
    for(auto& surf : s)
        sl.push_back(surf.second);
    int l = gmsh::model::geo::addSurfaceLoop(sl);
    gmsh::model::geo::addVolume({l});

    gmsh::model::geo::synchronize();

    // Mesh step
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", std::to_string(h));
    gmsh::model::mesh::field::setAsBackgroundMesh(f);
    // Create mesh
    gmsh::model::mesh::generate(3);

    // Get nodes data from mesh
    std::vector<double> nodesCoord;
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // Get tetrahedron data from mesh
    std::vector<std::size_t>* tetrsNodesTags = nullptr;
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);
    for(unsigned int i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] != GMSH_TETR_CODE)
            continue;
        tetrsNodesTags = &elementNodeTags[i];
    }

    // Check tetra data
    if(tetrsNodesTags == nullptr) {
        cout << "Can not find tetra data. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "The model has " <<  nodeTags.size() << " nodes and " << tetrsNodesTags->size() / 4 << " tetrs." << endl;

    // Check nodes data
    for(int i = 0; i < nodeTags.size(); ++i) {
        // Indexation starts with 1, not 0
        assert(i == nodeTags[i] - 1);
    }
    assert(tetrsNodesTags->size() % 4 == 0);

    // Initialize mesh
    CalcMesh mesh(nodesCoord, *tetrsNodesTags);

    gmsh::finalize();

    // Get 100 frames to make video
    mesh.doSnapshot(0);
    unsigned int N = 100;
    for (unsigned int i = 1; i < N; ++i) {
        mesh.doTimeStep(tau);
        mesh.doSnapshot(i);
    }
    return 0;
}