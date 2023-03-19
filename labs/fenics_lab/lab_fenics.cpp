#include <dolfin.h>
#include "Poisson.h"
#include <mshr.h>

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression {
    void eval(Array<double>& values, const Array<double>& x) const {
        values[0] = -exp(pow(-x[0], 2) + pow(-x[1], 2));
    }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression {
    void eval(Array<double>& values, const Array<double>& x) const {
        values[0] = 0;
    }
};

// Sub domain for Dirichlet boundary condition
class Boundary : public SubDomain {
    bool inside(const Array<double>& x, bool on_boundary) const {
        return (pow(x[0], 2) + pow(x[1], 2) > 4 - DOLFIN_EPS) && on_boundary;
    }
};

int main() {
    // Create mesh and function space
    auto circle = std::make_shared<mshr::Circle> (Point(0., 0.), 2.2, 32);
    auto mesh = mshr::generate_mesh(circle, 100);
    auto V = std::make_shared<Poisson::FunctionSpace>(mesh);
    // Define boundary condition
    auto u_b = std::make_shared<Constant>(0.0);
    auto boundary = std::make_shared<Boundary>();
    DirichletBC bc(V, u_b, boundary);

    // Define variational forms
    Poisson::BilinearForm a(V, V);
    Poisson::LinearForm L(V);
    auto f = std::make_shared<Source>();
    auto g = std::make_shared<dUdN>();
    L.f = f;
    L.g = g;

    // Compute solution
    Function u(V);
    solve(a == L, u, bc);

    // Save solution in VTK format
    File file("poisson_test.pvd");
    file << u;

    return 0;
}
