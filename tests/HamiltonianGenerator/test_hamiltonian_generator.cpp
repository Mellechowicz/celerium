#include <hamiltonian_generator.h>

using namespace celerium;

int main(int argc, char *argv[])
{

  cubacores(1, 1);

  if (argc != 2) {
    std::cerr << "Wrong number of command line arguments.";
    return -1;
  }

  // Load pseudopotential.
  LocalPotential potential_o("../../demo/CuO2Plane/O.UPF");
  LocalPotential potential_cu("../../demo/CuO2Plane/Cu.UPF");

  Element oxygen("O",            // Element name.
                 potential_o,    // Potential.
                 {});            // Do not add any orbital classes manually.

  Element copper("Cu",            // Element name.
                 potential_cu,    // Potential.
                 {});             // Do not add any orbital classes manually.
  
  oxygen.AddOrbitalClass(2, 0);    // 2s      
  oxygen.AddOrbitalClass(2, 1);    // 2p      

  copper.AddOrbitalClass(4, 0);    // 4s
  copper.AddOrbitalClass(4, 1);    // 4p
  copper.AddOrbitalClass(3, 2);    // 3d
  
  // Define basis vectors for CuO2 plane.
  double a {3.81380};
  double c {13.22490};
  std::array<ArithmeticVector, 3> basis;
  basis[0] = {{a, 0.0, 0.0}};
  basis[1] = {{0.0, a, 0.0}};
  basis[2] = {{0.0, 0.0, 0*1000.0*a + c}};

  ElementaryCell elementary_cell(basis);
  
  elementary_cell.AddSite("Cu", copper, {{0, 0, 0}});
  elementary_cell.AddSite("O(1)", oxygen, {{a*0.5, 0.0, 0.0}});
  elementary_cell.AddSite("O(2)", oxygen, {{0.0, a*0.5, 0.0}});
  elementary_cell.AddSite("O(3)", oxygen, {{0.0, 0.0, c*0.18623}});
  elementary_cell.AddSite("O(4)", oxygen, {{0.0, 0.0, -c*0.18623}});      

  elementary_cell.SetCrystalPotentialCutoff(40.0);
  Lattice lattice(elementary_cell);
  
  // Uncomment the below region to recalculate Wannier coefficients for Cr.

  
  /*
    cuba::Cuba engine(1e7,1e6,1e-3);
    engine.parameters.epsrel = 0;
    engine.parameters.epsabs = 1e-3;
    std::vector<std::pair<double, double>> integration_limits(3, {-15, 15});

    lattice.CalculateWannierData(
      {{2, 2, 0}},            // Range of orbital overlaps (in the units of lattice parameters).
      {{6, 6, 0}},            // Maximum range of orbitals contributing to the Wannier functions.
      {{0, 0, 0}},            // Wannier located at points (0, 0, 0)  will be caluclated.
      integration_limits,     // Intergation limits for orbital overlaps.
      engine,                 // Integration engine.
      0.001,                  // Lower cutoff for the Wannier coefficients.
      true,                   // Suppress coefficiens below numerical accuracy treshold.
      true);                  // Verbose = true. Pass info about progress to std::cerr
                              // to allow for progress monitoring.
  
  lattice.SaveWannierDataToFile("wanniers.dat"); // Save Wannier data to file for later use.
  return 0;
  */
   
  lattice.LoadWannierDataFromFile("wanniers.dat");

  HamiltonianGenerator hamiltonian_generator(3.0, lattice);

  cuba::Cuba engine(1e5, 1e5,0.0);
  engine.parameters.epsrel = 0;
  engine.parameters.epsabs = 1e-2;
  engine.parameters.maxeval = 1e6;
  engine.parameters.mineval = 1e6;
  engine.parameters.flatness = 100;
  engine.parameters.maxpass = 1e3;
  
  int hamiltonian_term_index = std::atoi(argv[1]);

  if (hamiltonian_term_index == -1) {
    std::cout << hamiltonian_generator.NTerms();
    return 0;
  }

  if (hamiltonian_term_index == -2) {
    for (size_t i  = 0; i < elementary_cell.NOrbitals(); ++i) {
      std::cout << elementary_cell.GetOrbitalDescription(i) << " ";
      std::cout << elementary_cell.GetOrbitalPositionInElementaryCell(i)[0] << " ";
      std::cout << elementary_cell.GetOrbitalPositionInElementaryCell(i)[1] << " ";
      std::cout << elementary_cell.GetOrbitalPositionInElementaryCell(i)[2] << "\n";
    }
    return 0;
  }


  double value, error;
  int steps;
  hamiltonian_generator.ComputeTerm(hamiltonian_term_index,
                                    engine,
                                    10.0,
                                    value,
                                    error,
                                    steps);

  std::cout <<
      hamiltonian_generator.GetTerms()[hamiltonian_term_index].wannier1_index << " ";
  std::cout <<
      hamiltonian_generator.GetTerms()[hamiltonian_term_index].wannier2_index << " ";
  std::cout <<
      hamiltonian_generator.GetTerms()[hamiltonian_term_index].cell_distance[0] << " ";
  std::cout <<
      hamiltonian_generator.GetTerms()[hamiltonian_term_index].cell_distance[1] << " ";
  std::cout <<
      hamiltonian_generator.GetTerms()[hamiltonian_term_index].cell_distance[2] << " ";

  if (hamiltonian_generator.GetTerms()[hamiltonian_term_index].type == E0) std::cout << "E0 ";
  if (hamiltonian_generator.GetTerms()[hamiltonian_term_index].type == T) std::cout << "T ";
  if (hamiltonian_generator.GetTerms()[hamiltonian_term_index].type == U) std::cout << "U ";
  if (hamiltonian_generator.GetTerms()[hamiltonian_term_index].type == Up) std::cout << "Up ";
  if (hamiltonian_generator.GetTerms()[hamiltonian_term_index].type == J) std::cout << "J ";
  
  std::cout << value << " " << error << " " << steps << "\n";
  
  return 0;
}
