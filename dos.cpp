#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

#include <omp.h>

const double sigma = 0.05;
const double padding_factor = 0.1;
const size_t num_points = 1000;

constexpr double PI = 3.14159265358979323844;

std::vector<double> read_eigenvalues(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return {};
  }

  std::vector<double> eigenvalues;
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    double value;
    if (iss >> value) {
      eigenvalues.push_back(value);
    } else {
      std::cerr << "Error reading eigenvalue from line: " << line << std::endl;
    }
  }

  return eigenvalues;
}

void print_dos(const std::vector<double>& eigenvalues) {
  if (eigenvalues.empty()) {
    std::cerr << "Error: Eigenvalue list is empty." << std::endl;
    return;
  }

  double min_val = eigenvalues[0];
  double max_val = eigenvalues[0];
  for (double val : eigenvalues) {
    min_val = std::min(min_val, val);
    max_val = std::max(max_val, val);
  }

  double padding = (max_val - min_val) * padding_factor;
  min_val -= padding;
  max_val += padding;

  double dE = (max_val - min_val) / (num_points - 1);
  double norm;

#if defined(LORENTZ_KERNEL)
  norm = 1.0 / PI / sigma;
#else
  norm = 1.0 / std::sqrt(2.0 * PI) / sigma;
#endif

  std::vector<double> dos_values(num_points, 0.0);
  double integral = 0.0;

#pragma omp parallel for reduction(+:integral)
  for (size_t i = 0; i < num_points; ++i) {
    double E = min_val + i * dE;
    for (double eigenvalue : eigenvalues) {
      double delta_E = (E - eigenvalue) / sigma;
#if defined(LORENTZ_KERNEL)
      dos_values[i] += 1.0 / (1.0 + delta_E * delta_E);
#else
      dos_values[i] += std::exp(-0.5 * delta_E * delta_E);
#endif
    }
    dos_values[i] *= norm;
    integral += dos_values[i] * dE;
  }

  for (size_t i = 0; i < num_points; ++i) {
    double E = min_val + i * dE;
    std::cout << E << " " << dos_values[i] / integral << std::endl;
  }
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  std::string filename = argv[1];
  std::vector<double> eigenvalues = read_eigenvalues(filename);
  if (!eigenvalues.empty()) {
    print_dos(eigenvalues);
  }

  return 0;
}
