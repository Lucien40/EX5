#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "ConfigFile.tpp"

using namespace std;

double puissance(vector<vector<double> > const& T, double const& kappa,
                 double const& h, double const& x1, double const& x2,
                 double const& y1, double const& y2);

int main(int argc, char* argv[]) {
  string inputPath("configuration.in");  // Fichier d'input par defaut
  if (argc > 1)  // Fichier d'input specifie par l'utilisateur ("./Exercice5
                 // config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath);  // Les parametres sont lus et stockes dans
                                     // une "map" de strings.

  for (int i(2); i < argc; ++i)  // Input complementaires ("./Exercice5
                                 // config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Geometrie:
  double L = configFile.get<double>("L");
  double xa = configFile.get<double>("xa");
  double xb = configFile.get<double>("xb");
  double xc = configFile.get<double>("xc");
  double xd = configFile.get<double>("xd");
  double ya = configFile.get<double>("ya");
  double yb = configFile.get<double>("yb");

  // Temperatures:
  double Tc = configFile.get<double>("Tc");
  double Tf = configFile.get<double>("Tf");
  double Tb = configFile.get<double>("Tb");

  // Physical coefficients
  double kappa = configFile.get<double>("kappa");
  double rho = configFile.get<double>("rho");
  double Ccoef = configFile.get<double>("Ccoef");
  double Dcoef = kappa / (rho * Ccoef);

  // Duree de la simulation:
  double tfin = configFile.get<double>("tfin");
  double eps =
      configFile.get<double>("eps");  // Condition d'arret si etat stationnaire

  // Discretisation:
  int N =
      configFile.get<int>("N");  // Nombre d'intervalles dans chaque dimension
  double dt = configFile.get<double>("dt");
  double h = L / N;
  double alpha = Dcoef / (h * h) * dt;

  // Fichiers de sortie:
  string output = configFile.get<string>("output");
  ofstream output_T((output + "_T.out").c_str());  // Temperature au temps final
  ofstream output_P(
      (output + "_P.out").c_str());  // Puissance au cours du temps
  output_T.precision(15);
  output_P.precision(15);

  // Tableaux:
  vector<vector<bool> > flag(N + 1, vector<bool>(N + 1));
  vector<vector<double> > T(N + 1, vector<double>(N + 1));
  vector<vector<double> > TOld(N + 1, vector<double>(N + 1));

  // TODO: Initialisation des tableaux
  //////////////////////////////////////

  for (size_t i = 0; i < flag.size(); i++) {
    for (size_t j = 0; j < flag.size(); j++) {
      if (i * h <= yb && i * h >= ya) {
        if (j * h <= xb && j * h >= xa) {
          flag[i][j] = true;
          T[i][j] = Tc;
          TOld[i][j] = Tc;
        } else if (j * h <= xd && j * h >= xc) {
          flag[i][j] = true;
          T[i][j] = Tf;
          TOld[i][j] = Tf;
        }
      } else {
        if (i == N || j == N || j == 0 || i == 0) {
          flag[i][j] = true;
        } else {
          flag[i][j] = false;
        }
        T[i][j] = Tb;
        TOld[i][j] = Tb;
      }
    }
  }

  double max(11);
  // Iterations:
  //////////////////////////////////////
  // TODO: Modifier la condition de sortie de la boucle temporelle pour tester
  // si l'etat stationnaire est atteint.
  for (int iter = 0; iter * dt < tfin && max > eps; ++iter) {
    // TODO: Schema a 2 niveaux et calcul de max(|dT/dt|)

    for (size_t i = 0; i < N + 1; i++) {
      for (size_t j = 0; j < N + 1; j++) {
        if (!flag[i][j]) {
          T[i][j] = alpha * (TOld[i + 1][j] + TOld[i - 1][j] + TOld[i][j + 1] +
                             TOld[i][j - 1] - 4 * TOld[i][j]) +
                    TOld[i][j];
        }
      }
    }

    max = (T[1][1] - TOld[1][1]) / dt;
    double deriv(max);
    for (size_t i = 0; i < N + 1; i++) {
      for (size_t j = 0; j < N + 1; j++) {
        if (!flag[i][j]) {
          deriv = (T[i][j] - TOld[i][j]) / dt;
          if (max < abs(deriv)) max = abs(deriv);
          TOld[i][j] = T[i][j];
        }
      }
    }

    // Diagnostiques:
    output_P << iter * dt << " " << puissance(T, kappa, h, xa, xb, ya, yb)
             << " " << puissance(T, kappa, h, xc, xd, ya, yb) << " "
             << puissance(T, kappa, h, xa, xd, ya, yb) << endl;
  }
  output_P.close();

  // Ecriture de la temperature finale:
  for (int i(0); i < N + 1; ++i)
    for (int j(0); j < N + 1; ++j)
      output_T << j * h << " " << i * h << " " << T[i][j] << endl;
  output_T.close();
  return 0;
}

// TODO: Calculer la puissance calorifique emise/recue par le rectangle allant
// de (x1,y1) a (x2,y2)
double puissance(vector<vector<double> > const& T, double const& kappa,
                 double const& h, double const& x1, double const& x2,
                 double const& y1, double const& y2) {
  return 0;
}
