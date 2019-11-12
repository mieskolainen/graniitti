// PDG class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <map>
#include <regex>

// Own
#include "Graniitti/MPDG.h"

// Libraries
#include "rang.hpp"


namespace gra {
// Read PDG particle data in
// Input as the full path to PDG .mcd file
void MPDG::ReadParticleData(const std::string &filepath) {
  std::cout << "MPDG::ReadParticleData: Reading in PDG tables: ";
  std::ifstream infile(filepath);
  std::string   line;

  // Clear particle data
  PDG_table.clear();

  while (std::getline(infile, line)) {
    std::istringstream iss(line);

    // Check * (comment) lines
    std::string first;
    iss >> first;
    std::string str = first;
    str.erase(str.begin() + 1, str.end());

    // Skip lines with *
    if (str.compare("*") == 0) continue;

    //  ******************** DEFINITION FROM PDG *******************
    //       1 -  8 \ Monte Carlo particle numbers as described in
    //       the "Review of
    //       9 - 16 | Particle Physics". Charge states appear, as
    //       appropriate,
    //      17 - 24 | from left-to-right in the order -, 0, +, ++.
    //      25 - 32 /
    //           33   blank
    //      34 - 51   central value of the mass (double precision)
    //           52   blank
    //      53 - 60   positive error
    //           61   blank
    //      62 - 69   negative error
    //           70   blank
    //      71 - 88   central value of the width (double precision)
    //           89   blank
    //      90 - 97   positive error
    //           98   blank
    //      99 -106   negative error
    //          107   blank
    //     108 -128   particle name left-justified in the field and
    //                charge states right-justified in the field.
    //                This field is for ease of visual examination
    //                of the file and
    //                should not be taken as a standardized
    //                presentation of
    //                particle names.

    const int L         = 40;
    char      cid[L]    = {0};
    char      cmass[L]  = {0};
    char      cwidth[L] = {0};
    char      cname[L]  = {0};

    sscanf(line.c_str(), "%33c %35c %35c %21c", cid, cmass, cwidth, cname);
    // printf("%s\t%s\t%s\t%s \n", cid, cmass, cwidth, cname);

    std::string sid(cid, L);
    std::string smass(cmass, L);
    std::string swidth(cwidth, L);
    std::string sname(cname, L);

    std::string temp_str;

    // -----------------------------------------------------
    // ID
    // First find out how many charges
    std::vector<int>   id;
    std::istringstream ss(sid);

    while (ss >> temp_str) {
      // Convert string to int
      int value = 0;
      std::stringstream(temp_str) >> value;
      if (value != 0) { id.push_back(value); }
    }

    // -----------------------------------------------------
    // MASS
    // Mass and errors
    std::vector<double> mass;
    std::stringstream   ss2(smass);

    while (ss2 >> temp_str) {
      // Convert string to int
      double value = 0;
      std::stringstream(temp_str) >> value;
      mass.push_back(value);
    }

    // -----------------------------------------------------
    // WIDTH
    // Width and errors
    std::vector<double> width;
    std::stringstream   ss3(swidth);

    // If the particle WIDTH columns were empty
    // -> name column string will be empty by reading logic. This
    // checks that.
    bool empty = true;
    if ((sname.find("0") != std::string::npos) || (sname.find("+") != std::string::npos) ||
        (sname.find("-") != std::string::npos) || (sname.find("++") != std::string::npos)) {
      empty = false;
    }

    if (!empty) {
      while (ss3 >> temp_str) {
        // Convert string to double
        double value = 0;
        std::stringstream(temp_str) >> value;
        width.push_back(value);
      }
    }

    // -----------------------------------------------------
    // NAME and Charges
    std::stringstream ss4;
    if (empty == true) {
      ss4 = std::stringstream(swidth);
    } else {
      ss4 = std::stringstream(sname);
    }
    std::string name;
    ss4 >> name;

    // Split by comma
    std::vector<int> chargeX3;

    while (ss4.good()) {
      std::string substr;
      std::getline(ss4, substr, ',');

      // Remove extra whitespace
      substr = std::regex_replace(substr, std::regex("^ +| +$|( ) +"), "$1");

      // Identify charge (ORDER is important here!)
      if (substr.find("++") != std::string::npos) {
        chargeX3.push_back(6);
      } else if (substr.find("--") != std::string::npos) {
        chargeX3.push_back(-6);
      } else if (substr.find("-1/3") != std::string::npos) {
        chargeX3.push_back(-1);
      } else if (substr.find("+2/3") != std::string::npos) {
        chargeX3.push_back(2);
      } else if (substr.find("+") != std::string::npos) {
        chargeX3.push_back(3);
      } else if (substr.find("-") != std::string::npos) {
        chargeX3.push_back(-3);
      } else if (substr.find("0") != std::string::npos) {
        chargeX3.push_back(0);
      }
    }

    // Loop over different charge assignments
    for (std::size_t k = 0; k < id.size(); ++k) {
      // New particle
      gra::MParticle p;

      // Add properties
      p.name = name;

      p.pdg = id[k];
      if (mass.size() > 0) { p.mass = mass[0]; }
      if (width.size() > 0) {  // Unstable particles have width
        p.width = width[0];
      }
      p.chargeX3 = chargeX3[k];
      p.tau      = PDG::hbar / p.width;  // mean lifetime in the rest frame

      int lastdigit = p.pdg % 10;       // Get last digit
      p.spinX2      = (lastdigit - 1);  // Get 2J
      p.wcut        = 0;                // Off-shell mass (width) cut

      // Neutral mesons/baryons get 0 for their name
      if (std::abs(p.chargeX3) == 0 && p.pdg > 100) {
        p.name = p.name + "0";
      } else if ((p.chargeX3 != 0) &&
                 !((p.pdg >= 1 && p.pdg <= 6) ||
                   (p.pdg == 12 || p.pdg == 14 || p.pdg == 16))) {  // quarks & neutrinos

        std::string signstr = p.chargeX3 > 0 ? std::string(std::abs(p.chargeX3 / 3), '+')
                                             : std::string(std::abs(p.chargeX3 / 3), '-');

        p.name = name + signstr;
      }

      // -------------------------------------------------
      // SPECIAL CASES (not following spin numbering)

      // SM bosons
      if (p.pdg == 21) {  // gluon
        p.spinX2 = 2;
        p.P      = -1;
        p.C      = 0;  // Not defined
      }
      if (p.pdg == 22) {  // gamma
        p.spinX2 = 2;
        p.P      = -1;
        p.C      = -1;
      }
      if (std::abs(p.pdg) == 24) {  // W+-
        p.spinX2 = 2;
        p.P      = 0;  // Not defined
        p.C      = 0;  // Not defined
      }
      if (p.pdg == 23) {  // Z
        p.spinX2 = 2;
        p.P      = 0;  // Not defined
        p.C      = 0;  // Not defined
      }
      if (p.pdg == 25) {  // H
        p.spinX2 = 0;
        p.P      = 1;  // Scalar SM higgs
        p.C      = 1;  // Scalar SM higgs
      }

      // SM fermions
      if (p.pdg >= 11 && p.pdg <= 16) {  // leptons and neutrinos
        p.spinX2 = 1;
        p.P      = 1;
        p.C      = 0;  // Not defined
      }
      if (p.pdg >= 1 && p.pdg <= 6) {  // quarks
        p.spinX2 = 1;
        p.P      = 1;
        p.C      = 0;  // Not defined
      }

      // -------------------------------------------------
      // Meson J^PC (L) assingments

      if (p.spinX2 == 0 || p.spinX2 == 2 || p.spinX2 == 4 || p.spinX2 == 8 ||
          p.spinX2 == 10) {  // Is a boson

        if (std::to_string(p.pdg).length() == 3 || std::to_string(p.pdg).length() == 5 ||
            std::to_string(p.pdg).length() == 6 || std::to_string(p.pdg).length() == 7) {
          unsigned int N = std::to_string(p.pdg).length();
          unsigned int m = 0;

          if (N == 5) { m = 0; }
          if (N == 6) { m = 1; }
          if (N == 7) { m = 2; }

          // [code   JPC   L]

          // L = J-1, S = 1
          // -----------------

          if (N == 3 || (std::to_string(p.pdg)[m] == '0' && std::to_string(p.pdg)[m + 1] == '0')) {
            // 00qq3   1--   0
            if (std::to_string(p.pdg)[N - 1] == '3') p.setPCL(-1, -1, 0);

            // 00qq5   2++   1
            if (std::to_string(p.pdg)[N - 1] == '5') p.setPCL(1, 1, 1);

            // 00qq7   3--   2
            if (std::to_string(p.pdg)[N - 1] == '7') p.setPCL(-1, -1, 2);

            // 00qq9   4++   3
            if (std::to_string(p.pdg)[N - 1] == '9') p.setPCL(1, 1, 3);
          }

          // L = J, S = 0
          // -----------------

          // 00qq1   0-+   0
          if (N == 3 || (std::to_string(p.pdg)[m] == '0' && std::to_string(p.pdg)[m + 1] == '0')) {
            if (std::to_string(p.pdg)[N - 1] == '1') p.setPCL(-1, 1, 0);
          }

          if (N >= 5) {
            // 10qq3   1+-   1
            if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[N - 1] == '3')
              p.setPCL(1, -1, 1);

            // 10qq5   2-+   2
            if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[N - 1] == '5')
              p.setPCL(-1, 1, 2);

            // 10qq7   3+-   3
            if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[N - 1] == '7')
              p.setPCL(+1, -1, 3);

            // 10qq9   4-+   4
            if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[N - 1] == '9')
              p.setPCL(-1, 1, 4);

            // L = J, S = 1
            // -----------------

            // 20qq3   1++   1
            if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[N - 1] == '3')
              p.setPCL(1, 1, 1);

            // 20qq5   2--   2
            if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[N - 1] == '5')
              p.setPCL(-1, -1, 2);

            // 20qq7   3++   3
            if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[N - 1] == '7')
              p.setPCL(1, 1, 3);

            // 20qq9   4--   4
            if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[N - 1] == '9')
              p.setPCL(-1, -1, 4);

            // L = J+1, S = 1
            // -----------------

            // 10qq1   0++   1
            if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[N - 1] == '1')
              p.setPCL(1, 1, 1);

            // 30qq3   1--   2
            if (std::to_string(p.pdg)[m] == '3' && std::to_string(p.pdg)[N - 1] == '3')
              p.setPCL(-1, -1, 2);

            // 30qq5   2++   3
            if (std::to_string(p.pdg)[m] == '3' && std::to_string(p.pdg)[N - 1] == '5')
              p.setPCL(1, 1, 3);

            // 30qq7   3--   4
            if (std::to_string(p.pdg)[m] == '3' && std::to_string(p.pdg)[N - 1] == '7')
              p.setPCL(-1, -1, 4);

            // 30qq9   4++   5
            if (std::to_string(p.pdg)[m] == '3' && std::to_string(p.pdg)[N - 1] == '9')
              p.setPCL(1, 1, 5);
          }
        }
      }

      // -------------------------------------------------
      // Baryon J^PC (L) assingments

      if (p.spinX2 == 1 || p.spinX2 == 3 || p.spinX2 == 5 || p.spinX2 == 7 ||
          p.spinX2 == 9) {  // Is a fermion

        if (std::to_string(p.pdg).length() == 4 || std::to_string(p.pdg).length() == 5 ||
            std::to_string(p.pdg).length() == 6 || std::to_string(p.pdg).length() == 7) {
          // [code   J^P]
          unsigned int N = std::to_string(p.pdg).length();
          unsigned int m = 0;

          if (N == 6) { m = 0; }
          if (N == 7) { m = 1; }

          // 00qqq2
          if (N == 4 || (std::to_string(p.pdg)[m] == '0' && std::to_string(p.pdg)[m + 1] == '0')) {
            if (std::to_string(p.pdg)[N - 1] == '2') p.setPCL(1, 0, 0);
          }
          // 20qqq2
          if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[m + 1] == '0') {
            if (std::to_string(p.pdg)[N - 1] == '2') p.setPCL(1, 0, 0);
          }
          // 21qqq2
          if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[m + 1] == '1') {
            if (std::to_string(p.pdg)[N - 1] == '2') p.setPCL(1, 0, 0);
          }
          // 10qqq2
          if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[m + 1] == '0') {
            if (std::to_string(p.pdg)[N - 1] == '2') p.setPCL(-1, 0, 1);
          }

          // 00qq4
          if (N == 4 || (std::to_string(p.pdg)[m] == '0' && std::to_string(p.pdg)[m + 1] == '0')) {
            if (std::to_string(p.pdg)[N - 1] == '4') p.setPCL(1, 0, 0);
          }
          // 20qqq4
          if (std::to_string(p.pdg)[m] == '2' && std::to_string(p.pdg)[m + 1] == '0') {
            if (std::to_string(p.pdg)[N - 1] == '4') p.setPCL(1, 0, 0);
          }
          // 11qqq2
          if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[m + 1] == '1') {
            if (std::to_string(p.pdg)[N - 1] == '2') p.setPCL(1, 0, 1);
          }
          // 12qqq4
          if (std::to_string(p.pdg)[m] == '1' && std::to_string(p.pdg)[m + 1] == '2') {
            if (std::to_string(p.pdg)[N - 1] == '4') p.setPCL(-1, 0, 1);
          }

          // Check this case, PDG seems ambiguous here
          if (N == 5) { p.setPCL(1, 0, 0); }
        }
      }

      // -------------------------------------------------
      // Finally, ADD TO THE TABLE
      PDG_table.insert(std::make_pair(p.pdg, p));

      // ANTIPARTICLE (id.size < 3 because 0,+,++)
      if (((int)std::abs(p.spinX2 / 2.0) != (double)std::abs(p.spinX2 / 2.0) && id.size() < 3) ||
          (std::abs(p.chargeX3) > 0 && id.size() < 3)) {
        gra::MParticle antip = p;

        // antiquarks & antineutrinos get tilde
        if ((1 <= p.pdg && p.pdg <= 6) || (p.pdg == 12 || p.pdg == 14 || p.pdg == 16)) {
          antip.name = name + "~";
        } else if (p.chargeX3 != 0) {
          std::string signstr = p.chargeX3 > 0 ? std::string(std::abs(p.chargeX3 / 3), '-')
                                               : std::string(std::abs(p.chargeX3 / 3), '+');

          antip.name = name + signstr;
        } else if (p.chargeX3 == 0) {  // e.g. anti K0(S)
          antip.name = name + "0~";
        }

        antip.pdg      = -p.pdg;
        antip.chargeX3 = -p.chargeX3;

        // Flip the parity for fermions
        if (p.spinX2 == 1 || p.spinX2 == 3 || p.spinX2 == 5 || p.spinX2 == 7 || p.spinX2 == 9) {
          antip.P = -p.P;
        }
        // Do not flip for the bosons

        // ADD TO THE TABLE
        PDG_table.insert(std::make_pair(antip.pdg, antip));
      }
    }
  }  // PDG file line while loop

  // Process now here fast variable setup >>
  std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;
}

// Recursive function to read the decay process string
void MPDG::TokenizeProcess(const std::string &str, int depth,
                           std::vector<gra::MDecayBranch> &branches) const {
  std::string mother;
  int         aM = 0;

  for (std::string::size_type i = 0; i < str.size(); ++i) {
    // printf("[%d,%d] \n", depth, i);

    // If no subdecays left, extract particles directly
    std::string substr = str.substr(i, str.size() - i);
    if (!IsDecay(substr)) {
      std::vector<std::string> particles = gra::aux::Extract(substr);
      for (std::size_t k = 0; k < particles.size(); ++k) {
        // std::cout << "L" << depth << " " <<
        // particles[k] << std::endl;

        // **************
        // Create decay branch
        gra::MDecayBranch branch;

        branch.p     = FindByPDGName(particles[k]);
        branch.depth = depth;
        branches.push_back(branch);
        // **************
      }
      return;
    }

    // Found sub decay '>', then find closing brackets
    if (str[i] == '>') {
      mother = str.substr(aM, i - aM);

      // Check if mother is actually several particles,
      // only the last is mother, others before the last are
      // not
      // (say e+ e- rho0 > (pi+ pi-) ), mother is rho0
      std::vector<std::string> particles = gra::aux::Extract(mother);
      for (std::size_t k = 0; k < particles.size(); ++k) {
        std::cout << "L" << depth << " " << particles[k] << std::endl;

        // **************
        // Create decay branch
        gra::MDecayBranch branch;

        branch.p     = FindByPDGName(particles[k]);
        branch.depth = depth;
        branches.push_back(branch);
        // **************
      }

      mother = particles[particles.size() - 1];  // Choose the last one

      // Find brackets (  )
      std::vector<std::string::size_type> L;
      std::vector<std::string::size_type> R;

      // Now find the
      unsigned int a = 0;
      unsigned int b = 0;
      for (std::string::size_type j = i + 1; j < str.size(); ++j) {
        if (str[j] == '{') { L.push_back(j); }
        if (str[j] == '}') { R.push_back(j); }
        if ((L.size() == R.size()) && L.size() > 0) {  // Found outer closing brackets
          a = L[0];
          b = R[R.size() - 1];

          // Take the token
          std::string token = str.substr(a + 1, b - a - 1);  // pos, len
          // printf("a = %d, b = %d \n", a,b );
          TokenizeProcess(token, depth + 1, branches[branches.size() - 1].legs);  // Recursion
          i  = b + 1;
          aM = b + 1;
          break;
        }
      }
    }
  }
}

// Check if string contains '>'
bool MPDG::IsDecay(const std::string &str) const {
  bool found = false;
  for (std::string::size_type i = 0; i < str.size(); ++i) {
    if (str[i] == '>') { found = true; }
  }
  return found;
}

// Print out PDG table
void MPDG::PrintPDGTable() const {
  std::cout << "MPDG::PrintPDGTable:" << std::endl << std::endl;
  printf("\t\tpdg\tmass\t\twidth\t\tcharge\tJ^PC\tname\n");

  std::map<int, gra::MParticle>::const_iterator it = PDG_table.begin();

  unsigned int counter = 0;
  while (it != PDG_table.end()) {
    const gra::MParticle p = it->second;

    printf("%d\t%10d\t%0.6f\t%0.6f\t%2s\t%s%s%s\t%s \n", ++counter, p.pdg, p.mass, p.width,
           gra::aux::Charge3XtoString(p.chargeX3).c_str(),
           gra::aux::Spin2XtoString(p.spinX2).c_str(), gra::aux::ParityToString(p.P).c_str(),
           gra::aux::ParityToString(p.C).c_str(), p.name.c_str());
    ++it;
  }
}

// Find particle by PDG ID
const gra::MParticle &MPDG::FindByPDG(int pdgcode) const {
  std::map<int, gra::MParticle>::const_iterator it = PDG_table.find(pdgcode);
  if (it != PDG_table.end()) { return it->second; }

  // Throw a fatal error, did not found the particle
  PrintPDGTable();
  std::string str = "MProcess::FindByPDG: Unknown PDG ID: " + std::to_string(pdgcode);
  throw std::invalid_argument(str);
}

// Find particle by PDG name
const gra::MParticle &MPDG::FindByPDGName(const std::string &pdgname) const {
  std::map<int, gra::MParticle>::const_iterator it = PDG_table.begin();
  while (it != PDG_table.end()) {
    if (it->second.name.compare(pdgname) == 0) { return it->second; }
    ++it;
  }

  // Check if we have PDG number input (we allow that too)
  if (gra::aux::IsIntegerDigits(pdgname)) {
    const int pdgcode = std::stoi(pdgname);
    return FindByPDG(pdgcode);
  }

  // Throw a fatal error, did not found the particle
  PrintPDGTable();
  std::string str = "MProcess::FindByPDGName: Unknown PDG name: " + pdgname;
  throw std::invalid_argument(str);
}

}  // namespace gra
