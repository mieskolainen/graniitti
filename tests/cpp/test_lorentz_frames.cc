// Unit tests for Lorentz frame transformations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Kinematics/M4Vec.h"

using namespace gra;

using gra::math::pow2;


TEST_CASE("generic Lorentz frames preserve the physical two-body state",
          "[Lorentz Frames][Kinematics]") {
  constexpr unsigned int trials = 128;
  constexpr double tolerance = 1e-8;
  constexpr double proton_mass = 0.938;
  constexpr double pion_mass = 0.14;

  MRandom rng;
  rng.SetSeed(123456);
  std::vector<M4Vec> final(2);

  for (unsigned int trial = 0; trial < trials; ++trial) {
    const double mother_mass = rng.U(2.0 * pion_mass + 0.01, 100.0);
    M4Vec mother;
    mother.SetPxPyPzM(rng.U(-100.0, 100.0), rng.U(-100.0, 100.0),
                      rng.U(-100.0, 100.0), mother_mass);

    const auto weight = gra::kinematics::TwoBodyPhaseSpace(
        mother, mother_mass, {pion_mass, pion_mass}, final, rng);
    REQUIRE(weight.GetW() > 0.0);
    REQUIRE(gra::math::CheckEMC(mother - final[0] - final[1], tolerance));
    REQUIRE(final[0].M2() == Approx(pion_mass * pion_mass).margin(tolerance));
    REQUIRE(final[1].M2() == Approx(pion_mass * pion_mass).margin(tolerance));

    const double beam_pz = rng.U(50.0, 7000.0);
    M4Vec beam1;
    M4Vec beam2;
    beam1.SetPxPyPzM(0.0, 0.0, beam_pz, proton_mass);
    beam2.SetPxPyPzM(0.0, 0.0, -beam_pz, proton_mass);

    const M4Vec system = final[0] + final[1];
    M4Vec beam1_rest;
    M4Vec beam2_rest;
    std::vector<M4Vec> final_rest;
    gra::kinematics::LorentFramePrepare(
        final, system, beam1, beam2, beam1_rest, beam2_rest, final_rest);
    REQUIRE((final_rest[0] + final_rest[1]).P3mod() ==
            Approx(0.0).margin(tolerance));
    REQUIRE((final_rest[0] + final_rest[1]).E() ==
            Approx(mother_mass).margin(tolerance));

    for (const std::string frame : {"CS", "HX", "AH", "PG", "CM"}) {
      CAPTURE(trial, frame);
      std::vector<M4Vec> transformed;
      gra::kinematics::LorentzFrame(
          transformed, beam1_rest, beam2_rest, final_rest, frame, -1);
      REQUIRE(transformed.size() == final.size());
      REQUIRE((transformed[0] + transformed[1]).P3mod() ==
              Approx(0.0).margin(tolerance));
      REQUIRE((transformed[0] + transformed[1]).E() ==
              Approx(mother_mass).margin(tolerance));

      for (std::size_t i = 0; i < transformed.size(); ++i) {
        REQUIRE(transformed[i].M2() ==
                Approx(final[i].M2()).margin(tolerance));
        REQUIRE(transformed[i].E() ==
                Approx(final_rest[i].E()).margin(tolerance));
        REQUIRE(transformed[i].P3mod2() ==
                Approx(final_rest[i].P3mod2()).margin(tolerance));
      }

      std::vector<M4Vec> standalone = final;
      if (frame == "CS") {
        gra::kinematics::CSframe(standalone, system, beam1, beam2);
      } else if (frame == "HX") {
        gra::kinematics::HXframe(standalone, system);
      } else if (frame == "AH") {
        gra::kinematics::AHframe(standalone, system, beam1, beam2);
      } else if (frame == "PG") {
        gra::kinematics::PGframe(standalone, system, -1, beam1, beam2);
      } else {
        standalone = final_rest;
      }

      for (std::size_t i = 0; i < transformed.size(); ++i) {
        REQUIRE(gra::math::CheckEMC(
            transformed[i] - standalone[i], tolerance));
      }
    }
  }
}

TEST_CASE("GJframe aligns the selected exchange in the central rest frame",
          "[Lorentz Frames][GJ][Kinematics]") {
  constexpr unsigned int trials = 128;
  constexpr double tolerance = 1e-8;
  constexpr double proton_mass = 0.938;
  constexpr double pion_mass = 0.14;

  MRandom rng;
  rng.SetSeed(789012);
  std::vector<M4Vec> final(4);

  for (unsigned int trial = 0; trial < trials; ++trial) {
    const double beam_pz = rng.U(50.0, 700.0);
    M4Vec beam1;
    M4Vec beam2;
    beam1.SetPxPyPzM(0.0, 0.0, beam_pz, proton_mass);
    beam2.SetPxPyPzM(0.0, 0.0, -beam_pz, proton_mass);
    const M4Vec initial = beam1 + beam2;

    const auto weight = gra::kinematics::NBodyPhaseSpace(
        initial, initial.M(), {proton_mass, pion_mass, pion_mass, proton_mass},
        final, false, rng);
    REQUIRE(weight.GetW() > 0.0);
    M4Vec final_sum;
    for (const auto& momentum : final) {
      final_sum += momentum;
    }
    REQUIRE(gra::math::CheckEMC(initial - final_sum, tolerance));

    const M4Vec system = final[1] + final[2];
    const M4Vec exchange1 = beam1 - final[0];
    const M4Vec exchange2 = beam2 - final[3];
    REQUIRE(gra::math::CheckEMC(
        system - exchange1 - exchange2, tolerance));

    for (const int direction : {-1, 1}) {
      CAPTURE(trial, direction);
      std::vector<M4Vec> transformed = {exchange1, exchange2};
      gra::kinematics::GJframe(
          transformed, system, direction, exchange1, exchange2);

      REQUIRE(transformed[0].M2() ==
              Approx(exchange1.M2()).margin(tolerance));
      REQUIRE(transformed[1].M2() ==
              Approx(exchange2.M2()).margin(tolerance));
      REQUIRE((transformed[0] + transformed[1]).P3mod() ==
              Approx(0.0).margin(tolerance));
      REQUIRE((transformed[0] + transformed[1]).E() ==
              Approx(system.M()).margin(tolerance));
      REQUIRE(transformed[0].Pz() + transformed[1].Pz() ==
              Approx(0.0).margin(tolerance));

      const std::size_t selected = direction == -1 ? 0 : 1;
      REQUIRE(transformed[selected].Px() ==
              Approx(0.0).margin(tolerance));
      REQUIRE(transformed[selected].Py() ==
              Approx(0.0).margin(tolerance));
    }
  }
}

TEST_CASE("gra::kinematics frame helpers reject invalid directions", "[Lorentz Frames]") {
	const bool DEBUG = false;
	const double mproton = 0.938;
	const double mpion = 0.14;

	MRandom rng;
	rng.SetSeed(123456);

	const double PZ = 6500.0;
	const M4Vec p1(0, 0, PZ, std::sqrt(pow2(mproton) + pow2(PZ)));
	const M4Vec p2(0, 0, -PZ, std::sqrt(pow2(mproton) + pow2(PZ)));

	std::vector<M4Vec> pf4(4);
	const auto weight = gra::kinematics::NBodyPhaseSpace(
	    p1 + p2, (p1 + p2).M(), {mproton, mpion, mpion, mproton}, pf4,
	    false, rng);
	REQUIRE(weight.GetW() > 0.0);
	M4Vec final_sum;
	for (const auto &momentum : pf4) { final_sum += momentum; }
	REQUIRE(gra::math::CheckEMC(p1 + p2 - final_sum, 1e-8));

	const M4Vec X  = pf4[1] + pf4[2];
	const M4Vec q1 = p1 - pf4[0];
	const M4Vec q2 = p2 - pf4[3];

	SECTION("GJframe throws on invalid direction") {
		std::vector<M4Vec> out = {q1, q2};
		REQUIRE_THROWS_AS(gra::kinematics::GJframe(out, X, 0, q1, q2, DEBUG),
		                  std::invalid_argument);
	}

	SECTION("PGframe throws on invalid direction") {
		std::vector<M4Vec> out = {pf4[1], pf4[2]};
		REQUIRE_THROWS_AS(gra::kinematics::PGframe(out, X, 0, p1, p2, DEBUG),
		                  std::invalid_argument);
	}

	SECTION("LorentzFrame throws on invalid PG direction") {
		M4Vec pb1boost;
		M4Vec pb2boost;
		std::vector<M4Vec> pfboost;
		gra::kinematics::LorentFramePrepare({pf4[1], pf4[2]}, X, p1, p2, pb1boost, pb2boost, pfboost);

		std::vector<M4Vec> out;
		REQUIRE_THROWS_AS(
		    gra::kinematics::LorentzFrame(out, pb1boost, pb2boost, pfboost, "PG", 0),
		    std::invalid_argument);
	}
}

// Require a finite four-vector with invariant mass preserved by a frame rotation
void RequireFiniteFrameVector(const M4Vec &input, const M4Vec &output) {
	REQUIRE(std::isfinite(output.Px()));
	REQUIRE(std::isfinite(output.Py()));
	REQUIRE(std::isfinite(output.Pz()));
	REQUIRE(std::isfinite(output.E()));
	REQUIRE(output.M2() == Approx(input.M2()).margin(1e-12));
	REQUIRE(output.P3mod2() == Approx(input.P3mod2()).margin(1e-12));
}

TEST_CASE("Lorentz frame axes are deterministic for exactly collinear beams",
          "[Lorentz Frames][Collinear]") {
	const std::vector<M4Vec> particles = {
	    M4Vec(0.7, -0.2, 1.1, 1.6), M4Vec(-0.4, 0.5, -0.8, 1.3)};
	const std::vector<std::string> frames = {"CS", "AH", "HX", "PG", "CM"};
	const std::vector<std::pair<M4Vec, M4Vec>> beam_pairs = {
	    {M4Vec(0.0, 0.0, 5.0, 5.2), M4Vec(0.0, 0.0, -4.0, 4.3)},
	    {M4Vec(0.0, 0.0, 5.0, 5.2), M4Vec(0.0, 0.0, 4.0, 4.3)}};

	for (const auto &beams : beam_pairs) {
		for (const auto &frame : frames) {
			std::vector<M4Vec> first;
			std::vector<M4Vec> second;
			gra::kinematics::LorentzFrame(
			    first, beams.first, beams.second, particles, frame, -1);
			gra::kinematics::LorentzFrame(
			    second, beams.first, beams.second, particles, frame, -1);
			REQUIRE(first.size() == particles.size());
			for (std::size_t i = 0; i < particles.size(); ++i) {
				RequireFiniteFrameVector(particles[i], first[i]);
				REQUIRE(first[i].Px() == Approx(second[i].Px()).margin(1e-15));
				REQUIRE(first[i].Py() == Approx(second[i].Py()).margin(1e-15));
				REQUIRE(first[i].Pz() == Approx(second[i].Pz()).margin(1e-15));
			}
		}
	}
}

TEST_CASE("standalone Collins-Soper and anti-helicity frames handle the collinear limit",
          "[Lorentz Frames][Collinear]") {
	const double beam_energy = 10.0;
	const M4Vec beam1(0.0, 0.0, 9.0, beam_energy);
	const M4Vec beam2(0.0, 0.0, -9.0, beam_energy);
	const M4Vec system(0.0, 0.0, 0.0, 4.0);
	const std::vector<M4Vec> particles = {
	    M4Vec(0.3, 0.4, 1.2, 2.0), M4Vec(-0.3, -0.4, -1.2, 2.0)};

	for (const std::string frame : {"CS", "AH"}) {
		std::vector<M4Vec> generic;
		gra::kinematics::LorentzFrame(generic, beam1, beam2, particles, frame, -1);
		std::vector<M4Vec> standalone = particles;
		if (frame == "CS") {
			gra::kinematics::CSframe(standalone, system, beam1, beam2);
		} else {
			gra::kinematics::AHframe(standalone, system, beam1, beam2);
		}
		for (std::size_t i = 0; i < particles.size(); ++i) {
			RequireFiniteFrameVector(particles[i], standalone[i]);
			REQUIRE(generic[i].Px() == Approx(standalone[i].Px()).margin(1e-12));
			REQUIRE(generic[i].Py() == Approx(standalone[i].Py()).margin(1e-12));
			REQUIRE(generic[i].Pz() == Approx(standalone[i].Pz()).margin(1e-12));
		}
	}
}
