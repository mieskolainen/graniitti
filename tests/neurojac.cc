// GRANIITTI - Monte Carlo event generator for high energy diffraction
// https://github.com/mieskolainen/graniitti
// 
// <NeuroJacobian test>
//
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


// C++ includes
#include <iostream>
#include <future>

// Own
#include <Graniitti/MNeuroJacobian.h>
#include <Graniitti/MRandom.h>
#include <Graniitti/MAux.h>
#include <Graniitti/MH2.h>
#include <Graniitti/MMath.h>
#include <Graniitti/MRandom.h>


// Main function
int main() {

	/*
	MRandom random;
	
	gra::neurojac::MNeuroJacobian neurojac;
	gra::neurojac::BATCHSIZE = 300;

    // Set network layer dimensions [first, ..., output]
    const int D = 2;
	gra::neurojac::par.D = D; // Integrand dimension

	gra::neurojac::par.L.push_back(gra::neurojac::Layer(3,D)); // Input
	gra::neurojac::par.L.push_back(gra::neurojac::Layer(3,3));
	gra::neurojac::par.L.push_back(gra::neurojac::Layer(3,3));
	gra::neurojac::par.L.push_back(gra::neurojac::Layer(D,3)); // Output

	neurojac.Optimize();

	
	//if (ALGO == "GRAD") {
	//FLAT_TARGET = true;
	//NaiveGradDesc(w, 30, 1e-3);

	//FLAT_TARGET = false;
	//NaiveGradDesc(w, 30, 1e-3);
	//}
	

    const int N_MC  = 1e5;
    const int N_evt = 1e2;

    auto NeuroSample = [&] (std::vector<double>& u)  {

		// Prior p(z) distribution sampling
		VectorXdual z(u.size());
		for (std::size_t i = 0; i < D; ++i) {
			z[i] = random.G(0,1);
		}
	    const double p = val(gra::neurojac::gaussprob(z,0,1));
	    
		// Evaluate network map
	    VectorXdual u_ = gra::neurojac::G_net(z);

	    // Evaluate the Jacobian matrix du/dz
	    MatrixXd dudz = jacobian(gra::neurojac::G_net, u_, z);
	    
	    // Abs Jacobian determinant and inverse prior
	    const double jacweight = abs(dudz.determinant()) / p;
	    
	    // Total event weight
	    for (std::size_t i = 0; i < D; ++i) {
	    	u[i] = val(u_[i]);
	    }
	    const double weight = gra::neurojac::func(u) * jacweight;

	    return weight;
    };


	{
		MH2 hist(30, 0.0, 1.0, 30, 0.0, 1.0, "NEUROJAC");

	    uint   samples = 0;
	    double maxweight = 0.0;

	    std::vector<double> u(D);

	    // Integration
		while (true) {

			const double weight = NeuroSample(u);
		    maxweight = weight > maxweight ? weight : maxweight;

		    hist.Fill(u[0], u[1], weight);
		    ++samples;

		    if (samples == N_MC) break;
		}

		// Event generation
		samples = 0;
		uint trials = 0;
		while (true) {

			const double weight = NeuroSample(u);
		    ++trials;

		    // Acceptance-Rejection
		    if (random.U(0,1) < (weight / maxweight)) {
		    	++samples;
		    }
		    if (samples == N_evt) break;
		}
		hist.Print();
		printf("Event generation efficiency:: %0.3E \n ", N_evt / (double)trials);
	}

 	{
		MH2 hist(30, 0.0, 1.0, 30, 0.0, 1.0, "FLAT");
	    double    sumw   = 0.0;
	    double   sumw2   = 0.0;
	    uint   samples   = 0;
	    double maxweight = 0.0;
		while (true) {

			// Fully uniform sampling
			VectorXdual z(D);
			for (std::size_t i = 0; i < D; ++i) {
				z[i] = random.U(0,1);
			}
		    std::vector<double> z_ = {val(z[0]), val(z[1])};
		    const double weight = gra::neurojac::func(z_);
		    maxweight = weight > maxweight ? weight : maxweight;

		    sumw  += weight;
		    sumw2 += weight*weight;

		    hist.Fill(val(z[0]), val(z[1]), weight);
		    ++samples;
		    if (samples == N_MC) break;
		}

		samples = 0;
		uint trials = 0;
		while (true) {

			// Fully uniform sampling
			VectorXdual z(D);
			for (std::size_t i = 0; i < D; ++i) {
				z[i] = random.U(0,1);
			}
		    std::vector<double> z_ = {val(z[0]), val(z[1])};
		    const double weight = gra::neurojac::func(z_);

		    ++trials;
		    if (random.U(0,1) < (weight / maxweight)){
		    	++samples;
		    }

		    if (samples == N_MC) break;
		}
		hist.Print();
		printf("Event generation efficiency:: %0.3E \n ", N_evt / (double)trials);
	}
	*/
	
	return 0;
}
