//============================================================================
// Name        : Dogfight.cpp
// Author      : Jacob Chisausky
// Version     :
// Copyright   : Your copyright notice
// Description : Dogfight game simulation
//============================================================================

//Paper notes: 'adversary' and 'pilot' in bold are confusing.

//Questions: What if US is not sure about eavesdropping?
	//c < 1.0

#include <iostream>
#include <vector>
#include <random>
#include "rndutils.hpp"

class US {
public:
	//double Pursue; //Chance of pursuing. Otherwise return
	double hP; //P if 'hit'
	double mP; //P if 'miss'

	double k;  //say 'miss' if hit
	double j;  //say 'hit' if miss

	double fitness = 0.0;

	US(double hP_, double mP_, double k_, double j_){
		hP = hP_;
		mP = mP_;
		k = k_;
		j = j_;
	}
};

class Sov {
public:
	//double Engage; //Chance of Engaging. Otherwise flee
	double hE;	// engage if 'hit'
	double mE;  // engage if 'miss'
	double nE;	// engage if no signal (which happens with chance 1-c)

	double fitness = 0.0;

	Sov(double hE_, double mE_, double nE_){
		hE = hE_;
		mE = mE_;
		nE = nE_;
	}
};


int main() {

	//Parameters
	double p = .5; //Chance of US plane being hit
	double c = 1.0;	//Chance of eavesdropping occurring

	double a = 0.5;
	double b = 0.25;

	int seed = 123456789;

	int G = 100000;
	int N = 1000;
	int partners = 5; //how many times each player plays the game with a different partner (US vs Sov)
	int replicates = 1;

	double mu_US = 0.01; //hP and mP
	double mu_k = 0.0;
	double mu_j = 0.0;
	double mu_Sov = 0.00; //hE and mE and nE. Make sure to ignore nE if c = 1

	double muSize = 0.01; //For now- just one mutation size. SD of normal dist with mean 0

	double init_hP = 0.5;
	double init_mP = 0.5;
	double init_k = 0.0;  //hit, say 'miss'
	double init_j = 0.0;  //miss, say 'hit

	double init_hE = 1.0;
	double init_mE = 0.0;
	double init_nE = 0.5;

	int reportFreq = 1000;	//Report data every reportFreq generations

	int k = 2; //For k-selection tournament

	//Payoff: parenthesis show values in text
	//If US plane NOT hit: US payoff
	//				Soviet
	//	|		|Engage	|Flee
	//US|Pursue	| 1		| 0
	//	|Return	| -b	| 0

	//If US plane IS hit: US payoff
	//				Soviet
	//	|		|Engage	|Flee
	//US|Pursue	| -a 	| 0
	//	|Return	| -b	| 0


	//Vector used for drawing random numbers without replacement
	std::vector<int> nullVec;
	for (int i = 0; i < N; i++){
		nullVec.push_back(i);
	}

	auto rng = std::default_random_engine {seed};
	std::uniform_real_distribution<double> prob(0,1);
	std::uniform_int_distribution<int> randN(0,N-1);

	//Mutation distributions
	auto muStepDist = std::normal_distribution<double>(0.0, std::abs(muSize));
	std::bernoulli_distribution mu_US_dist(mu_US);
	std::bernoulli_distribution mu_k_dist(mu_k);
	std::bernoulli_distribution mu_j_dist(mu_j);
	std::bernoulli_distribution mu_Sov_dist(mu_Sov);



	//Start replicate loop
	for (int rep = 1; rep <= replicates; rep++){


		//Initialize population
		std::vector<US> US_vec;
		std::vector<Sov> Sov_vec;

		std::vector<US> US_vec_offspring;
		std::vector<Sov> Sov_vec_offspring;

		//First initialize population
		for (int i = 0; i < N; i++){
			US_vec.push_back(US(init_hP,init_mP,init_k,init_j));
			Sov_vec.push_back(Sov(init_hE,init_mE,init_nE));
			US_vec_offspring.push_back(US(init_hP,init_mP,init_k,init_j));
			Sov_vec_offspring.push_back(Sov(init_hE,init_mE,init_nE));
		}

		for (int g = 1; g <= G; g++){

			for (int game = 0; game < partners; game++){//Playing multiple games with different partners

				//Shuffle population vectors. Pair US 1 with Sov 1, US 2 with Sov 2, etc...
				std::shuffle(std::begin(US_vec), std::end(US_vec), rng);
				std::shuffle(std::begin(Sov_vec), std::end(Sov_vec), rng);

				//Play the game for each n pair of players
				for (int n = 0; n < N; n++){

					bool hit = (prob(rng) < p);	//US plane hit or not?

					bool signal = hit; //0 = miss. 1 = hit

					if (hit == 1){
						//say hit, or lie with prob k
						if (prob(rng) < US_vec[n].k){
							signal = 0;
						}
					} else {
						//say miss, or lie with prob j
						if (prob(rng) < US_vec[n].j){
							signal = 1;
						}
					}
					//Signal sent by controller is determined now

					//does Sov see signal or no?
					bool SovSignal = 1; //Yes, signal seen
					if (prob(rng) > c){
						SovSignal = 0;  //No, signal not seen
					}

					//US and Sov actions based on signal
					int US_Action = 1;   //1 = pursue, 0 = return
					int Sov_Action = 1;  //1 = engage, 0 = flee

					if (signal == 1){ //Signal = 1 ('hit')
						if (prob(rng) > US_vec[n].hP){
							US_Action = 0; //US returns
						}
						if (SovSignal == 1){//Sov sees signal
							if (prob(rng) > Sov_vec[n].hE){
								Sov_Action = 0;
							}
						} else { //Sov doesn't see signal
							if (prob(rng) > Sov_vec[n].nE){
								Sov_Action = 0;
							}
						}
					} else { //Signal = 0 ('miss')
						if (prob(rng) > US_vec[n].mP){
							US_Action = 0; //US returns
						}
						if (SovSignal == 1){//Sov sees signal
							if (prob(rng) > Sov_vec[n].mE){
								Sov_Action = 0;
							}
						} else { //Sov doesn't see signal
							if (prob(rng) > Sov_vec[n].nE){
								Sov_Action = 0;
							}
						}
					}

				//	std::cout<<hit<<" "<<signal<<" "<<US_Action<<" "<<Sov_Action<<"\n";

					//US and Sov actions are decided. Determine Payoffs
					double US_Payoff = 0.0;
					if (hit == 0){ //hit = 0, miss
						if (US_Action == 0){
							if (Sov_Action == 0){
								//miss, 0, 0
								//miss, return, flee
								US_Payoff = 0.0;
							} else {
								//miss, 0, 1
								//miss, return, engage
								US_Payoff = -b;
							}
						} else {
							if (Sov_Action == 0){
								//miss, 1, 0
								//miss, pursue, flee
								US_Payoff = 0.0;
							} else {
								//miss, 1, 1
								//miss, pursue, engage
								US_Payoff = 1.0;
							}
						}
					} else { //hit = 1, hit
						if (US_Action == 0){
							if (Sov_Action == 0){
								//hit, return, flee
								US_Payoff = 0.0;
							} else {
								//hit, return, engage
								US_Payoff = -b;
							}
						} else {
							if (Sov_Action == 0){
								//hit, pursue, flee
								US_Payoff = 0.0;

							} else {
								//hit, pursue, engage
								US_Payoff = -a;
							}
						}
					}
					//Payoffs are now determined.
					//Add fitnesses to agents
					US_vec[n].fitness += US_Payoff;
					Sov_vec[n].fitness -= US_Payoff;

				}
			}

			//All games have been played. Fitnesses are accumulated for all players.
			//Reproduction step.

			//k-selection tournament
			for (int n = 0; n < N; n++){

				int id = randN(rng);
				int id_best = id;
				double maxFit = US_vec[id].fitness;
				for (int i = 0; i < k; i++){
					id = randN(rng);
					if (US_vec[id].fitness > maxFit){
						id_best = id;
					}
				}
				US_vec_offspring[n] = US_vec[id];
				US_vec_offspring[n].fitness = 0.0;
				//Mutation of US offspring
				//mutation to hP
				if (mu_US_dist(rng)){
					double mut = muStepDist(rng);
					US_vec_offspring[n].hP = std::min(std::max(US_vec_offspring[n].hP + mut, 0.0),1.0);
				}

				//mutation to mP
				if (mu_US_dist(rng)){
					double mut = muStepDist(rng);
					US_vec_offspring[n].mP = std::min(std::max(US_vec_offspring[n].mP + mut, 0.0),1.0);
				}

				//mutation to j
				if (mu_j_dist(rng)){
					double mut = muStepDist(rng);
					US_vec_offspring[n].j = std::min(std::max(US_vec_offspring[n].j + mut, 0.0),1.0);
				}

				//mutation to k
				if (mu_k_dist(rng)){
					double mut = muStepDist(rng);
					US_vec_offspring[n].k = std::min(std::max(US_vec_offspring[n].k + mut, 0.0),1.0);
				}


				id = randN(rng);
				id_best = id;
				maxFit = Sov_vec[id].fitness;
				for (int i = 0; i < k; i++){
					id = randN(rng);
					if (Sov_vec[id].fitness > maxFit){
						id_best = id;
					}
				}
				Sov_vec_offspring[n] = Sov_vec[id];
				Sov_vec_offspring[n].fitness = 0.0;
				//Mutation of Sov offspring

				//mutation to hE
				if (c > 0.0){ // only if chance of eavesdropping is > 0
					if (mu_Sov_dist(rng)){
						double mut = muStepDist(rng);
						Sov_vec_offspring[n].hE = std::min(std::max(Sov_vec_offspring[n].hE + mut, 0.0),1.0);
					}

					//mutation to mE
					if (mu_Sov_dist(rng)){
						double mut = muStepDist(rng);
						Sov_vec_offspring[n].mE = std::min(std::max(Sov_vec_offspring[n].mE + mut, 0.0),1.0);
					}
				}

				//mutation to nE
				if (c < 1.0){ //only if chance of eavesdropping is < 1
					if (mu_Sov_dist(rng)){
						double mut = muStepDist(rng);
						Sov_vec_offspring[n].nE = std::min(std::max(Sov_vec_offspring[n].nE + mut, 0.0),1.0);
					}
				}

			}

			//Reproduction and mutation is complete.

			//Report data
			if (g % reportFreq == 0 | g == G){
				long double sum_hP = 0.0;
				long double sum_mP = 0.0;
				long double sum_j = 0.0;
				long double sum_k = 0.0;
				long double sum_hE = 0.0;
				long double sum_mE = 0.0;

				for (int i = 0; i < N; i++){
					//std::cout << g << "  hP: "<< US_vec[i].hP<< " mP: "<< US_vec[i].mP << " j: "<< US_vec[i].j<< " k: " << US_vec[i].k << "\n";
					sum_hP += US_vec[i].hP;
					sum_mP += US_vec[i].mP;
					sum_j += US_vec[i].j;
					sum_k += US_vec[i].k;

					sum_hE += Sov_vec[i].hE;
					sum_mE += Sov_vec[i].mE;
				}

				std::cout << g << "  hP: " << sum_hP/double(N) << "  mP: " << sum_mP/double(N) << "  j: " << sum_j/double(N) << "  k: " << sum_k/double(N) << "  hE: " << sum_hE/double(N) << "  mE: " << sum_mE/double(N) << "\n";


			}


			//Replace Parents with Offspring
			US_vec.swap(US_vec_offspring);
			Sov_vec.swap(Sov_vec_offspring);


		}


	}
	return 0;
}
