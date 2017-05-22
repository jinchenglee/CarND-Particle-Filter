#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>

#include "particle_filter.h"

#ifdef DEBUG
void ParticleFilter::particles_print() {
	for (int i=0; i< num_particles; i++) {
	  std::cout<< "id " << particles[i].id
	     << ", x " << particles[i].x
	     << ", y " << particles[i].y
	     << ", theta " << particles[i].theta
	     << ", weight " << particles[i].weight << std::endl;
	}
}
#endif

double ParticleFilter::normpdf(double x, double mu, double std) {
    return (1.0/sqrt(2*3.1415926)/std)*exp(-0.5*((x-mu)/std)*((x-mu)/std));
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.

	// Gaussian distribution data generator
	std::default_random_engine gen;
	std::normal_distribution<double> pos_x_init(x, std[0]);
	std::normal_distribution<double> pos_y_init(y, std[1]);
	std::normal_distribution<double> pos_theta_init(theta, std[2]);

	// TODO: fine-tune number of particles.
	//   Experiments showed 300 is on the boundary. Passing sometimes, failing others.
	//   Set it to 500 to be on the safe side.
	num_particles = 500;

	Particle particle_tmp;
	for (int i=0; i< num_particles; i++) {
	     particle_tmp.id = i;
	     particle_tmp.x = pos_x_init(gen);
	     particle_tmp.y = pos_y_init(gen);
	     particle_tmp.theta = pos_theta_init(gen);
	     particle_tmp.weight = 1.0;
	     // Push generated particle into list
	     particles.push_back(particle_tmp);
	}
#ifdef DEBUG
        std::cout<<"PF init(): Init " << num_particles << " particles." << std::endl;
        particles_print();
#endif

        is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Add measurements to each particle and add random Gaussian noise.
	// Useful c++ std functions when adding noise 
    //  - std::normal_distribution and std::default_random_engine useful.
    // 
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Gaussian distribution data generator
	std::default_random_engine gen;

    // Prediction
    // x_{t+1} = x_{t} + (v/\dot{\theta}) \left [  sin(\theta_{t} + \dot{\theta}dt) - sin\theta_{t} \right ]
    // y_{t+1} = y_{t} + (v/\dot{\theta}) \left [  cos\theta_{t} - cos(\theta_{t} + \dot{\theta}dt) \right ]
    // \theta_{t+1} = \theta_{t} + \dot\theta dt
	for (int i=0; i< num_particles; i++) {

	    std::normal_distribution<double> pos_x_noise (0, std_pos[0]);
	    std::normal_distribution<double> pos_y_noise (0, std_pos[1]);
	    std::normal_distribution<double> pos_theta_noise (0, std_pos[2]);

        double tmp = yaw_rate * delta_t;
	    particles[i].x += velocity * (sin(particles[i].theta+tmp) - sin(particles[i].theta)) / yaw_rate + pos_x_noise(gen);
	    particles[i].y += velocity * (-cos(particles[i].theta+tmp) + cos(particles[i].theta)) / yaw_rate + pos_y_noise(gen);
        particles[i].theta += tmp + pos_theta_noise(gen);
	}
#ifdef DEBUG
    std::cout<<"PF prediction(): " << num_particles << " particles." << std::endl;
    particles_print();
#endif
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.

    // - implemented in a different way. No use of this function.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. 	
    //
    // The observations are given in the VEHICLE'S coordinate system while particles are located
	//   according to the MAP'S coordinate system. Need to transform between the two systems.
    //
	//   The following is a good resource for the rotation/translation theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   http://planning.cs.uiuc.edu/node99.html

    double weight_sum = 0.0;

    // Go through all particles, one after another
	for (int i=0; i< num_particles; i++) {

    
        //
        // Convert observations in vehicle coordinates to Map coordinates and only
        //  record those within senor_range. O(M) 
        //
        // Not doing the other way (convert world to robot coordinates) because that 
        //  means converting every landmark (counts of landmarks are more than observations), 
        //  inefficient compute wise. O(N)
        //
        // M observations, N landmarks. N>M.
        //
#ifdef DEBUG
        std::cout << "updateWeights(): obs_in_world conversion." << std::endl;
#endif
        Map obs_in_world; // Observation in world coordinates
    	for (int j=0; j<observations.size(); j++) {
            Map::single_landmark_s  single_landmark_tmp;
            single_landmark_tmp.id_i = j;
            single_landmark_tmp.x_f = particles[i].x + observations[j].x * cos(particles[i].theta)
                                        - observations[j].y * sin(particles[i].theta);
            single_landmark_tmp.y_f = particles[i].y + observations[j].x * sin(particles[i].theta)
                                        + observations[j].y * cos(particles[i].theta);;
            // Save into list
            obs_in_world.landmark_list.push_back(single_landmark_tmp);
#ifdef DEBUG
            std::cout << single_landmark_tmp.id_i << " " 
                    << single_landmark_tmp.x_f << " "
                    << single_landmark_tmp.y_f << std::endl;
#endif
        }



        //
        // Calculate distance between car pos (predicted) vs. landmark.
        // Maintain a list of landmarks of distance <= sensor_range. Only
        //  the landmark index is required to be recorded.
        //
        // O(N)
        //
#ifdef DEBUG
        std::cout << "updateWeights(): finding lm_idx_in_range." << std::endl;
#endif
        std::vector<int> lm_idx_in_range; // Indice ofLandmarks within sensor detection range
    	for (int j=0; j<map_landmarks.landmark_list.size(); j++) {
            double distance = dist(map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f,
                                particles[i].x, particles[i].y);
            if (distance <= sensor_range) 
                lm_idx_in_range.push_back(j);
        }
#ifdef DEBUG
    	for (int j=0; j<lm_idx_in_range.size(); j++) {
            std::cout << lm_idx_in_range[j] << " "
                    << map_landmarks.landmark_list[lm_idx_in_range[j]].x_f << " "
                    << map_landmarks.landmark_list[lm_idx_in_range[j]].y_f << std::endl;
        }
#endif



        //
        // Data association - in world coordinates
        //   J observations
        //   K landmarks within sensor range
        //   O(J*K) - loop over
        //
#ifdef DEBUG
        std::cout << "updateWeights(): dataAssociation." << std::endl;
#endif

        double min_distance = 10000.0; // The init value means nothing.
        double min_distance_x = 0.0;
        double min_distance_y = 0.0;

        for (int j=0; j<obs_in_world.landmark_list.size();j++) {
            for (int k=0; k<lm_idx_in_range.size();k++) {
                // Calculate Euclidean distance (nearest neighbour)
                double distance = dist( obs_in_world.landmark_list[j].x_f, 
                            obs_in_world.landmark_list[j].y_f,
                            map_landmarks.landmark_list[lm_idx_in_range[k]].x_f, 
                            map_landmarks.landmark_list[lm_idx_in_range[k]].y_f);
                if ( (k==0) // Start of a certain observed landmark
                    || (distance < min_distance) ) { // Better candidate found
                    min_distance = distance;
                    min_distance_x = std::fabs(map_landmarks.landmark_list[lm_idx_in_range[k]].x_f 
                                - obs_in_world.landmark_list[j].x_f);
                    min_distance_y = std::fabs(map_landmarks.landmark_list[lm_idx_in_range[k]].y_f 
                                - obs_in_world.landmark_list[j].y_f);
                }
            }
            // If obs matches lm perfectly, dist == 0. Thus expectated mean = 0 in normpdf().
            particles[i].weight *= normpdf(min_distance_x, 0.0, std_landmark[0]) * 
                            normpdf(min_distance_y, 0.0, std_landmark[1]);
#ifdef DEBUG
            std::cout.precision(17);
            std::cout << " WIP: particle " << i << " weight updated = " << std::fixed << particles[i].weight << ", x_normpdf " << normpdf(min_distance_x,0.0,std_landmark[0]) << ", y_normpdf " << normpdf(min_distance_y,0.0,std_landmark[1]) << std::endl;
#endif
        } // <end> go through each observations


        // Accumulate for next step of weight normalization
        weight_sum += particles[i].weight;

#ifdef DEBUG
        std::cout.precision(17);
        std::cout << " particle " << i << " weight updated = " << particles[i].weight << std::endl;
#endif

    } // <end> Go through each particles



    //
    // Weight normalization and update
    //
	for (int i=0; i< num_particles; i++) {
        particles[i].weight /= weight_sum;
#ifdef DEBUG
        std::cout << " particle " << i << " weight normalized = " << particles[i].weight << std::endl;
#endif
    } // <end> weight normalization

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    // Collect all the weights into a list.
    std::vector<double> weights;
	for (int i=0; i< num_particles; i++) {
        weights.push_back(particles[i].weight);
    }
    //Find the max_weight 
    double max_weight =  *std::max_element(std::begin(weights), std::end(weights));
#ifdef DEBUG
    std::cout << "Resampling(): max_weight = " << max_weight << std::endl;
#endif
    
    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    
    //Uniform random distributions for index and 2*max weight
    std::uniform_int_distribution<int> dist_num_parts(0, num_particles);
    std::uniform_real_distribution<double> weights_gen(0, 2*max_weight);
    
    std::vector<Particle> resampled_particles;

    int index = dist_num_parts(gen);
    double beta = 0.0;
    
    // Draw the wheel using Sebastian's method
    for (int i =0; i< num_particles; i++)
    {
        beta += weights_gen(gen);
        while (beta > weights[index])
        {
            beta -= weights[index];
            index = (index +1)% num_particles;
        }

        // No need to introduce noise here in newly gen'ed particle's x/y/z, as
        //   noise will be introduced anyways in prediction() step.
        resampled_particles.push_back(particles[index]);
    }

    particles = resampled_particles;
#ifdef DEBUG
    std::cout<<"Resampling(): " << num_particles << " particles." << std::endl;
    particles_print();
#endif
}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
