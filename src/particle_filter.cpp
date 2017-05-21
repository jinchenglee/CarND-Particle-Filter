/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.

	// Gaussian distribution data generator
	std::default_random_engine gen;
	std::normal_distribution<double> pos_x_init(x, std[0]);
	std::normal_distribution<double> pos_y_init(y, std[1]);
	std::normal_distribution<double> pos_theta_init(theta, std[2]);

	num_particles = 3; //TODO: check whether 100 particles is enough.

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

	    std::normal_distribution<double> pos_x_noise (particles[i].x, std_pos[0]);
	    std::normal_distribution<double> pos_y_noise (particles[i].y, std_pos[1]);
	    std::normal_distribution<double> pos_theta_noise (particles[i].theta, std_pos[2]);

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
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html

    // Go through all particles, one after another
	for (int i=0; i< num_particles; i++) {
    
        // Convert *all* landmark within sensor_range in Map coordinates into 
        //   observations in vehicle coordinates 
        LandmarkObs lm_obs_tmp;
        std::vector<LandmarkObs> lm_in_car_xy;

    	for (int i=0; i< map_landmarks.size(); i++) {

            Map::single_landmark_s  single_landmark_tmp = map_landmarks[i];

            // TODO: Calculate distance between car pos (predicted) vs. landmark
            // Only continue if distance <= sensor_range
            distance = ...

            if (distance <= sensor_range) {
                //TODO: the formula is wrong below, placeholders only for now
                lm_obs_tmp.id_i = i;
                lm_obs_tmp.x = single_landmark_tmp.x_f;
                lm_obs_tmp.y = single_landmark_tmp.y_f;

                lm_in_car_xy.push_back(lm_obs_tmp);
            }
        }

        // Data association processing
        // Clean up the lm_in_car_xy list to contain only those associated landmarks

        // Weight calculation and update
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
