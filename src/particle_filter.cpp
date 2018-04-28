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
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// Set num particles
	num_particles = 150;

	// Create distributions and random generator
	std::default_random_engine generator;
  std::normal_distribution<double> x_dist(x, std[0]);
	std::normal_distribution<double> y_dist(y, std[1]);
	std::normal_distribution<double> theta_dist(theta, std[2]);
	
	for (int i = 0; i < num_particles; ++i) {
		// Initialize particle
		Particle particle;
		particle.id           = i;
		particle.x            = x_dist(generator);
		particle.y            = y_dist(generator);
		particle.theta        = theta_dist(generator);
		particle.weight       = 1.0;
		// No need to initialize these variables
		//particle.sense_x      = 0;
		//particle.sense_y      = 0;
		//particle.associations = 0;

		// Angle normalization
		//while (particle.theta >  M_PI) particle.theta -= 2. * M_PI;
		//while (particle.theta < -M_PI) particle.theta += 2. * M_PI;

		// Add particle to list of particles
		particles.push_back(particle);
	}
	
	// Complete initialization
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	// Create random generator and distributions
	std::default_random_engine generator;
	std::normal_distribution<double> x_noise(0, std_pos[0]);
	std::normal_distribution<double> y_noise(0, std_pos[1]);
	std::normal_distribution<double> theta_noise(0, std_pos[2]);

	for (int i = 0; i < num_particles; ++i) {
		// Get old positions
		double x     = particles[i].x;
		double y     = particles[i].y;
		double theta = particles[i].theta;

		// Get new positions
		double x_new     = x;
		double y_new     = y;
		double theta_new = theta;
		if (abs(yaw_rate) < 1e-6) {
			x_new     = x     + velocity * delta_t * cos(theta);
			y_new     = y     + velocity * delta_t * sin(theta);
		}
		else {
			x_new     = x     + velocity/yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
			y_new     = y     + velocity/yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
			theta_new = theta + yaw_rate * delta_t;
		}

		// Update particle
		particles[i].x            = x_new     + x_noise(generator);
		particles[i].y            = y_new     + y_noise(generator);
		particles[i].theta        = theta_new + theta_noise(generator);

		// Angle normalization
		//while (particles[i].theta >  M_PI) particles[i].theta -= 2. * M_PI;
		//while (particles[i].theta < -M_PI) particles[i].theta += 2. * M_PI;
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Loop through each observed landmark
	for (unsigned int i = 0; i < observations.size(); ++i) {
		// Initialize default values for searching
		double closestLandmarkDistance = 1000000;
		int closestLandmarkID = -1;
		
		// Loop through each predicted landmark
		for (unsigned int j = 0; j < predicted.size(); ++j) {
			// Calculate distance between observed and predicted landmark
			double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			
			// If it is the closer than any other landmark...
			if (distance < closestLandmarkDistance) {
				closestLandmarkDistance = distance;
				closestLandmarkID       = predicted[j].id;
			}
		}

		// Set observation ID to closest landmark ID
		observations[i].id = closestLandmarkID;
	}

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int i = 0; i < num_particles; ++i) {
		// Get particle positions
		double x_part     = particles[i].x;
		double y_part     = particles[i].y;
		double theta_part = particles[i].theta;

		// Convert observations to map coordinates
		std::vector<LandmarkObs> translatedObservations;
		for (unsigned int j = 0; j < observations.size(); ++j) {
			LandmarkObs translatedObs;

			// Get observation positions
			double x_obs     = observations[j].x;
			double y_obs     = observations[j].y;

			// Calculate observation coordinates in map system
			translatedObs.x  = x_part + cos(theta_part)*x_obs - sin(theta_part)*y_obs;
			translatedObs.y  = y_part + sin(theta_part)*x_obs + cos(theta_part)*y_obs;
			translatedObs.id = -1;

			// Add to list
			translatedObservations.push_back(translatedObs);
		}

		cout << "num observations" << translatedObservations.size() << endl;

		// Find landmarks in sensor range
		std::vector<LandmarkObs> LandmarksInRange;
		for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
			LandmarkObs landmark;

			// Get landmark positions
			double x_lm     = map_landmarks.landmark_list[j].x_f;
			double y_lm     = map_landmarks.landmark_list[j].y_f;

			if (dist(x_part, y_part, x_lm, y_lm) < 2*sensor_range) {
				// Calculate observation coordinates in map system
				landmark.x  = x_lm;
				landmark.y  = y_lm;
				landmark.id = map_landmarks.landmark_list[j].id_i;

				// Add to list
				LandmarksInRange.push_back(landmark);
			}
		}

		cout << "num landmarks" << LandmarksInRange.size() << endl;

		// Compare observed landmarks to actual landmarks
		dataAssociation(LandmarksInRange, translatedObservations);

		// Update weights
		particles[i].weight = 1.0;
		for (unsigned int j = 0; j < translatedObservations.size(); ++j) {
			double x_obs     = observations[j].x;
			double y_obs     = observations[j].y;
			double x_lm      = 0;
			double y_lm      = 0;
			for (unsigned int k = 0; k < LandmarksInRange.size(); ++k) {
				if (LandmarksInRange[k].id == translatedObservations[j].id) {
					x_lm         = LandmarksInRange[k].x;
					y_lm         = LandmarksInRange[k].y;
				}
			}
			particles[i].weight *= 1 / exp((x_obs-x_lm)*(x_obs-x_lm)/(2*std_landmark[0]*std_landmark[0]));
			particles[i].weight *= 1 / exp((y_obs-y_lm)*(y_obs-y_lm)/(2*std_landmark[1]*std_landmark[1]));
			particles[i].weight *= 1 / (2*M_PI*std_landmark[0]*std_landmark[1]);
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	// Get vector of weights and maxWeight
	std::vector<double> weights;
	double maxWeight = 0.0;
	for (int i = 0; i < num_particles; ++i) {
		weights.push_back(particles[i].weight);
		if (particles[i].weight > maxWeight) {
			maxWeight = particles[i].weight;
		}
	}

	// Create random generator and distribution
	std::default_random_engine generator;
	std::discrete_distribution<> beta_dist(weights.begin(), weights.end());

	// Get random initial index
	uniform_int_distribution<int> index_dist(0, num_particles - 1);
	int index = index_dist(generator);

	// Initialize new particles list
	std::vector<Particle> resampledParticles;

	// Generate resampled particles list
	double beta = 0.0;
	for (int i = 0; i < num_particles; ++i) {
		beta += beta_dist(generator);
		while (beta > weights[index]) {
			beta -= weights[index];
			index += 1;
			if (index == num_particles) {
				index = 0;
			}
    }
		resampledParticles.push_back(particles[index]);
	}

	// Overwrite particle list with resampled particles
	particles = resampledParticles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

		return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
