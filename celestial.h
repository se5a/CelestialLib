//=====================================================================
//
//   ####  #####  ##      #####   ####  ######  ##    ###    ##
//  ##     ##     ##      ##     ##       ##    ##   ## ##   ##
//  ##     #####  ##      #####   ###     ##    ##  ##   ##  ##
//  ##     ##     ##      ##        ##    ##    ##  #######  ##
//   ####  #####  ######  #####  ####     ##    ##  ##   ##  ######
//
//=====================================================================

#ifndef CELESTIAL_C_CELESTIAL_H
#define CELESTIAL_C_CELESTIAL_H

//time_t
#include <time.h>

typedef int bool;
#define true 1
#define false 0

//
// vector
//

typedef struct vector_t *vector;

// Memory
vector vector_create();
void vector_free(vector vec);
vector vector_clone(vector vec);

//ctor
vector vector_init(double x, double y, double z);

// maths
double vector_distance_to(vector a, vector b);
void vector_abs(vector vec);
double vector_length(vector vec);
double vector_normalize(vector vec);
double vector_dot(vector a, vector b);
vector vector_cross(vector a, vector b);

// operators
vector vector_add(vector a, vector b);
vector vector_sub(vector a, vector b);
vector vector_div(vector a, vector b);
vector vector_mul(vector a, vector b);

// Util
void vector_print(vector vec);

//
// orbital_elements
//

typedef struct orbital_elements_t *orbital_elements;

typedef struct state_vectors_t *state_vectors;

//Simple malloc and free
orbital_elements orbital_elements_create();
void orbital_elements_free(orbital_elements elem);

state_vectors state_vectors_create();
void state_vectors_free(state_vectors vecs);

//so things can be pulled in from a planet
void orbital_elements_init_from_major_planet(
    orbital_elements,
    double parentMass,
    double myMass,
    double semiMajorAxis,
    double eccentricity,
    double inclination,
    double longitudeOfAscendingNode,
    double longitudeOfPeriapsis,
    double meanLongitude,
    time_t epoch);

void orbital_elements_calculate_extended_parameters(orbital_elements elem);

void orbital_elements_init_from_vector(
    orbital_elements elem,
    state_vectors state_vecs);

state_vectors orbital_elements_get_state_vectors(orbital_elements elem);


#endif //CELESTIAL_C_CELESTIAL_H
