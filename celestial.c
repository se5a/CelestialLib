//=====================================================================
//
//   ####  #####  ##      #####   ####  ######  ##    ###    ##
//  ##     ##     ##      ##     ##       ##    ##   ## ##   ##
//  ##     #####  ##      #####   ###     ##    ##  ##   ##  ##
//  ##     ##     ##      ##        ##    ##    ##  #######  ##
//   ####  #####  ######  #####  ####     ##    ##  ##   ##  ######
//
//=====================================================================

/*
    -Myvar:
    Coding conventions:
    1. Please add the thing your using in an include, with a comment
            ->
                //malloc, free
                #include <stdlib.h>

    2. Use Comment with one empty one above an one empty below
            -> 
                //
                // Some Comment
                //

    3. Use Sections to make code navagation easyer
            -> 
                //
                // [Helpers]
                //
                inline double rad_to_deg(double rad) { ... }

    4. use @ tags to mark things: (ensure tag index is upto date)
            -> 
                // @Hack !!! this code is hacky
                // @Cleanup Code is messy the names need to be updated
                // @Speed !! malloc a big block insted of meany small malocs for cacheing
  
    5. prefix all procedures with a structs name if its related to that struct
            -> void orbital_elements_init(orbital_elements);
                insted of just init

 */

//
// Tag Index
//

/*
    * @Todo
    * @Speed
    * @Cleanup
    * @Refactoring
    * @Hack
    * @bug
 */

//
// Todo/Ideas
//

/*
    * Orbital Mechanics
    * Stellar classification
*/

//malloc, free
#include <stdlib.h>

//printf
#include <stdio.h>

//time_t
#include <time.h>

//M_PI, sqrt, abs, round
#include <math.h>

#include "celestial.h"

#ifdef TEST

//assert
#include <assert.h>

#endif

//
// [Implmentation]
//

//
// [Helpers]
//
static inline double rad_to_deg(double rad)
{
    // @Cleanup: vs code things M_PI does not exists but it does wtf
    return rad * (180.0 / M_PI);
}

static inline double deg_to_rad(double deg)
{
    // @Cleanup: vs code things M_PI does not exists but it does wtf
    return deg * (M_PI / 180.0);
}

// @Refactoring, @Cleanup : this can easly be confused with sqrT maby rename it
static inline double sqr(double n)
{
    return n * n;
}

//
// [Math]
//

//
// [Vector]
//

struct vector_t
{
    double x;
    double y;
    double z;
};

vector vector_create()
{
    return malloc(sizeof(struct vector_t));
}

void vector_free(vector vec)
{
    free(vec);
}

vector clone(vector vec)
{
    vector a = vector_create();
    a->x = vec->x;
    a->y = vec->y;
    a->z = vec->z;

    return a;
}

vector vector_init(double x, double y, double z)
{
    vector vec = vector_create();

    vec->x = x;
    vec->y = y;
    vec->z = z;

    return vec;
}

// maths
double vector_distance_to(vector a, vector b)
{
    double x = a->x - b->x;
    double y = a->y - b->y;
    double z = a->z - b->z;

    return sqrt(x * x + y * y + z * z);
}

void vector_abs(vector vec)
{
    vec->x = abs(vec->x);
    vec->y = abs(vec->y);
    vec->z = abs(vec->z);
}

double vector_length(vector vec)
{
    //
    return sqrt(sqr(vec->x) + sqr(vec->y) + sqr(vec->z));
}

double vector_normalize(vector vec)
{
    double l = vector_length(vec);

    vec->x /= l;
    vec->y /= l;
    vec->z /= l;

    return l;
}

double vector_dot(vector a, vector b)
{
    return a->x * b->x + a->y * b->y + a->z * b->z;
}

vector vector_cross(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->y * b->z - a->z * b->y;
    vec->y = a->z * b->x - a->x * b->z;
    vec->z = a->x * b->y - a->y * b->x;

    return vec;
}

// operators
vector vector_add(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x + b->x;
    vec->y = a->y + b->y;
    vec->z = a->z + b->z;

    return vec;
}

vector vector_sub(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x - b->x;
    vec->y = a->y - b->y;
    vec->z = a->z - b->z;

    return vec;
}

vector vector_div(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x / b->x;
    vec->y = a->y / b->y;
    vec->z = a->z / b->z;

    return vec;
}

vector vector_mul(vector a, vector b)
{
    vector vec = vector_create();

    vec->x = a->x * b->x;
    vec->y = a->y * b->y;
    vec->z = a->z * b->z;

    return vec;
}

void vector_print(vector vec)
{
    printf("[X: %f, Y: %f, Z: %f]\n", vec->x, vec->y, vec->z);
}

//
// [Orbital_Mechanics]
//
struct orbital_elements_t
{
    double semi_major_axis_km;
    double eccentricity;
    double inclination_km;
    double longitude_of_ascending_node_rad;
    double argument_of_periapsis_rad;
    double mean_anomaly_rad;
    time_t epoch;
    double gravitational_parameter;
    time_t orbital_period;
    long mean_motion_rad;
    double apoapsis_km;
    double periapsis_km;
    double mass_kg;
    double parent_mass_kg;
};

// @Speed ! mallocing a big struct like this piecemeal is slow, it would be way faster to use an Array for Struct, with block mallocs
orbital_elements orbital_elements_create()
{
    return malloc(sizeof(struct orbital_elements_t));
}

void orbital_elements_free(orbital_elements elem) { free(elem); }

void orbital_elements_init_from_major_planet(
    orbital_elements elem,
    double parentMass,
    double myMass,
    double semiMajorAxis,
    double eccentricity,
    double inclination,
    double longitudeOfAscendingNode,
    double longitudeOfPeriapsis,
    double meanLongitude,
    time_t epoch)
{
}

void orbital_elements_calculate_extended_parameters(orbital_elements elem)
{
}
//
// Testing
//

#ifdef TEST

void vector_test()
{
    printf("Starting Vector Tests:\n");
    printf("--------------------\n");

    //init
    vector a = vector_init(10, 5, 7);
    vector b = vector_init(4, 3, 2);
    vector tmp = vector_init(-11, 5, 7);

    // vector_print(a);
    // vector_print(b);
    // vector_print(tmp);

    printf("vector_distance_to: ");
    assert(round(vector_distance_to(a, b)) == round(8.062258));
    printf("PASSED\n");

    printf("vector_abs: ");
    vector_abs(tmp);
    assert(tmp->x == 11);
    printf("PASSED\n");

    printf("vector_distance_to: ");
    assert(round(vector_length(a)) == round(13.190905));
    printf("PASSED\n");

    printf("vector_normalize: ");
    vector_normalize(tmp);
    assert(round(tmp->x) == round(0.787726));
    printf("PASSED\n");

    printf("vector_dot: ");
    assert(vector_dot(a, b) == 69);
    printf("PASSED\n");

    printf("vector_cross: ");
    vector_free(tmp);
    tmp = vector_cross(a, b);
    assert(tmp->y == 8);
    printf("PASSED\n");

    printf("vector_add: ");
    vector_free(tmp);
    tmp = vector_add(a, b);
    assert(tmp->x == a->x + b->x);
    printf("PASSED\n");

    printf("vector_sub: ");
    vector_free(tmp);
    tmp = vector_sub(a, b);
    assert(tmp->x == a->x - b->x);
    printf("PASSED\n");

    printf("vector_div: ");
    vector_free(tmp);
    tmp = vector_div(a, b);
    assert(tmp->x == a->x / b->x);
    printf("PASSED\n");

    printf("vector_mul: ");
    vector_free(tmp);
    tmp = vector_mul(a, b);
    assert(tmp->x == a->x * b->x);
    printf("PASSED\n");

    //free
    vector_free(a);
    vector_free(b);
    vector_free(tmp);
}

void orbital_elements_test()
{
    orbital_elements elm;
    elm = orbital_elements_create();

    orbital_elements_free(elm);
}

int main()
{
    vector_test();

    orbital_elements_test();

    return 0;
}

#endif