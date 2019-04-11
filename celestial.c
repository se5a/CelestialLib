/*
    -Myvar:
    Coding conventions:
    1. Please add the thing your using in an include, with a comment
        ->
            //malloc, free
            #include <stdlib.h>

    2. Use Comment with one empty one above an one empty below
        -> //
           // Some Comment
           //

    3. Use Sections to make code navagation easyer
        -> //
           // [Helpers]
           //
           inline double rad_to_deg(double rad) { ... }
           
    4. use @ tags to mark things: (ensure tag index is upto date)
        -> // @Hack !!! this code is hacky
           // @Cleanup Code is messy the names need to be updated
           // @Speed !! malloc a big block insted of meany small malocs for cacheing
 */

//
// Tag Index
//

/*
 * @Speed
 * @Cleanup
 * @Hack
 */

//malloc, free
#include <stdlib.h>

//printf
#include <stdio.h>

//time_t
#include <time.h>

//M_PI
#include <math.h>

#include "celestial.h"

/*
    TODO:
    * Orbital Mechanics
    * Stellar classification
*/

//
// [Implmentation]
//

//
// [Helpers]
//
inline double rad_to_deg(double rad)
{
    return rad * (180.0 / M_PI);
}

inline double deg_to_rad(double deg)
{
    return deg * (M_PI / 180.0);
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
};

// @Speed ! mallocing a big struct like this piecemeal is slow, it would be way faster to use an Array for Struct, with block mallocs
orbital_elements create_orbital_elements()
{
    return malloc(sizeof(struct orbital_elements_t));
}

void free_orbital_elements(orbital_elements h) { free(h); }

//
// Testing
//

#ifdef TEST

int main()
{
    orbital_elements elm;
    elm = create_orbital_elements();

    free_orbital_elements(elm);

    printf("Hello, World! \n");
    return 0;
}

#endif