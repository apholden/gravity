#include <cmath>
#include <iostream>
#include <vector>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <boost/multiprecision/cpp_dec_float.hpp>
using boost::multiprecision::cpp_dec_float_50;

typedef cpp_dec_float_50 Real;

const Real G = 6.67430e-11;  // m^3 kg^-1 s^-2
const int total_steps = 1000000;
const int steps_per_print = total_steps / 100;

/**
 * An entity in a 2D space.
 * @note The entity has a mass, position, velocity, and acceleration.
 * @note The entity is affected by gravity.
 */
class Entity
{
public:
    Entity(Real mass, Real x, Real y)
        : mass(mass), x(x), y(y) {}

    /**
     * Update the entity's position and velocity.
     * @param dt The time step.
     * @note The entity's position is updated using the formula:
     *     x = x + v_x * dt
     *     y = y + v_y * dt
     * @note The entity's velocity is updated using the formula:
     *    v_x = v_x + a_x * dt
     *    v_y = v_y + a_y * dt
     * @note The entity's acceleration is reset to zero after updating the velocity.
     */
    void Update(Real dt)
    {
        v_x += a_x * dt;
        v_y += a_y * dt;
        x += v_x * dt;
        y += v_y * dt;

        a_x = 0;
        a_y = 0;
    }

    Real const mass;
    Real x;
    Real y;
    Real v_x = 0;
    Real v_y = 0;
    Real a_x = 0;
    Real a_y = 0;
};

/**
 * Output an entity.
 * @param os The output stream.
 * @param e The entity.
 * @return The output stream.
 * @note The entity is output in the format:
 *    p<x>,<y>, v<speed>∠<angle>
 * @note The angle is in degrees.
 * @note The angle is measured counter-clockwise from the x-axis.
 */
std::ostream &operator<<(std::ostream &os, Entity const &e)
{
    auto speed = sqrt(e.v_x * e.v_x + e.v_y * e.v_y);
    auto angle = atan2(e.v_y, e.v_x) * 180 / pi<Real>();
    os << 'p' << e.x << ',' << e.y
        << ", v" << speed << "∠" << angle;
    return os;
}

/**
 * Output a vector of entities.
 * @param os The output stream.
 * @param entities The vector of entities.
 * @return The output stream.
 * @note The entities are output in the format:
 *     p<x>,<y>, v<speed>∠<angle>  p<x>,<y>, v<speed>∠<angle>  ...
 * @note The entities are separated by two spaces.
 */
std::ostream &operator<<(std::ostream &os, std::vector<Entity> const &entities)
{
    auto space = "";
    for(auto const &e : entities)
    {
        os << space << e;
        space = "  ";
    }
    return os;
}

/**
 * Calculate the force of attraction between two entities due to gravity.
 * @param a The first entity.
 * @param b The second entity.
 * @return The force of attraction between the two entities.
 * @note The force is always attractive.
 * @note The force is calculated using Newton's law of universal gravitation.
 * @note The force is calculated using the formula:
 *      F = G * m1 * m2 / r^2
 */
Real Get_attraction(Entity const &a, Entity const &b)
{
    auto dx = b.x - a.x;
    auto dy = b.y - a.y;
    auto d_squared = dx * dx + dy * dy;
    auto f = G * a.mass * b.mass / d_squared;
    return f;
}

/**
 * Add the force of attraction between two entities due to gravity.
 * @param entity_1 The first entity.
 * @param entity_2 The second entity.
 * @param force The force of attraction between the two entities.
 * @note The force is always attractive.
 * @note The force is added to the acceleration of each entity.
 * @note The force is added using Newton's second law of motion.
 * @note The force is added using the formula:
 *      F = m * a
 *      a = F / m
 */
void Add_force(Entity &entity_1, Entity &entity_2, Real force)
{
    auto dx = entity_1.x - entity_2.x;
    auto dy = entity_1.y - entity_2.y;

    auto theta = atan2(dy, dx);
    auto f_x = force * cos(theta);
    auto f_y = force * sin(theta);

    entity_1.a_x -= f_x / entity_1.mass;
    entity_1.a_y -= f_y / entity_1.mass;
    entity_2.a_x += f_x / entity_2.mass;
    entity_2.a_y += f_y / entity_2.mass;
}

int main(int argc, char *argv[])
{
    std::vector<Entity> entities;
    entities.push_back(Entity(1, -1, 0));
    entities.push_back(Entity(1, 1, 0));
    entities.push_back(Entity(2, 1, 1));

    std::cout << entities << std::endl;

    for (int i = 0; i < total_steps; ++i)
    {
        // Loop over all pairs of entities and calculate the force between them due to gravity.
        // Then add the force to each entity.
        for(auto outer = entities.begin(); outer != entities.end(); ++outer)
        {
            for(auto inner = entities.begin(); inner != outer; ++inner)
            {
                auto force = Get_attraction(*outer, *inner);
                Add_force(*inner, *outer, force);
            }
        }

        for(auto &e : entities)
        {
            e.Update(0.01);
        }

        if (i % steps_per_print == 0)
        {
            std::cout << entities << std::endl;
        }
    }

    std::cout << entities << std::endl;
    return 0;
}


