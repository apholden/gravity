#include <cmath>
#include <iostream>
#include <vector>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

#include <boost/multiprecision/cpp_dec_float.hpp>

typedef boost::multiprecision::cpp_dec_float_50 Real;

const Real G = 6.67430e-11;  // m^3 kg^-1 s^-2
const int total_steps = 1000000;
const int steps_per_print = total_steps / 100;

/**
 * An entity in a 2D space.
 * @note The entity has a mass, position, and velocity.
 * @note The entity is affected by gravity.
 */
class Entity
{
public:
    Entity(Real mass, Real x, Real y)
        : mass(mass), x(x), y(y) {}

    Real const mass;
    Real x;
    Real y;
    Real v_x = 0;
    Real v_y = 0;
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
 * Compute gravitational accelerations for all entities at the given positions.
 * @param entities The entities (used for masses).
 * @param pos_x The x-coordinates of each entity.
 * @param pos_y The y-coordinates of each entity.
 * @return A vector of (a_x, a_y) pairs for each entity.
 */
std::vector<std::pair<Real,Real>> compute_accelerations(
    std::vector<Entity> const &entities,
    std::vector<Real> const &pos_x,
    std::vector<Real> const &pos_y)
{
    size_t n = entities.size();
    std::vector<std::pair<Real,Real>> acc(n, {Real(0), Real(0)});

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            auto dx = pos_x[i] - pos_x[j];
            auto dy = pos_y[i] - pos_y[j];
            auto d_squared = dx * dx + dy * dy;
            auto f = G * entities[i].mass * entities[j].mass / d_squared;

            auto theta = atan2(dy, dx);
            auto f_x = f * cos(theta);
            auto f_y = f * sin(theta);

            acc[i].first -= f_x / entities[i].mass;
            acc[i].second -= f_y / entities[i].mass;
            acc[j].first += f_x / entities[j].mass;
            acc[j].second += f_y / entities[j].mass;
        }
    }

    return acc;
}

/**
 * Advance all entities by one time step using the 4th-order Runge-Kutta method.
 * @param entities The entities to update.
 * @param dt The time step.
 * @note For each entity the state is (x, y, v_x, v_y) and the derivative is
 *       (v_x, v_y, a_x, a_y) where accelerations come from gravitational interactions.
 */
void rk4_step(std::vector<Entity> &entities, Real dt)
{
    size_t n = entities.size();

    // Save initial state
    std::vector<Real> x0(n), y0(n), vx0(n), vy0(n);
    for (size_t i = 0; i < n; ++i)
    {
        x0[i] = entities[i].x;
        y0[i] = entities[i].y;
        vx0[i] = entities[i].v_x;
        vy0[i] = entities[i].v_y;
    }

    // k1: derivatives at the current state
    auto a1 = compute_accelerations(entities, x0, y0);

    // State at t + dt/2 using k1
    std::vector<Real> x1(n), y1(n), vx1(n), vy1(n);
    for (size_t i = 0; i < n; ++i)
    {
        x1[i] = x0[i] + dt/2 * vx0[i];
        y1[i] = y0[i] + dt/2 * vy0[i];
        vx1[i] = vx0[i] + dt/2 * a1[i].first;
        vy1[i] = vy0[i] + dt/2 * a1[i].second;
    }

    // k2: derivatives at the midpoint using k1
    auto a2 = compute_accelerations(entities, x1, y1);

    // State at t + dt/2 using k2
    std::vector<Real> x2(n), y2(n), vx2(n), vy2(n);
    for (size_t i = 0; i < n; ++i)
    {
        x2[i] = x0[i] + dt/2 * vx1[i];
        y2[i] = y0[i] + dt/2 * vy1[i];
        vx2[i] = vx0[i] + dt/2 * a2[i].first;
        vy2[i] = vy0[i] + dt/2 * a2[i].second;
    }

    // k3: derivatives at the midpoint using k2
    auto a3 = compute_accelerations(entities, x2, y2);

    // State at t + dt using k3
    std::vector<Real> x3(n), y3(n), vx3(n), vy3(n);
    for (size_t i = 0; i < n; ++i)
    {
        x3[i] = x0[i] + dt * vx2[i];
        y3[i] = y0[i] + dt * vy2[i];
        vx3[i] = vx0[i] + dt * a3[i].first;
        vy3[i] = vy0[i] + dt * a3[i].second;
    }

    // k4: derivatives at the endpoint using k3
    auto a4 = compute_accelerations(entities, x3, y3);

    // Final RK4 combination: y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
    for (size_t i = 0; i < n; ++i)
    {
        entities[i].x = x0[i] + dt/6 * (vx0[i] + 2*vx1[i] + 2*vx2[i] + vx3[i]);
        entities[i].y = y0[i] + dt/6 * (vy0[i] + 2*vy1[i] + 2*vy2[i] + vy3[i]);
        entities[i].v_x = vx0[i] + dt/6 * (a1[i].first + 2*a2[i].first + 2*a3[i].first + a4[i].first);
        entities[i].v_y = vy0[i] + dt/6 * (a1[i].second + 2*a2[i].second + 2*a3[i].second + a4[i].second);
    }
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
        rk4_step(entities, Real("0.01"));

        if (i % steps_per_print == 0)
        {
            std::cout << entities << std::endl;
        }
    }

    std::cout << entities << std::endl;
    return 0;
}
