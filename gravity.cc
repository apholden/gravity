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

class Entity
{
public:
    Entity(Real mass, Real x, Real y)
        : mass(mass), x(x), y(y) {}

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

std::ostream &operator<<(std::ostream &os, Entity const &e)
{
    auto speed = sqrt(e.v_x * e.v_x + e.v_y * e.v_y);
    auto angle = atan2(e.v_y, e.v_x) * 180 / pi<Real>();
    os << 'p' << e.x << ',' << e.y
        << ", v" << speed << "∠" << angle;
    return os;
}

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

Real Get_attraction(Entity const &a, Entity const &b)
{
    auto dx = b.x - a.x;
    auto dy = b.y - a.y;
    auto d_squared = dx * dx + dy * dy;
    auto f = G * a.mass * b.mass / d_squared;
    return f;
}

void Add_force(Entity &inner, Entity &outer, Real force)
{
    auto dx = inner.x - outer.x;
    auto dy = inner.y - outer.y;
    
    auto theta = atan2(dy, dx);
    auto f_x = force * cos(theta);
    auto f_y = force * sin(theta);

    outer.a_x += f_x / outer.mass;
    outer.a_y += f_y / outer.mass;
    inner.a_x -= f_x / inner.mass;
    inner.a_y -= f_y / inner.mass;
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


