# gravity
A simple physics simulation to observe the effect of Newtonian gravity on a multi-body simulation

In its current form, it places 3 stationary masses in a 2 x 1 meter space and
simulates Newtonian gravity to accelerate them toward each other.  The simulation
takes 10-millisecond steps and runs for 1,000,000 steps.  It prints the position
and velocity of the objects every 10,000 steps.

Using this, I was able to observe the velocities change direction over time,
implying that 3 bodies, placed in an asymetric pattern, is enough to introduce
angular momentum.

To compile this code, it depends on the [Boost Libraries](https://www.boost.org/),
specifically the math constants and multiprecision libraries.