#import "../lib.typ": *

#show: university-assignment.with(
  title: "Dynamics Homework",
  subtitle: "P01. ODE",
  author: "Suraj Dayma",
)

= Problem
Set up (define system, draw FBD, write ODEs) a particle problem. Just one particle. 2D or 3D, your choice. Use a force, or forces that you like (gravity, spring, air friction). Any example of interest. Find a numerical solution. Graph it. Animate it. Try to make an interesting observation.

= ODE

Random ODE solved via custom Euler solver.

```matlab
F1 = sin(t) * i_hat + cos(t) * j_hat;
F2 = cos(t) * i_hat + sin(t) * j_hat;

xdot = v;
vdot = F1 + F2;
```

= Trajectory

For the initial conditions of `x,y,vx,vy = z0 = [0, 0, 0.5, -0.5];`

#figure(
  image("imgs/States.png"),
  caption: [Variables over time],
)

#figure(
  image("imgs/PhasePlot.png"),
  caption: [Phase Plot],
)

#figure(
  image("imgs/animation.gif"),
  caption: [Animation of Trajectory],
)



