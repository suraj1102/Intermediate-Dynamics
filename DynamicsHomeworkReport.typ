#import "lib.typ": *

#show: university-assignment.with(
  title: "Dynamics Homework",
  subtitle: "Dynamics & Simulation",
  author: "Suraj Dayma",
)

// ---------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------

// Problem heading, e.g. "P13. Statics vs Dynamics".
#let problem(num, title) = heading(level: 1)[P#num. #title]

// Short statement of the question.
#let question(body) = block(above: 0.4em, below: 1em, text(style: "italic", body))

// Original report content: headings nest one level under the problem heading.
#let report(body) = {
  set heading(offset: 1)
  body
}

// Two figures side by side.
#let pair(a, b) = grid(
  columns: (1fr, 1fr),
  column-gutter: 1em,
  align: horizon,
  a, b,
)

// Embed PDF pages two per row to save space.
#let pdf-pages(path, pages) = grid(
  columns: (1fr, 1fr),
  column-gutter: 0.8em,
  row-gutter: 0.8em,
  ..range(1, pages + 1).map(p => image(path, page: p, width: 100%)),
)

// lib.typ's page numbering takes (current, total); the outline only passes
// one value, so outline entries print the plain page number instead.
#show outline.entry: it => {
  show link: l => l.body
  link(
    it.element.location(),
    it.indented(it.prefix(), [
      #it.body()
      #box(width: 1fr, repeat[.#h(4pt)])
      #str(it.element.location().page())
    ]),
  )
}

#set heading(numbering: none)

#outline(title: [Contents], depth: 2, indent: auto)
#pagebreak()


#problem("01", "ODE — A Particle Problem")
#report[
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
    image("01_ODE/imgs/States.png"),
    caption: [Variables over time],
  )

  #pair(
    figure(image("01_ODE/imgs/PhasePlot.png", width: 100%), caption: [Phase Plot]),
    figure(image("01_ODE/imgs/animation.gif", width: 100%), caption: [Animation of Trajectory]),
  )
]

#problem("02", "Vectors in 3D")
#question[
  Given four random points $bup(r)_A, bup(r)_B, bup(r)_C, bup(r)_D$ and a force
  $bup(F)$, write single vector formulas for: (a) the minimum distance between
  lines $A B$ and $C D$; (b) the end points of that shortest segment; (c) the
  volume of tetrahedron $A B C D$ (derive "⅓ base × height"); (d) the tension in
  bar $A D$ when $A, B, C$ are fixed and connected by rods to $D$, where
  $bup(F)$ is applied.
]

#problem("03", "More ODE & Animation Practice")
#question[
  Pick a simple ODE system (harmonic oscillator, pendulum, Lorenz …) and solve it
  three ways: (a) a single script with a plain Euler loop; (b) your own Euler
  solver function calling a RHS function; (c) ODE45. Then (d) study accuracy vs.
  step size with Euler and (e) vary `RelTol`/`AbsTol` in ODE45 to judge whether
  MATLAB estimates its own accuracy well.
]

#problem("06", "Spring and Mass (2D)")
#question[
  A mass $m$ on a spring $(k, L_0)$ pinned at the origin, under gravity, starting
  from $bup(r)_0, bup(v)_0$. (a) Find the EoM; (b) plot and animate the
  trajectory; (c) list and carry out as many sanity checks as possible
  (e.g. $k = 0$ → parabola; $x_0 = v_(x 0) = 0$ → stays on the $y$-axis;
  $g = 0, bup(v)_0 = 0$ → radial harmonic oscillation; $L_0 = 0$ → …).
]

#problem("07", "Simple Animation of a Shape")
#question[
  Draw an object (face, house, …) and make it move smoothly and interestingly on
  screen using only rigid motions — translations and rotations, no distortion.
]

#problem("08", "Simplest Dynamics in Polar Coordinates")
#question[
  A force-free particle on a plane, posed in polar coordinates: (a) write
  $bup(a) = bup(0)$ as ODEs in $(r, theta)$; (b) solve numerically for various
  ICs; (c) check the path is a straight line at constant speed; (d) define a
  quantitative "straightness" measure and evaluate it; (e) does it improve as
  tolerances are tightened?
]

#problem("09", "Cannon Ball (Linear Drag)")
#question[
  A cannon ball launched at angle $theta$ and speed $v_0$ with gravity and linear
  viscous drag $|c v|$. (a) Analytical $x(t), y(t)$; (b) numerical solution for
  $theta = pi slash 4$, $v_0 = g = m = c = 1$; (c) compare — error at $t = 2$ vs.
  tolerances; (d) trajectories as $v_0 -> infinity$ — what shape?; (e) where does
  the range-maximising angle $theta^*$ tend as $v_0 -> infinity$?
]

#problem("11", "Central Force")
#question[
  Only force is central, $bup(F) = -F(r) hat(bup(e))_r$. (a) Find $F(r)$ such
  that, across circular orbits, the speed is independent of radius.
  (b) Numerically (with `fsolve` over ICs and period) find a periodic orbit that
  is neither circular nor straight, for a force law other than a linear spring or
  inverse-square gravity; comment on chaotic-looking failed searches and on the
  minimum number of loops for power laws.
]
#report[
  #figure(image("11_CentralForce/imgs/P11_PartA.jpeg", width: 40%))

  = Periodic Solutions
  == For $F(r) = 1$
  Found using fsolve root finding
  with initial conditions - \[ 4.837094489229557 0.113739688509354 0.481812323218859 1.333913034846417 \]
  time period = 49.578757824567518

  Using jacobian analysis, we find this is an unstable periodic motion.
  eigvals of jacobian = \[ 0.786327175145648 0.786327175145648 1.187326408237210 0.842228377131789\]

  #figure(
    image("11_CentralForce/imgs/f1_1_unstable.png", width: 80%),
    caption: [Trajectories],
  )
]

#problem("13", "Statics vs Dynamics — Two Spring Damper System")
#question[
  A 2D point mass under gravity, connected to anchors $A$ and $B$ by springs
  $(k_A, L_A), (k_B, L_B)$ with parallel dashpots $(c_A, c_B)$ and air drag $c_d$.
  (a) Find an equilibrium; (b) simulate the motion; (c) find _all_ equilibria and
  the most/least number possible as parameters vary; (d) determine stability both
  by simulation near the equilibrium and by eigenvalues of the linearised system.
]
#report[
  = Intro

  System animation, fixed point analysis and clustering, and linear stability analysis for the *Two Spring and Two Damper System* with drag.

  = Problem and Dynamics:

  #figure(image("13_SpringAndMass/imgs/equations.png", width: 70%), caption: [Force Balance])

  = Fixed Points

  The fixed points were found using `fsolve`, and clustering of fixed points was performed by using the following algorithms:
  1. KMeans: within cluster sum of squares was found was varying K and the best K was determined by the elbow method.
  2. DBScan: the distance parameter was kept very small (0.005) which converged correctly for every scenario I tried (except the infinite fixed point cases).

  == Fixed Points = 2


  #grid(
    columns: (35%, 1fr),
    column-gutter: 1em,
    figure(
      ```matlab
      % System Constants
      p.m = 1;
      p.g = 10;

      p.b1 = 1;
      p.b2 = 1;

      p.d = 0.6; % drag coef

      p.k1 = 1;
      p.k2 = 10;
      p.l1_0 = 1;
      p.l2_0 = 1;

      p.r1 = [-1; 0];
      p.r2 = [1; 0];
      ```,
    ),
    figure(
      image("13_SpringAndMass/imgs/2FP.png", width: 100%),
      caption: [Fixed Points],
    ),
  )

  Both KMeans and DBScan provided the right answer for this.

  == Fixed Points = 5

  #grid(
    columns: (35%, 1fr),
    column-gutter: 1em,
    figure(
      ```matlab
      % System Constants
      p.m = 1;
      p.g = 0;

      p.b1 = 1;
      p.b2 = 1;

      p.d = 1; % drag coef

      p.k1 = 10;
      p.k2 = 10;
      p.l1_0 = 1.5;
      p.l2_0 = 1.5;

      p.r1 = [-1; 0];
      p.r2 = [1; 0];
      ```,
    ),
    figure(
      image("13_SpringAndMass/imgs/clustering/5FP.png", width: 100%),
      caption: [Fixed Points],
    ),
  )


  #pair(
    figure(image("13_SpringAndMass/imgs/clustering/5FP-KMeans.png", width: 100%), caption: [KMean - Best K = 3]),
    figure(
      image("13_SpringAndMass/imgs/clustering/5FP-DB.png", width: 100%),
      caption: [DBScan - Correct Identification],
    ),
  )

  == Fixed Points = 1

  #grid(
    columns: (35%, 1fr),
    column-gutter: 1em,
    figure(``` matlab
    % System Constants
    p.m = 1;
    p.g = 10;

    p.b1 = 1;
    p.b2 = 1;;

    p.d = 1; % drag coef

    p.k1 = 10;
    p.k2 = 10;
    p.l1_0 = 1;
    p.l2_0 = 1;

    p.r1 = [-1; 0];
    p.r2 = [1; 0];
    ```),
    figure(
      image("13_SpringAndMass/imgs/1FP.png", width: 100%),
      caption: [Fixed Points],
    ),
  )

  Both KMeans and DBScan provided the right answer for this.

  == Fixed Points = 4, Two of them are very close

  #grid(
    columns: (35%, 1fr),
    column-gutter: 1em,
    figure(``` matlab
    % System Constants
    p.m = 1;
    p.g = 0;

    p.b1 = 1;
    p.b2 = 1;;

    p.d = 1; % drag coef

    p.k1 = 10;
    p.k2 = 10;
    p.l1_0 = 1;
    p.l2_0 = 1;

    p.r1 = [-1; 0];
    p.r2 = [1; 0];
    ```),
    figure(
      image("13_SpringAndMass/imgs/clustering/closeFP.png", width: 100%),
      caption: [Fixed Points],
    ),
  )

  #figure(image("13_SpringAndMass/imgs/clustering/closeFP-zoomed.png", width: 50%), caption: [2 Fixed Points Closeby])


  #pair(
    figure(
      image("13_SpringAndMass/imgs/clustering/closeFP-KMeans.png", width: 100%),
      caption: [KMean - Best K = 3: WRONG],
    ),
    figure(
      image("13_SpringAndMass/imgs/clustering/closeFP-DB.png", width: 100%),
      caption: [DBScan - Correct Identification],
    ),
  )

  == Fixed Points = 3


  #grid(
    columns: (35%, 1fr),
    column-gutter: 1em,
    figure(
      ``` matlab
      % System Constants
      p.m = 1;
      p.g = 1;

      p.b1 = 1;
      p.b2 = 1;

      p.d = 1; % drag coef

      p.k1 = 1;
      p.k2 = 1;
      p.l1_0 = 1;
      p.l2_0 = 2;

      p.r1 = [-1; 0];
      p.r2 = [1; 0];
      ```,
    ),
    figure(
      image("13_SpringAndMass/imgs/3FP.png", width: 100%),
      caption: [Fixed Points],
    ),
  )

  DBScan was correct all the times, KMeans sometimes said best K was 2, while other times said it was 3.


  == Fixed Points = Infinite

  When gravity is off, and either of the spring constants are zero, all points on the free length circle of the active spring are valid unique fixed points.

  #figure(
    image("13_SpringAndMass/imgs/infFP.png", width: 50%),
    caption: [g = 0, k1 = 0 => inf FP],
  )

  Also when gravity is off, and both the springs are active, but the anchor positions are the same, we get infinite fixed points.

  #figure(
    image("13_SpringAndMass/imgs/infFP2.png", width: 50%),
    caption: [g = 0, rA=rB => inf FP],
  )


  = Linearised Motion

  *Task:*
  - Using the Jacobian supplied by fsolve at one of the stable fixed points, solve the linearized equations for motion for an initial condition near the fixed point.
  - Using the same initial condition, use ODE45 to find the motion.
  - Repeat for an unstable fixed point

  *System Used:*
  The same parameters that resulted in 3 fixed points were used.

  Stable Fixed Points:
  - [-1; 0]
  - [-0.5323; -1.7316]

  Unstable Fixed Points:
  - [1.2412; 0.2169]

  == Response near stable fixed point
  #figure(
    image("13_SpringAndMass/imgs/linResponseStable.png", width: 50%),
  )

  == Response near unstable fixed point
  #figure(
    image("13_SpringAndMass/imgs/linResponseUnstable.png", width: 50%),
  )

  == Response When L0 = 0, b = 0
  #pair(
    figure(image("13_SpringAndMass/imgs/linResponseAllZero.png", width: 100%)),
    figure(image("13_SpringAndMass/imgs/linResponseAllZero2.png", width: 100%)),
  )
]

#pagebreak(weak: true)
#problem("14", "What Means \"Rate of Change of Angular Momentum\"?")
#question[
  For a moving particle $P$ and a moving reference point $C$, decide for which
  definition of $bup(H)_(slash C)$ the equation $bup(M)_C = dot(bup(H))_(slash C)$
  holds (i.e. is consistent with $bup(F) = m bup(a)$): (a)
  $bup(r)_(P slash C') times bup(v)_(P slash 0) m$ with $C'$ fixed and
  instantaneously coincident with $C$; (b) $bup(r)_(P slash C) times bup(v)_(P slash 0) m$;
  (c) $bup(r)_(P slash C) times bup(v)_(P slash C) m$ — in general, or only in
  special cases?
]

#pdf-pages("14_AngularMomentum/P14_AngularMomentum.pdf", 2)

#problem("15", "Potential Energy")
#question[
  Find the potential-energy function for a set of 3D force fields: constant
  force, zero-rest-length spring, spring with rest length, inverse-square central
  force, general central force $f(r)$, spring and inverse-square attraction
  between two points, constant force along $hat(bup(lambda))$, and an $x$-only
  force $f(x)$ — verifying the given answers with multivariable calculus.
]

#pdf-pages("15_PotentialFunction/P15_PotentialFunction.pdf", 3)
