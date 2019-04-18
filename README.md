# indie-genre

The [`chippr`](https://github.com/aimalz/chippr) test suite is built on a forward model that deserves its own repo, the N-dimensional Generative Model (indie genre for short).

The core of the code uses [`pomegranate`](https://github.com/jmschrei/pomegranate) for one-dimensional PDFs that can be composed arbitrarily into a generic N-dimensional probability space from which points can be drawn and marginal PDFs can be evaluated.
One-dimensional marginal PDFs can be re-parameterized by [`qp`](https://github.com/aimalz/qp).
