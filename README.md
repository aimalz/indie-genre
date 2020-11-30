# `zeppole` (backronym TBD)

The [`chippr`](https://github.com/aimalz/chippr) test suite is built on a forward model that deserves its own repo to facilitate more generic applications of forward-modeled photo-z posteriors.

The core of the code uses [`pomegranate`](https://github.com/jmschrei/pomegranate) for one-dimensional PDFs that can be composed generically into any arbitrary 2-dimensional probability space of true and estimated scalars from which points can be drawn and marginal PDFs can be evaluated.
One-dimensional marginal PDFs can be re-parameterized by [`qp`](https://github.com/aimalz/qp).

At some point I'd like to adapt this to higher-dimensional quantities, so 2N dimensional probability spaces, but the scope is currently limited to photometric redshift PDFs needing a quick and dirty solution.
