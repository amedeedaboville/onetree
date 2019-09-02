# Onetree: a fun little TSP solver.

More is covered in this post (future link here).

It was a pretty much line for line copy of the answer in (this post)[https://stackoverflow.com/questions/7159259/optimized-tsp-algorithms].

This finds exact solutions to smallish (~150 cities) TSP problems.

It uses branch and bound and the Held-Karp approximation.

It is entirely single threaded, though amenable to Rayon in a few places.

It was a toy project of mine; do whatever you'd like with it :)