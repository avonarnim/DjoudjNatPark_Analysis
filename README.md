# Djoudj National Park Waterbirds

## Log

Thursday, April 27

- We created a model for the case where a disease breaks out amongst a migrating predator population.

Tuesday, April 25

- We slightly modified terms in `djoudjTouristImpactOnFish.m` and `djoudjTouristImpactOnFishSimplified.m` so that we could account for how much fishing depends on the tourist population.

Sunday, April 23

- We implemented the initial paper's model, finding slightly different results than the original paper. We used Euler's method to approximate the progress of the system.
- One issue with the paper is that it does not identify what the value of `D` (the half saturation constant) should be. This value is supposed to indicate/control the prey population at which the growth rate will be half of the maximum rate.
