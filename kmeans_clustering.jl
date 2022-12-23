using Clustering, Plots
using Random

# initialize plot
gr(size = (400, 400))
# generate random points
Random.seed!(1)

f1 = rand(100)
f2 = rand(100)

p_rand = scatter(f1, f2)

X = [f1 f2]'

k = 2
itr = 100

Random.seed!(1)

result = kmeans(X, k; maxiter = itr, display = :iter)

a = assignments(result)
c = counts(a)
mu = result.centers

p_kmeans_demo = scatter(f1, f2, group = a)

scatter!(mu[1, :], mu[2, :])

## Application ##
using RDatasets

cats = dataset("boot", "catsM")

vscodedisplay(cats)

p_cats = scatter(cats.BWt, cats.HWt)

f1 = cats.BWt
f2 = cats.HWt

f1_min = minimum(f1)
f2_min = minimum(f2)

f1_max = maximum(f1)
f2_max = maximum(f2)

f1_n = (f1 .- f1_min) ./ (f1_max - f1_min)
f2_n = (f2 .- f2_min) ./ (f2_max - f2_min)

X = [f1_n f2_n]'

p_cats = scatter(f1_n, f2_n)

k = 3

itr = 200000

Random.seed!(1)

result = kmeans(X, k; maxiter = itr, display = :iter)

a = assignments(result)
c = counts(a)
mu = result.centers

p_cats_n = scatter(f1_n, f2_n, group = a)

scatter!(mu[1, :], mu[2, :], color = :yellow, markersize = 10)
 
