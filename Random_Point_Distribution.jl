# cd(@__DIR__)

using Random       
using Statistics   
using CairoMakie  

cd(@__DIR__)

# Squared distance between two points. Note, this is computationally easier than finding distance of two points, but 
# is harmless for searching; squaring is a monotonic function. I report Euclidian distance by taking square root later.
function dist2(p::Tuple{Float64,Float64}, q::Tuple{Float64,Float64})
    dx = p[1] - q[1]
    dy = p[2] - q[2]
    return dx * dx + dy * dy
end

# ---------------------------
# Divide & Conquer / Convex Hull helpers
# ---------------------------

# Recursive function for the divide and conquer closest pair.
# Px and Py are the points sorted by x and y, respectively.
function closest_pair_dc(Px::Vector{Tuple{Float64,Float64}}, Py::Vector{Tuple{Float64,Float64}})
    n = length(Px)
    if n ≤ 3
        best = Inf
        best_pair = (Px[1], Px[2])
        for i in 1:(n-1)
            for j in i+1:n
                d = dist2(Px[i], Px[j])
                if d < best
                    best = d
                    best_pair = (Px[i], Px[j])
                end
            end
        end
        return best, best_pair
    end

    mid = n ÷ 2
    midPoint = Px[mid]
    Qx = Px[1:mid]
    Rx = Px[mid+1:end]

    Qy = Tuple{Float64,Float64}[]
    Ry = Tuple{Float64,Float64}[]
    for p in Py
        if p[1] ≤ midPoint[1]
            push!(Qy, p)
        else
            push!(Ry, p)
        end
    end

    d_left, pair_left = closest_pair_dc(Qx, Qy)
    d_right, pair_right = closest_pair_dc(Rx, Ry)
    if d_left < d_right
        d_min = d_left
        best_pair = pair_left
    else
        d_min = d_right
        best_pair = pair_right
    end

    delta = sqrt(d_min)
    strip = Tuple{Float64,Float64}[]
    for p in Py
        if abs(p[1] - midPoint[1]) < delta
            push!(strip, p)
        end
    end

    m = length(strip)
    for i in 1:(m-1)
        # Check at most the next 7 points (or until the y-distance is ≥ delta)
        for j in i+1:min(i+7, m)
            if (strip[j][2] - strip[i][2]) ≥ delta
                break
            end
            d = dist2(strip[i], strip[j])
            if d < d_min
                d_min = d
                best_pair = (strip[i], strip[j])
                delta = sqrt(d_min)
            end
        end
    end

    return d_min, best_pair
end

# Compute the closest pair using divide and conquer.
function closest_pair(points::Vector{Tuple{Float64,Float64}})
    # Sort points by x and by y.
    Px = sort(points, by = p -> p[1])
    Py = sort(points, by = p -> p[2])
    d2, pair = closest_pair_dc(Px, Py)
    return sqrt(d2), pair
end

# Cross product of OA and OB.
function cross(o::Tuple{Float64,Float64}, a::Tuple{Float64,Float64}, b::Tuple{Float64,Float64})
    return (a[1] - o[1]) * (b[2] - o[2]) - (a[2] - o[2]) * (b[1] - o[1])
end

# Compute the convex hull using Andrew’s monotone chain.
function convex_hull(points::Vector{Tuple{Float64,Float64}})
    sorted_points = sort(points, by = p -> (p[1], p[2]))
    if length(sorted_points) ≤ 1
        return sorted_points
    end

    lower = Tuple{Float64,Float64}[]
    for p in sorted_points
        while length(lower) ≥ 2 && cross(lower[end-1], lower[end], p) ≤ 0
            pop!(lower)
        end
        push!(lower, p)
    end

    upper = Tuple{Float64,Float64}[]
    for p in reverse(sorted_points)
        while length(upper) ≥ 2 && cross(upper[end-1], upper[end], p) ≤ 0
            pop!(upper)
        end
        push!(upper, p)
    end

    # Remove the last point of each list to avoid duplication.
    pop!(lower)
    pop!(upper)
    return vcat(lower, upper)
end

# Compute farthest pair of points on convex hull (brute force).
function farthest_pair(points::Vector{Tuple{Float64,Float64}})
    hull = convex_hull(points)
    m = length(hull)
    @assert m ≥ 2 "Convex hull has fewer than 2 points."
    max_d2 = -Inf
    best_pair = (hull[1], hull[2])
    for i in 1:(m-1)
        for j in i+1:m
            d = dist2(hull[i], hull[j])
            if d > max_d2
                max_d2 = d
                best_pair = (hull[i], hull[j])
            end
        end
    end
    return sqrt(max_d2), best_pair
end

###############################################################################
# Implementations for find_closest_and_farthest
###############################################################################

# Brute-force implementation (O(n²) for both closest and farthest pairs).
function find_closest_and_farthest_brute(points::Vector{Tuple{Float64,Float64}})
    n = length(points)
    @assert n ≥ 2 "Need at least 2 points."
    min_dist2 = Inf
    max_dist2 = -Inf
    min_i, min_j = 1, 2
    max_i, max_j = 1, 2

    for i in 1:(n-1)
        p1 = points[i]
        for j in i+1:n
            p2 = points[j]
            d2 = dist2(p1, p2)
            if d2 < min_dist2
                min_dist2 = d2
                min_i, min_j = i, j
            end
            if d2 > max_dist2
                max_dist2 = d2
                max_i, max_j = i, j
            end
        end
    end

    return sqrt(min_dist2), (points[min_i], points[min_j]),
           sqrt(max_dist2), (points[max_i], points[max_j])
end

# Divide & Conquer / Convex Hull implementation.
function find_closest_and_farthest_dc(points::Vector{Tuple{Float64,Float64}})
    @assert length(points) ≥ 2 "Need at least 2 points."
    # Closest pair using divide & conquer.
    min_dist, closest_pair_points = closest_pair(points)
    # Farthest pair using convex hull and brute force.
    max_dist, farthest_pair_points = farthest_pair(points)
    return min_dist, closest_pair_points, max_dist, farthest_pair_points
end

###############################################################################
# Trial functions (run_trials)
###############################################################################

# Runs trials using the brute force implementation.
function run_trials_brute(n::Int, r::Int)
    closest_points  = Vector{Tuple{Float64,Float64}}(undef, 2*r)
    farthest_points = Vector{Tuple{Float64,Float64}}(undef, 2*r)
    closest_dists   = Vector{Float64}(undef, r)
    farthest_dists  = Vector{Float64}(undef, r)

    for trial in 1:r
        points = [(rand(), rand()) for _ in 1:n]
        min_dist, (cp1, cp2), max_dist, (fp1, fp2) = find_closest_and_farthest_brute(points)
        closest_points[2*trial - 1] = cp1
        closest_points[2*trial    ] = cp2
        farthest_points[2*trial - 1] = fp1
        farthest_points[2*trial    ] = fp2
        closest_dists[trial] = min_dist
        farthest_dists[trial] = max_dist
    end

    return closest_points, farthest_points, closest_dists, farthest_dists
end

# Runs trials using the divide & conquer / convex hull implementation.
function run_trials_dc(n::Int, r::Int)
    closest_points  = Vector{Tuple{Float64,Float64}}(undef, 2*r)
    farthest_points = Vector{Tuple{Float64,Float64}}(undef, 2*r)
    closest_dists   = Vector{Float64}(undef, r)
    farthest_dists  = Vector{Float64}(undef, r)

    for trial in 1:r
        points = [(rand(), rand()) for _ in 1:n]
        min_dist, (cp1, cp2), max_dist, (fp1, fp2) = find_closest_and_farthest_dc(points)
        closest_points[2*trial - 1] = cp1
        closest_points[2*trial    ] = cp2
        farthest_points[2*trial - 1] = fp1
        farthest_points[2*trial    ] = fp2
        closest_dists[trial] = min_dist
        farthest_dists[trial] = max_dist
    end

    return closest_points, farthest_points, closest_dists, farthest_dists
end


###############################################################################
# Benchmarking
###############################################################################

"""
Runs the simulation for n = n_min, n_min+n_step, …, n_max.
For each n, it runs the simulation num_median_trials times,
measures time for each run, and then uses the median.
It performs this for both implementations (Brute Force and Divide & Conquer/Convex Hull).
Then I create corresponding plot.
"""

function run_benchmark(r::Int, num_median_trials::Int, n_min::Int, n_max::Int, n_step::Int, seed::Int)
    ns = collect(n_min:n_step:n_max)
    median_times_brute = Float64[]
    median_times_dc    = Float64[]
    
    for n in ns
        times_brute = Float64[]
        times_dc    = Float64[]
        for trial in 1:num_median_trials
            Random.seed!(seed)
            t_brute = @elapsed run_trials_brute(n, r)
            push!(times_brute, t_brute)
            
            Random.seed!(seed)
            t_dc = @elapsed run_trials_dc(n, r)
            push!(times_dc, t_dc)
        end
        push!(median_times_brute, median(times_brute))
        push!(median_times_dc, median(times_dc))
    end
    # Plot
    fig = Figure(size = (700, 500))
    ax = Axis(fig[1, 1],
    xlabel = "Number of Points (n)",
    ylabel = "Time (seconds)",
    title  = "Benchmark: Execution Time",
    xticks = 200:200:1200  
)

    # Plot Brute Force results (blue) and Divide & Conquer/Convex Hull results (red).
    lines!(ax, ns, median_times_brute, color = :blue, linestyle = :solid, label = "Brute Force")
    scatter!(ax, ns, median_times_brute, marker = :circle, color = :blue)
    lines!(ax, ns, median_times_dc, color = :red, linestyle = :solid, label = "Divide & Conquer / Convex Hull")
    scatter!(ax, ns, median_times_dc, marker = :circle, color = :red)

    # Legend
    axislegend(ax, position = (:left, :top))
   
     # Save
    filename = "benchmark_plot.png"
    CairoMakie.save(filename, fig)
    println("Saved benchmark plot as $filename")
    display(fig)
end

###############################################################################
# Main
###############################################################################

function main()
    seed = 42                        # Set seed 
    num_median_trials = 100         # Number of simulation runs over which the median is computed.
    r = 50                         # Number of internal trials per run.

    # Define grid parameters for n.
    n_min = 50
    n_max = 1200
    n_step = 50

    run_benchmark(r, num_median_trials, n_min, n_max, n_step, seed)
end;