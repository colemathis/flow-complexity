#==============================================================================#
# Residence Time Analysis: Lattice vs. Randomized Topology
#
# Addresses Reviewer 1 Major 4 / Reviewer 3 HL3:
#   Does topology matter because of structural properties, or simply because
#   randomization creates shortcuts to the outflow, reducing residence time?
#
# Two analyses:
#   1. Random walk first-passage time (FPT) from inflow (node 1) to outflow (node 25)
#   2. Calibration: can the FPT shift explain the observed α divergence?
#
# Comparison: 5×5 bidirectional lattice vs. 50 randomized realizations
#==============================================================================#

using Graphs
using StatsBase
using Random
using Statistics
using Printf
using CSV
using DataFrames

#==============================================================================#
# GRAPH CONSTRUCTION (from Ensemble.jl)
#==============================================================================#

function lattice_digraph_bidirectionnal(N)
    """
    Create a bidirectional directed lattice graph with N nodes.
    """
    L = Graphs.SimpleDiGraph(N)
    n = Int(sqrt(N))

    for i in 1:n
        for j in 0:(n-1)
            current_node = i + n * j
            right_neighbor = current_node + 1
            down_neighbor = current_node + n

            if i < n
                Graphs.add_edge!(L, (current_node, right_neighbor))
                Graphs.add_edge!(L, (right_neighbor, current_node))
            end
            if current_node + n <= N
                Graphs.add_edge!(L, (current_node, down_neighbor))
                Graphs.add_edge!(L, (down_neighbor, current_node))
            end
        end
    end

    return L
end

function randomize_graph!(g::Graphs.SimpleDiGraph, nswap::Int)
    """
    Randomize the directed graph while preserving degree distribution.
    Swaps pairs of bidirectional edges.
    """
    edges = collect(Graphs.edges(g))
    swaps = 0
    tries = 0

    while swaps < nswap && tries < 100 * nswap
        tries += 1

        e1, e2 = rand(edges, 2)
        u1, v1 = Graphs.src(e1), Graphs.dst(e1)
        u2, v2 = Graphs.src(e2), Graphs.dst(e2)

        # Skip if any overlap between the four nodes
        if length(Set([u1, v1, u2, v2])) < 4
            continue
        end

        # Check that mirrored edges exist
        if !(Graphs.has_edge(g, v1, u1) && Graphs.has_edge(g, v2, u2))
            continue
        end

        # Proposed new edges and their mirrors must not already exist
        if Graphs.has_edge(g, u1, v2) || Graphs.has_edge(g, u2, v1) ||
           Graphs.has_edge(g, v2, u1) || Graphs.has_edge(g, v1, u2) ||
           u1 == v2 || u2 == v1
            continue
        end

        # Perform the swap for both edges and their mirrors
        Graphs.rem_edge!(g, u1, v1)
        Graphs.rem_edge!(g, v1, u1)
        Graphs.rem_edge!(g, u2, v2)
        Graphs.rem_edge!(g, v2, u2)

        Graphs.add_edge!(g, u1, v2)
        Graphs.add_edge!(g, v2, u1)
        Graphs.add_edge!(g, u2, v1)
        Graphs.add_edge!(g, v1, u2)

        edges = collect(Graphs.edges(g))
        swaps += 1
    end

    return g
end

#==============================================================================#
# ANALYSIS FUNCTIONS
#==============================================================================#

function random_walk_fpt(g::Graphs.SimpleDiGraph, source::Int, target::Int)
    """
    Perform a random walk from source to target on graph g.
    At each step, uniformly pick a random out-neighbor.
    Returns the number of steps to reach target.
    """
    current = source
    steps = 0
    while current != target
        neighbors = Graphs.outneighbors(g, current)
        current = rand(neighbors)
        steps += 1
    end
    return steps
end

function run_fpt_trials(g::Graphs.SimpleDiGraph, source::Int, target::Int, n_trials::Int)
    """
    Run n_trials random walks and return vector of first-passage times.
    """
    fpts = Vector{Int}(undef, n_trials)
    for i in 1:n_trials
        fpts[i] = random_walk_fpt(g, source, target)
    end
    return fpts
end

function fpt_summary(fpts::Vector{Int})
    """
    Compute summary statistics for a vector of first-passage times.
    """
    return (
        mean   = mean(fpts),
        median = median(fpts),
        std    = std(fpts),
        min    = minimum(fpts),
        max    = maximum(fpts),
        q25    = quantile(fpts, 0.25),
        q75    = quantile(fpts, 0.75),
        q95    = quantile(fpts, 0.95),
    )
end

#==============================================================================#
# MAIN
#==============================================================================#

function main()
    t_start = time()
    Random.seed!(42)

    # Parameters
    N = 25              # 5×5 lattice
    source = 1          # top-left (inflow)
    target = N          # bottom-right (outflow)
    n_trials = 100_000  # random walks per graph
    n_random = 50       # number of randomized realizations

    println("=" ^ 70)
    println("  Residence Time Analysis: Lattice vs. Randomized Topology")
    println("  N = $N nodes, $n_trials random walks, $n_random randomizations")
    println("=" ^ 70)

    #=== Phase 1: Build lattice ===#
    println("\n[Phase 1] Building 5×5 bidirectional lattice...")
    lattice = lattice_digraph_bidirectionnal(N)
    println("  Nodes: $(nv(lattice)), Edges: $(ne(lattice))")

    #=== Phase 2: Lattice FPT ===#
    println("\n[Phase 2] Running $n_trials random walks on lattice...")
    fpts_lattice = run_fpt_trials(lattice, source, target, n_trials)
    stats_lattice = fpt_summary(fpts_lattice)
    L_min_lattice = length(Graphs.a_star(lattice, source, target))
    println("  Shortest path (lattice): $L_min_lattice")
    println("  Mean FPT: $(round(stats_lattice.mean, digits=1))")
    println("  Median FPT: $(round(stats_lattice.median, digits=1))")

    #=== Phase 3: Randomized graphs ===#
    println("\n[Phase 3] Running analysis on $n_random randomized graphs...")

    # Storage for randomized results
    rand_mean_fpts = Vector{Float64}(undef, n_random)
    rand_median_fpts = Vector{Float64}(undef, n_random)
    rand_std_fpts = Vector{Float64}(undef, n_random)
    rand_shortest_paths = Vector{Int}(undef, n_random)

    for r in 1:n_random
        # Randomize a copy of the lattice
        g = copy(lattice)
        randomize_graph!(g, 10 * ne(g))

        # Verify connectivity
        sp = Graphs.a_star(g, source, target)
        if isempty(sp)
            println("  WARNING: Realization $r is disconnected (source→target), skipping.")
            rand_mean_fpts[r] = NaN
            rand_median_fpts[r] = NaN
            rand_std_fpts[r] = NaN
            rand_shortest_paths[r] = -1
            continue
        end
        rand_shortest_paths[r] = length(sp)

        # FPT analysis
        fpts_r = run_fpt_trials(g, source, target, n_trials)
        s = fpt_summary(fpts_r)
        rand_mean_fpts[r] = s.mean
        rand_median_fpts[r] = s.median
        rand_std_fpts[r] = s.std

        if r % 10 == 0
            print("  Completed $r / $n_random realizations\r")
        end
    end
    println("  Completed $n_random / $n_random realizations")

    # Filter valid realizations
    valid = .!isnan.(rand_mean_fpts)

    #=== Phase 4: Summary table ===#
    println("\n" * "=" ^ 70)
    println("  SUMMARY TABLE")
    println("=" ^ 70)
    @printf("  %-30s %15s %20s\n", "Metric", "Lattice", "Randomized (mean±std)")
    println("  " * "-" ^ 65)
    @printf("  %-30s %15d %13.1f ± %.1f\n", "Shortest path (1→25)",
            L_min_lattice,
            mean(rand_shortest_paths[valid]),
            std(rand_shortest_paths[valid]))
    @printf("  %-30s %15.1f %13.1f ± %.1f\n", "Mean FPT",
            stats_lattice.mean,
            mean(rand_mean_fpts[valid]),
            std(rand_mean_fpts[valid]))
    @printf("  %-30s %15.1f %13.1f ± %.1f\n", "Median FPT",
            stats_lattice.median,
            mean(rand_median_fpts[valid]),
            std(rand_median_fpts[valid]))
    @printf("  %-30s %15.1f %13.1f ± %.1f\n", "Std FPT",
            stats_lattice.std,
            mean(rand_std_fpts[valid]),
            std(rand_std_fpts[valid]))
    @printf("  %-30s %15d %13.1f ± %.1f\n", "Min FPT",
            stats_lattice.min,
            mean(rand_shortest_paths[valid]),
            std(rand_shortest_paths[valid]))
    @printf("  %-30s %15.1f\n", "95th percentile FPT (lattice)",
            stats_lattice.q95)
    println("=" ^ 70)

    #=== Phase 5: Export CSVs ===#
    println("\n[Phase 5] Exporting data to CSV...")
    output_dir = joinpath(@__DIR__, "data")

    # --- CSV 1: FPT summary per graph ---
    fpt_rows = DataFrame(
        topology   = vcat(["lattice"], ["randomized_$(lpad(i,2,'0'))" for i in 1:n_random]),
        mean_fpt   = vcat([stats_lattice.mean], rand_mean_fpts),
        median_fpt = vcat([stats_lattice.median], rand_median_fpts),
        std_fpt    = vcat([stats_lattice.std], rand_std_fpts),
        shortest_path = vcat([L_min_lattice], rand_shortest_paths),
    )
    CSV.write(joinpath(output_dir, "fpt_summary.csv"), fpt_rows)
    println("  Saved fpt_summary.csv")

    #=== Phase 6: Calibration plot (FPT shift vs. topology effect) ===#
    println("\n[Phase 6] Generating calibration data...")
    fig6_dir = joinpath(@__DIR__, "fig6")

    # Load simulation parameters and power-law fits
    params = CSV.read(joinpath(fig6_dir, "params.csv"), DataFrame)
    lattice_fits = CSV.read(joinpath(fig6_dir, "lattice_powerlaw_fits.csv"), DataFrame)
    randomized_fits = CSV.read(joinpath(fig6_dir, "randomized_powerlaw_fits.csv"), DataFrame)

    # Filter to outflow node (chemostat_id == 25) and join with k_d
    lat_out = innerjoin(
        filter(r -> r.chemostat_id == 25, lattice_fits),
        select(params, :sim_number, :diffusion_rate),
        on = :sim_number
    )
    rand_out = innerjoin(
        filter(r -> r.chemostat_id == 25, randomized_fits),
        select(params, :sim_number, :diffusion_rate),
        on = :sim_number
    )

    # Sort by diffusion rate
    sort!(lat_out, :diffusion_rate)
    sort!(rand_out, :diffusion_rate)

    # FPT shift factor: 26% reduction in FPT ≈ multiplying k_d by 1.26
    fpt_shift_factor = stats_lattice.mean / mean(rand_mean_fpts[valid])

    # Create shifted lattice: interpolate lattice exponent at k_d * shift_factor
    # (this is what the lattice exponent would be if we only accounted for FPT change)
    lat_kd = lat_out.diffusion_rate
    lat_alpha = lat_out.slope
    shifted_kd = lat_kd .* fpt_shift_factor
    # Interpolate: for each shifted k_d, find the lattice α at that k_d
    shifted_alpha = zeros(length(lat_kd))
    for i in eachindex(lat_kd)
        target_kd = lat_kd[i] * fpt_shift_factor
        # Find bracketing indices in the original lattice data
        idx = searchsortedlast(lat_kd, target_kd)
        if idx < 1
            shifted_alpha[i] = lat_alpha[1]
        elseif idx >= length(lat_kd)
            shifted_alpha[i] = lat_alpha[end]
        else
            # Linear interpolation in log space
            t = (log10(target_kd) - log10(lat_kd[idx])) /
                (log10(lat_kd[idx+1]) - log10(lat_kd[idx]))
            shifted_alpha[i] = lat_alpha[idx] + t * (lat_alpha[idx+1] - lat_alpha[idx])
        end
    end

    # Compute gaps
    n_sims = length(lat_kd)
    rand_alpha = rand_out.slope[1:n_sims]
    total_gap = lat_alpha .- rand_alpha              # lattice α - randomized α (positive = lattice less steep)
    fpt_gap = lat_alpha .- shifted_alpha             # portion explained by FPT shift

    # Smooth with moving average (window = 7) to reduce noise
    function moving_avg(x, w)
        n = length(x)
        out = similar(x)
        for i in 1:n
            lo = max(1, i - w ÷ 2)
            hi = min(n, i + w ÷ 2)
            out[i] = mean(x[lo:hi])
        end
        return out
    end
    smooth_total = moving_avg(total_gap, 7)
    smooth_fpt = moving_avg(fpt_gap, 7)
    log_kd = log10.(lat_kd)

    # --- CSV 2: Calibration data ---
    cal_df = DataFrame(
        log_kd           = log_kd,
        lat_alpha        = lat_alpha,
        rand_alpha       = rand_alpha,
        shifted_alpha    = shifted_alpha,
        total_gap        = total_gap,
        fpt_gap          = fpt_gap,
        total_gap_smooth = smooth_total,
        fpt_gap_smooth   = smooth_fpt,
    )
    CSV.write(joinpath(output_dir, "calibration.csv"), cal_df)
    println("  Saved calibration.csv")

    # --- CSV 3: Metadata ---
    meta_df = DataFrame(
        fpt_shift_factor = [fpt_shift_factor],
        lattice_mean_fpt = [stats_lattice.mean],
        randomized_mean_fpt = [mean(rand_mean_fpts[valid])],
    )
    CSV.write(joinpath(output_dir, "meta.csv"), meta_df)
    println("  Saved meta.csv")

    t_elapsed = time() - t_start
    println("\n[Done] Total runtime: $(round(t_elapsed, digits=1)) seconds")
end

main()
