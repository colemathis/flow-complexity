using DrWatson
@quickactivate :FlowComplexity
const FC = FlowComplexity

using DataFrames, Statistics, Plots, LaTeXStrings, JLD2
using Printf # for @sprintf
using StatsBase # for fit, histogram
using LsqFit # for curvefit

# using RCall
using CSV

################################################################################

function plot_timeseries()

    # create a time series
    
    # only keep the first chemostat
    ts = copy(timeseries)
    filter!(row -> row.sim_number == 1, ts)
    
    # take average over simulations
    gts = groupby(ts, [:time, :integer])
    gts = combine(gts, :frequency => mean => :frequency)
    
    p = plot(title="Time series for sim 1, avg over all chemostats", xlabel="Time", ylabel="Count",
             legendtitle="Integers", legend=:outerright)
    
    for i in 1:10
        f = filter(row -> row.integer == i, gts)
        plot!(p, f.time, f.frequency, label="$i")
    end
    
    save_fig("time-series.pdf")

end

################################################################################

function plot_multiplot()

    # create a multiplot for all chemostats of sim 1
    
    ts_multi = copy(timeseries)
    sim_no = 1
    filter!(row -> row.sim_number == sim_no, ts_multi)
    
    # p = plot(layout = grid(5, 5), size = (900, 900), legend=nothing, link=:both, ylim = (0, 1000))
    p = plot(layout = grid(5, 5), size = (900, 900), legend=nothing, link=:both, ylim = (0, 10000))
    
    nspecies = 10
    nchem = sim_array[sim_no].params[:N_reactors]
    for i in 1:nchem
        f = filter(row -> row.chemostat_id == i, ts_multi)
        for j in 1:nspecies
            g = filter(row -> row.integer == j, f)
            x = g.time
            y = g.frequency
            plot!(p, x, y, subplot = i)
        end
    end
    
    save_fig("multi.pdf")

end

################################################################################

function plot_histogram()

    # create a histogram    
    
    ts = copy(timeseries) # sim_number, time, chemostat_id, integer, frequency
    filter!(row -> row.time == 100.0, ts)
    filter!(row -> 0 < row.integer < 2500, ts)
    
    gts = groupby(ts, [:sim_number, :integer])
    gts = combine(gts, :frequency => mean => :frequency)
    sort!(gts, :integer)

    p = plot(title="Parameter sweep: in-flow", subtitle="Average populations over all chemostats",
             xlabel="Integers", ylabel="Frequency", ylim=(-5,200),
             legendtitle="In-flow", legend=:topright)
    
    for i in 1:2:10
        
        f = filter(row -> row.sim_number == i, gts)
        filter!(row -> row.frequency > 1, f)
        expanded_data = [x for (x, y) in zip(f.integer, f.frequency) for _ in 1:y]
        bin_counts, bin_edges = calculate_histogram(data = expanded_data,
                                                    xmin = 0,
                                                    xmax = maximum(f.integer),
                                                    nbins = 20)
        
        # Replace all 0 values in bin_counts with 1
        # bin_counts = max.(bin_counts, repeat([0.1], nbins-1))
        
        x = bin_edges
        y = bin_counts
        
        x = x[2:end]
        y = y[2:end]
        model, fitted_params = fit_power_law(x, y)
        
        l = sim_array[i].params[:inflow_mols]
        label = @sprintf("%.1e", l)
        
        scatter!(p, x, y, label="$label (data)")
        plot!(p, x -> model(x, fitted_params), label="$label (fit)")

    end
    
    save_fig("hist_all.pdf")
    
end

################################################################################

function calculate_histogram(; data, xmin, xmax, nbins)

    # linear bins
    linBins = LinRange(xmin, xmax, nbins)
    
    # exp bins
    l = log10(xmax)
    c = ceil(l)
    e = LinRange(0, c, nbins)
    expBins = 10.0 .^ e
    
    hist = fit(Histogram, data, linBins)
    bin_counts = hist.weights
    bin_edges = hist.edges[1][1:end-1]

    return bin_counts, bin_edges

end

################################################################################

function fit_power_law(x, y)
    
    model(x, p) = p[1] * (p[2] * x) .^ (-p[3])
    
    p0 = [50000.0, 1.0, 1.0]
    
    f = curve_fit(model, x, y, p0, maxIter=100000)
    fitted_params = coef(f)

    return model, fitted_params

end

################################################################################

function fit_exp_law(x, y)

    # model(x, p) = p[1] * exp.(-p[2] * x) .+ p[3]

    # ...

    return

end

################################################################################

function save_fig(name)

    mkpath("figs")
    fn = "figs/$name"
    savefig(fn)
    run(`imgcat $fn`)

end

################################################################################

params = DataFrame(CSV.File("data/params.csv"))
timeseries = DataFrame(CSV.File("data/timeseries.csv"))
@load "data/sim_array.jld2" sim_array

mkpath("figs")

plot_timeseries()
plot_multiplot()
plot_histogram()
