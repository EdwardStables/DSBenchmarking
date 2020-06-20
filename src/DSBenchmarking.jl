module DSBenchmarking
using PyCall
using DirectSearch
using Dates
using Tables
using CSV

export RunSingleBenchmark, Benchmark, setup, run, csv_write, dims

const ROOT = "/home/eps116/FYP_Benchmarking/DS/testing/"

mutable struct Benchmark
    p::PyObject
    d::DSProblem

    times::Vector{Float64}
    test::String

    known_opt::Float64

    write::Bool
    file::String

    function Benchmark(test::String;file=false)  
        b = new()
        b.test = test

        b.times = Vector{Float64}()

        if file == false
            write = false
        else
            write = true
            file=file
        end

        return b
    end
end

function RunSingleBenchmark(test::String, poll::String, max_unsucc::Int, opt; limit=100 ,size=nothing)
    bm = Benchmark(test)

    #println("Setup")
    setup(bm, poll, max_unsucc, opt; size=size, limit=limit)
    #println("Run")
    run(bm)
    return bm
end

function result(bm::Benchmark)
    println("Done:")
    println("$(bm.test)")
    println("Result: $(bm.d.x_cost)")
    println("Iterations: $(bm.d.iteration)")
end

function setup(bm::Benchmark, poll::String, max_unsucc::Int, known_opt; size=nothing, limit=100)

    P = pyimport("cutest_interface")
    if size != nothing
        bm.p = P.problem(bm.test, user_dims=size)
    else
        bm.p = P.problem(bm.test)
    end



    bm.known_opt = known_opt
    
    if poll == "LT"
        bm.d = DSProblem(bm.p.n; max_unsuccessful=max_unsucc)
    elseif poll == "O"
        bm.d = DSProblem(bm.p.n, poll=OrthoMADS(bm.p.n); max_unsuccessful=max_unsucc)
    else
        error("gimme a real poll")
    end

    eq_cons = bm.p.is_eq_cons
    lb = bm.p.lb
    ub = bm.p.ub

    for i in 1:bm.p.m
        #if Bool(eq_cons[i])
        #    AddProgressiveConstraint(bm.d, x -> abs(bm.p.eval_con(x, i-1)))
        #else
            AddExtremeConstraint(bm.d, x -> lb[i] <= bm.p.eval_con(x, i-1) <= ub[i])
        #end
    end

    DS.BenchmarkKnownOptimum(bm.d, bm.known_opt)
    SetInitialPoint(bm.d, bm.p.x0)
    SetIterationLimit(bm.d, limit)
    SetObjective(bm.d, x -> bm.p.eval_obj(x))
end

function dims(bm::Benchmark)
    return bm.p.dims()
end

function run(bm::Benchmark)
    push!(bm.times, time())
    Optimize!(bm.d) 
    push!(bm.times, time())
end

function csv_write(bm)
    cols = ["Test Time", "Poll", "Test", "Variables", "Iterations", "Total Function Evaluations", "Function Evaluations Before Optimum", "Time", "Optimal Value"]

    
    runtime = bm.times[end]-bm.times[end-1]

    poll = typeof(bm.d.poll) == LTMADS{Float64} ? "LTMADS" : "OrthoMADS"

    row = [Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") poll bm.test bm.p.n bm.d.iteration bm.d.function_evaluations bm.d.optimum_found_after_evals runtime bm.d.x_cost]

    file = ROOT*"results/DS.csv"
    open(file, isfile(file) ? "a" : "w") do f
        CSV.write(f, Tables.table(row), header=cols, append=isfile(file))
    end
end

end # module
