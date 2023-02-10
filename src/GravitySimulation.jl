module GravitySimulation

include("UsingPackages.jl")

function GravityStep!(dPS::Vector{Float64}, PS::Vector{Float64}, Parameters::Vector{Float64}, t)
    M1, M2, M3, = Parameters

    dPS[1:9] = PS[10:18] # difference of position is velocity.

    R1 = PS[[1, 2, 3]]
    R2 = PS[[4, 5, 6]]
    R3 = PS[[7, 8, 9]]

    S12 = -1 / sum(abs2.(R1 - R2))^(3 / 2)
    S21 = S12
    S13 = -1 / sum(abs2.(R1 - R3))^(3 / 2)
    S31 = S13
    S23 = -1 / sum(abs2.(R2 - R3))^(3 / 2)
    S32 = S23

    dPS[10:12] = (S12 * M2) * (R1 - R2) .+ (S13 * M3) * (R1 - R3) # force on 1
    dPS[13:15] = (S21 * M1) * (R2 - R1) .+ (S23 * M3) * (R2 - R3) # force on 2
    dPS[16:18] = (S31 * M1) * (R3 - R1) .+ (S32 * M2) * (R3 - R2) # force on 3

    return dPS
end

function SolveTheProblem()

    Init_PS = zeros(Float64, 18)

    Init_PS[1] = 0.0
    Init_PS[4] = 10.0
    Init_PS[7] = 10.05

    Init_PS[11] = 0.0
    Init_PS[14] = 10.0
    Init_PS[17] = 2.0
    Init_PS[18] = 0.0

    M1, M2, M3 = 1000.0, 2.0, 0.0005
    Parameters = [M1, M2, M3]

    TimeAxis = collect(range(0, 10, 500))

    TargetProblem = ODEProblem(GravityStep!, Init_PS, (TimeAxis[begin], TimeAxis[end]), Parameters)
    Solution = solve(TargetProblem, Tsit5(), tstops = TimeAxis)

    Fig = Figure(resolution = (1600, 1600))
    Axe = Axis3(Fig[1,1])

    X = zeros(Float64, 3, length(Solution.t))
    Y = zeros(Float64, 3, length(Solution.t))
    Z = zeros(Float64, 3, length(Solution.t))

    for i in eachindex(Solution.t)
        PS = Solution.u[i]
        X[:, i] = PS[[1, 4, 7]]
        Y[:, i] = PS[[2, 5, 8]]
        Z[:, i] = PS[[3, 6, 9]]
    end

    for i in 1:3
        lines!(Axe, X[i, :], Y[i, :], Z[i, :])
    end

    wait(display(Fig))
end

end # module GravitySimulation
