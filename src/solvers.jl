#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# src/solvers.jl
# Handles conversion of the JuMP Model into a format that can be passed
# through the MathProgBase interface to solvers, and ongoing updating of
# that representation if supported by the solver.
#############################################################################

# Analyze a JuMP Model to determine its traits, and thus what solvers can
# be used to solve the problem
immutable ProblemTraits
    int::Bool  # has integer variables
    lin::Bool  # has only linear objectives and constraints
    qp ::Bool  # has a quadratic objective function
    qc ::Bool  # has a quadratic constraint
    nlp::Bool  # has general nonlinear objective or constraints
    soc::Bool  # has a second-order cone constraint
    sdp::Bool  # has an SDP constraint (or SDP variable bounds)
    sos::Bool  # has an SOS constraint
    conic::Bool  # has an SDP or SOC constraint
end
function ProblemTraits(m::Model)
    int = any(c-> !(c == :Cont || c == :Fixed), m.colCat)
    qp = !isempty(m.obj.qvars1)
    qc = !isempty(m.quadconstr)
    nlp = m.nlpdata !== nothing
    soc = !isempty(m.socconstr)
    # will need to change this when we add support for arbitrary variable cones
    sdp = !isempty(m.sdpconstr) || !isempty(m.varCones)
    sos = !isempty(m.sosconstr)
    ProblemTraits(int, !(qp|qc|nlp|soc|sdp|sos), qp, qc, nlp, soc, sdp, sos, soc|sdp)
end
function default_solver(traits::ProblemTraits)
    if traits.int | traits.sos
        MathProgBase.defaultMIPsolver
    elseif traits.sdp
        MathProgBase.defaultSDPsolver
    elseif traits.conic
        MathProgBase.defaultConicsolver
    elseif traits.qp | traits.qc
        MathProgBase.defaultQPsolver
    elseif traits.nlp
        MathProgBase.defaultNLPsolver
    else
        MathProgBase.defaultLPsolver
    end
end


function solve(m::Model; suppress_warnings=false,
                ignore_solve_hook=(m.solvehook===nothing), kwargs...)
    # If the user or an extension has provided a solve hook, call
    # that instead of solving the model ourselves
    if !ignore_solve_hook
        return m.solvehook(m; suppress_warnings=suppress_warnings, kwargs...)
    end

    # Clear warning counters
    m.getvalue_counter = 0
    m.operator_counter = 0

    # Remember if the solver was initially unset so we can restore
    # it to be unset later
    unset = m.solver == UnsetSolver()

    # Analyze the problems traits to determine what solvers we can use
    traits = ProblemTraits(m)

    # Build the MathProgBase model from the JuMP model
    buildInternalModel(m, traits, suppress_warnings=suppress_warnings)

    # If the model is a general nonlinear, use different logic in
    # nlp.jl to solve the problem
    traits.nlp && return solvenlp(m, traits, suppress_warnings=suppress_warnings)

    # Solve the problem
    MathProgBase.optimize!(m.internalModel)
    stat = MathProgBase.status(m.internalModel)

    # Extract solution from the solver
    numRows, numCols = length(m.linconstr), m.numCols
    m.objVal = NaN
    m.colVal = fill(NaN, numCols)
    m.linconstrDuals = Array(Float64, 0)

    if stat == :Optimal
        # If we think dual information might be available, try to get it
        # If not, return an array of the correct length
        # TODO: support conic duals
        if !(traits.int | traits.sos | traits.conic)
            m.redCosts = try
                MathProgBase.getreducedcosts(m.internalModel)[1:numCols]
            catch
                fill(NaN, numCols)
            end

            m.linconstrDuals = try
                MathProgBase.getconstrduals(m.internalModel)[1:numRows]
            catch
                fill(NaN, numRows)
            end
        end
    else
        # Problem was not solved to optimality, attempt to extract useful
        # information anyway
        suppress_warnings || warn("Not solved to optimality, status: $stat")
        # Some solvers provide infeasibility rays (dual) or unbounded
        # rays (primal) for linear problems. Store these as the solution
        # if the exist.
        if traits.lin
            if stat == :Infeasible
                m.linconstrDuals = try
                    infray = MathProgBase.getinfeasibilityray(m.internalModel)
                    @assert length(infray) == numRows
                    infray
                catch
                    suppress_warnings || warn("Infeasibility ray (Farkas proof) not available")
                    fill(NaN, numRows)
                end
            elseif stat == :Unbounded
                m.colVal = try
                    unbdray = MathProgBase.getunboundedray(m.internalModel)
                    @assert length(unbdray) == numCols
                    unbdray
                catch
                    suppress_warnings || warn("Unbounded ray not available")
                    fill(NaN, numCols)
                end
            end
        end
    end

    # If the problem was solved, or if it terminated prematurely, try
    # to extract a solution anyway. This commonly occurs when a time
    # limit or tolerance is set (:UserLimit)
    if !(stat == :Infeasible || stat == :Unbounded)
        try
            objVal = MathProgBase.getobjval(m.internalModel) + m.obj.aff.constant
            colVal = MathProgBase.getsolution(m.internalModel)[1:numCols]
            # Don't corrupt the answers if one of the above two calls fails
            m.objVal = objVal
            m.colVal = colVal
        end
    end

    # The MathProgBase interface defines a conic problem to always be
    # a minimization problem, so we need to flip the objective before
    # reporting it to the user
    if traits.conic && m.objSense == :Max
        m.objVal *= -1
    end

    # If the solver was initially not set, we will restore this status
    # and drop the internal MPB model. This is important for the case
    # where the solver used changes between solves because the user
    # has changed the problem class (e.g. LP to MILP)
    if unset
        m.solver = UnsetSolver()
        if traits.int
            m.internalModelLoaded = false
        end
    end

    # Return the solve status
    stat
end

# Converts the JuMP Model into a MathProgBase model based on the
# traits of the model
function buildInternalModel(m::Model, traits=ProblemTraits(m);
                            suppress_warnings=false)
    # Set solver based on the model's traits if it hasn't provided
    if isa(m.solver, UnsetSolver)
        m.solver = default_solver(traits)
    end

    # If the model is nonlinear, use different logic in nlp.jl
    # to build the problem
    traits.nlp && return _buildInternalModel_nlp(m, traits)

    if traits.conic
        # If the problem is conic then use only the objective
        # coefficients from prepProblemBounds
        f,_,_ = prepProblemBounds(m)

        # The conic MPB interface defines conic problems as
        # always being minimization problems, so flip if needed
        m.objSense == :Max && scale!(f, -1.0)

        # Obtain a fresh MPB model for the solver
        # If the problem is conic, we rebuild the problem from
        # scratch every time
        m.internalModel = MathProgBase.model(m.solver)

        # Build up the LHS, RHS and cones from the JuMP Model...
        A, b, var_cones, con_cones = conicconstraintdata(m)
        # ... and pass to the solver
        MathProgBase.loadconicproblem!(m.internalModel, f, A, b, con_cones, var_cones)
    else
        # Extract objective coefficients and linear constraint bounds
        f, rowlb, rowub = prepProblemBounds(m)
        # If we already have an MPB model for the solver...
        if m.internalModelLoaded
            # ... and if the solver supports updating bounds/objective
            if applicable(MathProgBase.setvarLB!, m.internalModel, m.colLower) &&
               applicable(MathProgBase.setvarUB!, m.internalModel, m.colUpper) &&
               applicable(MathProgBase.setconstrLB!, m.internalModel, rowlb) &&
               applicable(MathProgBase.setconstrUB!, m.internalModel, rowub) &&
               applicable(MathProgBase.setobj!, m.internalModel, f) &&
               applicable(MathProgBase.setsense!, m.internalModel, m.objSense)
                MathProgBase.setvarLB!(m.internalModel, m.colLower)
                MathProgBase.setvarUB!(m.internalModel, m.colUpper)
                MathProgBase.setconstrLB!(m.internalModel, rowlb)
                MathProgBase.setconstrUB!(m.internalModel, rowub)
                MathProgBase.setobj!(m.internalModel, f)
                MathProgBase.setsense!(m.internalModel, m.objSense)
            else
                # The solver doesn't support changing bounds/objective
                # We need to build the model from scratch
                if !suppress_warnings
                    Base.warn_once("Solver does not appear to support hot-starts. Model will be built from scratch.")
                end
                m.internalModelLoaded = false
            end
        end
        # If we don't already have a MPB model
        if !m.internalModelLoaded
            # Obtain a fresh MPB model for the solver
            m.internalModel = MathProgBase.model(m.solver)
            # Construct a LHS matrix from the linear constraints
            A = prepConstrMatrix(m)

            # If we have either:
            #   1) A solver that does not support the loadproblem! interface, or
            #   2) A QCP and a solver that does not support the addquadconstr! interface,
            # wrap everything in a ConicSolverWrapper
            if !applicable(MathProgBase.loadproblem!, m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense) ||
                ( applicable(MathProgBase.supportedcones, m.solver) && # feel like this should have a && traits.qc as well...
                  !method_exists(MathProgBase.addquadconstr!, (typeof(m.internalModel), Vector{Int}, Vector{Float64}, Vector{Int}, Vector{Int}, Vector{Float64}, Char, Float64)) &&
                  :SOC in MathProgBase.supportedcones(m.solver) )

                m.internalModel = MathProgBase.model(MathProgBase.ConicSolverWrapper(m.solver))
            end

            # Load the problem data into the model...
            MathProgBase.loadproblem!(m.internalModel, A, m.colLower, m.colUpper, f, rowlb, rowub, m.objSense)
            # ... and add quadratic and SOS constraints separately
            addQuadratics(m)
            addSOS(m)
        end
        # Update solver callbacks, if any
        registercallbacks(m)
    end

    # Update the type of each variable
    if applicable(MathProgBase.setvartype!, m.internalModel, Symbol[])
        colCats = vartypes_without_fixed(m)
        MathProgBase.setvartype!(m.internalModel, colCats)
    elseif traits.int
        # Solver that do not implement anything other than continuous
        # variables do not need to implement this method, so throw an
        # error if the model has anything but continuous
        error("Solver does not support discrete variables")
    end

    # Provide a primal solution to the solve, if the problem is integer
    # and the user has provided one.
    # TODO: change this so you can warm start continuous problems?
    if traits.int && !all(isnan(m.colVal))
        if applicable(MathProgBase.setwarmstart!, m.internalModel, m.colVal)
            MathProgBase.setwarmstart!(m.internalModel, m.colVal)
        else
            suppress_warnings || Base.warn_once("Solver does not appear to support providing initial feasible solutions.")
        end
    end

    # Some solvers need to have an explicit "update" phase, e.g. Gurobi
    if applicable(MathProgBase.updatemodel!, m.internalModel)
        MathProgBase.updatemodel!(m.internalModel)
    end

    # Record that we have a MPB model constructed
    m.internalModelLoaded = true
    nothing
end

# Add the quadratic part of the objective and all quadratic constraints
# to the internal MPB model
function addQuadratics(m::Model)
    # The objective function is always a quadratic expression, but
    # may have no quadratic terms (i.e. be just affine)
    if length(m.obj.qvars1) != 0
        # Check that no coefficients are NaN/Inf
        assert_isfinite(m.obj)
        # Check that quadratic term variables belong to this model
        # Affine portion is checked in prepProblemBounds
        if !(verify_ownership(m, m.obj.qvars1) &&
                verify_ownership(m, m.obj.qvars2))
            error("Variable not owned by model present in objective")
        end
        # Check for solver support for quadratic objectives happens in MPB
        MathProgBase.setquadobjterms!(m.internalModel,
            Cint[v.col for v in m.obj.qvars1],
            Cint[v.col for v in m.obj.qvars2], m.obj.qcoeffs)
    end

    # Add quadratic constraint to solver
    const sensemap = @compat Dict(:(<=) => '<', :(>=) => '>', :(==) => '=')
    for k in 1:length(m.quadconstr)
        qconstr = m.quadconstr[k]::QuadConstraint
        if !haskey(sensemap, qconstr.sense)
            error("Invalid sense for quadratic constraint")
        end
        s = sensemap[qconstr.sense]

        terms::QuadExpr = qconstr.terms
        # Check that no coefficients are NaN/Inf
        assert_isfinite(terms)
        # Check that quadratic and affine term variables belong to this model
        if !(verify_ownership(m, terms.qvars1) &&
                verify_ownership(m, terms.qvars2) &&
                verify_ownership(m, terms.aff.vars))
            error("Variable not owned by model present in quadratic constraint")
        end
        # Extract indices for MPB, and add the constraint (if we can)
        affidx  = Cint[v.col for v in terms.aff.vars]
        var1idx = Cint[v.col for v in terms.qvars1]
        var2idx = Cint[v.col for v in terms.qvars2]
        if applicable(MathProgBase.addquadconstr!, m.internalModel, affidx, terms.aff.coeffs, var1idx, var2idx, terms.qcoeffs, s, -terms.aff.constant)
            MathProgBase.addquadconstr!(m.internalModel,
                affidx, terms.aff.coeffs,           # aᵀx +
                var1idx, var2idx, terms.qcoeffs,    # xᵀQx
                s, -terms.aff.constant)             # ≤/≥ b
        else
            error("Solver does not support quadratic constraints")
        end
    end
    nothing
end

function addSOS(m::Model)
    for i in 1:length(m.sosconstr)
        sos = m.sosconstr[i]
        indices = Int[v.col for v in sos.terms]
        if sos.sostype == :SOS1
            if applicable(MathProgBase.addsos1!, m.internalModel, indices, sos.weights)
                MathProgBase.addsos1!(m.internalModel, indices, sos.weights)
            else
                error("Solver does not support SOS constraints")
            end
        elseif sos.sostype == :SOS2
            if applicable(MathProgBase.addsos2!, m.internalModel, indices, sos.weights)
                MathProgBase.addsos2!(m.internalModel, indices, sos.weights)
            else
                error("Solver does not support SOS constraints")
            end
        end
    end
end

# Returns coefficients for the affine part of the objective and the
# affine constraint lower and upper bounds, all as dense vectors
function prepProblemBounds(m::Model)

    # Create dense objective vector
    objaff::AffExpr = m.obj.aff
    # Check that no coefficients are NaN/Inf
    assert_isfinite(objaff)
    if !verify_ownership(m, objaff.vars)
        error("Variable not owned by model present in objective")
    end
    f = zeros(m.numCols)
    @inbounds for ind in 1:length(objaff.vars)
        f[objaff.vars[ind].col] += objaff.coeffs[ind]
    end

    # Create dense affine constraint bound vectors
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # -Inf means no lower bound, +Inf means no upper bound
    rowlb = fill(-Inf, numRows)
    rowub = fill(+Inf, numRows)
    @inbounds for ind in 1:numRows
        rowlb[ind] = linconstr[ind].lb
        rowub[ind] = linconstr[ind].ub
    end

    return f, rowlb, rowub
end

# Convert all the affine constraints into a sparse column-wise
# matrix of coefficients.
function prepConstrMatrix(m::Model)

    # We want a column-wise sparse A matrix, but we have the
    # data in a row-wise format. The solution is to build up
    # a row-wise sparse matrix, then use the built-in transpose
    # operator to convert to a column-wise matrix.
    linconstr = m.linconstr::Vector{LinearConstraint}
    numRows = length(linconstr)
    # Calculate the maximum number of nonzeros
    # The actual number may be less because of cancelling or
    # zero-coefficient terms
    nnz = 0
    for c in 1:numRows
        nnz += length(linconstr[c].terms.coeffs)
    end
    # The non-zero values
    rownzval = Array(Float64,nnz)
    # The column for each non-zero
    colval = Array(Int,nnz)
    # The index of the beginning of each row in rownzval
    rowptr = Array(Int,numRows+1)

    # Fill it up!
    # Number of nonzeros seen so far
    nnz = 0
    # Maintain a data structure for collapsing down duplicate
    # terms in each constraint
    tmprow = IndexedVector(Float64,m.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for c in 1:numRows
        rowptr[c] = nnz + 1
        # Check that no coefficients are NaN/Inf
        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars   = linconstr[c].terms.vars
        # Check that variables belong to this model
        if !verify_ownership(m, vars)
            error("Variable not owned by model present in a constraint")
        end
        # Add all terms into the IndexedVector, which will
        # combine any duplicate terms
        @inbounds for ind in 1:length(coeffs)
            addelt!(tmprow, vars[ind].col, coeffs[ind])
        end
        # Extract all the terms from the IndexedVector
        # Each variable will appear at most once
        @inbounds for i in 1:tmprow.nnz
            idx = tmpnzidx[i]   # The column index
            elt = tmpelts[idx]  # The coefficient
            # Do not pass zero coefficients through to the solver
            elt == 0.0 && continue
            # Store this nonzero value in the row-wise matrix
            nnz += 1
            colval[nnz] = idx
            rownzval[nnz] = elt
        end
        # Reset the IndexedVector for the next constraint
        empty!(tmprow)
    end
    # The last value in rowptr is the 1 after the last term added
    rowptr[numRows+1] = nnz + 1

    # Build the row-wise sparse matrix (although we pretend it is
    # a column-wise sparse matrix, it doesn't make a difference)
    rowmat = SparseMatrixCSC(m.numCols, numRows, rowptr, colval, rownzval)
    # Note that rowmat doesn't have sorted indices, so technically doesn't
    # follow SparseMatrixCSC format. But it's safe to take the transpose.
    A = rowmat'
end

function vartypes_without_fixed(m::Model)
    colCats = copy(m.colCat)
    for i in 1:length(colCats)
        if colCats[i] == :Fixed
            @assert m.colLower[i] == m.colUpper[i]
            colCats[i] = :Cont
        end
    end
    return colCats
end

function collect_expr!(m, tmprow, terms::AffExpr)
    empty!(tmprow)
    assert_isfinite(terms)
    coeffs = terms.coeffs
    vars = terms.vars
    # collect duplicates
    for ind in 1:length(coeffs)
        if !is(vars[ind].m, m)
            error("Variable not owned by model present in constraints")
        end
        addelt!(tmprow,vars[ind].col, coeffs[ind])
    end
    tmprow
end

function conicconstraintdata(m::Model)
    var_cones = Any[cone for cone in m.varCones]
    con_cones = Any[]
    nnz = 0

    # find starting column indices for sdp matrices
    numSDPRows = 0
    numSymRows = 0
    for c in m.sdpconstr
        n = size(c.terms,1)
        @assert n == size(c.terms,2)
        @assert ndims(c.terms) == 2
        numSDPRows += convert(Int, n*(n+1)/2)
        for i in 1:n, j in i:n
            nnz += length(c.terms[i,j].coeffs)
        end
        if !issym(c.terms)
            # symmetry constraints
            numSymRows += convert(Int, n*(n-1)/2)
        end
    end

    soc_cones  = Any[]
    rsoc_cones = Any[]
    numQuadRows = 0
    for qconstr in m.quadconstr
        q = copy(qconstr.terms)
        if qconstr.sense == :(>=)
            q *= -1
        end
        if !(isempty(q.aff.vars) && q.aff.constant == 0)
            error("Quadratic constraint $qconstr must be in second-order cone form")
        end
        n_pos_on_diag = 0
        off_diag_idx  = 0
        neg_diag_idx  = 0
        n = length(q.qvars1)
        for i in 1:n
            if q.qvars1[i].col == q.qvars2[i].col
                if q.qcoeffs[i] == 1
                    n_pos_on_diag += 1
                elseif q.qcoeffs[i] == -1
                    neg_diag_idx == off_diag_idx == 0 || error("Invalid SOC constraint $qconstr")
                    neg_diag_idx = i
                end
            else
                if q.qcoeffs[i] == -1
                    neg_diag_idx == off_diag_idx == 0 || error("Invalid rotated SOC constraint $qconstr")
                    off_diag_idx = i
                end
            end
        end
        cone = Array(Int, n)
        if n_pos_on_diag == n-1 && neg_diag_idx > 0
            cone[1] = q.qvars1[neg_diag_idx].col
            for i in 1:(neg_diag_idx-1); cone[i+1] = q.qvars1[i].col; end
            for i in (neg_diag_idx+1):n; cone[i]   = q.qvars1[i].col; end
            push!(soc_cones, cone)
        elseif n_pos_on_diag == n-1 && off_diag_idx > 0
            cone[1] = q.qvars1[off_diag_idx].col
            cone[2] = q.qvars2[off_diag_idx].col
            for i in 1:(off_diag_idx-1); cone[i+2] = q.qvars1[i].col; end
            for i in (off_diag_idx+1):n; cone[i+1] = q.qvars1[i].col; end
            push!(rsoc_cones, cone)
        else
            error("Quadratic constraint $qconstr is not conic representable")
        end
        numQuadRows += length(cone)
    end

    linconstr = m.linconstr::Vector{LinearConstraint}
    numLinRows = length(linconstr)
    numBounds = 0
    nonNeg  = Int[]
    nonPos  = Int[]
    free    = Int[]
    zeroVar = Int[]
    for i in 1:m.numCols
        seen = false
        lb, ub = m.colLower[i], m.colUpper[i]
        for (_,cone) in m.varCones
            if i in cone
                seen = true
                @assert lb == -Inf && ub == Inf
                break
            end
        end

        if !seen
            if !(lb == 0 || lb == -Inf)
                numBounds += 1
            end
            if !(ub == 0 || ub == Inf)
                numBounds += 1
            end
            if lb == 0 && ub == 0
                push!(zeroVar, i)
            elseif lb == 0
                push!(nonNeg, i)
            elseif ub == 0
                push!(nonPos, i)
            else
                push!(free, i)
            end
        end
    end

    if !isempty(zeroVar)
        push!(var_cones, (:Zero,zeroVar))
    end
    if !isempty(nonNeg)
        push!(var_cones, (:NonNeg,nonNeg))
    end
    if !isempty(nonPos)
        push!(var_cones, (:NonPos,nonPos))
    end
    if !isempty(free)
        push!(var_cones, (:Free,free))
    end

    nnz += numBounds
    for c in 1:numLinRows
        nnz += length(linconstr[c].terms.coeffs)
    end

    numSOCRows = 0
    for con in m.socconstr
        numSOCRows += length(con.normexpr.norm.terms) + 1
    end
    numRows = numLinRows + numBounds + numQuadRows + numSOCRows + numSDPRows + numSymRows

    b = Array(Float64, numRows)

    I = Int[]
    J = Int[]
    V = Float64[]
    @compat sizehint!(I, nnz)
    @compat sizehint!(J, nnz)
    @compat sizehint!(V, nnz)

    # Fill it up
    nnz = 0
    tmprow = IndexedVector(Float64,m.numCols)
    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    nonneg_rows = Int[]
    nonpos_rows = Int[]
    eq_rows     = Int[]
    for c in 1:numLinRows
        if linconstr[c].lb == -Inf
            b[c] = linconstr[c].ub
            push!(nonneg_rows, c)
        elseif linconstr[c].ub == Inf
            b[c] = linconstr[c].lb
            push!(nonpos_rows, c)
        elseif linconstr[c].lb == linconstr[c].ub
            b[c] = linconstr[c].lb
            push!(eq_rows, c)
        else
            error("We currently do not support ranged constraints with conic solvers")
        end

        assert_isfinite(linconstr[c].terms)
        coeffs = linconstr[c].terms.coeffs
        vars = linconstr[c].terms.vars
        # collect duplicates
        for ind in 1:length(coeffs)
            if !is(vars[ind].m, m)
                error("Variable not owned by model present in constraints")
            end
            addelt!(tmprow,vars[ind].col, coeffs[ind])
        end
        nnz = tmprow.nnz
        append!(I, fill(c, nnz))
        indices = tmpnzidx[1:nnz]
        append!(J, indices)
        append!(V, tmpelts[indices])
        empty!(tmprow)
    end

    c = numLinRows
    for idx in 1:m.numCols
        lb = m.colLower[idx]
        if !(lb == 0 || lb == -Inf)
            nnz += 1
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = lb
            push!(nonpos_rows, c)
        end
        ub = m.colUpper[idx]
        if !(ub == 0 || ub == Inf)
            c   += 1
            push!(I, c)
            push!(J, idx)
            push!(V, 1.0)
            b[c] = ub
            push!(nonneg_rows, c)
        end
    end

    if !isempty(nonneg_rows)
        push!(con_cones, (:NonNeg,nonneg_rows))
    end
    if !isempty(nonpos_rows)
        push!(con_cones, (:NonPos,nonpos_rows))
    end
    if !isempty(eq_rows)
        push!(con_cones, (:Zero,eq_rows))
    end
    @assert c == numLinRows + numBounds

    for cone in soc_cones
        n = length(cone)
        rng = (c+1):(c+n)
        append!(I, rng)
        append!(J, copy(cone))
        append!(V, [-1.0; ones(n-1)])
        push!(con_cones, (:SOC,rng))
        b[rng] = 0
        c += n
    end
    for cone in rsoc_cones
        n = length(cone)
        rng = (c+1):(c+n)
        append!(I, rng)
        append!(J, copy(cone))
        append!(V, [-1.0; -1.0; ones(n-2)])
        push!(con_cones, (:SOCRotated,rng))
        b[rng] = 0
        c += n
    end
    @assert c == numLinRows + numBounds + numQuadRows

    tmpelts = tmprow.elts
    tmpnzidx = tmprow.nzidx
    for con in m.socconstr
        expr = con.normexpr
        c += 1
        soc_start = c
        collect_expr!(m, tmprow, expr.aff)
        nnz = tmprow.nnz
        indices = tmpnzidx[1:nnz]
        append!(I, fill(c, nnz))
        append!(J, indices)
        append!(V, tmpelts[indices])
        b[c] = -expr.aff.constant
        for term in expr.norm.terms
            c += 1
            collect_expr!(m, tmprow, term)
            nnz = tmprow.nnz
            indices = tmpnzidx[1:nnz]
            append!(I, fill(c, nnz))
            append!(J, indices)
            append!(V, -expr.coeff*tmpelts[indices])
            b[c] = expr.coeff*term.constant
        end
        push!(con_cones, (:SOC, soc_start:c))
    end
    @assert c == numLinRows + numBounds + numQuadRows + numSOCRows

    for con in m.sdpconstr
        sdp_start = c + 1
        n = size(con.terms,1)
        for i in 1:n, j in i:n
            c += 1
            terms::AffExpr = con.terms[i,j]
            collect_expr!(m, tmprow, terms)
            nnz = tmprow.nnz
            indices = tmpnzidx[1:nnz]
            append!(I, fill(c, nnz))
            append!(J, indices)
            append!(V, -tmpelts[indices])
            b[c] = terms.constant
        end
        push!(con_cones, (:SDP, sdp_start:c))
        if !issym(con.terms)
            sym_start = c + 1
            # add linear symmetry constraints
            for i in 1:n, j in 1:(i-1)
                c += 1
                collect_expr!(m, tmprow, con.terms[i,j] - con.terms[j,i])
                nnz = tmprow.nnz
                indices = tmpnzidx[1:nnz]
                append!(I, fill(c, nnz))
                append!(J, indices)
                append!(V, tmpelts[indices])
                b[c] = 0
            end
            push!(con_cones, (:Zero, sym_start:c))
        end
    end
    @assert c == numRows

    A = sparse(I, J, V, numRows, m.numCols)
    # @show full(A), b
    # @show var_cones, con_cones

    # TODO: uncomment these lines when they work with Mosek
    # supported = MathProgBase.supportedcones(m.internalModel)
    # @assert (:NonNeg in supported) && (:NonPos in supported) && (:Free in supported) && (:SDP in supported)
    A, b, var_cones, con_cones
end


# returns (unsorted) column indices and coefficient terms for merged vector
# assume that v is zero'd
function merge_duplicates{CoefType,IntType<:Integer}(::Type{IntType},aff::GenericAffExpr{CoefType,Variable}, v::IndexedVector{CoefType}, m::Model)
    resize!(v, m.numCols)
    for ind in 1:length(aff.coeffs)
        var = aff.vars[ind]
        is(var.m, m) || error("Variable does not belong to this model")
        addelt!(v, aff.vars[ind].col, aff.coeffs[ind])
    end
    indices = Array(IntType,v.nnz)
    coeffs = Array(CoefType,v.nnz)
    for i in 1:v.nnz
        idx = v.nzidx[i]
        indices[i] = idx
        coeffs[i] = v.elts[idx]
    end
    empty!(v)

    return indices, coeffs

end
