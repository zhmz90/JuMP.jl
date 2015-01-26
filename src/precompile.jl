# JuMP.jl
precompile(Model, ())
precompile(Variable, (Model,Float64,Float64,Symbol,ASCIIString,Float64))
precompile(Variable, (Model,Float64,Float64,Symbol,UTF8String,Float64))
precompile(Variable, (Model,Int))
precompile(isequal, (Variable,Variable))
precompile(LinearConstraint, (AffExpr,Float64,Float64))
precompile(QuadConstraint, (QuadExpr,Symbol))

precompile(addConstraint, (Model, LinearConstraint))
precompile(addConstraint, (Model, QuadConstraint))
precompile(setObjectiveSense, (Model, Symbol))
precompile(setObjective, (Model, Symbol, Variable))
precompile(setObjective, (Model, Symbol, AffExpr))

precompile(verify_ownership, (Model, Vector{Variable}))

precompile(AffExpr, ())
precompile(push!, (AffExpr, Float64, Variable))
precompile(append!, (AffExpr, AffExpr))


# macros
if VERSION < v"0.4.0-"
    precompile(timescoef, (Expr,))
    precompile(timesvar, (Expr,))
    precompile(addToExpression, (AffExpr,Float64,Variable))
    precompile(addToExpression, (AffExpr,Int,Variable))
    precompile(addToExpression, (AffExpr,Float64,Float64))
    precompile(addToExpression, (AffExpr,Variable,Variable))
    precompile(addToExpression, (AffExpr,Float64,AffExpr))
    precompile(addToExpression, (AffExpr,Int,AffExpr))
    precompile(addToExpression, (AffExpr,AffExpr,Variable))
    precompile(addToExpression, (AffExpr,Variable,AffExpr))
    precompile(addToExpression, (AffExpr,Float64,QuadExpr))
    precompile(addToExpression, (AffExpr,Int,QuadExpr))
    precompile(addToExpression, (QuadExpr,Float64,Variable))
    precompile(addToExpression, (QuadExpr,Int,Variable))
    precompile(addToExpression, (QuadExpr,Float64,Float64))
    precompile(addToExpression, (QuadExpr,Float64,AffExpr))
    precompile(addToExpression, (QuadExpr,Int,AffExpr))
    precompile(addToExpression, (QuadExpr,QuadExpr,Number))
end

# operators
for sgn in [<=,==,>=]
    operands = [AffExpr,QuadExpr,Variable,Float64,Int]
    for o1 in operands
        for o2 in operands
            precompile(sgn, (o1,o2))
        end
    end
end
precompile(*, (Variable, Variable))
precompile(^, (Variable, Int))

# solvers.jl
precompile(solve, (Model,))
precompile(addQuadratics, (Model,))
precompile(addSOS, (Model,))
precompile(prepProblemBounds, (Model,))
precompile(prepConstrMatrix, (Model,))
precompile(solveLP, (Model,))
precompile(solveMIP, (Model,))
precompile(buildInternalModel, (Model,))

# utils.jl
precompile(length, (IndexedVector{Float64},))
