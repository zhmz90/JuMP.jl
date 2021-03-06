#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# iis.jl
#
# Extract an Irreducible Inconsistent Subsystem (IIS)
# from an infeasible model.
# This code is specific to Gurobi.
#############################################################################

using JuMP, Gurobi
import MathProgBase

function print_iis_gurobi(m::Model)

    grb = MathProgBase.getrawsolver(getInternalModel(m))
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    println("Irreducible Inconsistent Subsystem (IIS)")
    println("Variable bounds:")
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(getLower(v), " <= ", getName(v), " <= ", getUpper(v))
        elseif iislb[i] != 0
            println(getName(v), " >= ", getLower(v))
        elseif iisub[i] != 0
            println(getName(v), " <= ", getUpper(v))
        end
    end


    println("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(ConstraintRef{LinearConstraint}(m, i))
        end
    end
    println("End of IIS")
end

# The function below produces the following output:
#=
Irreducible Inconsistent Subsystem (IIS)
Variable bounds:
x <= 1.0
y <= 1.0
Constraints:
x + y >= 3
End of IIS
=#

function iis_example()
    m = Model(solver=GurobiSolver())

    @defVar(m, 0 <= x <= 1)
    @defVar(m, 0 <= y <= 1)

    @addConstraint(m, x+y >= 3)

    status = solve(m)
    @assert status == :Infeasible

    print_iis_gurobi(m)
end
