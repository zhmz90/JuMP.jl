#  Copyright 2015, Iain Dunning, Joey Huchette, Miles Lubin, and contributors
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# test/print.jl
# Testing for all pretty-printing-related functionality
#############################################################################
using JuMP, FactCheck
import JuMP.REPLMode, JuMP.IJuliaMode

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @fact sprint(print, obj) --> exp_str
        repl != :print && @fact sprint(show,  obj) --> exp_str
    else
        @fact sprint(writemime, "text/latex", obj) --> "\$\$ $exp_str \$\$"
    end
end


facts("[print] JuMPContainer{Variable}") do
    le, ge = JuMP.repl[:leq], JuMP.repl[:geq]
    m = Model()

    #------------------------------------------------------------------
    # Test bound printing
    context("bound printing") do
    @defVar(m,      bnd_free[2:5])
    @defVar(m,      bnd_lowb[2:5] >= 2)
    @defVar(m,      bnd_high[2:5] <= 5)
    @defVar(m, 2 <= bnd_both[2:5] <= 5)
    @defVar(m,      bnd_difflo[i=2:5] >= i)
    @defVar(m,      bnd_diffup[i=2:5] <= i)
    @defVar(m, i <= bnd_diffbo[i=2:5] <= 2i)
    @defVar(m, i <= bnd_difflo_with_up[i=2:5] <= 5)
    @defVar(m, 2 <= bnd_diffup_with_lo[i=2:5] <= i)

    io_test(REPLMode, bnd_free, "bnd_free[i] free for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_lowb, "bnd_lowb[i] $ge 2 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_high, "bnd_high[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_both, "2 $le bnd_both[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo, "bnd_difflo[i] $ge .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup, "bnd_diffup[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffbo, ".. $le bnd_diffbo[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo_with_up, ".. $le bnd_difflo_with_up[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup_with_lo, "2 $le bnd_diffup_with_lo[i] $le .. for all i in {2,3,4,5}")

    io_test(IJuliaMode, bnd_free, "bnd_free_{i} free \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_lowb, "bnd_lowb_{i} \\geq 2 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_high, "bnd_high_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_both, "2 \\leq bnd_both_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo, "bnd_difflo_{i} \\geq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup, "bnd_diffup_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffbo, ".. \\leq bnd_diffbo_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_difflo_with_up, ".. \\leq bnd_difflo_with_up_{i} \\leq 5 \\quad\\forall i \\in \\{2,3,4,5\\}")
    io_test(IJuliaMode, bnd_diffup_with_lo, "2 \\leq bnd_diffup_with_lo_{i} \\leq .. \\quad\\forall i \\in \\{2,3,4,5\\}")
    end

    #------------------------------------------------------------------
    # Test index set printing
    context("index set printing") do
    @defVar(m, rng_unit1[1:10])  # Array{Variable}
    @defVar(m, rng_unit2[-2:3])  # JuMPArray
    @defVar(m, rng_unit3[[1:10;]])  # JuMPDict
    @defVar(m, rng_step1[1:2:10])
    @defVar(m, rng_step2[-2:5:10])
    @defVar(m, rng_step3[1:5:3])
    @defVar(m, rng_step4[0:2:2])
    @defVar(m, arr_1[[:a,:b,:c]])
    @defVar(m, arr_2[[:a,1,"test"]])
    @defVar(m, arr_3[[:apple,:banana,:carrot,:diamonds]])
    @defVar(m, rng2_1[1:10,[:a,:b,:c]])
    @defVar(m, tri_1[i=1:3,j=i:3])
    @defVar(m, tri_2[i=1:3,j=-i])
    @defVar(m, tri_3[(i,j)=[(i,i+2) for i in 1:5],k=i:j])

    io_test(REPLMode, rng_unit1, "rng_unit1[i] free for all i in {1,2..9,10}")
    io_test(REPLMode, rng_unit2, "rng_unit2[i] free for all i in {-2,-1..2,3}")
    io_test(REPLMode, rng_unit3, "rng_unit3[i] free for all i in {1,2..9,10}")
    io_test(REPLMode, rng_step1, "rng_step1[i] free for all i in {1,3..7,9}")
    io_test(REPLMode, rng_step2, "rng_step2[i] free for all i in {-2,3,8}")
    io_test(REPLMode, rng_step3, "rng_step3[i] free for all i in {1}")
    io_test(REPLMode, rng_step4, "rng_step4[i] free for all i in {0,2}")
    io_test(REPLMode, arr_1, "arr_1[i] free for all i in {a,b,c}")
    io_test(REPLMode, arr_2, "arr_2[i] free for all i in {a,1,test}")
    io_test(REPLMode, arr_3, "arr_3[i] free for all i in {apple,banana,carrot,diamonds}")
    io_test(REPLMode, rng2_1, "rng2_1[i,j] free for all i in {1,2..9,10}, j in {a,b,c}")
    io_test(REPLMode, tri_1, "tri_1[i,j] free for all i in {1,2,3}, j in {..}")
    io_test(REPLMode, tri_2, "tri_2[i,j] free for all i in {1,2,3}, j in {..}")
    io_test(REPLMode, tri_3, "tri_3[(i,j),k] free for all (i,j) in {(1,3),(2,4)..(4,6),(5,7)}, k in {..}")

    io_test(IJuliaMode, rng_unit1, "rng_unit1_{i} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}")
    io_test(IJuliaMode, rng_unit2, "rng_unit2_{i} free \\quad\\forall i \\in \\{-2,-1,\\dots,2,3\\}")
    io_test(IJuliaMode, rng_unit3, "rng_unit3_{i} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}")
    io_test(IJuliaMode, rng_step1, "rng_step1_{i} free \\quad\\forall i \\in \\{1,3,\\dots,7,9\\}")
    io_test(IJuliaMode, rng_step2, "rng_step2_{i} free \\quad\\forall i \\in \\{-2,3,8\\}")
    io_test(IJuliaMode, rng_step3, "rng_step3_{i} free \\quad\\forall i \\in \\{1\\}")
    io_test(IJuliaMode, rng_step4, "rng_step4_{i} free \\quad\\forall i \\in \\{0,2\\}")
    io_test(IJuliaMode, arr_1, "arr_1_{i} free \\quad\\forall i \\in \\{a,b,c\\}")
    io_test(IJuliaMode, arr_2, "arr_2_{i} free \\quad\\forall i \\in \\{a,1,test\\}")
    io_test(IJuliaMode, arr_3, "arr_3_{i} free \\quad\\forall i \\in \\{apple,banana,carrot,diamonds\\}")
    io_test(IJuliaMode, rng2_1, "rng2_1_{i,j} free \\quad\\forall i \\in \\{1,2,\\dots,9,10\\}, j \\in \\{a,b,c\\}")
    io_test(IJuliaMode, tri_1, "tri_1_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{..\\}")
    io_test(IJuliaMode, tri_2, "tri_2_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{..\\}")
    io_test(IJuliaMode, tri_3, "tri_3_{(i,j),k} free \\quad\\forall (i,j) \\in \\{(1,3),(2,4),\\dots,(4,6),(5,7)\\}, k \\in \\{..\\}")
    end

    #------------------------------------------------------------------
    # Test category printing
    context("category printing") do
    @defVar(m, cat_bin[1:3], Bin)
    @defVar(m, 2 <= cat_int[1:3] <= 5, Int)
    @defVar(m, cat_semiint_both[2:3] >= 2, SemiInt)
    @defVar(m, i <= cat_semiint_difflow[i=2:3] <= 4, SemiInt)
    @defVar(m, 2 <= cat_semiint_diffup[i=2:3] <= i, SemiInt)
    @defVar(m, i <= cat_semiint_none[i=2:3] <= 2i, SemiInt)
    @defVar(m, cat_semicont_both[2:3] >= 2, SemiCont)
    @defVar(m, i <= cat_semicont_difflow[i=2:3] <= 4, SemiCont)
    @defVar(m, 2 <= cat_semicont_diffup[i=2:3] <= i, SemiCont)
    @defVar(m, i <= cat_semicont_none[i=2:3] <= 2i, SemiCont)
    @defVar(m, fixed_var[i=2:3] == i)

    io_test(REPLMode, cat_bin, "cat_bin[i] in {0,1} for all i in {1,2,3}")
    io_test(REPLMode, cat_int, "2 $le cat_int[i] $le 5, integer, for all i in {1,2,3}")
    io_test(REPLMode, cat_semiint_both, "cat_semiint_both[i] in {2..Inf} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_difflow, "cat_semiint_difflow[i] in {....4} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_diffup, "cat_semiint_diffup[i] in {2....} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semiint_none, "cat_semiint_none[i] in {......} or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_both, "cat_semicont_both[i] in [2,Inf] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_difflow, "cat_semicont_difflow[i] in [..,4] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_diffup, "cat_semicont_diffup[i] in [2,..] or {0} for all i in {2,3}")
    io_test(REPLMode, cat_semicont_none, "cat_semicont_none[i] in [..,..] or {0} for all i in {2,3}")
    io_test(REPLMode, fixed_var, "fixed_var[i] = .. for all i in {2,3}")

    io_test(IJuliaMode, cat_bin, "cat_bin_{i} \\in \\{0,1\\} \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_int, "2 \\leq cat_int_{i} \\leq 5, \\in \\mathbb{Z}, \\quad\\forall i \\in \\{1,2,3\\}")
    io_test(IJuliaMode, cat_semiint_both, "cat_semiint_both_{i} \\in \\{2,\\dots,\\intfy\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_difflow, "cat_semiint_difflow_{i} \\in \\{..,\\dots,4\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_diffup, "cat_semiint_diffup_{i} \\in \\{2,\\dots,..\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semiint_none, "cat_semiint_none_{i} \\in \\{..,\\dots,..\\} \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_both, "cat_semicont_both_{i} \\in \\[2,\\intfy\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_difflow, "cat_semicont_difflow_{i} \\in \\[..,4\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_diffup, "cat_semicont_diffup_{i} \\in \\[2,..\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, cat_semicont_none, "cat_semicont_none_{i} \\in \\[..,..\\] \\cup \\{0\\} \\quad\\forall i \\in \\{2,3\\}")
    io_test(IJuliaMode, fixed_var, "fixed_var_{i} = .. \\quad\\forall i \\in \\{2,3\\}")
    end

    #------------------------------------------------------------------
    # Tests for particular issues
    context("[print] Empty JuMPContainer printing (#124)") do
    @defVar(m, empty_free[1:0])
    io_test(REPLMode, empty_free, "Empty Array{Variable} (no indices)")
    io_test(IJuliaMode, empty_free, "Empty Array{Variable} (no indices)")
    @defVar(m, empty_set[[]])
    io_test(REPLMode, empty_set, "empty_set (no indices)")
    io_test(IJuliaMode, empty_set, "empty_set (no indices)")
    end
end



facts("[print] JuMPContainer{Number}") do
    # The same output for REPL and IJulia, so only testing one
    mod = Model()
    @defVar(mod, i*j <= w[i=9:10, [:Apple,5,:Banana], j=-1:+1] <= i*j)
    @defVar(mod, i*j*k <= x[i=9:11,j=99:101,k=3:4] <= i*j*k)
    @defVar(mod, i*j <= y[i=9:11,j=i:11] <= i*j)
    @defVar(mod, j <= z[i=[:a,'b'],j=1:3] <= j)
    solve(mod)

    # Deal with hashing variations
    if hash(5) < hash(:Apple)
        io_test(REPLMode, getValue(w), """
w: 3 dimensions, 18 entries:
 [ 9,     5,-1] = -9.0
 [ 9,     5, 0] = 0.0
 [ 9,     5, 1] = 9.0
 [ 9, Apple,-1] = -9.0
 [ 9, Apple, 0] = 0.0
 [ 9, Apple, 1] = 9.0
 [ 9,Banana,-1] = -9.0
 [ 9,Banana, 0] = 0.0
 [ 9,Banana, 1] = 9.0
 [10,     5,-1] = -10.0
 [10,     5, 0] = 0.0
 [10,     5, 1] = 10.0
 [10, Apple,-1] = -10.0
 [10, Apple, 0] = 0.0
 [10, Apple, 1] = 10.0
 [10,Banana,-1] = -10.0
 [10,Banana, 0] = 0.0
 [10,Banana, 1] = 10.0""", repl=:print)
    else
        io_test(REPLMode, getValue(w), """
w: 3 dimensions, 18 entries:
 [ 9, Apple,-1] = -9.0
 [ 9, Apple, 0] = 0.0
 [ 9, Apple, 1] = 9.0
 [ 9,Banana,-1] = -9.0
 [ 9,Banana, 0] = 0.0
 [ 9,Banana, 1] = 9.0
 [ 9,     5,-1] = -9.0
 [ 9,     5, 0] = 0.0
 [ 9,     5, 1] = 9.0
 [10, Apple,-1] = -10.0
 [10, Apple, 0] = 0.0
 [10, Apple, 1] = 10.0
 [10,Banana,-1] = -10.0
 [10,Banana, 0] = 0.0
 [10,Banana, 1] = 10.0
 [10,     5,-1] = -10.0
 [10,     5, 0] = 0.0
 [10,     5, 1] = 10.0""", repl=:print)
    end

    io_test(REPLMode, getValue(x), """
x: 3 dimensions:
[ 9,:,:]
  [ 9, 99,:]
    [ 9, 99,3] = 2673.0
    [ 9, 99,4] = 3564.0
  [ 9,100,:]
    [ 9,100,3] = 2700.0
    [ 9,100,4] = 3600.0
  [ 9,101,:]
    [ 9,101,3] = 2727.0
    [ 9,101,4] = 3636.0
[10,:,:]
  [10, 99,:]
    [10, 99,3] = 2970.0
    [10, 99,4] = 3960.0
  [10,100,:]
    [10,100,3] = 3000.0
    [10,100,4] = 4000.0
  [10,101,:]
    [10,101,3] = 3030.0
    [10,101,4] = 4040.0
[11,:,:]
  [11, 99,:]
    [11, 99,3] = 3267.0
    [11, 99,4] = 4356.0
  [11,100,:]
    [11,100,3] = 3300.0
    [11,100,4] = 4400.0
  [11,101,:]
    [11,101,3] = 3333.0
    [11,101,4] = 4444.0
""", repl=:print)

    io_test(REPLMode, getValue(y), """
y: 2 dimensions, 6 entries:
 [ 9, 9] = 81.0
 [ 9,10] = 90.0
 [ 9,11] = 99.0
 [10,10] = 100.0
 [10,11] = 110.0
 [11,11] = 121.0""")

    # Deal with hashing variations
    first_hash  = hash(:a) < hash('b') ? "a" : "b"
    second_hash = first_hash == "a" ? "b" : "a"
    io_test(REPLMode, getValue(z), """
z: 2 dimensions, 6 entries:
 [$first_hash,1] = 1.0
 [$first_hash,2] = 2.0
 [$first_hash,3] = 3.0
 [$second_hash,1] = 1.0
 [$second_hash,2] = 2.0
 [$second_hash,3] = 3.0""")

end



facts("[print] SOS constraints") do
    modS = Model()
    a = [1,2,3]
    @defVar(modS, x[1:3], Bin)
    addSOS1(modS, [a[i]x[i] for i in 1:3])
    s1 = JuMP.SOSConstraint([x[i] for i in 1:3],
                            [a[i] for i in 1:3], :SOS1)
    io_test(REPLMode, s1, "SOS1: {1 x[1], 2 x[2], 3 x[3]}")
    io_test(IJuliaMode, s1, "SOS1: \\{1 x[1], 2 x[2], 3 x[3]\\}")

    b = [5,4,7,2,1]
    @defVar(modS, y[1:5], Bin)
    s2 = JuMP.SOSConstraint([y[i] for i in 1:5],
                            [b[i] for i in 1:5], :SOS2)
    io_test(REPLMode, s2, "SOS2: {5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]}")
    io_test(IJuliaMode, s2, "SOS2: \\{5 y[1], 4 y[2], 7 y[3], 2 y[4], 1 y[5]\\}")
end



facts("[print] Model") do
    le, ge = JuMP.repl[:leq], JuMP.repl[:geq]

    #------------------------------------------------------------------

    mod_1 = Model()
    @defVar(mod_1, a>=1)
    @defVar(mod_1, b<=1)
    @defVar(mod_1, -1<=c<=1)
    @defVar(mod_1, a1>=1,Int)
    @defVar(mod_1, b1<=1,Int)
    @defVar(mod_1, -1<=c1<=1,Int)
    @defVar(mod_1, x, Bin)
    @defVar(mod_1, y)
    @defVar(mod_1, z, Int)
    @defVar(mod_1, sos[1:3], Bin)
    @defVar(mod_1, 2 <= si <= 3, SemiInt)
    @defVar(mod_1, 2 <= sc <= 3, SemiCont)
    @defVar(mod_1, fi == 9)
    @setObjective(mod_1, Max, a - b + 2a1 - 10x)
    @addConstraint(mod_1, a + b - 10c - 2x + c1 <= 1)
    @addConstraint(mod_1, a*b <= 2)
    addSOS1(mod_1, [i*sos[i] for i in 1:3])

    io_test(REPLMode, mod_1, """
Max a - b + 2 a1 - 10 x
Subject to
 a + b - 10 c - 2 x + c1 $le 1
 a*b - 2 $le 0
 SOS1: {1 sos[1], 2 sos[2], 3 sos[3]}
 sos[i] in {0,1} for all i in {1,2,3}
 a $ge 1
 b $le 1
 -1 $le c $le 1
 a1 $ge 1, integer
 b1 $le 1, integer
 -1 $le c1 $le 1, integer
 x in {0,1}
 y free
 z free, integer
 si in {2..3} or {0}
 sc in [2,3] or {0}
 fi = 9
""", repl=:print)

    io_test(IJuliaMode, mod_1, """
\\begin{alignat*}{1}\\max\\quad & a - b + 2 a1 - 10 x\\\\
\\text{Subject to} \\quad & a + b - 10 c - 2 x + c1 \\leq 1\\\\
 & a\\timesb - 2 \\leq 0\\\\
 & SOS1: \\{1 sos[1], 2 sos[2], 3 sos[3]\\}\\\\
 & sos_{i} \\in \\{0,1\\} \\quad\\forall i \\in \\{1,2,3\\}\\\\
 & a \\geq 1\\\\
 & b \\leq 1\\\\
 & -1 \\leq c \\leq 1\\\\
 & a1 \\geq 1, \\in \\mathbb{Z}\\\\
 & b1 \\leq 1, \\in \\mathbb{Z}\\\\
 & -1 \\leq c1 \\leq 1, \\in \\mathbb{Z}\\\\
 & x \\in \\{0,1\\}\\\\
 & y free\\\\
 & z free, \\in \\mathbb{Z}\\\\
 & si \\in \\{2,\\dots,3\\} \\cup \\{0\\}\\\\
 & sc \\in \\[2,3\\] \\cup \\{0\\}\\\\
 & fi = 9\\\\
\\end{alignat*}
""")

    #------------------------------------------------------------------

    mod_2 = Model()
    @defVar(mod_2, x, Bin)
    @defVar(mod_2, y, Int)
    @addConstraint(mod_2, x*y <= 1)

    io_test(REPLMode, mod_2, """
Feasibility problem with:
 * 0 linear constraints
 * 1 quadratic constraint
 * 2 variables: 1 binary, 1 integer
Solver set to Default""", repl=:show)

    mod_2 = Model()
    @defVar(mod_2, x)
    @addConstraint(mod_2, x <= 3)

    io_test(REPLMode, mod_2, """
Feasibility problem with:
 * 1 linear constraint
 * 1 variable
Solver set to Default""", repl=:show)

    #------------------------------------------------------------------

    mod_3 = Model()

    @defVar(mod_3, x[1:5])
    @addNLConstraint(mod_3, x[1]*x[2] == 1)
    @addNLConstraint(mod_3, x[3]*x[4] == 1)
    @addNLConstraint(mod_3, x[5]*x[1] == 1)
    @setNLObjective(mod_3, Min, x[1]*x[3])

    io_test(REPLMode, mod_3, """
Min (nonlinear expression)
Subject to
 3 nonlinear constraints
 x[i] free for all i in {1,2..4,5}
""", repl=:print)
    io_test(REPLMode, mod_3, """
Minimization problem with:
 * 0 linear constraints
 * 3 nonlinear constraints
 * 5 variables
Solver set to Default""", repl=:show)
    io_test(IJuliaMode, mod_3, """
\\begin{alignat*}{1}\\min\\quad & (nonlinear expression)\\\\
\\text{Subject to} \\quad & 3 nonlinear constraints\\\\
 & x_{i} free \\quad\\forall i \\in \\{1,2,\\dots,4,5\\}\\\\
\\end{alignat*}
""", repl=:print)
end

facts("[print] changing variable categories") do
    le, ge = JuMP.repl[:leq], JuMP.repl[:geq]
    mod = Model()
    @defVar(mod, x[1:3])
    @defVar(mod, y[i=1:3,i:3])
    setCategory(x[3], :SemiCont)
    setCategory(y[1,3], :Int)

    io_test(REPLMode, mod, """
Min 0
Subject to
 x[i] free for all i in {1,2,3}
 y[i,j] free for all i in {1,2,3}, j in {..}
 x[1] free
 x[2] free
 x[3] in [-Inf,Inf] or {0}
 y[1,1] free
 y[1,2] free
 y[1,3] free, integer
 y[2,2] free
 y[2,3] free
 y[3,3] free
""", repl=:print)

    io_test(IJuliaMode, mod, """
\\begin{alignat*}{1}\\min\\quad & 0\\\\
\\text{Subject to} \\quad & x_{i} free \\quad\\forall i \\in \\{1,2,3\\}\\\\
 & y_{i,j} free \\quad\\forall i \\in \\{1,2,3\\}, j \\in \\{..\\}\\\\
 & x_{1} free\\\\
 & x_{2} free\\\\
 & x_{3} \\in \\[-Inf,Inf\\] \\cup \\{0\\}\\\\
 & y_{1,1} free\\\\
 & y_{1,2} free\\\\
 & y_{1,3} free, \\in \\mathbb{Z}\\\\
 & y_{2,2} free\\\\
 & y_{2,3} free\\\\
 & y_{3,3} free\\\\
\\end{alignat*}
""")
end

facts("[print] expressions") do
    # Most of the expression logic is well covered by test/operator.jl
    # This is really just to check IJulia printing for expressions
    le, ge = JuMP.repl[:leq], JuMP.repl[:geq]

    #------------------------------------------------------------------
    mod = Model()
    @defVar(mod, x[1:5])
    @defVar(mod, y[i=2:4,j=i:5])
    @defVar(mod, z)

    @addConstraint(mod, x[1] + 2*y[2,3] <= 3)
    io_test(REPLMode, mod.linconstr[end], "x[1] + 2 y[2,3] $le 3")
    io_test(IJuliaMode, mod.linconstr[end], "x_{1} + 2 y_{2,3} \\leq 3")

    if VERSION > v"0.4.0-"
        @addConstraint(mod, (x[1]+x[2])*(y[2,2]+3.0) <= 1)
        io_test(REPLMode, mod.quadconstr[end], "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2] - 1 $le 0")
        io_test(IJuliaMode, mod.quadconstr[end], "x_{1}\\timesy_{2,2} + x_{2}\\timesy_{2,2} + 3 x_{1} + 3 x_{2} - 1 \\leq 0")
    end

    @addConstraint(mod, (y[2,2]+3.0)*(x[1]+x[2]) <= 1)
    io_test(REPLMode, mod.quadconstr[end], "x[1]*y[2,2] + x[2]*y[2,2] + 3 x[1] + 3 x[2] - 1 $le 0")
    io_test(IJuliaMode, mod.quadconstr[end], "x_{1}\\timesy_{2,2} + x_{2}\\timesy_{2,2} + 3 x_{1} + 3 x_{2} - 1 \\leq 0")
end



facts("[print] Variable") do
    m = Model()
    @defVar(m, 0 <= x <= 2, inconstraints=ConstraintRef{LinearConstraint}[], objective=0.0, coefficients=Float64[] )

    @fact    getName(x) --> "x"
    io_test(REPLMode,   x, "x")
    io_test(IJuliaMode, x, "x")

    setName(x, "x2")
    @fact    getName(x) --> "x2"
    io_test(REPLMode,   x, "x2")
    io_test(IJuliaMode, x, "x2")

    setName(x, "")
    @fact    getName(x) --> "col_1"
    io_test(REPLMode,   x, "col_1")
    io_test(IJuliaMode, x, "col_1")

    @defVar(m, z[1:2,3:5])
    @fact       getName(z[1,3]) --> "z[1,3]"
    io_test(REPLMode,   z[1,3],    "z[1,3]")
    io_test(IJuliaMode, z[1,3],    "z_{1,3}")
    @fact       getName(z[2,4]) --> "z[2,4]"
    io_test(REPLMode,   z[2,4],    "z[2,4]")
    io_test(IJuliaMode, z[2,4],    "z_{2,4}")
    @fact       getName(z[2,5]) --> "z[2,5]"
    io_test(REPLMode,   z[2,5],    "z[2,5]")
    io_test(IJuliaMode, z[2,5],    "z_{2,5}")

    @defVar(m, w[3:9,["red","blue","green"]])
    @fact    getName(w[7,"green"]) --> "w[7,green]"
    io_test(REPLMode,   w[7,"green"], "w[7,green]")
    io_test(IJuliaMode, w[7,"green"], "w_{7,green}")

    rng = 2:5
    @defVar(m, v[rng,rng,rng,rng,rng,rng,rng])
    a_v = v[4,5,2,3,2,2,4]
    @fact    getName(a_v) --> "v[4,5,2,3,2,2,4]"
    io_test(REPLMode,   a_v, "v[4,5,2,3,2,2,4]")
    io_test(IJuliaMode, a_v, "v_{4,5,2,3,2,2,4}")
end
