using Printf
"""
    convergence_rate(f, h, r=3)

Construct a LaTex table containing convergence data pertaining to function f.

# Arguments
- `f::Function`: a function depending on mesh size h which returns a numerical error value.
- `h::Real`: an initial grid mesh size.
- `r::Integer`: the number of mesh refinements to perform.

# Examples
```jldoctest
julia> h = 0.25
julia> f(x) = x
julia> r = 4
julia> convergence_rate(f, h, r)
{\\tabulinesep=1.5mm
\\begin{tabu}{|| c | c | c ||} 
\\hline \\hline
      \$h\$ & \$ \\left \\lVert u-u_h \\right \\rVert \$ & \$ \\log_2 \\big (e_h /e_{\\frac{h}{2}} \\big ) \$ \\\\  \\hline
  0.25000 &            0.250000 &       \$ \\emptyset \$ \\\\ \\hline
  0.12500 &            0.125000 &             1.00000 \\\\ \\hline
  0.06250 &            0.062500 &             1.00000 \\\\ \\hline
  0.03125 &            0.031250 &             1.00000 \\\\ \\hline
\\hline
\\end{tabu}}

```
"""
function convergence_rates(f, h, num_refinements=3)

    error = zeros(num_refinements)
    rate = zeros(size(error))
    @printf("{\\tabulinesep=1.5mm \n")
    @printf("\\begin{tabu}{|| c | c | c ||} \n")
    @printf("\\hline \\hline \n")
    @printf("%9s & %19s & %19s \\\\  \\hline \n",
            "\$h\$",
            "\$ \\left \\lVert u-u_h \\right \\rVert \$",
            "\$ \\log_2 \\big (e_h /e_{\\frac{h}{2}} \\big ) \$")

    for i = 1:num_refinements
        if i > 1
            error[i] = f(h)
            rate[i] = log2(error[i-1]/error[i])

            @printf("%9.5f & %19.f & %19.5f \\\\ \\hline  \n", h, error[i], rate[i])
        elseif i ==1
            error[i] = f(h)
            null = "\$ \\emptyset \$"

            @printf("%9.5f & %19.f & %19s \\\\ \\hline  \n", h, error[i], null)
        end


        h = h/2


    end
    @printf("\\hline \n")
    @printf("\\end{tabu}}")


end
