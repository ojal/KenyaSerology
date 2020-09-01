using CSV,DataFrames,Plots,Dates,Statistics,JLD2,Polynomials
using Plots.PlotMeasures
plotlyjs()
pyplot()
kenya_owid = DataFrame!(CSV.File("kenya_owid.csv"))

scatter(Date.(kenya_owid.date,DateFormat("dd/mm/yyyy")),kenya_owid.new_tests,
        lab = "",
        ylabel = "Number of new daily swab tests",
        title = "OWID testing data for Kenya",
        guidefont = 20,
        tickfont = 10,
        titlefont = 18,
        legendfont = 10,
        size = (700,500),dpi=250)

savefig("Kenya_testing_data.csv")

scatter(Date.(kenya_owid.date,DateFormat("dd/mm/yyyy")),kenya_owid.new_cases./kenya_owid.new_tests,
        lab = "",
        ylabel = "Number of new daily swab tests",
        title = "OWID testing data for Kenya")


ys = kenya_owid.new_tests[24:end]
xs = collect(0:127)[.~ismissing.(ys)]
_ys = Float64.(ys[.~ismissing.(ys)])

p = Polynomial(fit(xs,_ys./mean(skipmissing(ys)),1).coeffs)

pyplot()
dates = Date.(kenya_owid.date,DateFormat("dd/mm/yyyy"))
scatter(dates,kenya_owid.new_tests./mean(skipmissing(ys)),
        lab = "",
        ylabel = "Daily swab testing rate relative to mean",
        title = "Relative testing rate (Kenya-wide)",
        guidefont = 14,
        tickfont = 10,
        titlefont = 18,
        legendfont = 10,
        size = (700,500),dpi=250)
plot!([dates[24] + Day(Int64(x)) for x in xs],[p(x) for x in xs],
        lab = "Linear fit: relative testing rate increasing by $(round(p.coeffs[2],digits = 3)*100)% per day")

savefig("Kenya_testing_data_with_linear_fit.pdf")

date = collect(Date(2020,2,21):Day(1):Date(2021,1,1))
relative_testing_rate = [min(max(p((d - dates[24]).value),0.),4.) for d in date]
plot(date,relative_testing_rate)
relative_testing_rate_nairobi = (date = date, relative_testing_rate=relative_testing_rate)
@save("relative_testing_rate.jld2",relative_testing_rate_nairobi)
