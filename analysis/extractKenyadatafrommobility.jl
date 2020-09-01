using CSV,DataFrames,Plots,Dates,Statistics,JLD2,Polynomials
using Plots.PlotMeasures
pyplot()
## Extract Kenya data from the world data
mobilitydata = DataFrame!(CSV.File("Global_Mobility_Report_25082020.csv"))
mobilitydata_kenya = mobilitydata[mobilitydata.country_region .== "Kenya",:]
CSV.write("Kenya_Mobility_Report_25082020.csv",mobilitydata_kenya)
Kenyawide = mobilitydata_kenya[(ismissing.(mobilitydata_kenya.sub_region_1)).&(ismissing.(mobilitydata_kenya.metro_area)),:]
nairobidata = mobilitydata_kenya[(mobilitydata_kenya.sub_region_1 .=="Nairobi County").&(.~ismissing.(mobilitydata_kenya.sub_region_1)),:]
mombasadata = mobilitydata_kenya[(mobilitydata_kenya.sub_region_1 .=="Mombasa County").&(.~ismissing.(mobilitydata_kenya.sub_region_1)),:]
# ugishdata = mobilitydata_kenya[(mobilitydata_kenya.sub_region_1 .=="Uasin Gishu County").&(.~ismissing.(mobilitydata_kenya.sub_region_1)),:]


function sma(a::Array, n::Int)
    vals = zeros(size(a,1) - (n-1), size(a,2))

    for i in 1:size(a,1) - (n-1)
        for j in 1:size(a,2)
            vals[i,j] = mean(a[i:i+(n-1),j])
        end
    end

    vals
end
plot(mean(Matrix(Kenyawide[:,[9,10,12,13]]),dims = 2))
plot!(sma(Matrix(Kenyawide[:,9:end]),7))


## Plot some mobility trends

function riskycontacts_sma(df,n)
    areaname = ""
    if !ismissing(df.sub_region_1[1])
        areaname = df.sub_region_1[1]
    else
        areaname = "Kenyawide"
    end
    return (contactrate = 1 .+ sma(mean(Matrix(df[:,[9,10,12,13]]),dims = 2),n)./100, date = df.date[n:end], area = areaname)
end

function plot_mobilitytrends(df,title_str,n)
    plt = plot(df.date[n:end],sma(Matrix(df[:,9:end]),n),
        lab = ["Retail and recreation" "Grocery and pharmacy" "Parks" "Transit" "Workplace" "Residential"],
        ylims = (-60,80),
        lw = 1,
        ylabel = "% change from baseline",
        title = title_str*", $(n)-day moving av.",
        size = (700,500),dpi = 250)
    plot!(plt,df.date[n:end],sma(mean(Matrix(df[:,[9,10,12,13]]),dims = 2),n),
        lab = "Risky contacts",
        ylims = (-60,80),
        lw = 2,
        color = :black,
        ls = :dash,
        ylabel = "% change from baseline",
        title = title_str*", $(n)-day moving av.",
        size = (700,500),dpi = 250)
    return plt
end

plt_kenyawide = plot_mobilitytrends(Kenyawide,"Kenya-wide mobility trend (Google)",7)
plt_nairobi = plot_mobilitytrends(nairobidata,"Nairobi mobility trend (Google)",7)
plt_mombasa = plot_mobilitytrends(mombasadata,"Mombasa mobility trend (Google)",7)
# plt_ugish = plot_mobilitytrends(ugishdata,"Uasin Gishu mobility trend (Google)",7)

savefig(plt_kenyawide,"Kenyawidemobility.png")
savefig(plt_nairobi,"Nairobimobility.png")
savefig(plt_mombasa,"Mombasamobility.png")


contactrate_kenya = riskycontacts_sma(Kenyawide,7)
contactrate_nairobi = riskycontacts_sma(nairobidata,7)
contactrate_mombasa = riskycontacts_sma(mombasadata,7)

@save("contact_data_25082020.jld2",contactrate_kenya,contactrate_nairobi,contactrate_mombasa)


function projectcontactrate(contactrate_data,fittingdate::Date,projectiondate::Date)
    f = findfirst(contactrate_nairobi.date .== fittingdate)
    ys = contactrate_data.contactrate[f:end]
    xs = collect(1:length(ys)) .-1
    fit1 = fit(xs,ys,1)
    p = Polynomial(fit1.coeffs)
    _xs = collect((length(ys)+1):((projectiondate - contactrate_nairobi.date[end]).value))
    _ys = [min(1,p(x)) for x in _xs]
    _date = [contactrate_nairobi.date[end] + Day(t) for (t,d) in enumerate(_xs)]
     # vcat(contactrate_data.contactrate,_ys)
     #
     #  vcat(contactrate_data.date,_date)
    return (contactrate = vcat(contactrate_data.contactrate,_ys) , date = vcat(contactrate_data.date,_date), area = contactrate_data.area)
end

projected_contactrate_nairobi = projectcontactrate(contactrate_nairobi,Date(2020,7,1),Date(2021,3,1))
projected_contactrate_mombasa = projectcontactrate(contactrate_mombasa,Date(2020,7,1),Date(2021,3,1))
projected_contactrate_kenya = projectcontactrate(contactrate_kenya,Date(2020,7,1),Date(2021,3,1))

@save("projected_contact_data_25082020.jld2",projected_contactrate_nairobi,projected_contactrate_mombasa,projected_contactrate_kenya)

plot!(contactrate_nairobi.date,contactrate_nairobi.contactrate,lw=2)
plt_mob = plot(projected_contactrate_nairobi.date,projected_contactrate_nairobi.contactrate,
    ls = :dash,
    lab = "",
    color = :red,
    ylims = (0.5,1.3),
    size = (700,500),dpi = 250,
    right_margin = 1cm,
    title = "Contact rates: Google mobility trends and projections",
    ylabel = "Mobility compared to baseline")
plot!(plt_mob,projected_contactrate_mombasa.date,projected_contactrate_mombasa.contactrate,
    ls = :dash,
    lab = "",
    color = :green)
plot!(plt_mob,projected_contactrate_kenya.date,projected_contactrate_kenya.contactrate,
    ls = :dash,
    lab = "",
    color = :blue)

plot!(plt_mob,contactrate_nairobi.date,contactrate_nairobi.contactrate,
    lw = 2,
    lab = "Nairobi",
    color = :red)
plot!(plt_mob,contactrate_mombasa.date,contactrate_mombasa.contactrate,
    lw = 2,
    lab = "Mombasa",
    color = :green)
plot!(plt_mob,contactrate_kenya.date,contactrate_kenya.contactrate,
    lw = 2,
    lab = "Kenya (overall)",
    color = :blue)

savefig(plt_mob,"Kenyawidemobility_with_projections.png")
