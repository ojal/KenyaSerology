## Load PATH and relevant packages (most of these aren't required, but I run everything in Main module
# then put into KenyaSerology module)

push!(LOAD_PATH, joinpath(homedir(),"GitHub/KenyaSerologyPrivate/src"))
using Distributions,Plots,Dates,JLD2,TransformVariables,Optim,FileIO,DataFrames,CSV
using Revise
import KenyaSerology
## Choose plotting backend. I currently like pyplot (which uses the matplotlib python backend) because
#it renders Greek symbols
# pyplot()
# plotlyjs()
gr()


## Load data and analysis pipeline
@load("data/case_data_with_pos_neg_21feb_to_30sept.jld2")
@load("data/serologydata_21feb_6aug.jld2")

#Zero fill serodata to be same length as case data
zerofill = zeros(Int64,size(case_data_with_pos_neg.cases,1)-size(sero_data.serodata,1),30,2)
zerofilledserodata = vcat(sero_data.serodata,zerofill)
filleddates = vcat(sero_data.dates,[sero_data.dates[end] + Day(k+1) for k = 1:size(zerofill,1)])

sero_data = (serodata =zerofilledserodata, areas = sero_data.areas,dates = filleddates)

@load("data/relative_testing_rate01102020.jld2")
@load("data/projected_contact_data_09102020.jld2")
@load("data/death_data_by_area_21feb_to_17oct.jld2")
@load("data/p_ID.jld2")
df_kenya = DataFrame!(CSV.File("data/2019-population_census-report-per-county.csv"))
df_kenya.County
converted_countynames = copy(df_kenya.County)
converted_countynames[converted_countynames.=="Elgeyo-Marakwet"] .= "Elgeyo Marakwet"
converted_countynames[converted_countynames.=="Murang'a"] .= "Muranga"
converted_countynames[converted_countynames.=="Tharaka-Nithi"] .= "Tharaka Nithi"
converted_countynames[converted_countynames.=="Homa Bay"] .= "HomaBay"

# scatter(case_data_with_pos_neg.dates,case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .=="NAIROBI",1])
# scatter!(case_data_with_pos_neg.dates,case_data_with_pos_neg.cases[:,case_data_with_pos_neg.areas .=="NAIROBI",2])
# scatter(death_data.dates,death_data.deaths[:,death_data.areas .== "NAIROBI"],lab="",
#         color = :black, title = "Nairobi deaths")
#
# scatter(death_data.dates,sum(death_data.deaths,dims =2 ),lab="",
#         color = :black, title = "Kenyan deaths")

## Put in the variable transformation
trans_Peff = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))
trans_neg_PCR = as((R = asℝ₊,E₀ = asℝ₊,I₀ = asℝ₊,χ = asℝ₊,M_PCR = asℝ₊,α = asℝ₊,p_test = as𝕀,P_eff = as𝕀))




## Alphabetic list of counties
countynames = sort(case_data_with_pos_neg.areas)
nairobiname = ["Nairobi"]
mombasaname = ["Mombasa"]
neighbouringnairobi = ["Kajiado","Kiambu","Machakos"]
semiurbancounties = ["Kajiado", "Kiambu", "Kilifi", "Kisumu", "Machakos", "Uasin Gishu", "Vihiga", "Isiolo" , "Nakuru" ]
ruralcounties = ["Bungoma", "Busia", "Garissa", "HomaBay", "Kericho", "Kisii", "Kitui", "Kwale", "Lamu", "Makueni", "Marsabit",
                "Migori", "Narok","Nyamira", "Nyeri", "Siaya", "Taita Taveta", "Tana River", "Trans Nzoia", "Turkana", "Wajir",
                "Baringo", "Bomet","Elgeyo Marakwet", "Embu", "Kakamega", "Kirinyaga", "Laikipia", "Mandera", "Meru", "Muranga",
                "Nandi","Nyandarua","Samburu","Tharaka Nithi", "West Pokot"]

county_using_coast_linelist = ["Mombasa","Kilifi","Kwale","Taita Taveta", "Tana River","Lamu"]


length(vcat(nairobiname,mombasaname,semiurbancounties,ruralcounties))
countynames = vcat(nairobiname,mombasaname,semiurbancounties,ruralcounties)
# countynames = setdiff(vcat(nairobiname,mombasaname,semiurbancounties,ruralcounties),county_using_coast_linelist)


countyprior_pairs = [(nairobiname,KenyaSerology.nairobi_denominator_prior),
                        (mombasaname,KenyaSerology.nairobi_denominator_prior),#Now assuming same prior for Mombasa and Nairobi
                        (semiurbancounties,KenyaSerology.semiurban_denominator_prior),
                        (ruralcounties,KenyaSerology.rural_denominator_prior)]

countycontactrate_pairs = [(nairobiname,projected_contactrate_nairobi),
                        (mombasaname,projected_contactrate_mombasa),
                        (semiurbancounties,projected_contactrate_kenya),
                        (ruralcounties,projected_contactrate_kenya)]

countytestingrate_pairs = [(nairobiname,relative_testing_rate),
                        (neighbouringnairobi,relative_testing_rate),
                        (mombasaname,ones(length(relative_testing_rate))),
                        (setdiff(semiurbancounties,neighbouringnairobi),ones(length(relative_testing_rate))),
                        (ruralcounties,ones(length(relative_testing_rate)))]

countypopulation_pairs = [converted_countynames df_kenya.Total_Population19]


function findpairing(county,countypairs)
        countynamegroups = [pairs[1] for pairs in countypairs]
        return findfirst(in.(county,countynamegroups))
end
# findpairing("Garissa",countytestingrate_pairs)
function get_target_data(county,countypairs)
        f = findpairing(county,countypairs)
        return countypairs[f][2]
end

# plot(get_target_data("Nairobi",countytestingrate_pairs))
# plot!(get_target_data("Garissa",countytestingrate_pairs))
#

## Create function to loop over fitting counties
function loopoverfitting(listofcounties,countyprior_pairs,countycontactrate_pairs,countytestingrate_pairs,countypopulationsize_pairs,case_data,sero_data;
                        enddate,n,parameter_transformation,loglikelihood_func,savedir)
        countieswherefitfailed = []
        for areaname in listofcounties
                popsize = countypopulationsize_pairs[countypopulationsize_pairs[:,1].==areaname,2][1]
                logprior = get_target_data(areaname,countyprior_pairs)
                contactrate_data = get_target_data(areaname,countycontactrate_pairs)
                relative_testing_rate = get_target_data(areaname,countytestingrate_pairs)
                try
                        model = KenyaSerology.pipeline_for_fit(areaname,popsize,logprior,case_data,sero_data,contactrate_data,relative_testing_rate,enddate;
                                                n=n,parameter_transformation=parameter_transformation,loglikelihood_func = loglikelihood_func)
                        @save(joinpath(savedir,areaname*"_model.jld2"),model)
                catch
                        push!(countieswherefitfailed,areaname)
                end
        end
        return countieswherefitfailed
end


countieswherefitfailed = loopoverfitting(["Nakuru"],countyprior_pairs,countycontactrate_pairs,countytestingrate_pairs,countypopulation_pairs,case_data_with_pos_neg,sero_data,
                enddate = Date(2020,10,1),n=10000,parameter_transformation=trans_neg_PCR,loglikelihood_func = KenyaSerology.loglikelihood_with_negPCR,savedir = "modelfits_pos_neg/")
findfirst(countynames .== "Nairobi")
