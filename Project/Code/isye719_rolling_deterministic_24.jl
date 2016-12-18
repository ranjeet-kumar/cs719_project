using JuMP
using PyPlot
using Gurobi

include("helper_functions.jl")
#############################################################################################################
# Read load and price data from csv files
loadExceldata = readcsv("DATARCoast.csv");
rtmpriceExceldata = readcsv("AggregatedData_RTM_ALTA31GT_7_B1.csv")
dampriceExceldata = readcsv("AggregatedData_DAM_ALTA31GT_7_B1.csv")

generate_new_scenarios = 0; # 0: don't generate new scenarios for load profiles; 1: generate new scenarios
generate_new_sample_paths = 0; # 0: don't generate new sample paths for loads; 1: generate new sample paths
participate_rtm = 1; # 0: Don't participate in RTM; 1: participate in RTM
participate_dam = 1; # 0: Don't participate in DAM; 1: participate in DAM
makeplots = 0; # 0: Don't generate plots, 1: Generate plots

# Defining time parameters for the model and planning schedule
dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day-ahead market (hourly) intervals in a day
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtdam/dtrtm);		#No. of real time intervals in each hour
nrtmpoints = ndam*nrtm;		    #Total number of points in RTM data in one day
ndampoints = ndam;			  #Total number of points in hourly data in one day
weekly_ndays = 7;         #Number of days in every week
ndays = 365;              #Number of days data is available for
ndays_planning = 7;       #Number of days you want to plan for
nweeks_planning = Int64(ceil((ndays_planning/weekly_ndays)));
nhours_planning = ndays_planning*ndam;  #Number of hours we will plan the policy for
nrtm_planning = nhours_planning*nrtm;   #Number of rtm intervals we will plan the policy for
ndays_horizon = 1; # Number of days in horizon at every step of receding horizon scheme
nhours_horizon = ndays_horizon*ndam; # Number of hours in horizon at every step of receding horizon scheme
nrtm_horizon = nhours_horizon*nrtm; # Number of real-time intervals in horizon at every step of receding horizon scheme

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
ebat0 = ebat_max;		   #Initial State of charge
rampmin = -0.5*P_max;	          #Lower bound for ramp discharge, MW/5min
rampmax = 0.5*P_max;  	  #Upper bound for ramp discharge, MW/5min
NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;

# Price data
eprrtm = rtmpriceExceldata[1:nrtmpoints*ndays,4];	 #Real Time Market price, $/MWh
eprdam = dampriceExceldata[1:ndampoints*ndays,4];	 #Day Ahead Market price, $/MWh

# Generate new scenarios for loads
loadNSdatafile = "loads_scenarios_month.csv";
if generate_new_scenarios == 1
  load = Matrix{Float64}(loadExceldata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW
  loadNSdata = generate_weekly_loads_scenarios(load,1:52,ndays_planning,NS,loadNSdatafile);
end
loadNSdata = readcsv(loadNSdatafile)
ndays_data = (nweeks_planning+1)*weekly_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

#Reshape the data to matrices to be used in the model
rtmepr = reshape(eprrtm,nrtm,ndam,ndays);
damepr = reshape(eprdam,ndam,ndays);
load = loadNSplanningdata[:,:,1:ndays_planning,:]/1000; #MW

#Define sets to be used in the JuMP model
rtm = 1:nrtm; # {1,2,...,12}
dam = 1:ndam; # {1,2,...,24}
day = 1:ndays_planning; # {1,2,...,7}

if generate_new_sample_paths == 1
# Generate NS sample paths for realizations of loads in 7 days at hourly intervals
  (paths,loadperm) = generate_sample_paths(load,NS,"samplepaths.csv","sampleloadperm.csv");
end
# Take the NS sample paths for loads generated earlier
paths = Matrix{Int64}(readcsv("samplepaths.csv"));
loadperm = zeros(nrtm,ndam,ndays_planning,NS);
for s in S
  j = 1;
  for l in day
    for k in dam
        loadperm[:,k,l,s] = load[:,k,l,paths[j,s]];  #MW
        j = j+1;
    end
  end
end


load1 = loadNSplanningdata/1000; #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);

load1_mv = mean(load1,4);
loaddata1_mv = reshape(load1_mv,nrtm*ndam*ndays_data);

######################################################



dailyprofit = zeros(ndays);
loadplot = zeros(ndays*nrtmpoints);
rtmeprplot = zeros(ndays*nrtmpoints);
dameprplot = zeros(ndays*ndampoints);
ebatplot = zeros(ndays*nrtmpoints);
Prtmplot = zeros(ndays*nrtmpoints);
Pdamplot = zeros(ndays*ndampoints);
rtmeprofitplot = zeros(ndays*nrtmpoints);
dameprofitplot = zeros(ndays*ndampoints);

totalpower = zeros(nrtm,ndam,ndays);
totalpowerplot = zeros(nrtm,ndam,ndays);




Pdam_mv_rol = zeros(nhours_planning);

Prtm_realized = zeros(nrtm,nhours_planning);
unmetload_realized = zeros(nrtm,nhours_planning)
unmetcost_realized = zeros(nrtm,nhours_planning);
profitErtm_realized = zeros(nrtm,nhours_planning);
profitEdam_realized = zeros(nhours_planning);
profitE_realized = zeros(nhours_planning);
profittotal_realized = zeros(nhours_planning);
netobjective_realized = zeros(nhours_planning);
obj_dt_rh_NS = Vector()

mv_rol = nothing; m_rold = nothing;

for k in S[1] # Loop to evaluate cost along each scenario

realized_sequence = paths[:,k];

ebat0_mv = ebat_max;		  #Initial State of charge for mean-value problem, 100 means fully charged
ebat0 = ebat_max;               #Initial State of charge for all scenarios, 100 means fully charged

realized_sequence = Vector{Int64}(k*ones(nhours_planning));

Prtm_realized = zeros(nrtm,nhours_planning);
unmetload_realized = zeros(nrtm,nhours_planning)
unmetcost_realized = zeros(nrtm,nhours_planning);
profitErtm_realized = zeros(nrtm,nhours_planning);
profitEdam_realized = zeros(nhours_planning);
profitE_realized = zeros(nhours_planning);
profittotal_realized = zeros(nhours_planning);
netobjective_realized = zeros(nhours_planning);


j=1;
tic()
for p in 1:nhours_planning # Starting rolling horizon for mean-value problem


    #Load and price data
    load_mv = loaddata1_mv[(p-1)*nrtm+(1:nrtm_horizon)];	#Load, MW

    eprrtm = rtmpricedata[(p-1)*nrtm+(1:nrtm_horizon),4];	    	#Real Time Market price, $/MWh
    eprdam = dampricedata[(p-1)+(1:nhours_horizon),4];	    	#Day Ahead Market Selling price, $/MWh

    #Reshape the data to matrices
    rtmepr = reshape(eprrtm,nrtm,nhours_horizon);
    damepr = reshape(eprdam,nhours_horizon);
    load_mv = reshape(load_mv,nrtm,nhours_horizon);

    #Define sets to be used in the model defVar and addConstraint
    rtm = 1:nrtm;
    dam = 1:nhours_horizon;

    ################ Mean-value Model ##################

    mv_rol = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
        @variable(mv_rol, -P_max <= Prtm[rtm,dam] <= P_max)	                #Power sold to the real time market, kW
        @variable(mv_rol, -P_max <= Pdam[dam] <= P_max)    	                #Power sold to the day ahead market, kW
        @variable(mv_rol, 0 <= ebat[rtm,dam] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
        @variable(mv_rol, suppliedload[rtm,dam] >= 0)
        @variable(mv_rol, unmetload[rtm,dam] >= 0)
    		@expression(mv_rol, Pnet[i in rtm,k in dam], Prtm[i,k] + Pdam[k] + suppliedload[i,k])    #Net power discharged from battery in all 5-min interval, kW
    		@variable(mv_rol, profitErtm[rtm,dam])				        #Profit from the real time market, USD
        @variable(mv_rol, profitEdam[dam])	        			#Profit from the day ahead market, USD
        @variable(mv_rol, profittotal)		        	#Total profit in the day, USD
        @variable(mv_rol, unmetcost[rtm,dam])

        @constraint(mv_rol, InitialEnergy, ebat[1,1] == ebat0 - Pnet[1,1]*dtrtm)	#Inital energy in the battery
        @constraint(mv_rol, rtmEBalance[i in rtm[2:end],k in dam], ebat[i,k] == ebat[i-1,k] - Pnet[i,k]*dtrtm)	#Dynamics constraint
        @constraint(mv_rol, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end]], ebat[i,k] == ebat[iend,k-1] - Pnet[i,k]*dtrtm)	#Dynamics constraint
        @constraint(mv_rol, UnmetLoad[i in rtm,k in dam], suppliedload[i,k] + unmetload[i,k] >=  load_mv[i,k])
        @constraint(mv_rol, BoundSupplied[i in rtm,k in dam], suppliedload[i,k] <= load_mv[i,k])
        @constraint(mv_rol, BoundUnmet[i in rtm,k in dam], unmetload[i,k] <= load_mv[i,k])
    #=    @constraint(mv_rol, RTMRamp1[i in rtm[2:end],k in dam], Pnet[i,k,s]  - Pnet[i-1,k] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
    		@constraint(mv_rol, RTMRamp2[i in rtm[2:end],k in dam], Pnet[i,k,s]  - Pnet[i-1,k] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
    		@constraint(mv_rol, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end]], Pnet[i,k] - Pnet[iend,k-1] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
    		@constraint(mv_rol, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end]], Pnet[i,k] - Pnet[iend,k-1] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
    =#
    		@constraint(mv_rol, RTMEProfits[i in rtm,k in dam], profitErtm[i,k] == rtmepr[i,k]*Prtm[i,k]*dtrtm)	#Economic calculation
        @constraint(mv_rol, DAMEProfits[k in dam], profitEdam[k] == damepr[k]*Pdam[k]*dtdam)	#Economic calculation
        @constraint(mv_rol, TotalProfit, profittotal ==
                            sum{profitErtm[i,k], i in rtm, k in dam} + sum{profitEdam[k], k in dam})
        @constraint(mv_rol, UnmetCost[i in rtm, k in dam], unmetcost[i,k] == rtmepr[i,k]*unmetload[i,k])
        @constraint(mv_rol, NetDischarge1[i in rtm,k in dam], Pnet[i,k] <= P_max)
        @constraint(mv_rol, NetDischarge2[i in rtm,k in dam], Pnet[i,k] >= -P_max)
        @objective(mv_rol, Min, -profittotal + sum{unmetcost[i,k],i in rtm, k in dam})
        status = solve(mv_rol)

##########################################################################

    ebat0_mv = getvalue(getvariable(mv_rol,:ebat))[rtm[end],dam[1]];

    # Store first-stage sollution to be implemented at current step
    Pdam_mv_rol[p] = getvalue(getvariable(mv_rol,:Pdam))[1];


   ################ Resolving Stochastic Model For all scenarios to get second-stage ##################

    #Load data for all scenarios
    load = loaddata1[(p-1)*nrtm+(1:nrtm_horizon),:];	#Load, MW
    load = reshape(load,nrtm,nhours_horizon,NS);

    m_rold = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m_rold, -P_max <= Prtm[rtm,dam,S] <= P_max)	                #Power sold to the real time market, kW
    @variable(m_rold, -P_max <= Pdam[dam,S] <= P_max)    	                #Power sold to the day ahead market, kW
    @variable(m_rold, 0 <= ebat[rtm,dam,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
    @variable(m_rold, suppliedload[rtm,dam,S] >= 0)
    @variable(m_rold, unmetload[rtm,dam,S] >= 0)
		@expression(m_rold, Pnet[i in rtm,k in dam,s in S], Prtm[i,k,s] + Pdam[k,s] + suppliedload[i,k,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m_rold, profitErtm[rtm,dam,S])				        #Profit from the real time market, USD
    @variable(m_rold, profitEdam[dam,S])	        			#Profit from the day ahead market, USD
    @variable(m_rold, profittotal[S])		        	#Total profit in the day, USD
    @variable(m_rold, unmetcost[rtm,dam,S])

    @constraint(m_rold, InitialEnergy[s in S], ebat[1,1,s] == ebat0 - 1/eff*Pnet[1,1,s]*dtrtm)	#Inital energy in the battery
    @constraint(m_rold, rtmEBalance[i in rtm[2:end],k in dam,s in S], ebat[i,k,s] == ebat[i-1,k,s] - 1/eff*Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rold, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],s in S], ebat[i,k,s] == ebat[iend,k-1,s] - 1/eff*Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rold, UnmetLoad[i in rtm,k in dam,s in S], suppliedload[i,k,s] + unmetload[i,k,s] >=  load[i,k,s])
    @constraint(m_rold, BoundSupplied[i in rtm,k in dam,s in S], suppliedload[i,k,s] <= load[i,k,s])
    @constraint(m_rold, BoundUnmet[i in rtm,k in dam,s in S], unmetload[i,k,s] <= load[i,k,s])
#=    @constraint(m_rold, RTMRamp1[i in rtm[2:end],k in dam,s in S], Pnet[i,k,s]  - Pnet[i-1,k,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rold, RTMRamp2[i in rtm[2:end],k in dam,s in S], Pnet[i,k,s]  - Pnet[i-1,k,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rold, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],s in S], Pnet[i,k,s] - Pnet[iend,k-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rold, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],s in S], Pnet[i,k,s] - Pnet[iend,k-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
=#
		@constraint(m_rold, RTMEProfits[i in rtm,k in dam,s in S], profitErtm[i,k,s] == rtmepr[i,k]*Prtm[i,k,s]*dtrtm)	#Economic calculation
    @constraint(m_rold, DAMEProfits[k in dam,s in S], profitEdam[k,s] == damepr[k]*Pdam[k,s]*dtdam)	#Economic calculation
    @constraint(m_rold, TotalProfit[s in S], profittotal[s] ==
                        sum{profitErtm[i,k,s], i in rtm, k in dam} + sum{profitEdam[k,s], k in dam})
    @constraint(m_rold, UnmetCost[i in rtm, k in dam, s in S], unmetcost[i,k,s] == rtmepr[i,k]*unmetload[i,k,s])
    @constraint(m_rold, NetDischarge1[i in rtm,k in dam,s in S], Pnet[i,k,s] <= P_max)
    @constraint(m_rold, NetDischarge2[i in rtm,k in dam,s in S], Pnet[i,k,s] >= -P_max)
    # Fixing first stage variables at solutions from rolling horizon with mean-value
    @constraint(m_rold, Fix_PDAM[k in dam[1],s in S], Pdam[k,s] == Pdam_mv_rol[p])
    @objective(m_rold, Min, (1/NS)*sum{-profittotal[s] + sum{unmetcost[i,k,s],i in rtm, k in dam}, s in S})
    tic()
    status = solve(m_rold)
    time_taken_st_rolling = toc();
##########################################################################

    ebat0 = getvalue(getvariable(m_rold,:ebat))[rtm[end],dam[1],realized_sequence[p]];
    Prtm_realized[:,p] = getvalue(getvariable(m_rold,:Prtm))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetload_realized[:,p] = getvalue(getvariable(m_rold,:unmetload))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetcost_realized[:,p] = getvalue(getvariable(m_rold,:unmetcost))[1:rtm[end],dam[1],realized_sequence[p]];
    profitErtm_realized[:,p] = getvalue(getvariable(m_rold,:profitErtm))[1:rtm[end],dam[1],realized_sequence[p]];
    profitEdam_realized[p] = getvalue(getvariable(m_rold,:profitEdam))[dam[1],realized_sequence[p]];
    profittotal_realized[p] = sum(profitErtm_realized[:,p]) + profitEdam_realized[p];
    netobjective_realized[p] = sum(unmetcost_realized[:,p]) - profittotal_realized[p];
    println("Scenario $k, Step $p, netobjective_realized = $(netobjective_realized[p]), time = $(time_taken_st_rolling)")





#=

    println("\nTotal Profits on day ", p, ": ", getValue(profittotal),"\n")
    dailyprofit[j] = getValue(profittotal);

    n2 = [nrtm,ndam]

    loadplot[(j-1)*nrtm+1:j*nrtm] = reshape(load,nrtmpoints);
    rtmeprplot[(j-1)*nrtm+1:j*nrtm] = reshape(rtmepr,nrtmpoints);
    dameprplot[j] = reshape(damepr,ndampoints);
    damregupprplot[j] = reshape(damreguppr,ndampoints);
    damregdownprplot[j] = reshape(damregdownpr,ndampoints);
    ebatplot[(j-1)*nrtm+1:j*nrtm] = reshape(convertToArray2(ebat,n2),nrtmpoints);
    socplot[(j-1)*nrtm+1:j*nrtm] = reshape(convertToArray2(soc,n2),nrtmpoints);
    Prtmplot[(j-1)*nrtm+1:j*nrtm] = reshape(convertToArray2(Prtm,n2),nrtmpoints);
    Pdamplot[j] = reshape(convertToArray(Pdam),ndampoints);
    regupdamplot[j] = reshape(convertToArray(regupdam),ndampoints);
    regdowndamplot[j] = reshape(convertToArray(regdowndam),ndampoints);
    rtmeprofitplot[(j-1)*nrtm+1:j*nrtm] = reshape(convertToArray2(profitErtm,n2),nrtmpoints);
    dameprofitplot[j] = reshape(convertToArray(profitEdam),ndampoints);
    damregprofitplot[j] = reshape(convertToArray(profitregupdam),ndampoints) + reshape(convertToArray(profitregdowndam),ndampoints);


Prtmarray = convertToArray2(Prtm,n2);
Pdamarray = convertToArray(Pdam);
regupdamarray = convertToArray(regupdam);
regdowndamarray = convertToArray(regdowndam);


for i in rtm
    for k in dam
        totalpower[i,k,j] = Prtmarray[i,k] + Pdamarray[k];
        totalregup[i,k,j] = regupdamarray[k];
        totalregdown[i,k,j] = regdowndamarray[k];
        upperband[i,k,j] = totalpower[i,k,j] + totalregup[i,k,j];
        lowerband[i,k,j] = totalpower[i,k,j] - totalregdown[i,k,j];
    end
end

=#

j = j+1;

end # End rolling horizon mean-value problem

totalcost_after_rolling_dt = sum(netobjective_realized);

obj_dt_rh = totalcost_after_rolling_dt;


push!(obj_dt_rh_NS,obj_dt_rh)

end


expected_obj_dt_rh = mean(obj_dt_rh_NS);



#= Start comment here


println("Total Profits: ", sum(dailyprofit),"\n")



totalpowerplot = reshape(totalpower,nrtmpoints*ndays);
upperbandplot = reshape(upperband,nrtmpoints*ndays);
lowerbandplot = reshape(lowerband,nrtmpoints*ndays);





#################### Plotting #####################

timertm = 0:dtrtm:ndam*ndays;
timedam = 0:dtdam:ndam*ndays;

loadplot = push!(loadplot,loadplot[end]);
rtmeprplot = push!(rtmeprplot,rtmeprplot[end]);
dameprplot = push!(dameprplot,dameprplot[end]);
damregupprplot = push!(damregupprplot,damregupprplot[end]);
damregdownprplot = push!(damregdownprplot,damregdownprplot[end]);
#ebatplot = push!(ebatplot,ebatplot[end]);
ebatplot = [ebat_max,ebatplot];
#socplot = push!(socplot,socplot[end]);
socplot = [soc0,socplot];
Prtmplot = push!(Prtmplot,Prtmplot[end]);
Pdamplot = push!(Pdamplot,Pdamplot[end]);
regupdamplot = push!(regupdamplot,regupdamplot[end]);
regdowndamplot = push!(regdowndamplot,regdowndamplot[end]);
rtmeprofitplot = push!(rtmeprofitplot,rtmeprofitplot[end]);
dameprofitplot = push!(dameprofitplot,dameprofitplot[end]);
damregprofitplot = push!(damregprofitplot,damregprofitplot[end]);



totalpowerplot = push!(totalpowerplot,totalpowerplot[end]);
upperbandplot = push!(upperbandplot,upperbandplot[end]);
lowerbandplot = push!(lowerbandplot,lowerbandplot[end]);




function cumul(A)
    C = zeros(length(A)-1);
    for i in 1:length(A)-1
        C[i] = sum(A[1:i]);
    end
    C = push!(C,C[end]);
    return C;
end

function damtortm(A)
    A = A[1:ndampoints*ndays]
    B = zeros(nrtmpoints*ndays);
    for i in 1:ndampoints*ndays
        B[(i-1)*12+1:i*12] = repmat([A[i]],12);
    end
    B = push!(B,B[end])
    return B
end

rtmeprofitplot = cumul(rtmeprofitplot);
dameprofitplot = cumul(dameprofitplot);
totaleprofitplot = rtmeprofitplot + damtortm(dameprofitplot);
damregprofitplot = cumul(damregprofitplot);
totalregprofitplot = damtortm(damregprofitplot);
totalprofitplot = totaleprofitplot + totalregprofitplot;


figure()
plot(timertm,loadplot,drawstyle="steps-post",label="Load (MW)",LineWidth=1.5)
#hold(true)
#plot(timertm,[loadmeanplot,loadmeanplot[end]])
xlabel("Time (hr)",size=20)
ylabel("Load (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timertm[(1:nrtmpoints+1)],rtmeprplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="RTM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("RTM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("RTM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timertm[(1:nrtmpoints+1)],Prtmplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam[(1:ndampoints+1)],dameprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="DAM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam[(1:ndampoints+1)],Pdamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam[(1:ndampoints+1)],damregupprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timedam[(1:ndampoints+1)],damregdownprplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam[(1:ndampoints+1)],regupdamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in DAM",LineWidth=1.5)
hold(true)
plot(timedam[(1:ndampoints+1)],regdowndamplot[91*ndampoints+(1:ndampoints+1)],drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in DAM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm[(1:nrtmpoints+1)],socplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",label="SOC",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("SOC (%)",size=20)
grid()
tick_params(labelsize=18)
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")


figure()
plot(timertm[(1:nrtmpoints+1)],totalpowerplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="blue",label="Total Power (MW)",LineWidth=1.5)
hold(true)
plot(timertm[(1:nrtmpoints+1)],lowerbandplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="red",label="Lower band (MW)",LineWidth=1.5)
plot(timertm[(1:nrtmpoints+1)],upperbandplot[91*nrtmpoints+(1:nrtmpoints+1)],drawstyle="steps-post",color="green",label="Upper band (MW)",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("Power (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")




figure()
plot(timertm/24,rtmeprofitplot,drawstyle="steps-post",color="blue",label="RTM",LineWidth=1.5)
xlabel("Time (day)",size=20), ylabel("Revenues from Energy (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,dameprofitplot,drawstyle="steps-post",color="cyan",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totaleprofitplot,drawstyle="steps-post",color="red",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
xlabel("Time (day)",size=20), ylabel("Revenues from Reg (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,damregprofitplot,drawstyle="steps-post",color="cyan",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totalregprofitplot,drawstyle="steps-post",color="red",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm/24,totalprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
xlabel("Time (day)",size=14), ylabel("Total Revenue (\$)",size=20)
grid();tick_params(labelsize=18);#xlim(0,370)



#End comment here
=#







#=







figure()
plot(timertm,loadplot,drawstyle="steps-post",label="Load (MW)",LineWidth=1.5)
#hold(true)
#plot(timertm,[loadmeanplot,loadmeanplot[end]])
xlabel("Time (hr)",size=20)
ylabel("Load (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Load (MW)")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timertm,rtmeprplot,drawstyle="steps-post",color="blue",label="RTM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("RTM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("RTM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timertm,Prtmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("RTM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in RTM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam,dameprplot,drawstyle="steps-post",color="blue",label="DAM Energy price (USD/MWh)",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Energy prices (\$/MWh)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM energy prices")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam,edamplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM",LineWidth=1.5)
ylim(-1.1*P_max,1.1*P_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Power sold (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("Power sold in DAM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam,damregupprplot,drawstyle="steps-post",color="blue",label="Reg Up",LineWidth=1.5)
hold(true)
plot(timedam,damregdownprplot,drawstyle="steps-post",color="red",label="Reg Down",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg prices (\$/MW)",size=20)
grid()
tick_params(labelsize=18)
#title("DAM RegUp prices")
legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam,regupdamplot,drawstyle="steps-post",color="blue",label="Regulation Up committed (MW) in DAM",LineWidth=1.5)
hold(true)
plot(timedam,regdowndamplot,drawstyle="steps-post",color="red",label="Regulation Down committed (MW) in DAM",LineWidth=1.5)
ylim(-0.1,1.1*regup_max)
xlabel("Time (hr)",size=20)
ylabel("DAM Reg (MW)",size=20)
grid()
tick_params(labelsize=18)
#title("RegUp in DAM")
#legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm,socplot,drawstyle="steps-post",label="SOC",LineWidth=1.5)
xlabel("Time (hr)",size=20)
ylabel("SOC (%)",size=20)
grid()
tick_params(labelsize=18)
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")



figure()
plot(timertm/24,rtmeprofitplot,drawstyle="steps-post",label="RTM",LineWidth=1.5)
xlabel("Time (day)",size=20), ylabel("Revenues from Energy (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,dameprofitplot,drawstyle="steps-post",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totaleprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
xlabel("Time (day)",size=20), ylabel("Revenues from Reg (\$)",size=20)
#xlim(0,370)
grid();tick_params(labelsize=18)
hold(true)
plot(timedam/24,damregprofitplot,drawstyle="steps-post",label="DAM",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
plot(timertm/24,totalregprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
#xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
#xlim(0,24)
#grid();tick_params(labelsize=18)
legend(loc=2,fancybox="True", shadow="True")

figure()
plot(timertm/24,totalprofitplot,drawstyle="steps-post",label="Total",LineWidth=1.5)
xlabel("Time (day)",size=14), ylabel("Total Revenue (\$)",size=20)
grid();tick_params(labelsize=18);#xlim(0,370)

=#
