
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
makeplots = 1; # 0: Don't generate plots, 1: Generate plots

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
paths = readcsv("samplepaths.csv");
loadperm = zeros(nrtm,ndam,ndays_planning,NS);
for s in S
  j = 1;
  for l in day
    for k in dam
        loadperm[:,k,l,s] = load[:,k,l,paths[j,s]];
        j = j+1;
    end
  end
end
loadperm_mv = mean(loadperm,4); # Mean loads for mean-value problem

################ Mean Value Model ##################
tic()
mv = Model(solver = GurobiSolver(Threads=2))
		@variable(mv, -P_max <= Prtm[rtm,dam,day] <= P_max)	                #Power sold to the real time market, kW
		@variable(mv, -P_max <= Pdam[dam,day] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(mv, 0 <= ebat[rtm,dam,day] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(mv, suppliedload[rtm,dam,day] >= 0)
		@variable(mv, unmetload[rtm,dam,day] >= 0)
		@expression(mv, Pnet[i in rtm,k in dam,l in day], Prtm[i,k,l] + Pdam[k,l] + suppliedload[i,k,l])    #Net power discharged from battery in all 5-min interval, kW
		@variable(mv, profitErtm[rtm,dam,day])				        #Profit from the real time market, USD
		@variable(mv, profitEdam[dam,day])	        			#Profit from the day ahead market, USD
		@variable(mv, profittotal)	        	#Total profit in the day, USD
		@variable(mv, unmetcost)

		@constraint(mv, InitialEnergy[s in S], ebat[1,1,1] == ebat0 - Pnet[1,1,1]*dtrtm)	#Inital energy in the battery
		@constraint(mv, rtmEBalance[i in rtm[2:end],k in dam,l in day], ebat[i,k,l] == ebat[i-1,k,l] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day], ebat[i,k,l] == ebat[iend,k-1,l] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], ebat[i,k,l] == ebat[iend,kend,l-1] - Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, UnmetLoad[i in rtm,k in dam,l in day], suppliedload[i,k,l] + unmetload[i,k,l] >=  loadperm_mv[i,k,l])
		@constraint(mv, BoundSupplied[i in rtm,k in dam,l in day], suppliedload[i,k,l] <= loadperm_mv[i,k,l])
		@constraint(mv, BoundUnmet[i in rtm,k in dam,l in day], unmetload[i,k,l] <= loadperm_mv[i,k,l])
		#=    @constraint(mv, RTMRamp1[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, RTMRamp2[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(mv, RTMEProfits[i in rtm,k in dam,l in day], profitErtm[i,k,l] == rtmepr[i,k,l]*Prtm[i,k,l]*dtrtm)	#Economic calculation
		@constraint(mv, DAMEProfits[k in dam,l in day], profitEdam[k,l] == damepr[k,l]*Pdam[k,l]*dtdam)	#Economic calculation
		@constraint(mv, TotalProfit, profittotal ==
												sum{profitErtm[i,k,l], i in rtm, k in dam, l in day} + sum{profitEdam[k,l], k in dam, l in day})
		@constraint(mv, UnmetCost, unmetcost == sum{rtmepr[i,k,l]*unmetload[i,k,l], i in rtm, k in dam, l in day})
		@objective(mv, Min, -profittotal + unmetcost)
		if participate_dam == 0
      @constraint(mv, NoDAM[k in dam,l in day, s in S], Pdam[k,l,s] == 0)
    end
    if participate_rtm == 0
      @constraint(mv, NoRTM[i in rtm,k in dam,l in day, s in S], Prtm[i,k,l,s] == 0)
    end
status = solve(mv);
Pdam_mv = getvalue(getvariable(mv,:Pdam));

######################### Mean value Model End ###########################

################ Evaluating costs in all scenarios with mv-solution ##################

m_d = Model(solver = GurobiSolver(Threads=2))
		@variable(m_d, -P_max <= Prtm[rtm,dam,day,S] <= P_max)	                #Power sold to the real time market, kW
		@variable(m_d, -P_max <= Pdam[dam,day,S] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(m_d, 0 <= ebat[rtm,dam,day,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(m_d, suppliedload[rtm,dam,day,S] >= 0)
		@variable(m_d, unmetload[rtm,dam,day,S] >= 0)
		@expression(m_d, Pnet[i in rtm,k in dam,l in day,s in S], Prtm[i,k,l,s] + Pdam[k,l,s] + suppliedload[i,k,l,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m_d, profitErtm[rtm,dam,day,S])				        #Profit from the real time market, USD
		@variable(m_d, profitEdam[dam,day,S])	        			#Profit from the day ahead market, USD
		@variable(m_d, profittotal[S])		        	#Total profit in the day, USD
		@variable(m_d, unmetcost[S])

		@constraint(m_d, InitialEnergy[s in S], ebat[1,1,1,s] == ebat0 - Pnet[1,1,1,s]*dtrtm)	#Inital energy in the battery
		@constraint(m_d, rtmEBalance[i in rtm[2:end],k in dam,l in day,s in S], ebat[i,k,l,s] == ebat[i-1,k,l,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], ebat[i,k,l,s] == ebat[iend,k-1,l,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], ebat[i,k,l,s] == ebat[iend,kend,l-1,s] - Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m_d, UnmetLoad[i in rtm,k in dam,l in day, s in S], suppliedload[i,k,l,s] + unmetload[i,k,l,s] >=  loadperm[i,k,l,s])
		@constraint(m_d, BoundSupplied[i in rtm,k in dam,l in day,s in S], suppliedload[i,k,l,s] <= loadperm[i,k,l,s])
		@constraint(m_d, BoundUnmet[i in rtm,k in dam,l in day,s in S], unmetload[i,k,l,s] <= loadperm[i,k,l,s])
		#=    @constraint(m_d, RTMRamp1[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, RTMRamp2[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_d, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(m_d, RTMEProfits[i in rtm,k in dam,l in day,s in S], profitErtm[i,k,l,s] == rtmepr[i,k,l]*Prtm[i,k,l,s]*dtrtm)	#Economic calculation
		@constraint(m_d, DAMEProfits[k in dam,l in day,s in S], profitEdam[k,l,s] == damepr[k,l]*Pdam[k,l,s]*dtdam)	#Economic calculation
		@constraint(m_d, TotalProfit[s in S], profittotal[s] ==
												sum{profitErtm[i,k,l,s], i in rtm, k in dam, l in day} + sum{profitEdam[k,l,s], k in dam, l in day})
		@constraint(m_d, UnmetCost[s in S], unmetcost[s] == sum{rtmepr[i,k,l]*unmetload[i,k,l,s], i in rtm, k in dam, l in day})
    # Fixing first stage variables at mean-value solution
		@constraint(m_d, Fix_PDAM[k in dam,l in day,s in S], Pdam[k,l,s] == Pdam_mv[k,l])
    @objective(m_d, Min, (1/NS)*sum{-profittotal[s] + unmetcost[s], s in S})
		if participate_dam == 0
      @constraint(m_d, NoDAM[k in dam,l in day, s in S], Pdam[k,l,s] == 0)
    end
    if participate_rtm == 0
      @constraint(m_d, NoRTM[i in rtm,k in dam,l in day, s in S], Prtm[i,k,l,s] == 0)
    end
		status = solve(m_d)

time_taken_dt_fullproblem = toc();
###############################################################

println("\nExpected Objective with Mean Value first stage solution ", getobjectivevalue(m_d),"\n" )

obj_dt_fp = getobjectivevalue(m_d);


if makeplots == 1
################# PLOTTING #################
#=
# Plot of Scenarios of loads
xplot = 0:dtrtm:dtrtm*nrtm_planning
loadNSplot = loadNSdata[1:nrtm_planning,:];
loadNSplot = [loadNSplot;loadNSplot[end,:]];
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S
plot(xplot,loadNSplot[:,s], color="grey", drawstyle="steps-post");
end
plot(xplot,mean(loadNSplot,2), color="blue", drawstyle="steps-post",label="Mean scenario");
grid()
xlim(0,nhours_planning)
ylabel("Loads (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/loads_scenarios.pdf"))
close("all")
=#

r = 5;
# Plot of Pdam
n3 = [ndam,ndays_planning,NS];
eprdamplot = eprdam[1:nhours_planning];
push!(eprdamplot,eprdamplot[end]);
Pdamarray = convertToArray3(getvalue(getvariable(m_d,:Pdam)),n3);
Pdamplot = reshape(Pdamarray[:,:,1],nhours_planning);
push!(Pdamplot,Pdamplot[end]);
xplot = 0:dtdam:dtdam*nhours_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
subplot(2,1,1)
hold(true)
plot(xplot,eprdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Energy price (\$/MWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
plot(xplot,Pdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("DAM Power (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Pdam_fp_dt.pdf"))
close("all")

# Plot of Prtm
n4 = [nrtm,ndam,ndays_planning,NS];
eprrtmplot = eprdam[1:nrtm_planning];
push!(eprrtmplot,eprrtmplot[end]);
Prtmarray = convertToArray4(getvalue(getvariable(m_d,:Prtm)),n4);
Prtmplot = reshape(Prtmarray,nrtm_planning,NS);
Prtmplot = [Prtmplot;Prtmplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
subplot(2,1,1)
hold(true)
plot(xplot,eprrtmplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Energy price (\$/MWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
for s in S[r]
	plot(xplot,Prtmplot[:,s],  drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("RTM Power (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Prtm_fp_dt.pdf"))
close("all")

# Plot of supplied load and unmet load
n4 = [nrtm,ndam,ndays_planning,NS];
suppliedloadarray = convertToArray4(getvalue(getvariable(m_d,:suppliedload)),n4);
suppliedloadplot = reshape(suppliedloadarray,nrtm_planning,NS);
suppliedloadplot = [suppliedloadplot;suppliedloadplot[end,:]];
unmetloadarray = convertToArray4(getvalue(getvariable(m_d,:unmetload)),n4);
unmetloadplot = reshape(unmetloadarray,nrtm_planning,NS);
unmetloadplot = [unmetloadplot;suppliedloadplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,suppliedloadplot[:,s], drawstyle="steps-post", color = "green", label="Supplied");
	plot(xplot,unmetloadplot[:,s], drawstyle="steps-post", color = "red", label="Unmet");
end
grid()
xlim(0,nhours_planning)
ylabel("Supplied & unmet loads (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/supp_unmet_fp_dt.pdf"))
close("all")


# Plot of SOC
n4 = [nrtm,ndam,ndays_planning,NS];
socarray = convertToArray4(getvalue(getvariable(m_d,:ebat)),n4)/ebat_max*100;
socplot = reshape(socarray,nrtm_planning,NS);
socplot = [socplot;socplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,socplot[:,s], drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylim(-0.02,100.02)
ylabel("State of charge (%)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/soc_fp_dt.pdf"))
close("all")

# Pnet, regulation bands calculation
netpower = zeros(nrtm,ndam,ndays_planning,NS);
for i in rtm
        for k in dam
            for l in day
							for s in S
                netpower[i,k,l,s] = Prtmarray[i,k,l,s] + Pdamarray[k,l,s] + suppliedloadarray[i,k,l,s];
							end
            end
        end
end
netpowerplot = reshape(netpower,nrtm_planning,NS);
netpowerplot = [netpowerplot;netpowerplot[end,:]];

# Plot of netpower
n4 = [nrtm,ndam,ndays_planning,NS];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[r]
	plot(xplot,netpowerplot[:,s], drawstyle="steps-post", color = "blue", label="Net discharge");
end
grid()
xlim(0,nhours_planning)
ylim(-1.02,1.02)
ylabel("Net discharge (MW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/netpower_fp_dt.pdf"))
close("all")

# Plot of Profits from Pdam and Prtm and Total
n3 = [ndam,ndays_planning,NS];
profitEdamarray = convertToArray3(getvalue(getvariable(m_d,:profitEdam)),n3);
profitEdamplot = reshape(profitEdamarray[:,:,1],nhours_planning);
cumulprofitEdamplot = cumul(profitEdamplot);
cumulprofitEdamtortmplot = damtortm(cumulprofitEdamplot);

n4 = [nrtm,ndam,ndays_planning,NS];
profitErtmarray = convertToArray4(getvalue(getvariable(m_d,:profitErtm)),n4);
profitErtmplot = reshape(profitErtmarray,nrtm_planning,NS);
profitErtmplotmean = mean(profitErtmplot,2);
cumulprofitErtmplot = zeros(nrtm_planning,NS);
for s in S
  cumulprofitErtmplot[:,s] = cumul(profitErtmplot[:,s]);
end
cumulprofitErtmplotmean = cumul(profitErtmplotmean);


n4 = [nrtm,ndam,ndays_planning,NS];
eprrtmplot = eprrtm[1:nrtm_planning];
unmetloadarray = convertToArray4(getvalue(getvariable(m_d,:unmetload)),n4);
unmetloadplot = reshape(unmetloadarray,nrtm_planning,NS);
costunmetloadplot = zeros(nrtm_planning,NS);
cumulcostunmetloadplot = zeros(nrtm_planning,NS);
for s in S
  costunmetloadplot[:,s] = eprrtmplot.*unmetloadplot[:,s];
  cumulcostunmetloadplot[:,s] = cumul(costunmetloadplot[:,s]);
end
cumulcostunmetloadplotmean = mean(cumulcostunmetloadplot,2);
cumultotalprofitplot = repmat(cumulprofitEdamtortmplot,1,NS) + cumulprofitErtmplot - cumulcostunmetloadplot;
cumultotalprofitplotmean = mean(cumultotalprofitplot,2);

cumulprofitEdamtortmplot = [cumulprofitEdamtortmplot;cumulprofitEdamtortmplot[end]];
cumulprofitErtmplotmean = [cumulprofitErtmplotmean;cumulprofitErtmplotmean[end]];
cumulcostunmetloadplotmean = [cumulcostunmetloadplotmean;cumulcostunmetloadplotmean[end]];
cumultotalprofitplotmean = [cumultotalprofitplotmean;cumultotalprofitplotmean[end]];
cumulprofitErtmplot = [cumulprofitErtmplot;cumulprofitErtmplot[end,:]];
cumulcostunmetloadplot = [cumulcostunmetloadplot;cumulcostunmetloadplot[end,:]];
cumultotalprofitplot = [cumultotalprofitplot;cumultotalprofitplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
plot(xplot,cumulprofitEdamtortmplot, drawstyle="steps-post", color = "blue", label="Day-ahead market");
plot(xplot,cumulprofitErtmplotmean, drawstyle="steps-post", color = "green", label="Real-time market");
plot(xplot,cumulcostunmetloadplotmean, drawstyle="steps-post", color = "red", label="Unmet load");
plot(xplot,cumultotalprofitplotmean, drawstyle="steps-post", color = "indigo", label="Total revenue");
grid()
xlim(0,nhours_planning)
#ylim(-1.02,1.02)
ylabel("Expected Cumulative Revenue (\$)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=20)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/cumulative_rev_fp_dt.pdf"))
close("all")

end #END IF










#    println("\nTotal Profits ", getvalue(profittotal),"\n" )

#= #Start commenting here

dailyprofit = zeros(ndays_planning);
loadplot = zeros(ndays_planning*nrtmpoints);
rtmeprplot = zeros(ndays_planning*nrtmpoints);
dameprplot = zeros(ndays_planning*ndampoints);
damregupprplot = zeros(ndays_planning*ndampoints);
damregdownprplot = zeros(ndays_planning*ndampoints);
ebatplot = zeros(ndays_planning*nrtmpoints);
socplot = zeros(ndays_planning*nrtmpoints);
Prtmplot = zeros(ndays_planning*nrtmpoints);
Pdamplot = zeros(ndays_planning*ndampoints);
regupdamplot = zeros(ndays_planning*ndampoints);
regdowndamplot = zeros(ndays_planning*ndampoints);







n3 = [nrtm,ndam,ndays_planning]
n2 = [ndam,ndays_planning]

loadplot = reshape(load,nrtmpoints*ndays);
rtmeprplot = reshape(rtmepr,nrtmpoints*ndays);
dameprplot = reshape(damepr,ndampoints*ndays);
damregupprplot = reshape(damreguppr,ndampoints*ndays);
damregdownprplot = reshape(damregdownpr,ndampoints*ndays);
ebatplot = reshape(convertToArray3(ebat,n3),nrtmpoints*ndays);
socplot = reshape(convertToArray3(soc,n3),nrtmpoints*ndays);
Prtmplot = reshape(convertToArray3(Prtm,n3),nrtmpoints*ndays);
Pdamplot = reshape(convertToArray2(Pdam,n2),ndampoints*ndays);
regupdamplot = reshape(convertToArray2(regupdam,n2),ndampoints*ndays);
regdowndamplot = reshape(convertToArray2(regdowndam,n2),ndampoints*ndays);
rtmeprofitplot = reshape(convertToArray3(profitErtm,n3),nrtmpoints*ndays);
dameprofitplot = reshape(convertToArray2(profitEdam,n2),ndampoints*ndays);
damregprofitplot = reshape(convertToArray2(profitregupdam,n2),ndampoints*ndays) + reshape(convertToArray2(profitregdowndam,n2),ndampoints*ndays);


totalpower = zeros(nrtm,ndam,ndays);
totalregup = zeros(nrtm,ndam,ndays);
totalregdown = zeros(nrtm,ndam,ndays);
upperband = zeros(nrtm,ndam,ndays);
lowerband = zeros(nrtm,ndam,ndays);
totalpowerplot = zeros(nrtm,ndam,ndays);
upperbandplot = zeros(nrtm,ndam,ndays);
lowerbandplot = zeros(nrtm,ndam,ndays);

Prtmarray = convertToArray3(Prtm,n3);
Pdamarray = convertToArray2(Pdam,n2);
regupdamarray = convertToArray2(regupdam,n2);
regdowndamarray = convertToArray2(regdowndam,n2);

for i in rtm
        for k in dam
            for l in day
                totalpower[i,k,l] = Prtmarray[i,k,l] + Pdamarray[k,l];
                totalregup[i,k,l] = regupdamarray[k,l];
                totalregdown[i,k,l] = regdowndamarray[k,l];
                upperband[i,k,l] = totalpower[i,k,l] + totalregup[i,k,l];
                lowerband[i,k,l] = totalpower[i,k,l] - totalregdown[i,k,l];
            end
        end
end


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

=# #End commenting here
