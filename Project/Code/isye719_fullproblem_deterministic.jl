
using JuMP
using PyPlot
using Gurobi

close("all")

############# Funtions to convert JuMP returned variables to arrays ################

function convertToArray(x)
	y = getvalue(x)
	n = length(y)
	a = zeros(n)
	for i = 1:n
		a[i] = y[i]
	end
	return a
end

function convertToArray2(AA,n)
	AA = getvalue(A)
	m = (n[1],n[2])
	B = zeros(m)
	for i in 1:n[1]
		for j in 1:n[2]
			B[i,j] = AA[i,j]
		end
	end
	return B
end

function convertToArray3(AA,n)
	m = (n[1],n[2],n[3])
	B = zeros(m)
	for i in 1:n[1]
		for j in 1:n[2]
			for k in 1:n[3]
				B[i,j,k] = AA[i,j,k]
			end
		end
	end
	return B
end

function convertToArray4(AA,n)
	m = (n[1],n[2],n[3],n[4])
	B = zeros(m)
	for i in 1:n[1]
		for j in 1:n[2]
			for k in 1:n[3]
				for l in 1:n[4]
					B[i,j,k,l] = AA[i,j,k,l]
				end
			end
		end
	end
	return B
end

#############################################################################################################

loaddata = readcsv("DATARCoast.csv");
rtmpricedata = readcsv("AggregatedData_RTM_ALTA31GT_7_B1.csv")
dampricedata = readcsv("AggregatedData_DAM_ALTA31GT_7_B1.csv")

dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day ahead market (hourly) intervals in a day
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtdam/dtrtm);		#No. of real time intervals in each hour

nrtmpoints = ndam*nrtm;		        #Total number of points in RTM data in one day
ndampoints = ndam;			#Total number of points in hourly data in one day

#Model Parameters
ebat_max = 0.5;	          #Battery capacity, MWh
P_max = 1;	          #Maximum power, MW
regup_max = P_max;    #Regulation Up Capacity, MW
regdown_max = P_max;  #Regulation Up Capacity, MW
rampmin = -0.5*P_max;	          #Lower bound for ramp discharge, MW/5min
rampmax = 0.5*P_max;  	  #Upper bound for ramp discharge, MW/5min
eff = 1;                  #Discharging Efficiency of battery
ebat0 = ebat_max;		   #Initial State of charge
ebatend = ebat_max;		  #State of charge at the end of the day
ndays = 365;              #Number of days data is available for
ndays_planning = 7;       #Number of days you want to plan for
nhours_planning = ndays_planning*ndam;
nrtm_planning = nhours_planning*nrtm;

NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;




#Load and price data
load = Matrix{Float64}(loaddata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW

eprrtm = rtmpricedata[1:nrtmpoints*ndays,4];	    	#Real Time Market price, $/MWh
eprdam = dampricedata[1:ndampoints*ndays,4];	    	#Day Ahead Market Selling price, $/MWh
regupprdam = dampricedata[1:ndampoints*ndays,5];	    	#Day Ahead Market Regulation up price, $/MWh
regdownprdam = dampricedata[1:ndampoints*ndays,6]; 		#Day Ahead Market Regulation down price, $/MWh


#Reshape the data to matrices
rtmepr = reshape(eprrtm,nrtm,ndam,ndays);
damepr = reshape(eprdam,ndam,ndays);
damreguppr = reshape(regupprdam,ndam,ndays);
damregdownpr = reshape(regdownprdam,ndam,ndays);
load = reshape(load,nrtm,ndam,ndays);
loadvec = vec(load);



# Reshaping the data as weekly profiles
reshape_ndays = 7;
reshape_nrows = reshape_ndays*nrtmpoints;
reshape_ncolumns = Int64(floor(length(load)/reshape_nrows));
load_estimationdata = loadvec[1:reshape_nrows*reshape_ncolumns];
load_weekly = reshape(load_estimationdata,reshape_nrows,reshape_ncolumns);



# Loading the NS scenarios for weekly load profiles in kW generated from the fullproblem_stochastic code
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));

loadNSdata = readcsv("loads_scenarios_month.csv")


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

################################################################


load = loadNSplanningdata[:,:,1:ndays_planning,:]/1000; #MW

load_mv = mean(load,4);

#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
dam = 1:ndam;
day = 1:ndays_planning;




################ Mean Value Model ##################
tic()
mv = Model(solver = GurobiSolver(Threads=2))
		@variable(mv, -P_max <= Prtm[rtm,dam,day] <= P_max)	                #Power sold to the real time market, kW
		@variable(mv, -P_max <= Pdam[dam,day] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(mv, 0 <= ebat[rtm,dam,day] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(mv, suppliedload[rtm,dam,day] >= 0)
		@variable(mv, unmetload[rtm,dam,day] >= 0)
		@expression(mv, Pnet[i in rtm,k in dam,l in day], Prtm[i,k,l] + Pdam[k,l] + suppliedload[i,k,l])    #Net power discharged from battery in all 5-min interval, kW
		@variable(mv, profitErtm[rtm,dam,day])# >= 0)				        #Profit from the real time market, USD
		@variable(mv, profitEdam[dam,day])# >= 0)	        			#Profit from the day ahead market, USD
		@variable(mv, profittotal)# >= 0)		        	#Total profit in the day, USD
		@variable(mv, unmetcost)

		@constraint(mv, InitialEnergy[s in S], ebat[1,1,1] == ebat0 - 1/eff*Pnet[1,1,1]*dtrtm)	#Inital energy in the battery
		#    @constraint(mv, EndSOC[i in rtm,k in dam,l in day,s in S], soc[i,k,l,s] >= socend)		#Constraint on SOC at the end of the day
		@constraint(mv, rtmEBalance[i in rtm[2:end],k in dam,l in day], ebat[i,k,l] == ebat[i-1,k,l] - 1/eff*Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day], ebat[i,k,l] == ebat[iend,k-1,l] - 1/eff*Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], ebat[i,k,l] == ebat[iend,kend,l-1] - 1/eff*Pnet[i,k,l]*dtrtm)	#Dynamics constraint
		@constraint(mv, UnmetLoad[i in rtm,k in dam,l in day], suppliedload[i,k,l] + unmetload[i,k,l] >=  load_mv[i,k,l])
		@constraint(mv, BoundSupplied[i in rtm,k in dam,l in day], suppliedload[i,k,l] <= load_mv[i,k,l])
		@constraint(mv, BoundUnmet[i in rtm,k in dam,l in day], unmetload[i,k,l] <= load_mv[i,k,l])
		#=    @constraint(mv, RTMRamp1[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, RTMRamp2[i in rtm[2:end],k in dam,l in day], Pnet[i,k,l]  - Pnet[i-1,k,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day], Pnet[i,k,l] - Pnet[iend,k-1,l] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(mv, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end]], Pnet[i,k,l] - Pnet[iend,kend,l-1] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(mv, RTMEProfits[i in rtm,k in dam,l in day], profitErtm[i,k,l,s] == rtmepr[i,k,l]*Prtm[i,k,l,s]*dtrtm)	#Economic calculation
		@constraint(mv, DAMEProfits[k in dam,l in day], profitEdam[k,l] == damepr[k,l]*Pdam[k,l]*dtdam)	#Economic calculation
		@constraint(mv, TotalProfit, profittotal ==
												sum{profitErtm[i,k,l], i in rtm, k in dam, l in day} + sum{profitEdam[k,l], k in dam, l in day})
		@constraint(mv, UnmetCost, unmetcost == sum{rtmepr[i,k,l]*unmetload[i,k,l], i in rtm, k in dam, l in day})
		# Non-anticipativity constraints for first stage variables
		@constraint(mv, Nonant_PDAM[k in dam,l in day], Pdam[k,l] == Pdam[k,l])
		@objective(mv, Min, -profittotal + unmetcost[s])
#    print(mv)
status = solve(mv)
Pdam_mv = getvalue(getvariable(mv,:Pdam))


###############################################################

######################### Mean value Model End ###########################

################ Evaluating costs in all scenarios with mv-solution ##################

m = Model(solver = GurobiSolver(Threads=2))
		@variable(m, -P_max <= Prtm[rtm,dam,day,S] <= P_max)	                #Power sold to the real time market, kW
		@variable(m, -P_max <= Pdam[dam,day,S] <= P_max)    	                #Power sold to the day ahead market, kW
		@variable(m, 0 <= ebat[rtm,dam,day,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
		@variable(m, suppliedload[rtm,dam,day,S] >= 0)
		@variable(m, unmetload[rtm,dam,day,S] >= 0)
		@expression(m, Pnet[i in rtm,k in dam,l in day,s in S], Prtm[i,k,l,s] + Pdam[k,l,s] + suppliedload[i,k,l,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m, profitErtm[rtm,dam,day,S])# >= 0)				        #Profit from the real time market, USD
		@variable(m, profitEdam[dam,day,S])# >= 0)	        			#Profit from the day ahead market, USD
		@variable(m, profittotal[S])# >= 0)		        	#Total profit in the day, USD
		@variable(m, unmetcost[S])

		@constraint(m, InitialEnergy[s in S], ebat[1,1,1,s] == ebat0 - 1/eff*Pnet[1,1,1,s]*dtrtm)	#Inital energy in the battery
		#    @constraint(m, EndSOC[i in rtm,k in dam,l in day,s in S], soc[i,k,l,s] >= socend)		#Constraint on SOC at the end of the day
		@constraint(m, rtmEBalance[i in rtm[2:end],k in dam,l in day,s in S], ebat[i,k,l,s] == ebat[i-1,k,l,s] - 1/eff*Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], ebat[i,k,l,s] == ebat[iend,k-1,l,s] - 1/eff*Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m, dayEBalance[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], ebat[i,k,l,s] == ebat[iend,kend,l-1,s] - 1/eff*Pnet[i,k,l,s]*dtrtm)	#Dynamics constraint
		@constraint(m, UnmetLoad[i in rtm,k in dam,l in day, s in S], suppliedload[i,k,l,s] + unmetload[i,k,l,s] >=  load[i,k,l,s])
		@constraint(m, BoundSupplied[i in rtm,k in dam,l in day,s in S], suppliedload[i,k,l,s] <= load[i,k,l,s])
		@constraint(m, BoundUnmet[i in rtm,k in dam,l in day,s in S], unmetload[i,k,l,s] <= load[i,k,l,s])
		#=    @constraint(m, RTMRamp1[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, RTMRamp2[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		=#
		@constraint(m, RTMEProfits[i in rtm,k in dam,l in day,s in S], profitErtm[i,k,l,s] == rtmepr[i,k,l]*Prtm[i,k,l,s]*dtrtm)	#Economic calculation
		@constraint(m, DAMEProfits[k in dam,l in day,s in S], profitEdam[k,l,s] == damepr[k,l]*Pdam[k,l,s]*dtdam)	#Economic calculation
		@constraint(m, TotalProfit[s in S], profittotal[s] ==
												sum{profitErtm[i,k,l,s], i in rtm, k in dam, l in day} + sum{profitEdam[k,l,s], k in dam, l in day})
		@constraint(m, UnmetCost[s in S], unmetcost[s] == sum{rtmepr[i,k,l]*unmetload[i,k,l,s], i in rtm, k in dam, l in day})
    # Fixing first stage variables at mean-value solution
		@constraint(m, Fix_PDAM[k in dam,l in day,s in S], Pdam[k,l,s] == Pdam_mv[k,l])
    @objective(m, Min, (1/NS)*sum{-profittotal[s] + unmetcost[s], s in S})
#    print(m)
    status = solve(m)

time_taken_dt_fullproblem = toc();
###############################################################

println("\nExpected Objective with Mean Value first stage solution ", getobjectivevalue(m),"\n" )

obj_dt_fp = getobjectivevalue(m);


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

# Plot of Pdam
n3 = [ndam,ndays_planning,NS];
eprdamplot = eprdam[1:nhours_planning];
push!(eprdamplot,eprdamplot[end]);
Pdamarray = convertToArray3(getvalue(getvariable(m,:Pdam)),n3);
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
ylabel("Energy price (\$/kWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
plot(xplot,Pdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Net Power (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Pdam_fp_dt.pdf"))
close("all")

# Plot of Prtm
n4 = [nrtm,ndam,ndays_planning,NS];
eprrtmplot = eprdam[1:nrtm_planning];
push!(eprrtmplot,eprrtmplot[end]);
Prtmarray = convertToArray4(getvalue(getvariable(m,:Prtm)),n4);
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
ylabel("Energy price (\$/kWh)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
subplot(2,1,2)
hold(true)
for s in S
	plot(xplot,Prtmplot[:,s],  drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylabel("Net Power (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/Prtm_fp_dt.pdf"))
close("all")

# Plot of supplied load and unmet load
n4 = [nrtm,ndam,ndays_planning,NS];
suppliedloadarray = convertToArray4(getvalue(getvariable(m,:suppliedload)),n4);
suppliedloadplot = reshape(suppliedloadarray,nrtm_planning,NS);
suppliedloadplot = [suppliedloadplot;suppliedloadplot[end,:]];
unmetloadarray = convertToArray4(getvalue(getvariable(m,:unmetload)),n4);
unmetloadplot = reshape(unmetloadarray,nrtm_planning,NS);
unmetloadplot = [unmetloadplot;suppliedloadplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[1]
	plot(xplot,suppliedloadplot[:,s], drawstyle="steps-post", color = "green", label="Supplied");
	plot(xplot,unmetloadplot[:,s], drawstyle="steps-post", color = "red", label="Unmet");
end
grid()
xlim(0,nhours_planning)
ylabel("Supplied & unmet loads (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/supp_unmet_fp_dt.pdf"))
close("all")


# Plot of SOC
n4 = [nrtm,ndam,ndays_planning,NS];
socarray = convertToArray4(getvalue(getvariable(m,:ebat)),n4)/ebat_max;
socplot = reshape(socarray,nrtm_planning,NS);
socplot = [socplot;socplot[end,:]];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[1]
	plot(xplot,socplot[:,s], drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylabel("State of charge",size = 24)
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
netpowerplot = [netpowerplot;netpowerplot[end,:]]

# Plot of netpower
n4 = [nrtm,ndam,ndays_planning,NS];
xplot = 0:dtrtm:dtrtm*nrtm_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
hold(true)
for s in S[1]
	plot(xplot,netpowerplot[:,s], drawstyle="steps-post", color = "blue", label="Net discharge");
end
grid()
xlim(0,nhours_planning)
ylabel("Net discharge (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
#legend(loc="upper right",fancybox="True", shadow="True", fontsize = 15)
savefig(string("cs719figures/netpower_fp_dt.pdf"))
close("all")












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
