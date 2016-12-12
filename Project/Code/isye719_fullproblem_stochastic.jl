
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
regup_max = 0.5*P_max;    #Regulation Up Capacity, MW
regdown_max = 0.5*P_max;  #Regulation Up Capacity, MW
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


#=

# Using Ledoit-Wolfe Sample Covariance Estimator
(p,n) = size(load_weekly);
load_weeklymean = mean(load_weekly,2);
X = load_weekly-repmat(load_weeklymean,1,n); # Columns are 168*1 random vectors with mean 0 and covariance Sigma
Sn = X*X'/n;
mn = trace(Sn*eye(size(Sn)[1])')/p;
dn2 = trace((Sn-mn*eye(size(Sn)[1]))*(Sn-mn*eye(size(Sn)[1]))')/p;
bnbar2 = 0;
for k = 1:n
    bnbar2 = bnbar2 + trace((X[:,k]*X[:,k]'-Sn)*(X[:,k]*X[:,k]'-Sn)')/p;
end
bnbar2 = bnbar2/n^2;
bn2 = min(bnbar2,dn2);
an2 = dn2 - bn2;
Snstar = bn2/dn2*eye(p) + an2/dn2*Sn; # Estimator of variance Sigma


# Generating NS scenarios for weekly load profiles in kW
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));
loadNSdata = zeros(reshape_nrows*(nweeks_planning+1),NS);
R = chol(Snstar);
for j in 1:nweeks_planning+1
    loadNSdata[(j-1)*reshape_nrows+(1:reshape_nrows),:]  = repmat(load_weeklymean,1,NS) + R'*randn(p,NS);
end
loadNSdata[loadNSdata.<=0] = minimum(loadvec);
writecsv("loads_scenarios_month.csv",loadNSdata)


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

################################################################


load = loadNSplanningdata[:,:,1:ndays_planning,:]/1000; #MW

=#


# Loading the NS scenarios for weekly load profiles in kW generated from the fullproblem_stochastic code
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));

loadNSdata = readcsv("loads_scenarios_month.csv")


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

################################################################


load = loadNSplanningdata[:,:,1:ndays_planning,:]/1000; #MW



#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
dam = 1:ndam;
day = 1:ndays_planning;

################ Model ##################
tic()
m = Model(solver = GurobiSolver(Threads=2))

    @variable(m, -P_max <= Prtm[rtm,dam,day,S] <= P_max)	                #Power sold to the real time market, kW
    @variable(m, -P_max <= Pdam[dam,day,S] <= P_max)    	                #Power sold to the day ahead market, kW
    @variable(m, 0 <= ebat[rtm,dam,day,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
    @variable(m, 0 <= regupdam[dam,day,S] <= regup_max)                 #Amount of regulation up, kW
    @variable(m, 0 <= regdowndam[dam,day,S] <= regdown_max)             #Amount of regulation down, kW
    @variable(m, suppliedload[rtm,dam,day,S] >= 0)
    @variable(m, unmetload[rtm,dam,day,S] >= 0)
		@expression(m, Pnet[i in rtm,k in dam,l in day,s in S], Prtm[i,k,l,s] + Pdam[k,l,s] + suppliedload[i,k,l,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m, profitErtm[rtm,dam,day,S])# >= 0)				        #Profit from the real time market, USD
    @variable(m, profitEdam[dam,day,S])# >= 0)	        			#Profit from the day ahead market, USD
    @variable(m, profitregupdam[dam,day,S])# >= 0)			        #Profit from the day ahead market, USD
    @variable(m, profitregdowndam[dam,day,S])# >= 0)	        		#Profit from the day ahead market, USD
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

    @constraint(m, RTMRamp1[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, RTMRamp2[i in rtm[2:end],k in dam,l in day,s in S], Pnet[i,k,l,s]  - Pnet[i-1,k,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time

		@constraint(m, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],l in day,s in S], Pnet[i,k,l,s] - Pnet[iend,k-1,l,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time

		@constraint(m, DAYRamp1[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m, DAYRamp2[i=rtm[1],k=dam[1],iend=rtm[end],kend=dam[end],l in day[2:end],s in S], Pnet[i,k,l,s] - Pnet[iend,kend,l-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time

		@constraint(m, RegUp[i in rtm,k in dam,l in day,s in S], Pnet[i,k,l,s] + regupdam[k,l,s] <= P_max)	#Constraint on total power

    @constraint(m, RegDown[i in rtm,k in dam,l in day,s in S], Pnet[i,k,l,s] - regdowndam[k,l,s] >= -P_max)	#Constraint on total power

		@constraint(m, EnsureRegUp[i in rtm,k in dam,l in day,s in S], ebat[i,k,l,s] - dtrtm*(Pnet[i,k,l,s] + regupdam[k,l,s]) >= 0)
		@constraint(m, EnsureRegDown[i in rtm,k in dam,l in day,s in S], ebat[i,k,l,s] - dtrtm*(Pnet[i,k,l,s] - regupdam[k,l,s]) <= ebat_max)

    @constraint(m, RTMEProfits[i in rtm,k in dam,l in day,s in S], profitErtm[i,k,l,s] == rtmepr[i,k,l]*Prtm[i,k,l,s]*dtrtm)	#Economic calculation
    @constraint(m, DAMEProfits[k in dam,l in day,s in S], profitEdam[k,l,s] == damepr[k,l]*Pdam[k,l,s]*dtdam)	#Economic calculation

    @constraint(m, DAMregupProfits[k in dam,l in day,s in S], profitregupdam[k,l,s] == damreguppr[k,l]*regupdam[k,l,s])
    @constraint(m, DAMregdownProfits[k in dam,l in day,s in S], profitregdowndam[k,l,s] == damregdownpr[k,l]*regdowndam[k,l,s])

    @constraint(m, TotalProfit[s in S], profittotal[s] ==
                        sum{profitErtm[i,k,l,s], i in rtm, k in dam, l in day} + sum{profitEdam[k,l,s], k in dam, l in day}
                        + sum{profitregupdam[k,l,s], k in dam, l in day} + sum{profitregdowndam[k,l,s], k in dam, l in day})

    @constraint(m, UnmetCost[s in S], unmetcost[s] == sum{rtmepr[i,k,l]*unmetload[i,k,l,s], i in rtm, k in dam, l in day})


    # Non-anticipativity constraints for first stage variables
    @constraint(m, Nonant_PDAM[k in dam,l in day,s in S], Pdam[k,l,s] == (1/NS)*sum{Pdam[k,l,s], s in S})
    @constraint(m, Nonant_DAMregup[k in dam,l in day,s in S], regupdam[k,l,s] == (1/NS)*sum{regupdam[k,l,s], s in S})
    @constraint(m, Nonant_DAMregdown[k in dam,l in day,s in S], regdowndam[k,l,s] == (1/NS)*sum{regdowndam[k,l,s], s in S})

    @objective(m, Min, (1/NS)*sum{-profittotal[s] + unmetcost[s], s in S})

#    print(m)

    status = solve(m)

time_taken_st_fullproblem = toc();

###############################################################

obj_st_fp = getobjectivevalue(m);

#    println("\nTotal Profits ", getvalue(profittotal),"\n" )


################# PLOTTING #################
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
savefig(string("cs719figures/loads_scenarios.pdf"))
close("all")

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
subplot(2,1,2)
hold(true)
plot(xplot,Pdamplot, color="blue", drawstyle="steps-post");
grid()
xlim(0,nhours_planning)
ylabel("Net Power (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
savefig(string("cs719figures/Pdam_fp_st.pdf"))
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
savefig(string("cs719figures/Prtm_fp_st.pdf"))
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
	plot(xplot,socplot[:,s],  drawstyle="steps-post");
end
grid()
xlim(0,nhours_planning)
ylabel("State of charge",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
savefig(string("cs719figures/soc_fp_st.pdf"))
close("all")

# Plot of regup and regdown
n3 = [ndam,ndays_planning,NS];
regupprdamplot = regupprdam[1:nhours_planning];
regdownprdamplot = regdownprdam[1:nhours_planning];
push!(regupprdamplot,regupprdamplot[end]);
push!(regdownprdamplot,regdownprdamplot[end]);
regupdamarray = convertToArray3(getvalue(getvariable(m,:regupdam)),n3);
regupdamplot = reshape(regupdamarray[:,:,1],nhours_planning);
push!(regupdamplot,regupdamplot[end]);
regdowndamarray = convertToArray3(getvalue(getvariable(m,:regdowndam)),n3);
regdowndamplot = reshape(regdowndamarray[:,:,1],nhours_planning);
push!(regdowndamplot,regdowndamplot[end]);
xplot = 0:dtdam:dtdam*nhours_planning;
figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
subplot(2,1,1)
hold(true)
plot(xplot,regupprdamplot, color="blue", drawstyle="steps-post", label="Reg up");
plot(xplot,regdownprdamplot, color="red", drawstyle="steps-post", label="Reg down");
grid()
xlim(0,nhours_planning)
ylabel("Regulation prices (\$/kWh)",size = 24)
tick_params(labelsize=14)
subplot(2,1,2)
hold(true)
plot(xplot,regupdamplot, color="blue", drawstyle="steps-post", label="Reg up");
plot(xplot,regupdamplot, color="blue", drawstyle="steps-post", label="Reg down");
grid()
xlim(0,nhours_planning)
ylabel("Regulation capacity (kW)",size = 24)
xlabel("Time (hours)",size = 24)
tick_params(labelsize=14)
savefig(string("cs719figures/reg_fp_st.pdf"))
close("all")

# Pnet, regulation bands calculation
netpower = zeros(nrtm,ndam,ndays_planning,NS);
upband = zeros(nrtm,ndam,ndays_planning,NS);
downband = zeros(nrtm,ndam,ndays_planning,NS);
for i in rtm
        for k in dam
            for l in day
							for s in S
                netpower[i,k,l,s] = Prtmarray[i,k,l,s] + Pdamarray[k,l,s];
                upband[i,k,l,s] = netpower[i,k,l,s] + regupdamarray[k,l,s];
                downband[i,k,l,s] = netpower[i,k,l,s] - regdowndamarray[k,l,s];
							end
            end
        end
end
netpowerplot = reshape(netpower,nrtm_planning,NS);
upbandplot = reshape(upband,nrtm_planning,NS);
downbandplot = reshape(downband,nrtm_planning,NS);
netpowerplot = [netpowerplot;netpowerplot[end,:]]
upbandplot = [upbandplot;upbandplot[end,:]]
downbandplot = [downbandplot;downbandplot[end,:]]

#=

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

=#



#= #Start commenting here

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
