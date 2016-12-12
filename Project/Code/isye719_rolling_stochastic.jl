using JuMP
using Gurobi
using PyPlot
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

function convertToArray2(A,n)
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

function convertToArray3(A,n)
    AA = getvalue(A)
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
rampmin = -0.5*P_max;	          #Lower bound for ramp discharge, MW/5min
rampmax = 0.5*P_max;  	  #Upper bound for ramp discharge, MW/5min
eff = 1;                  #Discharging Efficiency of battery
ebat0 = ebat_max;		   #Initial State of charge
ebatend = ebat_max;		  #State of charge at the end of the day
ndays = 365;              #Number of days data is available for
ndays_planning = 7;       #Number of days you want to plan for
ndays_horizon = 1;        #Number of days in the horizon at every step
nhours_planning = ndays_planning*ndam;
nrtm_planning = nhours_planning*nrtm;
nhours_horizon = ndays_horizon*ndam;
nrtm_horizon = nhours_horizon*nrtm;

NS = 50; # Number of scenarios you want to sample from the distrbution
S = 1:NS;




load = Matrix{Float64}(loaddata[2:nrtmpoints+1,2+(1:ndays)]);	#Load, MW
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

######################################################

load1 = loadNSplanningdata/1000;   #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);

=#


# Loading the NS scenarios for weekly load profiles in kW generated from the fullproblem_stochastic code
nweeks_planning = Int64(ceil((ndays_planning/reshape_ndays)));

loadNSdata = readcsv("loads_scenarios_month.csv")


ndays_data = (nweeks_planning+1)*reshape_ndays;
loadNSplanningdata = reshape(loadNSdata,nrtm,ndam,ndays_data,NS);   #kW

################################################################


load1 = loadNSplanningdata/1000; #MW
loaddata1 = reshape(load1,nrtm*ndam*ndays_data,NS);





dailyprofit = zeros(ndays);
loadplot = zeros(ndays*nrtmpoints);
rtmeprplot = zeros(ndays*nrtmpoints);
dameprplot = zeros(ndays*ndampoints);
ebatplot = zeros(ndays*nrtmpoints);
socplot = zeros(ndays*nrtmpoints);
Prtmplot = zeros(ndays*nrtmpoints);
Pdamplot = zeros(ndays*ndampoints);
rtmeprofitplot = zeros(ndays*nrtmpoints);
dameprofitplot = zeros(ndays*ndampoints);

totalpower = zeros(nrtm,ndam,ndays);
totalpowerplot = zeros(nrtm,ndam,ndays);



#=
realized_sequence = rand(S,nhours_planning);
writecsv("realized_sequence.csv",realized_sequence)
=#

# realized_sequence = readcsv("realized_sequence.csv")
# realized_sequence = Vector{Int64}(ones(nhours_planning));

m_rol = nothing;

obj_st_rh_NS = Vector()

for k in S # Loop to evaluate cost along each scenario

ebat0 = ebat_max;		  #Initial State of charge, 100 means fully charged


realized_sequence = Vector{Int64}(k*ones(nhours_planning));


Prtm_realized = zeros(nrtm,nhours_planning);
unmetload_realized = zeros(nrtm,nhours_planning)
unmetcost_realized = zeros(nhours_planning);
profitErtm_realized = zeros(nrtm,nhours_planning);
profitEdam_realized = zeros(nhours_planning);
profitE_realized = zeros(nhours_planning);
profittotal_realized = zeros(nhours_planning);
netobjective_realized = zeros(nhours_planning);

tic()
j=1;
for p in 1:nhours_planning
    println("Scenario $k, Step $p")

    #Load and price data
    load = loaddata1[(p-1)*nrtm+(1:nrtm_horizon),:];	#Load, MW

    eprrtm = rtmpricedata[(p-1)*nrtm+(1:nrtm_horizon),4];	    	#Real Time Market price, $/MWh
    eprdam = dampricedata[(p-1)+(1:nhours_horizon),4];	    	#Day Ahead Market Selling price, $/MWh

    #Reshape the data to matrices
    rtmepr = reshape(eprrtm,nrtm,ndam);
    damepr = reshape(eprdam,ndam);
    load = reshape(load,nrtm,ndam,NS);

    #Define sets to be used in the model defVar and addConstraint
    rtm = 1:nrtm;
    dam = 1:nhours_horizon;

    ################ Model ##################

    m_rol = Model(solver = GurobiSolver(Threads = 2, OutputFlag = 0))
    @variable(m_rol, -P_max <= Prtm[rtm,dam,S] <= P_max)	                #Power sold to the real time market, kW
    @variable(m_rol, -P_max <= Pdam[dam,S] <= P_max)    	                #Power sold to the day ahead market, kW
    @variable(m_rol, 0 <= ebat[rtm,dam,S] <= ebat_max)      	#Energy stored in the battery at the end of each real time interval, kWh
    @variable(m_rol, suppliedload[rtm,dam,S] >= 0)
    @variable(m_rol, unmetload[rtm,dam,S] >= 0)
		@expression(m_rol, Pnet[i in rtm,k in dam,s in S], Prtm[i,k,s] + Pdam[k,s] + suppliedload[i,k,s])    #Net power discharged from battery in all 5-min interval, kW
		@variable(m_rol, profitErtm[rtm,dam,S])# >= 0)				        #Profit from the real time market, USD
    @variable(m_rol, profitEdam[dam,S])# >= 0)	        			#Profit from the day ahead market, USD
    @variable(m_rol, profittotal[S])# >= 0)		        	#Total profit in the day, USD
    @variable(m_rol, unmetcost[S])

    @constraint(m_rol, InitialEnergy[s in S], ebat[1,1,s] == ebat0 - 1/eff*Pnet[1,1,s]*dtrtm)	#Inital energy in the battery
#    @constraint(m_rol, EndSOC[i in rtm,k in dam,s in S], soc[i,k,s] >= socend)		#Constraint on SOC at the end of the day
    @constraint(m_rol, rtmEBalance[i in rtm[2:end],k in dam,s in S], ebat[i,k,s] == ebat[i-1,k,s] - 1/eff*Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rol, damEBalance[i=rtm[1],k in dam[2:end],iend=rtm[end],s in S], ebat[i,k,s] == ebat[iend,k-1,s] - 1/eff*Pnet[i,k,s]*dtrtm)	#Dynamics constraint
    @constraint(m_rol, UnmetLoad[i in rtm,k in dam,s in S], suppliedload[i,k,s] + unmetload[i,k,s] >=  load[i,k,s])
    @constraint(m_rol, BoundSupplied[i in rtm,k in dam,s in S], suppliedload[i,k,s] <= load[i,k,s])
    @constraint(m_rol, BoundUnmet[i in rtm,k in dam,s in S], unmetload[i,k,s] <= load[i,k,s])
#=    @constraint(m_rol, RTMRamp1[i in rtm[2:end],k in dam,s in S], Pnet[i,k,s]  - Pnet[i-1,k,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rol, RTMRamp2[i in rtm[2:end],k in dam,s in S], Pnet[i,k,s]  - Pnet[i-1,k,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rol, DAMRamp1[i in rtm[1],k in dam[2:end],iend=rtm[end],s in S], Pnet[i,k,s] - Pnet[iend,k-1,s] <= rampmax*dtrtm)   #Ramp discharge constraint at each time
		@constraint(m_rol, DAMRamp2[i in rtm[1],k in dam[2:end],iend=rtm[end],s in S], Pnet[i,k,s] - Pnet[iend,k-1,s] >= -rampmax*dtrtm)   #Ramp discharge constraint at each time
=#
		@constraint(m_rol, RTMEProfits[i in rtm,k in dam,s in S], profitErtm[i,k,s] == rtmepr[i,k]*Prtm[i,k,s]*dtrtm)	#Economic calculation
    @constraint(m_rol, DAMEProfits[k in dam,s in S], profitEdam[k,s] == damepr[k]*Pdam[k,s]*dtdam)	#Economic calculation
    @constraint(m_rol, TotalProfit[s in S], profittotal[s] ==
                        sum{profitErtm[i,k,s], i in rtm, k in dam} + sum{profitEdam[k,s], k in dam})
    @constraint(m_rol, UnmetCost[s in S], unmetcost[s] == sum{rtmepr[i,k]*unmetload[i,k,s], i in rtm, k in dam})
    # Non-anticipativity constraints for first stage variables
    @constraint(m_rol, Nonant_PDAM[k in dam,s in S], Pdam[k,s] == (1/NS)*sum{Pdam[k,s], s in S})
    @objective(m_rol, Min, (1/NS)*sum{-profittotal[s] + unmetcost[s], s in S})
    #    print(m)
    status = solve(m)

##########################################################################

    ebat0 = getvalue(getvariable(m_rol,:ebat))[rtm[end],dam[1],realized_sequence[p]];
    Prtm_realized[:,p] = getvalue(getvariable(m_rol,:Prtm))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetload_realized[:,p] = getvalue(getvariable(m_rol,:unmetload))[1:rtm[end],dam[1],realized_sequence[p]];
    unmetcost_realized[p] = sum(unmetload_realized.*eprrtm[1:rtm[end]]);
    profitErtm_realized[:,p] = getvalue(getvariable(m_rol,:profitErtm))[1:rtm[end],dam[1],realized_sequence[p]];
    profitEdam_realized[p] = getvalue(getvariable(m_rol,:profitEdam))[dam[1],realized_sequence[p]];
    profittotal_realized[p] = sum(profitErtm_realized[:,p]) + profitEdam_realized[p];
    netobjective_realized[p] = unmetcost_realized[p] - profittotal_realized[p];

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

end # End rolling horizon


totalcost_after_rolling_st = sum(netobjective_realized);

obj_st_rh = totalcost_after_rolling_st;

push!(obj_st_rh_NS,obj_st_rh)

end

time_taken_st_rolling = toc();

expected_obj_st_rh = mean(obj_st_rh_NS);




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
