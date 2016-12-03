
using JuMP
using PyPlot
using Gurobi

close("all")

loaddata = readcsv("DATARCoast.csv");
rtmpricedata = readcsv("AggregatedData_RTM.csv")
qhmpricedata = readcsv("AggregatedData_HASP.csv")
dampricedata = readcsv("AggregatedData_DAM.csv")

dtdam = 1; 				#time interval of each day ahead market (hourly) intervals [hour]
ndam = 24;				#No. of day ahead market (hourly) intervals in a day
dtqhm = 15/60;				#time interval of each qhour [hour]
nqhm = Int64(dtdam/dtqhm);		#No. of quarter hours in an hour
dtrtm = 5/60;				#time interval of each real time interval [hour]
nrtm = Int64(dtqhm/dtrtm);		#No. of real time intervals in each quarter-hour [hour]

nrtmpoints = ndam*nqhm*nrtm;		#Total number of points in RTM data
nqhmpoints = ndam*nqhm;			#Total number of points in quarter hourly data
ndampoints = ndam;			#Total number of points in hourly data	

#Model Parameters
ebat_max = 500;	        #Battery capacity, kWh
ebuy_max = 1000;	#Maximum energy purchase rate, kW
esell_max = 1000;	#Minumum energy selling rate, kW
rampmin = -100;	        #Lower bound for ramp discharge, kWh
rampmax = 100;  	#Upper bound for ramp discharge, kWh
eff = 0.9;		#Efficiency of battery
soc0 = 100;		#Initial State of charge, 100 means fully charged
socend = 100;		#State of charge at the end of the day

#Load and price data
load = loaddata[2:nrtmpoints+1,3];	#Load, kW
#= For day 2 use this
sellprrtm = 0.001*rtmpricedata[nrtmpoints+[1:nrtmpoints],4];	    	#Real Time Market Selling price, $/kWh (currently taking only day-1 price data)
buyprrtm = 0.001*0.5*rtmpricedata[nrtmpoints+[1:nrtmpoints],4]; 	#Real Time Market Purchase price, $/kWh (currently taking only day-1 price data)

sellprqhm = 0.001*qhmpricedata[nqhmpoints+[1:nqhmpoints],4];	    	#Quarter Hourly Market Selling price, $/kWh (currently taking only day-1 price data)
buyprqhm = 0.001*0.5*qhmpricedata[nqhmpoints+[1:nqhmpoints],4];    	#Quarter Hourly Market Purchase price, $/kWh (currently taking only day-1 price data)

sellprdam = 0.001*dampricedata[ndampoints+[1:ndampoints],4];	    	#Day Ahead Market Selling price, $/kWh (currently taking only day-1 price data)
buyprdam = 0.001*0.5*dampricedata[ndampoints+[1:ndampoints],4]; 	#Day Ahead Market Purchase price, $/kWh (currently taking only day-1 price data)
=#
sellprrtm = 0.001*rtmpricedata[[1:nrtmpoints],4];	    	#Real Time Market Selling price, $/kWh (currently taking only day-1 price data)
buyprrtm = 0.001*0.5*rtmpricedata[[1:nrtmpoints],4]; 		#Real Time Market Purchase price, $/kWh (currently taking only day-1 price data)
buyprrtm = buyprrtm[end:-1:1];

sellprqhm = 0.001*qhmpricedata[[1:nqhmpoints],4];	    	#Quarter Hourly Market Selling price, $/kWh (currently taking only day-1 price data)
buyprqhm = 0.001*0.5*qhmpricedata[[1:nqhmpoints],4];    	#Quarter Hourly Market Purchase price, $/kWh (currently taking only day-1 price data)
buyprqhm = buyprqhm[end:-1:1];

sellprdam = 0.001*dampricedata[[1:ndampoints],4];	    	#Day Ahead Market Selling price, $/kWh (currently taking only day-1 price data)
buyprdam = 0.001*0.5*dampricedata[[1:ndampoints],4]; 		#Day Ahead Market Purchase price, $/kWh (currently taking only day-1 price data)
buyprdam = buyprdam[end:-1:1];

regupprrtm = 0.001*rtmpricedata[[1:nrtmpoints],6];	    	#Real Time Market Selling price, $/kWh (currently taking only day-1 price data)
regdownprrtm = 0.001*rtmpricedata[[1:nrtmpoints],7]; 		#Real Time Market Purchase price, $/kWh (currently taking only day-1 price data)

regupprqhm = 0.001*qhmpricedata[[1:nqhmpoints],6];	    	#Quarter Hourly Market Selling price, $/kWh (currently taking only day-1 price data)
regdownprqhm = 0.001*qhmpricedata[[1:nqhmpoints],7];    	#Quarter Hourly Market Purchase price, $/kWh (currently taking only day-1 price data)

regupprdam = 0.001*dampricedata[[1:ndampoints],6];	    	#Day Ahead Market Selling price, $/kWh (currently taking only day-1 price data)
regdownprdam = 0.001*dampricedata[[1:ndampoints],7]; 		#Day Ahead Market Purchase price, $/kWh (currently taking only day-1 price data)


#Reshape the data to matrices
rtmsellpr = reshape(sellprrtm,nrtm,nqhm,ndam);
rtmbuypr = reshape(buyprrtm,nrtm,nqhm,ndam);

qhmsellpr = reshape(sellprqhm,nqhm,ndam);
qhmbuypr = reshape(buyprqhm,nqhm,ndam);

damsellpr = reshape(sellprdam,ndam);
dambuypr = reshape(buyprdam,ndam);

damreguppr = reshape(regupprdam,ndam);
damregdownpr = reshape(regdownprdam,ndam);

rtmreguppr = reshape(regupprrtm,nrtm,nqhm,ndam);
rtmregdownpr = reshape(regdownprrtm,nrtm,nqhm,ndam);

qhmreguppr = reshape(regupprqhm,nqhm,ndam);
qhmregdownpr = reshape(regdownprqhm,nqhm,ndam);

load = reshape(load,nrtm,nqhm,ndam);

loadmean = 300*ones(nrtm,nqhm,ndam);    #Mean load

#Define sets to be used in the model defVar and addConstraint
rtm = 1:nrtm;
qhm = 1:nqhm;
dam = 1:ndam;

################ Model ##################

m = Model(solver = GurobiSolver())

    @defVar(m, 0 <= ebuyrtm[rtm,qhm,dam] <= ebuy_max)	#Energy bought from the real time market, kW
    @defVar(m, 0 <= ebuyqhm[qhm,dam] <= ebuy_max)	#Energy bought from the quarter hourly market, kW
    @defVar(m, 0 <= ebuydam[dam] <= ebuy_max)		#Energy bought from the day ahead market, kW
    @defVar(m, 0 <= esellrtm[rtm,qhm,dam] <= esell_max)	#Energy bought from the real time market, kW
    @defVar(m, 0 <= esellqhm[qhm,dam] <= esell_max)	#Energy bought from the real time market, kW
    @defVar(m, 0 <= eselldam[dam] <= esell_max)    	#Energy bought from the real time market, kW
    @defVar(m, 0 <= ebat[rtm,qhm,dam] <= ebat_max)	#Energy stored in the battery at the end of each real time interval
    @defVar(m, 0 <= soc[rtm,qhm,dam] <= 100)		#SOC of the battery at the end of each real time interval      
    @defVar(m, regup[rtm,qhm,dam] >= 0)                 #Amount of regulation up, kW
    @defVar(m, regdown[rtm,qhm,dam] >= 0)               #Amount of regulation down, kW
    @defVar(m, profitrtm >= 0)				#Profit from the real time market, USD
    @defVar(m, profitqhm >= 0)				#Profit from the quarter hourly market, USD
    @defVar(m, profitdam >= 0)				#Profit from the day ahead market, USD		
    @defVar(m, profitreguprtm[rtm,qhm,dam] >= 0)			#Profit from the real time market, USD
    @defVar(m, profitregupqhm[qhm,dam] >= 0)			#Profit from the quarter hourly market, USD
    @defVar(m, profitregupdam[dam] >= 0)			#Profit from the day ahead market, USD
    @defVar(m, profitregdownrtm[rtm,qhm,dam] >= 0)			#Profit from the real time market, USD
    @defVar(m, profitregdownqhm[qhm,dam] >= 0)			#Profit from the quarter hourly market, USD
    @defVar(m, profitregdowndam[dam] >= 0)			#Profit from the day ahead market, USD
    @defVar(m, profittotal >= 0)			#Total profit in the day, USD		
#    @defVar(m, minimum(load[rtm,qhm,dam]) <= loadmean[rtm,qhm,dam] <= maximum(load[rtm,qhm,dam]))

    @addConstraint(m, InitialEnergy, ebat[1,1,1] == soc0/100*ebat_max + (eff*(ebuyrtm[1,1,1] + ebuyqhm[1,1] + ebuydam[1]) - (esellrtm[1,1,1] + esellqhm[1,1] + eselldam[1]) - load[1,1,1])*dtrtm)	#Inital energy in the battery
        
    @addConstraint(m, DefSOC[i=rtm,j=qhm,k=dam], soc[i,j,k] == ebat[i,j,k]/ebat_max*100)			#Define SOC

    @addConstraint(m, EndSOC[i=rtm[end],j=qhm[end],k=dam[end]], soc[i,j,k] >= socend)		#Constraint on SOC at the end of the day
    
    @addConstraint(m, rtmEBalance[i=rtm[2:end],j=qhm,k=dam], ebat[i,j,k] == ebat[i-1,j,k] + (eff*(ebuyrtm[i,j,k] + ebuyqhm[j,k] + ebuydam[k]) - (esellrtm[i,j,k] + esellqhm[j,k] + eselldam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint
    
    @addConstraint(m, qhmEBalance[i=rtm[1],j=qhm[2:end],k=dam,iend=rtm[end]], ebat[i,j,k] == ebat[iend,j-1,k] + (eff*(ebuyrtm[i,j,k] + ebuyqhm[j,k] + ebuydam[k]) - (esellrtm[i,j,k] + esellqhm[j,k] + eselldam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint
    
    @addConstraint(m, damEBalance[i=rtm[1],j=qhm[1],k=dam[2:end],iend=rtm[end],jend=qhm[end]], ebat[i,j,k] == ebat[iend,jend,k-1] + (eff*(ebuyrtm[i,j,k] + ebuyqhm[j,k] + ebuydam[k]) - (esellrtm[i,j,k] + esellqhm[j,k] + eselldam[k]) - load[i,j,k])*dtrtm)	#Dynamics constraint

    @addConstraint(m, RegUp[i=rtm,j=qhm,k=dam], regup[i,j,k] == max((load[i,j,k] - loadmean[i,j,k]),0))
    @addConstraint(m, RegDown[i=rtm,j=qhm,k=dam], regdown[i,j,k] == max(-(load[i,j,k] - loadmean[i,j,k]),0))

   # @addConstraint(m, RTMRamp[i=rtm[2:end],j=qhm,k=dam], rampmin*dtrtm <= ebat[i,j,k] - ebat[i-1,j,k] <= rampmax*dtrtm)   #Ramp discharge constraint at each time

   # @addConstraint(m, QHMRamp[i=rtm[1],j=qhm[2:end],k=dam,iend=rtm[end]], rampmin*dtrtm <= ebat[i,j,k] - ebat[iend,j-1,k] <= rampmax*dtrtm)   #Ramp charging rate constraint at each time
    
   # @addConstraint(m, DAMRamp[i=rtm[1],j=qhm[1],k=dam[2:end],iend=rtm[end],jend=qhm[end]], rampmin*dtrtm <= ebat[i,j,k] - ebat[iend,jend,k-1] <= rampmax*dtrtm)   #Ramp discharge rate constraint at each time
    
    @addConstraint(m, RTMProfits, profitrtm == sum{(rtmsellpr[i,j,k]*esellrtm[i,j,k] - rtmbuypr[i,j,k]*ebuyrtm[i,j,k])*dtrtm,i=rtm,j=qhm,k=dam})	#Economic calculation
    
    @addConstraint(m, QHMProfits, profitqhm == sum{(qhmsellpr[j,k]*esellqhm[j,k] - qhmbuypr[j,k]*ebuyqhm[j,k])*dtqhm,j=qhm,k=dam})	#Economic calculation
    
    @addConstraint(m, DAMProfits, profitdam == sum{(damsellpr[k]*eselldam[k] - dambuypr[k]*ebuydam[k])*dtdam,k=dam})	#Economic calculation
  
    @addConstraint(m, RTMregupProfits[i=rtm,j=qhm,k=dam], profitreguprtm[i,j,k] == rtmreguppr[i,j,k]*dtrtm*regup[i,j,k])
    @addConstraint(m, QHMregupProfits[j=qhm,k=dam], profitregupqhm[j,k] == qhmreguppr[j,k]*sum{dtrtm*regup[i,j,k],i = rtm})
    @addConstraint(m, DAMregupProfits[k=dam], profitregupdam[k] == damreguppr[k]*sum{dtrtm*regup[i,j,k], i=rtm,j=qhm})

    @addConstraint(m, RTMregdownProfits[i=rtm,j=qhm,k=dam], profitregdownrtm[i,j,k] == rtmregdownpr[i,j,k]*dtrtm*regdown[i,j,k])
    @addConstraint(m, QHMregdownProfits[j=qhm,k=dam], profitregdownqhm[j,k] == qhmregdownpr[j,k]*sum{dtrtm*regdown[i,j,k], i=rtm})
    @addConstraint(m, DAMregdownProfits[k=dam], profitregdowndam[k] == damregdownpr[k]*sum{dtrtm*regdown[i,j,k], i=rtm,j=qhm})
     
    @addConstraint(m, TotalProfit, profittotal == profitrtm + profitqhm + profitdam + sum(profitreguprtm) + sum(profitregupqhm) + sum(profitregupdam) + sum(profitregdownrtm) + sum(profitregdownqhm) + sum(profitregdowndam))
    
    @setObjective(m, Max, profittotal);
    
    #print(m)

    status = solve(m)

    #println("EBuyRTM = \n", getValue(ebuyrtm))
    #println("ESellRTM = \n", getValue(esellrtm))
    #println("EBuyQHM = \n", getValue(ebuyqhm))
    #println("ESellQHM = \n", getValue(esellqhm))
    #println("EBuyDAM = \n", getValue(ebuydam))
    #println("ESellDAM = \n", getValue(eselldam))
    
    #println("EBat = \n", getValue(ebat))
    println("\nObjective value: ", getObjectiveValue(m),"\n")


############# Funtions to convert JuMP returned variables to arrays ################

function convertToArray(x)
	y = getValue(x)
	n = length(y)
	a = zeros(n)
	for i = 1:n
		a[i] = y[i]
	end
	return a
end

function convertToArray2(A)	
	AA = getValue(A)	
	n = [4,24]
	m = (n[1],n[2])
	B = zeros(m)	
	for i in 1:n[1]
            for j in 1:n[2]
		B[i,j] = AA[i,j]
	    end
	end	
	return B
end

function convertToArray3(A)
	AA = getValue(A)	
	n = [3,4,24]
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

function qhmtortm(A) 
    A = A[1:nqhmpoints]
    B = zeros(nrtmpoints);
    for i in 1:nqhmpoints
        B[(i-1)*3+1:i*3] = repmat([A[i]],3);
    end
    B = push!(B,B[end])
    return B
end    

function damtortm(A) 
    A = A[1:ndampoints]
    B = zeros(nrtmpoints);
    for i in 1:ndampoints
        B[(i-1)*12+1:i*12] = repmat([A[i]],12);
    end
    B = push!(B,B[end])
    return B
end    


#################### Plotting #####################

timertm = 0:dtrtm:ndam;
timeqhm = 0:dtqhm:ndam;
timedam = 0:dtdam:ndam;

loadsplot = reshape(load,nrtmpoints)/1000;
loadsplot = push!(loadsplot,loadsplot[end]);

rtmsellprplot = reshape(rtmsellpr,nrtmpoints)*1000;
rtmsellprplot = push!(rtmsellprplot,rtmsellprplot[end]);
rtmbuyprplot = reshape(rtmbuypr,nrtmpoints)*1000;
rtmbuyprplot = push!(rtmbuyprplot,rtmbuyprplot[end]);

qhmsellprplot = reshape(qhmsellpr,nqhmpoints)*1000;
qhmsellprplot = push!(qhmsellprplot,qhmsellprplot[end]);
qhmbuyprplot = reshape(qhmbuypr,nqhmpoints)*1000;
qhmbuyprplot = push!(qhmbuyprplot,qhmbuyprplot[end]);

damsellprplot = reshape(damsellpr,ndampoints)*1000;
damsellprplot = push!(damsellprplot,damsellprplot[end]);
dambuyprplot = reshape(dambuypr,ndampoints)*1000;
dambuyprplot = push!(dambuyprplot,dambuyprplot[end]);

qhmregupprplot = reshape(qhmreguppr,nqhmpoints)*1000;
qhmregupprplot = push!(qhmregupprplot,qhmregupprplot[end]);

damregupprplot = reshape(damreguppr,ndampoints)*1000;
damregupprplot = push!(damregupprplot,damregupprplot[end]);

qhmregdownprplot = reshape(qhmregdownpr,nqhmpoints)*1000;
qhmregdownprplot = push!(qhmregdownprplot,qhmregdownprplot[end]);

damregdownprplot = reshape(damregdownpr,ndampoints)*1000;
damregdownprplot = push!(damregdownprplot,damregdownprplot[end]);

ebatplot = reshape(convertToArray3(ebat),nrtmpoints);
#ebatplot = push!(ebatplot,ebatplot[end]);
ebatplot = [ebat_max,ebatplot]/1000;

socplot = reshape(convertToArray3(soc),nrtmpoints);
#socplot = push!(socplot,socplot[end]);
socplot = [soc0,socplot];

esellrtmplot = reshape(convertToArray3(esellrtm),nrtmpoints);
esellrtmplot = push!(esellrtmplot,esellrtmplot[end])/1000;

ebuyrtmplot = reshape(convertToArray3(ebuyrtm),nrtmpoints);
ebuyrtmplot = push!(ebuyrtmplot,ebuyrtmplot[end])/1000;

esellqhmplot = reshape(convertToArray2(esellqhm),nqhmpoints);
esellqhmplot = push!(esellqhmplot,esellqhmplot[end])/1000;

ebuyqhmplot = reshape(convertToArray2(ebuyqhm),nqhmpoints);
ebuyqhmplot = push!(ebuyqhmplot,ebuyqhmplot[end])/1000;

eselldamplot = reshape(convertToArray(eselldam),ndampoints);
eselldamplot = push!(eselldamplot,eselldamplot[end])/1000;

ebuydamplot = reshape(convertToArray(ebuydam),ndampoints);
ebuydamplot = push!(ebuydamplot,ebuydamplot[end])/1000;

esellrtmplot = reshape(convertToArray3(esellrtm),nrtmpoints);
esellrtmplot = push!(esellrtmplot,esellrtmplot[end])/1000;

rtmregprofitplot = reshape(convertToArray3(profitreguprtm),nrtmpoints) + reshape(convertToArray3(profitregdownrtm),nrtmpoints);
rtmregprofitplot = push!(rtmregprofitplot,rtmregprofitplot[end]);
qhmregprofitplot = reshape(convertToArray2(profitregupqhm),nqhmpoints) + reshape(convertToArray2(profitregdownqhm),nqhmpoints);
qhmregprofitplot = push!(qhmregprofitplot,qhmregprofitplot[end]);
damregprofitplot = reshape(convertToArray(profitregupdam),ndampoints) + reshape(convertToArray(profitregdowndam),ndampoints);
damregprofitplot = push!(damregprofitplot,damregprofitplot[end]);

rtmprofitplot = rtmsellprplot.*esellrtmplot*dtrtm - rtmbuyprplot.*ebuyrtmplot*dtrtm + rtmregprofitplot;
qhmprofitplot = qhmsellprplot.*esellqhmplot*dtqhm - qhmbuyprplot.*ebuyqhmplot*dtqhm + qhmregprofitplot;
damprofitplot = damsellprplot.*eselldamplot*dtdam - dambuyprplot.*ebuydamplot*dtdam + damregprofitplot;

figure()
plot(timertm,loadsplot,drawstyle="steps-post",label="Load (MW)")
hold(true)
plot(timertm,[vec(loadmean),loadmean[end]]/1000)
xlabel("Time (hr)",size=14), ylabel("Load (MW)",size=14)
xlim(0,24)
grid()
#title("Load (kW)")
#legend(loc=2,fancybox="True", shadow="True")

figure()
subplot(211)
plot(timertm,rtmsellprplot,drawstyle="steps-post",color="blue",label="RTM Energy selling price (\$/MWh)")
hold(true)
plot(timertm,rtmbuyprplot,drawstyle="steps-post",color="green",label="RTM Energy purchase price (\$/MWh)")
xlabel("Time (hr)",size=14)
ylabel("RTM Energy prices (\$/MWh)",size=14)
xlim(0,24)
grid()

#title("RTM Prices of energy purchase and selling")
#legend(loc=2,fancybox="True", shadow="True")
subplot(212)
plot(timertm,esellrtmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in RTM")
ylim(-0.100,1.1*esell_max/1000)
xlim(0,24)
hold(true)
plot(timertm,ebuyrtmplot,drawstyle="steps-post",color="green",label="Power purchased (kW) in RTM")
ylim(-0.100,1.1*ebuy_max/1000)
xlabel("Time (hr)",size=14), ylabel("Power (MW)",size=14)
xlim(0,24)
grid()
#title("Rate of energy purchase and selling in RTM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timeqhm,qhmsellprplot,drawstyle="steps-post",color="blue",label="QHM Energy selling price (\$/kWh)")
hold(true)
plot(timeqhm,qhmbuyprplot,drawstyle="steps-post",color="green",label="QHM Energy purchase price (\$/kWh)")
xlabel("Time (hr)",size=20)
ylabel("QHM Energy prices (\$/MWh)",size=20)
xlim(0,24)
grid()
#title("QHM Prices of energy purchase and selling")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timeqhm,esellqhmplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in QHM")
ylim(-0.100,1.1*esell_max/1000)
hold(true)
plot(timeqhm,ebuyqhmplot,drawstyle="steps-post",color="green",label="Power purchased (MW) in QHM")
ylim(-0.100,1.1*ebuy_max/1000)
xlabel("Time (hr)",size=14), ylabel("Power (MW)",size=14)
xlim(0,24)
grid()
#title("Rate of energy purchase and selling in QHM")
#legend(loc=2,fancybox="True", shadow="True")


figure()
subplot(211)
plot(timedam,damsellprplot,drawstyle="steps-post",color="blue",label="DAM Energy selling price (\$/kWh)")
hold(true)
plot(timedam,dambuyprplot,drawstyle="steps-post",color="green",label="DAM Energy purchase price (\$/kWh)")
xlabel("Time (hr)",size=20)
ylabel("DAM Energy prices (\$/MWh)",size=20)
xlim(0,24)
grid()
#title("DAM Prices of energy purchase and selling")
#legend(loc=2,fancybox="True", shadow="True")

subplot(212)
plot(timedam,eselldamplot,drawstyle="steps-post",color="blue",label="Power sold (MW) in DAM")
ylim(-0.100,1.1*esell_max/1000)
hold(true)
plot(timedam,ebuydamplot,drawstyle="steps-post",color="green",label="Power purchased (MW) in DAM")
ylim(-0.100,1.1*ebuy_max/1000)
xlabel("Time (hr)",size=14), ylabel("Power (MW)",size=14)
xlim(0,24)
grid()
#title("Rate of energy purchase and selling in DAM")
#legend(loc=2,fancybox="True", shadow="True")


#figure()
#subplot(211)
#plot(timertm,ebatplot[1:length(ebatplot)],drawstyle="steps-post",label="Energy of battery (kWh)")
#xlabel("Time (hr)")
#ylabel("Energy (kWh)")
#xlim(0,24)
#title("Energy of Battery (kWh)")

figure()
plot(timertm,socplot,drawstyle="steps-post",label="SOC")
xlabel("Time (hr)",size=20), ylabel("SOC (%)",size=20)
xlim(0,24)
grid()
#title("SOC")
#legend(loc=2,fancybox="True", shadow="True")

function cumul(A)
    C = zeros(length(A)-1);
    for i in 1:length(A)-1
        C[i] = sum(A[1:i]);
    end
    C = push!(C,C[end]);
    return C;
end
rtmprofitplot = cumul(rtmprofitplot);
qhmprofitplot = cumul(qhmprofitplot);
damprofitplot = cumul(damprofitplot);

totalprofitplot = rtmprofitplot + qhmtortm(qhmprofitplot) + damtortm(damprofitplot);

figure()
plot(timertm,rtmprofitplot,drawstyle="steps-post",label="RTM Profits")
xlabel("Time (hr)",size=14), ylabel("RTM Revenue (\$)",size=14)
xlim(0,24)
grid()

figure()
plot(timeqhm,qhmprofitplot,drawstyle="steps-post",label="QHM Profits")
xlabel("Time (hr)",size=14), ylabel("QHM Revenue (\$)",size=14)
xlim(0,24)
grid()

figure()
plot(timedam,damprofitplot,drawstyle="steps-post",label="DAM Profits")
xlabel("Time (hr)",size=14), ylabel("DAM Revenue (\$)",size=14)
xlim(0,24)
grid()

figure()
plot(timertm,totalprofitplot,drawstyle="steps-post",label="Total Profits")
xlabel("Time (hr)",size=14), ylabel("Total Revenue (\$)",size=14)
xlim(0,24)
grid()

figure()
plot(timeqhm,qhmregupprplot,drawstyle="steps-post",color="blue",label="Reg Up")
hold(true)
plot(timeqhm,qhmregdownprplot,drawstyle="steps-post",color="green",label="Reg Down")
xlabel("Time (hr)",size=20)
ylabel("QHM Regulation prices (\$/MW)",size=20)
xlim(0,24)
grid()
#title("QHM Prices of energy purchase and selling")
legend(loc=2,fancybox="True", shadow="True")


figure()
plot(timedam,damregupprplot,drawstyle="steps-post",color="blue",label="Reg Up")
hold(true)
plot(timedam,damregdownprplot,drawstyle="steps-post",color="green",label="Reg Down")
xlabel("Time (hr)",size=20)
ylabel("DAM Regulation prices (\$/MW)",size=20)
xlim(0,24)
grid()
#title("DAM Prices of energy purchase and selling")
legend(loc=2,fancybox="True", shadow="True")
