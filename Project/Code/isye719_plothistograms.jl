using PyPlot

# PLOTTING HISTOGRAM OF ST COSTS
PyPlot.rc("ytick",labelsize=14)
PyPlot.rc("xtick",labelsize=14)
nbins = 20;
fig = figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
sub1 = subplot(4,1,1)
grid()
hold(true)
hist1 = plt[:hist](committedTotalCost_ST_VECTOR,nbins,weights=ones(length(committedTotalCost_ST_VECTOR))/length(committedTotalCost_ST_VECTOR),color="blue",label="Stochastic");
# PLOTTING HISTOGRAM OF DET COSTS
sub2 = subplot(3,1,2)
grid()
hold(true)
hist2 = plt[:hist](committedTotalCost_DET_VECTOR,nbins,weights=ones(length(committedTotalCost_DET_VECTOR))/length(committedTotalCost_DET_VECTOR),color="red",label="Deterministic");
sub3 = subplot(3,1,3)
grid()
hold(true)
hist3 = plt[:hist](committedTotalCost_PI_VECTOR,nbins,weights=ones(length(committedTotalCost_PI_VECTOR))/length(committedTotalCost_PI_VECTOR),color="green",label="Perfect Information");
# Setting xlims and ylims based on default lims of the two subplots
xlims = [collect(sub1[:get_xlim]()) collect(sub2[:get_xlim]()) collect(sub3[:get_xlim]())];
xlims = [minimum(xlims),maximum(xlims)];
ylims = [collect(sub1[:get_ylim]()) collect(sub2[:get_ylim]()) collect(sub3[:get_ylim]())];
ylims = [minimum(ylims),maximum(ylims)];
sub1[:set_xlim](xlims)
sub1[:set_ylim](ylims)
sub2[:set_xlim](xlims)
sub2[:set_ylim](ylims)
sub3[:set_xlim](xlims)
sub3[:set_ylim](ylims)
subplot(sub1)
plot(Expected_Total_Cost_ST*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Mean Cost") # Plot vertical line at expected value
ylabel("Probability",size = 14)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
subplot(sub2)
plot(Expected_Total_Cost_DET*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Mean Cost") # Plot vertical line at expected value
ylabel("Probability",size = 14)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
subplot(sub3)
plot(Expected_Total_Cost_PI*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Mean Cost") # Plot vertical line at expected value
ylabel("Probability",size = 14)
xlabel("Total Cost (\$/month)",size = 24)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
savefig(string(FiguresDirectory,"/histogram_costs_$(hor)day$(scaledf[scale]).pdf"))
close("all")
