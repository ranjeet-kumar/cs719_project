using PyPlot

obj_st_fp_all = -obj_st_fp_all;
obj_dt_fp_all = -obj_dt_fp_all;
obj_pi_fp_all = -obj_pi_fp_all;
obj_st_fp_socfix_all = -obj_st_fp_socfix_all;
# PLOTTING HISTOGRAM OF ST COSTS
PyPlot.rc("ytick",labelsize=14)
PyPlot.rc("xtick",labelsize=14)
nbins = 20;
fig = figure()
plt[:get_current_fig_manager]()[:full_screen_toggle]()
sub1 = subplot(4,1,1)
grid()
hold(true)
hist1 = plt[:hist](obj_dt_fp_all,nbins*2,weights=ones(length(obj_dt_fp_all))/length(obj_dt_fp_all),color="indigo",label="Mean-value problem");
# PLOTTING HISTOGRAM OF DET COSTS
sub2 = subplot(4,1,2)
grid()
hold(true)
hist2 = plt[:hist](obj_st_fp_socfix_all,nbins,weights=ones(length(obj_st_fp_socfix_all))/length(obj_st_fp_socfix_all),color="red",label="Restriction on states");
# PLOTTING HISTOGRAM OF PI COSTS
sub3 = subplot(4,1,3)
grid()
hold(true)
hist3 = plt[:hist](obj_st_fp_all,nbins,weights=ones(length(obj_st_fp_all))/length(obj_st_fp_all),color="blue",label="Stochastic");
# PLOTTING HISTOGRAM OF PI COSTS
sub4 = subplot(4,1,4)
grid()
hold(true)
hist4 = plt[:hist](obj_pi_fp_all,nbins,weights=ones(length(obj_pi_fp_all))/length(obj_pi_fp_all),color="green",label="Perfect Information");
# Setting xlims and ylims based on default lims of the two subplots
xlims = [collect(sub1[:get_xlim]()) collect(sub2[:get_xlim]()) collect(sub3[:get_xlim]()) collect(sub4[:get_xlim]())];
xlims = [minimum(xlims),maximum(xlims)];
ylims = [collect(sub1[:get_ylim]()) collect(sub2[:get_ylim]()) collect(sub3[:get_ylim]()) collect(sub4[:get_ylim]())];
ylims = [minimum(ylims),maximum(ylims)];
sub1[:set_xlim](xlims)
sub1[:set_ylim](ylims)
sub2[:set_xlim](xlims)
sub2[:set_ylim](ylims)
sub3[:set_xlim](xlims)
sub3[:set_ylim](ylims)
sub4[:set_xlim](xlims)
sub4[:set_ylim](ylims)
subplot(sub1)
plot(mean(obj_dt_fp_all)*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Expected Revenue") # Plot vertical line at expected value
ylabel("Probability",size = 14)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
subplot(sub2)
plot(mean(obj_st_fp_socfix_all)*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Expected Revenue") # Plot vertical line at expected value
ylabel("Probability",size = 14)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
subplot(sub3)
plot(mean(obj_st_fp_all)*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Expected Revenue") # Plot vertical line at expected value
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
ylabel("Probability",size = 14)
subplot(sub4)
plot(mean(obj_pi_fp_all)*ones(2),ylims,LineWidth=3,color="black",LineStyle="--",label="Expected Revenue") # Plot vertical line at expected value
ylabel("Probability",size = 14)
xlabel("Revenue (\$)",size = 24)
legend(loc="upper left",fancybox="True", shadow="True", fontsize = 14)
savefig(string(figuredirectory,"/histogram_costs.pdf"))
close("all")
