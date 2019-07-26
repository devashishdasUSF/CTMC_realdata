cpp("queue_sim.cpp")


KS0_intensity = KS_train()
save(KS0_intensity, file = "KS0_intensity.RData")

# load("KS0_intensity.RData")
x = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.0)
cut_KS  = quantile(x[,1], .90)
y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.1)
sum(y[,1]<cut_KS)/length(y[,1])
y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.2)
sum(y[,1]<cut_KS)/length(y[,1])
y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.3)
sum(y[,1]<cut_KS)/length(y[,1])
y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.5)
sum(y[,1]<cut_KS)/length(y[,1])
y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 2.0)
sum(y[,1]<cut_KS)/length(y[,1])


result = simulate()
intensity = (result$N)/(result$Y)
save(intensity, file = "Beta0.RData")

# load("Beta0.RData")
x  = test_stat(Beta0 = intensity)

cut  = quantile(x[,1], .90)
y = test_stat(Beta0 = intensity, slow = 1.0);
sum(y[,1]<cut)/length(y[,1])
y = test_stat(Beta0 = intensity, slow = 1.1);
sum(y[,1]<cut)/length(y[,1])
y = test_stat(Beta0 = intensity, slow = 1.2);
sum(y[,1]<cut)/length(y[,1])
y = test_stat(Beta0 = intensity, slow = 1.3);
sum(y[,1]<cut)/length(y[,1])
y = test_stat(Beta0 = intensity, slow = 1.5);
sum(y[,1]<cut)/length(y[,1])
y = test_stat(Beta0 = intensity, slow = 2.0);
sum(y[,1]<cut)/length(y[,1])

# ==========================================================================
#  RESULTS
# ==========================================================================


# > x  = test_stat(Beta0 = intensity)
# >
# > cut  = quantile(x[,1], .90)
# > y = test_stat(Beta0 = intensity, slow = 1.0);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.90161
# > y = test_stat(Beta0 = intensity, slow = 1.1);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.86198
# > y = test_stat(Beta0 = intensity, slow = 1.2);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.81568
# > y = test_stat(Beta0 = intensity, slow = 1.3);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.76586
# > y = test_stat(Beta0 = intensity, slow = 1.5);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.65191
# > y = test_stat(Beta0 = intensity, slow = 2.0);
# > sum(y[,1]<cut)/length(y[,1])
# [1] 0.39273


# > load("KS0_intensity.RData")
# > x = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.0)
# > cut_KS  = quantile(x[,1], .90)
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.1)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.89262
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.2)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.88392
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.3)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.87169
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.5)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.84644
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 2.0)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.76957


#
# CORRECTED RESUTS
#
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.1)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.89155
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.2)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.88304
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.3)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.87085
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 1.5)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.8442
# > y = KS_test_stat(KS0_intensity = KS0_intensity, slow = 2.0)
# > sum(y[,1]<cut_KS)/length(y[,1])
# [1] 0.767

