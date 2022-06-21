
lam_vals <- c(0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,30,30,30,30,30,30) #dummy array of first timepoint values

lam_vals_used <- vector()
low_sig_lim_list <- vector()

counter = 0

for (p in 1:length(lam_vals)) {
  
  low_sig_lim <- 2*lam_vals[p]
    
  if (counter == 0) {
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
  }
  
  #else if (lam_vals_used[!(lam_vals_used %in% lam_vals[p])] == lam_vals_used) { #this means the lam_vals value is not already in lam_vals_used. So append to list
  else if (setequal(lam_vals_used[!(lam_vals_used %in% lam_vals[p])], lam_vals_used)) {
    print(paste("lam_vals_used[!(lam_vals_used %in% lam_vals[p])] == lam_vals_used",lam_vals_used[!(lam_vals_used %in% lam_vals[p])] == lam_vals_used))
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
  }
  
  lam_vals_used <- append(lam_vals_used, lam_vals[p])
  counter = counter + 1
  
  
  print(paste("lam_vals_used",lam_vals_used))
  print(paste("low_sig_lim_list",low_sig_lim_list))
}
