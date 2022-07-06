
lam_vals <- c(0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,30,30,30,30,30,30) #dummy array of first timepoint values
lam_vals <- c(10,3,3,2,5,8,7,10,0,2)

lam_vals_used <- vector()
low_sig_lim_list <- vector()

#counter = 0

for (p in 1:length(lam_vals)) {
  
  low_sig_lim <- 2*lam_vals[p]

  if (p == 1) {  #the vector is of length zero in the beginning so we have to append in the first instnace
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
  }
  #here we are checking if lam_vals[p] appears in lam_vals_used. If it does then we do not do the else if loop
  else if (setequal(lam_vals_used[!(lam_vals_used %in% lam_vals[p])], lam_vals_used)) { 
    print(paste("lam_vals_used[!(lam_vals_used %in% lam_vals[p])] == lam_vals_used",lam_vals_used[!(lam_vals_used %in% lam_vals[p])] == lam_vals_used))
    low_sig_lim_list <- append(low_sig_lim_list, low_sig_lim)
  }
  
  lam_vals_used <- append(lam_vals_used, lam_vals[p])

  print(paste("lam_vals_used",lam_vals_used))
  print(paste("low_sig_lim_list",low_sig_lim_list))
}



