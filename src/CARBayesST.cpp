#include <Rcpp.h>
using namespace Rcpp;



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// General helper functions   //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector linpredcompute(NumericMatrix X, const int nsites, const int p, 
                          NumericVector beta, NumericVector offset)
{
//Create new objects
// Compute the linear predictor
NumericVector linpred(nsites);
double temp; 


//  Compute the linear predictor via a double for loop
     for(int j = 0; j < nsites; j++)
     {
     temp = 0;
      
          for(int l = 0; l < p; l++) temp = temp + X(j,l) * beta[l];     
          
     linpred[j] = temp + offset[j];  
     }


// Return the result
return linpred;
}




// [[Rcpp::export]]
double quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                    NumericVector phi, NumericVector theta, double rho)
{
// Compute a quadratic form for the random effects
// Create new objects 
double tau2_posteriorscale;
double tau2_quadform = 0, tau2_phisq = 0;
int row, col;
   
   
// Compute the off diagonal elements of the quadratic form
     for(int l = 0; l < n_triplet; l++)
     {
     row = Wtriplet(l,0) - 1;
     col = Wtriplet(l,1) - 1;
     tau2_quadform = tau2_quadform + phi[(Wtriplet(l,0) - 1)] * theta[(Wtriplet(l,1) - 1)] * Wtriplet(l,2); 
     }
 
 
 // Compute the diagonal elements of the quadratic form          
     for(int l = 0; l < nsites; l++)
     {
     tau2_phisq = tau2_phisq + phi[l] * theta[l] * (rho * Wtripletsum[l] + 1 - rho);    
     }
           
     
// Compute the quadratic form
tau2_posteriorscale = 0.5 * (tau2_phisq - rho * tau2_quadform);

 
// Return the simulated value
return tau2_posteriorscale;
}




// [[Rcpp::export]]
List gammaquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, 
                          const int n_triplet, const int nsites, const int ntime, NumericMatrix phi, double rho)
{    
  NumericVector phi_t(nsites), phi_tminus1(nsites);
  double num=0, den=0;
  
  // Compute the sum of quadratic forms for updating gamma
  for(int t = 1; t < ntime; t++)
  {
    phi_t = phi(_,t);
    phi_tminus1 = phi(_,(t-1));    
    num = num + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_t, phi_tminus1, rho);
    den = den + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus1, phi_tminus1, rho);    
  }
  
  
  List out(2);
  out[0] = num;
  out[1] = den;
  
  return out;
}




// [[Rcpp::export]]
List alphaquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, 
                          const int n_triplet, const int nsites, const int ntime, NumericMatrix phi, double rho, double tau2)
{    
  NumericVector phi_t(nsites), phi_tminus1(nsites), phi_tminus2(nsites);
  double B11=0, B22=0, B12=0, M1=0, M2=0;
  
  // Compute the sum of quadratic forms for updating alpha
  for(int t = 2; t < ntime; t++)
  {
    // Create the three elements phi_t, phi_t-1, phi_t-2
    phi_t = phi(_,t);
    phi_tminus1 = phi(_,(t-1));  
    phi_tminus2 = phi(_,(t-2));  
    
    // Add to the quadratic form
    B11 = B11 + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus1, phi_tminus1, rho);
    B22 = B22 + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus2, phi_tminus2, rho);
    B12 = B12 + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_tminus1, phi_tminus2, rho);
    M1 = M1 + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_t, phi_tminus1, rho);
    M2 = M2 + 2 * quadform(Wtriplet, Wtripletsum, n_triplet, nsites, phi_t, phi_tminus2, rho);
  }
  
  
  List out(5);
  out[0] = B11 / tau2;
  out[1] = B22 / tau2;
  out[2] = B12 / tau2;
  out[3] = M1 / tau2;
  out[4] = M2 / tau2;
  
  return out;
}



// [[Rcpp::export]]
double tauquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                          const int nsites, const int ntime, NumericMatrix phi, double rho, double gamma)
{    
  NumericVector temp(nsites);
  double num=0;
  
  // Compute the sum of quadratic forms for updating tau
  temp = phi(_,0);
  num = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
  
  for(int t = 1; t < ntime; t++)
  {
    temp = phi(_,t) - gamma * phi(_,(t-1));  
    num = num + quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
  }
  
  
  return num;
}




// [[Rcpp::export]]
double tauquadformcomputear2(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                             const int nsites, const int ntime, NumericMatrix phi, double rho, double alpha1,
                             double alpha2)
{    
  NumericVector temp(nsites);
  double num=0, num1=0, num2=0;
  
  // Compute the sum of quadratic forms for updating tau
  // First two quadratic forms
  temp = phi(_,0);
  num1 = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
  temp = phi(_,1);
  num2 = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
  num = num1 + num2;
  
  // Quadratic form over AR2 part
  for(int t = 2; t < ntime; t++)
  {
    temp = phi(_,t) - alpha1 * phi(_,(t-1)) - alpha2 * phi(_,(t-2));  
    num = num + quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
  }
  
  // Return the result
  return num;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Regression parameter updates   //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                           NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                           NumericVector prior_varbeta, const int nblock, double beta_tune, 
                           List block_list)
{
    // Compute the acceptance probability for beta
    //Create new objects
    int accept=0;
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    double acceptance;
    NumericVector lp_current(nsites), lp_proposal(nsites), mala_temp1(nsites);
    
    // Create two beta vectors
    NumericVector beta_old(p);
    NumericVector beta_new(p);
    for(int g=0; g<p; g++)
    {
        beta_old[g] = beta[g];
        beta_new[g] = beta[g];
    }
    
    // Update each block in turn
    for(int r=0; r<nblock; r++)
    {
        // Determine the block to update
        IntegerVector idx = block_list[r];
        int len = block_list[(nblock+r)];
        
        // Propose a value
        lp_current = linpredcompute(X, nsites, p, beta_old, offset);
        mala_temp1 = y - exp(lp_current);
        NumericVector mala_temp2(len), mala_old(len);
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
        }
        
        // Compute the acceptance ratio - full conditionals 
        oldlikebit = 0;
        newlikebit=0;
        lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
        for(int j = 0; j < nsites; j++)     
        {
            oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
            newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
        }
        likebit = newlikebit - oldlikebit;
        
        for(int g = 0; g < len; g++)     
        {
            priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
        }
        
        // Compute the acceptance ratio - proposal distributions
        mala_temp1 = y - exp(lp_proposal);
        NumericVector mala_new(len);
        double prop_accept=0;
        for(int g=0; g<len; g++)
        {
            mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
            mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
            prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
        }
        
        // Accept or reject the proposal      
        acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
        if(runif(1)[0] <= acceptance) 
        {
            for(int g=0; g<len; g++)
            {
                beta_old[idx[g]] = beta_new[idx[g]];  
            }
            accept = accept + 1;
        }
        else
        { 
            for(int g=0; g<len; g++)
            {
                beta_new[idx[g]] = beta_old[idx[g]];  
            }   
        }
    }
    
    
    
    // Compute the acceptance probability and return the value
    //acceptance = exp(likebit + priorbit);
    List out(2);
    out[0] = beta_new;
    out[1] = accept;
    return out;    
}




// [[Rcpp::export]]
List poissonbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                           NumericVector offset, NumericVector y, NumericVector prior_meanbeta, 
                           NumericVector prior_varbeta, const int nblock, double beta_tune, 
                           List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }
  
  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
    beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    

    // Compute the acceptance ratio - likelihood part
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    oldlikebit = 0;
    newlikebit=0;
    for(int j = 0; j < nsites; j++)     
    {
      oldlikebit = oldlikebit + y[j] * lp_current[j] - exp(lp_current[j]);
      newlikebit = newlikebit + y[j] * lp_proposal[j] - exp(lp_proposal[j]);
    }
    likebit = newlikebit - oldlikebit;

    // Compute the acceptance ratio - prior part
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }

    // Accept or reject the proposal    
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }
 
  // Return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}




// [[Rcpp::export]]
List binomialbetaupdateMALA(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                            NumericVector offset, NumericVector y,  NumericVector failures,
                            NumericVector trials, NumericVector prior_meanbeta, 
                            NumericVector prior_varbeta, const int nblock,double beta_tune, 
                            List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), mala_temp1(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }
  
  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    mala_temp1 = y -  trials * exp(lp_current) / (1 + exp(lp_current));
    NumericVector mala_temp2(len), mala_old(len);
    for(int g=0; g<len; g++)
    {
      mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
      mala_old[g] = beta_old[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_old[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
      beta_new[idx[g]] = rnorm(1, mala_old[g], beta_tune)[0];
    }
    
    // Compute the acceptance ratio - full conditionals  
    oldlikebit = 0;
    newlikebit=0;
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }
    
    // Compute the acceptance ratio - proposal distributions
    mala_temp1 = y -  trials * exp(lp_proposal) / (1 + exp(lp_proposal));
    NumericVector mala_new(len);
    double prop_accept=0;
    for(int g=0; g<len; g++)
    {
      mala_temp2[g] = sum(X(_,idx[g]) * mala_temp1);
      mala_new[g] = beta_new[idx[g]] + 0.5 * pow(beta_tune,2) * (-(beta_new[idx[g]] - prior_meanbeta[idx[g]]) / prior_varbeta[idx[g]] + mala_temp2[g]); 
      prop_accept = prop_accept +   pow((beta_new[idx[g]] - mala_old[g]), 2) -  pow((beta_old[idx[g]] - mala_new[g]), 2); 
    }
    
    // Accept or reject hte proposal      
    acceptance = exp(0.5 * prop_accept / pow(beta_tune,2) + likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }
  
  
  // Compute the acceptance probability and return the value
  //acceptance = exp(likebit + priorbit);
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}


// [[Rcpp::export]]
List binomialbetaupdateRW(NumericMatrix X, const int nsites, const int p, NumericVector beta, 
                          NumericVector offset, NumericVector y,  NumericVector failures,
                          NumericVector prior_meanbeta, NumericVector prior_varbeta, 
                          const int nblock,double beta_tune, List block_list)
{
  // Compute the acceptance probability for beta
  //Create new objects
  int accept=0;
  double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
  double acceptance;
  NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
  
  // Create two beta vectors
  NumericVector beta_old(p);
  NumericVector beta_new(p);
  for(int g=0; g<p; g++)
  {
    beta_old[g] = beta[g];
    beta_new[g] = beta[g];
  }
  
  // Update each block in turn
  for(int r=0; r<nblock; r++)
  {
    // Determine the block to update
    IntegerVector idx = block_list[r];
    int len = block_list[(nblock+r)];
    
    // Propose a value
    for(int g=0; g<len; g++)
    {
      beta_new[idx[g]] = rnorm(1, beta_old[idx[g]], beta_tune)[0];
    }
    
    
    // Compute the acceptance ratio - full conditionals  
    oldlikebit = 0;
    newlikebit=0;
    lp_current = linpredcompute(X, nsites, p, beta_old, offset);
    lp_proposal = linpredcompute(X, nsites, p, beta_new, offset);     
    for(int j = 0; j < nsites; j++)     
    {
      p_current[j] = exp(lp_current[j]) / (1 + exp(lp_current[j]));
      p_proposal[j] = exp(lp_proposal[j]) / (1 + exp(lp_proposal[j]));
      oldlikebit = oldlikebit +  y[j] * log(p_current[j]) + failures[j] * log((1-p_current[j]));
      newlikebit = newlikebit +  y[j] * log(p_proposal[j]) + failures[j] * log((1-p_proposal[j]));
    }
    likebit = newlikebit - oldlikebit;
    
    priorbit = 0;
    for(int g = 0; g < len; g++)     
    {
      priorbit = priorbit + 0.5 * pow((beta_old[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]] - 0.5 * pow((beta_new[idx[g]]-prior_meanbeta[idx[g]]),2) / prior_varbeta[idx[g]];
    }
    
    
    // Accept or reject hte proposal      
    acceptance = exp(likebit + priorbit);
    if(runif(1)[0] <= acceptance) 
    {
      for(int g=0; g<len; g++)
      {
        beta_old[idx[g]] = beta_new[idx[g]];  
      }
      accept = accept + 1;
    }
    else
    { 
      for(int g=0; g<len; g++)
      {
        beta_new[idx[g]] = beta_old[idx[g]];  
      }   
    }
  }
  
  
  // Compute the acceptance probability and return the value
  List out(2);
  out[0] = beta_new;
  out[1] = accept;
  return out;    
}






////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Spatial CAR and indep updates   /////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissoncarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                      NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                      double tau2, const NumericMatrix y, const double phi_tune, 
                      double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        // Accept or reject it
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        oldlikebit = 0;
        newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
            lpold = mult_offset[i] * phinew[j] + offset(j, i);
            lpnew = mult_offset[i] * propphi + offset(j, i); 
            oldlikebit = oldlikebit + y(j,i) * lpold - exp(lpold);
            newlikebit = newlikebit + y(j,i) * lpnew - exp(lpnew);
        }       
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List poissonindepupdateRW(const int nsites, NumericVector theta,double tau2, 
                        const NumericVector y, const double theta_tune, NumericVector offset)
{
    // Update the spatially independent random effects 
    //Create new objects
    int accept=0;
    double acceptance, proptheta, lpold, lpnew;
    double priorbit, oldlikebit, newlikebit;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // propose a value  
        proptheta = rnorm(1, thetanew[j], sqrt(theta_tune))[0];
        
        // Accept or reject it
        priorbit = (0.5/tau2) * (pow(thetanew[j], 2) - pow(proptheta, 2));
        lpold = thetanew[j] + offset[j];
        lpnew = proptheta + offset[j];
        oldlikebit = lpold * y[j]  - exp(lpold);
        newlikebit = lpnew * y[j]  - exp(lpnew);
        acceptance = exp(priorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            thetanew[j] = proptheta;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List binomialindepupdateRW(const int nsites, NumericVector theta, double tau2, 
                         const NumericVector y, const NumericVector failures, 
                         const double theta_tune, NumericVector offset)
{
    // Update the spatially independent random effects 
    //Create new objects
    int accept=0;
    double acceptance, proptheta, lpold, lpnew, pold, pnew;
    double priorbit, oldlikebit, newlikebit;
    NumericVector thetanew(nsites);
    
    
    //  Update each random effect in turn
    thetanew = theta;
    for(int j = 0; j < nsites; j++)
    {
        // propose a value  
        proptheta = rnorm(1, thetanew[j], sqrt(theta_tune))[0];
        
        // Accept or reject it
        priorbit = (0.5/tau2) * (pow(thetanew[j], 2) - pow(proptheta, 2));
        lpold = thetanew[j] + offset[j];
        lpnew = proptheta + offset[j];
        pold = exp(lpold) / (1 + exp(lpold));
        pnew = exp(lpnew) / (1 + exp(lpnew));
        
        oldlikebit = y[j] * log(pold)  + failures[j] * log((1-pold));
        newlikebit = y[j] * log(pnew)  + failures[j] * log((1-pnew));
        acceptance = exp(priorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            thetanew[j] = proptheta;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = thetanew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
List binomialcarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                       NumericVector Wtripletsum, const int nsites, NumericVector phi, 
                       double tau2, const NumericMatrix y,  const NumericMatrix failures, const double phi_tune, 
                       double rho, NumericMatrix offset, const int ntime, NumericVector mult_offset)
{
    // Update the spatially correlated random effects 
    //Create new objects
    int accept=0,rowstart=0, rowend=0;
    double acceptance, sumphi;
    double oldpriorbit, newpriorbit, oldlikebit, newlikebit;
    double priorvardenom, priormean, priorvar;
    double propphi, lpold, lpnew, pold, pnew;
    NumericVector phinew(nsites);
    
    
    //  Update each random effect in turn
    phinew = phi;
    
    for(int j = 0; j < nsites; j++)
    {
        // Calculate prior variance
        priorvardenom = rho * Wtripletsum[j] + 1 - rho;
        priorvar = tau2 / priorvardenom;
        
        // Calculate the prior mean
        rowstart = Wbegfin(j,0) - 1;
        rowend = Wbegfin(j,1);
        sumphi = 0;
        for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
        priormean = rho * sumphi / priorvardenom; 
        
        // propose a value  
        propphi = rnorm(1, phinew[j], sqrt(priorvar*phi_tune))[0];
        
        // Accept or reject it
        newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
        oldpriorbit = (0.5/priorvar) * pow((phinew[j] - priormean), 2);
        
        oldlikebit = 0;
        newlikebit = 0;
        for(int i=0; i < ntime; i++)
        {
            lpold = mult_offset[i] * phinew[j] + offset(j, i);
            lpnew = mult_offset[i] * propphi + offset(j, i); 
            pold = exp(lpold) / (1 + exp(lpold));
            pnew = exp(lpnew) / (1 + exp(lpnew));        
            
            oldlikebit = oldlikebit + y(j,i) * log(pold) + failures(j,i) * log((1-pold));
            newlikebit = newlikebit + y(j,i) * log(pnew) + failures(j,i) * log((1-pnew));
        }
        
        acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
        if(runif(1)[0] <= acceptance) 
        {
            phinew[j] = propphi;
            accept = accept + 1;
        }
        else
        { 
        }
    }
    
    
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}




// [[Rcpp::export]]
NumericVector gaussiancarupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
     NumericVector Wtripletsum, const int nsites, NumericVector phi, 
     double tau2, double nu2, const NumericVector offset, double rho, double ntime)
{
// Update the spatially correlated random effects 
//Create new objects
int rowstart=0, rowend=0;
double sumphi;
double fcmean, fcvar;
double priorvardenom, priormean, priorvar;
NumericVector phinew(nsites);

   
//  Update each random effect in turn
phinew = phi;

     for(int j = 0; j < nsites; j++)
    {
     // Calculate prior variance
     priorvardenom = rho * Wtripletsum[j] + 1 - rho;
     priorvar = tau2 / priorvardenom;
     
     // Calculate the prior mean
     rowstart = Wbegfin(j,0) - 1;
     rowend = Wbegfin(j,1);
     sumphi = 0;
          for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew[(Wtriplet(l,1) - 1)];
     priormean = rho * sumphi / priorvardenom; 
     
      // compute the full conditional  
      fcvar = 1 / (1 / priorvar + ntime / nu2);
      fcmean = fcvar * (priormean / priorvar +  offset[j]/ nu2);
      phinew[j] = rnorm(1, fcmean, sqrt(fcvar))[0];
     }
        
return phinew;
}









////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// ARCAR modl updates   ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
List poissonar1carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                           NumericVector Wtripletsum, const int nsites, const int ntime,
                           NumericMatrix phi, double tau2, double gamma, double rho, 
                           const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
                           NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * phinew(j,1);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * phinew(row,1));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    oldlikebit = ymat(j,0) * lpold - exp(lpold);
    newlikebit = ymat(j,0) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,0) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(gamma,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp2 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
        priormeantemp2 = priormeantemp2 + temp; 
      }
      priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      oldlikebit = ymat(j,t) * lpold - exp(lpold);
      newlikebit = ymat(j,t) * lpnew - exp(lpnew);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    oldlikebit = ymat(j,(ntime-1)) * lpold - exp(lpold);
    newlikebit = ymat(j,(ntime-1)) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-1)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}




// [[Rcpp::export]]
List poissonar2carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                           NumericVector Wtripletsum, const int nsites, const int ntime,
                           NumericMatrix phi, double tau2, double alpha1, double alpha2, double rho, 
                           const NumericMatrix ymat, const double phi_tune, NumericMatrix offset,
                           NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priormeantemp3, priormeantemp4, priormeantemp5, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(1));
    priormeantemp2 = -alpha2 * denoffset[j] * phinew(j,(2));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha2,2)) * phinew(row,(0)) + alpha1 * alpha2 * phinew(row,(1)) - alpha2 * phinew(row,(2)));   
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;
    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    oldlikebit = ymat(j,0) * lpold - exp(lpold);
    newlikebit = ymat(j,0) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,0) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(0));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(2));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(3));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (alpha1 * alpha2 * phinew(row,(0)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(2)) - alpha2 * phinew(row,(3)));   
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;
    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,1), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,1) - priormean), 2);
    lpold = phinew(j,1) + offset(j, 1);
    lpnew = propphi + offset(j, 1); 
    oldlikebit = ymat(j,1) * lpold - exp(lpold);
    newlikebit = ymat(j,1) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,1) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = -alpha2 * denoffset[j] * phinew(j,(t-2));
      priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t-1));
      priormeantemp3 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t+1));
      priormeantemp4 = -alpha2 * denoffset[j] * phinew(j,(t+2));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * (- alpha2 * phinew(row,(t-2)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t-1)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(t)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t+1)) - alpha2 * phinew(row,(t+2)));
        priormeantemp5 = priormeantemp5 + temp; 
      }
      priormean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / priorvardenom;
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      oldlikebit = ymat(j,t) * lpold - exp(lpold);
      newlikebit = ymat(j,t) * lpnew - exp(lpnew);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  
  
  ////////////////////////////////////////////////
  // Update the random effects at time N-1 in turn
  ////////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * phinew(j,(ntime-1));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(ntime-3));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(ntime-4));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha1,2)) *  phinew(row,(ntime-2))  - alpha1 * phinew(row,(ntime-1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(ntime-3)) - alpha2 * phinew(row,(ntime-4)));
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-2)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-2)) - priormean), 2);
    lpold = phinew(j,(ntime-2)) + offset(j, (ntime-2));
    lpnew = propphi + offset(j, (ntime-2)); 
    oldlikebit = ymat(j,(ntime-2)) * lpold - exp(lpold);
    newlikebit = ymat(j,(ntime-2)) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-2)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * (phinew(j,(ntime-2)));
    priormeantemp2 = -alpha2 * denoffset[j] * (phinew(j,(ntime-3)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - alpha1 * (phinew(row,(ntime-2))) - alpha2 * (phinew(row,(ntime-3))));
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    oldlikebit = ymat(j,(ntime-1)) * lpold - exp(lpold);
    newlikebit = ymat(j,(ntime-1)) * lpnew - exp(lpnew);
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-1)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}




// [[Rcpp::export]]
List binomialar1carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                            NumericVector Wtripletsum, const int nsites, const int ntime,
                            NumericMatrix phi, double tau2, double gamma, double rho, 
                            const NumericMatrix ymat, const NumericMatrix failuresmat,
                            const double phi_tune, NumericMatrix offset,NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold));
    newlikebit = ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,0) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(gamma,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp2 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
        priormeantemp2 = priormeantemp2 + temp; 
      }
      priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      pold = exp(lpold) / (1 + exp(lpold));
      pnew = exp(lpnew) / (1 + exp(lpnew));        
      oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
      newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold));
    newlikebit = ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-1)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}




// [[Rcpp::export]]
List binomialar2carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                            NumericVector Wtripletsum, const int nsites, const int ntime,
                            NumericMatrix phi, double tau2, double alpha1, double alpha2, double rho, 
                            const NumericMatrix ymat, const NumericMatrix failuresmat, const double phi_tune, 
                            NumericMatrix offset, NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priormeantemp3, priormeantemp4, priormeantemp5, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(1));
    priormeantemp2 = -alpha2 * denoffset[j] * phinew(j,(2));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha2,2)) * phinew(row,(0)) + alpha1 * alpha2 * phinew(row,(1)) - alpha2 * phinew(row,(2)));   
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;
    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,0), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,0) - priormean), 2);
    lpold = phinew(j,0) + offset(j, 0);
    lpnew = propphi + offset(j, 0); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,0) * log(pold) + failuresmat(j,0) * log((1-pold));
    newlikebit = ymat(j,0) * log(pnew) + failuresmat(j,0) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,0) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(0));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(2));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(3));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (alpha1 * alpha2 * phinew(row,(0)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(2)) - alpha2 * phinew(row,(3)));   
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;
    
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,1), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,1) - priormean), 2);
    lpold = phinew(j,1) + offset(j, 1);
    lpnew = propphi + offset(j, 1); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,1) * log(pold) + failuresmat(j,1) * log((1-pold));
    newlikebit = ymat(j,1) * log(pnew) + failuresmat(j,1) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,1) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = -alpha2 * denoffset[j] * phinew(j,(t-2));
      priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t-1));
      priormeantemp3 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t+1));
      priormeantemp4 = -alpha2 * denoffset[j] * phinew(j,(t+2));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * (- alpha2 * phinew(row,(t-2)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t-1)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(t)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t+1)) - alpha2 * phinew(row,(t+2)));
        priormeantemp5 = priormeantemp5 + temp; 
      }
      priormean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / priorvardenom;
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      pold = exp(lpold) / (1 + exp(lpold));
      pnew = exp(lpnew) / (1 + exp(lpnew));        
      oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
      newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  
  
  ////////////////////////////////////////////////
  // Update the random effects at time N-1 in turn
  ////////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * phinew(j,(ntime-1));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(ntime-3));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(ntime-4));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha1,2)) *  phinew(row,(ntime-2))  - alpha1 * phinew(row,(ntime-1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(ntime-3)) - alpha2 * phinew(row,(ntime-4)));
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-2)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-2)) - priormean), 2);
    lpold = phinew(j,(ntime-2)) + offset(j, (ntime-2));
    lpnew = propphi + offset(j, (ntime-2)); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,(ntime-2)) * log(pold) + failuresmat(j,(ntime-2)) * log((1-pold));
    newlikebit = ymat(j,(ntime-2)) * log(pnew) + failuresmat(j,(ntime-2)) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-2)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * (phinew(j,(ntime-2)));
    priormeantemp2 = -alpha2 * denoffset[j] * (phinew(j,(ntime-3)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - alpha1 * (phinew(row,(ntime-2))) - alpha2 * (phinew(row,(ntime-3))));
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;   
    
    // Propose a value and calculate the acceptance probability
    propphi = rnorm(1, phinew(j,(ntime-1)), sqrt(priorvar*phi_tune))[0];
    newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
    oldpriorbit = (0.5/priorvar) * pow((phinew(j,(ntime-1)) - priormean), 2);
    lpold = phinew(j,(ntime-1)) + offset(j, (ntime-1));
    lpnew = propphi + offset(j, (ntime-1)); 
    pold = exp(lpold) / (1 + exp(lpold));
    pnew = exp(lpnew) / (1 + exp(lpnew));        
    oldlikebit = ymat(j,(ntime-1)) * log(pold) + failuresmat(j,(ntime-1)) * log((1-pold));
    newlikebit = ymat(j,(ntime-1)) * log(pnew) + failuresmat(j,(ntime-1)) * log((1-pnew));
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(j,(ntime-1)) = propphi;
      accept = accept + 1;
    }
    else
    { 
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}



// [[Rcpp::export]]
NumericMatrix gaussianar1carupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                   NumericVector Wtripletsum, const int nsites, const int ntime,
                                   NumericMatrix phi, double tau2, double nu2, double gamma, double rho, 
                                   NumericMatrix offset, NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priorvar,priorvardenom;
  double propphi, fcvar, fcmean; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(gamma,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,1));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,0) - gamma * (phinew(row,1)));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,0)/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,0) = propphi;
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(gamma,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = gamma * denoffset[j] * (phinew(j,(t-1)) + phinew(j, (t+1)));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp2 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * ((1 + pow(gamma,2)) * phinew(row,t) - gamma * (phinew(row,(t-1)) + phinew(row,(t+1))));   
        priormeantemp2 = priormeantemp2 + temp; 
      }
      priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
      
      // Compute the full conditional and update phi
      fcvar = 1 / (1 / priorvar + 1 / nu2);
      fcmean = fcvar * (priormean / priorvar +  offset(j,t)/ nu2);
      propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
      phinew(j,t) = propphi;
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = gamma * denoffset[j] * (phinew(j,(ntime-2)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - gamma * (phinew(row,(ntime-2))));   
      priormeantemp2 = priormeantemp2 + temp; 
    }
    priormean = (priormeantemp1 + rho * priormeantemp2) / priorvardenom; 
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,(ntime-1))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,(ntime-1)) = propphi;
  }
  
  return phinew;
}




// [[Rcpp::export]]
NumericMatrix  gaussianar2carupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                                    NumericVector Wtripletsum, const int nsites, const int ntime,
                                    NumericMatrix phi, double tau2, double nu2, double alpha1, double alpha2, double rho, 
                                    NumericMatrix offset, NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp1, priormeantemp2, priormeantemp3, priormeantemp4, priormeantemp5, priorvar, priorvardenom;
  double propphi, fcvar, fcmean; 
  NumericMatrix phinew(nsites,ntime);
  phinew = clone(phi);
  int row, rowstart, rowend;
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(1));
    priormeantemp2 = -alpha2 * denoffset[j] * phinew(j,(2));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha2,2)) * phinew(row,(0)) + alpha1 * alpha2 * phinew(row,(1)) - alpha2 * phinew(row,(2)));   
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,0) / nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,0) = propphi;
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] * phinew(j,(0));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(2));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(3));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (alpha1 * alpha2 * phinew(row,(0)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(2)) - alpha2 * phinew(row,(3)));   
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,1)/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,1) = propphi;
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j] * (1 + pow(alpha1,2) + pow(alpha2,2));
      priorvar = tau2 / priorvardenom;
      priormeantemp1 = -alpha2 * denoffset[j] * phinew(j,(t-2));
      priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t-1));
      priormeantemp3 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(t+1));
      priormeantemp4 = -alpha2 * denoffset[j] * phinew(j,(t+2));
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * (- alpha2 * phinew(row,(t-2)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t-1)) + (1 + pow(alpha1,2) + pow(alpha2,2)) * phinew(row,(t)) + (alpha1 * alpha2 - alpha1) * phinew(row,(t+1)) - alpha2 * phinew(row,(t+2)));
        priormeantemp5 = priormeantemp5 + temp; 
      }
      priormean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / priorvardenom;
      
      // Compute the full conditional and update phi
      fcvar = 1 / (1 / priorvar + 1 / nu2);
      fcmean = fcvar * (priormean / priorvar +  offset(j,t)/ nu2);
      propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
      phinew(j,t) = propphi;
    }
  }
  
  
  
  ////////////////////////////////////////////////
  // Update the random effects at time N-1 in turn
  ////////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j] * (1 + pow(alpha1,2));
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * phinew(j,(ntime-1));
    priormeantemp2 = (alpha1 * alpha2 - alpha1) *  denoffset[j] * phinew(j,(ntime-3));
    priormeantemp3 = -alpha2 * denoffset[j] * phinew(j,(ntime-4));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp4 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * ((1 + pow(alpha1,2)) *  phinew(row,(ntime-2))  - alpha1 * phinew(row,(ntime-1)) + (alpha1 * alpha2 - alpha1) * phinew(row,(ntime-3)) - alpha2 * phinew(row,(ntime-4)));
      priormeantemp4 = priormeantemp4 + temp; 
    }
    priormean = (rho * priormeantemp4 - priormeantemp1 - priormeantemp2 - priormeantemp3) / priorvardenom;   
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,(ntime-2))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,(ntime-2)) = propphi;
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  for(int j = 0; j < nsites; j++)
  {
    // calculate prior mean and variance
    priorvardenom = denoffset[j];
    priorvar = tau2 / priorvardenom;
    priormeantemp1 = -alpha1 * denoffset[j] * (phinew(j,(ntime-2)));
    priormeantemp2 = -alpha2 * denoffset[j] * (phinew(j,(ntime-3)));
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = 0;
    for(int l = rowstart; l < rowend; l++)
    {
      row = Wtriplet(l,1) - 1;
      temp = Wtriplet(l, 2) * (phinew(row,(ntime-1)) - alpha1 * (phinew(row,(ntime-2))) - alpha2 * (phinew(row,(ntime-3))));
      priormeantemp3 = priormeantemp3 + temp; 
    }
    priormean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / priorvardenom;   
    
    // Compute the full conditional and update phi
    fcvar = 1 / (1 / priorvar + 1 / nu2);
    fcmean = fcvar * (priormean / priorvar +  offset(j,(ntime-1))/ nu2);
    propphi = rnorm(1, fcmean, sqrt(fcvar))[0];
    phinew(j,(ntime-1)) = propphi;
  }
  
  return phinew;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// MVARCAR functions   /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List binomialmvar1carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                              NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                              NumericMatrix phi, double alpha, double rho, NumericMatrix Sigmainv,
                              const NumericMatrix ymat, const NumericMatrix failuresmat, NumericMatrix innovations, NumericMatrix offset,
                              NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), pold(nvar), pnew(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / ((1 + pow(alpha,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-1) * nsite + j;
      priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp2 = - alpha * denoffset[j] *  phinew(row1,_);
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp3 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * t;
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * (t+1);
        temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _) - alpha * phinew(row_phi3, _);
        priormeantemp3 = priormeantemp3 + temp; 
      }
      fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = 0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
      pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
      oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
      newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
      
      
      // Accept or reject the value
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {  
      }
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime-1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t-1);
      temp = phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / (denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}





// [[Rcpp::export]]
List binomialmvar2carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                              NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                              NumericMatrix phi, double alpha1, double alpha2, double rho, NumericMatrix Sigmainv,
                              const NumericMatrix ymat,  const NumericMatrix failuresmat, NumericMatrix innovations, NumericMatrix offset,
                              NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, row_phi4, row_phi5, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), priormeantemp4(nvar), priormeantemp5(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), pold(nvar), pnew(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp2 = - alpha2 * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      row_phi3 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha2,2)) * phinew(row_phi1,_) + alpha1 * alpha2 * phinew(row_phi2, _) - alpha2 * phinew(row_phi3, _);
      priormeantemp3 = priormeantemp3 + temp; 
    }
    fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  t = 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      row_phi5 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi5, _) + alpha1 * alpha2 * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-2) * nsite + j;
      priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
      row1 = (t-1) * nsite + j;
      priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+2) * nsite + j;
      priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
      
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * (t-2);
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * t;
        row_phi4 = row_W + nsite * (t+1);
        row_phi5 = row_W + nsite * (t+2);
        temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) - alpha2 * phinew(row_phi5, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
        priormeantemp5 = priormeantemp5 + temp; 
      }
      fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = 0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
      pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
      oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
      newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
      
      
      // Accept or reject the value
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {  
      }
    }
  }
  
  
  //////////////////////////////////////////////////
  // Update the random effects at time N - 1 in turn
  //////////////////////////////////////////////////
  t = ntime - 2;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha1,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) - alpha1 * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3) / ((1 + pow(alpha1,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime - 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      temp = phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _)  - alpha1 * phinew(row_phi2, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2) / denoffset[j]; 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    pold = exp(offset(row1,_) + phinew(row1,_)) / (1 + exp(offset(row1,_) + phinew(row1,_)));
    pnew = exp(offset(row1,_) + propphi) / (1 + exp(offset(row1,_) + propphi));
    oldlikebit = sum(ymat(row1,_) * log(pold) + failuresmat(row1,_) * log(1-pold));
    newlikebit = sum(ymat(row1,_) * log(pnew) + failuresmat(row1,_) * log(1-pnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}





// [[Rcpp::export]]
List poissonmvar1carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                             NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                             NumericMatrix phi, double alpha, double rho, NumericMatrix Sigmainv,
                             const NumericMatrix ymat, NumericMatrix innovations, NumericMatrix offset,
                             NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), lpold(nvar), lpnew(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / ((1 + pow(alpha,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-1) * nsite + j;
      priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp2 = - alpha * denoffset[j] *  phinew(row1,_);
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp3 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * t;
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * (t+1);
        temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _) - alpha * phinew(row_phi3, _);
        priormeantemp3 = priormeantemp3 + temp; 
      }
      fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = 0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      lpold = offset(row1,_) + phinew(row1,_);
      lpnew = offset(row1,_) + propphi;
      oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
      newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
      
      
      // Accept or reject the value
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {  
      }
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime-1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t-1);
      temp = phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / (denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}




// [[Rcpp::export]]
List poissonmvar2carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                             NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                             NumericMatrix phi, double alpha1, double alpha2, double rho, NumericMatrix Sigmainv,
                             const NumericMatrix ymat, NumericMatrix innovations, NumericMatrix offset,
                             NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, row_phi4, row_phi5, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), priormeantemp4(nvar), priormeantemp5(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), lpold(nvar), lpnew(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp2 = - alpha2 * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      row_phi3 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha2,2)) * phinew(row_phi1,_) + alpha1 * alpha2 * phinew(row_phi2, _) - alpha2 * phinew(row_phi3, _);
      priormeantemp3 = priormeantemp3 + temp; 
    }
    fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  t = 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      row_phi5 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi5, _) + alpha1 * alpha2 * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-2) * nsite + j;
      priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
      row1 = (t-1) * nsite + j;
      priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+2) * nsite + j;
      priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
      
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * (t-2);
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * t;
        row_phi4 = row_W + nsite * (t+1);
        row_phi5 = row_W + nsite * (t+2);
        temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) - alpha2 * phinew(row_phi5, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
        priormeantemp5 = priormeantemp5 + temp; 
      }
      fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = 0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      lpold = offset(row1,_) + phinew(row1,_);
      lpnew = offset(row1,_) + propphi;
      oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
      newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
      
      
      // Accept or reject the value
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {  
      }
    }
  }
  
  
  //////////////////////////////////////////////////
  // Update the random effects at time N - 1 in turn
  //////////////////////////////////////////////////
  t = ntime - 2;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha1,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) - alpha1 * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3) / ((1 + pow(alpha1,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime - 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      temp = phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _)  - alpha1 * phinew(row_phi2, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2) / denoffset[j]; 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = 0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = 0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    lpold = offset(row1,_) + phinew(row1,_);
    lpnew = offset(row1,_) + propphi;
    oldlikebit = sum(ymat(row1,_) * lpold - exp(lpold));
    newlikebit = sum(ymat(row1,_) * lpnew - exp(lpnew));
    
    
    // Accept or reject the value
    acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}



// [[Rcpp::export]]
List gaussianmvar1carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                              NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                              NumericMatrix phi, double alpha, double rho, NumericMatrix Sigmainv,
                              const NumericVector nu2, NumericMatrix innovations, NumericMatrix offset,
                              NumericVector denoffset)
{    
  //////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), diffold(nvar), diffnew(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / ((1 + pow(alpha,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 2 to N-1 in turn
  //////////////////////////////////////////////////////
  for(int t = 1; t < (ntime-1); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-1) * nsite + j;
      priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp2 = - alpha * denoffset[j] *  phinew(row1,_);
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp3 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * t;
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * (t+1);
        temp = (1 + pow(alpha,2)) * phinew(row_phi1,_) - alpha * phinew(row_phi2, _) - alpha * phinew(row_phi3, _);
        priormeantemp3 = priormeantemp3 + temp; 
      }
      fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = -0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
      diffnew = pow((offset(row1,_) - propphi),2) / nu2;
      oldlikebit = -0.5 * sum(diffold);
      newlikebit = -0.5 * sum(diffnew);
      
      // Accept or reject the value
      acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {   
      }
    }
  }
  
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime-1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp1 = - alpha * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp2 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t-1);
      temp = phinew(row_phi1,_) - alpha * phinew(row_phi2, _);
      priormeantemp2 = priormeantemp2 + temp; 
    }
    fcmean = (rho * priormeantemp2 - priormeantemp1) / (denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}




// [[Rcpp::export]]
List gaussianmvar2carupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, 
                              NumericVector Wtripletsum, const int nsite, const int ntime, const int nvar,
                              NumericMatrix phi, double alpha1, double alpha2, double rho, NumericMatrix Sigmainv,
                              const NumericVector nu2, NumericMatrix innovations, NumericMatrix offset,
                              NumericVector denoffset)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  int NK = ntime * nsite, accept=0;
  int row_W, row1, row_phi1, row_phi2, row_phi3, row_phi4, row_phi5, rowstart, rowend;
  NumericMatrix fcprec(nvar, nvar);
  NumericVector priormeantemp1(nvar), priormeantemp2(nvar), priormeantemp3(nvar), priormeantemp4(nvar), priormeantemp5(nvar), temp(nvar), fcmean(nvar), propphi(nvar);
  NumericVector quadcurrent(nvar), quadprop(nvar), diffcurrent(nvar), diffprop(nvar), diffnew(nvar), diffold(nvar);  
  NumericMatrix phinew(NK, nvar);
  double oldpriorbit, newpriorbit, oldlikebit, newlikebit, acceptance;
  phinew = clone(phi);
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 1 in turn
  //////////////////////////////////////////////
  int t = 0;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t+1) * nsite + j;
    priormeantemp1 = alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp2 = - alpha2 * denoffset[j] *  phinew(row1,_);
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp3 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * t;
      row_phi2 = row_W + nsite * (t+1);
      row_phi3 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha2,2)) * phinew(row_phi1,_) + alpha1 * alpha2 * phinew(row_phi2, _) - alpha2 * phinew(row_phi3, _);
      priormeantemp3 = priormeantemp3 + temp; 
    }
    fcmean = (rho * priormeantemp3 - priormeantemp1 - priormeantemp2) / ((1 + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time 2 in turn
  //////////////////////////////////////////////
  t = 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  alpha1 * alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+2) * nsite + j;
    priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      row_phi5 = row_W + nsite * (t+2);
      temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi5, _) + alpha1 * alpha2 * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at times 3 to N-2 in turn
  //////////////////////////////////////////////////////
  for(int t = 2; t < (ntime-2); t++)
  {
    for(int j = 0; j < nsite; j++)
    {
      // Calculate the prior precision
      for(int r=0; r<nvar; r++)
      {
        fcprec(_,r) = (1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j] * Sigmainv(_,r);  
      }
      
      // calculate prior mean
      row1 = (t-2) * nsite + j;
      priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
      row1 = (t-1) * nsite + j;
      priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+1) * nsite + j;
      priormeantemp3 = (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
      row1 = (t+2) * nsite + j;
      priormeantemp4 = - alpha2 * denoffset[j] *  phinew(row1,_);
      
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp5 = rep(0,nvar);
      for(int l = rowstart; l < rowend; l++)
      {
        row_W = Wtriplet(l,1) - 1;
        row_phi1 = row_W + nsite * (t-2);
        row_phi2 = row_W + nsite * (t-1);
        row_phi3 = row_W + nsite * t;
        row_phi4 = row_W + nsite * (t+1);
        row_phi5 = row_W + nsite * (t+2);
        temp = (1 + pow(alpha1,2) + pow(alpha2,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) - alpha2 * phinew(row_phi5, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi4, _);       
        priormeantemp5 = priormeantemp5 + temp; 
      }
      fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3 - priormeantemp4) / ((1 + pow(alpha1,2) + pow(alpha2,2)) * denoffset[j]); 
      
      // Generate the proposal distribution mean and propose a value
      row1 = t * nsite + j;
      for(int r=0; r<nvar; r++)
      {
        propphi[r] = phinew(row1,r) + innovations(row1, r);
      }
      
      // Compute the prior ratio
      diffcurrent = phinew(row1,_) - fcmean;
      diffprop = propphi - fcmean;
      for(int r=0; r<nvar; r++)
      {
        quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
        quadprop[r] = sum(diffprop * fcprec(_,r));  
      }
      oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
      newpriorbit = -0.5 * sum(quadprop * diffprop);      
      
      // Likelihood ratio
      diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
      diffnew = pow((offset(row1,_) - propphi),2) / nu2;
      oldlikebit = -0.5 * sum(diffold);
      newlikebit = -0.5 * sum(diffnew);
      
      // Accept or reject the value
      acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(row1,_) = propphi;
        accept = accept + 1;
      }
      else
      {  
      }
    }
  }
  
  
  //////////////////////////////////////////////////
  // Update the random effects at time N - 1 in turn
  //////////////////////////////////////////////////
  t = ntime - 2;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = (1 + pow(alpha1,2)) * denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  (alpha1 * alpha2 - alpha1) * denoffset[j] *  phinew(row1,_);
    row1 = (t+1) * nsite + j;
    priormeantemp3 = - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      row_phi4 = row_W + nsite * (t+1);
      temp = (1 + pow(alpha1,2))  * phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _) + (alpha1 * alpha2 - alpha1) * phinew(row_phi2, _) - alpha1 * phinew(row_phi4, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2 - priormeantemp3) / ((1 + pow(alpha1,2)) * denoffset[j]); 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  //////////////////////////////////////////////
  // Update the random effects at time N in turn
  //////////////////////////////////////////////
  t = ntime - 1;
  for(int j = 0; j < nsite; j++)
  {
    // Calculate the prior precision
    for(int r=0; r<nvar; r++)
    {
      fcprec(_,r) = denoffset[j] * Sigmainv(_,r);  
    }
    
    // calculate prior mean
    row1 = (t-2) * nsite + j;
    priormeantemp1 = - alpha2 * denoffset[j] *  phinew(row1,_);
    row1 = (t-1) * nsite + j;
    priormeantemp2 =  - alpha1 * denoffset[j] *  phinew(row1,_);
    
    rowstart = Wbegfin(j,0) - 1;
    rowend = Wbegfin(j,1);
    priormeantemp5 = rep(0,nvar);
    for(int l = rowstart; l < rowend; l++)
    {
      row_W = Wtriplet(l,1) - 1;
      row_phi1 = row_W + nsite * (t-2);
      row_phi2 = row_W + nsite * (t-1);
      row_phi3 = row_W + nsite * t;
      temp = phinew(row_phi3,_) - alpha2 * phinew(row_phi1, _)  - alpha1 * phinew(row_phi2, _);       
      priormeantemp5 = priormeantemp5 + temp; 
    }
    fcmean = (rho * priormeantemp5 - priormeantemp1 - priormeantemp2) / denoffset[j]; 
    
    // Generate the proposal distribution mean and propose a value
    row1 = t * nsite + j;
    for(int r=0; r<nvar; r++)
    {
      propphi[r] = phinew(row1,r) + innovations(row1, r);
    }
    
    // Compute the prior ratio
    diffcurrent = phinew(row1,_) - fcmean;
    diffprop = propphi - fcmean;
    for(int r=0; r<nvar; r++)
    {
      quadcurrent[r] = sum(diffcurrent * fcprec(_,r));  
      quadprop[r] = sum(diffprop * fcprec(_,r));  
    }
    oldpriorbit = -0.5 * sum(quadcurrent * diffcurrent);
    newpriorbit = -0.5 * sum(quadprop * diffprop);      
    
    // Likelihood ratio
    diffold = pow((offset(row1,_) - phinew(row1,_)),2) / nu2;
    diffnew = pow((offset(row1,_) - propphi),2) / nu2;
    oldlikebit = -0.5 * sum(diffold);
    newlikebit = -0.5 * sum(diffnew);
    
    // Accept or reject the value
    acceptance = exp(newpriorbit - oldpriorbit + newlikebit - oldlikebit);
    if(runif(1)[0] <= acceptance) 
    {
      phinew(row1,_) = propphi;
      accept = accept + 1;
    }
    else
    {  
    }
  }
  
  
  /////////////////////
  // Return the results
  /////////////////////
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}






// [[Rcpp::export]]
double MVSTquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, NumericVector den_offset,
                           const int nsite, const int nvar, NumericMatrix phit, NumericMatrix phij, double rho, NumericMatrix Sigmainv)
{    
  // Computes phi_t %*% Q kronecker Sigma.inv %*% phi_j
  double quadform_diag= 0, quadform_offdiag = 0, QF=0;
  int row, col;
  NumericVector temp(nvar);
  
  // Compute the off diagonal elements of the quadratic form
  for(int l = 0; l < n_triplet; l++)
  {
    row = Wtriplet(l,0) - 1;
    col = Wtriplet(l,1) - 1;
    for(int r=0; r<nvar; r++)
    {
      temp[r] = sum(phit(row,_) * Sigmainv(_,r));  
    }
    quadform_offdiag = quadform_offdiag + sum(temp * phij(col,_));
  }
  
  // Compute the diagonal elements of the quadratic form          
  for(int l = 0; l < nsite; l++)
  {
    for(int r=0; r<nvar; r++)
    {
      temp[r] = den_offset[l] * sum(phit(l,_) * Sigmainv(_,r));  
    }
    quadform_diag = quadform_diag + sum(temp * phij(l,_));
  }
  
  
  // Compute the final quadratic form
  QF = -rho * quadform_offdiag + quadform_diag;
  
  
  // Return the result
  return QF;
}






// [[Rcpp::export]]
List MVSTrhoTAR1compute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, NumericVector den_offset,
                        const int nsite, const int ntime, const int nvar, NumericMatrix phi, double rho, NumericMatrix Sigmainv)
{    
  // Computes the numerator and denominator for the temporal AR(1) parameter full conditional
  double num=0, denom=0;
  int range_min, range_max;
  NumericMatrix phi_A(nsite,nvar), phi_B(nsite,nvar);
  
  // Iterate over t=2,...,N-1 as C++ starts loops at zero.
  for(int t=1; t < ntime; t++)
  {
    // Extract phi_t-1
    range_min = (t-1) * nsite;
    range_max = t * nsite - 1;
    phi_A = phi(Range(range_min, range_max),_);
    
    // Extract phi_t
    range_min = t * nsite;
    range_max = (t+1) * nsite - 1;
    phi_B = phi(Range(range_min, range_max),_);
    
    // Compute the quadratic forms
    num = num + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_B, rho, Sigmainv);
    denom = denom + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_A, rho, Sigmainv);
  }
  
  
  // Save the results
  List out(2);
  out[0] = num;
  out[1] = denom;
  return out;
} 



// [[Rcpp::export]]
List MVSTrhoTAR2compute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, NumericVector den_offset,
                        const int nsite, const int ntime, const int nvar, NumericMatrix phi, double rho, NumericMatrix Sigmainv)
{    
  // Computes the 5 different quadratic forms for the full conditional of alpha
  double V11=0, V12=0, V22=0, C=0, D=0;
  int range_min, range_max;
  NumericMatrix phi_A(nsite,nvar), phi_B(nsite,nvar), phi_C(nsite,nvar);
  
  // Iterate over t=2,...,N-1 as C++ starts loops at zero.
  for(int t=2; t < ntime; t++)
  {
    // Extract phi_t
    range_min = t * nsite;
    range_max = (t+1) * nsite - 1;
    phi_A = phi(Range(range_min, range_max),_);
    
    // Extract phi_t-1
    range_min = (t-1) * nsite;
    range_max = t * nsite - 1;
    phi_B = phi(Range(range_min, range_max),_);
    
    // Extract phi_t-2
    range_min = (t-2) * nsite;
    range_max = (t-1) * nsite - 1;
    phi_C = phi(Range(range_min, range_max),_);    
    
    // Compute the quadratic forms
    V11 = V11 + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_B, phi_B, rho, Sigmainv);
    V12 = V12 + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_B, phi_C, rho, Sigmainv);
    V22 = V22 + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_C, phi_C, rho, Sigmainv);  
    C = C + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_B, rho, Sigmainv);
    D = D + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_C, rho, Sigmainv);
  }
  
  
  // Save the results
  List out(5);
  out[0] = V11;
  out[1] = V12;
  out[2] = V22;
  out[3] = C;
  out[4] = D;
  return out;
} 




// [[Rcpp::export]]
double MVSTrhoSAR1compute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, NumericVector den_offset,
                          const int nsite, const int ntime, const int nvar, NumericMatrix phi, double rho, double alpha, NumericMatrix Sigmainv)
{    
  // Computes the quadratic form for the rho update
  double QF=0;
  int range_min, range_max;
  NumericMatrix phi_A(nsite,nvar), phi_B(nsite,nvar), phi_diff(nsite, nvar);
  
  // Compute phi+_1^T %*% Q Kronecker Sigma.inv %*% phi_1
  int t=1;
  range_min = (t-1) * nsite;
  range_max = t * nsite - 1;
  phi_A = phi(Range(range_min, range_max),_);
  QF = MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_A, rho, Sigmainv);
  
  // Iterate over t=2,...,N-1 as C++ starts loops at zero.
  for(int t=1; t < ntime; t++)
  {
    // Extract phi_t-1
    range_min = (t-1) * nsite;
    range_max = t * nsite - 1;
    phi_A = phi(Range(range_min, range_max),_);
    
    // Extract phi_t
    range_min = t * nsite;
    range_max = (t+1) * nsite - 1;
    phi_B = phi(Range(range_min, range_max),_);
    
    // Compute the AR(1) difference
    for(int r=0; r<nvar; r++)
    {
      phi_diff(_,r) = phi_B(_,r) - alpha * phi_A(_,r); 
    }
    
    // Compute the QF
    QF = QF + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_diff, phi_diff, rho, Sigmainv);
  }
  
  
  // Save the results
  return QF;
} 




// [[Rcpp::export]]
double MVSTrhoSAR2compute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, NumericVector den_offset,
                          const int nsite, const int ntime, const int nvar, NumericMatrix phi, double rho, double alpha1, double alpha2, NumericMatrix Sigmainv)
{    
  // Computes the quadratic form for the rho update
  double QF=0;
  int range_min, range_max;
  NumericMatrix phi_A(nsite,nvar), phi_B(nsite,nvar), phi_C(nsite, nvar), phi_diff(nsite, nvar);
  
  // Compute phi+_1^T %*% Q Kronecker Sigma.inv %*% phi_1
  int t=1;
  range_min = (t-1) * nsite;
  range_max = t * nsite - 1;
  phi_A = phi(Range(range_min, range_max),_);
  QF = MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_A, rho, Sigmainv);
  
  // Compute phi+_2^T %*% Q Kronecker Sigma.inv %*% phi_2
  t=2;
  range_min = (t-1) * nsite;
  range_max = t * nsite - 1;
  phi_A = phi(Range(range_min, range_max),_);
  QF = QF + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_A, phi_A, rho, Sigmainv);
  
  // Iterate over t=2,...,N-1 as C++ starts loops at zero.
  for(int t=2; t < ntime; t++)
  {
    // Extract phi_t
    range_min = t * nsite;
    range_max = (t+1) * nsite - 1;
    phi_A = phi(Range(range_min, range_max),_);
    
    // Extract phi_t-1
    range_min = (t-1) * nsite;
    range_max = t * nsite - 1;
    phi_B = phi(Range(range_min, range_max),_);
    
    
    // Extract phi_t-2
    range_min = (t-2) * nsite;
    range_max = (t-1) * nsite - 1;
    phi_C = phi(Range(range_min, range_max),_);
    
    
    // Compute the AR(1) difference
    for(int r=0; r<nvar; r++)
    {
      phi_diff(_,r) = phi_A(_,r) - alpha1 * phi_B(_,r) - alpha2 * phi_C(_,r); 
    }
    
    // Compute the QF
    QF = QF + MVSTquadformcompute(Wtriplet, Wtripletsum, n_triplet, den_offset, nsite, nvar, phi_diff, phi_diff, rho, Sigmainv);
  }
  
  
  // Save the results
  return QF;
} 


////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Adaptive model functions   //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
double qform(NumericMatrix Qtrip, NumericVector phi){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi[Qtrip(i, 1) - 1];
  }
  return(Qform);
}




// [[Rcpp::export]]
double qform_asym(NumericMatrix Qtrip, NumericVector phi1, NumericVector phi2){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  for(int i = 0; i < nzero; i++){
    Qform += phi1[Qtrip(i, 0) - 1] * Qtrip(i, 2) * phi2[Qtrip(i, 1) - 1];
  }
  return(Qform);
}




// [[Rcpp::export]]
double qformSPACETIME(NumericMatrix Qtrip, NumericVector phi, const int ntime, const int nsite){ 
  int nzero = Qtrip.nrow();
  double Qform = 0;
  int spaceBlock = 0;
  for(int j = 0; j < ntime; j++){
    spaceBlock = j*nsite - 1;
    for(int i = 0; i < nzero; i++){
      Qform += phi[Qtrip(i, 0) + spaceBlock] * Qtrip(i, 2) * phi[Qtrip(i, 1) + spaceBlock];
    }
  }
  return(Qform);
}




// [[Rcpp::export]]
List SPTICARphiGaussian(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phi, NumericVector nneighbours, double tau, double lik_var,
                        const NumericVector y, double alpha, NumericVector XB){
  
  // stuff associated with phiVarb update
  int    n;
  double meanFullCon;
  double varFullCon;
  double sumphiVarb;
  double priorvardenom;
  double priormean;
  double priorvar;
  double futuresumphiVarb;
  double pastsumphiVarb;
  double asqOne = 1 + pow(alpha, 2); 
  NumericVector phiVarb = clone(phi);
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
        if(ntimes > 1) {
          if((0 < i) && (i < (ntimes-1))){
            priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past and the futuren
            for(int l = 0; l < nsites; l++)  {
              futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
              pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            }
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
            priorvar = tau/(priorvardenom*asqOne);
          } else if(i == 0) {
            priormean = phiVarb[((i+1)*nsites) + j];
            //calculate the sum of the neighbours of j in the future
            for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
            priorvar =  tau/(priorvardenom*asqOne);
          } else if(i == (ntimes-1)) {
            priormean = phiVarb[((i-1)*nsites) + j];
            //calculate the sum of the neighbours of j in the past
            for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
            priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
            priorvar = tau/(priorvardenom);
          }
        } else if(ntimes == 1){
          priorvar  = tau/priorvardenom;
          priormean = sumphiVarb/priorvardenom; 
        }
      
      // get the mean and variance of the full conditional distribution
      varFullCon  = 1/((1/priorvar) + (1/lik_var));
      meanFullCon = ((priormean/priorvar) + (y[n] - XB[n])/lik_var)*varFullCon;
      phiVarb[n]  = rnorm(1, meanFullCon, sqrt(varFullCon))[0];    
    }
  }
  List out(2);
  out[0] = phiVarb;
  return out;
}




// [[Rcpp::export]]
double qform_difference_ST(NumericMatrix Qtrip, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qtrip.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    if( !(Qtrip(i, 2) == 0) ){
      spRow  = Qtrip(i, 0);
      spCol  = Qtrip(i, 1);
      for(int j = 0; j < nrowTime; j++){
        tiRow  = Qtime(j, 0);
        tiCol  = Qtime(j, 1);
        stRow  = (tiRow - 1)*nsites + spRow;
        stCol  = (tiCol - 1)*nsites + spCol;
        Qform += phi[stCol - 1] * Qtrip(i, 2) * Qtime(j, 2) * phi[stRow - 1];
      }
    }
  }
  return(Qform);
}




// [[Rcpp::export]]
double qform_ST(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi[stRow - 1];
    }
  }
  return(Qform);
}




// [[Rcpp::export]]
double qform_ST_asym(NumericMatrix Qspace, NumericMatrix Qtime, NumericVector phi1, NumericVector phi2, int nsites){ 
  int nrowSpace = Qspace.nrow();
  int nrowTime  = Qtime.nrow();
  int spRow, spCol, tiRow, tiCol, stRow, stCol;
  double Qform = 0;
  for(int i = 0; i < nrowSpace; i++){
    spRow  = Qspace(i, 0);
    spCol  = Qspace(i, 1);
    for(int j = 0; j < nrowTime; j++){
      tiRow  = Qtime(j, 0);
      tiCol  = Qtime(j, 1);
      stRow  = (tiRow - 1)*nsites + spRow;
      stCol  = (tiCol - 1)*nsites + spCol;
      Qform += phi1[stCol - 1] * Qspace(i, 2) * Qtime(j, 2) * phi2[stRow - 1];
    }
  }
  return(Qform);
}




// [[Rcpp::export]]
NumericMatrix update_Qtime(NumericMatrix Qtime, double alpha, int rowNumberLastDiag){
  int nr = Qtime.nrow();
  double alphasq = alpha*alpha;
  NumericMatrix Qtime_new = clone(Qtime);
  for(int i = 0; i < nr; i++){
    if(Qtime(i, 0) == Qtime(i, 1))  Qtime_new(i, 2) = 1 + alphasq;
    if(!(Qtime(i, 0) == Qtime(i, 1)) ) Qtime_new(i, 2) = -alpha;
  }
  Qtime_new(rowNumberLastDiag,2) = 1;
  return Qtime_new;
}




// [[Rcpp::export]]
NumericMatrix updatetriplets_rho(NumericMatrix trips, int nsites, double rho_old, double rho_new, double fixedridge){    
  //   create a clone of the triplet matrix that will be used for output              
  NumericMatrix tripsnew = clone(trips);   
  int rows               = tripsnew.nrow();
  double one_rho_old     = 1 - rho_old;
  double one_rho_new     = 1 - rho_new;
  //   loop over all of the diagonal elements first.  To update the diagonal elements, 
  //   we first substract (1 - rho_old) and the fixed ridge and divide by rho_old.  Then multiply by rho_new
  //   and add 1 - rho_new.
  for(int i = 0; i < nsites; i++) {
    tripsnew(i, 2) = ((trips(i, 2) - one_rho_old - fixedridge)/rho_old)*rho_new + one_rho_new + fixedridge;
  }
  //   loop over all of the off-diagonal elements.  These are updated by simply 
  //   dividing by rho_old and multiplying by rho_new
  for(int j = nsites; j < rows; j++) {
    tripsnew(j, 2)     = (trips(j, 2)/rho_old)*rho_new;
  }
  //   output the updated triplets
  return tripsnew;
}




// [[Rcpp::export]]
List updatetripList2(NumericMatrix trips, NumericVector vold, 
                     NumericVector vnew, const int nedges, int nsites, IntegerVector block, int block_length, double rho, 
                     double fixedridge){                   
  //   create a clone of the triplet matrix                   
  NumericMatrix temporary = clone(trips); 
  NumericMatrix difference = clone(trips);
  for(int l = 0; l < temporary.nrow(); l++) difference(l, 2) = 0;
  //   stuff needed for intermediate calculations
  double oldoffdiag_binary, newoffdiag_binary;
  double newoffdiag_logit,  oldoffdiag_logit;
  int    diag_inds_1,       diag_inds_2, i;
  double olddiag_binary_1,  olddiag_binary_2;
  double newdiag_binary_1,  newdiag_binary_2;
  //   loop over all edges, stop only at elements where vold and vnew differ
  //   then perform and update of the corresponding triplets
  for(int j = 0; j < block_length; j++) {  
    i = block[j] - 1;
    //       this is the old off diagonal element (-ve to make +ve)
    oldoffdiag_binary = -temporary(i + nsites, 2)/rho;
    //       old and new v elements
    oldoffdiag_logit  = vold[i];
    newoffdiag_logit  = vnew[i];
    //       convert new v onto [0,1] scale
    newoffdiag_binary = -(1/(1 + exp(-newoffdiag_logit)));
    //       replace triplets with new v
    difference(i + nsites, 2) = rho*(temporary(i + nsites, 2) - newoffdiag_binary);
    temporary(i + nsites, 2)  = rho*newoffdiag_binary;
    difference(i + nsites + nedges, 2) = rho*(temporary(i + nsites + nedges, 2) - newoffdiag_binary);
    temporary(i + nsites + nedges, 2)  = rho*(newoffdiag_binary);
    
    //       now need to find x and y coords of these offdiags to update diags
    diag_inds_1 = temporary(i + nsites, 0) - 1;
    diag_inds_2 = temporary(i + nsites, 1) - 1;
    //       get the old binary elements
    olddiag_binary_1 = (temporary(diag_inds_1, 2) - fixedridge - (1 - rho))/rho;
    olddiag_binary_2 = (temporary(diag_inds_2, 2) - fixedridge - (1 - rho))/rho;
    //       calculate and replace with new ones
    newdiag_binary_1 = (olddiag_binary_1 - oldoffdiag_binary)  - newoffdiag_binary;
    newdiag_binary_2 = (olddiag_binary_2 - oldoffdiag_binary)  - newoffdiag_binary;
    temporary(diag_inds_1, 2)  = (newdiag_binary_1*rho) + fixedridge + (1 - rho);
    temporary(diag_inds_2, 2)  = (newdiag_binary_2*rho) + fixedridge + (1 - rho);
    difference(diag_inds_1, 2) = trips(diag_inds_1, 2) -  temporary(diag_inds_1, 2);
    difference(diag_inds_2, 2) = trips(diag_inds_2, 2) -  temporary(diag_inds_2, 2);
  }
  //   output to a list where the first element is the updated diagonal triplets
  //   second element is the offdiagonal triplets
  List out(2);
  out[0] = temporary;
  out[1] = difference;
  return out;
}






// [[Rcpp::export]]
List SPTICARphiBinomial(NumericMatrix W, const int nsites, const int ntimes, 
                        NumericVector phi, NumericVector nneighbours, double tau,
                        const NumericVector y, double alpha, NumericVector XB, 
                        const double phiVarb_tune, NumericVector trials){
  
  // stuff associated with phiVarb update
  int    n;
  double theta_prop;
  double theta;
  double l_prob_new;
  double l_prob_old;
  double sumphiVarb;
  double priorvardenom;
  double priormean;
  double priorvar;
  double phi_new;
  double futuresumphiVarb;
  double pastsumphiVarb;
  double acceptance;
  NumericVector accept(4);
  double asqOne = 1 + pow(alpha, 2); 
  NumericVector phiVarb = clone(phi);
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
      if(ntimes > 1) {
        if((0 < i) && (i < (ntimes-1))){
          priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past and the futuren
          for(int l = 0; l < nsites; l++)  {
            futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          }
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
          priorvar = tau/(priorvardenom*asqOne);
        } else if(i == 0) {
          priormean = phiVarb[((i+1)*nsites) + j];
          //calculate the sum of the neighbours of j in the future
          for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
          priorvar =  tau/(priorvardenom*asqOne);
        } else if(i == (ntimes-1)) {
          priormean = phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past
          for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
          priorvar = tau/(priorvardenom);
        }
      } else if(ntimes == 1){
        priorvar  = tau/priorvardenom;
        priormean = sumphiVarb/priorvardenom; 
      }
      
      // get the theta parameter in [0,1] of the binomial under proposal and old 
      phi_new      = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
      theta        = 1 + exp(-XB[n] - phiVarb[n]);
      theta_prop   = 1 + exp(-XB[n] - phi_new);
      l_prob_old   = -y[n]*log(theta) + (trials[n] - y[n])*log(1 - (1/theta)) - (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
      l_prob_new   = -y[n]*log(theta_prop) + (trials[n] - y[n])*log(1 - (1/theta_prop)) - (0.5/priorvar) * pow(phi_new - priormean, 2);
      acceptance   = exp(l_prob_new - l_prob_old);
      if(acceptance >= 1){
        phiVarb[n] = phi_new;
        accept[1]++;
      } else {
        if(runif(1)[0] <= acceptance) {
          phiVarb[n] = phi_new;
          accept[1]++;
        }
      }
    }
  }
  List out(2);
  out[0] = accept;
  out[1] = phiVarb;
  return out;
}




// [[Rcpp::export]]
List SPTICARphiVarb(NumericMatrix W, const int nsites, const int ntimes, 
                    NumericVector phiVarb, NumericVector nneighbours, double tau, 
                    const NumericVector y, const NumericVector E, 
                    const double phiVarb_tune, double alpha, NumericVector XB){
  
  //stuff associated with model objects
  int n;
  int k = ntimes*nsites;
  
  //General MCMC things
  NumericVector accept(4);
  double acceptance;
  
  // stuff associated with phiVarb update
  double newpriorbit;
  double oldpriorbit;
  double sumphiVarb;
  double gubbins;
  double priorvardenom, priormean;
  double priorvar;
  double propphiVarb;
  double futuresumphiVarb, pastsumphiVarb;
  double asqOne = 1 + pow(alpha, 2); 
  
  // stuff associated with tau update
  NumericVector arc(k);
  //Function rtrunc("rtrunc");
  
  //  outer loop for each random effect
  for(int i = 0; i < ntimes; i++) { 
    for(int j = 0; j < nsites; j++){
      sumphiVarb = 0;
      futuresumphiVarb = 0;
      pastsumphiVarb = 0;
      n = (i*nsites) + j;
      //  calulate the sum of the neighbours of j at time i
      for(int l = 0; l < nsites; l++)  sumphiVarb +=  W(j, l) * phiVarb[i*nsites + l]; 
      // calculate prior variance
      priorvardenom = nneighbours[j];
      // calculate prior mean of phiVarb_jt: different for t=1, t=n and else
      if(ntimes > 1) {
        if((0 < i) && (i < (ntimes-1))){
          priormean = phiVarb[((i+1)*nsites) + j] + phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past and the futuren
          for(int l = 0; l < nsites; l++)  {
            futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
            pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          }
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb + alpha*pastsumphiVarb);
          priorvar = tau/(priorvardenom*asqOne);
        } else if(i == 0) {
          priormean = phiVarb[((i+1)*nsites) + j];
          //calculate the sum of the neighbours of j in the future
          for(int l = 0; l < nsites; l++)  futuresumphiVarb +=  W(j, l) * phiVarb[(i+1)*nsites + l];
          priormean = (alpha*priormean)/asqOne - (1/(priorvardenom*asqOne))*(alpha*futuresumphiVarb - asqOne*sumphiVarb);
          priorvar =  tau/(priorvardenom*asqOne);
        } else if(i == (ntimes-1)) {
          priormean = phiVarb[((i-1)*nsites) + j];
          //calculate the sum of the neighbours of j in the past
          for(int l = 0; l < nsites; l++)  pastsumphiVarb +=  W(j, l) * phiVarb[(i-1)*nsites + l];
          priormean = (alpha*priormean) - (1/(priorvardenom))*(alpha*pastsumphiVarb - sumphiVarb);
          priorvar = tau/(priorvardenom);
        }
      } else if(ntimes == 1){
        priorvar = tau/priorvardenom;
        priormean = 1*sumphiVarb/priorvardenom; 
      }
      
      // propose a value and accept or reject it 
      propphiVarb = rnorm(1, phiVarb[n], sqrt(priorvar*phiVarb_tune))[0];    
      newpriorbit = (0.5/priorvar) * pow(propphiVarb - priormean, 2); 
      oldpriorbit = (0.5/priorvar) * pow(phiVarb[n] - priormean, 2);
      gubbins = exp(E[n] + XB[n])*(exp(propphiVarb) - exp(phiVarb[n]));
      acceptance = exp(y[n]*(propphiVarb - phiVarb[n]) - newpriorbit + oldpriorbit - gubbins);
      arc[n] = acceptance;
      if(acceptance >= 1){
        phiVarb[n] = propphiVarb;
        accept[1]++;
      } else {
        if(runif(1)[0] <= acceptance) {
          phiVarb[n] = propphiVarb;
          accept[1]++;
        }
      }
    }
  }
  List out(3);
  out[0] = accept;
  out[1] = phiVarb;
  out[2] = arc;
  return out;
}






////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Localised model functions   //////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericMatrix Zupdatesqbin(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, NumericMatrix failures)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G), lp(G), prob(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,0) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,0) * log(prob) + failures(k,0) * log((1 - prob));
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }




//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        lp = Offset(k,j) + lambda;
        prob = exp(lp) / (1 + exp(lp));
        like1 = Y(k,j) * log(prob) + failures(k,j) * log((1 - prob));
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
    }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    lp = Offset(k,ntimeminus) + lambda;
    prob = exp(lp) / (1 + exp(lp));
    like1 = Y(k,ntimeminus) * log(prob) + failures(k,ntimeminus) * log((1 - prob));
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }


return Z;
}




// [[Rcpp::export]]
NumericMatrix Zupdatesqpoi(NumericMatrix Z, NumericMatrix Offset, NumericMatrix Y, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   


// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,0) * lambda - exp(lambda) * Offset(k,0);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = Y(k,j) * lambda - exp(lambda) * Offset(k,j);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = Y(k,ntimeminus) * lambda - exp(lambda) * Offset(k,ntimeminus);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }


return Z;
}





// [[Rcpp::export]]
NumericMatrix Zupdatesqgau(NumericMatrix Z, NumericMatrix Offset, const double delta, 
NumericVector lambda, const int nsites, const int ntime, const int G, NumericVector SS, NumericVector prioroffset,
const double Gstar, const double nu2)
{
// Quantites needed for the MCMC updating
NumericVector like1(G), prior1(G), prior2(G), posterior1(G), posterior2(G), posterior3(G);
NumericVector prior3(G), prior4(G);

// Elements to undertake the posterior sampling
double u;
double cs;
int test;
int counter;
int ntimeminus = ntime -1;     
int ntimeminus2 = ntime -2;   



// Update the elements in Z during time period 1
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,0) - lambda),2) / (2*nu2);
    prior1 = -delta * pow(SS - Z(k,1),2) - prioroffset;
    prior2 = -delta * pow(SS - Gstar,2);
    posterior1 = like1 + prior1 + prior2;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,0) = counter;
    }


//  Update the elements in Z during the remaining time periods except time ntime
    for(int j = 1; j < ntimeminus; j++)
    {
        for(int k = 0; k < nsites; k++) 
        {
        // Compute the full conditional
        like1 = -pow((Offset(k,j) - lambda),2) / (2*nu2);
        prior1 = -delta * pow(SS - Z(k,(j+1)),2) - prioroffset;
        prior2 = -delta * (pow(SS - Z(k,(j-1)),2) + pow(SS - Gstar,2));
        posterior1 = like1 + prior1 + prior2;
        posterior2 = posterior1 - max(posterior1);
        posterior3 = exp(posterior2) / sum(exp(posterior2));
    
        // Sample a new Z value
        u = runif(1)[0];          
        cs = posterior3[0];
        test = 1;
        counter = 1;
            while(test==1)
            {
                if(cs>u)
                {
                test = 0;     
                }else
                {
                counter = counter + 1;
                cs = cs + posterior3[(counter-1)];     
                }
            }
        Z(k,j) = counter;
        }
    }


// Update the elements in Z during time period ntime (the last one)
    for(int k = 0; k < nsites; k++) 
    {
    // Compute the full conditional
    like1 = -pow((Offset(k,ntimeminus) - lambda),2) / (2*nu2);
    prior1 = -delta * (pow(SS - Z(k,ntimeminus2),2) + pow(SS - Gstar,2));
    posterior1 = like1 + prior1;
    posterior2 = posterior1 - max(posterior1);
    posterior3 = exp(posterior2) / sum(exp(posterior2));
    
    // Sample a new Z value
    u = runif(1)[0];          
    cs = posterior3[0];
    test = 1;
    counter = 1;
        while(test==1)
        {
            if(cs>u)
            {
            test = 0;     
            }else
            {
            counter = counter + 1;
            cs = cs + posterior3[(counter-1)];     
            }
        }
    Z(k,ntimeminus) = counter;
    }




return Z;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// Sepspatial model functions   ////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector tau2compute(NumericVector tau2, NumericVector temp, const double tau2_shape, const double prior_tau2, const int N)
{
  NumericVector tau2new;
  tau2new = tau2;
  double tau2_scale;
  
  for (int t = 0; t < N; t++)
  {
    tau2_scale = temp[t] + prior_tau2;
    tau2new[t] = 1 / rgamma(1, tau2_shape, (1/tau2_scale))[0];
  }
  return tau2new;
}






// [[Rcpp::export]]
double rhoquadformcompute(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                          const int nsites, const int ntime, NumericMatrix phi, double rho, NumericVector tau2)
{    
  NumericVector temp(nsites);
  double num=0;
  
  for(int t = 0; t < ntime; t++)
  {
    temp = phi(_,t);  
    num = num + (quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho) / tau2[t]);
  }
  return num;
}





// [[Rcpp::export]]
List binomialsrecarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                         NumericMatrix phi, double rho, const NumericMatrix ymat, const NumericMatrix failuresmat,
                         const double phi_tune, NumericMatrix offset, NumericVector denoffset, NumericVector tau2)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew, pold, pnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  
  //////////////////////////////////////////////////////
  // Update the random effects at all time points
  //////////////////////////////////////////////////////
  for(int t = 0; t < ntime; t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j];
      priorvar = tau2[t] / priorvardenom;
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * phinew(row,t); 
        priormeantemp = priormeantemp + temp; 
      }
      priormean = (rho * priormeantemp) / priorvardenom;   
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      pold = exp(lpold) / (1 + exp(lpold));
      pnew = exp(lpnew) / (1 + exp(lpnew));        
      oldlikebit = ymat(j,t) * log(pold) + failuresmat(j,t) * log((1-pold));
      newlikebit = ymat(j,t) * log(pnew) + failuresmat(j,t) * log((1-pnew));
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}





// [[Rcpp::export]]
List poissonsrecarupdateRW(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntime,
                        NumericMatrix phi, double rho, const NumericMatrix ymat, const double phi_tune, NumericMatrix offset, NumericVector denoffset,
                        NumericVector tau2)
{    
  ///////////////////////////////////////////    
  // Specify variables needed in the function
  ///////////////////////////////////////////
  double temp, priormean, priormeantemp, priorvar, priorvardenom, acceptance;
  double propphi, oldpriorbit, newpriorbit, oldlikebit, newlikebit, lpold, lpnew; 
  NumericMatrix phinew(nsites,ntime);
  phinew = phi;
  int row, rowstart, rowend, accept=0;
  
  //////////////////////////////////////////////////////
  // Update the random effects at all time points
  //////////////////////////////////////////////////////
  for(int t = 0; t < ntime; t++)
  {
    for(int j = 0; j < nsites; j++)
    {
      // calculate prior mean and variance
      priorvardenom = denoffset[j];
      priorvar = tau2[t] / priorvardenom;
      rowstart = Wbegfin(j,0) - 1;
      rowend = Wbegfin(j,1);
      priormeantemp = 0;
      for(int l = rowstart; l < rowend; l++)
      {
        row = Wtriplet(l,1) - 1;
        temp = Wtriplet(l, 2) * phinew(row,t); 
        priormeantemp = priormeantemp + temp; 
      }
      priormean = (rho * priormeantemp) / priorvardenom;   
      
      // Propose a value and calculate the acceptance probability
      propphi = rnorm(1, phinew(j,t), sqrt(priorvar*phi_tune))[0];
      newpriorbit = (0.5/priorvar) * pow((propphi - priormean), 2); 
      oldpriorbit = (0.5/priorvar) * pow((phinew(j,t) - priormean), 2);
      lpold = phinew(j,t) + offset(j, t);
      lpnew = propphi + offset(j, t); 
      oldlikebit = ymat(j,t) * lpold - exp(lpold);
      newlikebit = ymat(j,t) * lpnew - exp(lpnew);
      acceptance = exp(oldpriorbit - newpriorbit - oldlikebit + newlikebit);
      if(runif(1)[0] <= acceptance) 
      {
        phinew(j,t) = propphi;
        accept = accept + 1;
      }
      else
      { 
      }
    }
  }
  
  List out(2);
  out[0] = phinew;
  out[1] = accept;
  return out;
}



// [[Rcpp::export]]
NumericVector tauquadformcompute2(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                          const int nsites, const int ntime, NumericMatrix phi, double rho)
{    
NumericVector temp(nsites), num(ntime);

// Compute the sum of quadratic forms for updating tau
for(int t = 0; t < ntime; t++)
{
  temp = phi(_,t);  
  num[t] = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rho);
}

return num;
}





////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Clusttrends functions   /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
NumericVector tempupdate(const int Nchains, double dt)
{
    //Create new objects
    NumericVector newtemps(Nchains);
    
    newtemps[0] = 1;
    
    for(int i = 1; i < Nchains; i++)
    {
        newtemps[i] = dt * newtemps[(i-1)];
    }
    // Return the result
    return newtemps;
}


// [[Rcpp::export]]
NumericMatrix matcomp(NumericMatrix X, NumericMatrix beta, NumericVector prop, const int p, const int Nchains)
{
    //Create new objects
    NumericVector newprop = clone(prop);
    NumericMatrix proposal(p, Nchains), newbeta = clone(beta), newX = clone(X);
    NumericVector gen(p);
    NumericVector matmult(p);
    
    Environment base("package:stats");
    Function sample = base["rnorm"];
    
    for(int i = 0; i < Nchains; i++)
    {
        gen = as<NumericVector>( rnorm(p, 0, 1) );
        
        for(int j = 0; j < p; j++)
        {
            matmult[j] = sum((sqrt(newprop[i]) * newX(j, _)) * gen);
        }
        proposal(_, i) = newbeta(_, i) +  matmult;
    }
    // Return the result
    return proposal;
}


// [[Rcpp::export]]
NumericMatrix offsetcompute(NumericMatrix w, NumericMatrix gamma, NumericMatrix time, const int Nchains, const int nsites, const int Ntrends, NumericVector begin)
{
    //Create new objects
    NumericMatrix wnew = clone(w), gammanew = clone(gamma), timenew = clone(time);
    NumericMatrix offset(nsites, Nchains);
    int rowstart;
    
    for(int i = 0; i < nsites; i++)
    {
        for(int j = 0; j < Ntrends; j++)
        {
            rowstart = begin[j] - 1;
            
            offset(i, _) = offset(i, _) + wnew((rowstart + i), _) * (gammanew((rowstart + i), _) * timenew((rowstart + i), _));
        }
        
    }
    
    return offset;
}


// [[Rcpp::export]]
NumericMatrix matN(NumericVector x, const int nsites, const int Nchains)
{
    //Create new objects
    NumericMatrix mat(nsites, Nchains);
    NumericVector newx = clone(x);
    
    for(int j = 0; j < nsites; j++)
    {
        mat(j, _) = newx; 
    }
    // Return the result
    return mat;
}  


// [[Rcpp::export]]
NumericMatrix linpredcomputeNchains(NumericMatrix X, const int nsites, const int p, NumericMatrix beta, const int Nchains)
{
    // Create new objects
    // Compute the linear predictor
    NumericMatrix linpred(nsites, Nchains);
    double temp;
    
    for(int j = 0; j < Nchains; j++)
    {
        // Compute the linear predictor via a double for loop
        for(int k = 0; k < nsites; k++)
        {
            temp = 0;
            
            for(int l = 0; l < p; l++) temp = temp + X(k, l) * beta(l, j);
            
            linpred(k, j) = temp;
        }
    }
    
    // Return the result
    return linpred;
}


// [[Rcpp::export]]
NumericVector gammaproposal(const int Nchains, NumericVector gamma, NumericVector gamma_tune, const int prior_vargamma, NumericVector Wareas, const int trend, const int knots)
{
    NumericVector proposal(Nchains), gammanew = clone(gamma);
    NumericVector Wareasnew = clone(Wareas);
    double inf = std::numeric_limits<double>::infinity();
    
    Environment base("package:truncdist");
    Function rtrunc = base["rtrunc"];
    
    for(int j = 0; j < Nchains; j++)
    {
        if(Wareasnew[j] == 0)
        {
            
            if((trend == 2) | (trend == 5) | (trend == 6) | ((trend >= 8) & (trend <= (8 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, 0, 0.01) );
            }else if((trend == 3) | (trend == 4) | (trend == 7) | ((trend >= (8 + knots + 1)) & (trend <= (8 + knots + 1 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, 0, 0.01) );
            }else
            {}
        }
        else
        {
            if((trend == 2) | (trend == 5) | (trend == 6) | ((trend >= 8) & (trend <= (8 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", -inf, 0, gammanew[j], gamma_tune[j]) );
            }else if((trend == 3) | (trend == 4) | (trend == 7) | ((trend >= (8 + knots + 1)) & (trend <= (8 + knots + 1 + knots))))
            {
                proposal[j] = as<double>( rtrunc(1, "norm", 0, inf, gammanew[j], gamma_tune[j]) );
            }else
            {}
        }
    }
    return proposal;
}


// [[Rcpp::export]]
NumericMatrix lambdaupdate(const int Nchains, NumericMatrix temp)
{
    Environment base("package:gtools");
    Function rdirichlet = base["rdirichlet"];
    //Create new objects
    NumericMatrix tempnew = clone(temp), lambdanew = clone(temp);
    
    for(int j = 0; j < Nchains; j++)
    {
        lambdanew(_, j) = as<NumericVector>( rdirichlet(1, tempnew(_, j)) );
    }
    return lambdanew;
}


// [[Rcpp::export]]
NumericVector tau2quadform(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, const int nsites, 
                           NumericMatrix phi, NumericMatrix theta, NumericVector rho, const int Nchains)
{
    // Compute a quadratic form for the random effects
    // Create new objects 
    NumericVector tau2_posteriorscale(Nchains), rhonew = clone(rho);
    NumericMatrix phinew = clone(phi), thetanew = clone(theta);
    double tau2_quadform = 0, tau2_phisq = 0;
    int row, col;
    
    for(int j = 0; j < Nchains; j++)
    {
        tau2_quadform = 0;
        tau2_phisq = 0;
        // Compute the off diagonal elements of the quadratic form
        for(int i = 0; i < n_triplet; i++)
        {
            row = Wtriplet(i, 0) - 1;
            col = Wtriplet(i, 1) - 1;
            tau2_quadform += phinew((Wtriplet(i, 0) - 1), j) * thetanew((Wtriplet(i, 1) - 1), j) * Wtriplet(i, 2); 
        }
        // Compute the diagonal elements of the quadratic form          
        for(int l = 0; l < nsites; l++)
        {
            tau2_phisq += phinew(l, j) * thetanew(l, j) * (rhonew[j] * Wtripletsum[l] + 1 - rhonew[j]);    
        }
        // Compute the quadratic form
        tau2_posteriorscale[j] = 0.5 * (tau2_phisq - rhonew[j] * tau2_quadform);
    }
    // Return the simulated values
    return tau2_posteriorscale;
}


// [[Rcpp::export]]
NumericVector tau2computeNchains(NumericVector temp, const double tau2_shape, const double prior_tau2, const int Nchains)
{
    NumericVector tau2new(Nchains);
    double tau2_scale;
    
    for (int j = 0; j < Nchains; j++)
    {
        tau2_scale = temp[j] + prior_tau2;
        tau2new[j] = 1 / rgamma(1, tau2_shape, (1/tau2_scale))[0];
    }
    return tau2new;
}


// [[Rcpp::export]]
NumericVector rhoquadformcomputeNchains(NumericMatrix Wtriplet, NumericVector Wtripletsum, const int n_triplet, 
                                        const int nsites, const int Nchains, NumericMatrix phi,
                                        NumericVector rho, NumericVector tau2)
{    
    NumericVector temp(nsites);
    NumericVector rho_quadform(Nchains);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix phinew = clone(phi);
    
    for(int j = 0; j < Nchains; j++)
    {
        temp = phinew(_, j);  
        rho_quadform[j] = quadform(Wtriplet, Wtripletsum, n_triplet, nsites, temp, temp, rhonew[j]) / tau2new[j];
    }
    return rho_quadform;
}


// [[Rcpp::export]]
NumericVector Qdet(const int Nchains, NumericVector rho, NumericVector Wstar_val)
{
    NumericVector detQ(Nchains);
    NumericVector rhonew = clone(rho);
    
    for(int j = 0; j < Nchains; j++)
    {
        detQ[j] = 0.5 * sum(log((rhonew[j] * Wstar_val + (1 - rhonew[j]))));
    }
    return detQ;
}


// [[Rcpp::export]]
List poissondevfit(NumericVector y, NumericMatrix fitted, const int nsites, const int Nchains)
{
    NumericVector newy = clone(y);
    NumericMatrix newfitted = clone(fitted), like_all(nsites, Nchains);
    NumericVector deviance(Nchains), deviance_all(nsites), fit(nsites);
    
    Environment base2("package:stats");
    Function dpois = base2["dpois"];
    
    for(int j = 0; j < Nchains; j++)
    {
        fit = newfitted(_, j);
        deviance_all = as<NumericVector>( dpois(newy, fit) );
        like_all(_, j) = deviance_all;
        deviance[j] = -2 * sum(log(deviance_all));
    }
    
    List out(2);
    out[0] = deviance;
    out[1] = like_all;
    return out;
}


// [[Rcpp::export]]
NumericVector poissonbetablockupdate(const int nsites, NumericMatrix beta, NumericMatrix betaprop, NumericMatrix lp_beta, NumericMatrix lp_betaprop,
                                     NumericMatrix offset, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                                     const int Nchains, NumericVector temps, const int p)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), acceptance(Nchains);
    NumericMatrix betanew = clone(beta), betapropnew = clone(betaprop);
    NumericMatrix newlp_beta = clone(lp_beta), newlp_betaprop = clone(lp_betaprop), newoffset = clone(offset);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        priorbit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newlp_beta(i, j) + newoffset(i, j);
            lp_proposal[i] = newlp_betaprop(i, j) + newoffset(i, j);
            
            p_current[i] = exp(lp_current[i]);
            p_proposal[i] = exp(lp_proposal[i]);
            
            oldlikebit += y[i] * lp_current[i] - p_current[i];
            newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
        }
        likebit = newlikebit - oldlikebit;
        
        // Create the prior acceptance component
        for(int k = 0; k < p; k++)
        {
            priorbit += 0.5 * pow((betanew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betapropnew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        }
        
        // Compute the acceptance probability and return the value
        acceptance[j] = exp((likebit + priorbit) * temps[j]);
    }
    
    return acceptance;
}


// [[Rcpp::export]]
List poissongammaupdate(const int nsites, NumericVector gamma, NumericVector proposal, NumericMatrix offset, NumericMatrix offset_proposal, NumericVector y,
                        double prior_meangamma, double prior_vargamma, const int Nchains, NumericVector temps)
{
    // Compute the acceptance probability for gamma
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
    NumericVector accept(Nchains);
    NumericVector gammanew = clone(gamma), proposalnew = clone(proposal);
    NumericMatrix newoffset = clone(offset), newproposal = clone(offset_proposal);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            
            lp_current[i] = newoffset(i, j);
            lp_proposal[i] = newproposal(i, j);
            
            p_current[i] = exp(lp_current[i]);
            p_proposal[i] = exp(lp_proposal[i]);
            
            oldlikebit += y[i] * lp_current[i] - p_current[i];
            newlikebit += y[i] * lp_proposal[i] - p_proposal[i];
        }
        likebit = newlikebit - oldlikebit;
        
        priorbit = 0.5 * pow((gammanew[j] - prior_meangamma), 2) / prior_vargamma - 0.5 * pow((proposalnew[j] - prior_meangamma), 2) / prior_vargamma;
        
        // Compute the acceptance probability and return the value
        acceptance = exp((likebit + priorbit) * temps[j]);
        if(runif(1)[0] <= acceptance)
        {
            gammanew[j] = proposalnew[j];
            accept[j] = accept[j] + 1;
        }
        else
        {
        }
    }
    
    List out(2);
    out[0] = gammanew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List poissonwupdate(const int nsites, const int ntimes, NumericMatrix w, NumericMatrix offset, NumericMatrix offset_proposal, NumericMatrix w_proposal,
                    NumericMatrix y, NumericMatrix lambda, const int Nchains, NumericVector temps, NumericVector begin, NumericVector regbegin, const int Ntrends)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix wnew = clone(w), wproposal = clone(w_proposal), lambdanew = clone(lambda), newoffset = clone(offset), newproposal = clone(offset_proposal);
    int rowstart, regstart;
    NumericVector proposal;
    
    for(int j = 0; j < Nchains; j++)
    {
        rowstart = begin[j] - 1;
        
        for(int k = 0; k < nsites; k++)
        {
            proposal = wproposal((rowstart + k), _);
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j);
                lp_proposal[t] = newproposal((regstart + k), j);
                
                p_current[t] = exp(lp_current[t]);
                p_proposal[t] = exp(lp_proposal[t]);
                
                oldlikebit += y(k, t) * lp_current[t] - p_current[t];
                newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
            }
            
            likebit = newlikebit - oldlikebit;
            
            for(int i = 0; i < Ntrends; i++)
            {
                priorbit += proposal[i] * log(lambdanew(i, j)) - wnew((rowstart + k), i) * log(lambdanew(i, j));
            }
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                wnew((rowstart + k), _) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = wnew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List poissonphiupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntimes,
                      NumericMatrix phi, NumericMatrix offset, NumericMatrix y, NumericVector tau2, NumericVector rho, const int Nchains,
                      NumericVector temps, NumericMatrix phi_tune, NumericVector regbegin)
{
    // Compute the acceptance probability
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix phinew = clone(phi), newoffset = clone(offset);
    double proposal;
    double priorvardenom, priormean, priorvar, sumphi;
    int rowstart=0, rowend=0;
    int regstart;
    
    for(int j = 0; j < Nchains; j++)
    {
        
        for(int k = 0; k < nsites; k++)
        {
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            // Calculate prior variance
            priorvardenom = rhonew[j] * Wtripletsum[k] + 1 - rhonew[j];
            priorvar = tau2new[j] / priorvardenom;
            
            // Calculate the prior mean
            rowstart = Wbegfin(k, 0) - 1;
            rowend = Wbegfin(k, 1);
            sumphi = 0;
            for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), j);
            priormean = rhonew[j] * sumphi / priorvardenom; 
            
            // propose a value 
            proposal = rnorm(1, phinew(k, j), sqrt(priorvar*phi_tune(k, j)))[0];
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j) + phinew(k, j);
                lp_proposal[t] = newoffset((regstart + k), j) + proposal;
                
                p_current[t] = exp(lp_current[t]);
                p_proposal[t] = exp(lp_proposal[t]);
                
                oldlikebit += y(k, t) * lp_current[t] - p_current[t];
                newlikebit += y(k, t) * lp_proposal[t] - p_proposal[t];
            }
            
            likebit = newlikebit - oldlikebit;
            priorbit = (0.5/priorvar) * pow((phinew(k, j) - priormean), 2) - (0.5/priorvar) * pow((proposal - priormean), 2);
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(k, j) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
int poissoncouplingAllupdate(const int nsites, const int K, const int p, NumericMatrix w, NumericMatrix offset, NumericMatrix beta, NumericMatrix gamma,
                             NumericMatrix lambda, NumericMatrix phi, NumericVector rho, NumericVector tau2, NumericVector Wtripletsum,
                             NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector y, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                             NumericVector prior_meantrends, NumericVector prior_vartrends, NumericVector prior_lambda, NumericVector prior_tau2,
                             NumericVector swap, NumericVector temps, NumericVector begin, const int Ntrends, const int TrendSel)
{
    // Compute the acceptance probability for swapping
    //Create new objects
    double acceptance, chain1oldlikebit=0, chain2oldlikebit=0, chain1newlikebit=0, chain2newlikebit=0, chain1likebit, chain2likebit, chain1priorbitbeta=0, chain2priorbitbeta=0;
    double wpriorchain1old=0, wpriorchain1new=0, wpriorchain2old=0, wpriorchain2new=0, chain1priorbitw, chain2priorbitw, chain1priorbitlambda, chain2priorbitlambda;
    double phichain1=0, phichain2=0;
    double chain1priorbitgamma=0, chain2priorbitgamma=0;
    double chain1priorbittau2, chain2priorbittau2, chain1priorbitphi, chain2priorbitphi;
    NumericVector lpchain1_current(nsites), lpchain1_proposal(nsites), pchain1_current(nsites), pchain1_proposal(nsites);
    NumericVector lpchain2_current(nsites), lpchain2_proposal(nsites), pchain2_current(nsites), pchain2_proposal(nsites);
    int accept = 0;
    NumericVector swapnew;
    swapnew = swap - 1;
    double swap1, swap2;
    swap1 = swapnew[0];
    swap2 = swapnew[1];
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix gammanew = clone(gamma), wnew = clone(w), newoffset = clone(offset), betanew = clone(beta), lambdanew = clone(lambda), phinew = clone(phi);
    int beginchain1 = begin[swap1] - 1;
    int beginchain2 = begin[swap2] - 1;
    double temp1, temp2;
    temp1 = temps[swap1];
    temp2 = temps[swap2];
    double priorvardenomchain1, priorvardenomchain2, phipriorvarchain1, phipriorvarchain2, sumphichain1=0, sumphichain2=0, phipriormeanchain1, phipriormeanchain2;
    int rowstart=0, rowend=0;
    
    Environment base("package:gtools");
    Function ddirichlet = base["ddirichlet"];
    
    for(int i = 0; i < nsites; i++)
    {
        
        lpchain1_current[i] = newoffset(i, swap1);
        pchain1_current[i] = exp(lpchain1_current[i]);
        
        lpchain2_current[i] = newoffset(i, swap2);
        pchain2_current[i] = exp(lpchain2_current[i]);
        
        lpchain1_proposal[i] = newoffset(i, swap2);
        pchain1_proposal[i] = exp(lpchain1_proposal[i]);
        
        lpchain2_proposal[i] = newoffset(i, swap1);
        pchain2_proposal[i] = exp(lpchain2_proposal[i]);
        
        chain1oldlikebit += y[i] * lpchain1_current[i] - pchain1_current[i];
        chain2oldlikebit += y[i] * lpchain2_current[i] - pchain2_current[i];
        
        chain1newlikebit += y[i] * lpchain1_proposal[i] - pchain1_proposal[i];
        chain2newlikebit += y[i] * lpchain2_proposal[i] - pchain2_proposal[i];
    }
    
    chain1likebit = chain1newlikebit - chain1oldlikebit;
    chain2likebit = chain2newlikebit - chain2oldlikebit;
    
    for(int j = 0; j < K; j++)
    {
        
        for(int l = 0; l < Ntrends; l++)
        {
            wpriorchain1old += wnew((beginchain1 + j), l) * log(lambdanew(l, swap1));
            wpriorchain1new += wnew((beginchain2 + j), l) * log(lambdanew(l, swap1));
            
            wpriorchain2old += wnew((beginchain2 + j), l) * log(lambdanew(l, swap2));
            wpriorchain2new += wnew((beginchain1 + j), l) * log(lambdanew(l, swap2));
        }
        
        priorvardenomchain1 = rhonew[swap1] * Wtripletsum[j] + 1 - rhonew[swap1];
        phipriorvarchain1 = tau2new[swap1] / priorvardenomchain1;
        
        priorvardenomchain2 = rhonew[swap2] * Wtripletsum[j] + 1 - rhonew[swap2];
        phipriorvarchain2 = tau2new[swap2] / priorvardenomchain2;
        
        rowstart = Wbegfin(j, 0) - 1;
        rowend = Wbegfin(j, 1);
        for(int l = rowstart; l < rowend; l++)
        {
            sumphichain1 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap1);
            sumphichain2 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap2);
        }
        phipriormeanchain1 = rhonew[swap1] * sumphichain1 / priorvardenomchain1; 
        phipriormeanchain2 = rhonew[swap2] * sumphichain2 / priorvardenomchain2; 
        
        phichain1 += (0.5/phipriorvarchain1) * pow((phinew(j, swap1) - phipriormeanchain1), 2) - (0.5/phipriorvarchain1) * pow((phinew(j, swap2) - phipriormeanchain1), 2);
        phichain2 += (0.5/phipriormeanchain2) * pow((phinew(j, swap2) - phipriormeanchain2), 2) - (0.5/phipriormeanchain2) * pow((phinew(j, swap1) - phipriormeanchain2), 2);
    }
    
    for(int k = 0; k < p; k++)     
    {
        chain1priorbitbeta += 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        chain2priorbitbeta += 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    }
    
    for(int l = 1; l < TrendSel; l++)
    {
        chain1priorbitgamma +=  0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l];
        chain2priorbitgamma +=  0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l];
    }
    
    chain1priorbitw = wpriorchain1new - wpriorchain1old;  
    chain2priorbitw = wpriorchain2new - wpriorchain2old;  
    
    chain1priorbitphi = phichain1;  
    chain2priorbitphi = phichain2; 
    
    chain1priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) ));
    chain2priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) ));
    
    chain1priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]) - (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]);
    chain2priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]) - (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]);
    
    // Compute the acceptance probability and return the value
    acceptance = exp(((chain1likebit + chain1priorbitbeta + chain1priorbitgamma + chain1priorbitw + chain1priorbitlambda + chain1priorbitphi + chain1priorbittau2) * temps[swap1]) + ((chain2likebit + chain2priorbitbeta + chain2priorbitgamma + chain2priorbitw + chain2priorbitlambda + chain2priorbitphi + chain2priorbittau2) * temps[swap2]));
    
    if(runif(1)[0] <= acceptance) 
    {
        accept = accept + 1;
    }
    else
    { 
    }
    return accept;
}


// [[Rcpp::export]]
List binomialdevfit(NumericVector y, NumericVector trials, NumericMatrix probs, const int nsites, const int Nchains)
{
    NumericVector newy = clone(y), newtrials = clone(trials);
    NumericMatrix newprobs = clone(probs), like_all(nsites, Nchains);
    NumericVector deviance(Nchains), deviance_all(nsites), binprob(nsites);
    
    Environment base2("package:stats");
    Function dbinom= base2["dbinom"];
    
    for(int j = 0; j < Nchains; j++)
    {
        binprob = newprobs(_, j);
        deviance_all = as<NumericVector>( dbinom(newy, newtrials, binprob) );
        like_all(_, j) = deviance_all; 
        deviance[j] = -2 * sum(log(deviance_all));
    }
    
    List out(2);
    out[0] = deviance;
    out[1] = like_all;
    return out;
    
    // return deviance;
}


// [[Rcpp::export]]
NumericVector binomialbetablockupdate(const int nsites, NumericMatrix beta, NumericMatrix betaprop, NumericMatrix lp_beta, NumericMatrix lp_betaprop,
                                      NumericMatrix offset, NumericVector y, NumericVector failures, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                                      const int Nchains, NumericVector temps, const int p)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double oldlikebit=0, newlikebit=0, likebit, priorbit=0;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites), acceptance(Nchains);
    NumericMatrix betanew = clone(beta), betapropnew = clone(betaprop);
    NumericMatrix newlp_beta = clone(lp_beta), newlp_betaprop = clone(lp_betaprop), newoffset = clone(offset);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        priorbit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newlp_beta(i, j) + newoffset(i, j);
            lp_proposal[i] = newlp_betaprop(i, j) + newoffset(i, j);
            
            p_current[i] = exp(lp_current[i]) / (1 + exp(lp_current[i]));
            p_proposal[i] = exp(lp_proposal[i]) / (1 + exp(lp_proposal[i]));
            
            // underflow
            if(p_current[i] == 1)
            {
                p_current[i] = 0.999;
            }else
            {}
            
            if(p_proposal[i] == 1)
            {
                p_proposal[i] = 0.999;
            }else
            {}
            
            if(p_current[i] == 0)
            {
                p_current[i] = 0.001;
            }else
            {}
            
            if(p_proposal[i] == 0)
            {
                p_proposal[i] = 0.001;
            }else
            {}
            //
            
            oldlikebit += y[i] * log(p_current[i]) + failures[i] * log((1 - p_current[i]));
            newlikebit += y[i] * log(p_proposal[i]) + failures[i] * log((1 - p_proposal[i]));
        }
        likebit = newlikebit - oldlikebit;
        
        // Create the prior acceptance component
        for(int k = 0; k < p; k++)
        {
            priorbit += 0.5 * pow((betanew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betapropnew(k, j) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        }
        
        // Compute the acceptance probability and return the value
        acceptance[j] = exp((likebit + priorbit) * temps[j]);
    }
    
    return acceptance;
}


// [[Rcpp::export]]
List binomialgammaupdate(const int nsites, NumericVector gamma, NumericVector proposal, NumericMatrix offset, NumericMatrix offset_proposal, NumericVector y,
                         NumericVector failures, double prior_meangamma, double prior_vargamma, const int Nchains, NumericVector temps)
{
    // Compute the acceptance probability for gamma
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(nsites), lp_proposal(nsites), p_current(nsites), p_proposal(nsites);
    NumericVector accept(Nchains);
    NumericVector gammanew = clone(gamma), proposalnew = clone(proposal);
    NumericMatrix newoffset = clone(offset), newproposal = clone(offset_proposal);
    
    for(int j = 0; j < Nchains; j++)
    {
        oldlikebit = 0;
        newlikebit = 0;
        
        for(int i = 0; i < nsites; i++)
        {
            lp_current[i] = newoffset(i, j);
            lp_proposal[i] = newproposal(i, j);
            
            p_current[i] = exp(lp_current[i]) / (1 + exp(lp_current[i]));
            p_proposal[i] = exp(lp_proposal[i]) / (1 + exp(lp_proposal[i]));
            
            oldlikebit += y[i] * log(p_current[i]) + failures[i] * log((1 - p_current[i]));
            newlikebit += y[i] * log(p_proposal[i]) + failures[i] * log((1 - p_proposal[i]));
        }
        likebit = newlikebit - oldlikebit;
        
        priorbit = 0.5 * pow((gammanew[j] - prior_meangamma), 2) / prior_vargamma - 0.5 * pow((proposalnew[j] - prior_meangamma), 2) / prior_vargamma;
        
        // Compute the acceptance probability and return the value
        acceptance = exp((likebit + priorbit) * temps[j]);
        if(runif(1)[0] <= acceptance)
        {
            gammanew[j] = proposalnew[j];
            accept[j] = accept[j] + 1;
        }
        else
        {
        }
    }
    
    List out(2);
    out[0] = gammanew;
    out[1] = accept;
    return out;
}

// [[Rcpp::export]]
List binomialwupdate(const int nsites, const int ntimes, NumericMatrix w, NumericMatrix offset, NumericMatrix offset_proposal, NumericMatrix w_proposal,
                     NumericMatrix y, NumericMatrix failures, NumericMatrix lambda, const int Nchains, NumericVector temps, NumericVector begin, NumericVector regbegin, const int Ntrends)
{
    // Compute the acceptance probability for beta
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix wnew = clone(w), wproposal = clone(w_proposal), lambdanew = clone(lambda), newoffset = clone(offset), newproposal = clone(offset_proposal);
    int rowstart, regstart;
    NumericVector proposal;
    
    for(int j = 0; j < Nchains; j++)
    {
        rowstart = begin[j] - 1;
        
        for(int k = 0; k < nsites; k++)
        {
            proposal = wproposal((rowstart + k), _);
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j);
                lp_proposal[t] = newproposal((regstart + k), j);
                
                p_current[t] = exp(lp_current[t]) / (1 + exp(lp_current[t]));
                p_proposal[t] = exp(lp_proposal[t]) / (1 + exp(lp_proposal[t]));
                
                oldlikebit += y(k, t) * log(p_current[t]) + failures(k, t) * log((1 - p_current[t]));
                newlikebit += y(k, t) * log(p_proposal[t]) + failures(k, t) * log((1 - p_proposal[t]));
            }
            
            likebit = newlikebit - oldlikebit;
            
            for(int i = 0; i < Ntrends; i++)
            {
                priorbit += proposal[i] * log(lambdanew(i, j)) - wnew((rowstart + k), i) * log(lambdanew(i, j));
            }
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                wnew((rowstart + k), _) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = wnew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
List binomialphiupdate(NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector Wtripletsum, const int nsites, const int ntimes,
                       NumericMatrix phi, NumericMatrix offset, NumericMatrix y, NumericMatrix failures, NumericVector tau2, NumericVector rho, const int Nchains,
                       NumericVector temps, NumericMatrix phi_tune, NumericVector regbegin)
{
    // Compute the acceptance probability
    //Create new objects
    double acceptance, oldlikebit=0, newlikebit=0, likebit, priorbit;
    NumericVector lp_current(ntimes), lp_proposal(ntimes), p_current(ntimes), p_proposal(ntimes);
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix accept(nsites, Nchains);
    NumericMatrix phinew = clone(phi), newoffset = clone(offset);
    double proposal;
    double priorvardenom, priormean, priorvar, sumphi;
    int rowstart=0, rowend=0;
    int regstart;
    
    for(int j = 0; j < Nchains; j++)
    {
        
        for(int k = 0; k < nsites; k++)
        {
            
            oldlikebit = 0;
            newlikebit = 0;
            priorbit = 0;
            
            // Calculate prior variance
            priorvardenom = rhonew[j] * Wtripletsum[k] + 1 - rhonew[j];
            priorvar = tau2new[j] / priorvardenom;
            
            // Calculate the prior mean
            rowstart = Wbegfin(k, 0) - 1;
            rowend = Wbegfin(k, 1);
            sumphi = 0;
            for(int l = rowstart; l < rowend; l++) sumphi += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), j);
            priormean = rhonew[j] * sumphi / priorvardenom; 
            
            // propose a value 
            proposal = rnorm(1, phinew(k, j), sqrt(priorvar*phi_tune(k, j)))[0];
            
            for(int t = 0; t < ntimes; t++)
            {
                regstart = regbegin[t] - 1;
                
                lp_current[t] = newoffset((regstart + k), j) + phinew(k, j);
                lp_proposal[t] = newoffset((regstart + k), j) + proposal;
                
                p_current[t] = exp(lp_current[t]) / (1 + exp(lp_current[t]));
                p_proposal[t] = exp(lp_proposal[t]) / (1 + exp(lp_proposal[t]));
                
                oldlikebit += y(k, t) * log(p_current[t]) + failures(k, t) * log((1 - p_current[t]));
                newlikebit += y(k, t) * log(p_proposal[t]) + failures(k, t) * log((1 - p_proposal[t]));
            }
            
            likebit = newlikebit - oldlikebit;
            priorbit = (0.5/priorvar) * pow((phinew(k, j) - priormean), 2) - (0.5/priorvar) * pow((proposal - priormean), 2);
            
            // Compute the acceptance probability and return the value
            acceptance = exp((likebit + priorbit) * temps[j]);
            if(runif(1)[0] <= acceptance) 
            {
                phinew(k, j) = proposal;
                accept(k, j) = accept(k, j) + 1;
            }
            else
            { 
            }
        }
    }
    List out(2);
    out[0] = phinew;
    out[1] = accept;
    return out;
}


// [[Rcpp::export]]
int binomialcouplingAllupdate(const int nsites, const int K, const int p, NumericMatrix w, NumericMatrix offset, NumericMatrix beta, NumericMatrix gamma,
                              NumericMatrix lambda, NumericMatrix phi, NumericVector rho, NumericVector tau2, NumericVector Wtripletsum,
                              NumericMatrix Wtriplet, NumericMatrix Wbegfin, NumericVector y, NumericVector failures, NumericVector prior_meanbeta, NumericVector prior_varbeta,
                              NumericVector prior_meantrends, NumericVector prior_vartrends, NumericVector prior_lambda, NumericVector prior_tau2,
                              NumericVector swap, NumericVector temps, NumericVector begin, const int Ntrends, const int TrendSel)
{
    // Compute the acceptance probability for swapping
    //Create new objects
    double acceptance, chain1oldlikebit=0, chain2oldlikebit=0, chain1newlikebit=0, chain2newlikebit=0, chain1likebit, chain2likebit, chain1priorbitbeta=0, chain2priorbitbeta=0;
    double wpriorchain1old=0, wpriorchain1new=0, wpriorchain2old=0, wpriorchain2new=0, chain1priorbitw, chain2priorbitw, chain1priorbitlambda, chain2priorbitlambda;
    double phichain1=0, phichain2=0;
    double chain1priorbitgamma=0, chain2priorbitgamma=0;
    double chain1priorbittau2, chain2priorbittau2, chain1priorbitphi, chain2priorbitphi;
    NumericVector lpchain1_current(nsites), lpchain1_proposal(nsites), pchain1_current(nsites), pchain1_proposal(nsites);
    NumericVector lpchain2_current(nsites), lpchain2_proposal(nsites), pchain2_current(nsites), pchain2_proposal(nsites);
    int accept = 0;
    NumericVector swapnew;
    swapnew = swap - 1;
    double swap1, swap2;
    swap1 = swapnew[0];
    swap2 = swapnew[1];
    NumericVector rhonew = clone(rho), tau2new = clone(tau2);
    NumericMatrix gammanew = clone(gamma), wnew = clone(w), newoffset = clone(offset), betanew = clone(beta), lambdanew = clone(lambda), phinew = clone(phi);
    int beginchain1 = begin[swap1] - 1;
    int beginchain2 = begin[swap2] - 1;
    double temp1, temp2;
    temp1 = temps[swap1];
    temp2 = temps[swap2];
    double priorvardenomchain1, priorvardenomchain2, phipriorvarchain1, phipriorvarchain2, sumphichain1=0, sumphichain2=0, phipriormeanchain1, phipriormeanchain2;
    int rowstart=0, rowend=0;
    
    Environment base("package:gtools");
    Function ddirichlet = base["ddirichlet"];
    
    for(int i = 0; i < nsites; i++)
    {
        lpchain1_current[i] = newoffset(i, swap1);
        pchain1_current[i] = exp(lpchain1_current[i]) / (1 + exp(lpchain1_current[i]));
        
        lpchain2_current[i] = newoffset(i, swap2);
        pchain2_current[i] = exp(lpchain2_current[i]) / (1 + exp(lpchain2_current[i]));
        
        lpchain1_proposal[i] = newoffset(i, swap2);
        pchain1_proposal[i] = exp(lpchain1_proposal[i]) / (1 + exp(lpchain1_proposal[i]));
        
        lpchain2_proposal[i] = newoffset(i, swap1);
        pchain2_proposal[i] = exp(lpchain2_proposal[i]) / (1 + exp(lpchain2_proposal[i]));
        
        chain1oldlikebit += y[i] * log(pchain1_current[i]) + failures[i] * log((1 - pchain1_current[i]));
        chain2oldlikebit += y[i] * log(pchain2_current[i]) + failures[i] * log((1 - pchain2_current[i]));
        
        chain1newlikebit += y[i] * log(pchain1_proposal[i]) + failures[i] * log((1 - pchain1_proposal[i]));
        chain2newlikebit += y[i] * log(pchain2_proposal[i]) + failures[i] * log((1 - pchain2_proposal[i]));
    }
    
    chain1likebit = chain1newlikebit - chain1oldlikebit;
    chain2likebit = chain2newlikebit - chain2oldlikebit;
    
    for(int j = 0; j < K; j++)
    {
        
        for(int l = 0; l < Ntrends; l++)
        {
            wpriorchain1old += wnew((beginchain1 + j), l) * log(lambdanew(l, swap1));
            wpriorchain1new += wnew((beginchain2 + j), l) * log(lambdanew(l, swap1));
            
            wpriorchain2old += wnew((beginchain2 + j), l) * log(lambdanew(l, swap2));
            wpriorchain2new += wnew((beginchain1 + j), l) * log(lambdanew(l, swap2));
        }
        
        priorvardenomchain1 = rhonew[swap1] * Wtripletsum[j] + 1 - rhonew[swap1];
        phipriorvarchain1 = tau2new[swap1] / priorvardenomchain1;
        
        priorvardenomchain2 = rhonew[swap2] * Wtripletsum[j] + 1 - rhonew[swap2];
        phipriorvarchain2 = tau2new[swap2] / priorvardenomchain2;
        
        rowstart = Wbegfin(j, 0) - 1;
        rowend = Wbegfin(j, 1);
        for(int l = rowstart; l < rowend; l++)
        {
            sumphichain1 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap1);
            sumphichain2 += Wtriplet(l, 2) * phinew((Wtriplet(l, 1) - 1), swap2);
        }
        phipriormeanchain1 = rhonew[swap1] * sumphichain1 / priorvardenomchain1; 
        phipriormeanchain2 = rhonew[swap2] * sumphichain2 / priorvardenomchain2; 
        
        phichain1 += (0.5/phipriorvarchain1) * pow((phinew(j, swap1) - phipriormeanchain1), 2) - (0.5/phipriorvarchain1) * pow((phinew(j, swap2) - phipriormeanchain1), 2);
        phichain2 += (0.5/phipriormeanchain2) * pow((phinew(j, swap2) - phipriormeanchain2), 2) - (0.5/phipriormeanchain2) * pow((phinew(j, swap1) - phipriormeanchain2), 2);
    }
    
    for(int k = 0; k < p; k++)     
    {
        chain1priorbitbeta += 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k];
        chain2priorbitbeta += 0.5 * pow((betanew(k, swap2) - prior_meanbeta[k]), 2) / prior_varbeta[k] - 0.5 * pow((betanew(k, swap1) - prior_meanbeta[k]), 2) / prior_varbeta[k];
    }
    
    for(int l = 1; l < TrendSel; l++)
    {
        chain1priorbitgamma +=  0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l];
        chain2priorbitgamma +=  0.5 * pow((gammanew(l, swap2) - prior_meantrends[l]), 2) / prior_vartrends[l] - 0.5 * pow((gammanew(l, swap1) - prior_meantrends[l]), 2) / prior_vartrends[l];
    }
    
    chain1priorbitw = wpriorchain1new - wpriorchain1old;  
    chain2priorbitw = wpriorchain2new - wpriorchain2old;  
    
    chain1priorbitphi = phichain1;  
    chain2priorbitphi = phichain2; 
    
    chain1priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) ));
    chain2priorbitlambda = log(as<double>( ddirichlet(lambdanew(_, swap2), prior_lambda) )) - log(as<double>( ddirichlet(lambdanew(_, swap1), prior_lambda) ));
    
    chain1priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]) - (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]);
    chain2priorbittau2 = (-prior_tau2[0] - 1) * log(tau2[swap2]) - (prior_tau2[1] / tau2[swap2]) - (-prior_tau2[0] - 1) * log(tau2[swap1]) - (prior_tau2[1] / tau2[swap1]);
    
    // Compute the acceptance probability and return the value
    acceptance = exp(((chain1likebit + chain1priorbitbeta + chain1priorbitgamma + chain1priorbitw + chain1priorbitlambda + chain1priorbitphi + chain1priorbittau2) * temps[swap1]) + ((chain2likebit + chain2priorbitbeta + chain2priorbitgamma + chain2priorbitw + chain2priorbitlambda + chain2priorbitphi + chain2priorbittau2) * temps[swap2]));
    
    if(runif(1)[0] <= acceptance) 
    {
        accept = accept + 1;
    }
    else
    { 
    }
    return accept;
}











////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Graph based optimisation algorithm  /////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
// This file contains the following functions:
// std::vector<std::vector<int>> optimise_graph(std::vector<std::vector<int>> adj, std::vector<double> data,
//    bool add=false, bool remove=true, bool remove_first=false) {
// Note that this function appears right at the end of this file.

// Needed for set_difference
#include <algorithm>
// shared_ptr
#include <memory>

// Create anonymous namespace to hold graph classes. Anything in an anonymous
// namespace is only available in this file.
namespace {


// Iterators to nicely generate powersets with minimal data copying
template<typename T>
class PowersetIterator {
private:
  std::vector<bool> counter_;
  const std::vector<T> & elts;
  bool end_;
public:
  typedef std::shared_ptr<const std::set<T>> value_type;
  typedef void difference_type;
  typedef value_type* pointer;
  typedef value_type& reference;
  typedef std::input_iterator_tag iterator_category;
  
  explicit PowersetIterator(const std::vector<T> & orig, bool end=false) : counter_(orig.size(), false), elts(orig), end_(end) { }
  std::shared_ptr<const std::set<T>> operator*() const {
    std::set<T> set;
    for(size_t i = 0; i < elts.size(); ++i) {
      if (counter_[i] == true) {
        set.insert(elts[i]);
      }
    }
    return std::make_shared<const std::set<T>>(set);
  }
  
  bool operator==(const PowersetIterator& other) const {
    // counter_ will be all false for both begin() and end(), so we use the
    // end_ member to note that we actually are at the end
    if ((end_ == true) && (other.end_ == true)) {
      return true;
    }
    if ((end_ == true) && (other.end_ == false)) {
      return false;
    }
    if ((end_ == false) && (other.end_ == true)) {
      return false;
    }
    return counter_ == other.counter_;
  }
  
  bool operator!=(const PowersetIterator& other) const { return !(*this == other); }
  
  PowersetIterator& operator++() {
    for(size_t i = 0; i < counter_.size(); ++i) {
      counter_[i] = !counter_[i];
      if (counter_[i] == true) {
        break;
      }
      // If we've flipped all bits, and not yet flipped anything
      // false->true, we must have just been at the point where everything
      // is true aka the last subset.
      if (i == (counter_.size() - 1)) {
        end_ = true;
      }
    }
    // If we are finding the powerset of an empty set, we find one empty
    // set and then reach the end.
    if (counter_.size() == 0) {
      end_ = true;
    }
    return *this;
  }
};

template<typename T>
class Powerset {
private:
  std::vector<T> elts;
public:
  Powerset(std::set<T> elts_) : elts(elts_.begin(), elts_.end()) { }
  PowersetIterator<T> begin() { return PowersetIterator<T>(elts); }
  PowersetIterator<T> end() { return PowersetIterator<T>(elts, true); }
};



// An edge in a graph
typedef std::tuple<int,int> Edge;

// A graph, which internally stores the adjacency matrix, and the current
// neighbourhood of each vertex.
class Graph {
public:
  Graph(std::vector<std::vector<int>> adj) {
    this->adj_ = adj;
    this->size_ = adj.size();
    this->nbs_.resize(this->size_);
    int row_id = 0;
    for(const auto row: this->adj_) {
      int col_id = 0;
      for(const auto entry: row) {
        if (entry == 1) {
          this->nbs_[row_id].insert(col_id);
        }
        col_id += 1;
      }
      row_id += 1;
    }
  }
  
  Graph(int size) : adj_(size, std::vector<int>(size)) {
    this->size_ = size;
    this->nbs_.resize(this->size_);
  }
  
  std::vector<std::vector<int>> adjMatrix() const {
    return std::vector<std::vector<int>>(this->adj_);
  }
  
  int size() const {
    return this->size_;
  }
  
  bool has_edge(int u, int v) {
    return this->adj_[u][v] == 1;
  }
  
  void add_edge(const Edge& edge) {
    int v0, v1;
    std::tie(v0, v1) = edge;
    this->adj_[v0][v1] = 1;
    this->adj_[v1][v0] = 1;
    this->nbs_[v0].insert(v1);
    this->nbs_[v1].insert(v0);
  }
  
  void add_edges(const std::list<Edge> & edges) {
    for(const auto edge: edges) {
      this->add_edge(edge);
    }
  }
  
  void remove_edge(const Edge& edge) {
    int v0, v1;
    std::tie(v0, v1) = edge;
    this->adj_[v0][v1] = 0;
    this->adj_[v1][v0] = 0;
    this->nbs_[v0].erase(v1);
    this->nbs_[v1].erase(v0);
  }
  void remove_edges(std::list<Edge> to_remove) {
    for(auto edge: to_remove) {
      this->remove_edge(edge);
    }
  }
  
  int degree(int v) const {
    return this->nbs_[v].size();
  }
  
  const std::set<int> nbs(int v) const {
    return this->nbs_[v];
  }
  
  const std::list<Edge> edges() const {
    std::list<Edge> result;
    int row_id = 0;
    for(const auto row: this->adj_) {
      int col_id = 0;
      for(const auto entry: row) {
        if (entry == 1) {
          result.push_back(std::make_tuple(col_id, row_id));
        }
        // Only check and return where col_id <= row_id
        if (col_id == row_id) {
          break;
        }
        col_id += 1;
      }
      row_id += 1;
    }
    return result;
  }
  
private:
  size_t size_;
  std::vector<std::vector<int>> adj_;
  std::vector<std::set<int>> nbs_;
};

// An "EmptyGraph" has a set of potential edges. This way we can restrict
// what edges we might add.
class EmptyGraph : public Graph {
  
public:
  EmptyGraph(std::vector<std::vector<int>> adj) : Graph(adj.size()), _possible(adj) {
  }
  
  const std::set<int> possNbs(int v) const {
    return this->_possible.nbs(v);
  }
  
  const std::list<Edge> origEdges() const {
    return this->_possible.edges();
  }
  
private:
  Graph _possible;
};

// A class to do our optimising. Using a class means we can store one central
// copy of the graph and the data, and still access them quickly
class Optimiser {
public:
  Optimiser(std::vector<std::vector<int>> adj, std::vector<double> data) : _graph(adj), _data(data) { }
  
  std::vector<std::vector<int>> adjMatrix() const { return this->_graph.adjMatrix(); }
  
  void iterative_opt(bool remove, bool add, bool remove_first) {
    if ((!add) || (remove && remove_first)) {
      this->_graph.add_edges(this->_graph.origEdges());
    } else {
      // Building from a mostly empty graph
      for(int v = 0; v < this->_graph.size(); ++v) {
        if (this->_graph.degree(v) >= 1) {
          continue;
        }
        signed int bestNb = -1;
        double minDiff = 0;
        for (auto u: this->_graph.possNbs(v)) {
          double diff = fabs(this->_data[v] - this->_data[u]);
          if ((bestNb == -1) || (diff < minDiff)) {
            bestNb = u;
            minDiff = diff;
          }
        }
        //Rprintf("Setup: Adding edge (%u, %u)\n", v, bestNb);
        this->_graph.add_edge(std::make_pair(v, bestNb));
      }
    }
    bool changed = true;
    while (changed) {
      changed = false;
      if (add && (! (remove_first && remove))) {
        double avsum = this->avsum();
        double oldscore = this->score();
        auto lastAdded = this->greedy_opt_add(avsum);
        double newscore = this->score();
        //Rprintf("adding\toldscore = %f\tnewscore = %f\tchange size is %u\n", oldscore, newscore, lastAdded.size());
        if (!lastAdded.empty()) {
          if (newscore > oldscore) {
            // Good change
            changed = true;
          } else {
            // Bad change, put the edges back
            this->_graph.remove_edges(lastAdded);
          }
        }
      }
      if (remove) {
        // Flip flag off so we add in next iteration.
        remove_first = false;
        double avsum = this->avsum();
        double oldscore = this->score();
        auto lastRemoved = this->greedy_opt_remove(avsum);
        double newscore = this->score();
        //Rprintf("removing\toldscore = %f\tnewscore = %f\tchange size is %u\n", oldscore, newscore, lastRemoved.size());
        if (!lastRemoved.empty()) {
          if (newscore > oldscore) {
            // Good change
            changed = true;
          } else {
            // Bad change, put the edges back
            this->_graph.add_edges(lastRemoved);
          }
        }
      }
    }
  }
private:
  EmptyGraph _graph;
  std::vector<double> _data;
  
  double disc(int v) const {
    double sum = 0;
    for(auto nb: this->_graph.nbs(v)) {
      sum += this->_data[nb];
    }
    return this->_graph.degree(v)*std::pow(this->_data[v] - sum/this->_graph.degree(v), 2);
  }
  
  double score() const {
    double first = 0;
    double second_inner = 0;
    for(int v = 0; v < this->_graph.size(); ++v) {
      first += std::log(this->_graph.degree(v));
      second_inner += this->disc(v);
    }
    return first/2.0f - (this->_graph.size()/2.0f)*std::log(second_inner);
  }
  
  double avsum() const {
    double sum = 0;
    for(int v = 0; v < this->_graph.size(); ++v) {
      sum += this->disc(v);
    }
    return sum;
  }
  
  double vertex_val(int v, const std::set<int> & fixed_nbs, const std::set<int> & opt_nbs, double avsum) const {
    double sum = 0;
    for(auto nb: fixed_nbs) {
      sum += this->_data[nb];
    }
    for(auto nb: opt_nbs) {
      sum += this->_data[nb];
    }
    int degree = fixed_nbs.size() + opt_nbs.size();
    double vx_av = degree * std::pow(this->_data[v] - sum/degree, 2);
    if (vx_av > avsum) {
      // The discrepancy at this vertex is already worse than the total
      // discrepancy in the old graph. This means this value should be bad,
      // but also breaks the calculations so we just return something very
      // negative.
      return -1e20;
    }
    return std::log(degree)/2.0f - (this->_graph.size()/2.0f)*std::log1p(vx_av / (avsum - vx_av));
  }
  
  // BEWARE: I define a < b to be true if a[1] > b[1] to get a reversed
  // multiset.
  class CondSortedEntry : public std::tuple<std::shared_ptr<const std::set<int>>, double> {
  public:
    using BaseType = std::tuple<std::shared_ptr<const std::set<int>>, double>;
    CondSortedEntry(const BaseType& other) : BaseType(other) { }
    bool operator<(const CondSortedEntry other) const {
      return std::get<1>(*this) > std::get<1>(other);
    }
  };
  
  typedef std::multiset<CondSortedEntry> CondSortedList;
  
  CondSortedList cond_sorted_nbs_remove(int v, const std::set<int> fixed_nbs, const std::set<int> options, double avsum) const {
    CondSortedList result;
    Powerset<int> pow(options);
    for(const auto option: pow) {
      // Must have degree >= 1
      if (fixed_nbs.size() + option->size() == 0) {
        continue;
      }
      result.insert(std::make_tuple(option, this->vertex_val(v, fixed_nbs, *option, avsum)));
    }
    return result;
  }
  
  CondSortedList cond_sorted_nbs_add(int v, const std::set<int> fixed, const std::set<int> options, double avsum) const {
    CondSortedList result;
    Powerset<int> pow(options);
    for(const auto option: pow) {
      // Must have degree >= 1
      if (fixed.size() + option->size() == 0) {
        continue;
      }
      result.insert(std::make_tuple(option, this->vertex_val(v, fixed, *option, avsum)));
    }
    return result;
  }
  
  double find_best_cond(const CondSortedList & list, signed int mustHave, signed int cantHave) const {
    for(const CondSortedEntry entry: list) {
      const std::set<int> & opts = * std::get<0>(entry);
      if ((mustHave != -1) && (opts.find((int)mustHave) == opts.end())) {
        continue;
      }
      if ((cantHave != -1) && (opts.find((int)cantHave) != opts.end())) {
        continue;
      }
      return std::get<1>(entry);
    }
    // TODO Error if we reach here?
    Rprintf("ERROR ---------------- ERROR\n");
    return 1;
  }
  
  std::list<Edge> greedy_opt_add(double avsum) {
    std::list<Edge> added;
    for(int v = 0; v < this->_graph.size(); v++) {
      std::set<int> nb_options;
      std::set<int> nb_fixed;
      for(int u: this->_graph.possNbs(v)) {
        // If we already have this edge, it is fixed for now.
        if (this->_graph.has_edge(u,v)) {
          nb_fixed.insert(u);
          continue;
        }
        // Only consider edges {v,u} where v < u, so this edge has not yet
        // been added, but has been considered.
        if (u < v) {
          continue;
        }
        nb_options.insert(u);
      }
      CondSortedList v_list = this->cond_sorted_nbs_add(v, nb_fixed, nb_options, avsum);
      std::list<Edge> add_here;
      for (auto u: nb_options) {
        // -1 indicates that we don't care.
        double vx_gain = this->find_best_cond(v_list, u, -1) - this->find_best_cond(v_list, -1, u);
        std::set<int> u_options;
        std::set<int> u_fixed;
        for(int u_nb: this->_graph.possNbs(u)) {
          if (this->_graph.has_edge(u, u_nb)) {
            u_fixed.insert(u_nb);
            continue;
          }
          if (u_nb < v) {
            continue;
          }
          u_options.insert(u_nb);
        }
        CondSortedList nbr_list = this->cond_sorted_nbs_add(u, u_fixed, u_options, avsum);
        double nb_gain = this->find_best_cond(nbr_list, v, -1) - this->find_best_cond(nbr_list, -1, v);
        //Rprintf("Adding (%u, %u), vx_gain is %f and nb_gain is %f\n", v, u, vx_gain, nb_gain);
        if (vx_gain + nb_gain > 0 ) {
          add_here.push_back(std::make_pair(v, u));
        }
      }
      this->_graph.add_edges(add_here);
      for(auto edge: add_here) {
        added.push_back(edge);
      }
    }
    return added;
  }
  
  std::list<Edge> greedy_opt_remove(double avsum) {
    std::list<Edge> removed;
    for(int v = 0; v < this->_graph.size(); v++) {
      if (this->_graph.degree(v) == 1) {
        continue;
      }
      std::set<int> nb_options;
      std::set<int> nb_fixed;
      for(int u: this->_graph.nbs(v)) {
        // Only consider edges {v,u} where v < u
        if ((u < v) || (this->_graph.degree(u) == 1)) {
          nb_fixed.insert(u);
        } else {
          nb_options.insert(u);
        }
      }
      CondSortedList v_list = this->cond_sorted_nbs_remove(v, nb_fixed, nb_options, avsum);
      std::list<Edge> remove_here;
      for (auto u: nb_options) {
        // Don't remove the last edge, just break
        if ((this->_graph.degree(v) - remove_here.size()) == 1) {
          break;
        }
        // -1 indicates that we don't care.
        double vx_gain = this->find_best_cond(v_list, u, -1) - this->find_best_cond(v_list, -1, u);
        std::set<int> u_fixed;
        std::set<int> u_options;
        for(int u_nb: this->_graph.nbs(u)) {
          if ((u_nb < v) || (this->_graph.degree(u_nb) == 1)) {
            u_fixed.insert(u_nb);
          } else {
            u_options.insert(u_nb);
          }
        }
        CondSortedList nbr_list = this->cond_sorted_nbs_remove(u, u_fixed, u_options, avsum);
        double nb_gain = this->find_best_cond(nbr_list, v, -1) - this->find_best_cond(nbr_list, -1, v);
        //Rprintf("Dropping (%u, %u), vx_gain is %f and nb_gain is %f\n", v, u, vx_gain, nb_gain);
        // the gains are from the edge being _in_ the graph, so we remove
        // if the gain is negative!
        if (vx_gain + nb_gain < 0 ) {
          remove_here.push_back(std::make_pair(v, u));
        }
      }
      this->_graph.remove_edges(remove_here);
      for(auto edge: remove_here) {
        removed.push_back(edge);
      }
    }
    return removed;
  }
};
} // namespace

// [[Rcpp::export]]
std::vector<std::vector<int>> optimise_graph(const IntegerMatrix& adj, const NumericVector& data,
                                             bool add=false, bool remove=true, bool remove_first=false) {
  //std::vector<double> dd = as<std::vector<double>>(data);
  std::vector<std::vector<int>> adj_;
  for(int r = 0; r < adj.nrow(); ++r) {
    std::vector<int> row;
    for(int c = 0; c < adj.ncol(); ++c) {
      row.push_back(adj(r,c));
    }
    adj_.push_back(row);
  }
  //Optimiser opt(as<std::vector<std::vector<int>>>(adj), dd);
  //Rcout << "adj is " << adj << "\n";
  //Rcout << "data is " << data << "\n";
  Optimiser opt(adj_, as<std::vector<double>>(data));
  opt.iterative_opt(remove, add, remove_first);
  return opt.adjMatrix();
}



