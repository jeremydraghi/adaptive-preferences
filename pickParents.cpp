#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix pickParents(int N, NumericMatrix pop, NumericMatrix fits, double p, double c, double r, NumericVector x) {
  NumericVector fitness(N);
  NumericMatrix retMat(N, 3);
  NumericVector intervals(N);
  IntegerVector parents(2);
  double p1;
  double p2;
  int tracker = 0;
  int cursor;
  double draw;
  int crossA;
  int crossB;
  int parentsNeeded;
  double tFit = 0;
  double mFit;
  
  for(int i = 0; i < N; i++)
  {
  	 p1 = (1 - p) / (1 - p * pop(i,0) * (1 - c));
  	 p2 = p * (1 - pop(i,0)) / (1 - p * pop(i,0) * (1 - c));
  	 if(x[tracker] > (p1 + p2))
  	 {
  	 	fitness[i] = 0;
  	 }
  	 else
  	 {
  	 	if(x[tracker] <= p1)
  	 	{
  	 		fitness[i] = fits((int)(pop(i,1))-1, 0);
  	 	}
  	 	else
  	 	{
  	 		fitness[i] = fits((int)(pop(i,2))-1, 1);
  	 	}
  	 }
  	 tracker++;
  }

  for(int i = 0; i < N; i++)
  {
  	 tFit += fitness[i];
  	 intervals[i] = tFit;
  }
  mFit = tFit / N;
  

  for(int i = 0; i < N; i++)
  {
  	if(x[tracker] < r)
  	{
  		crossA = 1;
  	}
  	else
  	{
  		crossA = 0;
  	}
  	tracker++;
  	if(x[tracker] < r)
  	{
  		crossB = 1;
  	}
  	else
  	{
  		crossB = 0;
  	}
  	tracker++;
  	parentsNeeded = 1;
  	if(crossA == 1 || crossB == 1) parentsNeeded++;
  	
  	for(int j = 0; j < parentsNeeded; j++)
  	{
  		draw = x[tracker] * tFit;
  		cursor = int(draw / mFit);
  		while(draw >= intervals[cursor]) cursor++;
  		while(cursor > 0 && draw < intervals[cursor-1]) cursor--;
  		parents[j] = cursor;
  		tracker++;
  	}
  		
  	retMat(i,0) = pop(parents[0],0);
  	if(crossA == 1)
  	{
  		retMat(i,1) = pop(parents[1],1);
  		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  	}
  	else
  	{
  		retMat(i,1) = pop(parents[0],1);
		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  	}
  } 
  return retMat;
}

// [[Rcpp::export]]
NumericMatrix pickParentsFreePref(int N, NumericMatrix pop, NumericMatrix fits, double p, double c, double r, NumericVector x) {
  NumericVector fitness(N);
  NumericMatrix retMat(N, 3);
  NumericVector intervals(N);
  IntegerVector parents(2);
  double p1;
  double p2;
  int tracker = 0;
  int cursor;
  double draw;
  int crossA;
  int crossB;
  int parentsNeeded;
  double tFit = 0;
  double mFit;
  
  for(int i = 0; i < N; i++)
  {
  	 p1 = (1 - p) / (1 - p * pop(i,0) * (1 - c));
  	 p2 = p * (1 - pop(i,0)) / (1 - p * pop(i,0) * (1 - c));
  	 if(x[tracker] > (p1 + p2))
  	 {
  	 	fitness[i] = 0;
  	 }
  	 else
  	 {
  	 	if(x[tracker] <= p1)
  	 	{
  	 		fitness[i] = fits((int)(pop(i,1))-1, 0);
  	 	}
  	 	else
  	 	{
  	 		fitness[i] = fits((int)(pop(i,2))-1, 1);
  	 	}
  	 }
  	 tracker++;
  }

  for(int i = 0; i < N; i++)
  {
  	 tFit += fitness[i];
  	 intervals[i] = tFit;
  }
  mFit = tFit / N;
  

  for(int i = 0; i < N; i++)
  {
  	if(x[tracker] < 0.5)
  	{
  		crossA = 1;
  	}
  	else
  	{
  		crossA = 0;
  	}
  	tracker++;
  	if(x[tracker] < r)
  	{
  		crossB = 1;
  	}
  	else
  	{
  		crossB = 0;
  	}
  	tracker++;
  	parentsNeeded = 1;
  	if(crossA == 1 || crossB == 1) parentsNeeded++;
  	
  	for(int j = 0; j < parentsNeeded; j++)
  	{
  		draw = x[tracker] * tFit;
  		cursor = int(draw / mFit);
  		while(draw >= intervals[cursor]) cursor++;
  		while(cursor > 0 && draw < intervals[cursor-1]) cursor--;
  		parents[j] = cursor;
  		tracker++;
  	}
  		
  	retMat(i,0) = pop(parents[0],0);
  	if(crossA == 1)
  	{
  		retMat(i,1) = pop(parents[1],1);
  		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  	}
  	else
  	{
  		retMat(i,1) = pop(parents[0],1);
		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  	}
  } 
  return retMat;
}


// [[Rcpp::export]]
NumericMatrix pickParentsCrowding(int N, NumericMatrix pop, NumericMatrix fits, double p, double c, double r, NumericVector x) {
  NumericVector fitness(N);
  NumericMatrix retMat(N, 3);
  NumericVector intervals(N);
  IntegerVector parents(2);
  double p1;
  double p2;
  int tracker = 0;
  int cursor;
  double draw;
  int crossA;
  int crossB;
  int parentsNeeded;
  double tFit = 0;
  double mFit;
  int countA = 0;
  IntegerVector env(N);
  
  for(int i = 0; i < N; i++)
  {
  	env[i] = 0;
  	 p1 = (1 - p) / (1 - p * pop(i,0) * (1 - c));
  	 p2 = p * (1 - pop(i,0)) / (1 - p * pop(i,0) * (1 - c));
  	 if(x[tracker] > (p1 + p2))
  	 {
  	 	fitness[i] = 0;
  	 }
  	 else
  	 {
  	 	if(x[tracker] <= p1)
  	 	{
  	 		fitness[i] = fits((int)(pop(i,1))-1, 0);
  	 		countA++;
  	 		env[i] = 1;
  	 	}
  	 	else
  	 	{
  	 		fitness[i] = fits((int)(pop(i,2))-1, 1);
  	 		env[i] = 2;
  	 	}
  	 }
  	 tracker++;
  }
  
  double crowdedA = (1 - p) * N / countA;
  double crowdedB = p * N / (N - countA);
  if(crowdedA > 1) crowdedA = 1;
  if(crowdedB > 1) crowdedB = 1;
  for(int i = 0; i < N; i++)
  {
  	 if(env[i] == 1) tFit += fitness[i] * crowdedA;
  	 if(env[i] == 2) tFit += fitness[i] * crowdedB;
  	 intervals[i] = tFit;
  }
  mFit = tFit / N;
  

  for(int i = 0; i < N; i++)
  {
  	if(x[tracker] < r)
  	{
  		crossA = 1;
  	}
  	else
  	{
  		crossA = 0;
  	}
  	tracker++;
  	if(x[tracker] < r)
  	{
  		crossB = 1;
  	}
  	else
  	{
  		crossB = 0;
  	}
  	tracker++;
  	parentsNeeded = 1;
  	if(crossA == 1 || crossB == 1) parentsNeeded++;
  	
  	for(int j = 0; j < parentsNeeded; j++)
  	{
  		draw = x[tracker] * tFit;
  		cursor = int(draw / mFit);
  		while(draw >= intervals[cursor]) cursor++;
  		while(cursor > 0 && draw < intervals[cursor-1]) cursor--;
  		parents[j] = cursor;
  		tracker++;
  	}
  		
  	retMat(i,0) = pop(parents[0],0);
  	if(crossA == 1)
  	{
  		retMat(i,1) = pop(parents[1],1);
  		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  	}
  	else
  	{
  		retMat(i,1) = pop(parents[0],1);
		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  	}
  } 
  return retMat;
}


// [[Rcpp::export]]
NumericMatrix pickParentsSoftSel(int N, NumericMatrix pop, NumericMatrix fits, double p, double c, double r, NumericVector x) {
  NumericVector fitness(N);
  NumericMatrix retMat(N, 3);
  NumericVector intervals(N);
  IntegerVector parents(2);
  double p1;
  double p2;
  int tracker = 0;
  int cursor;
  double draw;
  int crossA;
  int crossB;
  int parentsNeeded;
  double tFit = 0;
  double mFit;
  double meanFitA = 0;
  double meanFitB = 0;
  int countA = 0;
  IntegerVector env(N);
  
  for(int i = 0; i < N; i++)
  {
  	 env[i] = 0;
  	 p1 = (1 - p) / (1 - p * pop(i,0) * (1 - c));
  	 p2 = p * (1 - pop(i,0)) / (1 - p * pop(i,0) * (1 - c));
  	 if(x[tracker] > (p1 + p2))
  	 {
  	 	fitness[i] = 0;
  	 }
  	 else
  	 {
  	 	if(x[tracker] <= p1)
  	 	{
  	 		fitness[i] = fits((int)(pop(i,1))-1, 0);
  	 		meanFitA += fitness[i];
  	 		countA++;
  	 		env[i] = 1;
  	 	}
  	 	else
  	 	{
  	 		fitness[i] = fits((int)(pop(i,2))-1, 1);
  	 		meanFitB += fitness[i];
  	 		env[i] = 2;
  	 	}
  	 }
  	 tracker++;
  }
	
  meanFitA /= countA;
  meanFitB /= N - countA; 
  for(int i = 0; i < N; i++)
  {
  	 if(env[i] == 1) tFit += fitness[i] / meanFitA;
  	 if(env[i] == 2) tFit += fitness[i] / meanFitB;
  	 intervals[i] = tFit;
  }
  mFit = tFit / N;
  

  for(int i = 0; i < N; i++)
  {
  	if(x[tracker] < r)
  	{
  		crossA = 1;
  	}
  	else
  	{
  		crossA = 0;
  	}
  	tracker++;
  	if(x[tracker] < r)
  	{
  		crossB = 1;
  	}
  	else
  	{
  		crossB = 0;
  	}
  	tracker++;
  	parentsNeeded = 1;
  	if(crossA == 1 || crossB == 1) parentsNeeded++;
  	
  	for(int j = 0; j < parentsNeeded; j++)
  	{
  		draw = x[tracker] * tFit;
  		cursor = int(draw / mFit);
  		while(draw >= intervals[cursor]) cursor++;
  		while(cursor > 0 && draw < intervals[cursor-1]) cursor--;
  		parents[j] = cursor;
  		tracker++;
  	}
  		
  	retMat(i,0) = pop(parents[0],0);
  	if(crossA == 1)
  	{
  		retMat(i,1) = pop(parents[1],1);
  		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  	}
  	else
  	{
  		retMat(i,1) = pop(parents[0],1);
		if(crossB == 1)
  		{
  			retMat(i,2) = pop(parents[1],2);
  		}
  		else
  		{
  			retMat(i,2) = pop(parents[0],2);
  		}
  	}
  } 
  return retMat;
}