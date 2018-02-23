#' @title Optimal Exact Tests for Multiple Binary Endpoints
#'
#' @description Calculates global tests and multiple testing procedures to compare two groups with respect to multiple binary endpoints based on optimal rejection regions.
#'    
#' @param x a data frame of binary response vectors, or an array of numbers of failures and successes in the treatment group, or a list of marginal \emph{2 by 2} tables, see details.
#' @param y a vector of group allocations, or an array of numbers of failures and successes in the reference group, see details.
#' @param method a character variable indicating which optimization procedure to use.
#'	This can be one of \code{"alpha.greedy"}, \code{"alpha"}, \code{"number"}, \code{"power"} or \code{"bonferroni.greedy"}, see details.
#' @param alpha nominal significance level, the default is 0.025. Note that the test is one-sided.
#' @param p1 an array of assumed probabilities for failure and success in the treatment group, see details.
#' @param p0 an array of assumed probabilities for failure and success in the reference group, see details.
#' @param max.iter the maximal number of iterations in the branch and bound optimzation algorithm. Defaults to 10^5.
#' @param limit the value below which contributions to alpha are set to zero (and alpha is lowered accordingly) to speed up computation. Defaults to 0.
#' @param show.region logical, if \code{TRUE} a data frame indicating which possible outcome is element of the rejection region of the global test is added to the output. Defaults to \code{FALSE}.
#' @param closed.test logical, if \code{TRUE} adjusted p-values for the elementary null hypotheses are calculated by applying the specified test to all intersection hypotheses in a closed testing scheme. This can be 
#'	computer intensive, depending on the number of endpoints. 
#' @param consonant logical indicating if the global test should be constrained such that the resulting closed test is consonant. This option is only available for two endpoints. Note that the
#' Bonferroni greedy method is always consonant by construction.
#'
#' @details The null hypothesis for the global test is an identical multidimensional distribution of successes and failures in both groups.
#' The alternative hypothesis is a larger success proportion in the treatment group in at least one endpoint.
#'
#' \code{x} can be a data frame with one row per subject and one column for each endpoint. Only values of 0 or 1 are allowed,
#' with 0 indicating failure and 1 indicating success of the subject for the particular endpoint. In that case \code{y} needs to be a vector of group assignemnts with values 0 and 1,
#' where 0 is the reference group and 1 the treatment group.
#' Alternatively, \code{x} and \code{y} can be contingency tables in terms of \emph{2 by 2 by ... by 2} arrays. Each dimension of the array corresponds to one endpoint, the first coordinate position
#' in each dimension refers to failure in that endpoint, the second coordinate position refers to success. The array contains the number of subjects that were observed
#' for each possible combination of failures and successes. 
#' If \code{x} is a list of marginal \emph{2 by 2} tables, the Bonferroni greedy method is used. Matching the other input
#' variants, the \emph{2 by 2} tables are assumed to have the number of failures in the first row and the number of successes in the second row, and the first column to correspond to
#' the reference group, the second column to the treatment group.
#'
#' The methods \code{"alpha.greedy"}, \code{"alpha"}, \code{"number"} and \code{"power"} are based on the multivariate permutation distribution of the data conditional 
#' on the observed numbers of successes and failures across both groups. The method \code{"alpha.greedy"} uses a greedy algorithm aiming to exhaust the nominal significance level.
#' The methods \code{"alpha"}, \code{"number"} and \code{"power"} use a branch and bound algorithm to find rejection regions with, respectively,
#' maximal exhaustion of the nominal significance level, maximal number of elements or maximal power for the alternative given by \code{p1} and \code{p0}. 
#' The method \code{"bonferroni.greedy"} uses a greedy algorithm aiming to exhaust the nominal significance level of a weighted Bonferroni adjustment of multiple Fisher's exact tests.
#' See reference for further details.
#'
#' \code{p1} and \code{p0} are \emph{2 by 2 by ... by 2} arrays. Each dimension of the array corresponds to one endpoint, the first coordinate position
#' in each dimension refers to failure in that endpoint, the second coordinate position refers to success. 
#' The array contains the assumed true probabilities for each possible combination of failures and successes. 
#'
#' @return A list with class \code{multfisher} containing the following components:
#' \describe{
#'	\item{\code{call}}{the function call.}
#'	\item{\code{data}}{a data frame showing the aggregated input data. If \code{p1} and \code{p0} are provided they are included in vectorized form.}
#'	\item{\code{alpha}}{the value of \code{alpha}.}
#'	\item{\code{method}}{the chosen method as found by argument match to \code{method}.}
#'	\item{\code{statistic}}{the vector of test statistics, these are the marginal numbers of successes in the treatment group.}
#'	\item{\code{p.value}}{the p-value of the global test. See reference for details on the calculation.}
#'	\item{\code{conditional.properties}}{a list of the actual significance level, the number of elements and the power of the global test. The values are calculated from the permutation 
#'		distribution of the date and they are conditional on the observed total numbers of successes and failures. The power is calculated for the alternative defined through
#'		\code{p1} and \code{p0}. If \code{p1} and \code{p0} are not specified, the value for power is \code{NA}.}
#'	\item{\code{rej.region}}{Provided if \code{show.region} is \code{TRUE} and method is in \code{c("alpha","number","power","alpha.greedy")}. A data frame showing in the column rejection.region
#'		if a multidimensional test statistic, indicated by the previous columns, is element of the rejection region (value of 1) or not (value of 0) for the global level alpha test.
#'		The column alpha gives the probability of  observing the particular vector of test statistics under the null hypothesis and conditional on the observed total numbers of
#'		successes and failures. Values of 0 occur if a combination of test statistics is not possible in the conditional distribution. The column power shows the conditional probability
#'		under the alternative defined through \code{p1} and \code{p0}. If \code{p1} and \code{p0} are not specified, the values for power are \code{NA}.}
#'	\item{\code{elementary.tests}}{a data frame showing for each endpoint the marginal odds ratio, the unadjusted one-sided p-value of Fisher's exact test and the adjusted 
#'		p-value resulting from application of the optimal exact test in a closed testing procedure.}
#'	\item{\code{closed.test}}{a data frame indicating all intersection hypotheses in the closed test and giving their p-values.}
#'	\item{\code{consonant.constraint}}{logical indicating whether the consonance constraint was used.}
#'	\item{\code{OPT}}{a list summarizing the optimization success, if applicable. The number of iterations of the branch and bound algorithm is given, as well as the 
#'		specified maximal iteration number and a logical variable indicating whether the optimization (in all steps of the closed test, if applicable) was finished.
#'		The number of iterations may be 0, which indicates that the optimization problem was solved in a pre-processing step.}
#' }
#'
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @references Robin Ristl, Dong Xi, Ekkehard Glimm, Martin Posch (2018), Optimal exact tests for multiple binary endpoints. 
#'     \emph{Computational Statistics and Data Analysis}, \strong{122}, 1-17. doi: 10.1016/j.csda.2018.01.001 (open access)
#'
#' @seealso \code{\link{print.multfisher}}, \code{\link{plot.multfisher}}
#'
#' @examples
#' ## Examples with two endpoints
#' data<-data.frame(endpoint1=c(0,0,1,1,1,0,0,0,0,1,1,1,1,1,1, 0,0,1,0,0,1,1,1,1,1,1,1,1,1,1),
#'			  endpoint2=c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1, 0,0,0,1,1,1,1,1,1,1,1,1,1,1,1),
#'			  group=rep(c(0,1),each=15))
#' ## maximal power under a specified alternative
#' p1<-matrix(c(0.1,0.2,0.2,0.5),2,2)
#' p0<-matrix(c(0.75,0.1,0.1,0.05),2,2)
#' rownames(p1)<-rownames(p0)<-c("EP1_failure","EP1_success")
#' colnames(p1)<-colnames(p0)<-c("EP2_failure","EP2_success")
#' testpower<-mfisher.test(x=data[,c(1:2)],y=data$group,method="power",
#'		p1=p1,p0=p0,closed.test=TRUE,show.region=TRUE)
#' print(testpower)
#' plot(testpower,cex=2)
#' str(testpower)
#'
#' ## maximal alpha with consonance constraint and using aggregated data as input
#' tab1<-table(data$endpoint1[data$group==1],data$endpoint2[data$group==1])
#' tab0<-table(data$endpoint1[data$group==0],data$endpoint2[data$group==0])
#' testalpha<-mfisher.test(x=tab1,y=tab0,method="alpha",closed.test=TRUE,
#'		show.region=TRUE,consonant=TRUE)
#' print(testalpha)
#' plot(testalpha,cex=2)
#'
#' ## Examples with three endpoints
#' data3EP<-data.frame(endpoint1=c(0,0,0,0,0,1,1,0,0,0, 0,0,0,0,1,1,1,1,1,1),
#'			     endpoint2=c(0,0,0,0,0,1,0,1,0,0, 0,0,1,1,1,1,1,1,1,1),
#'			     endpoint3=c(0,0,0,0,0,0,0,0,1,1, 0,0,0,1,1,1,1,1,1,1),
#'			     group=rep(c(0,1),each=10))
#'
#' ## greedy alpha exhaustion
#' testgreedy3EP<-mfisher.test(x=data3EP[,1:3],y=data3EP$group,method="alpha.greedy",
#'		show.region=TRUE,closed.test=TRUE)
#' print(testgreedy3EP)
#' par(mfrow=c(3,3))
#' for(i in 1:9) {
#'	plot(testgreedy3EP,dim=c(1,2),slice=list(T3=i),show.titles=FALSE,cex=2,xlim=c(0,8),ylim=c(0,10))
#'	title(paste("T3 =",i))
#' }
#'
#' ## Bonferroni greedy
#' mfisher.test(x=data3EP[,1:3],y=data3EP$group,method="bonferroni.greedy",closed.test=TRUE)
#' ## Bonferroni greedy with alternative input of marginal tables
#' mfisher.test(x=list(table(data3EP$endpoint1,data3EP$group),
#'		table(data3EP$endpoint2,data3EP$group),table(data3EP$endpoint3,data3EP$group)),
#'		method="bonferroni.greedy",closed.test=TRUE)
#'
#' @export
mfisher.test <- function(x,y=NULL,method=c("alpha.greedy","alpha","number","power","bonferroni.greedy"),alpha=0.025,p1=NULL,p0=NULL,max.iter=10^5,limit=0,show.region=FALSE,closed.test=FALSE,consonant=FALSE) {

	method<-match.arg(method,choices=c("alpha.greedy","alpha","number","power","bonferroni.greedy"))

	if(!method%in%c("alpha","number","power","alpha.greedy","bonferroni.greedy")) stop('method must be one of "alpha", "number", "power", "alpha.greedy", "bonferroni.greedy".')
	if(!is.logical(consonant)) stop('consonant must be logical.')
	
	if(is.list(x) & !is.data.frame(x)) {
		if(method!="bonferroni.greedy") stop(paste("Method",method,"not applicable to the provided data."))
		if(any(is.na(x))) stop("No missing values allowed if passing aggregated data.")
		anz.dim<-length(x)
		if(is.null(names(x)) | any(names(x)=="")) EPnames<-paste("EP",1:anz.dim,sep="") else EPnames<-names(x)
		dataTab<-x
		names(dataTab)<-EPnames	
	} else {
		if(is.vector(y) & !is.matrix(y) & !is.array(y)) {
			x<-as.matrix(x)
			if(length(y)!=dim(x)[1]) stop("the number of rows in x must be identical to the length of y.")
			keep<-rowSums(cbind(is.na(x),is.na(y)))==0
			if(any(keep==FALSE)) {
				x<-x[keep,]
				y<-y[keep]
				warning(paste(sum(!keep),"observation(s) removed due to missing values."))
			}
			if(any(x!=1 & x!=0)) stop("x must contain 0 and 1 only.")
			if(any(y!=1 & y!=0)) stop("y must contain 0 and 1 only.")
			anz.dim<-dim(x)[2]
			if(is.null(colnames(x))) EPnames<-paste("EP",1:anz.dim,sep="") else EPnames<-colnames(x)
			basis<-2^c(0:(anz.dim-1))
			category<-factor(x%*%basis,levels=0:(2^anz.dim-1))
			tab<-matrix(table(category,y),ncol=2)
		} else {
			if(any(is.na(x)) | any(is.na(y))) stop("No missing values allowed if passing aggregated data.")
			if(length(dim(x))!=length(dim(y))) stop("x and y must be of same dimension.")
			anz.dim<-length(dim(x))
			tab<-matrix(c(as.vector(y),as.vector(x)),ncol=2)
			x.names<-names(dimnames(x))
			y.names<-names(dimnames(y))
			#if(any(x.names!=y.names)) warning("Names of x and y do not match. Is the data correct?")
			#if(any(x.names=="")) 
			EPnames<-paste("EP",1:anz.dim,sep="") #else EPnames<-x.names
		}
		M<-namen.fn(anz.dim)
		dataTab<-cbind(M,tab)
		colnames(dataTab)<-c(paste("Success_",EPnames,sep=""),"CountGroup0","CountGroup1")
		rownames(tab)<-rownames(dataTab)<-apply(M,1,paste,collapse="")
		
		if(method=="bonferroni.greedy") {
			#marginal table list
			#x als Liste von 2x2 Tafeln. Schema: Zeilen: Failure, Success; Spalten: Ctr, Trt
			x<-vector(mode="list",length=anz.dim)
			for(i in 1:anz.dim) {
				set<-rep(FALSE,anz.dim)
				set[i]<-TRUE
				P<-subset.matrix.fn(set,anz.dim)
				x[[i]]<-P%*%tab
			}			
		}

		tab<-tab[,c(2,1)]
		
		#Elementary test statistics
		M<-t(namen.fn(anz.dim))
		T<-as.numeric(M%*%tab[,1])	

		if(method=="power") if(is.null(p1) | is.null(p0)) stop('p1 and p0 must be provided for method "power".')
		if(!is.null(p1) | !is.null(p0)) {
			#if(length(p1)!=(2^anz.dim) | length(p0)!=(2^anz.dim)) stop("Incorrect length of p1 or p0.")
			if( any( c(dim(p1),dim(p0))!=2) | length(dim(p1))!=anz.dim | length(dim(p0))!=anz.dim) stop("Incorrect dimension of p1 or p0.")
			p1<-as.vector(p1)
			p0<-as.vector(p0)
		}

		if(method%in%c("alpha","number","power","alpha.greedy")) {
			if(consonant==TRUE & anz.dim!=2) {
				consonant<-FALSE
				warning('Consonance constraint is only available for two endpoints.')
			}
			#Permutation distribution
			A<-mult.fisher.2J.fn(tab=tab,p1=p1,p2=p0,calc.marg.p=consonant)
		}
	}

	#Rejection region
	if(method%in%c("alpha","number","power")) {
		R<-optimal.2J.fn(V=A$vek.out,objective=method,alpha=alpha,algorithm="BB",limit=limit,maxiter=max.iter,return.array=FALSE,elementary.consonance=consonant)
		V<-p.val.fun(V=R$V,T=T,anz.dim=anz.dim,consonant=FALSE,alpha=alpha)
		if(R$opt.status==-1) {
			OPT<-list(iterations=0,maxiter=max.iter,optimization.finished=TRUE)
		} else {
			OPT<-list(iterations=R$OPT$iterations,maxiter=R$OPT$maxiter,optimization.finished=R$OPT$optimization.finished)
		}
	}
	if(method=="alpha.greedy") {
		R<-list(V=A$vek.out,anz.dim=anz.dim)
		if(consonant==FALSE) {
			V<-p.val.fun(V=R$V,T=T,anz.dim=anz.dim,all.p.values=show.region,alpha=alpha)
			if(show.region) V$V$solution[!is.na(V$V$p.val) & V$V$p.val<=alpha]<-1
		} else {
			Vtemp<-p.val.fun(V=R$V,T=T,anz.dim=anz.dim,all.p.values=TRUE,consonant=TRUE,alpha=alpha)
			Vtemp$V$solution[!is.na(Vtemp$V$p.val) & Vtemp$V$p.val<=alpha & (Vtemp$V$pval1<=alpha | Vtemp$V$pval2<=alpha)]<-1
			V<-p.val.fun(V=Vtemp$V,T=T,anz.dim=anz.dim,consonant=FALSE,alpha=alpha)
			#if(show.region) V$V$solution[!is.na(V$V$p.val) & V$V$p.val<=alpha]<-1
		}
		OPT<-list(iterations=NA,maxiter=NA,optimization.finished=NA)
	}
	if(method=="bonferroni.greedy") {
		test<-bonferroni.greedy.fn(x)
		V<-list(p.value=test$p)
		OPT<-list(iterations=NA,maxiter=NA,optimization.finished=NA)
		T<-test$stat
	}

	#closed test
	if(closed.test==TRUE & anz.dim==1) {
		warning("Just one endpoint defined. No closed test is calculated.")
		closed.test<-FALSE	
	}
	if(closed.test==TRUE) {
		S<-namen.fn(anz.dim)[-1,] == 1
		k<-dim(S)[1]
		p.intersect<-rep(NA,k)
		opt.intersect<-rep(FALSE,k)
		if(method=="bonferroni.greedy") OR<-test$OR else OR<-rep(NA,anz.dim)
		
		p.intersect[k]<-V$p.value
		opt.intersect[k]<-OPT$optimization.finished

		#layer<-1+rowSums(!S) #there are as many layers as dimensions (Endpoints)
		#with k we get the full table and hence the global intersection test
		for(i in (k-1):1) {

			if(method%in%c("alpha","number","power","alpha.greedy")) {

				#Projection matrix
				P<-subset.matrix.fn(S[i,],anz.dim)
				tab.lokal<-P%*%tab
				if(dim(tab.lokal)[1]==2) {
					p<-fisher.test(tab.lokal[c(2,1),],alternative="greater")$p.value
					p.intersect[i]<-p
					opt.intersect[i]<-TRUE
					OR[(1:anz.dim)[S[i,]]]<-tab.lokal[2,1]/tab.lokal[1,1]*tab.lokal[1,2]/tab.lokal[2,2]
				} else {
					if(!is.null(p1)) p1.lokal<-as.numeric(P%*%p1) else p1.lokal<-NULL
					if(!is.null(p0)) p2.lokal<-as.numeric(P%*%p0) else p2.lokal<-NULL
					A.lokal<-mult.fisher.2J.fn(tab=tab.lokal,p1=p1.lokal,p2=p2.lokal)
					if(method%in%c("alpha","number","power")) {
						R.lokal<-optimal.2J.fn(V=A.lokal$vek.out,objective=method,alpha=alpha,algorithm="BB",limit=limit,maxiter=max.iter,return.array=FALSE)
						opt.intersect[i]<-R.lokal$opt.status
					} 
					if(method=="alpha.greedy") {
						R.lokal<-list(V=A.lokal$vek.out,anz.dim=A.lokal$anz.dim)
						opt.intersect[i]<-NA
					}
					V.lokal<-p.val.fun(V=R.lokal$V,T=T[S[i,]],anz.dim=R.lokal$anz.dim,alpha=alpha)
					p.intersect[i]<-V.lokal$p.value
				}
			}
			if(method=="bonferroni.greedy") {
				test.intersect<-bonferroni.greedy.fn(x,subset=(1:anz.dim)[S[i,]])
				p.intersect[i]<-test.intersect$p
				opt.intersect[i]<-NA
			}
		}
		
		#unadjusted p-values
		p<-p.intersect[rowSums(S)==1]
		#adjusted p-values
		p.adj<-rep(NA,anz.dim)
		for(i in 1:anz.dim) p.adj[i]<-max(p.intersect[S[,i]])
		out.CT<-data.frame(OR,p,p.adj)
		rownames(out.CT)<-EPnames
		CT.details<-as.data.frame(cbind(S,p.intersect))
		names(CT.details)[1:anz.dim]<-EPnames
		if(method%in%c("alpha","number","power")) OPT$optimization.finished<-all(opt.intersect!=0)
	} else {
		out.CT<-NA
		CT.details<-NA
	}

	if(method%in%c("alpha","number","power")) {
		if(OPT$optimization.finished==FALSE) warning("Maximal number of iterations reached, optimization was not finished.")
	}
	
	names(T)<-paste("T",1:anz.dim,sep="")
	if(method%in%c("alpha","number","power","alpha.greedy")) {
		dataTab<-as.data.frame(dataTab)
		dataTab$p0<-p0
		dataTab$p1<-p1
	}
	
	out<-list(call=match.call(),data=dataTab,alpha=alpha,method=method,statistic=T,p.value=V$p.value,conditional.properties=NA,rej.region=NA,elementary.tests=out.CT,closed.test=CT.details,consonance.constraint=consonant,OPT=OPT)
	#Properties of global test conditional on observed data
	if(method%in%c("alpha","number","power","alpha.greedy")) {
		out$conditional.properties<-list(alpha=sum(V$V$alpha[V$V$solution==1]),number=sum(V$V$solution),power=sum(V$V$power[V$V$solution==1]))
	} else {
		out$conditional.properties<-NULL
	}
	if(show.region) {
		if(method=="bonferroni.greedy") {
			warning("show.region ignored when method is bonferroni.greedy.")
		} else {
			out$rej.region<-V$V[c(paste("T",1:anz.dim,sep=""),"solution","alpha","power")]
			colnames(out$rej.region)[anz.dim+1]<-"rejection.region"
		}
	} else {
		out$rej.region<-NULL
	}
	if(closed.test==FALSE) {
		out$elementary.tests<-NULL
		out$closed.test<-NULL
	}
	class(out)<-"multfisher"
	return(out)
} 

#' @title Print Values from a \code{multfisher} Object
#'
#' @description Print the test results.
#'    
#' @param x an object of class \code{multfisher}
#' @param ... further arguments passed to other methods. Not used.
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{mfisher.test}}, \code{\link{plot.multfisher}}
#' @export
print.multfisher<-function(x,...) {
	cat("\n    Optimal exact test for multiple binary endpoints\n\n")
	cat("Aggregated data:\n")
	print(x$data)
	cat("\n")
	if(x$method%in%c("alpha","number","power")) {
		if(x$consonance.constraint==TRUE) constext<-"       with consonance constraint\n" else constext<-""
		cat(paste('Method: Optimization of',paste('"',x$method,'"',sep=""),"for nominal alpha =",x$alpha,"\n",constext))
	} else {
		if(x$consonance.constraint==TRUE) constext<-"with consonance constraint\n" else constext<-"\n"
		cat(paste('Method:',paste('"',x$method,'"',sep=""),constext))
	}
	cat("\n")
	cat(paste("Global test p-value =",round(x$p.value,4),"\n"))
	cat("\n")
	cat("Null hypothesis: Equal success distribution in both groups.\n")
	cat("Alternative hypothesis: Larger success proportion in the treatment group in at least one endpoint.\n")
	if(!is.null(x$elementary.tests)) {
		cat("\nMarginal odds ratios and adjusted p-values:\n")
		x$elementary.tests[,1]<-round(x$elementary.tests[,1],2)
		x$elementary.tests[,2:3]<-round(x$elementary.tests[,2:3],4)
		print(x$elementary.tests)
	}
	cat("\n")
}

#' @title Plot Rejection Region from a \code{multfisher} Object
#'
#' @description Plot two dimensions of the rejection region.
#'    
#' @param x an object of class \code{multfisher}
#' @param dim a vector of length two, indicating which two dimensions of the rejection region to plot. The default is \code{c(1,2)}.
#'	The argument is ignored if \code{x} was calculated for only one endpoint.
#' @param slice a named list of numeric values at which the test statistics for the dimensions (assumedly) not included in \code{dim} are held constant, see details.
#' @param show.titles logical indicating whether the plot title and explanatory subtitles are shown.
#' @param ... further arguments passed to the generic \code{\link[graphics]{plot}} function.
#' @details The function produces plots of the multivariate rejection regions calculated by \code{\link{mfisher.test}}. 
#'	For more than two dimensions, the default \code{slice=NULL} shows a
#'	projection of the rejection region on the two dimensions indicated in \code{dim}. \code{slice} may be specified to produce plots of slices through the 
#'	multivariate rejection region. In that case \code{slice} must be a named list of numeric objects. The names must be of the form \code{Ti}, where \code{i} is replaced
#'	by the number of the dimension to be indicated. The numeric value defines at which value the test statistic for the indicated dimension is held constant.
#'	(If there are dimensions that are neither indicated in \code{dim} nor in \code{slice}, the plot is still a projection.)
#' @author Robin Ristl, \email{robin.ristl@@meduniwien.ac.at}
#' @seealso \code{\link{mfisher.test}}, \code{\link{print.multfisher}}
#'	
#' @examples
#' ## Example with two endpoints
#' data<-data.frame(endpoint1=c(0,0,1,1,1,0,0,0,0,1,1,1,1,1,1, 0,0,1,0,0,1,1,1,1,1,1,1,1,1,1),
#'			endpoint2=  c(0,0,0,0,0,1,1,1,1,1,1,1,1,1,1, 0,0,0,1,1,1,1,1,1,1,1,1,1,1,1),
#'			group=rep(c(0,1),each=15))
#' plot(mfisher.test(x=data[,c(1:2)],y=data$group,show.region=TRUE),cex=2)
#'
#' ## Example with three endpoints
#' data3EP<-data.frame(endpoint1=  c(0,0,0,0,0,1,1,0,0,0, 0,0,0,0,1,1,1,1,1,1),
#'			  endpoint2=  c(0,0,0,0,0,1,0,1,0,0, 0,0,1,1,1,1,1,1,1,1),
#'			  endpoint3=  c(0,0,0,0,0,0,0,0,1,1, 0,0,0,1,1,1,1,1,1,1),
#'			  group=rep(c(0,1),each=10))
#' testgreedy3EP<-mfisher.test(x=data3EP[,1:3],y=data3EP$group,method="alpha.greedy",
#'		show.region=TRUE,closed.test=TRUE)
#' ## Projecion on the first two dimensions
#' plot(testgreedy3EP,dim=c(1,2),cex=2)
#' ## Slice at a value of 5 for the third dimension
#' plot(testgreedy3EP,dim=c(1,2),slice=list(T3=5),cex=2)
#' ## Show all slices through the third dimension
#' par(mfrow=c(3,3))
#' for(i in 1:9) {
#'	plot(testgreedy3EP,dim=c(1,2),slice=list(T3=i),show.titles=FALSE,cex=2,xlim=c(0,8),ylim=c(0,10))
#'	title(paste("T3 =",i))
#'}
#'
#' @export
plot.multfisher<-function(x,dim=c(1,2),slice=NULL,show.titles=TRUE,...) {
	if(is.null(x$rej.region)) stop("Rejection region required. Use mfisher.test with show.region=TRUE.")
	R<-x$rej.region[x$rej.region$alpha>0,]
	pch.region<-ifelse(R$rejection.region==1,19,1)
	if(length(x$statistic)>1) {
		if(length(dim)!=2) stop("dim must be a vector of length two.")
		if(max(dim)>length(x$statistic)) stop("Values in dim must be less or equal the number of dimensions of the rejection region.")
		if(!is.null(slice)) {
			if(!is.list(slice)) stop("slice must be a list.")
			set<-colSums(t(R[,names(slice),drop=FALSE])==unlist(slice))==length(slice)
			if(sum(set)==0) stop("Provided value for slice results in an empty set.")
			R<-R[set,]
			pch.region<-pch.region[set]
		}
		plot(R[,dim], xlab=paste("Statistic endpoint",dim[1]),ylab=paste("Statistic endpoint",dim[2]),pch=pch.region,...)	
			abline(v=unique(R[,dim[1]]),col="grey")
		abline(h=unique(R[,dim[2]]),col="grey")
		if(!is.null(slice) & show.titles==TRUE) mtext(side=1,text=paste("Other statistics held constant at:",paste(names(slice),slice,sep="=",collapse=", ")),line=4,adj=0)
	} else {
		plot(R[,1],rep(1,dim(R)[1]), xlab=paste("Statistic endpoint",dim[1]),ylab="",pch=pch.region,axes=FALSE,...)	
		box()
		axis(1)
			abline(v=unique(R[,dim[1]]),col="grey")
		abline(h=1,col="grey")
	}
	if(show.titles==TRUE) {
		title("Optimal exact test rejection region")
		mtext(side=3,text="Solid: rejection region. Empty: not in rejection region. No dot: not in permutation distribution.")
	}
}


C.ii.fn<-function(ko,sort.it=FALSE, keep.empty=FALSE) {
	#Creates cond(ii) Constraint matrix
	#set keep.empty=TRUE to keep these constraints. This is necessary, if my branch and bound algorithm is used, because
	#it uses the n x n constraint matrix as lookup table for cond(ii).
	#constraint matrix for criterion (ii) (no holes)
	if(!is.matrix(ko) & !is.data.frame(ko)) ko<-matrix(ko,ncol=1)
	if(sort.it==TRUE) for(i in 1:dim(ko)[2]) ko<-ko[order(ko[,i]),] #so it also works for one-dimensional problems
	n<-dim(ko)[1]
	C.ii<-matrix(0,n,n) 
	if(n > 0) {
		not.empty<-rep(TRUE,n)
		for(i in 1:n) {
			for(j in 1:n) if(i!=j & all(ko[j,]>=ko[i,])) C.ii[i,j]<- -1
			summe<- -sum(C.ii[i,])
			
			if(summe>0)	{
				#C.ii.t[,i]<-C.ii.t[,i]/summe
				#C.ii.t[i,i]<-1
				C.ii[i,i]<- summe

			} else {
				if(keep.empty==TRUE) C.ii[i,i]<- -1 else not.empty[i]<-FALSE
				# any number <= 0 is possible,
				#better to remove empty constraints
			}
		}
		C.out<-C.ii[not.empty,] 
	} else {
		C.out<-C.ii
	}
	if(!is.matrix(C.out)) C.out<-matrix(C.out,nrow=n)
	C.out
}

C.ii.t.fn<-function(ko,sort.it=FALSE, keep.empty=FALSE) {
	#set keep.empty=TRUE to keep these constraints. This is necessary, if my branch and bound algorithm is used, because
	#it uses the n x n constraint matrix as lookup table for cond(ii).
	#constraint matrix for criterion (ii) (no holes)
	if(!is.matrix(ko) & !is.data.frame(ko)) ko<-matrix(ko,ncol=1)
	if(sort.it==TRUE) for(i in 1:dim(ko)[2]) ko<-ko[order(ko[,i]),] #so it also works for one-dimensional problems
	n<-dim(ko)[1]
	C.ii.t<-matrix(0,n,n) 
	if(n > 0) {
		not.empty<-rep(TRUE,n)
		for(i in 1:n) {
			for(j in 1:n) if(i!=j & all(ko[j,]>=ko[i,])) C.ii.t[j,i]<- -1
			summe<- -sum(C.ii.t[,i])
			
			if(summe>0)	{
				#C.ii.t[,i]<-C.ii.t[,i]/summe
				#C.ii.t[i,i]<-1
				C.ii.t[i,i]<- summe

			} else {
				if(keep.empty==TRUE) C.ii.t[i,i]<- -1 else not.empty[i]<-FALSE
				# any number <= 0 is possible,
				#better remove empty constraints
			}
		}
		C.out<-C.ii.t[,not.empty] #we calculate and return the transposed constrained matrix 
	} else {
		C.out<-C.ii.t
	}
	if(!is.matrix(C.out)) C.out<-matrix(C.out,nrow=n)
	C.out
}

vek2array.fn<-function(vek,d,name.shift=0) { #previously shift=-1, but we don't need this option anymore
	anz.dim<-length(d)
	dim.namen<-vector(length=anz.dim,mode="list")
	for(i in 1:anz.dim) dim.namen[[i]]<-(1:d[i])+name.shift
	array(vek,dim=d,dimnames=dim.namen)
}

mult.fisher.2J.fn<-function(tab=NULL,n=NULL,m=NULL,M=NULL,p1=NULL,p2=NULL,margin.prob=FALSE,calc.marg.p=FALSE) {
	#eiter a full table must be provided via argument tab
	#or the marginal sums are provided as n and m
	if(!is.null(tab)) {
		m<-rowSums(tab)
		n<-colSums(tab)
	}
	#m...row margins
	#n...column margins (number of observations in the two groups)
	
	#Matrix M to get test statistics as sums of the approproate cell entries
	#M%*%tab[,1] gives the test statistics
	#e.g. for two dimensions M=matrix(c(0,1,0,1,0,0,1,1),byrow=TRUE,nrow=2)
	if(is.null(M)) {
		anz.dim<-round(log2(length(m))) #round just in case there are numerical issues
		M<-t(namen.fn(anz.dim))
	} else {
		anz.dim<-dim(M)[1]
	}
	
	k<-length(m)
	x<-rep(NA,k)
	N<-sum(m)
	cur.lim<-matrix(NA,nrow=k,ncol=2)
	stop<-FALSE
	j<-1:k
	u<-0 #so that we start with the full loop
	ws<-NULL

	#Upper bound for the test statistics
	#this should be improved, so that we look only at the possible range of values for the test statistics
	m.marg<-as.numeric(M%*%m) #marginal sum of successess (succ A, succ B, ...)
	#now use the usual formula for 2x2 tables:
	low.marg<-up.marg<-rep(NA,anz.dim)
	for(i in 1:anz.dim) {
		#lowest value the i-th test stat can take:
		low.marg[i]<-max(0,m.marg[i]-n[2])
		#largest value the i-th test stat can take:
		up.marg[i]<-min(m.marg[i],n[1]) #maximum
	}

	#Coordinates and test statistics have this relation:
	#T_i = X_i + low.marg_i - 1 
	#X_i = T_i + 1 - low.marg_i
	#Tmax<-as.numeric(M%*%m) #this is also marginal number of successess
	#d<-Tmax+1 #in dimension i, the array ranges from 0 to d[i]
	d<- up.marg-low.marg + 1 #this is the size of the nonzero elements of the distribution in each dimension
	#anz.dim<-length(d) #number of dimensions
	#prepare vectorized output
	WA.vek<-Z.vek<-F.vek<-rep(0,prod(d))   # array(0,dim=d) #+1 because 0 is ingeneral a possible value for the test stats
	cp<-cumprod(d[-length(d)]) #requried to set indices
	WS<-NA
	WS.vek<-NA
	if(margin.prob==TRUE) ws.obs.margins<-0 else ws.obs.margins<-NA #probability, to observe a table with the provided margins, this may be used when calculating unconditional power
	#log.total.perm<-lchoose(sum(n),n[1])
	total.perm<-choose(sum(n),n[1])

	while(!stop) {
		#for i>u: calculate cur.lims anew and set x[i] to the lower limit. x are the observed frequencies for group 1 in the current permutation
		if(u+1<=k-1) {
			for(i in (u+1):(k-1)) {
				x.sum<-ifelse(i==1,0,sum(x[1:(i-1)])) #x sum set so far #sum(c(0,x[0:(i-1)]))
				m.sum<-sum(m[(i+1):k]) #remaining margin sum
				low<-max(n[1]-x.sum-m.sum,0) #c(0,x[0])=0
				up<-min(n[1]-x.sum,m[i])
				cur.lim[i,]<-c(low,up)
				x[i]<-low
			}
		}
		x[k]<-n[1]-sum(x[1:(k-1)])
		cur.lim[k,]<-c(x[k],x[k]) #here, only one possibility remains
		#Calculate probabiliy for current permutation
		mgl<-prod(choose(m,x))
		#ws<-c(ws,mgl)

		#using logs also works, but it seems to be less exact. Is it faster?
		#log.mgl<-sum(lchoose(m,x))
		#ws<-c(ws,log.mgl)

		#Current Vector of test statistics (now not plus 1, because the coordinates are calc below
		T.cur<-M%*%x # + 1
		#Current vector of coordinates X_i = T_i + 1 - low.marg_i
		X.cur<-T.cur + 1 - low.marg
		
		#Coordinates in vectorized matrix (also im Vektor)
		ko<-sum(cp*(X.cur[-1]-1))+X.cur[1]
		
		Z.vek[ko]<-Z.vek[ko]+mgl
		#F.vek[ko]<-F.vek[ko]+exp(log.mgl-log.total.perm)
		#F.vek[ko]<-F.vek[ko]+mgl/total.perm

		#old version using matrix: #Z[T.A+1,T.B+1]<-Z[T.A+1,T.B+1]+mgl
		
		#WS, Probability under the alternative
			if(!is.null(p1) & !is.null(p2) ) {
				x2<-m-x
				p.u.log<-sum(x*(log(p1)-log(p2)))-sum(lfactorial(x))-sum(lfactorial(x2))
				WA.vek[ko]<-WA.vek[ko]+exp(p.u.log)

				#unconditional probability of these table margins under given alternative p1 and p2:
				if(margin.prob==TRUE) ws.obs.margins<-ws.obs.margins+exp(dmultinom(x=x,prob=p1,log=TRUE)+dmultinom(x=x2,prob=p2,log=TRUE))
			
			}

		#Stop or next iteration
		if(all(x==cur.lim[,2])) {
			stop<-TRUE
		} else {
			#increase x[i] by 1, where i=max{j=1,...,k : x[j]<cur.lim[j,2]}
			u<-max(j[x<cur.lim[,2]])
			x[u]<-x[u]+1
		}
	}

		
	#Distribution under H0
	F.vek<-Z.vek/sum(Z.vek)
	#F.vek<-Z.vek/total.perm

	#done above, using log/log.total

	#make dimnames, and in the same loop, make the coordinate vectors X1, X2, ... for the vectorized output
	dim.namen<-vector(length=anz.dim,mode="list")
	koord<-matrix(NA,ncol=anz.dim,nrow=prod(d))
	colnames(koord)<-paste("X",1:anz.dim,sep="")

	#caclulate marginal p-values (one sided), these are used for the p-value based rejection regions
	marg.pval<-matrix(NA,ncol=anz.dim,nrow=prod(d))
	colnames(marg.pval)<-paste("pval",1:anz.dim,sep="")

	for(i in 1:anz.dim) {
		dim.namen[[i]]<-low.marg[i]:up.marg[i]     #0:(d[i]-1)

		#Coordinates, use recycling of vectors in R
		koord[,i]<-rep(1:d[i],each=prod( c(1,d)[1:i] ) )
			#c(1,d), so that in the first dimension each entry is repeated once, and otherwise the number of repeats is the number of values that occur in the previous dimension
		
		#and the marginal p-values:
		if(calc.marg.p==TRUE) {		
			t.temp<-low.marg[i]:up.marg[i]    #0:Tmax[i] #ist 0:(d[i]-1), Tmax ist auch gleichzeitig marginal number of successes
			p.temp<-phyper(q=t.temp,m=n[1],n=n[2],k=m.marg[i],lower.tail=FALSE)+
					dhyper(x=t.temp,m=n[1],n=n[2],k=m.marg[i])
			marg.pval[,i]<-rep(p.temp,each=prod( c(1,d)[1:i] ) )   #c(1,d), as above
		}
	}

	#Distributions under H0 and alternative in array form:
	#null distribution
	F<-array(F.vek,dim=d,dimnames=dim.namen)
	#Joint distribution of T under the alternative p1,p2
	if(!is.null(p1) & !is.null(p2) ) {
		WS.vek<-WA.vek/sum(WA.vek)
		WS<-array(WS.vek,dim=d,dimnames=dim.namen)
	}
	
	#Vektorized output, the LP requires this type of input anyways and it works for arbitrary dimensions
	#Teststatistics from coordinates
	#T_i = X_i + low.marg_i - 1 
	T<-matrix(NA,nrow=dim(koord)[1],ncol=dim(koord)[2])
	for(i in 1:dim(koord)[2]) T[,i]<-koord[,i] + low.marg[i] -1

	colnames(T)<-paste("T",1:anz.dim,sep="")
	vek.out<-data.frame(koord,T,index=1:prod(d),alpha=F.vek,power=WS.vek,solution=0,E=0,marg.pval)
	vek.out	

	if(!is.null(tab)) obs.prob<-prod(choose(m,tab[,1])) else obs.prob<-NA #number of Permutations that give the observed table (if provided)
	list(vek.out=vek.out,d=d,anz.dim=anz.dim,M=M,n=n,m=m,marginal.succ=m.marg,F=F,WS=WS,p1=p1,p2=p2,ws.obs.margins=ws.obs.margins,total.perm=total.perm,sumZ=sum(Z.vek),sum(Z.vek)==total.perm,dim.out=dim(vek.out),dim.koord=dim(koord),l.F.vek=length(F.vek))
}

optimal.2J.fn<-function(V,preproc.results=NULL,report.pp=TRUE,elementary.consonance=FALSE,objective=c("alpha","power","area","specified")[3],obj.spec=NULL,true.power=NULL,anz.dim=sum(grepl("X",names(V))),alpha=0.025,limit=0,limit.obj=0,allowed=NULL,koo.cur=NULL,subset=NULL,return.array=FALSE,rescale.obj=10000,algorithm=c("LP","Rglpk","BB")[1],keep.empty=(algorithm=="BB"),presolution=NULL,maxiter=100000,...) {
	if(objective=="number") objective<-"area"
	#if preproc.results, which contains all results of the preprocessing, the function skips the preprocessing, 
	#and directly starts the optimization.
	#This is much faster, because the preprocessing can take a second or so.
	#Also any parts referring to "allowed" are skipped, as we assume that the supplied X3 and C.mat
	#already are matching any such restrictions.
	#Take care: the setting for limit must be the same as in the calculation of the preprocessing

	#rescale.obj is a factor by which the objective contributions are multiplied when the objective is
	#aplha or power. Otherwise the lpsolver sometimes runs into problems.
	#limit is the limit under which contributions to alpha are set to zero (and alpha is lowered accordingly) to speed up computation
	#limit.opt is the limit under which contributions to the objective function are considered too small to allow stable computation
	#all such values are set to the value of limit.opt
	
	#koo.cur ... current coordinates and solutions of lower level tests in the closed test
	#points that are not contained in any one-step-lower level rejection region are set to "not allowed",
	#this enforces consonance
	#is used after pre-processing step 1, to make it faster
	
	#With a specified objective, e.g. power other than that in V can be used.
	#We could also use true.power, but take care: If preproc.result is used, the V will be supplied from there
	#so V must match the optimization aim.
	if(objective=="specified") {
		if(is.null(obj.spec)) stop("Objective vector obj.spec is required.")
		V$specified<-obj.spec
	}

   if(is.null(preproc.results)) {
	preproc.results<-NA

	#check general argument "allowed"
 	if(is.null(allowed)) allowed<-rep(TRUE,dim(V)[1]) #"allowed" can be used to exclude points and enforce consonance
	if(length(allowed)!=dim(V)[1]) stop("Vector 'allowed' not of same length as data.")

	#Add information on the level of the subset test in the closed test
	if(!is.null(subset)) {
		dim.sub<-sum(subset)
		V$subset.level<-dim.sub
	}
	
	#Check consonance constraint due to solution of lower level tests (the locally optimal bottom up approach).
	#This overrides the argument "allowed".
	#if(!is.null(koo.cur) & !is.null(subset)) {
	#	anz.dim.total<-sum(grepl("X",names(koo.cur)))
	#	#these are the points in the relevant dimensions that are lower level solutions
	#	G<-koo.cur[dim.sub-koo.cur$subset.level==1 & koo.cur$solution==1 , (1:anz.dim.total)[subset]]
	#	allowed<-rep(FALSE,dim(V)[1])
	#	#check for each relevant dimension, which points have lower level solutions
	#	if(dim(G)[1]>0 & dim(V)[1]>0) {  #otherwise there are no lower level solution points or no points in V
	#		for(i in 1:dim(V)[1]) {
	#			for(j in 1:dim(G)[1]) {
	#				delta<-V[i,1:anz.dim] - G[j,]
	#				if( all(delta[!is.na(delta)] )==0 ) allowed[i]<-TRUE
	#			}
	#		}
	#	}
	#
	#	#for(i in 1:dim.sub) allowed<-allowed & (V[,i]%in%G[,i]  )  #all subdimension coordinates must match
	#}

	#Check if consonance due to elementary tests is demanded. This uses the pvalues supplied with V
	#This step does not overrides the argument koo.cur. It may used inaddition to koo.cur, because then X is smaller
	#when it comes to checking the consonance with koo.cur
	if(elementary.consonance==TRUE) {
		allowed<-rep(TRUE,dim(V)[1])
		#only works for EP>1, otherwise PVALS will be a vector, not a matrix
		PVALS<-V[,grepl("pval",names(V))]
		for(i in 1:dim(V)[1]) if(all(PVALS[i,]>alpha)) allowed[i]<-FALSE
	}
	
	V$allowed<-allowed
	#plot.muf.fn(F=A$F,region=vek2array.fn(V$allowed,d=A$d),gr=1.25,r=1,main="Optimal power cons",print.area="prob.area",gr.label=1,label.line=4,gr.ax=1.25)


	#1) remove lines with prob 0 and lines that are excluded because of some external constraint (e.g. to enforce consonance)
	keep<-V$alpha>0 & V$allowed==TRUE
	#all((V$alpha>0) == (V$alpha!=0)) #ok
	X<-V[keep,]
	
    #check, if X is empty. If so, skip remaining parts and report V. We need to check this twice, to avoid a conflict with the preproc.results jump
    if(dim(X)[1]>0) {

	#2) actual pre-processing, cond (i) und (ii)
	K<-dim(X)[1] #number of points (rows)
	#condition (i)
	for(i in 1:K) {
		set<-rep(TRUE,K)
		for(j in 1:anz.dim) set<-set&(X[,j]>=X[i,j]) #the coordinates have to be in the first columns of V
		alpha.min<-sum(X$alpha[set])
		if(alpha.min<=alpha) X$E[i]<-1 #eligible points, according to cond(i)
	}
	#only keep eligible
	X1<-X[X$E==1,]
	K1<-dim(X1)[1]
	

	#here is the better place for testing the consonance contraint, because this can be slow
	#so it is better to have it after pp step 1 (it is much faster this way)
	#Check consonance constraint due to solution of lower level tests (the locally optimal bottom up approach).
	#This is in addition to the argument "allowed".
	if(!is.null(koo.cur) & !is.null(subset)) {
		anz.dim.total<-sum(grepl("X",names(koo.cur)))
		#these are the points in the relevant dimensions that are lower level solutions
		G<-koo.cur[dim.sub-koo.cur$subset.level==1 & koo.cur$solution==1 , (1:anz.dim.total)[subset]]
		#allowed<-rep(FALSE,dim(V)[1])
		#check for each relevant dimension, which points have lower level solutions
		if(dim(G)[1]>0 & K1>0) {  #otherwise there are no lower level solution points or no points in X1
			allowed.X1<-rep(FALSE,K1)
			for(i in 1:K1) {
				for(j in 1:dim(G)[1]) {
					delta<-X1[i,1:anz.dim] - G[j,]
					if( all(delta[!is.na(delta)] )==0 ) allowed.X1[i]<-TRUE
				}
			}
		}
		X1<-X1[allowed.X1,]
		K1<-dim(X1)[1]
		#for(i in 1:dim.sub) allowed<-allowed & (V[,i]%in%G[,i]  )  #all subdimension coordinates must match
	}
	
	max.alpha<-sum(X1$alpha)

	#condition (ii)
	if(K1>0) {
		for(i in 1:K1) {
			set<-rep(TRUE,K1)
			for(j in 1:anz.dim) set<-set&(X1[,j]<=X1[i,j])
			max.level.ohne.i<-max.alpha-sum(X1$alpha[set])
			max.level.mit.i<-max.level.ohne.i + X1$alpha[i]
			if(max.level.mit.i<=alpha) X1$E[i]<-2 #points, which are contained in eversy optimal rej. region
		}
		X2<-X1[X1$E==1,]
	} else {
		X2<-X1
	}
	alpha.frei<-alpha-sum(X1$alpha[X1$E==2])

	
	#3) To speed up comp.: alpha very close to zero are treated as 0
	#In the null distribution, one can set to zero all probabilities that are very small (close to zero)
	#such that the sum of these probabilities is below some pre-specified threshold, e.g. 10^-4. These points are then removed from the search space and the optimization is done for a nominal level alpha minus the threshold value.
	#In the end all these "almost zero prob." points can be added for free, subject to the "no holes" constraint. 
	
	K2<-dim(X2)[1]
	anz.quasi.0<-0 #this is required later, even if limit==0
	if(limit>0 & alpha.frei>limit & K2>0) {
		rk<-rank(X2$alpha, ties.method = "first")
		cums<-cumsum(sort(X2$alpha))
		anz.quasi.0<-sum(cums<=limit)
		if(anz.quasi.0>0) {
			set.not.0<- rk>anz.quasi.0 #I.e, those with rank larger than the number of quasi 0 elements
			X3<-X2[set.not.0,]
			alpha.new<-alpha.frei-cums[anz.quasi.0]
		} else {
		set.not.0<-NA
		X3<-X2
		alpha.new<-alpha.frei
		}
	} else {
		set.not.0<-NA
		X3<-X2
		alpha.new<-alpha.frei
	}
	
	K3<-dim(X3)[1]

	#4) Make constraint matrix
	#Criterion (i) is alpha
	#criterion (ii) is no holes:

	#C.t<-C.ii.t.fn(ko=X3[,1:anz.dim]) 
	C.mat<-C.ii.fn(ko=X3[,1:anz.dim],keep.empty=keep.empty) 
	
	#save preproc results
	if(report.pp) preproc.results<-list(X=X,X1=X1,X2=X2,X3=X3,K=K,K3=K3,keep=keep,alpha.new=alpha.new,C.mat=C.mat,set.not.0=set.not.0,anz.quasi.0=anz.quasi.0)
    } else {  #brackets closed from first check if X is empty
	if(report.pp) preproc.results<-list(X=X,X1=NA,X2=NA,X3=NA,K=NA,K3=NA,keep=keep,alpha.new=NA,C.mat=NA,set.not.0=NA,anz.quasi.0=NA)
    }
	

   } else { #brackets closed from checking if the preprocessing results were supplied to the function
	#use preproc.results
	X<-preproc.results$X
	X1<-preproc.results$X1
	X2<-preproc.results$X2
	X3<-preproc.results$X3
	K<-preproc.results$K
	K3<-preproc.results$K3
	keep<-preproc.results$keep
	alpha.new<-preproc.results$alpha.new
	C.mat<-preproc.results$C.mat
	set.not.0<-preproc.results$set.not.0
	anz.quasi.0<-preproc.results$anz.quasi.0

   }

    opt.status<- -1 #this means, if the preprocessing does everything (and this is possible) we see -1. 0 means optimization successfull, 1 means optimization not successfull
	OPT<-NA #in case there is no optimization
    #check, if X is empty second time. If so, skip remaining parts and report V. We need to check this twice, to avoid a conflict with the preproc.results jump
    if(dim(X)[1]>0) {


	#5) Optimization with LP, Rglpk or branch and bound
	if(K3>0) {
		C.all<-rbind(X3$alpha,C.mat)
		#rhs
		rhs<-c(alpha.new,rep(0,dim(C.mat)[1]))
	
		X3$area<-1
		#this optimizes area, and if there is more than one solution picks that with best alpha exhaustion
		X3$area.and.alpha<- 1 + X3$alpha 
		if(objective=="specified" & !is.null(preproc.results)) X3$specified<-obj.spec[X3$index]
		objs<-X3[[objective]]
		if(limit.obj>0) objs[objs<=limit.obj]<-limit
		if(objective!="area" & objective!="area.and.alpha") objs<-objs*rescale.obj
			
	
		#Optimization
		#if(algorithm=="LP") OPT<-lp(direction="max",objective.in=objs,const.mat=C.all, const.dir="<=" ,const.rhs=rhs,transpose.constraints=TRUE,all.bin=TRUE,num.bin.solns=1)
		#if(algorithm=="Rglpk") OPT<-Rglpk_solve_LP(obj=objs, mat=C.all, dir=rep("<=",dim(C.all)[1]), rhs=rhs, types = "B", max = TRUE)
		if(algorithm=="BB") OPT<-bb.fn(X=X3,C.mat=C.mat,obj=objs,alpha.new=alpha.new,presolution=presolution,maxiter=maxiter)
		if(algorithm=="none") OPT<-list(solution=X3$solution,status=NA)
		X3$solution<-OPT$solution
		#opt.status<-OPT$status
		opt.status<-OPT$optimization.finished

	}

	#Write the solution to V and X
	added.by.lp<-X3$index[X3$solution==1]
	added.by.pp<-X1$index[X1$E==2] 
	total<-c(added.by.lp,added.by.pp) #they are all unique
	V$solution[total]<-1 #because the rownumbers "index" apply to the original table V
	X$solution<-V$solution[keep] #keep are the points with truly alpha>0 haben and which are "allowed"
	#alternatively we could use X$solution[X$index%in%total]<-1 #I don't know which is faster

	#Fill up quasi alpha 0 points (this has to be done, otherwise there can be holes)
	#I.e. all "larger" points have to be in the solution or they are quasi 0 points themselves
	if(limit>0 & anz.quasi.0>0) {
		V$quasi.0<-0
		V$quasi.0[X2$index[!set.not.0]]<-1
		X$quasi.0<-V$quasi.0[keep] #Work at level of X, because it has less rows than V
		for(i in 1:K) {
			set<-rep(TRUE,K)
			for(j in 1:anz.dim) set<-set&(X[,j]>=X[i,j])
			if(all(X$solution[set]==1 | X$quasi.0[set]==1)) X$solution[i]<-1
		}
		V$solution[X$index[X$solution==1]]<-1
	}

     } #bracket closed from parts that are only executed if X is not empty, second time
	#Output:
	RX<-X$solution==1
	level<-sum(X$alpha[RX])
	if(is.null(true.power)) power<-sum(X$power[RX]) else power<-sum(true.power[keep][RX])
	area<-sum(RX)
	R<-NA
	if(return.array==TRUE) R<-vek2array.fn(V$solution,...)
	if(report.pp==FALSE) preproc.results<-NA

	list(V=V,OPT=OPT,preproc.results=preproc.results,anz.dim=anz.dim,level=level,power=power,area=area,R=R,opt.status=opt.status,algorithm=algorithm)
}

dec2bin<-function(x,digits=8) {
	b<-rep(0,digits)
	for(k in 1:digits) {
		b[k]<-x%%2	
		x<-(x - b[k])/2
	}
	b
}

#namen.fn<-function(anz.dim) {
#	#create binary numbers from 0 to 2^anz.dim - 1
#	x<-0:(2^anz.dim - 1)
#	x
#	M<-t(sapply(x,dec2bin,digits=anz.dim))
#	#sort, so that we get 0,0,0 ; 1,0,0 ; 0,1,0 ; ...
#	if(anz.dim>1) out<-M[order(rowSums(M)),] else out<-t(M)  #because in this case the t() in sapply is not required
#	out 
#}

namen.fn<-function(anz.dim) {
	#create binary numbers from 0 to 2^anz.dim - 1
	x<-0:(2^anz.dim - 1)
	M<-t(sapply(x,dec2bin,digits=anz.dim))
	if(anz.dim>1) out<-M else out<-t(M)  #because in this case the t() in sapply is not required, no sorting
	out 
}



subset.matrix.fn<-function(subset,anz.dim) {  #subset is a logical vector
	M<-namen.fn(anz.dim)
	n.sub.dim<-sum(subset)
	Msub<-namen.fn(n.sub.dim)
	nr<-2^n.sub.dim
	nc<-2^anz.dim
	K<-matrix(NA,nrow=nr,ncol=nc)
	for(i in 1:nr) {
		for(j in 1:nc) K[i,j]<-as.numeric(all(Msub[i,]==M[j,subset]))
	}
	K
}

bb.fn<-function(X,C.mat,obj,alpha.new,presolution=NULL,maxiter=10000) {
	C.1<-(C.mat!=0)

	#X is a data.frame such as the output V of mult.fisher.2J.fn,
	#possibly after pre-processing
	#C.mat is the constraint matrix, C.t is the transposed constraint matrix (transposed, because for lpSolve this is
	#supposed to be faster and we want to use the code we already have).
	#presolution is a feasible (probably non-optimal) solution that is supplied.
	#It can be found be faster algorithms such as marg.fisher.opt.fn or pval.regions.fn.
	#obj is the vector of contributions to the objective function (obj%*%solution is the value of the objective function)

	#set length of the solution vector
	d.x<-dim(X)[1]
	#check presolution and set current best lower bound using the presolution
	if(!is.null(presolution)) {
		max.low<-sum(obj[presolution==1])
		if(sum(X$alpha[presolution==1])>alpha.new | !all(C.mat%*%presolution <=0)) {
			#all(t(presolution)%*%C.t <=0)
			warning("Presolution not feasible.")
			presolution<-NULL
		}
	}
	if(is.null(presolution)) {
		presolution<-rep(0,d.x)
		max.low<-0
	}
	
	#the function to add a node:
	index<-1:d.x
	add.fn<-function(n,z) {
		x<-n[-c(1:3)]
		#find first not yet defined entry
		i<-min(index[x == -1])
		#set new entry and all others according to cond(ii) and return low, up, branched out, x, if cond(i) is met
		if(z==0) {
			x[C.1[,i]]<-0
			return(c(sum(obj[x==1]),sum(obj[x!=0]),all(x>=0),x))
		}
		if(z==1) {
			x[C.1[i,]]<-1
			#check cond(i), this is only required for z=1.
			if(sum(X$alpha[x==1])>alpha.new) return(NULL) else return(c(sum(obj[x==1]),sum(obj[x!=0]),all(x>=0),x))
		}
	}	

	#prepare matrix of nodes
	#the columns are: lower bound, upper bound, fully branched, and then the solution vector
	#for a branched out solution, low=up.
	#for a solution 1 means included, 0 means not included and -1 means not yet determined

	N<-matrix(
		c(max.low,max.low,1,presolution,
		0,0,0,rep(-1,d.x)),byrow=TRUE,nrow=2)
	
	nc<-dim(N)[2]		
	active<-N[,3]==0
	N.act<-matrix(N[active,],ncol=nc)
	N.finished<-matrix(N[!active,],ncol=nc)

	iter<-0L
	#maxiter<-807
	while(dim(N.act)[1]>0 & iter<maxiter) {
		iter<-iter+1L
		#select of all not fully branched nodes that with largest lower bound 
		#as alternative we might alternatingly also select that with largest upper bound.
		#so they get removed before everything has branched out.
		#In the example, this takes the same number of iterations, but is slower, probably because N can gets larger.
		#index 1 or 2:
		#indi<-1+iter%%2
		indi<-1
		j<-which(N.act[,indi]==max(N.act[,indi]))[1]
		n<-N.act[j,]
		#add new nodes and remove the one just branched (j), rbind(NULL,n) gives a matrix, so it N always stays a matrix
		N<-rbind(N.finished,N.act[-j,],add.fn(n,0),add.fn(n,1))
		#update current best lower bound
		max.low<-max(N[,1])
		#remove all nodes with up<max.low, N can become 0 x (3+d.x) Matrix, but it will not disappear
		N<-matrix(N[N[,2]>=max.low,],ncol=nc)
		#select active nodes
		active<-N[,3]==0
		N.act<-matrix(N[active,],ncol=nc)
		N.finished<-matrix(N[!active,],ncol=nc)


	}
	optimization.finished<-dim(N.act)[1]==0
	if(optimization.finished) {
		solution<-N[1,-c(1:3)]
	} else {
		#else we still provide a feasible solution, namely that with the largest lower bound
		j<-which(N[,1]==max(N[,1]))[1]
		solution<-N[j,-c(1:3)]
		solution[solution==-1]<-0 #this is always feasible because of the filling up in add.fn
	}
	list(solution=solution,obj.value=sum(obj[solution==1]),nodes=N,n.solutions=dim(N)[1],iterations=iter,maxiter=maxiter,optimization.finished=optimization.finished,status=as.numeric(!optimization.finished))
}

check.cond.ii.fn<-function(ind,VV,anz.dim) {
	#in each direction, the neighboring point needs to be in the solution, or be of prob zero	
	set<-rep(TRUE,dim(VV)[1])
	for(i in 1:anz.dim) set<-set&(VV[,i]>=VV[VV$index==ind,i])
	cond.ii<-sum(VV$solution.temp[set]) == (sum(set)-1)  #-1, because the point itself is in the set, but not yet in the solution
	#cond.i<- (sum(V1.s$alpha[V1.s$solution==1])+V1.s$alpha[V1.s$index==ind])  <= alpha
	#cond.i & cond.ii
	cond.ii
}

check.cond.ii.remove.fn<-function(ind,VV,anz.dim) {
	#in each direction, the neighboring point needs to be in the solution, or be of prob zero	
	set<-rep(TRUE,dim(VV)[1])
	for(i in 1:anz.dim) set<-set&(VV[,i]<=VV[VV$index==ind,i])
	cond.ii<-sum(VV$solution.temp[set]) == 1  #1, because the point itself is in the set
	#cond.i<- (sum(V1.s$alpha[V1.s$solution==1])+V1.s$alpha[V1.s$index==ind])  <= alpha
	#cond.i & cond.ii
	cond.ii
}

p.val.fun<-function(V,T,anz.dim,all.p.values=FALSE,consonant=FALSE,alpha) {
	V$p.val<-NA
	if(consonant==TRUE) {
		V$allowed<-rep(TRUE,dim(V)[1])
		#only works for EP>1, otherwise PVALS will be a vector, not a matrix
		PVALS<-V[,grepl("pval",names(V))]
		for(i in 1:dim(V)[1]) if(all(PVALS[i,]>alpha)) V$allowed[i]<-FALSE
	}
	
	#V should contain the vector solution. If not, the p-value calculation just gives the greedy test.
	#remove points with prob zero
	V1<-V[V$alpha>0,]

	

	indT<-V$index[colSums(t(V[,(anz.dim+1):(2*anz.dim)])==T)==anz.dim]
	
	#start with points that are not in the rejection region
	ZZ0<- sum(1-V1$solution)
	if(V[indT,"solution"]==0 | (all.p.values==TRUE & ZZ0>0) ) {
		#sort 
		if(consonant==FALSE) {
			V.s<-V1[order(V1$alpha),] #sort by allowed, in order to account for the consonance constraint
		} else {
			V.s<-V1[order(!V1$allowed,V1$alpha),]
		}
		p.curr<-sum(V.s$alpha[V.s$solution==1])
		V.s$solution.temp<-V.s$solution
		indices<-V.s$index[V.s$solution==0]
		outer.stop<-FALSE
		zz<-0
		while(!outer.stop) {
			zz<-zz+1
			Z<-length(indices)
			stop<-FALSE
			z<-0
			#ind<-min(ind.s[V.s$solution.temp==0]) #sum(V.s$solution)+1
			while(!stop) {
				z<-z+1
				#V.s[ind,]
				if(check.cond.ii.fn(indices[z],V.s,anz.dim)==TRUE) {
					p.curr<-p.curr+V$alpha[indices[z]]
					V$p.val[indices[z]]<-p.curr
					V.s$solution.temp[V.s$index==indices[z]]<-1
					if(indices[z]==indT & all.p.values==FALSE) outer.stop<-TRUE
					indices<-indices[-z]
					stop<-TRUE
				}
				#should not be necessary:
				if(z==Z) stop<-TRUE
			}
		if(zz==ZZ0) outer.stop<-TRUE
		}
	}

	#now regard points that are in the rejection region
	ZZ1<-sum(V1$solution)
	if(V[indT,"solution"]==1 | (all.p.values==TRUE & ZZ1>0)) {
	#if(number.in.solution>0) {
		V.s<-V1[order(V1$alpha,decreasing=TRUE),]
		#V.s$counter<-NA
		p.curr<-sum(V.s$alpha[V.s$solution==1])

		V.s$solution.temp<-V.s$solution

		indices<-V.s$index[V.s$solution==1]
		#for(zz in 1:(sum(V.s$solution))) {
		outer.stop<-FALSE
		zz<-0
		
		while(!outer.stop) {
			Z<-length(indices)
			stop<-FALSE
			z<-0
			while(!stop) {
				z<-z+1
				#V.s[ind,]
				if(check.cond.ii.remove.fn(indices[z],V.s,anz.dim)==TRUE) {
					V$p.val[indices[z]]<-p.curr
					p.curr<-p.curr-V$alpha[indices[z]]
					V.s$solution.temp[V.s$index==indices[z]]<-0
					#V.s$counter[V.s$index==indices[z]]<-zz
					if(indices[z]==indT & all.p.values==FALSE) outer.stop<-TRUE
					indices<-indices[-z]
					stop<-TRUE
				}
				#should not be necessary:
				if(z==Z) stop<-TRUE
			}
		if(zz==ZZ1) outer.stop<-TRUE
		}
	}
	list(V=V,indT=indT,p.value=V$p.val[indT])
}

dhyper.fn<-function(n1.,n.1,n..) {
	n2.<-n..-n1.
	n.2<-n..-n.1
	start<-max(0,n.1-n2.)
	ende<-min(n1.,n.1)
	prob<-dhyper(x=start:ende,m=n1.,n=n..-n1.,k=n.1)	
	out<-data.frame(t=start:ende,prob=prob)
	out
}

bonferroni.greedy.fn<-function(x,subset=1:length(x)) {
	anz.dim<-length(subset)
	tabList<-vector(mode="list",length=anz.dim)
	T<-OR<-rep(NA,anz.dim)
	zeile<-rep(NA,anz.dim)
	
	for(i in 1:anz.dim) {
		subtab<-x[[subset[i]]][c(2,1),c(2,1)]
		T[i]<-subtab[1,1]
		OR[i]<-subtab[1,1]/subtab[2,1]*subtab[2,2]/subtab[1,2]
		tabList[[i]]<-dhyper.fn(n1.=sum(subtab[1,]),n.1=sum(subtab[,1]),n..=sum(subtab))
		zeile[i]<-dim(tabList[[i]])[1]
	}
	p<-0
	stop<-FALSE
	while(!stop) {
		eligible<-rep(NA,anz.dim)
		for(i in 1:anz.dim) eligible[i]<-tabList[[i]][zeile[i],2]
		j<-which(eligible==min(eligible))[1]
		p<-p+tabList[[j]][zeile[j],2]
		if(tabList[[j]][zeile[j],1]==T[j]) stop<-TRUE else zeile[j]<-zeile[j]-1
	}
	list(stat=T,OR=OR,p=p)
}


