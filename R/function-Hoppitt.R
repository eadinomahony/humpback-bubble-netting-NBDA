#A collection of functions for detecting social transmission from order of acquisition and time of acquisition data

#For fitting of Cox models
require(survival)
require(combinat)


setClass("oaData",representation(idname="vector",assMatrix="matrix",asoc="matrix",orderAcq="vector",ties="vector",trueTies="list",groupid="character",taskid="character",coxdata="data.frame",mldata="data.frame"));

setMethod("initialize",
          signature(.Object = "oaData"),
          function (.Object, idname=NULL, assMatrix,asoc,orderAcq,ties=NULL,trueTies=list(NULL),groupid="1",taskid="1",...) 
          {
            #Create default names vector if none is provided, if it is, convert to a factor
            if(is.null(idname)) idname<-(1:dim(assMatrix)[1]);
            
            
            time1<-vector();
            time2<-vector();
            status<-vector();
            
            
            
            id<-vector();
            identity<-vector();
            stMetric<-vector();
            totalMetric<-vector();
            learnMetric<-vector();
            totalAsoc<-vector();
            
            #Make sure the asoc is in matrix form
            asoc<-as.matrix(asoc);
            
            #Calculate the asocial learning variables for the learning individual at each step
            learnAsoc<-matrix(nrow=dim(asoc)[2],ncol=length(orderAcq))
            for(i in 1:length(orderAcq)) learnAsoc[,i]<-asoc[orderAcq[i]]
            
            #If there are no ties vector of zeroes
            if(is.null(ties)) ties<-rep(0,length(orderAcq));
            
            #Define functions for calculating total st Metric
            newFunc<-function(x) x*(1-statusTracker)
            newFunc2<-function(x) x*statusTracker2
            
            
            #Generate a variable for tracking the status of individuals in the group in question
            statusTracker<-rep(0,dim(assMatrix)[1]);
            statusTracker2<-rep(0,dim(assMatrix)[1]);
            
            #Loop through acquisition events
            for(i in 1:length(orderAcq)){
              
              if (ties[i]==0)statusTracker2<-statusTracker;
              
              #Loop through each individual
              for(j in 1:dim(assMatrix)[1]){
                
                #Only naive individuals considered
                if(statusTracker[j]==0){
                  
                  #Record variables for Cox model
                  time1<-c(time1,i-1);
                  time2<-c(time2,i);
                  identity<-c(identity,j);
                  stMetric<-c(stMetric,sum(assMatrix[j,]*statusTracker2));
                  id<-c(id,idname[j]);
                  
                  #Record status as one if individual acquires trait, zero otherwise for Cox model
                  if(j==orderAcq[i]){
                    status<-c(status,1);
                    
                    #Record the social transmission metric and asocial learning variables for the learning individual (ML model)
                    learnMetric<-c(learnMetric,sum(assMatrix[j,]*statusTracker2));
                    
                  }else{
                    status<-c(status,0);
                    
                  }
                }
              }
              
              #Calculate total st metric over all individuals in the group for each acquisition event (ML model)
              newMatrix<-apply(assMatrix,2,newFunc)
              totalMetric<-c(totalMetric,sum(apply(newMatrix,1,newFunc2)));
              
              #Set statusTracker to one if individual acquires trait
              statusTracker[orderAcq[i]]<-1;
            }
            
            #Record individual variables for each line of data in the Cox model
            indVar<-matrix(asoc[identity,],ncol=ncol(asoc),dimnames=dimnames(asoc));
            group<-rep(groupid,length(time1));
            task<-rep(taskid,length(time1));
            coxdata<-data.frame(id=as.factor(id),time1,time2,status,identity,stMetric,indVar,group,task);
            #And for the ML models
            group<-rep(groupid,length(orderAcq));
            task<-rep(taskid,length(orderAcq));
            mldata<-data.frame(group,task,orderAcq,learnMetric,totalMetric,asoc[orderAcq,]);
            
            
            
            
            
            
            
            callNextMethod(.Object,idname=idname,assMatrix=assMatrix,asoc=asoc,orderAcq=orderAcq,ties=ties,trueTies=trueTies,groupid=groupid,taskid=taskid,coxdata=coxdata,mldata=mldata)
            
          }
)

oaData<-function(idname=NULL, assMatrix,asoc,orderAcq,ties=NULL,trueTies=list(NULL),groupid="1",taskid="1"){
  
  
  new("oaData",idname=idname,assMatrix=assMatrix,asoc=asoc,orderAcq=orderAcq,ties=ties,trueTies=trueTies,groupid=groupid,taskid=taskid);
  
}

#Gets likelihood for a multiplicative model using a transformed cox model with specified value of s
multiCoxLikelihood<-function(s,oadata,formula=NULL,bounded=FALSE){
  #Fit transformed cox model
  
  #Calulate the appropriate value of s for each line of data
  if(is.null(oadata@coxdata$sParamIndex)){
    
    sVect<-rep(s,dim(oadata@coxdata)[1])
    
  }else{
    
    stemp<-1*(oadata@coxdata$sParamIndex>0)
    stemp[which(oadata@coxdata$sParamIndex>0)]<-s[oadata@coxdata$sParamIndex[oadata@coxdata$sParamIndex>0]]
    sVect<-stemp;
    
  }
  
  if(bounded==T){
    model<-coxph(Surv(time1,time2,status)~strata(group,task)+offset(log(sVect*stMetric+(1-sVect))),data=oadata@coxdata);
  }else	
  {
    model<-coxph(Surv(time1,time2,status)~strata(group,task)+offset(log(sVect*stMetric+1)),data=oadata@coxdata);
  }
  if(is.null(formula)){
    model1<-model
    return(-model1$loglik[1])
  }else{
    model1<-update(model,formula);
  }
  -model1$loglik[2];
  
}



#Combine the coxdata from two or more oaData objects, into an oaData object with NA in all other slots
combineOaCoxData<-function(oaNames,sParam=NULL){
  
  #Set the default sParam (same s for all diffusions) if it is not defined
  if(is.null(sParam)) sParam<-rep(1,length(oaNames));
  
  newOaObject<-eval(as.name(oaNames[1]));
  sParamIndex<-rep(sParam[1],dim(newOaObject@coxdata)[1])
  
  for(i in 2:length(oaNames)){
    
    
    newOaObject@coxdata<-rbind(newOaObject@coxdata,eval(as.name(oaNames[i]))@coxdata);
    sParamIndex<-c(sParamIndex,rep(sParam[i],dim(eval(as.name(oaNames[i]))@coxdata)[1]));
    
  }
  newOaObject@idname<-NA;
  newOaObject@assMatrix<-matrix(NA);
  newOaObject@asoc<-matrix(NA);
  newOaObject@orderAcq<-NA;
  newOaObject@groupid<-"NA";
  newOaObject@taskid<-"NA";
  newOaObject@mldata<-data.frame(NA);
  
  newOaObject@coxdata<-cbind(newOaObject@coxdata,sParamIndex)
  
  return(newOaObject);
}


#Define class of object for the fitted multiplicative model
setClass("multiCoxFit",representation(optimisation="list",sParam="numeric",bounded="logical",formula="formula",coxcall="call",coxdf="numeric",coef="matrix",conf.int="matrix",loglik="numeric",aic="numeric", aicc="numeric",nulllik="numeric",aicNull="numeric",aiccNull="numeric",LRTsocTransTS="numeric",LRTsocTransPV="numeric",fitted="numeric",time="numeric",status="numeric",task="factor",group="factor"));

#Method for initializing object- including model fitting
setMethod("initialize",
          signature(.Object = "multiCoxFit"),
          function (.Object, oadata,sParam=NULL,formula,startValue=startValue,bounded=bounded,interval=interval,method=method,...) 
          {
            
            #If the oadata if a character vector containing multiple oaData objects, combine into a single object first
            if(is.character(oadata)){		
              #Set the default sParam (same s for all diffusions) if it is not defined
              if(is.null(sParam)) sParam<-rep(1,length(oadata));
              oadata<-combineOaCoxData(oadata,sParam=sParam);
            }else{
              
              sParam<-1;	
            }
            
            sampleSize<-sum(oadata@coxdata$status);
            noSParam<-length(levels(as.factor(sParam[sParam>0])))
            
            #Set staring values if not specified by the user
            if(is.null(startValue)) 		startValue<-rep(0, noSParam);
            
            
            if(length(startValue)>1) method="nlminb";
            
            #Optimise for s
            #If bounded==T the social learning parameter is constrained between 0 and 1, otherwise between 0 and Inf
            if(method=="nlminb"){	
              if(bounded==T){
                fit1<-nlminb(start=startValue,objective=multiCoxLikelihood,lower=0,upper=1,oadata=oadata,formula=formula,bounded=T);
              }else{			
                fit1<-nlminb(start=startValue,objective=multiCoxLikelihood,lower=0,upper=Inf,oadata=oadata,formula=formula,bounded=F);
              }
              stFitted<-fit1$par;
            }
            if(method=="optimise"){	
              if(bounded==T){
                fit1<-optimise(f=multiCoxLikelihood,interval=c(0,1),oadata=oadata,formula=formula,bounded=T);
              }else{			
                fit1<-optimise(f=multiCoxLikelihood,interval=interval,oadata=oadata,formula=formula,bounded=F);
              }
              stFitted<-fit1$minimum;
            }
            
            
            #Calulate the appropriate value of s for each line of data
            if(is.null(oadata@coxdata$sParamIndex)){
              
              sVect<-rep(stFitted,dim(oadata@coxdata)[1])
              
            }else{
              
              stemp<-1*(oadata@coxdata$sParamIndex>0)
              stemp[which(oadata@coxdata$sParamIndex>0)]<-stFitted[oadata@coxdata$sParamIndex[oadata@coxdata$sParamIndex>0]]
              sVect<-stemp;
              
            }
            
            
            
            #Feed back into the cox model and record appropriate output
            if(bounded==F) model<-coxph(Surv(time1,time2,status)~strata(group,task)+offset(log(sVect*stMetric+1)),data=oadata@coxdata);
            if(bounded==T) model<-coxph(Surv(time1,time2,status)~strata(group,task)+offset(log(sVect*stMetric+(1-sVect))),data=oadata@coxdata);
            if(is.null(formula)){
              model1<-model
            }else{
              model1<-update(model,formula);
            }
            
            if(is.null(formula)){
              
              coef<-matrix(nrow=1,ncol=1);
              conf.int<-matrix(nrow=1,ncol=1);
            }else{
              coef<-summary(model1)$coef;
              conf.int<-summary(model1)$conf.int
            }
            
            #For a frailty model the summary does not output an object, the stats are printed to the screen, therefore we cannot store them in the object
            if(is.null(coef)) {
              coef<-matrix(nrow=1,ncol=1);
            }
            if(is.null(conf.int)) conf.int<-matrix(nrow=1,ncol=1);
            
            loglik<-model1$loglik;
            coxcall<-model1$call;
            nullmodel<-coxph(Surv(time1,time2,status)~strata(group,task),data=oadata@coxdata);
            if(is.null(formula)){
              nullmodel<-nullmodel
              nulllik<-nullmodel$loglik[1];		
            }else{
              nullmodel<-update(nullmodel,formula);
              nulllik<-nullmodel$loglik[2];
            }
            
            #Store appropriate formula if none is given
            if(is.null(formula)) {
              formula<-~.
              loglik<-rep(loglik,2)
            }
            
            #Record the total df for the cox model
            coxdf<-round(sum(model1$df));
            
            LRTsocTransTS<-2*(-nulllik+loglik[2])
            LRTsocTransPV<-1-pchisq(LRTsocTransTS,df= noSParam);
            
            
            #Record aic of fitted and null model. Also small sample AIC
            if(is.na(coef[1])){
              aic<-2*(length(stFitted)+coxdf)-2*loglik[2];
              aicc<-2*(length(stFitted)+coxdf)*(sampleSize/(sampleSize-length(stFitted)-coxdf-1))-2*loglik[2];
              aicNull<-2*coxdf-2*nulllik;
              aiccNull<-2*(coxdf)*(sampleSize/(sampleSize-coxdf-1))-2*nulllik;
            }else{
              aic<-2*(dim(coef)[1]+length(stFitted))-2*loglik[2];
              aicc<-2*(dim(coef)[1]+length(stFitted))*(sampleSize/(sampleSize-(dim(coef)[1]+length(stFitted))-1))-2*loglik[2];
              aicNull<-2*(dim(coef)[1])-2*nulllik;
              aiccNull<-2*(dim(coef)[1])*(sampleSize/(sampleSize-(dim(coef)[1])-1))-2*nulllik;
            }
            
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aic)|is.nan(aicc)){}else{
              if(aicc<aic) aicc<-NaN;
            }
            if(is.nan(aicNull)|is.nan(aiccNull)){}else{
              if(aiccNull<aicNull) aiccNull<-NaN;
            }
            
            #Record fitted values
            fitted<-exp(predict(model1));
            time<-oadata@coxdata$time2;
            status<-oadata@coxdata$status;
            task<-oadata@coxdata$task;
            group<-oadata@coxdata$group;       
            
            
            
            callNextMethod(.Object,bounded=bounded,optimisation=fit1,sParam=sParam,formula=formula,coxcall=coxcall,coxdf=coxdf,coef=coef,conf.int=conf.int,loglik=loglik,aic=aic,aicc=aicc,nulllik=nulllik,aicNull=aicNull,aiccNull=aiccNull,LRTsocTransTS=LRTsocTransTS,LRTsocTransPV=LRTsocTransPV,fitted=fitted,time=time,status=status,task=task,group=group,...)
            
            
            
          }
)

#Function for implementing the initialization
multiCoxFit<-function(oadata,sParam=NULL,formula=NULL,startValue=NULL,bounded=FALSE,interval=c(0,999),method="optimise"){
  
  new("multiCoxFit",oadata=oadata, sParam=sParam, formula=formula,startValue=startValue,bounded=bounded,interval=interval,method=method)	
}


setMethod("summary",
          signature(object = "multiCoxFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])));
            cat("Summary of Multiplicative Social Transmission Model\nOrder of acquisition data\n");
            
            
            if(is.na(object@coef[1])){
              
              if(object@bounded==T) {
                cat("Bounded parameterisation.");
                sumtable<-data.frame(Estimate=object@optimisation$minimum,Unbounded=object@optimisation$minimum/(1-object@optimisation$minimum), row.names="Social Transmission");
              }
              if(object@bounded==F) {
                cat("Unbounded parameterisation.");
                sumtable<-data.frame(Estimate=object@optimisation$minimum,Bounded=object@optimisation$minimum/(1+object@optimisation$minimum), row.names="Social Transmission");
              }
              
              
            }else{
              
              if(is.null(object@optimisation$par)){
                
                if(object@bounded==T) {
                  cat("Bounded parameterisation.");
                  sumtable<-data.frame(Estimate=c(object@optimisation$minimum,object@coef[,1]),Unbounded=c(object@optimisation$minimum/(1-object@optimisation$minimum),rep(NA,length(object@coef[,1]))),se=c(rep(NA,noSParam),object@coef[,3]),z=c(rep(NA,noSParam),object@coef[,4]),p=c(rep(NA,noSParam),object@coef[,5]), row.names=c(paste("Social transmission",1:noSParam),row.names(object@coef)));			}
                if(object@bounded==F) {
                  cat("Unbounded parameterisation.");
                  sumtable<-data.frame(Estimate=c(object@optimisation$minimum,object@coef[,1]),Bounded=c(object@optimisation$minimum/(1+object@optimisation$minimum),rep(NA,length(object@coef[,1]))),se=c(rep(NA,noSParam),object@coef[,3]),z=c(rep(NA,noSParam),object@coef[,4]),p=c(rep(NA,noSParam),object@coef[,5]), row.names=c(paste("Social transmission",1:noSParam),row.names(object@coef)));
                  
                }
              }else{
                
                if(object@bounded==T) {
                  cat("Bounded parameterisation.");
                  sumtable<-data.frame(Estimate=c(object@optimisation$par,object@coef[,1]),Unbounded=c(object@optimisation$par/(1-object@optimisation$par),rep(NA,length(object@coef[,1]))),se=c(rep(NA,noSParam),object@coef[,3]),z=c(rep(NA,noSParam),object@coef[,4]),p=c(rep(NA,noSParam),object@coef[,5]), row.names=c(paste("Social transmission",1:noSParam),row.names(object@coef)));			}
                if(object@bounded==F) {
                  cat("Unbounded parameterisation.");
                  sumtable<-data.frame(Estimate=c(object@optimisation$par,object@coef[,1]),Bounded=c(object@optimisation$par/(1+object@optimisation$par),rep(NA,length(object@coef[,1]))),se=c(rep(NA,noSParam),object@coef[,3]),z=c(rep(NA,noSParam),object@coef[,4]),p=c(rep(NA,noSParam),object@coef[,5]), row.names=c(paste("Social transmission",1:noSParam),row.names(object@coef)));
                }
              }
            }	
            cat("\n\nCoefficients:\n");
            print(sumtable);
            cat("\n\n")
            anova(object);
            
          }
)


setMethod("anova",
          signature(object = "multiCoxFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            if(is.na(object@coef[1])){
              
              table<-data.frame(Df=c(noSParam+object@coxdf,0+object@coxdf),LogLik=c(-object@loglik[2],-object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));			
              
            }else{
              
              table<-data.frame(Df=c((dim(object@coef)[1]+ noSParam),(dim(object@coef)[1])),LogLik=c(-object@loglik[2],-object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
            }
            
            atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine multiplicatively\n", class = c("anova", "data.frame"));
            
            
            return(atable);
            
          }
)


nbdaProfile<-function(data=NULL,model=NULL,range=NULL,range2=NULL, additive=TRUE, bounded=F,resolution=1000,confInt=c(0.95,0.99),ylim=NULL,param=1,otherParam=NULL,progress=F, stepLength=NULL){
  
  if(class(model)=="multiCoxFit"){
    additive<-FALSE
    formula<-model@formula;
    #If the oadata if a character vector containing multiple oaData objects, combine into a single object first
    if(is.character(data)){		
      if(additive==FALSE) data<-combineOaCoxData(data,sParam=model@sParam)
    }
  }
  
  if(class(model)=="tadaFit"){
    additive<-model@additive
    bounded<-F
    #If the oadata if a character vector containing multiple oaData objects, combine into a single object first
    if(is.character(data)){		
      data<-combineTadaData(data,sParam=model@sParam)
    }
  }	
  
  if(class(model)=="discreteTadaFit"){
    additive<-model@additive;
    bounded=F;
    
  }
  
  if(is.null(model)){
    fitted.value<-NULL;
  }else{
    if(is.null(model@optimisation$minimum)){ 
      fitted.value<-model@optimisation$par;
    }else{
      fitted.value<-model@optimisation$minimum;
    }
    
    #Convert parameters in model to bounded/unbounded as specified
    if(bounded==T){
      if(model@bounded==F) {
        if(class(model)=="multiCoxFit"){ 
          fitted.value<-fitted.value/(1+fitted.value);
        }else{
          if(is.null(model@asocialVar)){
            fitted.value<-fitted.value/(1+fitted.value);
          }else{
            fitted.value[1:(length(fitted.value)-length(model@asocialVar))]<-fitted.value[1:(length(fitted.value)-length(model@asocialVar))]/(1+fitted.value[1:(length(fitted.value)-length(model@asocialVar))]);	
          }	
        }
      }
    }
    if(bounded==F){
      if(model@bounded==T) {
        if(class(model)=="multiCoxFit"){ 
          fitted.value<-fitted.value/(1-fitted.value)
        }else{
          if(is.null(model@asocialVar)){
            fitted.value<-fitted.value/(1-fitted.value)
          }else{
            fitted.value[1:(length(fitted.value)-length(model@asocialVar))]<-fitted.value[1:(length(fitted.value)-length(model@asocialVar))]/(1+fitted.value[1:(length(fitted.value)-length(model@asocialVar))]);	
          }	
        }
      }
    }
  }
  
  
  #If values are specified for the other parameters, replace the fitted values with these
  if(is.null(otherParam)){}else{
    
    fitted.value[-param]<-otherParam
    
  }	
  
  #Set ranges of plot to 0,1 if bounded, user specified if unbounded with  default of 0,100, unless defined by user
  
  if(is.null(range)){
    if(bounded==T){
      range<-c(0,1);
      if(fitted.value[param[1]]>1) range<-c(0,2*fitted.value[param[1]]);
    }else{
      range<-c(0,2*fitted.value[param[1]]);
    }
  }
  
  if(is.null(range2)){
    if(bounded==T){
      range2<-c(0,1);
      if (is.na(fitted.value[param[2]])==FALSE){
        if(fitted.value[param[2]]>1) range2<-c(0,2*fitted.value[param[2]]);
      }
    }else{
      range2<-c(0,100);		
      if (is.na(fitted.value[param[2]])==FALSE){
        if(fitted.value[param[2]]>0) range2<-c(0,2*fitted.value[param[2]]);
        if(fitted.value[param[2]]<0) range2<-c(2*fitted.value[param[2]],0);
      }
    }
  }
  
  #If the user specifies only one parameter to be plotted, plot in 1D with other parameters set at fitted values, unless the user specifes values for the other parameters
  
  if(length(param)==1){
    
    #Create vector of x values
    xValues<-seq(range[1],range[2],length.out=resolution);
    #Cut out the last value (likelihood is 0) if bounded
    if(bounded==T) xValues<-xValues[-resolution];
    #Initialise y values
    yValues<-vector(mode="numeric",length=length(xValues));
    
    #If there is more than one parameter set the others to their fitted value
    if(is.null(model@optimisation$minimum)){
      
      newXValues<-t(matrix(rep(fitted.value,length(xValues)),nrow=length(fitted.value)));
      newXValues[,param]<-xValues;
      
    }else{
      
      newXValues<-as.matrix(xValues);
      
      
    }
    
    if(class(model)=="multiCoxFit") {
      asocialVar<-NULL;
    }else{
      asocialVar<-model@asocialVar;
      if(asocialVar=="NaN") asocialVar<-NULL;
    }
    
    #Calculate likelihood for each x value	
    for(i in 1:dim(newXValues)[1]){
      if(class(model)=="multiCoxFit") {
        yValues[i]<-multiCoxLikelihood(newXValues[i,],oadata=data,formula=model@formula,bounded=bounded);
        labels<-paste("Social transmission parameter", param);
      }else{
        
        if(class(model)=="tadaFit"){
          yValues[i]<-tadaLikelihood(l=newXValues[i,],tadata=data,asocialVar=asocialVar,bounded=bounded,additive=additive, task=model@task, group=model@group,sParam=model@sParam)
          labels<-model@varNames[param];
          
        }else{
          
          if(class(model)=="discreteTadaFit"){
            
            yValues[i]<-discreteTadaLikelihood(l=newXValues[i,],tadata=data,asocialVar=asocialVar,bounded=bounded,additive=additive, task=model@task, group=model@group, stepLength=stepLength,sParam=model@sParam);
            labels<-model@varNames[param];
            
            
            
          }else{
            
            yValues[i]<-addLikelihood(newXValues[i,],data=data,asocialVar=asocialVar,bounded=bounded);
            labels<-model@varNames[param];
          }
        }
      }
    }
    
    
    #Plot
    plot(xValues,yValues,type="l",xlab=labels[1],ylab="log Likelihood",ylim=ylim);
    
    
    if(is.null(fitted.value)){}else{
      
      if(class(model)=="multiCoxFit") {
        confIntValue<-(multiCoxLikelihood(fitted.value,oadata=data,formula=formula,bounded=bounded)+qchisq(confInt,1)/2);			}else{
          
          if(class(model)=="tadaFit"){
            confIntValue<-(tadaLikelihood(fitted.value,tadata=data,asocialVar=asocialVar,bounded=bounded, task=model@task, group=model@group, sParam=model@sParam)+qchisq(confInt,1)/2);
          }else{
            
            if(class(model)=="discreteTadaFit"){
              
              confIntValue<-(discreteTadaLikelihood(fitted.value,tadata=data,asocialVar=asocialVar,bounded=bounded, task=model@task, group=model@group, sParam=model@sParam, stepLength=stepLength))+(qchisq(confInt,1)/2);
              
            }else{
              
              confIntValue<-(addLikelihood(fitted.value,data=data,asocialVar=asocialVar,bounded=bounded, sParam=model@sParam)+qchisq(confInt,1)/2);
            }
          }
        }
      
      abline(h=confIntValue,lty=1+1:length(confInt));
      
      #Find confidence intervals. Replace with NA if on boundary (unless boundary is zero for ST paramter)
      STlCI<-vector();
      STuCI<-vector();
      for (i in 1:length(confInt)){
        STlCI<-c(STlCI,min(xValues[yValues<(confIntValue[i])]));
        STuCI<-c(STuCI,max(xValues[yValues<(confIntValue[i])]));
        if(STlCI[i]==0) {}else{
          if(STlCI[i]==min(xValues)) STlCI[i]<-NA;
        }
        if(STuCI[i]==max(xValues)) STuCI[i]<-NA;
        
      }
    }
    ciTable<-cbind(confInt,STlCI,STuCI);
    if(class(model)=="multiCoxFit") {
      dimnames(ciTable)[2]<-list(c("Confidence level",paste("Social Transmission",paste(rep(param,each=2),c("lower","upper")))));
    }else{
      dimnames(ciTable)[2]<-list(c("Confidence level",paste(rep(model@varNames[param],each=2),c("lower","upper"))));
    }
    return(ciTable)
    
    
  }else{
    #If the user specifies two parameters to be plotted, plot contours in 2D with other parameters set at fitted values, unless the user specifes values for the other parameters		
    
    #Create vectors of x and y values
    xValues<-seq(range[1],range[2],length.out=resolution);
    #Cut out the last value (likelihood is 0) if bounded
    if(bounded==T) xValues<-xValues[-resolution];
    yValues<-seq(range2[1],range2[2],length.out=resolution);
    #Cut out the last value (likelihood is 0) if bounded
    if(bounded==T) yValues<-yValues[-resolution];
    
    
    x<-matrix(nrow=length(xValues),ncol=length(yValues));
    y<-matrix(nrow=length(xValues),ncol=length(yValues));
    z<-matrix(nrow=length(xValues),ncol=length(yValues));
    
    
    asocialVar<-model@asocialVar;
    if(asocialVar=="NaN") asocialVar<-NULL;	
    
    #Generate grid of x y and z values
    for(i in 1:length(xValues)){
      
      for(j in 1:length(yValues)){
        
        if(progress==T){
          cat(i," ",j,"\n");
          flush.console();
        }
        
        x[i,j]<-xValues[i];
        y[i,j]<-yValues[j];
        
        input<-fitted.value;
        input[param[1]]<-xValues[i];
        input[param[2]]<-yValues[j];
        
        if(class(model)=="multiCoxFit") {
          z[i,j]<-multiCoxLikelihood(s=input,oadata=data,formula=model@formula,bounded=bounded);					labels<-paste("Social transmission parameter", param);
        }else{
          
          if(class(model)=="tadaFit"){
            z[i,j]<-tadaLikelihood(l=input,tadata=data,asocialVar=asocialVar,bounded=bounded,additive=additive,task=model@task, group=model@group,sParam=model@sParam);
            labels<-model@varNames[param];					}else{
              
              if(class(model)=="discreteTadaFit"){
                
                z[i,j]<-discreteTadaLikelihood(l=input,tadata=data,asocialVar=asocialVar,bounded=bounded,additive=additive,task=model@task, group=model@group,sParam=model@sParam,stepLength=stepLength);							labels<-model@varNames[param];					
                
              }else{
                
                z[i,j]<-addLikelihood(l=input,data=data,asocialVar=asocialVar,bounded=bounded,sParam=model@sParam);
                labels<-model@varNames[param];
              }
            }
        }
      }
    }
    
    if(is.null(model)){
      
      filled.contour(xValues,yValues,z,xlab=labels[1],ylab=labels[2]);
      
    }else{
      confIntValue<-model@optimisation$objective+c(0,(qchisq(confInt,1)/2));
      filled.contour(xValues,yValues,z,levels=confIntValue,xlab=labels[1],ylab=labels[2]);
      
      #Find confidence intervals. Replace with NA if on boundary (unless boundary is zero for ST paramter)
      STlCI<-vector();
      STuCI<-vector();
      ALlCI<-vector();
      ALuCI<-vector();
      for (i in 1:length(confInt)){
        STlCI<-c(STlCI,min(x[z<(confIntValue[i+1])]));
        STuCI<-c(STuCI,max(x[z<(confIntValue[i+1])]));
        if(STlCI[i]==0) {}else{
          if(STlCI[i]==min(x)) STlCI[i]<-NA;
        }
        if(STuCI[i]==max(x)) STuCI[i]<-NA;
        ALlCI<-c(ALlCI,min(y[z<(confIntValue[i+1])]));
        ALuCI<-c(ALuCI,max(y[z<(confIntValue[i+1])]));
        if(ALlCI[i]==min(y)) ALlCI[i]<-NA;
        if(ALuCI[i]==max(y)) ALuCI[i]<-NA;
      }
      ciTable<-cbind(confInt,STlCI,STuCI,ALlCI,ALuCI);
      if(class(model)=="multiCoxFit") {
        dimnames(ciTable)[2]<-list(c("Confidence level",paste("Social Transmission",paste(rep(param,each=2),c("lower","upper")))));
      }else{
        dimnames(ciTable)[2]<-list(c("Confidence level",paste(rep(model@varNames[param],each=2),c("lower","upper"))));
      }
      return(ciTable)
      
      
    }
  }
}


#Define class of object for the fitted nbda model
setClass("tadaFit",representation(optimisation="list",additive="logical",sParam="numeric",bounded="logical",loglik="numeric",aic="numeric", aicc="numeric",nulllik="numeric",aicNull="numeric",aiccNull="numeric",LRTsocTransTS="numeric",LRTsocTransPV="numeric",asocialVar="numeric",varNames="character",task="logical",group="logical")); 


#Method for initializing object- including model fitting
setMethod("initialize",
          signature(.Object = "tadaFit"),
          function (.Object, tadata,additive,sParam,startValue,asocialVar,bounded,task,group,...) 
          {
            
            #If the tadata if a character vector containing multiple oaData objects, combine into a single object first
            #Also set the default sParam if nothing is specified
            if(is.character(tadata)){		
              if(is.null(sParam)) sParam<-rep(1,length(tadata));			tadata<-combineTadaData(tadata,sParam=sParam);
              
            }else{
              sParam<-1;				
            }
            
            #Calculate the number of social transmission parameters to be fitted
            noSParam<-length(levels(as.factor(sParam[sParam>0])))
            
            sampleSize<-sum(tadata@nbdaData$status);
            
            extraParam<-0;
            
            if(task) extraParam<-extraParam+length(levels(tadata@nbdaData$task))-1;
            if(group) extraParam<-extraParam+length(levels(tadata@nbdaData$group))-1;
            
            #Set staring values if not specified by the user
            if(is.null(startValue)){
              startValue<-c(rep(0,noSParam),1,rep(0,length(asocialVar)+extraParam));
              startInd<-1;
            }else{
              
              startInd<-0;
              
            }
            
            if((length(asocialVar)+extraParam)==0){
              null<-optimise(nullTadaLikelihood,interval=c(0,max(tadata@nbdaData$time2)+1000),asocialVar=asocialVar,tadata=tadata, additive=additive,sParam=sParam);
              if(startInd==1) startValue[2]<-null$minimum;
            }else{
              null<-nlminb(start=startValue[-(1: noSParam)],nullTadaLikelihood,lower=c(0,rep(-Inf,length(asocialVar)+extraParam)),asocialVar=asocialVar,tadata=tadata,task=task,group=group,additive=additive,sParam=sParam);
              if(startInd==1) startValue[-(1: noSParam)]<-null$par;
            }
            
            
            
            fit1<-nlminb(start=startValue,tadaLikelihood,lower=c(rep(0,noSParam+1),rep(-Inf,length(asocialVar)+extraParam)),sParam=sParam,asocialVar=asocialVar,tadata=tadata,bounded=bounded,task=task,group=group,additive=additive);
            
            #Perform LRT for social transmission
            loglik<-fit1$objective;
            nulllik<-null$objective;
            LRTsocTransTS<-2*(nulllik-loglik)
            LRTsocTransPV<-1-pchisq(LRTsocTransTS,df=noSParam);
            
            if(is.null(asocialVar)) asocialVar<-NaN;
            
            
            #Calculate aic and for model without social transmission
            k<-length(fit1$par);
            aic<-(2*k)+(2*loglik);
            aicc<-aic+(2*k*(k+1)/(sampleSize-k-1));
            
            k<-length(null$par);
            aicNull<-(2*k)+(2*nulllik);
            aiccNull<-aicNull+(2*k*(k+1)/(sampleSize-k-1));
            
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aic)|is.nan(aicc)){}else{
              if(aicc<aic) aicc<-NaN;
            }
            if(is.nan(aicNull)|is.nan(aiccNull)){}else{
              if(aiccNull<aicNull) aiccNull<-NaN;
            }
            
            #Extract names of variables
            if(is.nan(asocialVar[1])){
              
              varNames<-c(paste("Social Transmission",1:noSParam),"1/Rate scaling");			
              
            }else{
              varNames<-c(paste("Social Transmission",1:noSParam),"1/Rate scaling",unlist(dimnames(tadata@nbdaData)[2])[asocialVar+6]);
            }
            if(task) varNames<-c(varNames,paste("Task",levels(tadata@nbdaData$task)[-1]));
            if(group) varNames<-c(varNames,paste("Group",levels(tadata@nbdaData$group)[-1]));
            
            
            
            
            
            callNextMethod(.Object,optimisation=fit1,additive=additive,sParam=sParam,bounded=bounded,loglik=loglik, aic=aic,aicc=aicc,nulllik=nulllik,aicNull=aicNull,aiccNull=aiccNull,LRTsocTransTS=LRTsocTransTS,LRTsocTransPV=LRTsocTransPV,asocialVar=asocialVar,varNames=varNames,task=task,group=group,...) 
            
            
            
            
          }
)


#Function for implementing the initialization
tadaFit<-function(tadata,sParam=NULL,asocialVar=NULL,startValue=NULL,bounded=FALSE,task=FALSE,group=FALSE,additive=TRUE){
  
  
  if(bounded==T) {
    
    cat("Bounded parameterisation not yet available for time of acquisition models");
    return("Bounded parameterisation not yet available for time of acquisition models");	
    
  }	
  
  new("tadaFit",tadata=tadata,sParam=sParam,asocialVar=asocialVar,startValue=startValue,bounded=bounded,task=task,group=group,additive=additive)	
}



#Gets likelihood for a multiplicative nbda model 
tadaLikelihood<-function(l,tadata,sParam=NULL,bounded=FALSE,asocialVar=NULL,task=FALSE,group=FALSE,additive=TRUE){
  
  
  
  
  if(is.null(sParam)){
    noSParam<-1
    
  }else{
    #Calculate the number of social transmission parameters to be fitted
    noSParam<-length(levels(as.factor(sParam[sParam>0])))
  }
  
  #Found that paramterising by 1/rate works better for optimisation
  l[noSParam+1]<-1/l[noSParam+1];
  
  if(is.null(asocialVar)){
  }else{
    #		if(max(asocialVar)>dim(tadata@oadata@asoc)[2]){
    
    #			return("Invalid asocial learning variables selected")
    
    
    #		}
  }
  
  #Cut down to the asocial variables specified by asocialVar
  asoc<-as.matrix(tadata@nbdaData[,asocialVar+6]);
  
  if(task){
    
    taskNames<-paste("task",levels(tadata@nbdaData$task),sep="");
    taskMatrix<-matrix(nrow=dim(tadata@nbdaData)[1],ncol=length(taskNames)-1,dimnames=list(NULL,taskNames[-1]))
    for(i in 2:length(taskNames)){
      
      taskMatrix[,i-1]<-(tadata@nbdaData$task==levels(tadata@nbdaData$task)[i])*1
      
    }
    asoc<-cbind(asoc,taskMatrix);
    
  }
  
  if(group){
    
    groupNames<-paste("group",levels(tadata@nbdaData$group),sep="");
    groupMatrix<-matrix(nrow=dim(tadata@nbdaData)[1],ncol=length(groupNames)-1,dimnames=list(NULL,groupNames[-1]))
    for(i in 2:length(groupNames)){
      
      groupMatrix[,i-1]<-(tadata@nbdaData$group==levels(tadata@nbdaData$group)[i])*1
      
    }
    asoc<-cbind(asoc,groupMatrix);		
  }
  
  
  #Calulate the appropriate value of s for each line of data
  if(is.null(tadata@nbdaData$sParamIndex)){
    
    sVect<-rep(l[1],dim(tadata@nbdaData)[1])
    
  }else{
    
    stemp<-1*(tadata@nbdaData$sParamIndex>0)
    stemp[which(tadata@nbdaData$sParamIndex>0)]<-l[tadata@nbdaData$sParamIndex[tadata@nbdaData$sParamIndex>0]]
    sVect<-stemp;
    
  }
  
  
  #Check if there are any asocial variables, set rate to one if not
  if(dim(asoc)[2]==0){
    lpRateAsoc<-rep(0,dim(tadata@nbdaData)[1]);
  }else{
    lpRateAsoc<-apply(t(matrix(rep(l[-(1:(noSParam+1))],dim(tadata@nbdaData)[1]),nrow=dim(asoc)[2]))*asoc,1,sum);
  }
  
  lpSocTrans<-sVect*tadata@nbdaData[6];
  
  if(bounded==FALSE){
    
    if(additive==TRUE){
      
      
      likelihood<-tadata@nbdaData$status*(log(l[noSParam+1])+log(lpSocTrans+exp(lpRateAsoc))-(l[noSParam+1]*(lpSocTrans+exp(lpRateAsoc)))*tadata@nbdaData$time)+(1-tadata@nbdaData$status)*((-l[noSParam+1]*(lpSocTrans+exp(lpRateAsoc)))*tadata@nbdaData$time);
      
    }else{
      
      likelihood<-tadata@nbdaData$status*(log(l[noSParam+1])+log(lpSocTrans+1)+lpRateAsoc-l[noSParam+1]*(lpSocTrans+1)*exp(lpRateAsoc)*tadata@nbdaData$time)+(1-tadata@nbdaData$status)*(-l[noSParam+1]*(lpSocTrans+1)*exp(lpRateAsoc)*tadata@nbdaData$time);
      
    }
    
    
  }
  
  if(bounded==TRUE){
    
    return("Bounded parameterisation not yet available for time of acquisition models");
    
  }
  
  
  return(-sum(likelihood));
  
}

nullTadaLikelihood<-function(l,tadata,sParam=NULL,asocialVar=NULL,task=FALSE,group=FALSE,additive=TRUE){
  
  if(is.null(sParam)){
    noSParam<-1
    
  }else{
    #Calculate the number of social transmission parameters to be fitted
    noSParam<-length(levels(as.factor(sParam[sParam>0])))
  }
  
  l<-c(rep(0,noSParam),l);
  tadaLikelihood(l,tadata,bounded=FALSE,asocialVar=asocialVar,task=task,group=group,additive=additive,sParam=sParam);
  
}

#Combine the coxdata from two or more oaData objects, into an oaData object with NA in all other slots
combineTadaData<-function(taNames,sParam=NULL){
  
  
  #Set the default sParam if nothing is specified
  if(is.null(sParam)) sParam<-rep(1,length(taNames));
  
  newTaObject<-eval(as.name(taNames[1]));
  sParamIndex<-rep(sParam[1],dim(newTaObject@nbdaData)[1])
  
  for(i in 2:length(taNames)){
    
    
    newTaObject@nbdaData<-rbind(newTaObject@nbdaData,eval(as.name(taNames[i]))@nbdaData);
    sParamIndex<-c(sParamIndex,rep(sParam[i],dim(eval(as.name(taNames[i]))@nbdaData)[1]));		
  }
  newTaObject@oadata@idname<-NA;
  newTaObject@oadata@assMatrix<-matrix(NA);
  newTaObject@oadata@asoc<-matrix(NA);
  newTaObject@oadata@orderAcq<-NA;
  newTaObject@oadata@groupid<-"NA";
  newTaObject@oadata@taskid<-"NA";
  newTaObject@oadata@mldata<-data.frame(NA);
  newTaObject@oadata@coxdata<-data.frame(NA);
  
  newTaObject@nbdaData<-cbind(newTaObject@nbdaData,sParamIndex)
  
  return(newTaObject);
}


setMethod("anova",
          signature(object = "tadaFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            if(is.null(object@optimisation$par[1])){
              
              table<-data.frame(Df=c(1,0),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
              
            }else{   
              
              table<-data.frame(Df=c(length(object@optimisation$par),length(object@optimisation$par)-noSParam),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
            }
            
            if(object@additive){atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine additively\n", class = c("anova", "data.frame"))}else{
              atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine multiplicatively\n", class = c("anova", "data.frame"))			
            }
            
            return(atable);
            
          }
)


setMethod("summary",
          signature(object = "tadaFit"),
          function (object, ...) 
          {
            
            if(object@additive==T)cat("Summary of Additive Social Transmission Model\nTime of Acquisition Data\n");
            if(object@additive==F)cat("Summary of Multiplicative Social Transmission Model\nTime of Acquisition Data\n");
            
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            # 		if(length(object@optimisation$par)<3){
            
            
            if(object@bounded==F) {
              cat("Unbounded parameterisation.");
              sumtable<-data.frame(Estimate=c(object@optimisation$par),Bounded=c(object@optimisation$par[1: noSParam]/(1+object@optimisation$par[1: noSParam]),rep("",length(object@optimisation$par)-noSParam)),row.names=object@varNames);   
              
            }		
            
            
            
            #   	}else{		
            
            #    	}
            
            
            cat("\n\nCoefficients\n");
            print(sumtable);
            cat("\n\n")
            anova(object);
            
          }
)


#Define class of object for the fitted additive model
setClass("addFit",representation(optimisation="list",sParam="numeric",bounded="logical",loglik="numeric",aic="numeric",aicc="numeric",nulllik="numeric",aicNull="numeric",aiccNull="numeric",LRTsocTransTS="numeric",LRTsocTransPV="numeric",asocialVar="numeric",varNames="character"));


#Method for initializing addFit object- including model fitting
setMethod("initialize",
          signature(.Object = "addFit"),
          function (.Object, oadata,sParam=NULL, startValue=startValue,bounded=bounded,method=method,interval=interval,asocialVar=NULL,...) 
          {
            
            if(class(oadata)=="oaData") {
              sParam<-1;
            }else{
              if(is.null(sParam)) sParam<-rep(1,length(oadata));
            }
            
            noSParam<-length(levels(as.factor(sParam[sParam>0])))
            
            #Set staring values if not specified by the user
            if(is.null(startValue)) 		startValue<-rep(0,length(asocialVar)+ noSParam);
            
            
            #Set vector of upper and lower values for each parameter. Unbounded for asocial learning variables
            lower<-c(rep(0,noSParam),rep(-Inf,length(asocialVar)));
            if (bounded==T)	upper<-c(rep(1,noSParam),rep(Inf,length(asocialVar)));
            if (bounded==F)	upper<-rep(Inf,length(asocialVar)+noSParam);
            
            if(is.null(interval)) interval<-c(0,999);
            
            #Optimise for s
            #If there are no asocial learning variables and a single st parameter we can use optimise. This is the default unless nlminb is specified. If there are asocial learning variables, nlminb is used
            if(length(startValue)==1){
              if(method=="nlminb"){	
                if(bounded==T){
                  fit1<-nlminb(start=startValue,objective=addLikelihood,lower=0,upper=1,data=oadata,bounded=T,asocialVar=NULL);
                  
                }else{			
                  fit1<-nlminb(start=startValue,objective=addLikelihood,lower=0,upper=Inf,data=oadata,bounded=F,asocialVar=NULL);
                }
              }else{
                if(bounded==T){
                  fit1<-optimise(f=addLikelihood,interval=c(0,1),data=oadata,bounded=T,asocialVar=NULL);
                }else{			
                  fit1<-optimise(f=addLikelihood,interval=interval,data=oadata,bounded=F,asocialVar=NULL);
                }					
              }
              nulllik<-nulladdLikelihood(0,oadata,bounded=bounded);
            }else{
              if(bounded==T){
                fit1<-nlminb(start=startValue,objective=addLikelihood,lower=lower,upper=upper,data=oadata,bounded=T,asocialVar=asocialVar,sParam=sParam);
                nullfit<-nlminb(start=startValue[-(1:noSParam)],objective=nulladdLikelihood,data=oadata,bounded=T,asocialVar=asocialVar);
              }else{			
                fit1<-nlminb(start=startValue,objective=addLikelihood,lower=lower,upper=upper,data=oadata,bounded=F,asocialVar=asocialVar,sParam=sParam);
                nullfit<-nlminb(start=startValue[-(1:noSParam)],objective=nulladdLikelihood,data=oadata,bounded=F,asocialVar=asocialVar);
                
              }
              nulllik<-nullfit$objective;
            }
            
            
            #Is the model being adjusted for true ties? CAUSES PROBLEMS FOR MULTIPLE DIFUSSIONS
            #		if (is.null(unlist(oadata@trueTies))) {
            #			trueTies<-FALSE;
            #		}else{
            #			trueTies<-TRUE;
            #		}
            
            
            
            #Perform LRT for social transmission
            loglik<-fit1$objective;
            LRTsocTransTS<-2*(nulllik-loglik)
            LRTsocTransPV<-1-pchisq(LRTsocTransTS,df=noSParam);
            
            if(is.null(asocialVar)) asocialVar<-NaN;
            
            sampleSize<-sampSizeExtract(oadata);
            
            #Calculate aic and for model without social transmission
            aic<-2*length(fit1$par)+2*loglik;
            aicc<-2*(length(fit1$par))*(sampleSize/(sampleSize-(length(fit1$par))-1))+2*loglik;
            
            if (is.nan(asocialVar[1])){
              aicNull<-2*nulllik;
              aiccNull<-aicNull;
              aic<-2+2*loglik;
              aicc<-2*(sampleSize/(sampleSize-2))+2*loglik;
            }else{
              aicNull<-2*length(asocialVar)+2*nulllik;
              aiccNull<-2*(length(asocialVar))*(sampleSize/(sampleSize-(length(asocialVar))-1))+2*nulllik;
            }	
            
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aic)|is.nan(aicc)){}else{
              if(aicc<aic) aicc<-NaN;
            }
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aicNull)|is.nan(aiccNull)){}else{
              if(aiccNull<aicNull) aiccNull<-NaN;
            }
            
            #Extract names of variables
            if(is.nan(asocialVar[1])){
              
              varNames<-c("Social Transmission");			
              
            }else{
              
              if(is.character(oadata)){
                
                varNames<-c(paste("Social transmission",1:noSParam),unlist(dimnames(eval(as.name(oadata[1]))@coxdata)[2])[asocialVar+6]);				
                
              }else{
                
                varNames<-c(paste("Social transmission",1:noSParam),unlist(dimnames(oadata@coxdata)[2])[asocialVar+6]);
              }
            }
            
            
            callNextMethod(.Object,bounded=bounded, optimisation=fit1,sParam=sParam,loglik=loglik, aic=aic,aicc=aicc,nulllik=nulllik,aicNull=aicNull,aiccNull=aiccNull,LRTsocTransTS=LRTsocTransTS,LRTsocTransPV=LRTsocTransPV,asocialVar=asocialVar,varNames=varNames,...)
            
            
            
          }
)		

#Function for implementing the initialization
addFit<-function(oadata,sParam=NULL,startValue=NULL,bounded=FALSE,interval=c(0,999),method="optimise",asocialVar=NULL){
  
  new("addFit",oadata=oadata,sParam=sParam,startValue=startValue,bounded=bounded,interval=interval,method=method,asocialVar=asocialVar)	
}


#Calculates likelihood for the additive model. If true ties are given then the likelihood is adjusted for these
addLikelihood<-function(l,data,asocialVar=NULL,bounded=FALSE, sParam=NULL){
  
  #Run appropriate function depending on which parameterisation for s has been specified
  if(bounded==T){
    
    likelihood<-addLikelihoodBounded(l=l,data=data,asocialVar=asocialVar, sParam= sParam);
    
  } else {
    
    likelihood<-addLikelihoodNotBounded(l=l,data=data,asocialVar=asocialVar, sParam= sParam);
    
  }
  
  if(is.character(data)){
    
    return(likelihood);
    
  }else{
    if(is.null(unlist(data@trueTies[1]))){
      return(likelihood);
    }else{
      
      for (i in 1:length(data@trueTies)){
        
        likelihood<-likelihood+addCorrectTrueTie(l=l,data=data,tiePosition=unlist(data@trueTies[i]),asocialVar=asocialVar,bounded=bounded);
        
      }
      return(likelihood);
    }
  }
}


#Calculates null likelihood for the additive model (s is fixed to zero, so only the other parameters are optimised by nlminb)
nulladdLikelihood<-function(l,data,asocialVar=NULL,bounded=FALSE){
  
  #Adds a zero to the parameter vector to correspond to no social transmission
  l<-c(0,l)
  
  likelihood<-addLikelihood(l=l,data=data,asocialVar=asocialVar,bounded=bounded)
  
  
  return(likelihood);
  
}


#Calculates likelihood for the additive model (single group) with s between 0 and 1
addLikelihoodBounded<-function(l,data,asocialVar=NULL, sParam=NULL){
  
  #Define required function
  sumWithoutNA<-function(x) sum(na.omit(x))
  
  #If the data if a character vector containing multiple oaData objects, make separate calls to addLikelihood and add together
  if(is.character(data)){		
    
    #If sParam is NULL this means only one s parameter is required
    if(is.null(sParam)) sParam<-rep(1,length(data));
    #Calculate the number of different s parameters
    noSParam<-length(levels(as.factor(sParam[sParam>0])))
    
    totalLikelihood<-0;
    
    for(i in 1:length(data)){
      
      
      if (sParam[i]==0){
        
        newl<-l[-(1:noSParam)];
        totalLikelihood<- totalLikelihood + nulladdLikelihood(l=newl,data=eval(as.name(data[i])),asocialVar=asocialVar,bounded=T);
        
      }else{
        #Build the parameter vector for the diffusion in question by selecting the correct st parameter and adding the others on the end
        newl<-c(l[sParam[i]],l[-(1:noSParam)]);
        totalLikelihood<-totalLikelihood+addLikelihoodBounded(l=newl,data=eval(as.name(data[i])),asocialVar=asocialVar);
      }
      
    }
    
    return(totalLikelihood);
    
  }else{
    
    
    #Cut down to the asocial variables specified by asocialVar
    asoc<-as.matrix(data@asoc[,asocialVar]);
    
    #Check if there are any asocial variables, set rate to one if not
    if(is.null(asocialVar)){
      solverAsocialRate<-rep(1,length(data@orderAcq));
    }else{
      #Multiply asocial variable by coefficients
      multiCoef<-as.matrix(asoc[data@orderAcq,]*t(matrix(l[-1],ncol=length(data@orderAcq),nrow=length(asocialVar))));
      #Add rows for the linear predictor for asocial learning for each solving individual
      solverLinearPred<-apply(multiCoef,1,sum);
      #Take exponentials to get rate of asocial learning
      solverAsocialRate<-exp(solverLinearPred);
    }
    
    #Now find the rate of social transmission by multiplying by the social transmission parameter
    solverSocialTransRate<-data@mldata$learnMetric*l[1];
    
    #Total rate for the solver at each stage
    solverTotalRate<-solverAsocialRate*(1-l[1])+solverSocialTransRate;
    
    #Take logs and add across acquisition events
    lComp1<-sum(log(solverTotalRate));
    
    
    #Generate a linear predictor for the rate of asocial acquisition for each ind at each time step	
    statusTracker<-rep(0,dim(data@assMatrix)[1])
    
    
    
    nsAsocialRate<-vector(length=length(data@orderAcq));
    
    for(i in 1:length(data@orderAcq)){
      
      #If there are no asocial learning variables set everyone to zero
      if(is.null(asocialVar)){
        nsLinearPred<-rep(0,dim(data@assMatrix)[1])
      }else{
        
        #Multiply asocial learning variables by coefficients
        nsMultiCoef<-as.matrix(asoc*t(matrix(l[-1],ncol=dim(data@assMatrix)[1],nrow=length(asocialVar))));
        #Add rows to get linear predictors for asocial learning
        nsLinearPred<-apply(nsMultiCoef,1,sum);
        
      }
      
      #Take exponentials and set =0 if for individuals who have already acquired the behaviour
      nsAsocialRate[i]<-sum((1-statusTracker)*exp(nsLinearPred));
      
      statusTracker[data@orderAcq[i]]<-1;
    }
    
    
    #Now find the rate of non-solver social transmission
    nsSocialTransRate<-data@mldata$totalMetric*l[1]
    #Add together
    nsTotalRate<-nsSocialTransRate+nsAsocialRate*(1-l[1]);
    
    
    #Sum the asocial and social rates, take logs and add across acquisition events
    lComp2<-sum(log(nsTotalRate));	
    #Return the likelihood for the additive model
    return(-lComp1+lComp2)
  }
}




#Calculates likelihood for the additive model with s between 0 and Inf
addLikelihoodNotBounded<-function(l,data,asocialVar=NULL,sParam=NULL){
  
  #Define required function
  sumWithoutNA<-function(x) sum(na.omit(x))
  
  #If the data is a character vector containing multiple oaData objects, make separate calls to addLikelihood and add together
  if(is.character(data)){		
    
    #If sParam is NULL this means only one s parameter is required
    if(is.null(sParam)) sParam<-rep(1,length(data));
    #Calculate the number of different s parameters
    noSParam<-length(levels(as.factor(sParam[sParam>0])))
    
    totalLikelihood<-0;
    
    for(i in 1:length(data)){
      
      
      if (sParam[i]==0){
        #The st is contstrained to zero for this diffusion
        newl<-l[-(1:noSParam)];
        totalLikelihood<- totalLikelihood + nulladdLikelihood(l=newl,data=eval(as.name(data[i])),asocialVar=asocialVar,bounded=F);
        
      }else{
        #Build the parameter vector for the diffusion in question by selecting the correct st parameter and adding the others on the end
        newl<-c(l[sParam[i]],l[-(1:noSParam)]);
        totalLikelihood<-totalLikelihood+addLikelihoodNotBounded(l=newl,data=eval(as.name(data[i])),asocialVar=asocialVar);
      }
      
    }
    
    return(totalLikelihood);
    
  }else{
    
    
    #Cut down to the asocial variables specified by asocialVar
    asoc<-as.matrix(data@asoc[,asocialVar]);
    
    #Check if there are any asocial variables, set rate to one if not
    if(is.null(asocialVar)){
      solverAsocialRate<-rep(1,length(data@orderAcq));
    }else{
      #Multiply asocial variable by coefficients
      multiCoef<-as.matrix(asoc[data@orderAcq,]*t(matrix(l[-1],ncol=length(data@orderAcq),nrow=length(asocialVar))));
      #Add rows for the linear predictor for asocial learning for each solving individual
      solverLinearPred<-apply(multiCoef,1,sum);
      #Take exponentials to get rate of asocial learning
      solverAsocialRate<-exp(solverLinearPred);
    }
    
    #Now find the rate of social transmission by multiplying by the social transmission parameter
    solverSocialTransRate<-data@mldata$learnMetric*l[1];
    
    #Total rate for the solver at each stage
    solverTotalRate<-solverAsocialRate+solverSocialTransRate;
    
    #Take logs and add across acquisition events
    lComp1<-sum(log(solverTotalRate));
    
    
    #Generate a linear predictor for the rate of asocial acquisition for each ind at each time step	
    statusTracker<-rep(0,dim(data@assMatrix)[1])
    
    
    
    nsAsocialRate<-vector(length=length(data@orderAcq));
    
    for(i in 1:length(data@orderAcq)){
      
      #If there are no asocial learning variables set everyone to zero
      if(is.null(asocialVar)){
        nsLinearPred<-rep(0,dim(data@assMatrix)[1])
      }else{
        
        #Multiply asocial learning variables by coefficients
        nsMultiCoef<-as.matrix(asoc*t(matrix(l[-1],ncol=dim(data@assMatrix)[1],nrow=length(asocialVar))));
        #Add rows to get linear predictors for asocial learning
        nsLinearPred<-apply(nsMultiCoef,1,sum);
        
      }
      
      #Take exponentials and set =0 if for individuals who have already acquired the behaviour
      nsAsocialRate[i]<-sum((1-statusTracker)*exp(nsLinearPred));
      
      statusTracker[data@orderAcq[i]]<-1;
    }
    
    
    #Now find the rate of non-solver social transmission
    nsSocialTransRate<-data@mldata$totalMetric*l[1]
    #Add together
    nsTotalRate<-nsSocialTransRate+nsAsocialRate;
    
    
    #Sum the asocial and social rates, take logs and add across acquisition events
    lComp2<-sum(log(nsTotalRate));	
    #Return the likelihood for the additive model
    return(-lComp1+lComp2)
  }
}

setMethod("summary",
          signature(object = "addFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])));
            
            
            cat("Summary of Additive Social Transmission Model\nOrder of acquisition data\n");
            if(is.null(object@optimisation$par[1])){
              if(object@bounded==T) {
                cat("Bounded parameterisation");
                sumtable<-data.frame(Estimate=object@optimisation$minimum,Unbounded=object@optimisation$minimum/(1-object@optimisation$minimum),row.names=object@varNames);
              }
              if(object@bounded==F) {
                cat("Unbounded parameterisation");
                sumtable<-data.frame(Estimate=object@optimisation$minimum,Bounded=object@optimisation$minimum/(1+object@optimisation$minimum),row.names=object@varNames);
              }
              
            }else{
              
              if(object@bounded==T) {
                cat("Bounded parameterisation");
                sumtable<-data.frame(Estimate=c(object@optimisation$par),Unbounded=c(object@optimisation$par[1:noSParam]/(1-object@optimisation$par[1:noSParam]),rep(NA,length(object@optimisation$par)-noSParam)),row.names=object@varNames);
              }
              if(object@bounded==F) {
                cat("Unbounded parameterisation");
                sumtable<-data.frame(Estimate=c(object@optimisation$par),Bounded=c(object@optimisation$par[1:noSParam]/(1+object@optimisation$par[1:noSParam]),rep(NA,length(object@optimisation$par)-noSParam)),row.names=object@varNames);
              }
              
              
            }
            
            
            if(object@bounded==T) 
              if(object@bounded==F) cat("Unbounded parameterisation");		
            
            cat("\n\nCoefficients\n");
            print(sumtable);
            cat("\n\n")
            anova(object);
            
          }
)

setMethod("anova",
          signature(object = "addFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            if(is.null(object@optimisation$par[1])){
              
              table<-data.frame(Df=c(1,0),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
              
            }else{   
              
              table<-data.frame(Df=c(length(object@optimisation$par),length(object@optimisation$par)-noSParam),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
            }
            
            
            atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine additively\n", class = c("anova", "data.frame"));
            
            
            return(atable);
            
          }
)


#Extracts the sample size for oaData in order for aicc to be calculated. As Burnham and Anderson (2002) point out, sample size is not always a straithforward issue. Here we take it to be the number of acquisition events
sampSizeExtract<-function(oadata){
  
  if(class(oadata)=="oaData"){
    
    return(length(oadata@orderAcq))
    
  }
  
  if(class(oadata)=="taData"){
    
    return(length(oadata@oadata@orderAcq))
    
  }	
  
  
  if(is.character(oadata)){		
    
    totalSS<-0;
    
    for(i in 1:length(oadata)){
      
      totalSS<-totalSS+sampSizeExtract(eval(as.name(oadata[i])));
      
    }
    
    return(totalSS);	
  }
  
}


#This function takes an oadata object and the specified position of a "true tie": where the order in which a specified subset of the diffusion chain is not known, and returns a correction to the likelihood calculated from the order given in the oadata object. If the number of tied individuals is too large there are too many permutations and the function will not work. This is the version for the additive model. There is not yet a version for the multiplicative model.
addCorrectTrueTie<-function(l,data,tiePosition,asocialVar=NULL,bounded=FALSE){
  
  status<-rep(0,length(data@orderAcq));
  if(min(tiePosition)!=1)status[data@orderAcq[1:(min(tiePosition)-1)]]<-1;
  
  #Cut down to the asocial variables specified by asocialVar
  asoc<-as.matrix(data@asoc[,asocialVar]);	
  
  #Get likelihood for the given order
  #Check if there are any asocial variables, set rate to one if not
  if(is.null(asocialVar)){
    solverAsocialRate<-rep(1,length(tiePosition));
  }else{
    #Multiply asocial variable by coefficients
    multiCoef<-as.matrix(asoc[data@orderAcq[tiePosition],]*t(matrix(l[-1],ncol=length(tiePosition),nrow=length(asocialVar))));
    #Add rows for the linear predictor for asocial learning for each solving individual
    solverLinearPred<-apply(multiCoef,1,sum);
    #Take exponentials to get rate of asocial learning
    solverAsocialRate<-exp(solverLinearPred);
  }
  
  #Now find the rate of social transmission by multiplying by the social transmission parameter
  solverSocialTransRate<-data@mldata$learnMetric[tiePosition]*l[1];
  
  #Total rate for the solver at each stage
  if (bounded==T )solverTotalRate<-solverAsocialRate*(1-l[1])+solverSocialTransRate;
  if (bounded==F )solverTotalRate<-solverAsocialRate+solverSocialTransRate;
  
  #Take logs and add across acquisition events
  lComp1<-sum(log(solverTotalRate));	
  #Generate a linear predictor for the rate of asocial acquisition for each ind at each time step	
  #Set status tracker to the start of the tied sequence
  statusTracker<-status;
  
  nsAsocialRate<-vector(length=length(tiePosition));
  
  for(i in tiePosition){
    
    #If there are no asocial learning variables set everyone to zero
    if(is.null(asocialVar)){
      nsLinearPred<-rep(0,dim(data@assMatrix)[1])
    }else{
      
      #Multiply asocial learning variables by coefficients
      nsMultiCoef<-as.matrix(asoc*t(matrix(l[-1],ncol=dim(data@assMatrix)[1],nrow=length(asocialVar))));
      #Add rows to get linear predictors for asocial learning
      nsLinearPred<-apply(nsMultiCoef,1,sum);
      
    }
    
    #Take exponentials and set =0 if for individuals who have already acquired the behaviour
    nsAsocialRate[i+1-tiePosition[1]]<-sum((1-statusTracker)*exp(nsLinearPred));
    
    statusTracker[data@orderAcq[i]]<-1;
  }
  
  
  #Now find the rate of non-solver social transmission
  nsSocialTransRate<-data@mldata$totalMetric[tiePosition]*l[1]
  #Add together
  if (bounded==T) nsTotalRate<-nsSocialTransRate+nsAsocialRate*(1-l[1]);
  if (bounded==F) nsTotalRate<-nsSocialTransRate+nsAsocialRate;	
  
  #Sum the asocial and social rates, take logs and add across acquisition events
  lComp2<-sum(log(nsTotalRate));	
  
  #Likelihood for the observed order of data within the tie
  givenOrderLogLik<--lComp1+lComp2
  
  #Now find the number of possible permutations within the tie
  
  perms<-permn(data@orderAcq[tiePosition]);
  
  
  #Now record the likelihood for each possible order within the tie
  
  likelihoodrecord<-vector(length=length(perms));
  
  for(i in 1:length(perms)){
    
    newOrderAcq<-data@orderAcq;
    newOrderAcq[tiePosition]<-unlist(perms[i]);
    
    #Get likelihood for each order
    #Check if there are any asocial variables, set rate to one if not
    if(is.null(asocialVar)){
      solverAsocialRate<-rep(1,length(tiePosition));
    }else{
      #Multiply asocial variable by coefficients
      multiCoef<-as.matrix(asoc[newOrderAcq[tiePosition],]*t(matrix(l[-1],ncol=length(tiePosition),nrow=length(asocialVar))));
      #Add rows for the linear predictor for asocial learning for each solving individual
      solverLinearPred<-apply(multiCoef,1,sum);
      #Take exponentials to get rate of asocial learning
      solverAsocialRate<-exp(solverLinearPred);
    }
    
    #Now find the rate of social transmission (for solver and total)
    
    #Define functions for calculating total st Metric
    newFunc<-function(x) x*(1-statusTracker)
    newFunc2<-function(x) x*statusTracker
    
    newLearnMetric<-vector();
    newtotalMetric<-vector();
    nsAsocialRate<-vector();
    
    
    #Set status tracker to the start of the tied sequence
    statusTracker<-status;
    
    for(j in tiePosition){
      
      newLearnMetric<-c(newLearnMetric,sum(data@assMatrix[newOrderAcq[j],]*statusTracker));
      #Calculate total st metric over all individuals in the group for each acquisition event (ML model)
      newMatrix<-apply(data@assMatrix,2,newFunc)
      newtotalMetric<-c(newtotalMetric,sum(apply(newMatrix,1,newFunc2)));
      
      #If there are no asocial learning variables set everyone to zero
      if(is.null(asocialVar)){
        nsLinearPred<-rep(0,dim(data@assMatrix)[1])
      }else{
        
        #Multiply asocial learning variables by coefficients
        nsMultiCoef<-as.matrix(asoc*t(matrix(l[-1],ncol=dim(data@assMatrix)[1],nrow=length(asocialVar))));
        #Add rows to get linear predictors for asocial learning
        nsLinearPred<-apply(nsMultiCoef,1,sum);
        
      }
      
      #Take exponentials and set =0 if for individuals who have already acquired the behaviour
      nsAsocialRate<-c(nsAsocialRate,sum((1-statusTracker)*exp(nsLinearPred)));
      
      
      statusTracker[newOrderAcq[j]]<-1;
      
    }
    
    solverSocialTransRate<-newLearnMetric*l[1];
    
    #Total rate for the solver at each stage
    if (bounded==T )solverTotalRate<-solverAsocialRate*(1-l[1])+solverSocialTransRate;
    if (bounded==F )solverTotalRate<-solverAsocialRate+solverSocialTransRate;
    
    #Take logs and add across acquisition events
    lComp1<-sum(log(solverTotalRate));	
    
    
    #Now find the rate of non-solver social transmission
    nsSocialTransRate<-newtotalMetric*l[1]
    #Add together
    if (bounded==T) nsTotalRate<-nsSocialTransRate+nsAsocialRate*(1-l[1]);
    if (bounded==F) nsTotalRate<-nsSocialTransRate+nsAsocialRate;	
    
    #Sum the asocial and social rates, take logs and add across acquisition events
    lComp2<-sum(log(nsTotalRate));	
    #Return the likelihood for the additive model
    likelihoodrecord[i]<--lComp1+lComp2		
  }
  
  #Calculate the total likelihood of the observed data within the tie
  logLikTie<--log(sum(exp(-likelihoodrecord)))
  
  #Calculate required adjustment by taking away the likelihood of the order given and adding the total likelihood of any order that results in the observed tie
  adjustment<--givenOrderLogLik+logLikTie;
  
  return(adjustment);
  
  
}


#A new class of data object for fitting time of acquisition models, similar to Franz & Nunn's NBDA

setClass("taData",representation(oadata="oaData",timeAcq="vector",endTime="numeric", nbdaData="data.frame"));

#Initialise method just takes an oaData object and an time of acquisition vector and makes a new taData object
setMethod("initialize",
          signature(.Object = "taData"),
          function (.Object, idname=NULL, oadata,timeAcq,endTime,...) 
          {
            
            timeVect<-c(0,timeAcq,endTime);
            
            nbdaData<-oadata@coxdata;
            nbdaData$time1<-timeVect[oadata@coxdata$time1+1];
            nbdaData$time2<-timeVect[oadata@coxdata$time1+2];
            
            
            temp<-oadata@coxdata[(oadata@coxdata$time2==length(oadata@orderAcq))&(oadata@coxdata$status==0),];	
            if(dim(temp)[1]>0){
              temp$time1<-timeVect[length(timeVect)-1];
              temp$time2<-timeVect[length(timeVect)];
              
              nbdaData<-rbind(nbdaData,temp);
              nbdaData$time<- nbdaData$time2-nbdaData$time1;
            }
            
            nbdaData$time<-nbdaData$time2-nbdaData$time1;
            
            #If there are tied data cut down the data appropriately
            if( sum(nbdaData$time==0)>0){
              #Cut down tied data appropriately
              rem<-(nbdaData$time==0)&nbdaData$status==1;
              statUpgrade<-rep(FALSE,dim(nbdaData)[1]);
              for(i in 1: length(nbdaData[rem,]$identity)){
                statUpgrade[as.numeric(rownames(nbdaData[nbdaData$identity==nbdaData[rem,]$identity[i],])[dim(nbdaData[nbdaData$identity==nbdaData[rem,]$identity[i],])[1]-1])]<-TRUE;
              }
              nbdaData$status[statUpgrade]<-1;
              nbdaData<-nbdaData[nbdaData$time>0,]	
            }
            
            
            
            
            callNextMethod(.Object,oadata=oadata,timeAcq=timeAcq,endTime=endTime,nbdaData=nbdaData)
            
          }
)


#Creates a new taData object taking either an oaData object and time of acquisition vector, or all components separately
taData<-function( timeAcq, endTime, idname=NULL,oadata=NULL,assMatrix=NULL,asoc=NULL,orderAcq=NULL,groupid=NULL,taskid=NULL){
  
  if(is.null(oadata)){
    
    if(is.null(taskid)) taskid<-"1";
    if(is.null(groupid)) groupid<-"1";
    oadata<-oaData(idname=idname, assMatrix=assMatrix,asoc=asoc,orderAcq=orderAcq,groupid=groupid,taskid=taskid);
    
    
  }
  
  
  
  new("taData",oadata=oadata,timeAcq=timeAcq, endTime=endTime)		
  
  
}


#Define class of object for the fitted discrete tada model
setClass("discreteTadaFit",representation(optimisation="list",additive="logical",sParam="numeric",bounded="logical",loglik="numeric",aic="numeric", aicc="numeric",nulllik="numeric",aicNull="numeric",aiccNull="numeric",LRTsocTransTS="numeric",LRTsocTransPV="numeric",asocialVar="numeric",varNames="character",task="logical",group="logical")); 


#Method for initializing object- including model fitting
setMethod("initialize",
          signature(.Object = "discreteTadaFit"),
          function (.Object, tadata,additive,sParam,startValue,asocialVar,bounded,task,group,stepLength,...) 
          {
            
            if(is.character(tadata)){		
              
              
              #If sParam is NULL this means only one s parameter is required
              if(is.null(sParam)) sParam<-rep(1,length(tadata));
              #Calculate the number of different s parameters
              noSParam<-length(levels(as.factor(sParam[sParam>0])))
              
              #For the purposes of working out how many task and group levels there are only:
              combinedData<-combineTadaData(tadata,sParam=sParam);
              sampleSize<-sum(combinedData@nbdaData$status);
              maxRate<-max(combinedData@nbdaData$time2)+1000
              
              noTasks<-length(levels(combinedData@nbdaData$task));
              noGroups<-length(levels(combinedData@nbdaData$group));
              #Calculate the number of social transmission parameters to be fitted
              noSParam<-length(levels(as.factor(sParam[sParam>0])))
              
              extraParam<-0;
              
              if(task) extraParam<-extraParam+length(levels(combinedData@nbdaData$task))-1;
              if(group) extraParam<-extraParam+length(levels(combinedData@nbdaData$group))-1;
              
            }else{
              
              sParam<-1;
              noSParam<-1;
              extraParam<-0;
              maxRate<-max(tadata@nbdaData$time2)+1000
              sampleSize<-sum(tadata@nbdaData$status);
              noTasks<-0;
              noGroups<-0;
            }
            
            #Set staring values if not specified by the user
            if(is.null(startValue)){
              startValue<-c(rep(0,noSParam),1,rep(0,length(asocialVar)+extraParam));
              startInd<-1;
            }else{
              
              startInd<-0;
              
            }
            
            if((length(asocialVar)+extraParam)==0){
              null<-optimise(nullDiscreteTadaLikelihood,interval=c(0,maxRate),asocialVar=asocialVar,tadata=tadata, additive=additive,sParam=sParam,group=group,task=task, stepLength=stepLength);
              if(startInd==1) startValue[noSParam+1]<-null$minimum;
            }else{
              null<-nlminb(start=startValue[-(1: noSParam)], nullDiscreteTadaLikelihood,lower=c(0,rep(-Inf,length(asocialVar)+extraParam)),asocialVar=asocialVar,tadata=tadata,additive=additive,sParam=sParam,group=group,task=task, stepLength=stepLength);
              if(startInd==1) startValue[-(1: noSParam)]<-null$par;
            }
            
            
            
            fit1<-nlminb(start=startValue,discreteTadaLikelihood,lower=c(rep(0,noSParam+1),rep(-Inf,length(asocialVar)+extraParam)),sParam=sParam,asocialVar=asocialVar,tadata=tadata,task=task,group=group,additive=additive, stepLength=stepLength);
            
            #Perform LRT for social transmission
            loglik<-fit1$objective;
            nulllik<-null$objective;
            LRTsocTransTS<-2*(nulllik-loglik)
            LRTsocTransPV<-1-pchisq(LRTsocTransTS,df=noSParam);
            
            if(is.null(asocialVar)) asocialVar<-NaN;
            
            
            #Calculate aic and for model without social transmission
            k<-length(fit1$par);
            aic<-(2*k)+(2*loglik);
            aicc<-aic+(2*k*(k+1)/(sampleSize-k-1));
            
            if(is.null(null$par)){k<-1}else{k<-length(null$par)};
            aicNull<-(2*k)+(2*nulllik);
            aiccNull<-aicNull+(2*k*(k+1)/(sampleSize-k-1));
            
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aic)|is.nan(aicc)|(aicc==Inf)){}else{
              if(aicc<aic) aicc<-NaN;
            }
            #To prevent a low AICc when there are more parameters than data!
            if(is.nan(aicNull)|is.nan(aiccNull)|(is.na(aiccNull))|(aiccNull==Inf)){}else{
              if(aiccNull<aicNull) aiccNull<-NaN;
            }
            
            #Extract names of variables
            if(is.nan(asocialVar[1])){
              
              varNames<-c(paste("Social Transmission",1:noSParam),"1/Rate scaling");			
              
            }else{
              varNames<-c(paste("Social Transmission",1:noSParam),"1/Rate scaling",unlist(dimnames(tadata@nbdaData)[2])[asocialVar+6]);
            }
            if(task) varNames<-c(varNames,paste("Task",levels(combinedData@nbdaData$task)[-1]));
            if(group) varNames<-c(varNames,paste("Group",levels(combinedData@nbdaData$group)[-1]));
            
            
            
            
            
            callNextMethod(.Object,optimisation=fit1,additive=additive,sParam=sParam,bounded=bounded,loglik=loglik, aic=aic,aicc=aicc,nulllik=nulllik,aicNull=aicNull,aiccNull=aiccNull,LRTsocTransTS=LRTsocTransTS,LRTsocTransPV=LRTsocTransPV,asocialVar=asocialVar,varNames=varNames,task=task,group=group,...) 
            
            
            
            
          }
)


#Function for implementing the initialization
discreteTadaFit<-function(tadata,sParam=NULL,asocialVar=NULL,startValue=NULL,bounded=FALSE,task=FALSE,group=FALSE,additive=TRUE,stepLength=1){
  
  
  if(bounded==T) {
    
    cat("Bounded parameterisation not yet available for time of acquisition models");
    return("Bounded parameterisation not yet available for time of acquisition models");	
    
  }	
  
  new("discreteTadaFit",tadata=tadata,sParam=sParam,asocialVar=asocialVar,startValue=startValue,bounded=bounded,task=task,group=group,additive=additive,stepLength=stepLength)	
}


nullDiscreteTadaLikelihood<-function(l,tadata,sParam=NULL,bounded=FALSE,asocialVar=NULL,task=FALSE,group=FALSE,additive=TRUE,stepLength=1){
  
  if(is.null(sParam)){noSParam<-1}else{noSParam<-length(levels(as.factor(sParam[sParam>0])))};
  newl<-c(rep(0,noSParam),l);
  return(discreteTadaLikelihood(newl,tadata,sParam=sParam,bounded=bounded,asocialVar=asocialVar,task=task,group=group, additive=additive,stepLength=stepLength));
}


#Gets likelihood for a multiplicative nbda model 
discreteTadaLikelihood<-function(l,tadata,sParam=NULL,bounded=FALSE,asocialVar=NULL,task=FALSE,group=FALSE,additive=TRUE,stepLength=1){
  
  #If the data is a character vector containing multiple oaData objects, make separate calls to addLikelihood and add together
  
  
  if(is.character(tadata)){		
    
    
    #If sParam is NULL this means only one s parameter is required
    if(is.null(sParam)) sParam<-rep(1,length(tadata));
    #Calculate the number of different s parameters
    noSParam<-length(levels(as.factor(sParam[sParam>0])))
    
    #For the purposes of working out how many task and group levels there are only:
    combinedData<-combineTadaData(tadata,sParam=sParam);
    sampleSize<-sum(combinedData@nbdaData$status);
    
    noTasks<-length(levels(combinedData@nbdaData$task));
    noGroups<-length(levels(combinedData@nbdaData$group));
    
    sParamVect<-l[1:noSParam];
    rateScale<-l[(noSParam+1)];
    if(is.null(asocialVar)){
      asocVarVect<-NULL;
      asocLength<-0;
    }else{
      asocLength<-length(asocialVar);
      asocVarVect<-l[(noSParam+2):(noSParam+2+ asocLength)];
    }
    lengthTask<-0
    taskParam<-NULL
    if(task){
      lengthTask<-noTasks-1
      if(lengthTask>0){
        taskParam<-c(0,l[(noSParam+2+asocLength):(noSParam+2+asocLength+lengthTask)])
      }
    }
    lengthGroup<-0
    groupParam<-NULL
    if(group){
      lengthGroup<-noGroups-1
      if(lengthGroup>0){
        groupParam<-c(0,l[(noSParam+2+asocLength+lengthTask):(noSParam+2+asocLength+lengthTask+lengthGroup-1)])
      }
    }
    
    totalLikelihood<-0;
    
    
    
    
    
    
    for(i in 1:length(tadata)){
      
      
      if(length(stepLength)==1) stepLength<-rep(stepLength,eval(as.name(tadata[i]))@endTime);
      
      if(dim(as.matrix(stepLength))[2]==1) {indStepLength<-as.matrix(stepLength)[,1]} else{indStepLength<-as.matrix(stepLength)[,i]}
      
      if (sParam[i]==0){
        #The st is contstrained to zero for this diffusion
        
        newl<-c(0,rateScale,asocVarVect, taskParam[which(eval(as.name(tadata[i]))@oadata@taskid==levels(combinedData@nbdaData$task))], groupParam[which(eval(as.name(tadata[i]))@oadata@groupid==levels(combinedData@nbdaData$group))])
        
        totalLikelihood<- totalLikelihood + discreteTadaLikelihood(l=newl,tadata=eval(as.name(tadata[i])),asocialVar=asocialVar, task=task, group=group, additive=additive, stepLength=stepLength);
        
      }else{
        #Build the parameter vector for the diffusion in question by selecting the correct st parameter and adding the others on the end
        newl<-c(l[sParam[i]],rateScale,asocVarVect, taskParam[which(eval(as.name(tadata[i]))@oadata@taskid==levels(combinedData@nbdaData$task))], groupParam[which(eval(as.name(tadata[i]))@oadata@groupid==levels(combinedData@nbdaData$group))])
        
        totalLikelihood<- totalLikelihood + discreteTadaLikelihood(l=newl,tadata=eval(as.name(tadata[i])),asocialVar=asocialVar, task=task, group=group, additive=additive, stepLength=stepLength);
      }
      
    }
    
    return(totalLikelihood);
    
  }else{
    
    
    #Cut down to the asocial variables specified by asocialVar
    asoc<-as.matrix(tadata@oadata@asoc[,asocialVar]);
    
    #If there are group and/or task effects, store them and take them off the end of the parameter vector
    if(group==TRUE) {
      groupEffect<-l[length(l)];
      l<-l[-length(l)];
    }else{groupEffect<-0}
    if(task==TRUE) {
      taskEffect<-l[length(l)];
      l<-l[-length(l)];
    }else{taskEffect<-0}
    
    
    #Check if there are any asocial variables, set rate to one if not
    if(is.null(asocialVar)){
      linearPred<-rep(0,dim(tadata@oadata@assMatrix)[1]);
    }else{
      #Multiply asocial variable by coefficients and sum for each individual
      linearPred<-apply(as.matrix(asoc*t(matrix(l[-(1:2)],ncol=dim(asoc)[1],nrow=length(asocialVar)))),1,sum);
    }
    #Take exponentials to get rate of asocial learning
    asocialRate<-exp(linearPred+groupEffect + taskEffect);
    
    timeSolve<-rep(tadata@endTime+1,length=dim(tadata@oadata@assMatrix)[1])
    
    for(i in 1:length(tadata@oadata@orderAcq)){
      
      timeSolve[tadata@oadata@orderAcq[i]]<-tadata@timeAcq[i]
      
    }
    
    
    if(length(stepLength)==1) stepLength<-rep(stepLength,tadata@endTime)
    
    logLikelihood<-0
    statusTracker<-rep(0,dim(tadata@oadata@assMatrix)[1])
    i<-0
    statusTracker[which(timeSolve==i)]<-1;	
    
    #Cycle through timeSteps, updating log likelihood and status of individuals
    for(i in 1:tadata@endTime){
      socialRate<-apply((tadata@oadata@assMatrix*statusTracker),2,sum)*l[1];
      if(additive){
        stepRate<-((1/l[2])*(socialRate+asocialRate)*stepLength[i]);
      }else{
        stepRate<-((1/l[2])*((socialRate+1)*asocialRate)*stepLength[i]);
      }
      logLikelihood<-logLikelihood +sum(((1-1*(timeSolve==i))*(stepRate)-log(1-exp(-stepRate))*(timeSolve==i))*(1-statusTracker));
      statusTracker[which(timeSolve==i)]<-1;
    }
    
    return(logLikelihood)
    
  }
}


setMethod("anova",
          signature(object = "discreteTadaFit"),
          function (object, ...) 
          {
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            if(is.null(object@optimisation$par[1])){
              
              table<-data.frame(Df=c(1,0),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
              
            }else{   
              
              table<-data.frame(Df=c(length(object@optimisation$par),length(object@optimisation$par)-noSParam),LogLik=c(object@loglik,object@nulllik),AIC=c(object@aic,object@aicNull),AICc=c(object@aicc,object@aiccNull),LR=c(object@LRTsocTransTS,NA),p=c(object@LRTsocTransPV,NA),row.names=c("With Social Transmission","Without Social Transmission"));
              
            }
            
            if(object@additive){atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine additively\n", class = c("anova", "data.frame"))}else{
              atable<-structure(table, heading = "Likelihood Ratio Test for Social Transmission:\n\nNull model includes all other specified variables\nSocial transmission and asocial learning assumed to combine multiplicatively\n", class = c("anova", "data.frame"))			
            }
            
            return(atable);
            
          }
)


setMethod("summary",
          signature(object = "discreteTadaFit"),
          function (object, ...) 
          {
            
            if(object@additive==T)cat("Summary of Additive Social Transmission Model\nDiscrete Time of Acquisition Data\n");
            if(object@additive==F)cat("Summary of Multiplicative Social Transmission Model\nDiscrete Time of Acquisition Data\n");
            
            
            noSParam<-length(levels(as.factor(object@sParam[object@sParam>0])))
            
            
            
            if(object@bounded==F) {
              cat("Unbounded parameterisation.");
              sumtable<-data.frame(Estimate=c(object@optimisation$par),Bounded=c(object@optimisation$par[1: noSParam]/(1+object@optimisation$par[1: noSParam]),rep("",length(object@optimisation$par)-noSParam)),row.names=object@varNames);   
              
            }		
            
            
            
            
            cat("\n\nCoefficients\n");
            print(sumtable);
            cat("\n\n")
            anova(object);
            
          }
)
