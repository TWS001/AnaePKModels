#' @useDynLib AnaePKModels PropSchnider PropSchnider2 PropEleveldFinal
#' @title Propofol Schnider
#'
#' @description V1..V3,Cl1..Cl3,ke0
#' @param age Age in years
#' @param wt Weight (kg)
#' @param ht Height(cm)
#' @param male Gender (male=1)
#' @export
#' @examples
#' propSchnider(20,70,172,1)


propSchnider <- function(age,wt,ht,male){
  .C("PropSchnider",
     as.double(age),
     as.double(wt),
     as.double(ht),
     as.logical(male),
    res=double(7))$res
}

#' @title Propofol Schnider (2)
#'
#' @description V1..V3,Cl1..Cl3,ke0
#' @param age Age in years
#' @param wt Weight (kg)
#' @param ht Height(cm)
#' @param lbm Lean body mass (James)
#' @keywords models
#' @export
#' @examples
#' propSchnider(20,70,172,50.6)

propSchnider2 <- function(age,wt,ht,lbm){
  .C("PropSchnider2",
     as.double(age),
     as.double(wt),
     as.double(ht),
     as.double(lbm),
     res=double(7))$res
}


#' @title Propofol, Diprifusor
#'
#' @description V1..V3,Cl1..Cl3,ke0
#' @param age Age in years
#' @param wt Weight
#' @keywords models
#' @export
#' @examples
#' diprifusor(20,76)

diprifusor <- function(age,wt){
  V1 = wt * 0.228
  if(age<16) V1 = 0.0
  V2 = V1*(0.114/0.055)
  V3 = V1*(0.0419/0.0033)
  Cl1 = V1*0.119
  Cl2 = V2*0.055
  Cl3 = V3*0.0033
  ke0 = 0.26;
  c(V1,V2,V3,Cl1,Cl2,Cl3,ke0)
}

#' @title Propofol, Marsh
#'
#' @description V1..V3,Cl1..Cl3, dummy ke0 of 10
#' @param wt Weight
#' @keywords models
#' @export
#' @examples
#' propMarsh(20)

propMarsh <- function(wt){
  V1 = wt * 0.228
  V2 = V1*(0.112/0.055)
  V3 = V1*(0.0419/0.0033)
  Cl1 = V1*0.119
  Cl2 = V2*0.055
  Cl3 = V3*0.0033
  ke0 = 20;
  c(V1,V2,V3,Cl1,Cl2,Cl3,ke0)
}



#' @title Remifentanil(incl. ke0), Minto et al. (Age, LBM) Anesthesiology 1997; 86:10-23.
#'
#' @description V1..V3,Cl1..Cl3,ke0
#' @param age Age in years
#' @param lbm Lean body mass (James)
#' @keywords models
#' @export
#' @examples
#' remiMinto(20,50.6)



# Remifentanil(incl. ke0), Minto et al. (Age, LBM) Anesthesiology 1997; 86:10-23.
remiMinto<-function(age, lbm)
    {
    V1 = 5.1-(0.0201 * (age-40))+(0.072*(lbm-55))
    V2 = 9.82-(0.0811 * (age-40))+(0.108*(lbm-55))
    V3 = 5.42
    Cl1= 2.6-(0.0162 * (age-40))+(0.0191 * (lbm-55))
    Cl2 = 2.05-(0.0301 * (age-40))
    Cl3 = 0.076-(0.00113 * (age-40))
    ke0 = 0.595-0.007 * (age-40)
    c(V1,V2,V3,Cl1,Cl2,Cl3,ke0)
    }

#' @title Propofol, Eleveld BJA 2018
#'
#' @description V1..V3,Cl1..Cl3,ke0
#' @param wt Weight
#' @param ht Height
#' @param age Age in years
#' @param male Gender, male?
#' @param pma postmenstrual age
#' @param venous Venous concentration?
#' @param opiates Are opiates used?
#' @param patient Is it a patient?
#' @param tpeak Time to peak effect site concentration (1.6 min)
#' @keywords models
#' @export
#' @examples
#' propEleveld(70,172,40,1)

propEleveld <- function(wt,ht,age,male, pma=40,venous=0,opiates=1,patient=1,tpeak = 0 ){
  ref<- "Eleveld 2018 to be completed"
  tempPK<-.C("PropEleveldFinal",
     as.double(wt),
     as.double(ht),
     as.double(pma),
     as.double(age),
     as.logical(male),
     as.logical(venous),
     as.logical(opiates),
     as.logical(patient),
     res=double(7))$res
     if(tpeak == 0)
       pk<-tempPK
    else{
     tempPK[7]<-PKPDTools::tpeak2ke0(PKPDTools::vcl2hyb(tempPK[1:6]),tpeak)
     pk<-tempPK
    }
  var<-NULL
  list(pk=pk,var=var,ref=ref)
}

#' @title Alfentanil, Maitre 1987 as tested 1988 Anesthesiology
#'
#' @description V1..V3,Cl1..Cl3
#' @param age Age in years
#' @param wt Weight
#' @param MALE Gender (MALE=1)
#' @keywords models
#' @export
#' @examples
#' alfentMaitre(20,76,1)

alfentMaitre <- function(age,wt,MALE){
  Cl1 = 0.356
  V1 = 0.111*wt
  k31 = 0.0126

  if(age > 40){
    Cl1 = Cl1-(0.00269*(age-40))
    k31= k31 -(0.000113*(age-40))
  }

  if(!MALE)
    V1= V1*1.15

  V2 = V1*1.545319
  V3 = V1*0.017/k31
  Cl2 = V2*0.0673
  Cl3 = V3*k31

  c(V1,V2,V3,Cl1,Cl2,Cl3)
}

#' @title Alfentanil, Scott/Stanski JPET 1987 with ke0
#' @description V1..V3,Cl1..Cl3, ke0
#' @keywords models
#' @export
#' @examples
#' alfentScott()


alfentScott <- function(){
  V1 = 2.19
  V2 = 5.915888
  V3 = 13.17647
  Cl1 = 0.195
  Cl2 = 1.266
  Cl3 = 0.224
  ke0 = 0.770164

  c(V1,V2,V3,Cl1,Cl2,Cl3,ke0)
}

#' @title Scott Fentanyl weight scaled with ke0 0.147
#'
#' @description V1..V3,Cl1..Cl3, ke0
#' @param wt Weight
#' @keywords models
#' @export
#' @examples
#' fentaScott(76)


fentaScott <- function(wt){
  V1 = wt*0.1814286
  V2 = wt*0.7049256
  V3 = wt*4.241187
  Cl1 = wt*0.01016
  Cl2 = wt*0.06767286
  Cl3 = wt*0.03265714
  ke0 = 0.147

  c(V1,V2,V3,Cl1,Cl2,Cl3,ke0)
}


