
#' @title Enoxaparin, Br J Clin Pharmacol 2003; 56, 96–103
#'
#' @description Reference time = hour!
#' @param wt Weight
#' @param lbm Lean body mass (James)
#' @keywords models
#' @export
#' @examples
#' enoxaparin(70,45)


enoxaparin <- function(wt,lbm){
  ref<- "Br J Clin Pharmacol 2003; 56, 96 - 103"
  V1 = 3.67*wt/70
  V2 = 13.6
  V3 = 0
  Cl1 = lbm*1.03/70
  Cl2 = 0.363
  Cl3 = 0
  kA = 0.195;
  pk<-c(V1,V2,V3,Cl1,Cl2,Cl3,kA)
  CVV1 = 24.5
  CVV2 = 34.4
  CVV3 = 0
  CVCl1 = 6.8
  CVCl2 = 33.3
  CVCl3 = 0
  CVka = 25.6
  sigma = 30
  var<-c(CVV1,CVV2,CVV3,CVCl1,CVCl2,CVCl3,sigma)
  list(pk=pk,var=var,ref=ref)
}

#' @title Enoxaparin, Br J Clin Pharmacol. 2007 Oct; 64(4): 428–438.
#'
#' @description Reference time = hour!
#' @param age Age
#' @param wt Weight
#' @param male Gender, male?
#' @param sCr Serum creatinin
#' @keywords models
#' @export
#' @examples
#' enoxaparin2(60,70,1,34)


enoxaparin2<-function(age,wt,male,sCr){
  ref<- "Br J Clin Pharmacol. 2007 Oct; 64(4): 428-438."
  V1 = 6.43 *(wt/65)^1.25
  V2= 8.18
  V3 = 0
  Cl1 = 0.7*(wt/65)^0.78 * (gfrMDRD(age,sCr,male)/69)^0.25
  Cl2 = 0.34
  Cl3 = 0
  kA = 0.63
  pk<-c(V1,V2,V3,Cl1,Cl2,Cl3,kA)
  list(pk=pk,var=NULL,ref=ref)
}

#' @title Enoxaparin,  Br J Clin Pharmacol. 2006 Aug;62(2):165-76.
#'
#' @description Reference time = hour!
#' @param age Age
#' @param wt Weight
#' @param IBW Ideal body weight
#' @param male Gender, male?
#' @param sCr Serum creatinin
#' @keywords models
#' @export
#' @examples
#' enoxaparin3(80,70,63,1,34)


enoxaparin3<-function(age,wt,IBW,male,sCr){
  ref<- "Br J Clin Pharmacol. 2006 Aug;62(2):165-76."
  V1 = wt*6.78/70
  V2= 6.19
  V3 = 0
  Cl1 = 0.229 + (0.744*(60*gfrCoGa(age,IBW,male,sCr)/4800))
  Cl2 = 0.429
  Cl3 = 0
  kA = 0.476
  pk<-c(V1,V2,V3,Cl1,Cl2,Cl3,kA)
  list(pk=pk,var=NULL,ref=ref)
}




# https://www.niddk.nih.gov/health-information/communication-programs/nkdep/laboratory-evaluation/glomerular-filtration-rate-calculators/mdrd-adults-conventional-units
gfrMDRD<-function(age,sCr,male,african = 0){
  # Creatinin in mg/dl -> Umrechnen von microMol/l
  sCr<-sCr/88.4
   tmp<-175 * (sCr^-1.154) * (age^-0.203)
  if(!male){
    tmp<-tmp*0.742
  }
  if(african){
    tmp<-tmp*1.212
  }
  return(tmp)
}
# https://www.researchgate.net/publication/21910265_Cockcroft_DW_Gault_MHPrediction_of_creatine_clearance_from_serum_creatinine_Nephron_16_31-41
gfrCoGa<-function(age,wt,male,sCr){
  # Creatinin in mg/dl -> Umrechnen von microMol/l
  sCr<-sCr/88.4
  f=1
  if(!male)f=0.85
  return ((f*(140-age)*wt)/(72*sCr))
}
