library(dplyr)
library(MASS)

set.seed(17022024)

# parameter ----
n.ind=200
ind_rep=5
n.groups=10

# Fixed effect
mu1 <- 50
mu2 <- 0

beta1 <- 1.7
beta2 <- 0.9
beta3 <- -0.8

# means
M_I<-c(0,0,0,0)
M_Gr<-c(0,0,0,0)

# variances
V_x1 <- 1 
V_int1   <- 0.3    # variance direct effect (mean behaviour)
V_slope1 <- 0.1
V_int2   <- 0.35    # variance direct effect (mean behaviour)
V_slope2 <- 0.12

# covariances
C_int1_slope1   <- 0.3*(sqrt(V_int1*V_slope1))         
C_int1_int2     <- -0.5*(sqrt(V_int1*V_int2)) 
C_int1_slope2   <- -0.3*(sqrt(V_int1*V_slope2))         
C_int2_slope1   <- 0*(sqrt(V_int2*V_slope1))
C_slope1_slope2 <- 0*(sqrt(V_slope1*V_slope2))
C_int2_slope2   <- -0.2*(sqrt(V_int2*V_slope2)) 

# Variance I-matrix
I<-matrix(NA,4,4)
I[1,1]<-V_int1
I[2,2]<-V_slope1
I[3,3]<-V_int2
I[4,4]<-V_slope2

I[1,2]<-I[2,1]<-C_int1_slope1
I[1,3]<-I[3,1]<-C_int1_int2
I[1,4]<-I[4,1]<-C_int1_slope2
I[2,3]<-I[3,2]<-C_int2_slope1
I[2,4]<-I[4,2]<-C_slope1_slope2
I[3,4]<-I[4,3]<-C_int2_slope2


## Group varcov
V_Gr_int1   <- 0.1    # variance direct effect (mean behaviour)
V_Gr_slope1 <- 0.05
V_Gr_int2   <- 0.1    # variance direct effect (mean behaviour)
V_Gr_slope2 <- 0.05

# covariances
C_int1_slope1   <- 0.6*(sqrt(V_Gr_int1*V_Gr_slope1))         
C_int1_int2     <- 0*(sqrt(V_Gr_int1*V_Gr_slope1)) 
C_int1_slope2   <- 0*(sqrt(V_Gr_int1*V_Gr_slope1))         
C_int2_slope1   <- 0*(sqrt(V_Gr_int1*V_Gr_slope1))
C_slope1_slope2 <- 0*(sqrt(V_Gr_int1*V_Gr_slope1))
C_int2_slope2   <- 0.2*(sqrt(V_Gr_int1*V_Gr_slope1)) 

# Variance I-matrix
GR<-matrix(NA,4,4)
GR[1,1]<-V_Gr_int1
GR[2,2]<-V_Gr_slope1
GR[3,3]<-V_Gr_int2
GR[4,4]<-V_Gr_slope2

GR[1,2]<-GR[2,1]<-C_int1_slope1
GR[1,3]<-GR[3,1]<-C_int1_int2
GR[1,4]<-GR[4,1]<-C_int1_slope2
GR[2,3]<-GR[3,2]<-C_int2_slope1
GR[2,4]<-GR[4,2]<-C_slope1_slope2
GR[3,4]<-GR[4,3]<-C_int2_slope2

Ve1 = 0.7
Ve2 = 0.3

#Random effects ftness
Vw_I = 0.5
Vw_e = 0.4

#Fixed effects on fitness
mu=5
betaw_int=1
betaw_int_q=-0.4

betaw_slope=0.5
betaw_slope_q=0

betaw_int_c=1

betaw_partner=0.1
betaw_partner_q=0
betaw_partner_c=-1.1


# study design ----
  IDi <- rep(1:n.ind, ind_rep)
  IDj <- rep(1:n.ind, ind_rep)+1
  IDj <- ifelse(IDj > n.ind, IDj-n.ind, IDj)
  df <- data.frame(IDi,IDj)
  df$Population <- rep(1:n.groups, 10)
  

# sim func ----
###Create phenotypes input needs output of samp.design
  blups_I=as.data.frame(mvrnorm(n.ind, M_I, I))
  colnames(blups_I)<-c("int1", "slope1", "int2", "slope2")
  blups_I$ID=1:n.ind
  
  df$int1=blups_I[match(df$IDi,blups_I$ID),"int1"]
  df$slope1=blups_I[match(df$IDi,blups_I$ID),"slope1"]
  df$int2=blups_I[match(df$IDi,blups_I$ID),"int1"]
  df$slope2=blups_I[match(df$IDi,blups_I$ID),"slope2"]
  
  df$int1_partner <- blups_I[match(df$IDj,blups_I$ID),"int1"]
  
  blups_Group=as.data.frame(mvrnorm(n.ind, M_Gr, GR))
  colnames(blups_Group)<-c("int1", "slope1", "int2", "slope2")
  blups_Group$Group=1:n.groups
  
  df$Gr_int1=blups_Group[match(df$Population,blups_Group$Group),"int1"]
  df$Gr_slope1=blups_Group[match(df$Population,blups_Group$Group),"slope1"]
  df$Gr_int2=blups_Group[match(df$Population,blups_Group$Group),"int1"]
  df$Gr_slope2=blups_Group[match(df$Population,blups_Group$Group),"slope2"]
  
  df$e1<-rnorm(nrow(df),0, sqrt(Ve1))
  df$e2<-rnorm(nrow(df),0, sqrt(Ve2))
  
  W_blups_I <-as.data.frame(rnorm(n.ind,0, sqrt(Vw_I)))
  colnames(W_blups_I) <- "W_int"
  W_blups_I$ID=1:n.ind
  df$W_int=W_blups_I[match(df$IDi,W_blups_I$ID),"W_int"]
  df$Vw_e<-rnorm(nrow(df),0, sqrt(Vw_e))
  
  df$x1 <- rnorm(nrow(df),10, sqrt(V_x1))
  df$x2 <- rep(-0.5:0.5, each=5, 100)
  
 ## Phenotypic Equation
  
 # Phenotype results from an individual's breeding value
 # and the population response psi to a known fixed phenotype (Xj) 
 # and the social breeding value of j
  df$z1=mu1 + beta1*scale(df$x1)[,] + beta2*df$x2 + beta3*scale(df$x1)[,]*df$x2 + df$int1 + df$slope1*df$x1 + df$Gr_int1 + df$Gr_slope1*scale(df$x1)[,] + df$e1
  df$z2=mu2 + beta1*scale(df$x1)[,] + beta2*df$x2 + beta3*scale(df$x1)[,] *df$x2 +  df$int2 + df$slope2*scale(df$x1)[,] + df$Gr_int2 + df$Gr_slope2*scale(df$x1)[,] + df$e2
  
  df$w1=mu + betaw_int*df$int1 + betaw_int_q*df$int1*df$int1 + betaw_slope*df$slope1 + 
        betaw_slope_q*df$slope1*df$slope1 + betaw_int_c*df$int1*df$slope1 +
        betaw_partner*df$int1_partner + betaw_partner_q*df$int1_partner +
        betaw_partner_c*df$int1*df$int1_partner
  
  df$w2=mu + df$W_int + betaw_int*df$int1 + betaw_int_q*df$int1*df$int1 + betaw_slope*df$slope1 + 
    betaw_slope_q*df$slope1*df$slope1 + betaw_int_c*df$int1*df$slope1 +
    betaw_partner*df$int1_partner + betaw_partner_q*df$int1_partner +
    betaw_partner_c*df$int1*df$int1_partner + df$Vw_e
  
## Rename to suit behavioural ecology group
data<-df[,c("IDi", "IDj", "Population", "x1", "x2", "z1", "z2", "w1","w2")]
data$x2 <- ifelse(data$x2==0.5, "M", "F")
colnames(data)<- c("Individual", "Partner", "Population", "Density", "Sex", "Exploration", "Aggression", "LRS", "Feeding_rate")
data$Density <- round(data$Density, 1)
data$Aggression_score <- rpois(1000, lambda = exp((data$Aggression - min(data$Aggression)) / 60))
data$Population <- as.factor(data$Population)
levels(data$Population)<- c("Bagley Wood", "Boshoek", "Donana", "Forstenrieder park", "Gotland", "Mont Ventoux", "Oulu", "Westerheide", "Whytham Woods", "Zvenigorod") 
data$Individual <- as.factor(data$Individual)
data$Individual <- paste0("C3F7", sprintf("%03d", data$Individual))
data$Partner <-paste0("C3F7", sprintf("%03d", data$Partner))

# write data
write.csv(data, "fake_SPI_birds.csv", row.names = F)

# Group by Individual and calculate means
summary_data <- data %>%
  group_by(Individual) %>% summarise(
    Mean_Exploration = mean(Exploration, na.rm = TRUE),
    Mean_Aggression = mean(Aggression, na.rm = TRUE),
    Mean_LRS = mean(LRS, na.rm = TRUE))
    

ggplot(summary_data, aes(x=Mean_Aggression, y=Mean_Exploration))+
  geom_point() +
  scale_x_continuous(name="Individual mean aggression") +
  scale_y_continuous(name="Individual mean exploration") +
  scale_color_discrete(guide="none") +
  theme_classic()

ggplot(summary_data, aes(x=Mean_Exploration, y=Mean_LRS))+
  geom_point() +
  scale_x_continuous(name="Individual mean exploration") +
  scale_y_continuous(name="Lifetime reproductive success") +
  scale_color_discrete(guide="none") +
  theme_classic()

ggplot(summary_data, aes(x=blups_I$int1, y=Mean_LRS))+
  geom_point() +
  scale_x_continuous(name="Individual mean exploration") +
  scale_y_continuous(name="Lifetime reproductive success") +
  scale_color_discrete(guide="none") +
  theme_classic()

ggplot(summary_data, aes(x=Mean_Exploration, y=Mean_LRS))+
  geom_point() +
  scale_x_continuous(name="Individual mean exploration") +
  scale_y_continuous(name="Lifetime reproductive success") +
  scale_color_discrete(guide="none") +
  theme_classic()
